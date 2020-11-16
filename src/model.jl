# for a single target tree
# for a single orthoDist
# for a single cell (bc::DataFrameRow) do:
# generate a 3-tuple of ligand::Vector,pathgene::Matrix,target::Matrix collections
function getPathSignals(g::AbstractMetaGraph,dist::Int,bc::DataFrameRow)
        ag = addExpression(g,bc,:orthoHGNC,:rnaSeq);
	colnames = [:vertex,:hgnc,:exp]
	coltypes = [Int64,String,Float64]
	expTable = DataFrame(coltypes,colnames)
        for v in 1:nv(ag)
		if haskey(props(ag,v),:orthoDist)
			odi = props(ag,v)[:orthoDist][dist][:members]
			for m in odi
		                if haskey(m,:rnaSeq)
		                        acc = m[:rnaSeq][:accession]
		        		val = m[:rnaSeq][:value]
		        		# println("vert:$v - gene:$acc is $val")
					data = tuple(v,acc,val)
					push!(expTable,initRow(colnames,colnames,data));
		                end
			end
		end
        end
	(ag=ag,exp=expTable)
end

# for a single target tree
# for multiple cells
# get tuple: lig, pathgene, target
function getCellTrees(targTree::AbstractMetaGraph,dist::Int,bcs::DataFrame)
	df = DataFrame();
	for bc in eachrow(bcs)
		bcode = bc.barcode
		hpf = bc.hpf
		bctype = bc.celltype
		bcnctype = bc.labels
		println("getting bc $bcode");
		ag = getPathSignals(targTree,dist,bc);
		insertcols!(ag.exp,1,(:hpf=>fill(hpf,size(ag.exp,1))))
		insertcols!(ag.exp,1,(:barcode=>fill(bcode,size(ag.exp,1))))
		insertcols!(ag.exp,1,(:subtype=>fill(bcnctype,size(ag.exp,1))))
		insertcols!(ag.exp,1,(:type=>fill(bctype,size(ag.exp,1))))
		append!(df,ag.exp);
	end
	df
end

# for a single target tree
# for multiple cells
# get tuple: lig, pathgene, target
function getCellTreesThread(targTree::AbstractMetaGraph,dist::Int,bcs::DataFrame)
	dat = Vector{Vector{DataFrame}}()
	gr = Vector{Vector{Dict}}()
	for i in 1:Threads.nthreads()
	  push!(dat,Vector{DataFrame}())
	  push!(gr,Vector{Dict}())
	end
	Threads.@threads for bc in eachrow(bcs)
		bcode = bc.barcode
		hpf = bc.hpf
		bctype = bc.celltype
		bcnctype = bc.labels
		println("getting bc $bcode");
		ag = getPathSignals(targTree,dist,bc);
		insertcols!(ag.exp,1,(:hpf=>fill(hpf,size(ag.exp,1))))
		insertcols!(ag.exp,1,(:barcode=>fill(bcode,size(ag.exp,1))))
		insertcols!(ag.exp,1,(:subtype=>fill(bcnctype,size(ag.exp,1))))
		insertcols!(ag.exp,1,(:type=>fill(bctype,size(ag.exp,1))))
		push!(dat[Threads.threadid()],ag.exp)
		push!(gr[Threads.threadid()],Dict(Symbol(bcode)=>ag.ag))
	end
	dat = reduce(vcat,reduce(vcat,dat))
	gr = reduce(vcat,reduce(vcat,gr))
	return (df=dat,g=gr)
end

# filter by cell type and a vector of attributes (each defines a separate output table)
# find |featInds| barcode/gene/expression matrices from
# a cell type/subtype filter dict
# variable number of featInds vectors of valid vertex indices
function filterTypes(df::DataFrame,T::Dict,featInds::Vector...)
	arrT = []
	for i in 1:length(featInds)
		if !ismissing(first(values(T))) # if the subtype is present, filter type and subtype
			println("filtering on type=",first(keys(T)),",subtype=",first(values(T)))
			df_type = filter(row->row.type==first(keys(T)) &&
						  row.subtype==first(values(T)),df)
		else # if the subtype is not present, filter type only
			println("filtering on type=",first(keys(T)))
			df_type = filter(row->row.type==first(keys(T)),df)
		end
		df_feat = filter(row->any(in.(row.vertex,featInds[i])),df_type)
		# make sure barcode/type/gene is unique (eg by ignoring when the same gene comes up in different complexes, generating non-unique expression values for barcode/type/vertex/exp )
		select!(df_feat, Not(:vertex))
		unique!(df_feat)
		push!(arrT,df_feat)
	end
	arrT
end

# get the product of all populations in arrT for the features in arrFeat
# for each type in arrT, a corresponding set of features in arrFeat is pulled
function getMultiObs(df_orig::DataFrame,arrT::Vector,arrFeat::Vector)
	df = df_orig[findall(x->x==0,nonunique(df_orig)),:]

	arrTF = []
	for t in 1:length(arrT)
		arrT_i_F = filterTypes(df,arrT[t],[arrFeat[t]])
		push!(arrTF,arrT_i_F...)
	end

	arrTFpiv = unstack.(arrTF,:barcode,:hgnc,:exp) # for each barcode, get genes as columns and expression as values
	arrTFpiv = map(df->select(df,Not(:barcode)),arrTFpiv) # remove the barcode column, leaving only a matrix of expression
	arrTFcross = crossjoin(arrTFpiv...,makeunique=true) # generate single cell samples of the cluster-wise state by computing the product
	# Matrix{Float64}(arrTFcross)
	# arrTFcross
end

# uniform sample from product of all populations in arrT for the features in arrFeat
# for each type in arrT, a corresponding set of features in arrFeat is pulled
function sampleMultiObs(df_orig::DataFrame,arrT::Vector,arrFeat::Vector,nobs=Int)
	df = df_orig[findall(x->x==0,nonunique(df_orig)),:]

	arrTF = []
	for t in 1:length(arrT)
		arrT_i_F = filterTypes(df,arrT[t],[arrFeat[t]])
		push!(arrTF,arrT_i_F...)
	end
	arrTFpiv = unstack.(arrTF,:barcode,:hgnc,:exp) # for each barcode, get genes as columns and expression as values
	arrTFpiv = map(df->select(df,Not(:barcode)),arrTFpiv) # remove the barcode column, leaving only a matrix of expression

	sizes = map(o->size(o)[1],arrTFpiv)
	indices = StatsBase.sample.(map(x->[1:x;],sizes),nobs)
	d = reduce(+,map(o->size(o)[2],arrTFpiv))
	mat = Array{Float64,2}(undef,nobs,d)
	for i in 1:nobs
		matrow = []
		for j in 1:length(arrTFpiv)
			r = indices[j][i]
			# println(arrTFpiv[j][r,:])
			matrow = vcat(matrow,convert(Array,arrTFpiv[j][r,:]))
			# println(length(matrow))

		end
		mat[i,:] .= vec(matrow)
	end
	mat
end

# condition and return a Turing.jl model on the data and training parameters
function singlePathPairedModel(targT1,ligT1,ligT2,pathT1,σ_pair,λ²_pair)
	yTrain = targT1'
	xTrain = [ligT1';pathT1';ligT2';pathT1']

	@model function grouped_lasso(y, X, σ,λ²)
		# number of observations and features
		p, nobs = size(X)
	    mk = div(p, 2)

		# set variance prior (shrinkage of the group-wise linear coefficients)
		τ² ~ filldist(Gamma((mk + 1) / 2, 2 / λ²), 2)

		# set the coefficient prior
		β ~ arraydist(MvNormal.(mk, σ .* sqrt.(τ²)))

		# set the target distribution
	    for i in 1:nobs
	        mu = view(X, :, i)' * vec(β)
			y[:, i] ~ MvNormal(fill(mu,size(y)[1]), Matrix(σ*I, size(y)[1], size(y)[1]))
		end
	end
	model = grouped_lasso(yTrain,xTrain,σ_pair,λ²_pair)
end

# run a turing model on all provided population combinations
#
function getSinglePathPaired(
			df::DataFrame,combos::Array{Tuple{String,String}},
			Td::Dict,Id::Dict,Sd::Dict)
	# get the combos for pairwise obs
	T1arr = Dict.(map(get(Td,:T1p,""),combos) .=> map(get(Td,:T1sp,""),combos))
	T2arr = Dict.(map(get(Td,:T2p,""),combos) .=> map(get(Td,:T2sp,""),combos))
	dat = []
	for i in 1:length(combos)
		println("obs combo $i")
		T1=T1arr[i]
		T2=T2arr[i]
		obspair = getPairwiseObs(df,Id[:targInd],Id[:ligInd],Id[:pathInd],
										Id[:ligInd],Id[:pathInd],T1,T2)
		println("inference combo $i")
		σ_pair = 1.0
		λ²_pair = 2.0
		model = singlePathPairedModel(obspair.targT1,obspair.ligT1,obspair.ligT2,obspair.pathT1,σ_pair,λ²_pair)
		chain = sample(model, NUTS(0.65), Sd[:iter]);
		push!(dat,(T1=T1,T2=T2,chain=chain))
	end
	dat
end

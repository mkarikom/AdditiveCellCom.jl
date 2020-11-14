# load single cell data into a data frame
function loadData()
end

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

# get combinations of cells for pairwise population signaling
# df is an expression matrix output by getCellTrees
# T1,2::Dict is df.type=>df.subtype
function getPairwiseObs(
			df_orig::DataFrame,
			T1targIdx::Vector,T1ligIdx::Vector,T1pathIdx::Vector,
			T2ligIdx::Vector,T2pathIdx::Vector,
			T1::Dict,T2::Dict)
	df = df_orig[findall(x->x==0,nonunique(df_orig)),:]
	if !ismissing(first(values(T1)))
		t1targ = filter(
					row->row.type==first(keys(T1)) &&
		    		row.subtype==first(values(T1)) &&
				    any(in.(row.vertex,T1targIdx)),df)
					select!(t1targ, Not(:vertex))
					t1targ = t1targ[findall(x->x==0,nonunique(t1targ)),:]
		t1lig = filter(
					row->row.type==first(keys(T1)) &&
		    		row.subtype==first(values(T1)) &&
				    any(in.(row.vertex,T1ligIdx)),df)
					select!(t1lig, Not(:vertex))
					t1lig = t1lig[findall(x->x==0,nonunique(t1lig)),:]
		t1path = filter(
					row->row.type==first(keys(T1)) &&
		    		row.subtype==first(values(T1)) &&
				    any(in.(row.vertex,T1pathIdx)),df)
					select!(t1path, Not(:vertex))
					t1path = t1path[findall(x->x==0,nonunique(t1path)),:]
	else
		t1targ = filter(
					row->row.type==first(keys(T1)) &&
				    any(in.(row.vertex,T1targIdx)),df)
					select!(t1targ, Not(:vertex))
					t1targ = t1targ[findall(x->x==0,nonunique(t1targ)),:]
		t1lig = filter(
					row->row.type==first(keys(T1)) &&
				    any(in.(row.vertex,T1ligIdx)),df)
					select!(t1lig, Not(:vertex))
					t1lig = t1lig[findall(x->x==0,nonunique(t1lig)),:]
		t1path = filter(
					row->row.type==first(keys(T1)) &&
				    any(in.(row.vertex,T1pathIdx)),df)
					select!(t1path, Not(:vertex))
					t1path = t1path[findall(x->x==0,nonunique(t1path)),:]
	end
	if !ismissing(first(values(T2)))
		t2lig = filter(
					row->row.type==first(keys(T2)) &&
		    		row.subtype==first(values(T2)) &&
				    any(in.(row.vertex,T2ligIdx)),df)
					select!(t2lig, Not(:vertex))
					t2lig = t2lig[findall(x->x==0,nonunique(t2lig)),:]
		t2path = filter(
					row->row.type==first(keys(T2)) &&
		    		row.subtype==first(values(T2)) &&
				    any(in.(row.vertex,T2pathIdx)),df)
					select!(t2path, Not(:vertex))
					t2path = t2path[findall(x->x==0,nonunique(t2path)),:]
	else
		t2lig = filter(
					row->row.type==first(keys(T2)) &&
				    any(in.(row.vertex,T2ligIdx)),df)
					select!(t2lig, Not(:vertex))
					t2lig = t2lig[findall(x->x==0,nonunique(t2lig)),:]
		t2path = filter(
					row->row.type==first(keys(T2)) &&
				    any(in.(row.vertex,T2pathIdx)),df)
					select!(t2path, Not(:vertex))
					t2path = t2path[findall(x->x==0,nonunique(t2path)),:]
	end
	t1bc = unique(t1lig.barcode)
	t2bc = unique(t2lig.barcode)
	n = length(t1bc)*length(t2bc) # the total number of observations

	targT1 = Array{Float64,2}(undef,n,length(unique(t1targ.hgnc)))
	ligT1 = Array{Float64,2}(undef,n,length(unique(t1lig.hgnc)))
	ligT2 = Array{Float64,2}(undef,n,length(unique(t2lig.hgnc)))
	pathT1 = Array{Float64,2}(undef,n,length(unique(t1path.hgnc)))
	pathT2 = Array{Float64,2}(undef,n,length(unique(t2path.hgnc)))

	idx = 1
	for bc1 in t1bc
		for bc2 in t2bc
			targT1[idx,:] .= filter(row->row.barcode==bc1,t1targ).exp
			ligT1[idx,:] .= filter(row->row.barcode==bc1,t1lig).exp
			pathT1[idx,:] .= filter(row->row.barcode==bc1,t1path).exp
			ligT2[idx,:] .= filter(row->row.barcode==bc2,t2lig).exp
			pathT2[idx,:] .= filter(row->row.barcode==bc2,t2path).exp
			idx += 1
		end
	end
	(targT1=targT1,ligT1=ligT1,ligT2=ligT2,pathT1=pathT1,pathT2=pathT2)
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


# # run a turing model on all provided population combinations
# #
# function getSinglePathPairedThread(
# 			df::DataFrame,combos::Array{Tuple{String,String}},
# 			Td::Dict,Id::Dict,Sd::Dict)
# 	# get the combos for pairwise obs
# 	T1arr = Dict.(map(get(Td,:T1p,""),combos) .=> map(get(Td,:T1sp,""),combos))
# 	T2arr = Dict.(map(get(Td,:T2p,""),combos) .=> map(get(Td,:T2sp,""),combos))
# 	dat = Vector{Vector{NamedTuple}}()
#
# 	for i in 1:Threads.nthreads()
# 	  push!(dat,Vector{NamedTuple}())
# 	end
# 	println("running threads: ",Threads.nthreads())
# 	Threads.@threads for i in 1:length(combos)
# 		println("obs combo $i")
# 		T1=T1arr[i]
# 		T2=T2arr[i]
# 		obspair = getPairwiseObs(df,Id[:targInd],Id[:ligInd],Id[:pathInd],
# 										Id[:ligInd],Id[:pathInd],T1,T2)
# 		println("inference combo $i")
# 		σ_pair = 1.0
# 		λ²_pair = 2.0
# 		model = GCom.singlePathPairedModel(obspair.targT1,obspair.ligT1,obspair.ligT2,obspair.pathT1,σ_pair,λ²_pair)
# 		chain = sample(model, NUTS(0.65), Sd[:iter]);
# 		push!(dat[Threads.threadid()],(T1=T1,T2=T2,chain=chain))
# 	end
# 	dat = reduce(vcat,reduce(vcat,dat))
# end

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

# for a single pathway
# for a single cell
# for multiple target trees
# get tuple: lig, pathgene, target
function getCellTrees(targTree::AbstractMetaGraph,dist::Int,bcs::DataFrame)
	df = DataFrame();
	for bc in eachrow(bcs)
		bcode = bc[:barcode]
		bctype = bc[:labels]
		println("getting bc $bcode");
		ag = getPathSignals(targTree,dist,bc);
		insertcols!(ag.exp,1,(:barcode=>fill(bcode,size(ag.exp,1))))
		insertcols!(ag.exp,1,(:label=>fill(bctype,size(ag.exp,1))))
		append!(df,ag.exp);
	end
	df
end

# for a single pathway
# for a single cell
# for multiple target trees
# get tuple: lig, pathgene, target
function getCellTreesThread(targTree::AbstractMetaGraph,dist::Int,bcs::DataFrame)
	dat = Vector{Vector{DataFrame}}()
	for i in 1:Threads.nthreads()
	  push!(dat,Vector{DataFrame}())
	end
	Threads.@threads for bc in eachrow(bcs)
		bcode = bc[:barcode]
		bctype = bc[:celltype]
		bcnctype = bc[:labels]
		println("getting bc $bcode");
		ag = getPathSignals(targTree,dist,bc);
		insertcols!(ag.exp,1,(:barcode=>fill(bcode,size(ag.exp,1))))
		insertcols!(ag.exp,1,(:subtype=>fill(bcnctype,size(ag.exp,1))))
		insertcols!(ag.exp,1,(:type=>fill(bctype,size(ag.exp,1))))
		push!(dat[Threads.threadid()],ag.exp)
	end
	dat = reduce(vcat,reduce(vcat,dat))
end


# for several pathways
# for a collection of cells do:
# getPathSignal(cell,pathway)
function getMultiPathModel()
	@model function grouped_lasso(ligT1, ligT2, pathT1, pathT2, targT1, targT2, σ², λ²)
		# set the scale prior

		# grouping priors: one prior per group per target
		τ²[1] ~ Gamma((2+1)/2,λ²/2) # prior on coefs for cell type 1 autonomous
		τ²[2] ~ Gamma((2+1)/2,λ²/2) # prior on coefs for cell type 1 non-autonomous
		τ²[3] ~ Gamma((2+1)/2,λ²/2) # prior on coefs for cell type 2 autonomous
		τ²[4] ~ Gamma((2+1)/2,λ²/2) # prior on coefs for cell type 2 non-autonomous

		npath = size(pathT1,2) # number of pathway genes
		nlig = size(ligT1,2)# number of ligand genes

		β[1] ~ MvNormal(npath+nlig,σ²*τ²[1]) # coef cell type 1 autonomous
		β[2] ~ MvNormal(npath+nlig,σ²*τ²[2]) # coef cell type 1 non autonomous
		β[3] ~ MvNormal(npath+nlig,σ²*τ²[3]) # coef cell type 2 autonomous
		β[4] ~ MvNormal(npath+nlig,σ²*τ²[4]) # coef cell type 2 non autonomous

		targT1 ~ MvNormal(β[1]*[ligT1 pathT1] + β[2]*[ligT2 pathT1]) # piecewise autonomous (T1) non-autonomous (T2)
		targT2 ~ MvNormal(β[1]*[ligT2 pathT2] + β[2]*[ligT1 pathT2]) # piecewise autonomous (T2) non-autonomous (T1)
	end
end

# compose and fit the linear model
function fitModel()
end

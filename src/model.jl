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
		println("getting bc $bcode");
		ag = getPathSignals(targTree,dist,bc);
		append!(df,ag.exp);
	end
	df
end




# for several pathways
# for a collection of cells do:
# getPathSignal(cell,pathway)
function getMultiPathModel()
end

# compose and fit the linear model
function fitModel()
end

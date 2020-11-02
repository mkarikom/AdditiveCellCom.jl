
function annotateGraphLRT!(g,dbParams)
	## annotate graph vertices where a given GO molecular function was defined in the literature
	# compose the fcnParams
	goterms = Dict("GO_0038023"=>"receptor","GO_0048018"=>"ligand")
	fcnParams = Dict(
		:goFilter=>delimitValues(collect(keys(goterms)),"http://nextprot.org/rdf/terminology/","<>"),
		:dbFilter=>x->x.entId[findall(db->db=="uniprot knowledgebase",skipmissing(x.entIdDb))],
		:annotationMapAncTerm=>x->map(t->get(goterms,split(t,"/")[end],""),skipmissing(x.goAncestor)),
		:vertKey=>:entId,
		:annotationKey=>:roleLR)
	[fcnParams[k] = v for (k,v) in dbParams]
	annLR = annotateGraphFcn!(fcnParams,g)

	## annotate graph vertices with the native (eg human for pathway commons) gene ids
	annGene = annotateGraphGene!(dbParams,g)
end


# # find the simple lig rec
# lrverts = filterVertices(g,:roleLR,v->true,p->get(p,:roleLR,missing)) # the vertices representing complexes
# # find the complexes which include a ligand, receptor or both, excluding simple lig rec

# traverse inward edges from a start node and get values of target features
# return a dict with a vector of verts for each key
function getLRtree(g::AbstractMetaGraph,dir::Symbol,
					start::Int,targetFeatures::Dict)
	# get the subgraph
	vertFound = Dict()
	ng = bfs_tree(g, start, dir=dir)
	connected_v = findall(x->x > 0,degree(ng))
	pc_sg,vmap = induced_subgraph(g,connected_v)
	for k in keys(targetFeatures)
		k_fnd = []
		for v in targetFeatures[k]
			res = filterVertices(pc_sg,k,p->p==v,p->get(p,k,missing)) # the vertices representing complexes
			vertFound[v]=sort(res.ind)
		end
	end
	(v=vertFound,g=pc_sg)
end

# identify biochemical reactions where the input is Dna and output is protein
function getTransTargs(g::AbstractMetaGraph)
    ctrl = ["http://www.biopax.org/release/biopax-level3.owl#Catalysis",
 			"http://www.biopax.org/release/biopax-level3.owl#Control"]

    colnames = [:ctrlInd,:geneInd,:protInd,
                :ctrlRef,:geneRef,:protRef]
    coltypes = [Union{Missing,Int64},Union{Missing,Int64},Union{Missing,Int64},
                Union{Missing,String},Union{Missing,String},Union{Missing,String}]
    edgeTable = DataFrame(coltypes,colnames)

	for v in vertices(g)
		# identify biochemical rxn
		if haskey(props(g,v),:intType)
			if props(g,v)[:intType] == "http://www.biopax.org/release/biopax-level3.owl#BiochemicalReaction"
				ns_o = outneighbors(g,v)
		        ns_i = inneighbors(g,v)
				println("out neigh = ",length(ns_o),", in neigh = ",length(ns_i))
				for n_o in ns_o
					for n_i in ns_i
						if haskey(props(g,n_i),:participantType)
							if props(g,n_i)[:participantType] == "http://www.biopax.org/release/biopax-level3.owl#Dna"
								if haskey(props(g,n_o),:participantType)
									if props(g,n_o)[:participantType] == "http://www.biopax.org/release/biopax-level3.owl#Protein"
										rxref = props(g,v)[:interaction]
										inref = props(g,n_i)[:participant]
										outref = props(g,n_o)[:participant]
										data = tuple(v,n_i,n_o,rxref,inref,outref)
										push!(edgeTable,initRow(colnames,colnames,data))
									end
								else
									println("no output participant $n_o")
								end
							end
						else
							println("no input participant $n_i")
						end
					end
				end
			end
		end
	end
	edgeTable
end

function getTxGraphs(g::AbstractMetaGraph,dir::Symbol,txList::DataFrame,targetFeatures::Dict)
	# pull the protein ind (should already have ENSG id for zebrafish)
	paths = []
	for i in 1:size(txList,1)
		protInd = txList[i,:protInd]
		lrt = getLRtree(g,:in,protInd,targetFeatures)

		# find the product protInd in the new lrt graph
		inds = []
		for v in 1:nv(lrt.g)
			prp = props(lrt.g,v)
			if prp == props(g,protInd)
				push!(inds,v)
			end
		end
		orig = reduce(vcat,inds)
		p = dijkstra_shortest_paths(reverse(lrt.g), orig);
		pathverts = []
		allverts = Vector{eltype(collect(vertices(g)))}()
		dstverts = []
		for ft in targetFeatures[:roleLR]
			println("processing target feature $ft")
			dst = filterVertices(lrt.g,:roleLR,f->f==ft)
			ep = enumerate_paths(p,dst)
			push!(pathverts,(ft,unique(reduce(vcat,ep))))
			push!(allverts,reduce(vcat,ep)...)
			push!(dstverts,(ft,unique(dst)))
		end
		allverts = unique(allverts)

		# construct the lig rec targ graph
		lrtg,vmap = induced_subgraph(lrt.g,allverts)
		recind = filterVertices(lrtg,:roleLR,f->f=="receptor",p->get(p,:ensId,missing))
		ligind = filterVertices(lrtg,:roleLR,f->f=="ligand",p->get(p,:ensId,missing))
		targgene = props(g,txList[i,:geneInd])[:entId]
		targind = filterVertices(lrtg,:ensId,f->f==targgene)
		recgenes = unique(recind.val)
		liggenes = unique(ligind.val)
		push!(paths,(recind=recind,
					 ligind=ligind,
					 targind=targind,
					 recgene=recgenes,
					 liggene=liggenes,
					 targgene=targgene,
					 graph=lrtg))
	end
	paths
end

# get a collection of subgraphs, one for each target tree
function getLigRecTargs(g)
	# get the transcription targets (via graph traversal on g)
	tt = getTransTargs(g)

	# get the graphs and special verts corresponding to lig rec
	targetFeatures = Dict(:roleLR=>["ligand","receptor"])
	paths = getTxGraphs(g,:in,tt,targetFeatures)
end

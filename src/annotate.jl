
function annotateGraphLRT!(g,dbParams)
	## annotate graph vertices where a given GO molecular function was defined in the literature
	# compose the fcnParams
	# goterms = Dict("GO_0038023"=>"receptor","GO_0005102"=>"ligand")
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
	annP = annotateGraphP!(dbParams,g)

	annG = annotateGraphG!(dbParams,g)
	(lr=annLR,p=annP,g=annG)
end


# traverse inward edges from a start node and get values of target features
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
	(v=vertFound,g=pc_sg,vmap=vmap)
end

# like getLRtree except that it is protein target agnostic to account for degenerate control reacions like https://reactome.org/content/detail/R-HSA-1980065
# traverse inward edges from the controller of a target gene, and append the graph to the gene as an out neighbor
# filter gene inneighbors so that only the target gene is included
function getLRgeneTree(g::AbstractMetaGraph,dir::Symbol,
					geneInd::Int,ctrlInd::Int,targetFeatures::Dict)
	# get the subgraph
	ng = bfs_tree(g, ctrlInd, dir=dir)
	connected_v = findall(x->x > 0,degree(ng))
	pc_sg,vmap = induced_subgraph(g,connected_v)

	# find all the extraneous inneighbor genes of the ctrl vertex (these are caused by degenerate ctrl rxn like https://reactome.org/content/detail/R-HSA-1980065)
	pc_sg_ctrl_in = inneighbors(pc_sg,findfirst(v->v==ctrlInd,vmap))
	pc_sg_genes = filterVertices(pc_sg,:entIdDb,n->n=="ensembl")
	pc_sg_gene_ind = findfirst(v->v==geneInd,vmap)
	nontarg = setdiff(pc_sg_genes,pc_sg_gene_ind)
	vremove = intersect(nontarg,pc_sg_ctrl_in)

	# remove extraneous genes and update the vmap for the induced subgraph
	for v in vremove
		rem_vertex!(pc_sg, v) # remove the vertex
		vmap[end] = vmap[v]
		vmap = vmap[1:end-1]
	end

	# get all target features
	vertFound = Dict()
	for k in keys(targetFeatures)
		k_fnd = []
		for v in targetFeatures[k]
			res = filterVertices(pc_sg,k,p->p==v,p->get(p,k,missing)) # the vertices representing complexes
			vertFound[v]=sort(res.ind)
		end
	end
	(v=vertFound,g=pc_sg,vmap=vmap)
end

# identify unique triples of:
#    1. ctrlInd (the vertex ind of the control physical entity that activates or suppresses transcription at geneInd))
#    2. geneInd (the vertex ind of gene whose transcription is controlled by "ctrlInd", in general ctrlInd->geneInd is one-to-many)
#    3. protInd (the vertex ind of protein product, in general geneInd->protInd is one-to-many)
# CAUTION: ctrl rxns like https://reactome.org/content/detail/R-HSA-1980065 will cause the geneInd->protInd to break
# In order to circumvent this only look at unique(edgeTable,[:ctrlInd,:geneInd]) for downstream analysis
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

# get a vector of single-target trees, one for each unique target gene in txList
function getTxGraphs(g::AbstractMetaGraph,dir::Symbol,txList::DataFrame,targetFeatures::Dict)
	paths = []
	for i in 1:size(txList,1)
		# get the subgraph of all reachable vertices from a transcriptional target
		geneInd = txList[i,:geneInd]
		ctrlInd = txList[i,:ctrlInd]
		lrt = getLRgeneTree(g,:in,geneInd,ctrlInd,targetFeatures)
		sgCtrlInd = findfirst(v->v==ctrlInd,lrt.vmap) # the vert ind of the ctrl entity in the new subgraph
		sgGeneInd = findfirst(v->v==geneInd,lrt.vmap) # the vert ind of the target gene in the new subgraph

		# trace all paths from the product to the target features, eg ligands and receptors
		p = dijkstra_shortest_paths(reverse(lrt.g), sgCtrlInd);
		pathverts = []
		allverts = Vector{eltype(collect(vertices(g)))}()
		dstverts = []
		for ft in targetFeatures[:roleLR]
			println("processing target feature $ft")
			dst = filterVertices(lrt.g,:roleLR,f->f==ft)
			ep = enumerate_paths(p,dst)
			push!(pathverts,(ft,unique(reduce(vcat,ep)))) # eg for targetFeatures[:ligand] push ("ligand", [46, 48, 49, 50, 51, 66, 67, 68, 45, 44  …  36, 26, 38, 37, 84, 88, 118, 122, 123, 127])
			push!(allverts,reduce(vcat,ep)...) # eg for ex above push all verts in path [46, 48, 49, 50, 51, 66, 67, 68, 45, 44  …  36, 26, 38, 37, 84, 88, 118, 122, 123, 127]...
			push!(dstverts,(ft,unique(dst))) # eg for ex above push just the feature verts
		end
		push!(allverts,sgGeneInd) # add the target gene to the vertex list
		allverts = unique(allverts)

		# construct the lig rec targ graph
		lrtg,vmap = induced_subgraph(lrt.g,allverts)
		recind = filterVertices(lrtg,:roleLR,f->f=="receptor",p->get(p,:ensId,missing))
		ligind = filterVertices(lrtg,:roleLR,f->f=="ligand",p->get(p,:ensId,missing))
		targgene = props(g,txList[i,:geneInd])[:entId]
		targind = findfirst(v->v==sgGeneInd,vmap)
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

# get a multi-target tree, representing the union of all single-target trees for each unique target protein in txList
# where target protein is a unique protein product of a transcriptional target gene
function getTxGraph(g::AbstractMetaGraph,dir::Symbol,txList::DataFrame,targetFeatures::Dict)
	# get the subgraph of all reachable vertices from a set of transcriptional targets

	pathverts = []
	allverts = Vector{eltype(collect(vertices(g)))}()
	dstverts = []
	for i in 1:size(txList,1)
		# get the subgraph of all reachable vertices from a transcriptional target
		geneInd = txList[i,:geneInd]
		ctrlInd = txList[i,:ctrlInd]
		lrt = getLRgeneTree(g,:in,geneInd,ctrlInd,targetFeatures)
		# make sure that receptors and ligands are reachable
		if all(length.([lrt.v[v] for v in first(values(targetFeatures))]) .> 0)
			sgCtrlInd = findfirst(v->v==ctrlInd,lrt.vmap) # the vert ind of the ctrl entity in the new subgraph
			sgGeneInd = findfirst(v->v==geneInd,lrt.vmap) # the vert ind of the target gene in the new subgraph

			# trace all paths from the product to the target features, eg ligands and receptors
			p = dijkstra_shortest_paths(reverse(lrt.g), sgCtrlInd);
			for ft in targetFeatures[:roleLR]
				println("processing target feature $ft")
				dst = filterVertices(lrt.g,:roleLR,f->f==ft)
				ep = enumerate_paths(p,dst)
				push!(pathverts,(targ=props(lrt.g,sgGeneInd)[:entId],feat=ft,verts=lrt.vmap[unique(reduce(vcat,ep))])) # eg for targetFeatures[:ligand] push ("ligand", [46, 48, 49, 50, 51, 66, 67, 68, 45, 44  …  36, 26, 38, 37, 84, 88, 118, 122, 123, 127])
				push!(allverts,lrt.vmap[unique(reduce(vcat,ep))]...) # eg for ex above push all verts in path [46, 48, 49, 50, 51, 66, 67, 68, 45, 44  …  36, 26, 38, 37, 84, 88, 118, 122, 123, 127]...
				push!(dstverts,(targ=props(lrt.g,sgGeneInd)[:entId],feat=ft,verts=lrt.vmap[unique(dst)])) # eg for targetFeatures[:ligand] push ("ligand", [46, 48, 49, 50, 51, 66, 67, 68, 45, 44  …  36, 26, 38, 37, 84, 88, 118, 122, 123, 127])
			end
			push!(allverts,geneInd) # add the target gene to the vertex list
		else
			tgene = map(x->props(g,x)[:displayName],txList.geneInd[i])
			for v in first(values(targetFeatures))
				lv = length(lrt.v[v])
				if lv <= 0
					println("$tgene unreachable from any $v")
				end
			end
		end
	end
	allverts = unique(allverts)

	# construct the lig rec targ graph
	lrtg,vmap = induced_subgraph(g,allverts)
	recinds = filterVertices(lrtg,:roleLR,f->f=="receptor",p->get(p,:ensId,missing))
	liginds = filterVertices(lrtg,:roleLR,f->f=="ligand",p->get(p,:ensId,missing))
	targinds = [findfirst(x->x==vert,vmap) for vert in txList.geneInd]
	targgenes = map(v->props(lrtg,v)[:entId],targinds)
	recgenes = unique(recinds.val)
	liggenes = unique(liginds.val)
	return (recinds=recinds,liginds=liginds,targinds=targinds,
			recgenes=recgenes,liggenes=liggenes,targgenes=targgenes,
			graph=lrtg)
end

# get a collection of subgraphs, one for each target tree
function getLigRecTargs(g)
	# get the transcription targets (via graph traversal on g)
	tt = getTransTargs(g)
	tt = unique(tt,[:ctrlInd,:geneInd])
	# get the graphs and special verts corresponding to lig rec
	targetFeatures = Dict(:roleLR=>["ligand","receptor"])
	paths = getTxGraphs(g,:in,tt,targetFeatures)
end

# get a subgraph of the union of single-target trees
function getLigRecTarg(g)
	# get the transcription targets (via graph traversal on g)
	tt = getTransTargs(g)
	tt = unique(tt,[:ctrlInd,:geneInd])
	# get the graphs and special verts corresponding to lig rec
	targetFeatures = Dict(:roleLR=>["ligand","receptor"])
	paths = getTxGraph(g,:in,tt,targetFeatures)
end

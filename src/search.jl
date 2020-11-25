## infer ligand/receptor pairs by search
function searchLR(g)
    # output
    res = []
    # find the lig/rec vertices
    lr = filterVertices(g,:roleLR,v->true)
    for v in lr
        if indegree(g,v) == 0
            # define the original lig or rec
            sParam = Dict(:role0=>props(g,v)[:roleLR], # the role of the vertex LR v0
                          :role1=>first(filter(r->r!=props(g,v)[:roleLR],("ligand","receptor"))), # the complement role
                          :v0=>v, # the index of v0
                          :v1=>[], # the indices of complement vertices
                          :vcurrent=>nothing, # State: current vertex position
                          :vexclude=>[v], # State: previously searched parents
                          :dec0=>nothing) # index of the closest descendent of v0 and v1
            println("searching $v")

            for otn in outneighbors(g,v)
                # start recursive search for complement in parent search mode
                # when stop:
                # 1) sParam[:v1] is populated with any complement vertices
                # 2) sParam[:dec0] is populated with nearest descendent vertex
                sParam[:vcurrent]=otn
                push!(sParam[:vexclude],v)
                sParam[:dec0]=v
                searchParents!(g,sParam)
            end
            if length(sParam[:v1]) > 0
                push!(res,sParam)
            end
        end
    end
    res
end

## recursive loop search for complement
# search for complement among parents
function searchParents!(g,sParam)
    push!(sParam[:vexclude],sParam[:vcurrent])
    # stop if no parents
    println("searching ",sParam[:vcurrent])
    if length(inneighbors(g,sParam[:vcurrent])) == 0
        pv = props(g,sParam[:vcurrent])
        if haskey(pv,:roleLR) && pv[:roleLR]==sParam[:role1]# filter L,
                push!(sParam[:v1],sParam[:vcurrent]) # add complement
                return
        end
    else # there are parents
        for v in inneighbors(g,sParam[:vcurrent]) # search all parents excluding previously searched ancestors
            if !any(in.(v,sParam[:vexclude]))
                sParam[:dec0]=sParam[:vcurrent]
                sParam[:vcurrent]=v
                searchParents!(g,sParam)
            end
        end
    end
end

# search for complement among children (run after parent mode comes up dry)
function searchChildren!(g,sParam)
    # Not Implmented
end


## infer LRT triples by searching from LR sParam[:dec0] to a reachable target
# only the target control vertex is generally accessible to transduction so add the index of the protein afterward
function searchT(g,lrset,targs)
    paths = []
    for t in 1:size(targs)[1]
        for l in lrset
            dsp = dijkstra_shortest_paths(g,l[:dec0])
            enum = enumerate_paths(dsp,targs.ctrlInd[t])
            if length(enum) > 0
                push!(paths,Dict(
                        :dec0=>l[:dec0],
                        :targCtrl=>targs.ctrlInd[t],
                        :targProt=>targs.protInd[t],
                        :paths=>enum,
                        :v0=>l[:v0],
                        :v1=>l[:v1]))
            end
        end
    end
    paths
end

## find all possible triples based on target paths from sParam[:d0]
function enumerateLRT(g,paths,drange)
    up0 = []
    ens0 = []
    hgnc0 = []
    upOrth = []
    hgncOrth = []
    for p in paths
        v0 = p[:v0]
        for v1 in p[:v1]
            pv0 = props(g,v0)
            pv1 = props(g,v1)
            ptarg = props(g,p[:targProt])
            push!(up0,Dict(
                    Symbol(pv0[:roleLR])=>pv0[:entId],
                    Symbol(pv1[:roleLR])=>pv1[:entId],
                    :target=>ptarg[:entId]))
            push!(ens0,Dict(
                    Symbol(pv0[:roleLR])=>pv0[:ensId],
                    Symbol(pv1[:roleLR])=>pv1[:ensId],
                    :target=>ptarg[:ensId]))
            push!(hgnc0,Dict(
                    Symbol(pv0[:roleLR])=>pv0[:gname],
                    Symbol(pv1[:roleLR])=>pv1[:gname],
                    :target=>ptarg[:gname]))

            orthv0hgnc = Dict()
            orthv0upid = Dict()
            for d in drange
                odv0 = getOrthDist(g,v0,d)
                odv1 = getOrthDist(g,v1,d)
                odt = getOrthDist(g,p[:targProt],d)
                push!(orthv0hgnc,d=>Dict(Symbol(pv0[:roleLR])=>odv0.hgnc,
                                         Symbol(pv1[:roleLR])=>odv1.hgnc,
                                         :target=>odt.hgnc))
                push!(orthv0upid,d=>Dict(Symbol(pv0[:roleLR])=>odv0.upid,
                                         Symbol(pv1[:roleLR])=>odv1.upid,
                                         :target=>odt.upid))
            end
            push!(hgncOrth,orthv0hgnc)
            push!(upOrth,orthv0upid)
        end
    end
    (up0=unique(up0),ens0=unique(ens0),hgnc0=unique(hgnc0),hgncOrth=hgncOrth,upOrth=upOrth)
end

# create a triples dataframe from output of enumerateLRT
function lrtDist(arr,dist)
    df = DataFrame(ligand=Vector{Vector{String}}(),receptor=Vector{Vector{String}}(),target=Vector{Vector{String}}())
    for a in arr
        push!(df,a[dist])
    end
    unique(df)
end

# get the max ortholog expression for a native gene at the given vertex at a given distance
function quantMax(g,v,dist)
    quantmax = missing
    if haskey(props(g,v),:orthoDist)
        vorth = props(g,v)[:orthoDist][dist][:members]
        if length(vorth) > 0
            println("length of vorth is ",length(vorth))
            if length(collect(skipmissing([haskey(vorth[o],:rnaSeq) ? vorth[o][:rnaSeq][:value] : missing for o in 1:length(vorth)]))) > 0
                quantmax = reduce(max,skipmissing([haskey(vorth[o],:rnaSeq) ? vorth[o][:rnaSeq][:value] : missing for o in 1:length(vorth)]))
                println(quantmax)
            else
                # println("no expression of ortholog")
            end
        else
            # println("vorth is empty")
        end
    else
        # println("orthoDist missing")
    end
    quantmax
end

# get all orthologs at the specified distance
# ensemble gene ids are not always available from orthodb so omit them here
function getOrthDist(g,v,dist)
    hgnc = []
    upid = []
    if haskey(props(g,v),:orthoDist)
        vorth = props(g,v)[:orthoDist][dist][:members]
        if length(vorth) > 0
            println("length of vorth is ",length(vorth))
            for o in vorth
                println("saving ids")
                push!(hgnc,o[:orthoHGNC])
                push!(upid,o[:orthoEntIds]...)
            end
        else
            # println("vorth is empty")
        end
    else
        # println("orthoDist missing")
    end
    (hgnc=hgnc,upid=upid)
end

# Make a prediction given an input vector. (from Turing.jl docs)
# in most cases eg NUTS n_adapt the "burn-in" samples have already been removed
function prediction(chain, x, startInd=1)
    p = get_params(chain[startInd:end, :, :])
    xb = reduce(hcat,map(β->reduce(vcat,β),p.β))*x
    intercept = reduce(hcat,map(i->reduce(vcat,i),p.intercept))
    targets = [xb .+ c for c in eachcol(intercept)]
    return mean(reshape(reduce(vcat,targets)',(size(x)[2],size(intercept)[2],:)),dims=3)[:,:,1]'
end

# make sure the data can be scaled and add noise if not
# this is a temporary fix, eventually we will just choose features that do not have this problem
function checkScaling!(data,varname)
    allsame = []
    rows = collect(eachrow(data))
    for r in 1:length(rows)
        if all(i->i==rows[r][1],rows[r])
            push!(allsame,r)
        end
    end
    if length(allsame) > 0
        println("$varname: zero variance detected in features: $allsame, adding noise")
        # data += rand(size(data)...)*1e-6
        data .= data + rand(size(data)...)*1e-6
    end
    nothing
end

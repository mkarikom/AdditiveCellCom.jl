# create a subdirectory of parent and store separate jld object for each dict
function archiveBarcodeGraphs(parentDir::String,dirName::String,data::Vector)
    drnm = joinpath(parentDir, dirName)
    mkpath(drnm)

    for d in data
        fname = String(first(keys(d)))
        gr = first(values(d))
        println("saving bc tree $fname")
        @save joinpath(drnm,"$fname.jld2") {compress=true} gr # @load joinpath(drnm,"wt.jld2") wt
    end
end

# create a subdirectory of parent and store the barcode graphs and summary expression dataframe (output of getCellTrees)
function archiveBarcodeExpression(parentDir::String,dirName::String,data::NamedTuple)
    drnm = joinpath(parentDir, dirName)
    mkpath(drnm)
    df = data.df
    @save joinpath(drnm,"expression.jld2") {compress=true} df  # @load joinpath(drnm,"wt.jld2") wt

    archiveBarcodeGraphs(drnm,"bc_trees",data.g)

end

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

function summarizeTrajectories(g,pathlist,human=true)
	dkeys = ["ligand","receptor","expression_target","expression_control_rxn","trajectory"]
	colnames = Symbol.(dkeys)
	coltypes = fill(Union{Missing,Any},length(colnames))
	df = DataFrame(colnames .=> [type[] for type in coltypes])
	for p in pathlist
		row_data = []
		# identify v0
		if props(g,p[:v0])[:roleLR]=="ligand"
			ligand = human ? props(g,p[:v0])[:displayName] : p[:v0]
		elseif props(g,p[:v0])[:roleLR]=="receptor"
			receptor = human ? props(g,p[:v0])[:displayName] : p[:v0]
		else
			print("v0 is: $(props(g,p[:v0]))")
		end


		# # identify v1
		for v1 in p[:v1]
			if props(g,v1)[:roleLR]=="ligand"
				ligand = human ? props(g,v1)[:displayName] : v1
			elseif props(g,v1)[:roleLR]=="receptor"
				receptor = human ? props(g,v1)[:displayName] : v1
			else
				print("v1 is: $(props(g,v1)[:roleLR])")
			end
			expression_target = human ? props(g,p[:targProt])[:displayName] : p[:targProt]
			expression_ctrl_rxn = human ? props(g,p[:targCtrl])[:interaction] : p[:targCtrl]
			trajectory = p[:paths]
			row_data = [ligand,receptor,expression_target,expression_ctrl_rxn,trajectory]
			push!(df,initRow(colnames,colnames,tuple(row_data...)))
		end 
	end
	unique(df)
end

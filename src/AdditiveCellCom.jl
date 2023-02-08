module AdditiveCellCom

using Turing
using DataFrames,Query
using Graphs,MetaGraphs
using PCquery
using JLD2
using LinearAlgebra
using Distributed
using StatsBase

export searchLR,searchT,enumerateLRT,quantMax,getOrthDist,lrtDist
export getPathSignals,getCellTrees,getCellTreesThread,getMultiObs,singlePathPairedModel,getSinglePathPaired,filterTypes
export archiveBarcodeGraphs,archiveBarcodeExpression
export annotateGraphLRT!,getLRtree,getLRgeneTree,getTransTargs,getTxGraphs,getTxGraph,getLigRecTargs,getLigRecTarg
export prediction, checkScaling!

include("search.jl") # search the graph for ligand/receptor combos

include("annotate.jl") # modification of graph

include("features.jl") # acquire features from the graph

include("model.jl") # fit signaling model

include("util.jl") # storage, etc
end

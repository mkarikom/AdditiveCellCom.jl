module GCom

using Turing
using DataFrames,LightGraphs,MetaGraphs
using PCquery
using JLD2
using LinearAlgebra
using Distributed
using StatsBase

export getPathSignals,getCellTrees,getCellTreesThread,getMultiObs,singlePathPairedModel,getSinglePathPaired,filterTypes
export archiveBarcodeGraphs,archiveBarcodeExpression
export annotateGraphLRT!,getLRtree,getLRgeneTree,getTransTargs,getTxGraphs,getTxGraph,getLigRecTargs,getLigRecTarg
export prediction, checkScaling!

include("annotate.jl") # modification of graph

include("features.jl") # acquire features from the graph

include("model.jl") # fit signaling model

include("util.jl") # storage, etc
end

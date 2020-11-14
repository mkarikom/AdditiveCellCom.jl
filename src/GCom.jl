module GCom

using Turing
using DataFrames,LightGraphs,MetaGraphs
using PCquery
using JLD2
using LinearAlgebra
using Distributed

export loadData,getPathSignal,getPairwiseObs
export getSinglePathPairedDistributed,getSinglePathPaired,singlePathPairedModelDistributed,singlePathPairedModel
export annotateGraphLRT!,getLigRecTargs,getTxGraphs,getLigRecTarg,getTxGraph

include("annotate.jl") # modification of graph

include("model.jl") # fit signaling model

include("util.jl") # storage, etc
end

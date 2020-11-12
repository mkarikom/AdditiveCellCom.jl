module GCom

using Turing
using DataFrames,LightGraphs,MetaGraphs
using PCquery

export loadData,getPathSignal,getMultiPathModel,fitModel

export annotateGraphLRT!,getLigRecTargs,getTxGraphs,getLigRecTarg,getTxGraph

include("annotate.jl") # modification of graph

include("model.jl") # fit signaling model
end

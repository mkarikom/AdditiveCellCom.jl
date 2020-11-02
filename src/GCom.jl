module GCom

using Turing
using DataFrames,LightGraphs,MetaGraphs
using PCquery,SeuratRDS

export loadData,getPathSignal,getMultiPathModel,fitModel

export annotateGraphLRT!,getLigRecTargs,getTxGraphs

include("annotate.jl") # modification of graph

include("model.jl") # fit signaling model
end

module ITPM_SimX

using NearestNeighbors, LinearAlgebra, SparseArrays, SuiteSparse

include("libs.jl")
include("funs_2D.jl")
include("funs_sat.jl")
include("linAlgLib.jl")

export make_gdm, make_well, make_sim,
       make_gdm_prop_sat,
       make_sim2f,
       calc_sat_step,
       calc_wtc

# Write your package code here.

end

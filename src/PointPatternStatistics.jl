module PointPatternStatistics
using LinearAlgebra: norm
using StaticArrays
using RecipesBase
using Requires


struct PointPattern{PT, WT}
    data::Vector{PT}
    window::WT
end
Base.length(pp::PointPattern) = npoints(pp)
npoints(pp) = length(pp.data)
window(pp) = pp.window

Base.iterate(pp::PointPattern) = iterate(pp.data)
Base.iterate(pp::PointPattern, state) = iterate(pp.data, state)
Base.eltype(pp::PointPattern) = eltype(pp.data)

getx(x::NamedTuple) = x.x
gety(x::NamedTuple) = x.y
getx(x) = x[1]
gety(x) = x[2]

include("window.jl")
include("weighted_distance_histogram.jl")
include("pair_correlation_function.jl")
include("K.jl")
include("kaplanmeier.jl")
#include("distance_transform.jl")
#include("bdist.jl")
#include("Fdt.jl")
include("Fnn.jl")
#include("F.jl")
include("G.jl")
include("envelope_erl.jl")

include("simple_sequential_inhibition.jl")

include("plotrecipes.jl")
Fest(pp, r, N=100) = Fkmnn(pp.data, window(pp), r, N)
const Gest = Gkm
const Kest = Ktrans
const globalenvelope = erlenvelope
export PointPattern
export Kest, pcf, Fest, Gest, globalenvelope, Lest, L12
export SimpleSequentialInhibition
export inside, window

function __init__()
    @require RCall="6f49c342-dc21-5d91-9882-a32aef131414" include("convertR.jl")
end

end # module

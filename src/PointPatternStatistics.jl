module PointPatternStatistics
using LinearAlgebra: norm
using StaticArrays


struct PointPattern{PT, WT}
    data::Vector{PT}
    window::WT
end
npoints(pp) = length(pp.data)
window(pp) = pp.window

getx(x::NamedTuple) = x.x
gety(y::NamedTuple) = x.y
getx(x::Tuple) = x[1]
gety(x::Tuple) = x[2]
getx(x::SVector) = x[1]
gety(x::SVector) = x[2]

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
Fest(pp, r, N=100) = Fkmnn(pp.data, window(pp), r, N)
Gest = Gkm
Kest = Ktrans
globalenvelope = erlenvelope
export PointPattern
export Kest, pcf, Fest, Gest, globalenvelope, Lest, L12
export SimpleSequentialInhibition
export inside

end # module

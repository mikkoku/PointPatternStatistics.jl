module PointPatternStatistics
using LinearAlgebra: norm
using StaticArrays

getxy(x::Tuple{<:Number, <:Number}) = x
getxy(x::SVector{2, <:Number}) = x
getxy(x::NamedTuple) = (x.x, x.y)

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
Fest(x, w, r, N=100) = Fkmnn(x, w, r, N)
Gest = Gkm
Kest = Ktrans
globalenvelope = erlenvelope
export Kest, pcf, Fest, Gest, globalenvelope, Lest, L12
export inside

end # module

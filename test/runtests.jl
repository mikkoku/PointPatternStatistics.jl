using PointPatternStatistics
using Test
include("reference_implementation.jl")

window = (x=(-1.3,1), y=(0,1.6))
randu((a, b)) = a + rand() * (b-a)
randw(window) = (randu(window.x), randu(window.y))
xy = [randw(window) for _ in 1:300]
r = range(0, 0.25, length=512)
est = pcf(xy, window, r)
ref = pcfref(xy, window, r)

@test sum(abs2, (est - ref)[2:end]) < 0.0001

est = Kest(xy, window, r)
ref = Kref(xy, window, r)
@test est ≈ ref

include("F.jl")
include("G.jl")
#include("envelope.jl")

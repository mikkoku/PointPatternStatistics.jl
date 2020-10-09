using Distances
using LinearAlgebra
using PointPatternStatistics: iscovered

window = (x=(0.0, 100.0), y=(30.0, 200.0))
x = rand(SimpleSequentialInhibition(window, 10.0))
d = pairwise(Euclidean(), reshape(reinterpret(Float64, x.data), 2, :), dims=2)
d[diagind(d)] .= Inf
@test minimum(d) >= 10.0
@test all(x -> inside(window)(x) >= 0, x.data)


@test !iscovered(0.0, 1.0, 0.0, 1.0, [(0.5, 0.5)], sqrt(0.5) - 0.00001)
@test iscovered(0.0, 1.0, 0.0, 1.0, [(0.5, 0.5)], sqrt(0.5) + 0.00001)
@test iscovered(0.0, 1.0, 0.0, 1.0, [(1.0, 0.5), (0.0,0.5)], sqrt(0.5)+0.00001)
@test !iscovered(0.0, 1.0, 0.0, 1.0, [(1.0, 0.5), (0.0,0.5)], sqrt(0.5)-0.00001)
@test !iscovered(0.0, 1.0, 0.0, 1.0, [(1.5, 0.5), (-0.5,0.5),(0.5,1.5),(0.5,-0.5)], 1-0.00001)
@test iscovered(0.0, 1.0, 0.0, 1.0, [(1.5, 0.5), (-0.5,0.5),(0.5,1.5),(0.5,-0.5)], 1+0.00001)

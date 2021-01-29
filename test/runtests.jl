using PointPatternStatistics
using Test
include("reference_implementation.jl")
@testset "PointPatternStatistics.jl" begin
    window = (x=(-1.3,1), y=(0,1.6))
    randu((a, b)) = a + rand() * (b-a)
    randw(window) = (randu(window.x), randu(window.y))
    xy = [randw(window) for _ in 1:300]
    r = range(0, 0.25, length=512)
    @testset "iterate" begin
        @test collect(PointPattern(xy, window)) == xy
    end

    @testset "pcf" begin
        est = pcf(PointPattern(xy, window), r)
        ref = pcfref(xy, window, r)
        @test sum(abs2, (est - ref)[2:end]) < 0.0001
    end

    @testset "K" begin
        est = Kest(PointPattern(xy, window), r)
        ref = Kref(xy, window, r)
        @test est ≈ ref
    end

    @testset "F" begin
        include("F.jl")
    end
    @testset "G" begin
        include("G.jl")
    end
    #include("envelope.jl")
    @testset "SSI" begin
        include("SSI.jl")
    end
end

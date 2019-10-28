using OffsetArrays
using FFTW

"""
    pcf(xy, window, r, bandwidth=nothing)

Return pair-correlation function estimate for point pattern xy using translation edge correction.

The kernel estimator only uses the r values provided plus padding.

Stoyan & Stoyan, 1994, p 284
"""
function pcf(xy, window, r, bandwidth=nothing)
    if length(xy) == 0
        return fill(0.0/0.0, length(r))
    end
    xy = getxy.(xy)
    λ = length(xy)/area(window)
    a = if bandwidth == nothing
        #min(0.15/sqrt(λ), last(r))
        0.15/sqrt(λ)
    else
        bandwidth
    end
    buf = 4*ceil(Int, a/step(r))
    A = pi * 2 * step(r) * r * length(xy) * λ
    h = disthisttrans(xy, step(r), length(r) + buf, window, true)
    p1 = vcat(fill(0.0, buf), h)
    rk = range(-step(r)*(div(length(r),2) + buf - 1), step=step(r), length=length(r) + 2buf)
    kernel = (x -> ifelse(-a < x < a, 3/4 * (1 - (x/a)^2)/a, 0.0)).(rk)
    #println(OffsetArray(rk, -div(length(kernel), 2))[0])
    kernel = OffsetArray(kernel, -div(length(kernel), 2))
    kernel ./= sum(kernel)
    real.(ifft(fft(p1) .* fft(kernel)))[(1+buf:end-buf) .- 1] ./ A
end

using OffsetArrays
using FFTW

bwstoyan(pp) = 0.15 / sqrt(npoints(pp)/area(window(pp)))
"""
    pcf(xy, window, r, bandwidth=bwstoyan)

Return pair-correlation function estimate for point pattern xy using translation edge correction.

The kernel estimator only uses the r values provided plus padding.

Bandwidth can be specified as a function of xy and window or as a number.

Stoyan & Stoyan, 1994, p 284
"""
pcf(pp, r::AbstractRange, bandwidth=bwstoyan) = pcf(pp, r, bandwidth(pp))
pcf(pp, r::AbstractRange, bandwidth::Number) = pcf(pp, pp, r, float(bandwidth))
pcf(pp1, pp2::PointPattern, r::AbstractRange, bandwidth=bwstoyan) = pcf(pp1, pp2, r, bandwidth(pp2))
# spatstat uses the bandwidth computed from the second pattern.
function pcf(pp1, pp2::PointPattern, r::AbstractRange, bandwidth::Float64)
    window(pp1) == window(pp2) || throw(ArgumentError("Both point patterns must be on the same window"))
    win = pp1.window
    if length(pp1) == 0 || length(pp2) == 0
        return fill(0.0/0.0, length(r))
    end
    #xy = getxy.(xy)
    λ = length(pp2)/area(win)
    a = bandwidth

    buf = 4*ceil(Int, a/step(r))
    A = pi * 2 * step(r) * r * length(pp1) * λ
    h = disthisttrans12(pp1.data, pp2.data, step(r), length(r) + buf, win, true)
    p1 = vcat(fill(0.0, buf), h)
    rk = range(-step(r)*(div(length(r),2) + buf - 1), step=step(r), length=length(r) + 2buf)
    kernel = (x -> ifelse(-a < x < a, 3/4 * (1 - (x/a)^2)/a, 0.0)).(rk)
    #println(OffsetArray(rk, -div(length(kernel), 2))[0])
    kernel = OffsetArray(kernel, -div(length(kernel), 2))
    kernel ./= sum(kernel)
    irfft(rfft(p1) .* rfft(kernel), length(kernel))[(1+buf:end-buf) .- 1] ./ A
end

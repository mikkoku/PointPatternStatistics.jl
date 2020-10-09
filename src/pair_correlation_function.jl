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
pcf(pp, r, bandwidth=bwstoyan) = pcf(pp, r, bandwidth(pp))
function pcf(pp, r, bandwidth::Number)
    xy = pp.data
    window = pp.window
    if length(xy) == 0
        return fill(0.0/0.0, length(r))
    end
    #xy = getxy.(xy)
    λ = length(xy)/area(window)
    a = bandwidth

    buf = 4*ceil(Int, a/step(r))
    A = pi * 2 * step(r) * r * length(xy) * λ
    h = disthisttrans(xy, step(r), length(r) + buf, window, true)
    p1 = vcat(fill(0.0, buf), h)
    rk = range(-step(r)*(div(length(r),2) + buf - 1), step=step(r), length=length(r) + 2buf)
    kernel = (x -> ifelse(-a < x < a, 3/4 * (1 - (x/a)^2)/a, 0.0)).(rk)
    #println(OffsetArray(rk, -div(length(kernel), 2))[0])
    kernel = OffsetArray(kernel, -div(length(kernel), 2))
    kernel ./= sum(kernel)
    irfft(rfft(p1) .* rfft(kernel), length(kernel))[(1+buf:end-buf) .- 1] ./ A
end

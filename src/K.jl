# Stoyan & Stoyan, 1994, 15.11

"""
    Ktrans(xy, r, window)

Return K-function estimate with translation edge correction.
"""
function Ktrans(pp, r)
  n = npoints(pp)
  cumsum(disthisttrans(pp.data, step(r), length(r), window(pp))) / (n*(n-1)/area(window(pp)))
end

"""
   Kmark(xy, r, window, f)

Return mark-weighted K-function with translation edge correction.

This is not properly normalised.
"""
function Kmark(pp, r, f)
  n = npoints(pp)
  cumsum(disthisttransf(pp.data, f, step(r), length(r), window(pp))) / (n*(n-1)/area(window(pp)))
end

function K12(x, y, r)
  window(x) == window(y) || throw(ArgumentError("Both point patterns must be on the same window"))
  n1 = npoints(x)
  n2 = npoints(y)
  cumsum(disthisttrans12(x.data, y.data, step(r), length(r), window(x))) / (n1*n2/area(window(x)))
end
function L12(x, y, r)
  sqrt.(K12(x, y, r) ./ pi)
end
Lest(pp, r) = sqrt.(Kest(pp, r) ./ pi)

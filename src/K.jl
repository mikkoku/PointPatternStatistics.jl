# Stoyan & Stoyan, 1994, 15.11

"""
    Ktrans(xy, r, window)

Return K-function estimate with translation edge correction.
"""
function Ktrans(xy, window, r)
  xy = getxy.(xy)
  n = length(xy)
  cumsum(disthisttrans(xy, step(r), length(r), window)) / (n*(n-1)/area(window))
end

"""
   Kmark(xy, r, window, f)

Return mark-weighted K-function with translation edge correction.

This is not properly normalised.
"""
function Kmark(xy, window, r, f)
  xy = getxy.(xy)
  n = length(xy)
  cumsum(disthisttransf(xy, f, step(r), length(r), window)) / (n*(n-1)/area(window))
end

function K12(x, y, window, r)
  n1 = length(x)
  n2 = length(y)
  cumsum(disthisttrans12(getxy.(x), getxy.(y), step(r), length(r), window)) / (n1*n2/area(window))
end
function L12(x, y, window, r)
  sqrt.(K12(x, y, window, r) ./ pi)
end
Lest(xy, window, r) = sqrt.(Kest(xy, window, r) ./ pi)

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

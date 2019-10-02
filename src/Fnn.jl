using NearestNeighbors
using StaticArrays
"""
    Fkmnn(xy, r, window, N)

Return empty space function based on Kaplan-Meier estimator.
"""
function Fkmnn(xy, window, r, N)
    xy = getxy.(xy)
    (x1,x2), (y1, y2) = window.x, window.y
    tree = KDTree(SVector.(xy))
    sx = (x2-x1)/N
    Ny = round(Int, N*(y2-y1)/(x2-x1))
    sy = (y2-y1)/Ny
    rx = range(x1 + 0.5sx, x2 - 0.5sx, length=N)
    ry = range(y1 + 0.5sy, y2 - 0.5sy, length=Ny)
    grid = vec([SVector(x, y) for x in rx, y in ry])
    _, D = knn(tree, grid, 1)
    D = first.(D)
    B = [bdist(window, x) for x in grid]
    ind = D .<= B
    T = min.(D, B)
    km(ind, T, r)
end

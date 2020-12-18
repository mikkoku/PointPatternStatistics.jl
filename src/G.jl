"""
    Gkm(xy, r, window)

Return nearest neighbour distance distribution based on Kaplam-Meier estimator.
"""
function Gkm(pp, r)
    xy = pp.data
    if length(xy) == 1
        return zeros(Float64, length(r))
    end
    svec = @. SVector(getx(xy), gety(xy))
    tree = KDTree(svec)
    # Take 2 nearest neighbors because, the first nearest neighbor is the point itself
    _, D = knn(tree, svec, 2, true)
    D = last.(D)
    B = [bdist(window(pp), x) for x in xy]
    ind = D .<= B
    T = min.(D, B)
    km(ind, T, r)
end

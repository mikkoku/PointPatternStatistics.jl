"""
    Gkm(xy, r, window)

Return nearest neighbour distance distribution based on Kaplam-Meier estimator.
"""
function Gkm(pp, r)
    xy = pp.data
    svec = @. SVector(getx(xy), gety(xy))
    tree = KDTree(svec)
    _, D = knn(tree, svec, 2, true)
    D = last.(D)
    B = [bdist(window(pp), x) for x in xy]
    ind = D .<= B
    T = min.(D, B)
    km(ind, T, r)
end

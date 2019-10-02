function km(ind, T, r)
    perm = sortperm(vec(T))
    lastT = T[perm[1]]
    n = 0
    d = length(perm)
    prod = 1.0
    Fr = ones(length(r))
    ri = 1
    for (i, p) in enumerate(perm)
        if T[p] != lastT
            prod *= (1 - n/d)
            while ri <= length(r) && r[ri] <= T[p]
                Fr[ri] = 1-prod
                ri += 1
            end
            n = 0
            d = length(perm) - i + 1
            lastT = T[p]
        end
        n += ind[p]
    end
    prod *= (1 - n/d)
    while ri <= length(r) && r[ri] <= lastT
        Fr[ri] = 1-prod
        ri += 1
    end
    Fr
end

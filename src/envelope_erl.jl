function ranks(x)
    ranks!(similar(x, Int), x)
end
function ranks!(out, x, temp=similar(x, Float64))
    perm = temp
    sortperm!(perm, x)
    b = out #zeros(length(x))#Float64.(a)
    i = 0
    while i < length(x)
        i += 1
        x1 = x[perm[i]]
        j = i
        while j+1 <= length(x) && x1 == x[perm[j+1]]
            j += 1
        end
        for k = i:j
            b[perm[k]] = (i+j)/2
        end
        i = j
    end
    b
end

"""
    erlenvelope(A, alpha=0.05)

Return envelope with global significance level alpha based on
Extreme Rank Length.

"""
function erlenvelope(A, alpha=0.05)
    N = length(A)
    ndims = length(size(A[1])) + 1
    colons = ntuple(i -> :, ndims-1)
    A1 = Array{eltype(A[1]), ndims}(undef, size(A[1])..., N)
    for (i, a) in enumerate(A)
        A1[colons..., i] = a
    end
    temp1 = similar(A, Int)
    temp2 = similar(A, Float64)
    Ar = mapslices(x -> ranks!(temp2, x, temp1), A1, dims=ndims)

    @. Ar = min(Ar, N - Ar + 1)
    Aerl = sort!(reshape(Ar, :, N), dims=1)
    n = round(Int, alpha*length(A))
    inds = partialsortperm(collect(eachcol(Aerl)), n+1:N)

    hi = similar(A[1])
    lo = similar(hi)
    for i in eachindex(CartesianIndices(hi))
        lo[i], hi[i] = extrema(view(A1, i, inds))
    end
    lo, hi
end

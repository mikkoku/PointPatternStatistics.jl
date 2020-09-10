struct SimpleSequentialInhibition{T <: Real}
    minx::T
    maxx::T
    miny::T
    maxy::T
    R::T
end
""" SimpleSequentialInhibition(window, R)

Define a simple sequential inhibition (a.k.a random sequential adsorbtion) point
process in `window` with minimum interpoint distance `R`.

The rand function has keyword argument `T` specifying the maximum number of
uniform proposal points as `T*area(window)/R^2`.

The algorithm used for sampling from the process is based on (Wang, 1994) with a
few modifications. For the skipped uniform proposals only the expected count is
used instead of a random variable.

J.-S. Wang, A fast algorithm for random sequential adsorption of discs,
International Journal of Modern Physics C 1994 05:04, 707-715.

```julia
rand(SimpleSequentialInhibition((x=(0.0, 100.0), y=(0.0, 100.0)), 10.0), T=10.0^20)
```
"""
SimpleSequentialInhibition(window, R) =
    SimpleSequentialInhibition{Float64}(window.x[1], window.x[2], window.y[1], window.y[2], R)
const RandomSequentialAdsorbtion = SimpleSequentialInhibition


# Algorithm 1
function tooclose(p, p0, R)
    x0, y0 = p0
    R2 = Float64(R^2)
    for i1 in eachindex(p)
        x1, y1 = p[i1]
        d2 = (x1-x0)^2 + (y1-y0)^2
        if d2 < R2
            return true
        end
    end
    false
end
function randSSIfillsimple!(p, (maxx, maxy), R, T)
    t = 0.0
    dt = R^2 / (maxx*maxy)
    while t <= T
        t += dt
        p1 = (rand()*maxx, rand()*maxy)
        if !tooclose(p, p1, R)
            push!(p, p1)
        end
    end
    p
end
function randSSIsimple(window, R, T)
    p = randp(window, 1)
    randSSIfillsimple!(p, window, R, T)
end

# Algorithm 2 with a grid based search for nearby points.
function tooclose!(closepoints::Array{<:Integer, 3}, pp, (x0,y0), R)
    xi = unsafe_trunc(Int, x0/R) + 2
    yi = unsafe_trunc(Int, y0/R) + 2
    R2 = Float64(R^2)
    for vyi in yi-1:yi+1
        for vxi in xi-1:xi+1
            ii = 1
            while ii <= 3 && (pi = closepoints[ii, vxi, vyi]) > 0
                @inbounds x1, y1 = pp[pi]
                d2 = (x1-x0)^2 + (y1-y0)^2
                if d2 < R2
                    return true
                end
                ii += 1
            end
        end
    end
    i = 1
    while closepoints[i, xi, yi] > 0
        i += 1
    end
    closepoints[i, xi, yi] = length(pp) + 1
    return false
end
function initializeclosepoints(pp, (maxx,maxy), R)
    closepoints = fill(-1, (3, Int(cld(maxx, R)+2), Int(cld(maxy, R))+2))
    for (pi, (x,y)) in enumerate(pp)
        xi = Int(fld(x, R)) + 2
        yi = Int(fld(y, R)) + 2
        i = 1
        while closepoints[i, xi, yi] > 0
            i += 1
        end
        closepoints[i, xi, yi] = pi
    end
    closepoints
end
# Part of alg 3
function iscovered(x0, x1, y0, y1, discs, closepoints::Array{Int,3}, R, tmp)
    relevantdiscs, tmpintervals = tmp
    xi = unsafe_trunc(Int, (x0 + x1)/(2R)) + 2
    yi = unsafe_trunc(Int, (y0 + y1)/(2R)) + 2
    empty!(relevantdiscs)
    for vyi in yi-1:yi+1
        for vxi in xi-1:xi+1
            ii = 1
            while ii <= 3 && (pi = closepoints[ii, vxi, vyi]) > 0
                x = discs[pi]
                if isrelevant(x0, x1, y0, y1, x, R)
                    push!(relevantdiscs, x)
                end
                ii += 1
            end
        end
    end

    iscovered(x0, x1, y0, y1, relevantdiscs, R, tmpintervals)
end

function randSSIfill!(pp, (maxx, maxy), R, T)
    window = (x=(0, maxx), y=(0, maxy))
    closepoints = initializeclosepoints(pp, (maxx, maxy), R)
    t = 0
    dt = R^2 / (maxx*maxy)
    while t <= T
        t += dt
        p1 = (rand()*maxx, rand()*maxy)

        if !tooclose!(closepoints, pp, p1, R)
            push!(pp, p1)
        end
    end
    pp
end


function randSSI((maxx, maxy), R, T)
    p = [(rand()*maxx, rand()*maxy)]
    randSSIfill!(p, (maxx, maxy), R, T)
end

# Algorithm 3 with occupied cells
function isrelevant(x0, x1, y0, y1, (x, y), R)
    x0 <= x <= x1 && y0-R <= y <= y0+R && return true
    y0 <= y <= y1 && x0-R <= x <= x0+R && return true
    R2 = R^2
    (x0-x)^2 + (y0-y)^2 <= R2 && return true
    (x0-x)^2 + (y1-y)^2 <= R2 && return true
    (x1-x)^2 + (y1-y)^2 <= R2 && return true
    (x1-x)^2 + (y0-y)^2 <= R2 && return true
    return false
end
function iscovered_tmp()
    [Vector{NTuple{2, Float64}}() for _ in 1:4]
end
function iscovered(x0, x1, y0, y1, discs, R)
    iscovered(x0, x1, y0, y1, discs, R, iscovered_tmp())
end
function iscovered(x0, x1, y0, y1, discs, R, tmp)
    intervalss = tmp
    # 1.
    #discs = filter(x -> isrelevant(x0, x1, y0, y1, x, R), discs)
    length(discs) == 0 && return false

    # 2. Covered by a single disc
    R2 = R^2
    for (x,y) in discs
        (x0-x)^2 + (y0-y)^2 <= R2 &&
        (x0-x)^2 + (y1-y)^2 <= R2 &&
        (x1-x)^2 + (y1-y)^2 <= R2 &&
        (x1-x)^2 + (y0-y)^2 <= R2 && return true
    end

    # 3. Not all corners covered
    any(discs) do (x,y)
        (x0-x)^2 + (y0-y)^2 <= R2  end || return false
    any(discs) do (x,y)
        (x0-x)^2 + (y1-y)^2 <= R2  end || return false
    any(discs) do (x,y)
        (x1-x)^2 + (y1-y)^2 <= R2  end || return false
    any(discs) do (x,y)
        (x1-x)^2 + (y0-y)^2 <= R2  end || return false

    # 4.
    # first edge x0,y0 - x1,y0
    # y0 < y1
    sx = x1-x0
    sy = y1-y0
    foreach(empty!, intervalss)

    addinterval(i,a,b) = push!(intervalss[i], (a,b))
    for (x,y) in discs
        dx2 = R^2 - (y-y0)^2
        if dx2 >= 0.0
            dx = sqrt(dx2)
            a0 = x-x0-dx
            a1 = x-x0+dx
            if a0 <= sx && a1 >= 0
                addinterval(1, a0, a1)
            end
        end
        dx2 = R^2 - (y-y1)^2
        if dx2 >= 0.0
            dx = sqrt(dx2)
            a0 = x-x0-dx
            a1 = x-x0+dx
            if a0 <= sx && a1 >= 0
                addinterval(2, a0, a1)
            end
        end

        dx2 = R^2 - (x-x0)^2
        if dx2 >= 0.0
            dx = sqrt(dx2)
            a0 = y-y0-dx
            a1 = y-y0+dx
            if a0 <= sy && a1 >= 0
                addinterval(3, a0, a1)
            end
        end
        dx2 = R^2 - (x-x1)^2
        if dx2 >= 0.0
            dx = sqrt(dx2)
            a0 = y-y0-dx
            a1 = y-y0+dx
            if a0 <= sy && a1 >= 0
                addinterval(4, a0, a1)
            end
        end
    end
    for i in eachindex(intervalss)
        intervals = intervalss[i]
        # Use in-place sorting
        sort!(intervals, alg=Base.Sort.InsertionSort, by=first)
        b = 0.0
        @assert intervals[1][1] <= 0.0
        for (a0, a1) in intervals
            a0 <= b || return false
            b = max(a1, b)
        end
        @assert b >= ifelse(i <= 2, sx, sy)
    end
    length(discs) == 2 && return true

    # 5.
    for i in 1:length(discs)
        for j in i+1:length(discs)
            ax, ay = discs[i]
            bx, by = discs[j]
            d2 = (ax-bx)^2 + (ay-by)^2
            d2 <= 4R2 || continue
            d = sqrt(d2)
            cx = 0.5*(ax + bx)
            cy = 0.5*(ay + by)
            l = sqrt(R^2 - d2/4)
            fx = cx + l*(ay-by)/d
            fy = cy - l*(ax-bx)/d
            gx = cx - l*(ay-by)/d
            gy = cy + l*(ax-bx)/d
            # if f or g is not in the box it doesn't need to be checked
            fcovered = !(x0 <= fx <= x1 && y0 <= fy <= y1)
            gcovered = !(x0 <= gx <= x1 && y0 <= gy <= y1)
            for k in 1:length(discs)
                if k != i && k != j
                    x,y = discs[k]
                    fcovered |= (fx-x)^2 + (fy-y)^2 <= R2
                    gcovered |= (gx-x)^2 + (gy-y)^2 <= R2
                end
            end
            (fcovered && gcovered) || return false
        end
    end
    return true
end

function randSSIfill2!(pp, (maxx, maxy), R, T)
    window = (x=(0, maxx), y=(0, maxy))
    closepoints = initializeclosepoints(pp, (maxx, maxy), R)
    tmp = (Vector{NTuple{2, Float64}}(), iscovered_tmp())
    covered = [iscovered(x, x+R, y, y+R, pp, closepoints, R, tmp) for x in 0:R:maxx, y in 0:R:maxy]
    t = 0.0
    dt = R^2 / (maxx*maxy)
    while t <= T
        t += dt
        (x,y) = (rand()*maxx, rand()*maxy)
        p1 = (x,y)
        xi = unsafe_trunc(Int, x/R) + 1
        yi = unsafe_trunc(Int, y/R) + 1
        if !covered[xi, yi]
            tc = tooclose!(closepoints, pp, p1, R)
            if !tc
                push!(pp, p1)
            end
        end
    end
    pp
end
function build_notcovered(closepoints, pp, (maxx, maxy), R, tmp)
    notcovered = Vector{NTuple{2, Float64}}()
    for x in 0:R:maxx, y in 0:R:maxy
        if x < maxx && y < maxy && !iscovered(x, min(x+R, maxx), y, min(y+R, maxy), pp, closepoints, R, tmp)
            push!(notcovered, (x,y))
        end
    end
    notcovered
end
function rebuild_notcovered!(notcovered, closepoints, pp, (maxx, maxy), R, cellsize, tmp)
    filter!(notcovered) do (x,y)
        x < maxx && y < maxy &&
        !iscovered(x, min(x+cellsize, maxx), y, min(y+cellsize, maxy), pp, closepoints, R, tmp)
    end
    # Here is maybe a possibility of infinite loop if cellsize is too small
    maxs = max(maxx, maxy)
    if maxs == maxs + 0.5cellsize
        throw(("Cellsize got too small", cellsize))
    end
    if  0 < length(notcovered) <= 100
        l = length(notcovered)
        resize!(notcovered, 4*l)
        for i in l:-1:1
            x0, y0 = notcovered[i]
            x1 = x0+0.5cellsize
            y1 = y0+0.5cellsize
            notcovered[4i] = (x1, y1)
            notcovered[4i-1] = (x1, y0)
            notcovered[4i-2] = (x0, y0)
            notcovered[4i-3] = (x0, y1)
        end
        rebuild_notcovered!(notcovered, closepoints, pp, (maxx, maxy), R, 0.5cellsize, tmp)
    else
        cellsize
    end
end
function randSSIfill3!(pp, (maxx, maxy), R, T)
    window = (x=(0, maxx), y=(0, maxy))
    closepoints = initializeclosepoints(pp, (maxx, maxy), R)
    tmp = (Vector{NTuple{2, Float64}}(), iscovered_tmp())
    cellsize = R
    notcovered = build_notcovered(closepoints, pp, (maxx, maxy), cellsize, tmp)

    # proposal_area = (length(notcovered)*cellsize^2)
    # The expected time used before getting a proposal in the not covered squares is
    # total area / proposal_area
    # Normally the time increment is R^2 / total area
    dt = 1/length(notcovered)
    t = 0.0
    proposals = 0
    total = 0
    acc = 0

    isempty(notcovered) && return (status=:jamm, cellsize=cellsize, time=t, totalproposals=total)
    while t <= T
        t += dt
        proposals += 1
        total += 1
        if total > 1_000_000
            println(("Something is wrong", length(pp), length(notcovered), (cellsize, t, dt)))
            throw("asdf")
        end
        if 3acc >= length(notcovered) || proposals >= 10*length(notcovered)
            cellsize = rebuild_notcovered!(notcovered, closepoints, pp, (maxx, maxy), R, cellsize, tmp)
            dt = R^2 /(length(notcovered) * cellsize^2)
            if isempty(notcovered)
                return (status=:jamm, cellsize=cellsize, time=t, totalproposals=total)
            end
            proposals = 0
            acc = 0
        end
        (x0, y0) = rand(notcovered)
        (x,y) = (x0+rand()*cellsize, y0+rand()*cellsize)
        p1 = (x,y)

        if x <= maxx && y <= maxy && !tooclose!(closepoints, pp, p1, R)
            acc += 1
            push!(pp, p1)
        end
    end
    (status=:ok, cellsize=cellsize, time=t, totalproposals=total)
end

function Base.rand(d::SimpleSequentialInhibition; T=10.0^20, alg=:gridrand, info=(info, p) -> nothing)
    w = d.maxx-d.minx
    h = d.maxy-d.miny
    status = :ok
    if alg == :simple
        return randSSIsimple((x=(d.minx, d.maxx), y=(d.miny, d.maxy)), d.R, T)
    elseif alg == :gridsearch
        p = randSSI((w, h), d.R, T)
    elseif alg == :gridrand
        # https://arxiv.org/pdf/cond-mat/9402066.pdf
        p = randSSI((w, h), d.R, min(16.0, T))
        p = [(rand()*w, rand()*h)]
        status = randSSIfill3!(p, (w, h), d.R, T)
    else
        throw(ArgumentError("$alg is not a valid algorithm for SimpleSequentialInhibition."))
    end
    for i in eachindex(p)
        x,y = p[i]
        p[i] = (x+d.minx, y+d.miny)
    end
    info(status, p)
    p
end

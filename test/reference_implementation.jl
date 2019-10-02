# Reference implementations for summary functions.

function area(window)
    -(window.x...) * -(window.y...)
end

function pairdist(xy)
    d((a1, a2), (b1, b2)) = sqrt((a1-b1)^2 + (a2-b2)^2)
    d.(reshape(xy, 1, :), xy)
end
epanechnikov(t, h) = if -h <= t <= h 3/4*(1-(t/h)^2)/h else 0.0 end
function sumpairs(dist, A, r, window, f)
    # Use TwicePrecision because adding a lot of small numbers causes some
    # round off error.
    # Doesn't really matter because there are roundoff errors in other implementations.
    s = Base.TwicePrecision(0.0)
    for j in axes(dist, 2)
        for i in axes(dist, 1)
            if i != j
                s += f(r, dist[i,j])/A[i,j]
            end
        end
    end
    Float64(s)
end
function Atrans(window, d)
    d1, d2 = abs.(d)
    x1, x2 = window.x
    y1, y2 = window.y
    ((x2-x1-d1)*(y2-y1-d2))
end
function Apairs(xy, window)
    ((a, b) -> Atrans(window, a .- b)).(reshape(xy, 1, :), xy)
end
function pcfref(xy, window, r)
    λ = length(xy)/area(window)
    h = 0.15/sqrt(λ)
    dist = pairdist(xy)
    A = Apairs(xy, window)
    (r -> sumpairs(dist, A, r, window, (r, d) -> epanechnikov(r - d, h)) / (2 * pi * r) / λ^2 ).(r)
end


function Kref(xy, window, r)
    dist = pairdist(xy)
    A = Apairs(xy, window)
    n = length(xy)
    λ² = n*(n-1)/area(window)^2
    (r -> sumpairs(dist, A, r, window, (r, d) -> d < r) / λ² ).(r)
end

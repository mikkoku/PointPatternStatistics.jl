
randu((a, b)) = a + rand() * (b-a)
randw(window) = (randu(window.x), randu(window.y))

function area(window::NamedTuple{(:x,:y)})
    -(window.x...) * -(window.y...)
end

function edgetrans(window::NamedTuple{(:x,:y)}, d)
    d1, d2 = abs.(d)
    x1, x2 = window.x
    y1, y2 = window.y
    (x2-x1)*(y2-y1) / ((x2-x1-d1)*(y2-y1-d2))
end

function bdist(window::NamedTuple{(:x,:y)}, xy)
    x, y = getx(xy), gety(xy)
    (x1,x2), (y1, y2) = window.x, window.y
    min(x-x1, y-y1, x2-x, y2-y)
end

""" inside(point, window)

Is point inside or on the boundary of window.
"""
function inside(T::Type{<:Union{Bool, Int}}, point, window::NamedTuple{(:x,:y)})
    x1, x2 = window.x
    y1, y2 = window.y
    x, y = getx(point), gety(point)
    if T==Bool
        x1 <= x <= x2 && y1 <= y <= y2
    else
        x1 < x < x2 && y1 < y < y2 && return 1
        x1 <= x <= x2 && y1 <= y <= y2 && return 0
        -1
    end
end
inside(point, window) = inside(Bool, point, window)
inside(window) = Base.Fix2(inside, window)

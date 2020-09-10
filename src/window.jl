
randu((a, b)) = a + rand() * (b-a)
randw(window) = (randu(window.x), randu(window.y))

function area(window)
    -(window.x...) * -(window.y...)
end

function edgetrans(window, d)
    d1, d2 = abs.(d)
    x1, x2 = window.x
    y1, y2 = window.y
    (x2-x1)*(y2-y1) / ((x2-x1-d1)*(y2-y1-d2))
end

function bdist(window, (x,y))
    (x1,x2), (y1, y2) = window.x, window.y
    min(x-x1, y-y1, x2-x, y2-y)
end

""" inside(point, window)
"""
function inside(point, window)
    x1, x2 = window.x
    y1, y2 = window.y
    x, y = point
    x1 <= x <= x2 && y1 <= y <= y2
end
inside(window) = Base.Fix2(inside, window)

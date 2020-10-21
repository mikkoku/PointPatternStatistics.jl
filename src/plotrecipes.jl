@recipe function plot(pp::PointPattern; window=true)
    if window && pp.window isa NamedTuple{(:x, :y)}
        @series begin
            seriestype := :path
            primary := false
            linecolor := :black
            markershape := :none
            x0, x1 = pp.window.x
            y0, y1 = pp.window.y
            ([x0, x0, x1, x1, x0], [y0, y1, y1, y0, y0])
        end
    end
    seriestype --> :scatter
    markershape --> :+
    markercolor --> :auto
    (getx.(pp.data), gety.(pp.data))
end


function disthisttrans(xy, step, N, window, pcf=false)
    hist = zeros(N)
    d2((x, y)) = Float64(x^2 + y^2)
    Rmax = (N-1)*step
    R2max = Rmax^2
    for i in 1:length(xy)-1, j in i+1:length(xy)
        dxy = xy[i] .- xy[j]
        dd = d2(dxy)
        if dd <= R2max
            d = sqrt(dd)/step
            i = ceil(Int, d)+1
            frac = i - 1 - d
            @assert i <= N
            val = 2*edgetrans(window, dxy)
            if pcf && i > 1
                hist[i-1] += frac*val
                hist[i] += (1-frac)*val
            else
                hist[i] += val
            end
        end
    end
    hist
end

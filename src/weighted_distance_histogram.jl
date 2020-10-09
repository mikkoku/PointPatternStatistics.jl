
disthisttrans(xy, step, N, window, pcf=false) =
    disthisttransf(xy, (_,_) -> 1.0, step, N, window, pcf)
disthisttransf(xy, f, step, N, window, pcf=false) =
    disthisttrans12f(xy, xy, f, step, N, window, pcf)
disthisttrans12(xy1, xy2, step, N, window, pcf=false) =
    disthisttrans12f(xy1, xy2, (_,_) -> 1.0, step, N, window, pcf)
function disthisttrans12f(xy1, xy2, f, step, N, window, pcf=false)
    hist = zeros(N)
    Rmax = (N-1)*step
    T = promote_type(typeof(getx(xy1[1])), typeof(step))
    R2max = T(Rmax^2)
    mult = xy1 === xy2 ? 2.0 : 1.0

    for i1 in 1:length(xy1)
        i2start = xy1 === xy2 ? i1+1 : 1
        for i2 in i2start:length(xy2)
            #dxy = xy1[i1] .- xy2[i2]
            dx = getx(xy1[i1]) - getx(xy2[i2])
            dy = gety(xy1[i1]) - gety(xy2[i2])
            dd = T(dx^2 + dy^2)
            # dd = d2(dxy)
            if dd <= R2max
                d = sqrt(dd)/step
                k = ceil(Int, d)+1
                frac = k - 1 - d
                @assert k <= N
                val = mult*edgetrans(window, (dx, dy))*f(xy1[i1], xy2[i2])
                if pcf && k > 1
                    hist[k-1] += frac*val
                    hist[k] += (1-frac)*val
                else
                    hist[k] += val
                end
            end
        end
    end
    hist
end

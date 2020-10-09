using .RCall
import .RCall: rcopy, rcopytype, RClass
function rcopy(::Type{PointPattern}, s::Ptr{VecSxp})
    x = rcopy(Vector, s[:x])
    y = rcopy(Vector, s[:y])
    window = rcopy(Dict{Symbol, Any}, s[:window])

    PointPattern(collect(zip(x,y)), (x=window[:xrange], y=window[:yrange]))
end

rcopytype(::Type{RClass{:ppp}}, ::Ptr{VecSxp}) = PointPattern

import .RCall: sexp, protect, unprotect, setclass!, RClass, sexpclass
function sexp(::Type{RClass{:ppp}}, pp::PointPattern)
    # win = protect(sexp(Dict(:type => "regtangle",
    #     :xrange => [pp.window.x...],
    #     :yrange => [pp.window.y...])))
    # setclass!(win, sexp("owin"))
    # r = protect(sexp(Dict(
    #     :n => length(pp.data),
    #     :x => getx.(pp.data),
    #     :y => gety.(pp.data),
    #     :window => win)))
    # setclass!(r, sexp("ppp"))
    # unprotect(2)
    r = R"spatstat::ppp"(getx.(pp.data), gety.(pp.data),
        window=R"spatstat::owin"([pp.window.x...], [pp.window.y...]))
    r
end

sexpclass(f::PointPattern) = RClass{:ppp}

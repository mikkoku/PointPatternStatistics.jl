# PointPatternStatistics.jl

A small collection of tools for model evaluation in point pattern statistics.

For data analysis there is an extensive R-package
[spatstat](http://spatstat.org).

Currently implemented
- K and pcf with "translate" edge correction
- F and G using Kaplan-Meier estimator
- Global envelope based on Extreme Rank Length

## TODO
- Point pattern type, window type
- Useful with logpdf and rand

```
using PointPatternStatistics
window = (x=(0,1), y=(0,1))
xy = [(rand(), rand()) for _ in 1:300]
r = range(0, 0.25, length=512)
pcf(xy, window, r)
Kest(xy, window, r)
Fest(xy, window, r)
Gest(xy, window, r)

xys = [[(rand(), rand()) for _ in 1:30] for _ in 1:2499]
A = [pcf(xy, window, r) for xy in xys]
globalenvelope(A)
```

## Bibliography

Stoyan, D., & Stoyan, H. (1994). Fractals, random shapes, and point fields: methods of geometrical statistics (Vol. 302). John Wiley & Sons Inc.

Adrian Baddeley, Ege Rubak, Rolf Turner (2015). Spatial Point Patterns: Methodology and Applications
with R. London: Chapman and Hall/CRC Press, 2015. URL
http://www.crcpress.com/Spatial-Point-Patterns-Methodology-and-Applications-with-R/Baddeley-Rubak-Turner/9781482210200/

Baddeley, A.J. and Gill, R.D. Kaplan-Meier estimators of interpoint distance distributions for spatial point processes. Annals of Statistics 25 (1997) 263-292.

Myllymäki, M., Mrkvicka, T., Grabarnik, P., Seijo, H. and Hahn, U. (2017). Global envelope tests for
  spatial processes. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 79:
  381-404. doi: 10.1111/rssb.12172

include("PolygonBase.jl")
include("silhouette&pieces.jl")
plotly()

P = tri_l
H = tri_s

plot(make_hole_GEOS(P, H))
include("PolygonBase.jl")
include("silhouette&pieces.jl")
plotly()

polygons = [tri_s]

display(polygons..., vertex=false)
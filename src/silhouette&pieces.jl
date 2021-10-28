tri_s = MYPolygon([0, 1, 1], [0, 0, 1])
tri_m = MYPolygon([0, √2, √2], [0, 0, √2])
tri_l = MYPolygon([0, 2, 2], [0, 0, 2])
square_s = MYPolygon([0, 1, 1, 0], [0, 0, 1, 1])
parallelogram = MYPolygon([0, √2, √(2)*3/2, √(2)/2], [0, 0, √(2)/2, √(2)/2])

# tri_m, square_s
house2 = MYPolygon([0, 1, 1, 1.5, 0.5, -0.5, 0], [0, 0, 1, 1, 2, 1, 1])

# tri_l, tri_m, tri_s, tri_s
house4 = MYPolygon([-0.5*√2, 0.5*√2, 0.5*√2, √2, 0, -√2, -0.5*√2], [0, 0, √2, √2, 2*√2, √2, √2])

# tri_s, tri_s, parallelogram
# tri_s, tri_s, square_s
rect_s = MYPolygon([0, 2, 2, 0], [0, 0, 1, 1])

#tri_s, tri_s, parallelogram
hurt_s = MYPolygon([0, 1, 1, 0, -1, -1], [0, 1, 2, 1, 2, 1])

# tri_m, tri_s, tri_s, parallelogram
hexagon_m = MYPolygon([0, √2, 2.12132, √2, 0, -0.7071068], [0, 0, 0.7071068, √2, √2, 0.7071068])

# tri_m, tri_s, tri_s
square_m = MYPolygon([0, √2, √2, 0], [0, 0, √2, √2])

# tri_l, tri_l, tri_m, tri_s, tri_s, square_s, parallelogram
square_l = MYPolygon([0, 2*√2, 2*√2, 0], [0, 0, 2*√2, 2*√2])

# tri_l, tri_l
betterfly = MYPolygon([0, √2, 2*√2, 2*√2, √2, 0], [0, √2, 0, 2*√2, √2, 2*√2])

# tri_s, square_s, parallelogram
note = MYPolygon([0, 0.5*√2, 0.5*√2, 1.70, 0, 0, -0.5*√2], [-0.5*√2, 0, 0.413, 0.413, 1.5*√2, 0.5*√2, 0])

# tri_l, tri_s, parallelogram
diamond3 = MYPolygon([0, √2, 0.5*√2, -0.5*√2, -√2], [-√2, 0, 0.5*√2, 0.5*√2, 0])

# tri_l, tri_s, tri_s, square_s
antena = MYPolygon([-√2, √2, 0, 0.5*√2, 0.5*√2, -√2, 0], [0, 0, √2, 1.5*√2, 3.53, √2, √2])

# tri_l, tri_l, tri_m, tri_s, tri_s, square_s, parallelogram
katana = MYPolygon([-1, 0, 1, 2, 1.22, 5.63, 5.63, 4.62, 4.62, 3.20, 3.20, 1.21, 1.21, 0, 0, -1], 
                   [0, 0, 1, 1, 1.79, 6.21, 7.21, 6.22, 5.21, 5.21, 3.79, 3.79, 1.79, 3, 2, 1])
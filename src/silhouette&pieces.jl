tri_s = MYPolygon([0, 1, 1], [0, 0, 1])
tri_m = MYPolygon([0, √2, √2], [0, 0, √2])
tri_l = MYPolygon([0, 2, 2], [0, 0, 2])
square_s = MYPolygon([0, 1, 1, 0], [0, 0, 1, 1])
parallelogram = MYPolygon([0, √2, √(2)*3/2, √(2)/2], [0, 0, √(2)/2, √(2)/2])

# tri_m, square_s
house2 = MYPolygon([0, 1, 1, 1.5, 0.5, -0.5, 0], [0, 0, 1, 1, 2, 1, 1])
p_house2 = ("house2", house2, [tri_m, square_s])

# tri_l, tri_m, tri_s, tri_s
house4 = MYPolygon([-0.5*√2, 0.5*√2, 0.5*√2, √2, 0, -√2, -0.5*√2], [0, 0, √2, √2, 2*√2, √2, √2])
p_house4 = ("house4", house4, [tri_l, tri_m, tri_s, tri_s])

# tri_l, square_s
arrow2 = MYPolygon([0, 3, 1, 1, 0], [0, 0, 2, 1, 1])
p_arrow2 = ("arrow2", arrow2, [tri_l, square_s])

# tri_s, parallelogram
trapezoid2_1 = MYPolygon([-√2, √2, 0.5*√2, -0.5*√2], [0, 0, 0.5*√2, 0.5*√2])
p_trapezoid2_1 = ("trapezoid2_1", trapezoid2_1, [tri_s, parallelogram])

# tri_s, parallelogram
trapezoid2_2 = MYPolygon([0, 0.5*√2, 0, -√2], [0, 0.5*√2, √2, 0])
p_trapezoid2_2 = ("trapezoid2_2", trapezoid2_2, [tri_s, parallelogram]) 

# tri_s, square_s
tri_above_square = MYPolygon([0, 1, 1, 0, 0, 1, 0, 0, 0], [0, 0, 1, 1, 2, 3, 3, 2, 1])
p_tri_above_square = ("tri_above_square", tri_above_square, [tri_s, square_s])

# tri_s, square_s
fish2 = MYPolygon([0, 0.5*√2, √2, 0.5*√2, 0, -0.5*√2, -0.5*√2], [0, -0.5*√2, 0, 0.5*√2, 0, 0.5*√2, -0.5*√2])
p_fish2 = ("fish2", fish2, [tri_s, square_s])

# tri_m, tri_s
small_stand = MYPolygon([0, 1, 1, √2, 0], [0, 0, 1, √2, √2])
p_small_stand = ("small_stand", small_stand, [tri_m, tri_s])

# tri_m, tri_s, tri_s
triforce3 = MYPolygon([0, -1, -0.5*√2 - 1, 0.5*√2 - 1, -1, 1, 1 - 0.5*√2, 1 + 0.5*√2, 1], [1, 0, -0.5*√2, -0.5*√2, 0, 0, -0.5*√2, -0.5*√2, 0])
p_triforce3 = ("triforce3", triforce3, [tri_m, tri_s, tri_s])

# tri_s, tri_s, parallelogram
# tri_s, tri_s, square_s
rect_s = MYPolygon([0, 2, 2, 0], [0, 0, 1, 1])
p_rect_s = ("rect_s", rect_s, [tri_s, tri_s, parallelogram])

# tri_s, tri_s, parallelogram
hurt3 = MYPolygon([0, 1, 1, 0, -1, -1], [0, 1, 2, 1, 2, 1])
p_hurt3 = ("hurt3", hurt3 , [tri_s, tri_s, parallelogram])

# tri_m, tri_s, tri_s, parallelogram
# tri_m, tri_s, tri_s, square_s
hexagon_m = MYPolygon([0, √2, 2.12132, √2, 0, -0.7071068], [0, 0, 0.7071068, √2, √2, 0.7071068])
p_hexagon_m1 = ("hexagon_m", hexagon_m, [tri_m, tri_s, tri_s, parallelogram])
p_hexagon_m2 = ("hexagon_m", hexagon_m, [tri_m, tri_s, tri_s, square_s])

# tri_m, tri_s, tri_s
square_3 = MYPolygon([0, √2, √2, 0], [0, 0, √2, √2])
p_square_3 = ("square_3", square_3, [tri_m, tri_s, tri_s])

# tri_l, tri_s, tri_s, square_s
# tri_l, tri_m, tri_s, square_s
# tri_l, tri_s, tri_s, parallelogram
square_4 = MYPolygon([0, 2, 2, 0], [0, 0, 2, 2])
p_square_4 = ("square_4", square_4, [tri_l, tri_s, tri_s, square_s])

# tri_l, tri_l, tri_m, tri_s, tri_s, square_s, parallelogram
square_7 = MYPolygon([0, 2*√2, 2*√2, 0], [0, 0, 2*√2, 2*√2])
p_square_7 = ("square_7", square_7, [tri_l, tri_l, tri_m, tri_s, tri_s, square_s, parallelogram])

# tri_m, tri_s, tri_s
butterfly_s = MYPolygon([0, -1, 1, 0, 1, -1], [0, -1, -1, 0, 1, 1])
p_butterfly_s = ("butterfly_s", butterfly_s, [tri_m, tri_s, tri_s])

# tri_l, tri_l
butterfly = MYPolygon([0, √2, 2*√2, 2*√2, √2, 0], [0, √2, 0, 2*√2, √2, 2*√2])
p_butterfly = ("butterfly", butterfly, [tri_l, tri_l])

# tri_l, tri_l
pesudo_butterfly = MYPolygon([-1.16, -0.25, 1.16, 1.16, 0.25, -1.16], [-1.16, -0.25, -1.67, 1.16, 0.25, 1.67])
p_pesudo_butterfly = ("pesudo_butterfly", pesudo_butterfly, [tri_l, tri_l])

# tri_l, tri_l, tri_s, tri_s
butterfly4 = MYPolygon([-√2, -0.5*√2, -0.5*√2, -√2, √2, 0.5*√2, 0.5*√2, √2], [√2, 0.5*√2, -0.5*√2, -√2, -√2, -0.5*√2, 0.5*√2, √2])
p_butterfly4 = ("butterfly4", butterfly4, [tri_l, tri_l, tri_s, tri_s])

# tri_s, square_s, parallelogram
note = MYPolygon([0, 0.5*√2, 0.5*√2, 1.70, 0, 0, -0.5*√2], [-0.5*√2, 0, 0.413, 0.413, 1.5*√2, 0.5*√2, 0])
p_note = ("note", note, [tri_s, square_s, parallelogram])

# tri_l, tri_s, parallelogram
diamond3 = MYPolygon([0, √2, 0.5*√2, -0.5*√2, -√2], [-√2, 0, 0.5*√2, 0.5*√2, 0])
p_diamond = ("diamond", diamond3, [tri_l, tri_s, parallelogram])

# tri_l, tri_s, tri_s, square_s
antena = MYPolygon([-√2, √2, 0, 0.5*√2, 0.5*√2, -√2, 0], [0, 0, √2, 1.5*√2, 3.53, √2, √2])
p_antena = ("antena", antena, [tri_l, tri_s, tri_s, square_s])

# tri_l, tri_m, tri_s, square_s
pencile4 = MYPolygon([0, √2, √2, 0.5*√2, 0], [0, 0, 2*√2, 2.5*√2, 2*√2])
p_pencile4 = ("pencile4", pencile4, [tri_l, tri_m, tri_s, square_s])

# tri_l, tri_s, tri_s
# tri_m, tri_s, tri_s, square_s
# tri_m, tri_s, tri_s, parallelogram
crown4 = MYPolygon([-√2, √2, √2, 0.5*√2, 0, -0.5*√2, -√2], [0, 0, √2, 0.5*√2, √2, 0.5*√2, √2])
p_crown = ("crown", crown4, [tri_l, tri_s, tri_s])

# tri_l, tri_l, tri_s, tri_s
snake = MYPolygon([-2, -1, -1, 1, 1, 0, 0, -2], [-2, -2, -1, -1, 1, 1, 0, 0])
p_snake = ("snake", snake, [tri_l, tri_l, tri_s, tri_s])

# tri_m, tri_s, tri_s, square_s, parallelogram
hexagon_5 = MYPolygon([-1.5, -0.5, 0.5, 1.5, 0.5, -0.5], [0, -1, -1, 0, 1, 1])
p_hexagon5 = ("hexagon5", hexagon_5, [tri_m, tri_s, tri_s, square_s, parallelogram])

# tri_s, tri_s, tri_s, tri_s, tri_s, tri_s
egu_force = MYPolygon([-1, 0, -1, 0, 0, 1, -1, 0, 0, 1, 1, 2, -1], [-2, -1, -1, 0, -1, 0, 0, 1, 0, 1, 0, 1, 1])
p_egu_force = ("egu_force", egu_force, [tri_s, tri_s, tri_s, tri_s, tri_s, tri_s])

# tri_l, tri_l, tri_m, tri_s, tri_s, square_s, parallelogram
faucet = MYPolygon([-2, -1, -1, 2, 2, 1, 1, 2, -1, 0, 0, -2], [-2, -2, -1, -1, 0, 0, 1, 2, 2, 1, 0, 0])
p_faucet = ("faucet", faucet, [tri_l, tri_l, tri_m, tri_s, tri_s, square_s, parallelogram])

# tri_l, tri_m, tri_s, square_s, parallelogram
Taiikusuwaring_human = MYPolygon([1, -0.5*√2, -√2, -0.5*√2, -2, -√2, -1.5*√2, -0.5*√2, 1-√2, 1], [0, 1+0.5*√2, 1, 0.291, -1, -1, -1-0.5*√2, -1-0.5*√2, -√2, √2])
p_Taiikusuwari = ("Taiikusuwari", Taiikusuwaring_human, [tri_l, tri_m, tri_s, square_s, parallelogram])

# tri_l, tri_l, tri_m, tri_s, tri_s, square_s, parallelogram
crab_scissor = MYPolygon([0, √2, √2, 1.5*√2, 0.5*√2, 0.5*√2, 0.5*√2-1, 0.5*√2-1, -1-0.5*√2, -0.5*√2, -0.5*√2, -√2, -0.5*√2], [0, 0, √2, 1.5*√2, 2.5*√2, 1+0.5*√2, 1+0.5*√2, 1+2.5*√2, 1+1.5*√2, 1.5*√2, 0.5*√2, 0, -0.5*√2])
p_crab_scissor = ("crab_scissor", crab_scissor, [tri_l, tri_l, tri_m, tri_s, tri_s, square_s, parallelogram])

# tri_l, tri_l, tri_m, tri_s, tri_s, square_s, parallelogram
candy7 = MYPolygon([-1.5-√2, -1.5, -0.5, 0.5, 1.5, 1.5+√2, 1.5+√2, 1.5, 0.5, -0.5, -1.5, -1.5-√2], [-√2, 0, -1, -1, 0, -√2, √2, 0, 1, 1, 0, √2])
p_candy7 = ("candy7", candy7, [tri_l, tri_l, tri_m, tri_s, tri_s, square_s, parallelogram])

# ミスアリ tri_l, tri_l, tri_m, tri_s, tri_s, square_s, parallelogram
flame6 = MYPolygon([0, -0.5, -0.5, -1, 0, -√2, -0.5*√2, -0.5*√2, -√2, √2, 0.5*√2, 0.5*√2, √2, 0, 1, 0.5, 0.5], [2, 2, 1, 1, 0, 0, -0.5*√2, -1.5*√2, -2*√2, -2*√2, -1.5*√2, -0.5*√2, 0, 0, 1, 1, 2])
# p_flame6 = ("flame6", flame6, [tri_l, tri_l, tri_m, tri_s, tri_s, square_s, parallelogram])

# tri_l, tri_l, tri_m, tri_s, tri_s, square_s, parallelogram
flame7 = MYPolygon([1, 0, 0, -0.5, -0.5, -1, 0, -√2, -0.5*√2, -0.5*√2, -√2, √2, 0.5*√2, 0.5*√2, √2, 0, 1, 0.5, 0.5, 0, 1], [4, 3, 2, 2, 1, 1, 0, 0, -0.5*√2, -1.5*√2, -2*√2, -2*√2, -1.5*√2, -0.5*√2, 0, 0, 1, 1, 2, 2, 3])
p_flame7 = ("flame7", flame7, [tri_l, tri_l, tri_m, tri_s, tri_s, square_s, parallelogram])

# tri_l, tri_l, tri_m, tri_s, tri_s, square_s, parallelogram
hurt7 = MYPolygon([-2, -1, -1, -0.294, 0, 2, √2, 0.5*√2, 1, 1, 0, 0, -0.5*√2, -√2, -0.5*√2, -√2], [0, -1, -1-√2, -1-0.5*√2, -2, 0, 0, 0.5*√2, 1, 2, 1, √2, 1.5*√2, √2, 0.5*√2, 0])
p_hurt7 = ("hurt7", hurt7, [tri_l, tri_l, tri_m, tri_s, tri_s, square_s, parallelogram])

# tri_l, tri_l, tri_m, tri_s, tri_s, square_s, parallelogram
katana = MYPolygon([-1, 0, 1, 2, 1.22, 5.63, 5.63, 4.62, 4.62, 3.20, 3.20, 1.21, 1.21, 0, 0, -1], 
                   [0, 0, 1, 1, 1.79, 6.21, 7.21, 6.22, 5.21, 5.21, 3.79, 3.79, 1.79, 3, 2, 1])
p_katana = ("katana", katana, [tri_l, tri_l, tri_m, tri_s, tri_s, square_s, parallelogram])
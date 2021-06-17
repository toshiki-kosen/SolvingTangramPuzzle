using LibGEOS
using Plots
gr() # 図が早くて綺麗に書けるおまじない

mutable struct MYPolygon{T}
    vertexes::Matrix{T}

    n::Int64
    center::Vector{T}
end
 
# 多角形の重心
function calc_center(X, Y)
    A = zero(X[1])
    for i in 1:length(X)-1
        A += X[i] * Y[i+1] - X[i+1] * Y[i]
    end
    A += X[end] * Y[1] - X[1] * Y[end]
    A /= 2

    Cx = Cy = 0
    for i in 1:length(X)-1
        Cx += (X[i] + X[i+1]) * (X[i] * Y[i+1] - X[i+1] * Y[i])
        Cy += (Y[i] + Y[i+1]) * (X[i] * Y[i+1] - X[i+1] * Y[i])
    end
    Cx += (X[end] + X[1]) * (X[end] * Y[1] - X[1] * Y[end])
    Cy += (Y[end] + Y[1]) * (X[end] * Y[1] - X[1] * Y[end])

    return [Cx/(6A), Cy/(6A)]
end

# Figure の宣言用関数
# X x座標の配列, Y y座標の配列
function MYPolygon(X::Array, Y::Array, T=Float64)
    MYPolygon{T}(
        hcat(X, Y)', 
        length(X),
        calc_center(X, Y)
    )
end

# M (N*2)行列
function MYPolygon(M::Matrix, T=Float64)
    return MYPolygon(M[1, :], M[2, :], T)
end

function Base.copy(P::MYPolygon{T}) where T
    return MYPolygon{T}(P.vertexes, P.n, P.center)
end

# Polygon を (x, y) ベクトルだけ動かす
function move!(P::MYPolygon{T}, x::T, y::T) where T
    P.vertexes .+= [x, y]
    P.center += [x, y]
    return nothing
end

function move(P::MYPolygon{T}, x::T, y::T) where T
    newP = deepcopy(P)
    move!(newP, x, y)
    return newP
end

# Polygon を θ だけ回転させる
function rotate!(P::MYPolygon{T}, θ::T) where T
    R = [cos(θ) -sin(θ); sin(θ) cos(θ)]
    
    P.vertexes.-= P.center
    P.vertexes = R * P.vertexes
    P.vertexes .+= P.center
    return nothing
end

function rotate(P::MYPolygon{T}, θ::T) where T
    newP = deepcopy(P)
    rotate!(newP, θ)
    return newP
end

# 図形を描画する
function display(Polygons::MYPolygon...; center=false, vertex=false)
    pl = plot(xlim = (-1, 6), ylim = (-1, 6), size=(400, 400))

    for P in Polygons
        # 図形が閉じるように終点を追加
        V = hcat(P.vertexes, P.vertexes[:, 1])
        plot!(pl, V[1, :], V[2, :], fill=0)
        if center != false
            scatter!(pl, [P.center[1]], [P.center[2]], marker=center)
        end
        if vertex != false
            scatter!(pl, V[1, :], V[2, :], marker=vertex)
        end
    end
    return pl
end

function MYPolygon2LibGEOS(P::MYPolygon)
    str = "POLYGON(("
    for i in 1:size(P.vertexes)[2]
        str *= string(P.vertexes[1, i])
        str *= " "
        str *= string(P.vertexes[2, i])
        str *= ","
    end
    str *= string(P.vertexes[1, 1])
    str *= " "
    str *= string(P.vertexes[2, 1])
    str *= "))"
    return readgeom(str)
end

# 図形と図形の共通部分
function intersectionArea(P1::MYPolygon, P2::MYPolygon)
    geP1 = MYPolygon2LibGEOS(P1)
    geP2 = MYPolygon2LibGEOS(P2)
    gePi = LibGEOS.intersection(geP1, geP2)
    return LibGEOS.area(gePi)
end

function unionArea(P1::MYPolygon, P2::MYPolygon)
    geP1 = MYPolygon2LibGEOS(P1)
    geP2 = MYPolygon2LibGEOS(P2)
    gePu = LibGEOS.union(geP1, geP2)
    return LibGEOS.area(gePu)
end
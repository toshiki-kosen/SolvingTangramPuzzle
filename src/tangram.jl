using Plots
gr() # 図が早くて綺麗に書けるおまじない

mutable struct Polygon{T}
    vertexes::Matrix{T}

    n::Int64
    center::Vector{T}
end

# 多角形の重心
function calc_center(X, Y)
    A = 0
    for i in 1:length(X)-1
        A += X[i]*Y[i+1] - X[i+1] * Y[i]
    end
    A += X[1] * Y[end] - X[end] * Y[1]
    A /= 2

    Cx = Cy = 0
    for i in 1:length(X)-1
        Cx += (X[i] + X[i+1]) * (X[i]*Y[i+1] - X[i+1]*Y[i])
        Cy += (Y[i] + Y[i+1]) * (X[i]*Y[i+1]- X[i+1]*Y[i])
    end
    Cx += (X[end] + X[1]) * (X[end]*Y[1] - X[1]*Y[end])
    Cy += (Y[end] + Y[1]) * (X[end]*Y[1] - X[1]*Y[end])

    return [Cx/(6A), Cy/(6A)]
end

# Figure の宣言用関数
# X x座標の配列, Y y座標の配列
function Polygon(X::Array, Y::Array, T=Float64)
    Polygon{T}(
        hcat(X, Y)', 
        length(X),
        calc_center(X, Y)
    )
end

# M (2 x N)行列
function Polygon(M::Matrix, T=Float64)
    return Polygon(M[1, :], M[2, :], T)
end

# Polygon を (x, y) ベクトルだけ動かす
function move!(P::Polygon{T}, x::T, y::T) where T
    P.vertexes .+= [x, y]
    P.center += [x, y]
end

# Polygon を θ だけ回転させる
function rotate!(P::Polygon{T}, θ::T) where T
    R = [cos(θ) -sin(θ); sin(θ) cos(θ)]
    
    P.vertexes.-= P.center
    P.vertexes = R * P.vertexes
    P.vertexes .+= P.center
end

# 図形を描画する
function display(P::Polygon)
    # 図形が閉じるように終点を追加
    V = hcat(P.vertexes, P.vertexes[:, 1])
    pl = plot(V[1, :], V[2, :], xlim = (-1, 3), ylim = (-1, 3), size=(500, 500))
    scatter!(pl, [P.center[1]], [P.center[2]])
    return pl
end

function line_intersect(A::Vector, B::Vector, C::Vector, D::Vector)
    AB = B - A
    AC = C - A
    AD = D - A
    CD = C - D
    CA = -AC
    CB = B - C
    BD = D - B

    s1 = AB[1] * AC[2] - AB[2] * AC[1] # AB×AC
    t1 = AB[1] * AD[2] - AB[2] * AD[1] # AB×AD

    s2 = CD[1] * CA[2] - CD[2] * CA[1] # CD×CA
    t2 = CD[1] * CB[2] - CD[2] * CB[1] # CD×CB

    p = [0.0, 0.0]

    if s1*t1 < 0 && s2*t2 < 0
        λ = (CD[1] * AD[2] - CD[2] * AD[1])/(CD[1] * AB[2] - CD[2] * AB[1])
        p = A + λ * AB
    end
    return (s1*t1 < 0 && s2*t2 < 0), p
end

function is_crossing_xray(P1::Polygon, P2::Polygon, i::Int64, j::Int64)
    if P2.vertexes[2, i] < P1.vertexes[2, 1] < P2.vertexes[2, j]
        return 1
    elseif P2.vertexes[2, i] > P1.vertexes[2, 1] > P2.vertexes[2, j]
        return 1
    elseif P2.vertexes[2, i] == P1.vertexes[2, 1] && P1.vertexes[2, 1] < P2.vertexes[2, j]
        return 1
    elseif P2.vertexes[2, j] == P1.vertexes[2, 1] && P1.vertexes[2, 1] < P2.vertexes[2, i]
        return 1
    end
end

# 図形と図形の共通部分
function intersect(P1::Polygon, P2::Polygon)
    # 多角形の交点をpoly1, poly2 にそれぞれ追加
    # N^2 なので改善したい
    # Intersections = [(P1における番号, P2における番号, 座標, 使用したかどうか), ...]
    Intersections = Array{Tuple{Int64, Int64, Vector, Bool}}()
    for i in 1:P1.n-1
        for j in 1:P2.n-1
            is_inter, intersection = line_intersect(P1.vertexes[i], P1.vertexes[i+1], P2.vertexes[j], P2.vertexes[j+1])
            if is_inter
                push!(Intersections, (i+length(Intersections)+1, j+length(Intersections)+1, intersection, false))
            end
        end
    end

    num_inter = length(Intersections)

    # 交点が追加された頂点の集合
    npoly1 = P1.vertexes[1:Intersections[1][1] - 1]
    npoly2 = P2.vertexes[1:Intersections[1][2] - 1]

    for i in 1:num_inter-1
        index1 = Intersections[i][1] - 1
        index2 = Intersections[i][2] - 1
        hcat!(npoly1, Intersections[i][3])
        hcat!(npoly2, Intersections[i][3])

        hcat!(npoly1, P1.vertexes[Intersections[i][1]-i:Intersections[i+1][1] - i - 1])
        hcat!(npoly2, P1.vertexes[Intersections[i][2]-i:Intersections[i+1][2] - i - 1])
    end
    hcat!(npoly1, Intersections[end][3])
    hcat!(npoly2, Intersections[end][3])

    hcat!(npoly1, P1.vertexes[Intersections[end][1]-num_inter:end])
    hcat!(npoly2, P1.vertexes[Intersections[end][2]-num_inter:end])

    # P1 の１つ目の頂点が P2 の内点か調べる
    isonEdge = false
    counter = 0
    for i in 1:P2.n-1
        counter += is_crossing_xray(P1, P2, i, i+1)
    end
    counter += is_crossing_xray(P1, P2, P2.n, 1)

    isonEdge = counter % 2 == 1 ? true : false

    # poly1 の１つ目の頂点から巡回しながら交点が来たらそこから poly2 を巡回する
    # (共通部分を構成していたら poly1 と poly2 を行ったり来たりして、そうでなかったそうしない)
    # TODO; 共通部分が連結でない場合があるので、連結成分を個別に出力できるようにする
    counter = 0
    num_inter = length(Intersections)
    onP1 = true
    index = 1


    polygons = Array{typeof(P1)}()
    vertexes = Matrix{(typeof(P1.vertexes[1]))}()
    while counter < num_inter
        if isonEdge
            if onP1
                hcat!(vertexes, poly1.vertexes[index])
            else
                hcat!(vertexes, poly2.vertexes[index])
            end
            counter += 1
        end

        if onP1
            index = index < length(poly1) ? index + 1 : 1
        else
            index = index < length(poly2) ? index + 1 : 1
        end

        is_polygon_complete = false
        if onP1
            for i in 1:length(vertexes[1, :])
                if poly1[index] == vertexes[:, i]
                    push!(polygons, Polygon(vertexes))
                    empty!(vertexes)
                    is_polygon_complete = true
                end
            end
        else
            for i in 1:length(vertexes[1, :])
                if poly2[index] == vertexes[:, i]
                    push!(polygons, Polygon(vertexes))
                    empty!(vertexes)
                    is_polygon_complete = true
                end
            end
        end

        # スタート地点を変える
        if is_polygon_complete; continue end

        intersection_index = 0
        onIntersection = false
        if onP1
            for i in Intersections
                if i[1] == index
                    onIntersection = true
                    intersection_index = i[2]
                    break
                end
            end
        else
            for i in Intersections
                if i[2] == index
                    onIntersection = true
                    intersection_index = i[1]
                    break
                end
            end
        end

        if isonEdge
            onP1 = !onP1
            index = intersection_index == 0 ? index : intersection_index
        else
            isonEdge= true
        end
    end

    return polygons
end
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

# M (N*2)行列
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
function display(Polygons::Polygon...; center=false, vertex=false)
    pl = plot(xlim = (-1, 3), ylim = (-1, 3), size=(400, 400))

    for P in Polygons
        # 図形が閉じるように終点を追加
        V = hcat(P.vertexes, P.vertexes[:, 1])
        plot!(pl, V[1, :], V[2, :])
        if center != false
            scatter!(pl, [P.center[1]], [P.center[2]], marker=center)
        end
        if vertex != false
            scatter!(pl, V[1, :], V[2, :], marker=vertex)
        end
    end
    return pl
end

function ×(v::Vector{T}, u::Vector{T}) where T
    return v[1] * u[2] - v[2] * u[1]
end

function segment_intersect(A::Vector, B::Vector, C::Vector, D::Vector)
    s1 = (B - A) × (C - A) # AB×AC
    t1 = (B - A) × (D - A)# AB×AD

    s2 = (D - C) × (A - C) # CD×CA
    t2 = (D - C) × (B - C) # CD×CB

    p = [-10.0, -10.0]

    if s1*t1 < 0 && s2*t2 < 0
        AB = B - A
        DC = C - D
        M = [AB[1] DC[1]
             AB[2] DC[2]]
        λ = inv(M) * (C - A)
        p = A + λ[1] * AB
    end
    return (s1*t1 < 0 && s2*t2 < 0), p
end

function is_crossing_XrayFromP1FirstPoint(p1::Vector, P2::Polygon, i::Int64, j::Int64)
    is_intersect, intersection = segment_intersect(p1, [2^13, p1[2]], P2.vertexes[:, i], P2.vertexes[:, j])
    if P2.vertexes[:, i] == intersection || P2.vertexes[:, j] == intersection
        if P2.vertexes[2, i] > P2.vertexes[2, j]
            return true
        end
    end

    return is_intersect
end

function next_walker(walker::Tuple{Int64, Int64, Bool}, P1::Polygon, P2::Polygon)
    if walker[1] == 1
        ind = walker[2] + 1 <= P1.n ? walker[2] + 1 : 1
    else
        ind = walker[2] + 1 <= P2.n ? walker[2] + 1 : 1
    end

    return (walker[1], ind, walker[3])
end

# 図形と図形の共通部分
function intersect(P1::Polygon, P2::Polygon)
    # 多角形の交点をpoly1, poly2 にそれぞれ追加

    # N^2 なので改善したい
    # Intersections = [(divP1における番号, 座標, 使用したかどうか), ...]
    P1Intersections = Array{Tuple{Int64, Vector{typeof(P1.vertexes[1])}, Bool}, 1}()
    for i in 1:P1.n
        for j in 1:P2.n
            if i < P1.n && j < P2.n
                is_inter, intersection = segment_intersect(P1.vertexes[:, i], P1.vertexes[:, i+1], 
                                                        P2.vertexes[:, j], P2.vertexes[:, j+1])
            elseif i < P1.n && j == P2.n
                is_inter, intersection = segment_intersect(P1.vertexes[:, i], P1.vertexes[:, i+1],
                                                        P2.vertexes[:, end], P2.vertexes[:, 1])
            elseif i == P1.n && j < P2.n
                is_inter, intersection = segment_intersect(P1.vertexes[:, end], P1.vertexes[:, 1],
                                                        P2.vertexes[:, j], P2.vertexes[:, j+1])
            else
                is_inter, intersection = segment_intersect(P1.vertexes[:, end], P1.vertexes[:, 1],
                                                        P2.vertexes[:, end], P2.vertexes[:, 1])
            end
            if is_inter
                # TODO; P2 のインデックスが正しくないので要修正
                push!(P1Intersections, (i + length(P1Intersections) + 1, intersection, false))
            end
        end
    end

    # Intersections = [(divP2における番号, 座標, 使用したかどうか), ...]
    P2Intersections = Array{Tuple{Int64, Vector{typeof(P1.vertexes[1])}, Bool}, 1}()
    for i in 1:P2.n
        for j in 1:P1.n
            if i < P2.n && j < P1.n
                is_inter, intersection = segment_intersect(P2.vertexes[:, i], P2.vertexes[:, i+1], 
                                                        P1.vertexes[:, j], P1.vertexes[:, j+1])
            elseif i < P2.n && j == P1.n
                is_inter, intersection = segment_intersect(P2.vertexes[:, i], P2.vertexes[:, i+1],
                                                        P1.vertexes[:, end], P1.vertexes[:, 1])
            elseif i == P2.n && j < P1.n
                is_inter, intersection = segment_intersect(P2.vertexes[:, end], P2.vertexes[:, 1],
                                                        P1.vertexes[:, j], P1.vertexes[:, j+1])
            else
                is_inter, intersection = segment_intersect(P2.vertexes[:, end], P2.vertexes[:, 1],
                                                        P1.vertexes[:, end], P1.vertexes[:, 1])
            end
            if is_inter
                push!(P2Intersections, (i + length(P2Intersections) + 1, intersection, false))
            end
        end
    end
    num_inter = length(P1Intersections)

    divP1TodivP2 = Dict{Int64, Int64}()
    divP2TodivP1 = Dict{Int64, Int64}()
    for inter1 in P1Intersections
        i2 = findfirst(x -> prod(x[2] .- 1e-2 .< inter1[2] .< x[2] .+ 1e-2), P2Intersections)
        divP2TodivP1[i2] = inter1[1]
        divP1TodivP2[inter1[1]] = i2
    end

    # 交点が追加された頂点の集合
    nV1 = P1.vertexes[:, 1:P1Intersections[1][1] - 1]
    nV2 = P2.vertexes[:, 1:P2Intersections[1][1] - 1]

    for i in 1:num_inter-1
        # index1 = Intersections[i][1] - 1
        # index2 = Intersections[i][2] - 1
        nV1 = hcat(nV1, P1Intersections[i][2])
        nV2 = hcat(nV2, P2Intersections[i][2])

        nV1 = hcat(nV1, P1.vertexes[:, P1Intersections[i][1]-i+1:P1Intersections[i+1][1] - i - 1])
        nV2 = hcat(nV2, P2.vertexes[:, P2Intersections[i][1]-i+1:P2Intersections[i+1][1] - i - 1])
    end

    nV1 = hcat(nV1, P1Intersections[end][2])
    nV2 = hcat(nV2, P2Intersections[end][2])

    # nV1 = hcat(nV1, P1.vertexes[:, P1Intersections[end][1]-num_inter+1:end])
    # nV2 = hcat(nV2, P2.vertexes[:, P2Intersections[end][1]-num_inter+1:end])

    nV1 = P1Intersections[end][1]-num_inter < P1.n ? hcat(nV1, P1.vertexes[:, P1Intersections[end][1]-num_inter+1:end]) : nV1
    nV2 = P2Intersections[end][1]-num_inter < P2.n ? hcat(nV2, P2.vertexes[:, P2Intersections[end][1]-num_inter+1:end]) : nV2

    # P1 の１つ目の頂点が P2 の内点か調べる
    counter = 0
    for i in 1:P2.n-1
        counter += is_crossing_XrayFromP1FirstPoint(P1.vertexes[:, 1], P2, i, i+1)
    end
    counter += is_crossing_XrayFromP1FirstPoint(P1.vertexes[:, 1], P2, P2.n, 1)

    # poly1 の１つ目の頂点から巡回しながら交点が来たらそこから poly2 を巡回する
    # (共通部分を構成していたら poly1 と poly2 を行ったり来たりして、そうでなかったそうしない)
    # TODO; 共通部分が連結でない場合があるので、連結成分を個別に出力できるようにする
    IntersectionsPoly = num_inter > 0 ? Array{typeof(P1), 1}() : false
    V = Array{typeof(P1.vertexes[1]), 1}()
 
    divP1 = Polygon(nV1)
    divP2 = Polygon(nV2)

    markP1 = falses(divP1.n)
    markP2 = falses(divP2.n)

    now_initial = true
    while (IntersectionsPoly) != false
        global walker, initWalker

        # 初期化
        # while の外だとエラー出まくり侍なのでここでやる
        if now_initial
            walker = (1, 1, counter % 2 == 1)
            initWalker = walker
            now_initial = false
        end

        # ここからループ
        # 現在いる点を記録する
        if walker[3]
            if isempty(V)
                V = walker[1] == 1 ? divP1.vertexes[:, walker[2]] : divP2.vertexes[:, walker[2]]
            else
                V = walker[1] == 1 ? hcat(V, divP1.vertexes[:, walker[2]]) : hcat(V, divP2.vertexes[:, walker[2]])
            end
        end

        # 次の頂点に移動
        walker = next_walker(walker, divP1, divP2)

        # 交点にいたら図形を変える
        if walker[3]
            if walker[1] == 1
                if walker[2] ∈ P1Intersections[:][1]
                    i1 = findfirst(x -> x == walker[2], P1Intersections[:][1])
                    i2 = findfirst(x -> x == divP1TodivP2[walker[2]], P2Intersections[:][1])
                    P1Intersections[i1][3] = true
                    P2Intersections[i2][3] = true
                    walker = (2, divP1TodivP2[walker[2]], true)
                end
            else
                if walker[2] ∈ P2Intersections[:][1]
                    i1 = findfirst(x -> x == divP2TodivP1[walker[2]], P1Intersections[:][1])
                    i2 = findfirst(x -> x == walker[2], P2Intersections[:][1])
                    P1Intersections[i1][3] = true
                    P2Intersections[i2][3] = true
                    walker = (2, divP2TodivP1[walker[2]], true)
                end
            end
        else
            if ((walker[1] == 1 && walker[2] ∈ P1Intersections[:][1]) || 
                (walker[1] == 2 && walker[2] ∈ P2Intersections[:][1]))
                walker = (walker[1], walker[2], true)
            end
        end
        
        # 現在地点と作成した図形の初点が一致するか？
        # 一致するなら図形が完成している
        if walker[1:2] == initWalker[1:2]
            push!(IntersectionsPoly, Polygon(V))

            # 未探索の交点はまだあるか？
            # ないなら、未探索の共通部分はなくループを抜ける
            is_uncomplete = P1Intersections[1][3]
            for p1i in P1Intersections[2:end]
                is_uncomplete *= p1i[3]
            end
            if is_uncomplete
                ini = findfirst(!(P1Intersections[:][3]))
                walker = (1, ini, true)
                initWalker = walker
            else
                break
            end
        end
    end

    return IntersectionsPoly
end

P1 = Polygon([0, 1, 2, 0], [1, 1, 2, 2])
P2 = Polygon([0.5, 1.5, 1], [0, 0, 3])

iPs = intersect(P1, P2);

# 面積計算はグリーンの定理で
function area(Poly::Polygon)

    return S
end

function test_segment_intersection(; xlim=(-5, 5), ylim=(-5, 5))
    X = (xlim[2] - xlim[1]) * rand(Float64, 4) .+ xlim[1]
    Y = (ylim[2] - ylim[1]) * rand(Float64, 4) .+ ylim[1]
    XY = hcat(X, Y)
    is_intersect, p = segment_intersect(XY[1, :], XY[2, :], XY[3, :], XY[4, :])

    plt = plot(X[1:2], Y[1:2], m = :circle, xlim = xlim, ylim=ylim, size=(300, 300))
    plot!(plt, X[3:4], Y[3:4], m = :circle)
    if is_intersect
        scatter!(plt, [p[1]], [p[2]])
    end
    return plt
end

function test_is_innerPoint(; xlim=(-5, 5), ylim=(-5, 5))
    n = rand(3:20)
    X = (xlim[2] - xlim[1]) * rand(Float64, n) .+ xlim[1]
    Y = (ylim[2] - ylim[1]) * rand(Float64, n) .+ ylim[1]
    Poly = Polygon(X, Y)

    p = [0.0, 0.0]
    p[1] = (xlim[2] - xlim[1]) * rand(Float64) .+ xlim[1]
    p[2] = (ylim[2] - ylim[1]) * rand(Float64) .+ ylim[1]

    counter = 0
    for i in 1:n-1
        counter += is_crossing_XrayFromP1FirstPoint(p, Poly, i, i+1)
    end
    counter += is_crossing_XrayFromP1FirstPoint(p, Poly, n, 1)

    println(counter)

    push!(X, X[1]); push!(Y, Y[1])
    plt = plot(X, Y, m = :circle, xlim = xlim, ylim=ylim, size=(500, 500))
    scatter!(plt, [p[1]], [p[2]])

    return plt
end
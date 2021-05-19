using Plots
gr() # 図が早くて綺麗に書けるおまじない

mutable struct Polygon{T}
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
    return nothing
end

# Polygon を θ だけ回転させる
function rotate!(P::Polygon{T}, θ::T) where T
    R = [cos(θ) -sin(θ); sin(θ) cos(θ)]
    
    P.vertexes.-= P.center
    P.vertexes = R * P.vertexes
    P.vertexes .+= P.center
    return nothing
end

# 図形を描画する
function display(Polygons::Polygon...; center=false, vertex=false)
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

function ×(v::Vector{T}, u::Vector{T}) where T
    return v[1] * u[2] - v[2] * u[1]
end

# 戻り値: (一点で交わっているか, 追加すべきか, 交点座標)
function segment_intersect(A::Vector, B::Vector, C::Vector, D::Vector)
    s1 = (B - A) × (C - A) # AB×AC
    t1 = (B - A) × (D - A)# AB×AD

    s2 = (D - C) × (A - C) # CD×CA
    t2 = (D - C) × (B - C) # CD×CB

    p = [-1000.0, -1000.0]

    if s1*t1 < 0 && s2*t2 < 0
        AB = B - A
        DC = C - D
        M = [AB[1] DC[1]
             AB[2] DC[2]]
        λ = inv(M) * (C - A)
        p = A + λ[1] * AB
        return true, true, p
    elseif s1*t1 == 0 && s2*t2 < 0
        if s1 == 0
            p = C
        else
            p = D
        end
        return true, true, p
    elseif s1*t1 < 0 && s2*t2 == 0
        if s2 == 0
            p = A
        else
            p = B
        end
        return true, false, p
    # elseif s1*t1 == 0 && s2*t2 == 0
    #     t = (D-A)[1]/(D-C)[1]
    #     if 0 <= t <= 1
    #         p = A
    #     else
    #         p = B
    #     end
    end
    return false, false, p
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
    # Intersections = [(divP1における番号, 座標, 使用したかどうか, 交点を追加するかどうか), ...]
    P1Intersections = Array{Tuple{Int64, Vector{typeof(P1.vertexes[1])}, Bool, Bool}, 1}()
    total_add = 0
    for i in 1:P1.n
        tmp_x = Array{typeof(P1.vertexes[1]), 1}()
        tmp_y = Array{typeof(P1.vertexes[1]), 1}()
        tmp_a = Array{Bool, 1}()

        for j in 1:P2.n
            if i < P1.n && j < P2.n
                is_inter, is_add, intersection = segment_intersect(P1.vertexes[:, i], P1.vertexes[:, i+1], 
                                                        P2.vertexes[:, j], P2.vertexes[:, j+1])
            elseif i < P1.n && j == P2.n
                is_inter, is_add, intersection = segment_intersect(P1.vertexes[:, i], P1.vertexes[:, i+1],
                                                        P2.vertexes[:, end], P2.vertexes[:, 1])
            elseif i == P1.n && j < P2.n
                is_inter, is_add, intersection = segment_intersect(P1.vertexes[:, end], P1.vertexes[:, 1],
                                                        P2.vertexes[:, j], P2.vertexes[:, j+1])
            else
                is_inter, is_add, intersection = segment_intersect(P1.vertexes[:, end], P1.vertexes[:, 1],
                                                        P2.vertexes[:, end], P2.vertexes[:, 1])
            end
            if is_inter
                push!(tmp_x, intersection[1])
                push!(tmp_y, intersection[2])
                push!(tmp_a, is_add)
            end
        end

        r = sortperm((tmp_x .- P1.vertexes[1, i]).^2 + (tmp_y .- P1.vertexes[2, i]).^2)
        for j in 1:length(r)
            # ここがダメらしい
            if j > 1 || tmp_a[r[j]]
                push!(P1Intersections, (i + total_add + 1, [tmp_x[r[j]], tmp_y[r[j]]], false, tmp_a[r[j]]))
                total_add += tmp_a[j]
            else
                push!(P1Intersections, (i + total_add, [tmp_x[r[j]], tmp_y[r[j]]], false, tmp_a[r[j]]))
            end
        end
    end

    # Intersections = [(divP2における番号, 座標, 使用したかどうか, 交点を追加するかどうか), ...]
    P2Intersections = Array{Tuple{Int64, Vector{typeof(P1.vertexes[1])}, Bool, Bool}, 1}()
    total_add = 0
    for i in 1:P2.n
        tmp_x = Array{typeof(P1.vertexes[1]), 1}()
        tmp_y = Array{typeof(P1.vertexes[1]), 1}()
        tmp_a = Array{Bool, 1}()

        for j in 1:P1.n
            if i < P2.n && j < P1.n
                is_inter, is_add, intersection = segment_intersect(P2.vertexes[:, i], P2.vertexes[:, i+1], 
                                                        P1.vertexes[:, j], P1.vertexes[:, j+1])
            elseif i < P2.n && j == P1.n
                is_inter, is_add, intersection = segment_intersect(P2.vertexes[:, i], P2.vertexes[:, i+1],
                                                        P1.vertexes[:, end], P1.vertexes[:, 1])
            elseif i == P2.n && j < P1.n
                is_inter, is_add, intersection = segment_intersect(P2.vertexes[:, end], P2.vertexes[:, 1],
                                                        P1.vertexes[:, j], P1.vertexes[:, j+1])
            else
                is_inter, is_add, intersection = segment_intersect(P2.vertexes[:, end], P2.vertexes[:, 1],
                                                        P1.vertexes[:, end], P1.vertexes[:, 1])
            end
            if is_inter
                push!(tmp_x, intersection[1])
                push!(tmp_y, intersection[2])
                push!(tmp_a, is_add)
            end
        end

        r = sortperm((tmp_x .- P2.vertexes[1, i]).^2 + (tmp_y .- P2.vertexes[2, i]).^2)
        for j in 1:length(r)
            # ここがダメらしい
            if j > 1 || tmp_a[r[j]]
                push!(P2Intersections, (i + total_add + 1, [tmp_x[r[j]], tmp_y[r[j]]], false, tmp_a[r[j]]))
                total_add += tmp_a[r[j]]
            else
                push!(P2Intersections, (i + total_add, [tmp_x[r[j]], tmp_y[r[j]]], false, tmp_a[r[j]]))
            end
        end
    end
    num_inter = length(P1Intersections)

    divP1TodivP2 = Dict{Int64, Int64}()
    divP2TodivP1 = Dict{Int64, Int64}()
    for inter1 in P1Intersections
        i2 = findfirst(x -> prod(abs.(inter1[2] .- x[2]) .<  .+ 1e-3), P2Intersections)
        # println(findall(x -> prod(x[2] .- 1e-2 .< inter1[2] .< x[2] .+ 1e-2), P2Intersections))
        #  i2 = findfirst(x -> x[2] == inter1[2], P2Intersections)
        divP2TodivP1[P2Intersections[i2][1]] = inter1[1]
        divP1TodivP2[inter1[1]] = P2Intersections[i2][1]
    end

    # 交点が追加された頂点の集合
    nV1 = num_inter > 0 ? P1.vertexes[:, 1:P1Intersections[1][1] - 1] : P1.vertexes
    nV2 = num_inter > 0 ? P2.vertexes[:, 1:P2Intersections[1][1] - 1] : P2.vertexes

    if num_inter > 0
        counter1 = 0
        counter2 = 0
        for i in 1:num_inter-1
            # index1 = Intersections[i][1] - 1
            # index2 = Intersections[i][2] - 1
            if P1Intersections[i][4]
                counter1 += 1
                nV1 = hcat(nV1, P1Intersections[i][2])
                nV1 = hcat(nV1, P1.vertexes[:, P1Intersections[i][1] - counter1 + 1:P1Intersections[i+1][1] - counter1 - 1])
            else
                nV1 = hcat(nV1, P1Intersections[i][2])
                nV1 = hcat(nV1, P1.vertexes[:, P1Intersections[i][1] - counter1 + 1:P1Intersections[i+1][1] - counter1 - 1])
            end
            if P2Intersections[i][4]
                counter2 += 1
                nV2 = hcat(nV2, P2Intersections[i][2])
                nV2 = hcat(nV2, P2.vertexes[:, P2Intersections[i][1] - counter2 + 1:P2Intersections[i+1][1] - counter2 - 1])
            else
                nV2 = hcat(nV2, P2Intersections[i][2])
                nV2 = hcat(nV2, P2.vertexes[:, P2Intersections[i][1] - counter2 + 1:P2Intersections[i+1][1] - counter2 - 1])
            end
        end

        nV1 = hcat(nV1, P1Intersections[end][2])
        nV2 = hcat(nV2, P2Intersections[end][2])

        counter1 += P1Intersections[end][4] ? 1 : 0
        counter2 += P2Intersections[end][4] ? 1 : 0

        # nV1 = hcat(nV1, P1.vertexes[:, P1Intersections[end][1]-num_inter+1:end])
        # nV2 = hcat(nV2, P2.vertexes[:, P2Intersections[end][1]-num_inter+1:end])

        nV1 = P1Intersections[end][1]-counter1 < P1.n ? hcat(nV1, P1.vertexes[:, P1Intersections[end][1]-counter1+1:end]) : nV1
        nV2 = P2Intersections[end][1]-counter2 < P2.n ? hcat(nV2, P2.vertexes[:, P2Intersections[end][1]-counter2+1:end]) : nV2
    end

    # P1 の１つ目の頂点が P2 の内点か調べる
    counter1 = 0
    for i in 1:P2.n-1
        counter1 += is_crossing_XrayFromP1FirstPoint(P1.vertexes[:, 1], P2, i, i+1)
    end
    counter1 += is_crossing_XrayFromP1FirstPoint(P1.vertexes[:, 1], P2, P2.n, 1)

    # P2 の１つ目の頂点が P1 の内点か調べる
    counter2 = 0
    for i in 1:P1.n-1
        counter2 += is_crossing_XrayFromP1FirstPoint(P2.vertexes[:, 1], P1, i, i+1)
    end
    counter2 += is_crossing_XrayFromP1FirstPoint(P2.vertexes[:, 1], P1, P1.n, 1)

    # poly1 の１つ目の頂点から巡回しながら交点が来たらそこから poly2 を巡回する
    # (共通部分を構成していたら poly1 と poly2 を行ったり来たりして、そうでなかったそうしない)
    # TODO; 共通部分が連結でない場合があるので、連結成分を個別に出力できるようにする
    IntersectionsPoly = num_inter > 0 ? Array{typeof(P1), 1}() : false

    # P1 が P2 に完全に含まれている場合
    if (IntersectionsPoly) == false && counter1 % 2 == 1
        return [P1]
    # P2 が P1 に完全に含まれている場合
    elseif (IntersectionsPoly) == false && counter2 % 2 == 1
        return [P2]
    end

    V = Array{typeof(P1.vertexes[1]), 1}()
 
    divP1 = Polygon(nV1)
    divP2 = Polygon(nV2)

    now_initial = true

    # viw = [10000.0, 10000.0]
    while (IntersectionsPoly) != false
        global walker, initWalker# , viw

        # 初期化
        # while の外だとエラー出まくり侍なのでここでやる
        if now_initial
            if P1Intersections[1][1] == 1
                walker = (1, 1, true)
            else
                walker = (1, 1, counter1 % 2 == 1)
            end
            if walker[3]
                initWalker = walker
            end
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
            # P1 について
            if walker[1] == 1
                # ここがおかしい
                on_intersection = false
                for P1I in P1Intersections
                    if walker[2] == P1I[1]
                        on_intersection = true
                        break
                    end
                end
                if on_intersection
                    i1 = findfirst(x -> x[1] == walker[2], P1Intersections)
                    i2 = findfirst(x -> x[1] == divP1TodivP2[walker[2]], P2Intersections)
                    P1Intersections[i1] = (P1Intersections[i1][1], P1Intersections[i1][2], true, P1Intersections[i1][4])
                    P2Intersections[i2] = (P2Intersections[i2][1], P2Intersections[i2][2], true, P2Intersections[i2][4])
                    walker = P1Intersections[i1][4] ? (2, P2Intersections[i2][1], true) : walker
                end
            # P2 について
            else
                # ここもおかしい
                on_intersection = false
                for P2I in P2Intersections
                    if walker[2] == P2I[1]
                        on_intersection = true
                        break
                    end
                end
                if on_intersection
                    i1 = findfirst(x -> x[1] == divP2TodivP1[walker[2]], P1Intersections)
                    i2 = findfirst(x -> x[1] == walker[2], P2Intersections)
                    P1Intersections[i1] = (P1Intersections[i1][1], P1Intersections[i1][2], true, P1Intersections[i1][4])
                    P2Intersections[i2] = (P2Intersections[i2][1], P2Intersections[i2][2], true, P2Intersections[i2][4])
                    walker = P2Intersections[i2][4] ? (1, P1Intersections[i1][1], true) : walker
                end
            end 
        else
            # ここもおかしくなりそう
            on_intersection = false
            if walker[1] == 1
                for P1I in P1Intersections
                    if walker[2] == P1I[1]
                        on_intersection = true # ここも P1I[4] にしないといけないか？
                        break
                    end
                end
            elseif walker[1] ==2 
                for P2I in P2Intersections
                    if walker[2] == P2I[1]
                        on_intersection = true # ここも P2I[4] にしないといけないか？
                        break
                    end
                end
            end

            if on_intersection
                walker = (walker[1], walker[2], true)
                initWalker = walker
                # viw =  walker[1] == 1 ? divP1.vertexes[:, walker[2]] : divP2.vertexes[:, walker[2]]
                now_initial = true
            end
        end

        # 現在地点と作成した図形の初点が一致するか？
        # 一致するなら図形が完成している
        is_comeback = false
        if walker[1] == 1 && walker[3]
            is_comeback = walker[2] == initWalker[2]
        elseif walker[1] == 2 && walker[3]
            is_comeback = walker[2] == divP1TodivP2[initWalker[2]]
        end
        if  is_comeback && !(now_initial)
            if size(V)[2] >= 3
                push!(IntersectionsPoly, Polygon(V))
            end

            # 未探索の交点はまだあるか？
            # ないなら、未探索の共通部分はなくループを抜ける
            is_complete = P1Intersections[1][3]
            for p1i in P1Intersections[2:end]
                is_complete *= p1i[3]
            end
            if is_complete
                break
            else
                ini1 = findfirst(x -> x[3] == false, P1Intersections)
                # ini2 = findfirst(x -> x[1] == divP1TodivP2[ini1], P2Intersections)
                # P1Intersections[ini1] = (P1Intersections[ini1][1], P1Intersections[ini1][2], false)
                # P2Intersections[ini2] = (P2Intersections[ini2][1], P2Intersections[ini2][2], false)
                walker = (1, P1Intersections[ini1][1]-1, false)
                # viw = [10000.0, 10000.0]
                V = Array{typeof(P1.vertexes[1]), 1}()
            end
        end
        now_initial = false
    end

    return isempty(IntersectionsPoly) ? false : IntersectionsPoly
end

# テスト用の図形
P1 = Polygon([0, 2, 1], [ 0, 0, 2])
P2 = Polygon([0, 2, 1], [0, 0, 2])

# P1 = Polygon( [])

move!(P2, 0.5, 2.0)

iPs = intersect(P1, P2)

# P1 = Polygon([1.5, 1.0, 0.5, 1.0], [1.0, 2.0, 1.0, 0.0])

# 面積計算はグリーンの定理で
function area(Poly::Polygon)
    S = zero(Poly.vertexes[1])
    for i in 1:Poly.n-1
        S += Poly.vertexes[1, i] * Poly.vertexes[2, i+1] - Poly.vertexes[1, i+1] * Poly.vertexes[2, i]
    end
    S += Poly.vertexes[1, end] * Poly.vertexes[2, 1] - Poly.vertexes[1, 1] * Poly.vertexes[2, end]

    return abs(S)/2
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
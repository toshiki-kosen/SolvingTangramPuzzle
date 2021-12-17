using LinearAlgebra: sqrt
using Printf: Threads
using Plots: Threads
using ProgressMeter
using Printf
include("PolygonBase.jl")
include("ESBase.jl")
include("silhouette&pieces.jl")

# X[1つめのピースのx座標, 1つめのピースのy座標, １つめの回転角Θ, 2つめのx座標, ...]
function loss_poly(X::Array, loss_args::Array{Float64, 1})
    global pices, silhouette

    # silhouetteとpieceの共通部分の面積を求める
    tmp_p = move(pieces[1], X[1], X[2])
    rotate!(tmp_p, X[3] * 2π)
    unioned = MYPolygon2LibGEOS(tmp_p)
    sum_pieces = LibGEOS.area(MYPolygon2LibGEOS(tmp_p))

    for i in 2:length(pieces)
        tmp_p = move(pieces[i], X[3i-2], X[3i-1])
        rotate!(tmp_p, X[3i] * 2π)
        unioned = LibGEOS.union(unioned, MYPolygon2LibGEOS(tmp_p))
        sum_pieces += MYPolygon2LibGEOS(tmp_p) |> LibGEOS.area
    end
    
    A = SymmetricDifference(unioned, silhouette) / LibGEOS.area(MYPolygon2LibGEOS(silhouette))

    # 頂点の誤差成分 (シルエットベース 緩い)
    Δv = 0.0
    for v_s in 1:silhouette.n
        Δmin = 100.0
        for j in 1:length(pieces)
            tmp_p = move(pieces[j], X[3j - 2], X[3j - 1])
            rotate!(tmp_p, X[3j] * 2π)
            for v_p in 1:tmp_p.n
                Δmin = min(Δmin, sum((silhouette.vertexes[:,v_s] - tmp_p.vertexes[:, v_p]).^2))
            end
        end
        Δv += Δmin / silhouette.n
        # Δv = max(Δv, Δmin)
    end
    
    # 頂点の誤差成分 (ピースベース 厳しい)
    """
    Δv = 0.0
    for i in 1:length(pieces)
        p = move(pieces[i], X[3i - 2], X[3i - 1])
        rotate!(p, X[3i] * 2π)
        Δmax = 0.0
        # 頂点ごとに誤差を見る
        for v_p in 1:p.n
            Δmin = 100.0
            # シルエットについて見る
            for v_s in 1:silhouette.n
                Δmin = min(Δmin, sum((silhouette.vertexes[:,v_s] - p.vertexes[:, v_p]).^2))
            end
            # その他のピースについて見る
            for j in 1:length(pieces)
                if j == i continue end
                q = move(pieces[j], X[3j - 2], X[3j - 1])
                rotate!(p, X[3j] * 2π)
                for v_q in 1:q.n
                    Δmin = min(Δmin, sum((q.vertexes[:,v_q] - p.vertexes[:, v_p]).^2))
                end
            end
            Δmax = max(Δmax, Δmin)
        end
        # Δv += 2tanh(2Δmin) / (π*silhouette.n)
        # Δv = max(Δv, Δmin)
        Δv += Δmax
    end
    """

    # 辺の誤差成分
    E = 0.0
    for i in 1:length(pieces)
        P = move(pieces[i], X[3i-2], X[3i-1])
        rotate!(P, X[3i] * 2π)

        cosθ = 0.0

        # シルエットの各頂点と距離が最小な点を探す
        min_d = 100.0
        for p in 0:P.n-1, s in 0:silhouette.n-1
            d = sum((P.vertexes[:, p + 1] - silhouette.vertexes[:, s + 1]).^2)
            # 最小の大きさの角度を探す
            if d < min_d
                eP = P.vertexes[:, mod(p, P.n) + 1] - P.vertexes[:, mod(p-1, P.n) + 1]
                eS = silhouette.vertexes[:, mod(s, silhouette.n) + 1] - silhouette.vertexes[:, mod(s-1, silhouette.n) + 1]
                cosθ1 = abs(eP' * eS / (norm(eP) * norm(eS) ))

                eP = P.vertexes[:, mod(p, P.n) + 1] - P.vertexes[:, mod(p-1, P.n) + 1]
                eS = silhouette.vertexes[:, mod(s+1, silhouette.n) + 1] - silhouette.vertexes[:, mod(s, silhouette.n) + 1]
                cosθ2 = abs(eP' * eS / (norm(eP) * norm(eS) ))

                eP = P.vertexes[:, mod(p+1, P.n) + 1] - P.vertexes[:, mod(p, P.n) + 1]
                eS = silhouette.vertexes[:, mod(s, silhouette.n) + 1] - silhouette.vertexes[:, mod(s-1, silhouette.n) + 1]
                cosθ3 = abs(eP' * eS / (norm(eP) * norm(eS) ))

                eP = P.vertexes[:, mod(p+1, P.n) + 1] - P.vertexes[:, mod(p, P.n) + 1]
                eS = silhouette.vertexes[:, mod(s+1, silhouette.n) + 1] - silhouette.vertexes[:, mod(s, silhouette.n) + 1]
                cosθ4 = abs(eP' * eS / (norm(eP) * norm(eS) ))

                cosθ = max(cosθ1, cosθ2, cosθ3, cosθ4)
            end
        end
        E += √(1 - cosθ^2) / length(pieces)
    end

    return loss_args[1] * A + loss_args[2] * Δv + loss_args[3] * E 
end

function loss_poly2(X::Array, loss_args::Array{Float64, 1})
    global pices, silhouette

    # silhouetteとpieceの共通部分の面積を求める
    tmp_p = move(pieces[1], X[1], X[2])
    rotate!(tmp_p, X[3] * 2π)
    unioned = MYPolygon2LibGEOS(tmp_p)
    sum_pieces = LibGEOS.area(MYPolygon2LibGEOS(tmp_p))

    for i in 2:length(pieces)
        tmp_p = move(pieces[i], X[3i-2], X[3i-1])
        rotate!(tmp_p, X[3i] * 2π)
        unioned = LibGEOS.union(unioned, MYPolygon2LibGEOS(tmp_p))
        sum_pieces += MYPolygon2LibGEOS(tmp_p) |> LibGEOS.area
    end
    F = LibGEOS.area(LibGEOS.intersection(unioned, silhouette)) / silhouette.area

    # 辺の誤差成分
    E = 0.0
    for P in pieces
        minimum_sinθ = 1
        for i in 1:P.n
            p1 = P.vertexes[:, i]
            p2 = i < n ? P.vertexes[:, i+1] : P.vertexes[:, 1]

            p21 = p2 - p1
            norm_p21 = norm(p21)

            # シルエットを見る
            for j in 1:silhouette.n
                s1 = silhouette.vertexes[:, j]
                s2 = j < silhouette.n ? silhouette.vertexes[:, j+1] : silhouette.vertexes[:, 1]
                sinθ = cross(p21, s2-s1)/(norm_p21*norm(s2-s1))
                minimum_sinθ = min(sinθ, minimum_sinθ)
            end
            # 他のピースを見る
            for Q in pieces
                if P == Q
                    continue
                end
                for j in 1:Q.n
                    q1 = Q.vertexes[j, :]
                    q2 = j < Q.n ? Q.vertexes[:, j+1] : Q.vertexes[:, 1]
                    sinθ = cross(p21, q2-q1)/(norm_p21*norm(q2-q1))
                    minimum_sinθ = min(sinθ, minimum_sinθ)
                end
            end
        end
        E += minimum_sinθ / length(pieces) * 2P.area
    end

    # 頂点の誤差成分
    V = 0.0
    for P in pieces
        v = 0.0
        for i in 1:P.n
            p = P.vertexes[:, i]
            
            min_d = 100.0
            # シルエットについて見る
            for j in 1:silhouette.n
                ps = norm(p - silhouette.vertexes[:, j])
                min_d = min(min_d, ps)
            end

            # 他のピースについて見る
            for Q in pieces
                if Q == P
                    continue
                end
                for j in 1:Q.n
                    pq = norm(p - Q.vertexes.[:, j])
                    min_d = min(min_d, pq)
                end
            end
            v += min_d / P.n
        end
        V += 2 * P.area * v / length(pieces)
    end

    return loss_args[1] * F + loss_args[2] * E + loss_args[3] * V 
end

# 面積ベースのピース批評
# シルエットから飛び出たやつしか扱えない
function reviewer_1(pieces::Array{MYPolygon{Float64}, 1}, X, silhouette::MYPolygon, threshold)
    # silhouetteと各pieceの共通部分の面積を求める
    score = zeros(Bool, length(pieces))

    for i in 1:length(pieces)
        tmp_p = move(pieces[i], X[3i-2], X[3i-1])
        rotate!(tmp_p, X[3i] * 2π)
        GEOS_p = MYPolygon2LibGEOS(tmp_p)
        score[i] = threshold < intersectionArea(tmp_p, silhouette) / LibGEOS.area(GEOS_p)
    end
    return score
end

# 頂点ベースの批評
# おそらく、あらゆるケースを扱えるが計算量がエグそう
function reviewer_2(pieces::Array{MYPolygon{Float64}, 1}, X, silhouette::MYPolygon, threshold::T) where T <: Number
    # 各ピースの各頂点の中でもっともシルエットや他のピースの辺から離れているものを探す
    score = trues(length(pieces))

    for i in 1:length(pieces)
        optied_piece = move(pieces[i], X[3i-2], X[3i-1])
        rotate!(optied_piece, X[3i] * 2π)
        for p in 1:optied_piece.n
            P = optied_piece.vertexes[:, p]
            vertex_check = false

            # シルエットについて
            for s in 1:silhouette.n
                A = silhouette.vertexes[:, s]
                if s < silhouette.n
                    B = silhouette.vertexes[:, s+1]
                else
                    B = silhouette.vertexes[:, 1]
                end

                if distance_seg2point(A, B, P) < threshold
                    vertex_check = true
                    break
                end
            end

            # 他のピースについて
            for j in 1:length(pieces)
                if vertex_check
                    break   
                elseif j == i
                    continue
                end
                optied_qiece = move(pieces[j], X[3j-2], X[3j-1])
                rotate!(optied_qiece, X[3j] * 2π)
                for q in 1:pieces[j].n
                    A = optied_qiece.vertexes[:, q]
                    if q < pieces[j].n
                        B = optied_qiece.vertexes[:, q+1]
                    else
                        B = optied_qiece.vertexes[:, 1]
                    end

                    if distance_seg2point(A, B, P) < threshold
                        vertex_check = true
                        break
                    end
                end
                if vertex_check
                    break
                end
            end

            # 一つでもダメな頂点があったらそのピースは終わり
            if vertex_check
            else
                score[i] = false
                break
            end
        end
    end
    return score
end

silhouette = antena
pieces = [tri_l, tri_s, tri_s, square_s]

loss_args = [100.0, 92.0, 64.0]

# 初期化
max_gen = 128
sample_num = 256
λ = 18
recycling_durability = 4
recycle_threshold = 0.1

# 本体
rm("outputs", force=true, recursive=true)
mkdir("outputs")

best_fitnesses = Array{Float64, 1}()
p = Progress(sample_num)
for t in 1:sample_num
    cmaes = init_CMAES(zeros(3 * length(pieces)), 0.8, λ)
    fitnesses = zeros(cmaes.dim)
    X = zeros((cmaes.dim, cmaes.λ))
    Y = zeros((0, cmaes.λ))
    rng = MersenneTwister(t)

    for recycle in 1:recycling_durability
        for gen in 1:max_gen
            # 個体生成
            X = samplePopulation(cmaes, rng=rng)

            # 回転角のパラメータを -π ~ π までに正規化
            for j in 1:size(X)[2]
                for i in 3:3:size(X)[1]
                    X[i, j] = mod(X[i, j], 1.0)
                end
            end

            # 評価値 
            loss(x) = loss_poly(x, loss_args)
            fitnesses = get_fitness(vcat(X, Y), loss)

            # 更新 
            update!(cmaes, X, fitnesses, gen)
        end
        # 正解っぽいピースと不正解っぽいピースに分ける
        Z = vcat(X, Y)
        best_arg = argmin(fitnesses)
        score = reviewer_2(pieces, Z[:, best_arg], silhouette, recycle_threshold)
        
        new_pieces = Array{MYPolygon{Float64}, 1}()
        y = zeros(Float64, 0)
        x = zeros(Float64, 0)
        bad_num = 0
        for i in 1:length(pieces)
            if score[i]
                push!(new_pieces, pieces[i])
                push!(y, Z[3i-2:3i, best_arg]...)
            else
                pushfirst!(new_pieces, pieces[i])
                pushfirst!(x, Z[3i-2:3i, best_arg]...)
                bad_num += 1
            end
        end
        if bad_num == 0 break end
        
        global pieces = deepcopy(new_pieces)
        Y = repeat(y, 1, λ)
        X = repeat(x, 1, λ)
        cmaes = init_CMAES(zeros(3 * bad_num), 1.0, λ)
    end
    push!(best_fitnesses, maximum(fitnesses))
    
    # 表示用
    best_arg = argmin(fitnesses)
    best_pieces = Array{MYPolygon, 1}()
    Z = vcat(X, Y)
    for i in 1:length(pieces)
        todisplay = move(pieces[i], Z[3i-2, best_arg], Z[3i-1, best_arg])
        rotate!(todisplay, Z[3i, best_arg] * 2π)
        push!(best_pieces, todisplay)
    end
    display(silhouette, best_pieces...)
    savefig(@sprintf "outputs\\%02.0f_%d_results" -best_fitnesses[end] t)

    next!(p)
end

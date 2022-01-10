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
    tmp = LibGEOS.intersection(unioned, MYPolygon2LibGEOS(silhouette))
    A = LibGEOS.area(tmp) / LibGEOS.area(MYPolygon2LibGEOS(silhouette))

    # 各2pieceの共通部分を求める
    inter = 1.0 - LibGEOS.area(unioned) / sum_pieces
    outer = 1.0 - LibGEOS.area(tmp)/sum_pieces

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
        E += cosθ / length(pieces)
    end

    return -loss_args[1] * A + loss_args[2] * Δv + loss_args[3] * outer + loss_args[4] * (E - 1)
end

# X[1つめのピースのx座標, 1つめのピースのy座標, １つめの回転角Θ, 2つめのx座標, ...]
function loss_poly_light(X::Array)
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
    tmp = LibGEOS.intersection(unioned, MYPolygon2LibGEOS(silhouette))
    A = LibGEOS.area(tmp) / LibGEOS.area(MYPolygon2LibGEOS(silhouette))

    return -A * 100
end

function loss_poly2(X::Array, loss_args::Array{Float64, 1})
    global pices, silhouette

    all_area = 0
    for P in pieces
        all_area += P.area
    end

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
    F = LibGEOS.area(LibGEOS.intersection(unioned, MYPolygon2LibGEOS(silhouette))) / silhouette.area

    # 辺の誤差成分
    E = 0.0
    for P in pieces
        sinθ_ave = 1
        for i in 1:P.n
            minimum_e = 1
            p1 = P.vertexes[:, i]
            p2 = i < P.n ? P.vertexes[:, i+1] : P.vertexes[:, 1]

            p21 = p2 - p1
            norm_p21 = norm(p21)

            # シルエットを見る
            for j in 1:silhouette.n
                s1 = silhouette.vertexes[:, j]
                s2 = j < silhouette.n ? silhouette.vertexes[:, j+1] : silhouette.vertexes[:, 1]
                sinθ = cross(p21, s2-s1)/(norm_p21*norm(s2-s1))
                minimum_e = min(sinθ, minimum_e)
            end
            # 他のピースを見る
            for Q in pieces
                if P == Q
                    continue
                end
                for j in 1:Q.n
                    q1 = Q.vertexes[:, j]
                    q2 = j < Q.n ? Q.vertexes[:, j+1] : Q.vertexes[:, 1]
                    sinθ = cross(p21, q2-q1)/(norm_p21*norm(q2-q1))
                    minimum_e = min(abs(sinθ), minimum_e)
                end
            end
           sinθ_ave += minimum_e / P.n
        end
        E += P.area *sinθ_ave / (length(pieces) * all_area)
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
                    pq = norm(p - Q.vertexes[:, j])
                    min_d = min(min_d, pq)
                end
            end
            v += min_d / P.n
        end
        V += P.area * v / (length(pieces) * all_area)
    end

    return -loss_args[1] * F + loss_args[2] * E + loss_args[3] * V 
end

# silhouette = hexagon_m
# pieces = [tri_s, tri_s, parallelogram, tri_m]
# silhouetteとpiecesとパズル名のトリプル
# p_Taiikusuwari は図形の頂点指定ミスでエラー起こる
# puzzles = [p_house2, p_house4, p_arrow2, p_trapezoid2_1, p_trapezoid2_2, p_tri_above_square, p_fish2, p_small_stand, p_triforce3, p_rect_s, p_hurt3, p_hexagon_m1, p_hexagon_m2, p_hexagon5, p_square_3, p_square_4, p_square_7, p_butterfly_s, p_butterfly, p_pesudo_butterfly, p_butterfly4, p_note, p_diamond, p_antena, p_pencile4, p_crown, p_snake, p_faucet, p_Taiikusuwari, p_crab_scissor, p_candy7, p_flame7, p_hurt7, p_katana]
loss_args = [100.0, 260.0, 0.0, 60.0]

# 初期化
max_gen = 128
sample_num = 128

for puzzle in puzzles
    global puzzle_name, silhouette, pieces = puzzle

    fpath = "outputs\\outputs_$puzzle_name"
    rm(fpath, force=true, recursive=true)
    mkdir(fpath)
 
    best_fitnesses = Array{Float64, 1}() 
    p = Progress(sample_num, desc="$puzzle_name: ")
    for t in 1:sample_num
        # local cmaes, rng
        cmaes = init_CMAES(zeros(3 * length(pieces)), 0.8, 18)
        fitnesses = zeros(cmaes.dim)
        X = zeros((cmaes.dim, cmaes.λ))
        rng = MersenneTwister(t)

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
            fitnesses = get_fitness(X, loss)

            # 更新
            update!(cmaes, X, fitnesses, gen)
        end
        push!(best_fitnesses, maximum(fitnesses))

        # 表示用
        best_arg = argmin(fitnesses)
        best_pieces = Array{MYPolygon, 1}()
        for i in 1:length(pieces)
            todisplay = move(pieces[i], X[3i-2, best_arg], X[3i-1, best_arg])
            rotate!(todisplay, X[3i, best_arg] * 2π)
            push!(best_pieces, todisplay)
        end
        display(silhouette, best_pieces...)
        savefig(@sprintf "outputs\\outputs_%s\\%02.0f_%d_results" puzzle_name -best_fitnesses[end] t)

        next!(p)
    end
end
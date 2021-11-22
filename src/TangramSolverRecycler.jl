using Plots: string
using Base: Float16, Float64
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
        rotate!(P, X[3i])

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

    return -loss_args[1] * A + loss_args[2] * Δv + loss_args[3] * outer + loss_args[4] * E
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

silhouette = small_stand
pieces = [tri_l, tri_s, tri_s]

# 初期化
λ = 18
cmaes = init_CMAES(zeros(3 * length(pieces)), 1.0, λ)
seed = MersenneTwister()
max_gen = 96

recycling_durability = 4
recycle_threshold = 0.1
bad_num = length(pieces)
X = zeros(Float64, 3length(pieces), λ)
Y = zeros(Float64, 0, λ)

println(seed)

fitnesses_ave = Array{Float64, 1}()
fitnesses_max = Array{Float64, 1}()

loss_args = [100.0, 92.0, 0.0, 64.0]

for recycle in 1:recycling_durability
    # best_argのグローバル化
    best_arg = 0
    # Z = vcat(X, Y)
    anim = @animate for gen in 1:max_gen
        # 個体生成 
        global X = samplePopulation(cmaes, rng=seed)

        # 回転角のパラメータを -π ~ π までに正規化
        for j in 1:size(X)[2]
            for i in 3:3:size(X)[1]
                X[i, j] = mod(X[i, j], 1.0)
            end
        end

        # 評価値
        Z = vcat(X, Y)
        loss(x) = loss_poly(x, loss_args)
        fitnesses = get_fitness(Z, loss)

        # 表示用
        best_arg = argmin(fitnesses)
        best_pieces = Array{MYPolygon, 1}()
        for i in 1:length(pieces)
            todisplay = move(pieces[i], Z[3i-2, best_arg], Z[3i-1, best_arg])
            rotate!(todisplay, Z[3i, best_arg] * 2π)
            push!(best_pieces, todisplay)
        end
        display(silhouette, best_pieces...)

        # 評価用
        push!(fitnesses_ave, sum(fitnesses)/size(Z)[2])
        push!(fitnesses_max, minimum(fitnesses))

        # 更新
        update!(cmaes, X, fitnesses, gen)
    end
    gif(anim, "best_pieces_$(recycle).gif", fps=10)

    # 正解っぽいピースと不正解っぽいピースに分ける
    Z = vcat(X, Y)
    score = reviewer_2(pieces, Z[:, best_arg], silhouette, recycle_threshold)
    println(score)
    
    new_pieces = Array{MYPolygon{Float64}, 1}()
    global Y = zeros(Float64, 0)
    global bad_num = 0
    for i in 1:length(pieces)
        if score[i]
            push!(new_pieces, pieces[i])
            push!(Y, Z[3i-2:3i, best_arg]...)
        else
            pushfirst!(new_pieces, pieces[i])
            global bad_num += 1
        end
    end
    if bad_num == 0 break end

    global pieces = deepcopy(new_pieces)
    global Y = repeat(Y, 1, λ)
    global cmaes = init_CMAES(zeros(3 * bad_num), 1.0, λ)
end

# pl = plot(1:max_gen, -fitnesses_max, lab="Average", xaxis="generation", yaxis="fitness [%]", legend = false, w=2.5)
# plot!(pl, 1:max_gen, -fitnesses_max, line=:dash, lab="Maximum", ylim=(minimum(-fitnesses_ave), 110), yticks = 100:-20:minimum(-fitnesses_ave))
# plot!(pl, xtickfontsize=13, ytickfontsize=13, xguidefontsize=13, yguidefontsize=13)

# savefig(pl, "fitnesses.png")
# Base.display(pl)

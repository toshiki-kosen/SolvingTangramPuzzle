using LinearAlgebra: sqrt
using Printf: Threads
using Plots: Threads
using ProgressMeter
using Printf
include("PolygonBase.jl")
include("ESBase.jl")
include("silhouette&pieces.jl")

# X[1つめのピースのx座標, 1つめのピースのy座標, １つめの回転角Θ, 2つめのx座標, ...]
function loss_poly(X::Array; Ks::Float64 = 1.0, Ki::Float64 = 0.0, Kv::Float64 = 0.8)
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

    # 頂点の誤差成分 (シルエットベース)
    """
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
        # Δv += 2tanh(2Δmin) / (π*silhouette.n)
        Δv = max(Δv, Δmin)
    end
    """
    # 頂点の誤差成分 (ピースベース)
    Δv = 0.0
    for i in 1:length(pieces)
        p = move(pieces[i], X[3i - 2], X[3i - 1])
        rotate!(p, X[3i] * 2π)
        Δmin = 100.0
        # シルエットについて見る
        for v_s in 1:silhouette.n
            for v_p in 1:p.n
                Δmin = min(Δmin, sum((silhouette.vertexes[:,v_s] - p.vertexes[:, v_p]).^2))
            end
        end
        # その他のピースについて見る
        for j in 1:length(pieces)
            if j == i continue end
            q = move(pieces[j], X[3j - 2], X[3j - 1])
            rotate!(p, X[3j] * 2π)
            for v_q in 1:q.n
                for v_p in 1:p.n
                    Δmin = min(Δmin, sum((q.vertexes[:,v_q] - p.vertexes[:, v_p]).^2))
                end
            end
        end
        # Δv += 2tanh(2Δmin) / (π*silhouette.n)
        Δv = max(Δv, Δmin)
    end

    return -100Ks * A + 100Kv * Δv + 100Ki * inter
end

silhouette = square_s
pieces = [tri_s, tri_s]

# 初期化
max_gen = 96
sample_num = 256

# 最終結果を保存するか否か
save_results = false
if save_results
    rm("outputs", force=true, recursive=true)
    mkdir("outputs")

    best_fitnesses = Array{Float64, 1}()
    correct_time = 0
    p = Progress(sample_num)
    for t in 1:sample_num
        global best_fitnesses, correct_time
        local cmaes, rng
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
            loss(x) = loss_poly(x,Ks = 1.0, Ki=0.5, Kv=0.0)
            fitnesses = get_fitness(X, loss)

            # 更新
            update!(cmaes, X, fitnesses, gen)
        end
        push!(best_fitnesses, maximum(fitnesses))

        if -best_fitnesses[end] > 95
            correct_time += 1
        end
        
        # 表示用
        best_arg = argmin(fitnesses)
        best_pieces = Array{MYPolygon, 1}()
        for i in 1:length(pieces)
            todisplay = move(pieces[i], X[3i-2, best_arg], X[3i-1, best_arg])
            rotate!(todisplay, X[3i, best_arg] * 2π)
            push!(best_pieces, todisplay)
        end
        display(silhouette, best_pieces...)
        savefig(@sprintf "outputs\\%02.0f_%d_results" -best_fitnesses[end] t)

        next!(p)
    end
else
    best_fitnesses = Array{Float64, 1}()
    correct_time = Threads.Atomic{Int}(0)
    p = Progress(sample_num÷Threads.nthreads())
    Threads.@threads for t in 1:sample_num
        global best_fitnesses, correct_time
        local cmaes
        cmaes = init_CMAES(zeros(3 * length(pieces)), 1.0, 0)
        fitnesses = zeros(cmaes.dim)
        X = zeros((cmaes.dim, cmaes.λ))
        rng = MersenneTwister(t)

        for gen in 1:max_gen
            # 個体生成
            X = samplePopulation(cmaes, rng=rng)

            # 回転角のパラメータを -π ~ π までに正規化
            for j in 1:size(X)[2]
                for i in 3:3:size(X)[1]
                    X[i, j] %= 1
                end
            end

            # 評価値
            loss(x) = loss_poly(x,Ks = 1.0, Ki=0.0 , Kv=0.0)
            fitnesses = get_fitness(X, loss)

            # 更新
            update!(cmaes, X, fitnesses, gen)
        end
        push!(best_fitnesses, maximum(fitnesses))

        if -best_fitnesses[end] > 94
            Threads.atomic_add!(correct_time, 1)
        end
        
        if Threads.threadid() == 1
            next!(p)
        end
    end
end

pl = histogram(-best_fitnesses, bins = 40:5:100, xlabel="Fitness", ylabel="The number of answers", lab=false)

println("correct prob: $(correct_time[])/$(sample_num) = $(100correct_time[]/sample_num)%")

savefig(pl, "fitness_distribusion_$(sample_num)samples_$(correct_time[])-corrects.png")
Base.display(pl)
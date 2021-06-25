using Base: Float64
using Printf
include("PolygonBase.jl")
include("ESBase.jl")
include("silhouette&pieces.jl")

levy(X::Array) = sin(3π*X[1])^2 + (X[1] - 1)^2 * (1 + sin(3π*X[2])^2) + (X[2] - 1)^2 * (1 + sin(2π * X[2])^2)
levy(x, y) = levy([x, y])

rastrigin(x, y) = 10 * 2 + x^2 - 10cos(2π*x) + y^2 - 10cos(2π*y)
rastrigin(X::Array) = rastrigin(X[1], X[2])

eggholder(x, y) = -(100*y + 47)*sin(sqrt(abs(50x + 100y + 47))) - 100x*sin(sqrt(abs(100x - 100y - 47)))
eggholder(X::Array) = eggholder(X[1], X[2])

# X[1つめのピースのx座標, 1つめのピースのy座標, １つめの回転角Θ, 2つめのx座標, ...]
function loss_poly(X::Array; Ki::Float64 = 0.0, Kv::Float64 = 0.8)
    global pices, silhouette

    # silhouetteとpieceの共通部分の面積を求める
    tmp_p = move(pieces[1], X[1], X[2])
    rotate!(tmp_p, X[3])
    unioned = MYPolygon2LibGEOS(tmp_p)
    sum_pieces = LibGEOS.area(MYPolygon2LibGEOS(tmp_p))

    for i in 2:length(pieces)
        tmp_p = move(pieces[i], X[3i-2], X[3i-1])
        rotate!(tmp_p, X[3i])
        unioned = LibGEOS.union(unioned, MYPolygon2LibGEOS(tmp_p))
        sum_pieces += MYPolygon2LibGEOS(tmp_p) |> LibGEOS.area
    end
    tmp = LibGEOS.intersection(unioned, MYPolygon2LibGEOS(silhouette))
    A = LibGEOS.area(tmp) / LibGEOS.area(MYPolygon2LibGEOS(silhouette))

    # 各2pieceの共通部分を求める
    inter = 1.0 - LibGEOS.area(unioned) / sum_pieces

    # 頂点の誤差成分 
    Δv = 0.0
    for i in 1:silhouette.n
        Δmin = 100.0
        for i in 1:length(pieces)
            tmp_p = move(pieces[i], X[3i - 2], X[3i - 1])
            rotate!(tmp_p, X[3i])
            for j in 1:tmp_p.n
                Δmin = min(Δmin, sum((silhouette.vertexes[:,i] - tmp_p.vertexes[:, j]).^2))
            end
        end
        Δv = max(Δv, Δmin)
    end

    return -100A + 100Kv * Δv + 100Ki * inter
end

silhouette = house2
pieces = [square_s, tri_m]

# 初期化
max_gen = 96
sample_num = 256

# 最終結果を保存するか否か
save_results = true
if save_results
    rm("outputs", force=true, recursive=true)
    mkdir("outputs")
end

best_fitnesses = Array{Float64, 1}()
correct_time = 0
for t in 1:sample_num
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
                X[i, j] %= 2π - π
            end
        end

        # 評価値
        loss(x) = loss_poly(x, Ki=0.0, Kv=0.8)
        fitnesses = get_fitness(X, loss)

        # 更新
        update!(cmaes, X, fitnesses, gen)
    end
    push!(best_fitnesses, maximum(fitnesses))
    correct_time += -best_fitnesses[end] > 97 ? 1 : 0
    if save_results
        # 表示用
        best_arg = argmin(fitnesses)
        best_pieces = Array{MYPolygon, 1}()
        for i in 1:length(pieces)
            todisplay = move(pieces[i], X[3i-2, best_arg], X[3i-1, best_arg])
            rotate!(todisplay, X[3i, best_arg])
            push!(best_pieces, todisplay)
        end
        display(silhouette, best_pieces...)
        savefig(@sprintf "outputs\\%02.0f_%d_results" -best_fitnesses[end] t)
    end
end

pl = histogram(-best_fitnesses, bins = 40:5:100)

println("correct prob: $(correct_time)/$(sample_num) = $(100correct_time/sample_num)%")

savefig(pl, "fitness_distribusion_$(sample_num)samples.png")
Base.display(pl)


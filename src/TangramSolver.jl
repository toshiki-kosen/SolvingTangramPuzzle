using Plots: length
using Base: Float64
include("PolygonBase.jl")
include("ESBase.jl")

levy(X::Array) = sin(3π*X[1])^2 + (X[1] - 1)^2 * (1 + sin(3π*X[2])^2) + (X[2] - 1)^2 * (1 + sin(2π * X[2])^2)
levy(x, y) = levy([x, y])

rastrigin(x, y) = 10 * 2 + x^2 - 10cos(2π*x) + y^2 - 10cos(2π*y)
rastrigin(X::Array) = rastrigin(X[1], X[2])

eggholder(x, y) = -(100*y + 47)*sin(sqrt(abs(50x + 100y + 47))) - 100x*sin(sqrt(abs(100x - 100y - 47)))
eggholder(X::Array) = eggholder(X[1], X[2])

silhouette = MYPolygon([0, 1, 1, 0], [0, 0, 1, 1])
pieces = [MYPolygon([-1, 0, 0, -1], [0.5, 0.5, 1.5, 1.5])]
rotate!(pieces[1], 0.5)

# X[1つめのピースのx座標, 1つめのピースのy座標, １つめの回転角Θ, 2つめのx座標, ...]
function loss_poly(X::Array)
    global pices, silhouette

    tmp_p = move(pieces[1], X[1], X[2])
    rotate!(tmp_p, X[3])
    unioned = MYPolygon2LibGEOS(tmp_p)

    for i in 2:length(pieces)
        tmp_p = move(pieces[i], X[3i-2], X[3i-1])
        rotate!(tmp_p, X[3i])
        unioned = LibGEOS.union(unioned, MYPolygon2LibGEOS(tmp_p))
    end

    tmp = LibGEOS.intersection(unioned, MYPolygon2LibGEOS(silhouette))

    return -100 * LibGEOS.area(tmp)/ LibGEOS.area(MYPolygon2LibGEOS(silhouette))
end

cmaes = init_CMAES([0.0, 0.0, 0.0], 2.0, 32)

fitnesses_ave = Array{Float64, 1}()
fitnesses_max = Array{Float64, 1}()

max_gen = 100

anim = @animate for gen in 1:max_gen
    global fitnesses_ave, fitnesses_max
    local X, fitnesses
    # 個体生成
    X = samplePopulation(cmaes)

    # 回転角のパラメータを -π ~ π までに正規化
    for j in 1:size(X)[2]
        for i in 3:3:size(X)[1]
            X[i, j] %= 2π - π
        end
    end

    # 評価値
    fitnesses = get_fitness(X, loss_poly)

    # 表示用
    best_arg = argmin(fitnesses)
    best_pieces = Array{MYPolygon, 1}()
    for i in 1:length(pieces)
        todisplay = move(pieces[i], X[3i-2, best_arg], X[3i-1, best_arg])
        rotate!(todisplay, X[3i, best_arg])
        push!(best_pieces, todisplay)
    end
    display(silhouette, best_pieces...)

    # 評価用
    push!(fitnesses_ave, sum(fitnesses)/size(X)[2])
    push!(fitnesses_max, minimum(fitnesses))

    # 更新
    update!(cmaes, X, fitnesses, gen)
end

gif(anim, "best_pieces.gif", fps=10)

plot(1:max_gen, -fitnesses_ave)
plot!(1:max_gen, -fitnesses_max, line=:dash)
plot!(1:max_gen, 100 .+ zeros(max_gen), line=:dot)


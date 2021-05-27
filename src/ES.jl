include("PolygonBase.jl")
using Random

f(X) = 3sin(X[1]/3.0) + 2cos(X[2]/2.0)

τ(n) = 1/√(10n)

# A の各列を引数としてFの最も高い1~n位のインデックスを求める
function best_n(A::Array, n::Int64, F::Function)
    score = zeros(size(A)[2])
    for i in 1:size(A)[2]
        score[i] = F(A[:, i])
    end
    
    return sortperm(score, rev=true)[1:n]
end

function test_commaES(;μ::Int64, ρ::Int64, λ::Int64, F::Function, t=100)
    if λ > μ
        throw(DomainError(λ, "λ must be larger than μ"))
    end
    # 親群の初期化 (x, y, σ) x 10
    P = 2 * rand(Float64, (3, μ)) .- 1
    
    for n in 1:t
        # 子供を作る
        perm = best_n(P[1:2, :], ρ, F)
        Pρ = copy(P[:, perm])

        # 組み替え & 突然変異
        C = zeros(3, λ)
        ave_σ = sum(Pρ[3, :])/ρ
        C[3, :] = ave_σ * exp.(τ(n) .* randn(λ))

        ave_y = (sum(Pρ[1, :])./ρ, sum(Pρ[2, :])./ρ)
        for i in 1:λ
            C[1:2, i] = ave_y .+ C[3, i] .* randn(2)
        end

        # 次の親を選ぶ
        perm = best_n(C[1:2, :], μ, F)
        P = copy(C[:, perm])
    end

    return P
end

function test_plusES(;μ::Int64, ρ::Int64, λ::Int64, F::Function, t=100)
    # 親群の初期化 (x, y, σ) x 10
    P = 100 * randn(Float64, (3, μ))
    
    for n in 1:t
        # 子供を作る
        perm = best_n(P[1:2, :], ρ, F)
        Pρ = copy(P[:, perm])

        # 組み替え & 突然変異
        C = zeros(3, λ)
        ave_σ = sum(Pρ[3, :])/ρ
        C[3, :] = ave_σ * exp.(τ(n) .* randn(λ))

        ave_y = (sum(Pρ[1, :])./ρ, sum(Pρ[2, :])./ρ)
        for i in 1:λ
            C[1:2, i] = ave_y .+ C[3, i] .* randn(2)
        end

        # 次の親を選ぶ
        candidate = hcat(C, P)
        perm = best_n(candidate[1:2, :], μ, F)
        P = copy(candidate[:, perm])
    end

    return P
end
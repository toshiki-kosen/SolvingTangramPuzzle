using Base: UInt16, Real
using LinearAlgebra: Matrix
using Random, LinearAlgebra

f(X) = 3sin(X[1]/3.0) + 2cos(X[2]/2.0)

τ(n) = 1/√(10n)

mutable struct CMAES{T}
    dim::Int
    λ::Int
    μ::Int

    centroid::Array{T, 1}
    c_m::T

    weights::Array{T, 1}
    μ_eff::T

    σ::T
    p_σ::Array{T, 1}
    c_σ::T
    d_σ::T

    C::Array{T, 2}
    p_c::Array{T, 1}
    c_c::T
    c_1::T
    c_μ::T
end

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

# size(X) == (λ, dim)
function get_fitness(X::Matrix{Float64}, f::Function)
    ToReturn = zeros(Float64, size(X)[2])

    for i in 1:size(X)[2]
        ToReturn[i] = f(X[:, i])
    end

    return ToReturn
end

# CMAESパラメータの初期化
function init_CMAES(center::Array{T, 1}, σ::T, λ::Int) where T <: Real
    ## 初期化フェーズ
    # 次元数
    dim = length(center)

    # 世代の総個体数 λ と エリート数 μ
    λ = λ > 0 ? λ : round(Int, 4 + 3 * log(dim))
    μ = λ÷2

    # 正規分布の中心と学習率
    centroid = T.(center)
    c_m = one(T)

    # 順位重み係数
    weights = [log(0.5(λ+1)) - log(i) for i in 1:μ]
    weights = T.(weights ./ sum(weights))
    μ_eff = 1 / sum(weights.^2)

    # ステップサイズ
    p_σ = zeros(T, dim)
    c_σ = (μ_eff + 2) / (dim + μ_eff + 5)
    d_σ = 1 + 2 * max(0, sqrt((μ_eff - 1)/(dim + 1)) - 1) + c_σ

    # 共分散行列 進化パス p と rank-μ, rank-one更新の学習率c
    C = Matrix{T}(I, dim, dim)
    p_c = zeros(T, dim)
    c_c = (4 + μ_eff / dim) / (dim + 4 + 2 * μ_eff/dim)
    c_1 = 2 / ((dim + 1.3f0)^2 + μ_eff)
    c_μ = min(1 - c_1, 2 * (μ_eff - 2 + 1/μ_eff) / ((dim+2)^2 + μ_eff))

    ToReturn = CMAES(dim, λ, μ, centroid, c_m, weights, μ_eff, σ, p_σ, c_σ, d_σ, 
                C, p_c, c_c, c_1, c_μ)
    return ToReturn
end

# 個体の生成
function samplePopulation(self::CMAES{T}; rng=MersenneTwister(123)) where T <: Real
    ## 探索フェーズ
    # z ~ N(0, I) なる個体を λ 個生成
    Z = randn(rng, (self.dim, self.λ))

    # C を固有値分解
    Ei = eigen(self.C)
    diagD = sqrt(Diagonal(Ei.values))
    B = Ei.vectors
    BD = B * diagD

    # y ~ N(0, C)
    Y = BD * Z
    # X ~ N(μ, σC)
    X = self.centroid .+ self.σ * Y

    return X
end

# 個体および行列 C の更新
function update!(self::CMAES{T}, X, fitnesses, gen) where T <: Real
    ### 更新フェーズ
    """ 正規分布パラメータの更新
    X : 個体群, shape == (λ, μ)
    fitnesses : 適合度
    gen : 現代の世代数
    """

    ## 1. Selection and recombination
    old_centroid = self.centroid
    old_σ = self.σ

    # fitnesses が上位 μ までのインデックスを抽出
    elite_indices = sortperm(fitnesses)[1:self.μ]

    X_elite = X[:, elite_indices]
    Y_elite = (X_elite .- old_centroid) ./ old_σ

    X_w = (X_elite * self.weights)
    Y_w = (Y_elite * self.weights)

    # 正規分布中心の更新
    self.centroid = (1 - self.c_m) * old_centroid .+ self.c_m * X_w

    ## 2. Step-size control
    Ei = eigen(self.C)
    diagD = sqrt(Diagonal(Ei.values))
    B = Ei.vectors
    inv_diagD = diagD^(-1)

    # Note. B*Z == C_ * Y
    C_ = B * inv_diagD * B'

    new_p_σ = (1 - self.c_σ) * self.p_σ
    new_p_σ += sqrt(self.c_σ * (2 - self.c_σ) * self.μ_eff) * C_ * Y_w
    self.p_σ = new_p_σ

    E_normal = sqrt(self.dim) * (1 - 1/(4*self.dim) + 1/(21 * self.dim ^ 2)) # 定数パラメータ
    self.σ = self.σ * exp((self.c_σ / self.d_σ) * (sqrt(sum(self.p_σ .^ 2))/E_normal - 1))

    # 3. Covariance matrix adaptation (CMA)
    # Note. h_σ(heaviside 関数) はステップサイズ σ が大きいときに C の更新を中断させるのに使う
    left = sqrt(sum(self.p_σ.^2)) / sqrt(1 - (1 - self.c_σ) ^ (2 * (gen + 1)))
    right = (1.4 + 2 / (self.dim + 1)) * E_normal
    hσ = left < right ? 1 : 0
    d_hσ = (1 - hσ) * self.c_c * (2 - self.c_c)

    # p_c の更新
    new_p_c = (1 - self.c_c) * self.p_c
    new_p_c += hσ * sqrt(self.c_c * (2 - self.c_c) * self.μ_eff) * Y_w
    self.p_c = new_p_c

    # C の更新
    new_C = (1 + self.c_1 * d_hσ - self.c_1 - self.c_μ) * self.C
    new_C += self.c_1 * [i * j for i in self.p_c, j in self.p_c]

    # 愚直な実装(スマートな実装はdeapのcma.pyを参照)
    wyy = zeros(T, (self.dim, self.dim))
    for i in 1:self.μ
        y_i = Y_elite[i]
        wyy .+= self.weights[i] * [y1 * y2 for y1 in y_i, y2 in y_i]
    end
    new_C += self.c_μ * wyy

    self.C = new_C
end
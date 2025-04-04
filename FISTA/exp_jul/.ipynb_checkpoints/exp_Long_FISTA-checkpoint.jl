module Long

export LongModel, 
    　 create_long_model, bond, breakingdown, Z_LP, long_df, solve, logit, long_armijo, Lyapunov

using SparseArrays
using LinearAlgebra
using Optim
using CSV
using DataFrames
using Plots

struct LongModel
    prm::Parameter
    algprm::AlgoParameter
    short::ShortProblem
end


function create_long_model(prm::Parameter, algprm::AlgoParameter, short::ShortProblem)
    return LongModel(prm, algprm, short)
end


function bond(long::LongModel, m::Vector{Float64}, n::Matrix{Float64})
    """企業・家計分布を一つの結合リストにする関数

    Parameters
    ----------
    m: Vector{Float64} (K,)
        企業分布
    n: Matrix{Float64} (K, K)
        家計分布

    Returns
    ----------
    mn: Vector{Float64} (K + K * K,)
        企業・家計分布の結合リスト
    """
    n_list = vec(n)  # 行列を1次元のベクトルに変換
    mn = vcat(m, n_list)
    return mn
end


function breakingdown(long::LongModel, mn::Vector{Float64})
    """結合リストを企業・家計分布に分解する関数

    Parameters
    ----------
    mn: Vector{Float64} (K + K * K,)
        企業・家計分布の結合リスト

    Returns
    ----------
    m: Vector{Float64} (K,)
        企業分布
    n: Matrix{Float64} (K, K)
        家計分布
    """
    K = long.prm.K
    m = mn[1:K]
    n_list = mn[K+1:end]
    n = reshape(n_list, K, K)
    return m, n
end

function Z_LP(long::LongModel, mn::Vector{Float64}, RW::Vector{Float64})
    """目的関数を計算する関数

    Parameters
    ----------
    mn: Vector{Float64} (K + K * K,)
        企業・家計分布の結合リスト

    Returns
    ----------
    F: Float64
        目的関数値
    """
    m, n = breakingdown(long, mn)
    Short = Z_SD(long.short, RW, m, n)

    F_value = -Short \
              - 0.5 * m' * long.prm.D * m \
              + (1 / long.prm.theta_firm) * (dot(m, log.(m ./ long.prm.M))) \
              + (1 / long.prm.theta_house) * sum(n .* log.(n ./ long.prm.N))

    return F_value
end

function long_df(long::LongModel, mn::Vector{Float64}, RW::Vector{Float64})
    """目的関数の勾配を計算する関数

    Parameters
    ----------
    mn: Vector{Float64} (K + K * K,)
        企業・家計分布の結合リスト
    RW: Vector{Float64} (2 * K,)
        価格変数の結合リスト

    Returns
    ----------
    dF: Vector{Float64} (K + K * K,)
        目的関数の勾配
    """
    m, n = breakingdown(long, mn)
    dF_m = -pi(long.short, RW[1:long.prm.K], RW[long.prm.K+1:end], m) + (1 / long.prm.theta_firm) * (log.(m ./ long.prm.M) .+ 1)
    dF_n = -v(long.short, RW[1:long.prm.K], RW[long.prm.K+1:end]) .+ (1 / long.prm.theta_house) * (log.(n ./ long.prm.N) .+ 1)

    dF = bond(long, dF_m, dF_n)
    return dF
end

function solve(long::LongModel, m0::Vector{Float64}, n0::Matrix{Float64}, RW_ini::Vector{Float64}, err_short::Float64, err_long::Float64, obj_corr::Float64, long_itr::Int)
    """長期均衡を解く関数

    Parameters
    ----------
    m0: Vector{Float64} (K,)
        企業の初期分布
    n0: Matrix{Float64} (K, K)
        家計の初期分布
    RW_ini: Vector{Float64} (2 * K,)
        価格変数の初期値
    err_short: Float64
        収束判定の閾値
    err_long: Float64
        長期収束判定の閾値
    obj_corr: Float64
        目的関数の正解値
    long_itr: Int
        最大反復回数

    Returns
    -------
    m_true: Vector{Float64} (K,)
        最終的な企業分布
    n_true: Matrix{Float64} (K, K)
        最終的な家計分布
    RW: Vector{Float64} (2 * K,)
        最終的な価格変数
    """
    mn0 = bond(long, m0, n0)
    mn_before = mn0
    obj_before = 0.0

    long_iteration = []
    obj_list = []
    obj_rel_list = []

    for k in 1:long_itr
        m_before, n_before = breakingdown(long, mn_before)

        R, W, iteration, short_obj_rel = short_solve(long.short, RW_ini, m_before, n_before, err_short, 100000, 0.0)
        RW = vcat(R, W)

        mn_d = logit(long.short, R, W, m_before)

        obj = Z_LP(long, mn_before, RW)
        push!(long_iteration, k)
        push!(obj_list, obj)

        obj_rel = abs(obj - obj_corr) / abs(obj_corr)
        push!(obj_rel_list, obj_rel)

        if maximum(abs.((mn_before .- mn_d) ./ mn_before)) < err_long
            break
        end

        alpha = 0.05
        mn = (1 - alpha) * mn_before .+ alpha * mn_d

        mn_before = mn
        obj_before = obj
    end

    m_true, n_true = breakingdown(long, mn)
    return m_true, n_true, RW, long_iteration, obj_list, obj_rel_list
end

function logit(long::LongModel, R::Vector{Float64}, W::Vector{Float64}, m::Vector{Float64})
    """長期均衡の探索方向を計算する関数

    Parameters
    ----------
    R: Vector{Float64} (K,)
        地代
    W: Vector{Float64} (K,)
        賃金
    m: Vector{Float64} (K,)
        企業分布

    Returns
    -------
    mn_d: Vector{Float64} (K + K * K,)
        企業分布の探索方向と家計分布の探索方向の結合リスト
    """
    pi = long.short.pi(R, W, m)
    v = long.short.v(R, W)

    m_d = long.prm.M .* exp.(long.prm.theta_firm * pi) ./ sum(exp.(long.prm.theta_firm * pi))

    n_d = long.prm.N .* exp.(long.prm.theta_house * v) ./ sum(exp.(long.prm.theta_house * v))

    mn_d = bond(long, m_d, n_d)

    return mn_d
end

function long_armijo(long::LongModel, Z_LP, dZ, mn::Vector{Float64}, mn_d::Vector{Float64}, RW::Vector{Float64}; c_1=0.5, beta=0.6)
    """Armijo条件を満たすステップサイズを探す関数

    Parameters
    ----------
    Z_LP: Function
        目的関数
    dZ: Function
        勾配関数
    mn: Vector{Float64} (K + K * K,)
        現在の企業・家計分布
    mn_d: Vector{Float64} (K + K * K,)
        探索方向
    RW: Vector{Float64} (2 * K,)
        価格変数の結合リスト
    c_1: Float64
        Armijo条件のパラメータ
    beta: Float64
        ステップサイズ縮小係数

    Returns
    -------
    t: Float64
        適切なステップサイズ
    """
    t = 1.0

    while true
        println("Armijo")
        if Z_LP(mn .+ t .* mn_d, RW) < Z_LP(mn, RW) + c_1 * t * (dZ(mn, RW)' * mn_d)
            break
        end
        t *= beta
    end

    return t
end

function Lyapunov(long::LongModel, m::Vector{Float64}, n::Matrix{Float64}, RW::Vector{Float64})
    """Lyapunov関数を計算する関数

    Parameters
    ----------
    m: Vector{Float64} (K,)
        企業分布
    n: Matrix{Float64} (K, K)
        家計分布
    RW: Vector{Float64} (2 * K,)
        価格変数の結合リスト

    Returns
    -------
    G: Float64
        Lyapunov関数の値
    """
    R = RW[1:long.prm.K]
    W = RW[long.prm.K+1:end]

    pi = long.short.pi(R, W, m)
    v = long.short.v(R, W)

    G = sum(pi .* m) + sum(v .* n) \
        - (long.prm.M * log(sum(exp.(long.prm.theta_firm * pi))) / long.prm.theta_firm) \
        - (long.prm.N * log(sum(exp.(long.prm.theta_house * v))) / long.prm.theta_house) \
        - (sum(m .* log.(m ./ long.prm.M)) / long.prm.theta_firm) \
        - (sum(n .* log.(n ./ long.prm.N)) / long.prm.theta_house)

    return G
end
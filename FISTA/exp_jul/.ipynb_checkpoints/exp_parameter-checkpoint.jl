module Parameter

export Parameter, Algo_Parameter, create_parameter, create_algo_parameter

using LinearAlgebra
using SparseArrays

# モデルのパラメータを保持するクラス
struct Parameter
    Col::Int
    K::Int
    tau::Float64
    t::Float64
    Scaling::Float64
    S_total::Float64
    S_bar::Float64
    theta_firm::Float64
    theta_house::Float64
    E::Float64
    RW_proj::Float64
    alpha_1::Float64
    alpha_2::Float64
    beta_1::Float64
    beta_2::Float64
    alter_T_num::Float64
    M::Float64
    N::Float64
    distance_matrix::Matrix{Float64}
    T::Matrix{Float64}
    D::Matrix{Float64}

    function create_parameter(
        Col::Int, distance_matrix::Matrix{Float64}, t::Float64, tau::Float64, Scaling::Float64, 
        S_total::Float64, S_bar::Float64, theta_firm::Float64, theta_house::Float64, 
        E::Float64, RW_proj::Float64, alter_T_num::Float64, alpha_1::Float64, alpha_2::Float64, beta_1::Float64, 
        beta_2::Float64, M::Float64, N::Float64
    )
        K = Col * Col
        T = max.(Scaling * t * alter_T_num, t * distance_matrix)
        D = exp.(-tau * distance_matrix)

        new(
            Col, K, tau, t, Scaling, S_total, S_bar, theta_firm, theta_house, 
            E, RW_proj, alpha_1, alpha_2, beta_1, beta_2, alter_T_num, M, N, 
            distance_matrix, T, D
        )
    end
end

# アルゴリズム中のパラメータを保持するクラス
struct Algo_Parameter
    L::Float64
    eta::Float64
    p_proj::::Float64

    function create_algo_parameter(L::Float64, eta::Float64, p_proj::Float64)
        new(L, eta, p_proj)
    end
end

end # module ShortAlgo
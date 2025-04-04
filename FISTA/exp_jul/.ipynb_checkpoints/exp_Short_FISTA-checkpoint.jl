using LinearAlgebra

export ShortModel,
       create_short_model, demand_sparse, 

struct ShortModel
    prm::Parameter
    algprm::AlgoParameter
    power_S_H_R::Float64
    power_S_H_W::Float64
    power_L_H_R::Float64
    power_L_H_W::Float64
    power_S_F_R::Float64
    power_S_F_W::Float64
    power_L_F_R::Float64
    power_L_F_W::Float64
    T_power::Float64
    alpha_S_R::Float64
    alpha_S_W::Float64
    alpha_L_R::Float64
    alpha_L_W::Float64
    beta_S_R::Float64
    beta_S_W::Float64
    beta_L_R::Float64
    beta_L_W::Float64
    T_power_v::Float64
    pi_beta::Float64
    Z_I::Array{Float64, 2}
    Z_J::Array{Float64, 2}
    one_K::Vector{Float64}
    one_2K::Vector{Float64}
end

function create_short_model(prm::Parameter, algprm::AlgoParameter)
    power_S_H_R = (1 - prm.alpha_2) / (1 - prm.alpha_1 - prm.alpha_2)
    power_S_H_W = prm.alpha_2 / (1 - prm.alpha_1 - prm.alpha_2)
    power_L_H_R = prm.alpha_1 / (1 - prm.alpha_1 - prm.alpha_2)
    power_L_H_W = (1 - prm.alpha_1) / (1 - prm.alpha_1 - prm.alpha_2)

    power_S_F_R = (1 - prm.beta_2) / (1 - prm.beta_1 - prm.beta_2)
    power_S_F_W = prm.beta_2 / (1 - prm.beta_1 - prm.beta_2)
    power_L_F_R = prm.beta_1 / (1 - prm.beta_1 - prm.beta_2)
    power_L_F_W = (1 - prm.beta_1) / (1 - prm.beta_1 - prm.beta_2)

    T_power = (1 / prm.T) ^ (1 / (1 - prm.alpha_1 - prm.alpha_2))
    alpha_S_R = prm.alpha_1 ^ power_S_H_R
    alpha_S_W = prm.alpha_2 ^ power_S_H_W
    alpha_L_R = prm.alpha_1 ^ power_L_H_R
    alpha_L_W = prm.alpha_2 ^ power_L_H_W

    beta_S_R = prm.beta_1 ^ power_S_F_R
    beta_S_W = prm.beta_2 ^ power_S_F_W
    beta_L_R = prm.beta_1 ^ power_L_F_R
    beta_L_W = prm.beta_2 ^ power_L_F_W

    T_power_v = (1 - prm.alpha_1 - prm.alpha_2) * T_power
    pi_beta = 1 - prm.beta_1 - prm.beta_2

    Z_I = zeros(Float64, prm.K, prm.K * prm.K)
    Z_J = zeros(Float64, prm.K, prm.K * prm.K)

    rows_I = collect(1:prm.K)
    cols_I = collect(0:(prm.K - 1)) * prm.K

    for i in rows_I
        for j in 1:prm.K
            Z_I[i, cols_I[i] + j] = 1.0
        end
    end

    for i in 1:prm.K
        for j in 1:prm.K
            Z_J[i, (i - 1) * prm.K + j] = 1.0
        end
    end

    one_K = ones(Float64, prm.K)
    one_2K = ones(Float64, 2 * prm.K)

    return ShortModel(prm, algprm, power_S_H_R, power_S_H_W, power_L_H_R, power_L_H_W,
                      power_S_F_R, power_S_F_W, power_L_F_R, power_L_F_W, T_power,
                      alpha_S_R, alpha_S_W, alpha_L_R, alpha_L_W,
                      beta_S_R, beta_S_W, beta_L_R, beta_L_W, T_power_v,
                      pi_beta, Z_I, Z_J, one_K, one_2K)
end

function demand_sparse(model::ShortModel, R::Vector{Float64}, W::Vector{Float64})
    S_H = model.T_power *
          (model.alpha_S_R * (R .^ -model.power_S_H_R)) * (model.alpha_S_W * (W .^ -model.power_S_H_W))

    L_H = model.T_power *
          (model.alpha_L_R * (R .^ -model.power_L_H_R)) * (model.alpha_L_W * (W .^ -model.power_L_H_W))

    S_F = model.beta_S_R * (R .^ -model.power_S_F_R) * model.beta_S_W * (W .^ -model.power_S_F_W)
    L_F = model.beta_L_R * (R .^ -model.power_L_F_R) * model.beta_L_W * (W .^ -model.power_L_F_W)

    return S_H, L_H, S_F, L_F
end

function pi(model::ShortModel, R::Vector{Float64}, W::Vector{Float64}, m::Vector{Float64})
    # 利潤関数を計算する関数

    pi = model.prm.D * m +
         (1 - model.prm.beta_1 - model.prm.beta_2) *
         (model.beta_L_R * (R .^ -model.power_L_F_R)) *
         (model.beta_S_W * (W .^ -model.power_S_F_W))

    return pi
end

function pi_noex(model::ShortModel, R::Vector{Float64}, W::Vector{Float64}, m::Vector{Float64})
    # 需要関数を設定する関数

    pi = model.pi_beta *
         (model.beta_L_R * (R .^ -model.power_L_F_R)) *
         (model.beta_S_W * (W .^ -model.power_S_F_W))

    return pi
end

function v(model::ShortModel, R::Vector{Float64}, W::Vector{Float64})
    # 間接効用関数を計算する関数

    v = model.prm.E * W .+
        model.T_power_v *
        (model.alpha_L_R * (R .^ -model.power_L_H_R)) * 
        (model.alpha_S_W * (W .^ -model.power_S_H_W))

    return v
end

function Z_SD(model::ShortModel, RW::Vector{Float64}, m::Vector{Float64}, n::Matrix{Float64})
    # 目的関数を計算する関数

    K = model.prm.K
    R = RW[1:K]
    W = RW[(K + 1):end]

    # 間接効用関数の計算
    v_value = calculate_v(model, R, W)
    
    # flatten() の代わりに ravel()
    F_value = sum(v_value .* n) +
              sum(calculate_pi_noex(model, R, W, m) .* m) +
              model.prm.S_bar * sum(R)

    return F_value
end

function short_dual_df(model::ShortModel, RW::Vector{Float64}, m::Vector{Float64}, n::Matrix{Float64})
    # 目的関数の勾配を計算する関数

    K = model.prm.K
    R = RW[1:K]
    W = RW[(K + 1):end]

    # 需要関数の計算
    S_H, L_H, S_F, L_F = demand_sparse(model, R, W)

    # 勾配の計算
    dR = model.prm.S_bar .- sum(S_H .* n, dims=2)[:] .- S_F .* m  # forall i
    dW = sum((model.prm.E .- L_H) .* n, dims=1)' .- L_F .* m  # forall j

    # 勾配を結合して返す
    dRW = vcat(dR, dW)

    return dRW
end

function short_solve(
    model::ShortModel,
    RW_ini::Vector{Float64},
    m_fixed::Vector{Float64},
    n_fixed::Matrix{Float64},
    err_short::Float64,
    short_itr::Int,
    rel::int
)
    # 短期均衡を解く関数

    RW_before = RW_ini .* model.one_2K
    p_bar_before = model.one_2K
    L_before = model.algprm.L
    t_before = 1.0
    max_value = 0
    obj_hist = Dict{Int, Float64}(0 => 1000.0)
    obj_rel_list = []

    R_ini = RW_before[1:model.prm.K]
    W_ini = RW_before[(model.prm.K + 1):end]

    if minimum(calculate_v(model, R_ini, W_ini)) < 0
        throw(ArgumentError("v must be non-negative"))
    end

    for k in 1:short_itr
        # Step 1: Backtracking
        L = backtracking(model, p_bar_before, L_before, m_fixed, n_fixed)
        dZ_SD = calculate_short_dual_df(model, p_bar_before, m_fixed, n_fixed)

        # Step 2: 解の更新
        RW = max.(model.prm.RW_proj, p_bar_before .- (dZ_SD ./ L))

        # Step 3: 収束判定
        if maximum(abs.((RW .- RW_before) ./ RW_before)) < err_short
            break
        end

        # Step 4: Adaptive restart
        if dot(calculate_short_dual_df(model, RW_before, m_fixed, n_fixed), RW .- RW_before) > 0
            t_before = 1.0
        end

        # Step 5: momentum項の計算
        t = (1.0 + sqrt(1.0 + 4.0 * t_before^2)) / 2.0
        p_bar = max.(model.prm.RW_proj, RW .+ ((t_before - 1.0) / t) * (RW .- RW_before))

        max_value = k
        RW_before = RW
        p_bar_before = p_bar
        L_before = L
        t_before = t
    end

    R_hist = RW[1:model.prm.K]
    W_hist = RW[(model.prm.K + 1):end]
    RW_hist = vcat(R_hist, W_hist)

    return R_hist, W_hist, max_value, obj_rel_list
end

function Q(
    model::ShortModel,
    p::Vector{Float64},
    p_bar::Vector{Float64},
    L_bar::Float64,
    m_fixed::Vector{Float64},
    n_fixed::Matrix{Float64}
)
    Q = calculate_Z_SD(model, p_bar, m_fixed, n_fixed) +
        dot(p .- p_bar, calculate_short_dual_df(model, p_bar, m_fixed, n_fixed)) +
        0.5 * L_bar * dot(p .- p_bar, p .- p_bar)
    return Q
end


function backtracking(
    model::ShortModel,
    p_bar::Vector{Float64},
    L::Float64,
    m_fixed::Vector{Float64},
    n_fixed::Matrix{Float64}
)
    i = 0
    L_bar = L
    dRdW = calculate_short_dual_df(model, p_bar, m_fixed, n_fixed)
    Z_SD_p_bar = calculate_Z_SD(model, p_bar, m_fixed, n_fixed)
    m_flat = vec(m_fixed)
    n_flat = vec(n_fixed)

    for k in 1:10^6
        p = max.(model.algprm.p_proj, p_bar .- dRdW ./ L_bar)
        p_R = p[1:model.prm.K]
        p_W = p[(model.prm.K + 1):end]

        v = model.prm.E * p_W .+
            model.T_power_v *
            (model.alpha_L_R .* p_R.^(-model.power_L_H_R)) *
            transpose(model.alpha_S_W .* p_W.^(-model.power_S_H_W))

        pi = model.pi_beta *
            (model.beta_L_R .* p_R.^(-model.power_L_F_R)) *
            (model.beta_S_W .* p_W.^(-model.power_S_F_W))

        Z = dot(vec(v), n_flat) +
            dot(vec(pi), m_flat) +
            model.prm.S_bar * sum(p_R)

        if Z - (Z_SD_p_bar + dot(p .- p_bar, dRdW) + 0.5 * L_bar * dot(p .- p_bar, p .- p_bar)) <= 0.0
            break
        end

        L_bar *= model.algprm.eta
        i += 1
    end

    return L_bar
end

function equilibrium(
    model::ShortModel,
    R::Vector{Float64},
    W::Vector{Float64},
    m_fixed::Vector{Float64},
    n_fixed::Matrix{Float64}
)
    RW = vcat(R, W)
    dRW = calculate_short_dual_df(model, RW, m_fixed, n_fixed)

    dR = dRW[1:model.prm.K]
    dW = dRW[(model.prm.K + 1):end]

    equ_R = maximum(R .* dR)
    equ_W = maximum(W .* dW)

    return equ_R, equ_W
end
module Short

    export demand_sparse, cal_pi, pi_noex, v_cal, Z_SD, short_dual_df, short_solve,
           Q, backtracking, equilibrium

    # include("exp_parameter.jl")

    using LinearAlgebra
    using Base.Threads
    # using .Parameter

    function demand_sparse(prm, R::Vector{Float64}, W::Vector{Float64})
    
        R_S_H = prm.alpha_S_R * (R .^ -prm.power_S_H_R)
        W_S_H = prm.alpha_S_W * (W .^ -prm.power_S_H_W)
        
        for i in 1:prm.K
            for j in 1:prm.K
                prm.RW_S_H[i, j] = R_S_H[i] * W_S_H[j]
            end
        end
    
        # RW_S_H = [R_S_H[i] * W_S_H[j] for i in 1:prm.K, j in 1:prm.K]

        R_L_H = prm.alpha_L_R * (R .^ -prm.power_L_H_R)
        W_L_H = prm.alpha_L_W * (W .^ -prm.power_L_H_W)
    
        for i in 1:prm.K
            for j in 1:prm.K
                prm.RW_L_H[i, j] = R_L_H[i] * W_L_H[j]
            end
        end
    
        # RW_L_H = [R_L_H[i] * W_L_H[j] for i in 1:prm.K, j in 1:prm.K]

        R_S_F = prm.beta_S_R * (R .^ -prm.power_S_F_R)
        W_S_F = prm.beta_S_W * (W .^ -prm.power_S_F_W)

        R_L_F = prm.beta_L_R * (R .^ -prm.power_L_F_R)
        W_L_F = prm.beta_L_W * (W .^ -prm.power_L_F_W)

        S_H = prm.T_power .* prm.RW_S_H
        L_H = prm.T_power .* prm.RW_L_H

        S_F = R_S_F .* W_S_F
        L_F = R_L_F .* W_L_F

        return S_H, L_H, S_F, L_F
    
    end


    function cal_pi(prm, R::Vector{Float64}, W::Vector{Float64}, m::Vector{Float64})
        # 利潤関数を計算する関数

        R_pi = prm.beta_L_R * (R .^ -prm.power_L_F_R)
        W_pi = prm.beta_S_W * (W .^ -prm.power_S_F_W)

        pi = prm.D * m +
             prm.pi_beta * (R_pi .* W_pi)

        return pi
    end

    function pi_noex(prm, R::Vector{Float64}, W::Vector{Float64}, m::Vector{Float64})
        # 需要関数を設定する関数

        R_pi = prm.beta_L_R * (R .^ -prm.power_L_F_R)
        W_pi = prm.beta_S_W * (W .^ -prm.power_S_F_W)

        pi = prm.pi_beta * (R_pi .* W_pi)

        return pi
    end

    function v_cal(prm, R::Vector{Float64}, W::Vector{Float64})
        # 間接効用関数を計算する関数

        R_v = prm.alpha_L_R * (R .^ -prm.power_L_H_R)
        W_v = prm.alpha_S_W * (W .^ -prm.power_S_H_W)
        RW_v = [R_v[i] * W_v[j] for i in 1:prm.K, j in 1:prm.K]

        v = prm.E * W .+
            prm.T_power_v .*
            RW_v

        return v
    end

    function Z_SD(prm, RW::Vector{Float64}, m::Vector{Float64}, n::Matrix{Float64})
    # 目的関数を計算する関数

        R = RW[1:prm.K]
        W = RW[(prm.K + 1):end]

        # 間接効用関数の計算
        v_value = v_cal(prm, R, W)

        F_value = sum(v_value .* n) +
                  dot(pi_noex(prm, R, W, m), m) +
                  prm.S_bar * sum(R)

        return F_value
    end

    function short_dual_df(prm, RW::Vector{Float64}, m::Vector{Float64}, n::Matrix{Float64})
        # 目的関数の勾配を計算する関数

        R = RW[1:prm.K]
        W = RW[(prm.K + 1):end]

        # 需要関数の計算
        S_H, L_H, S_F, L_F = demand_sparse(prm, R, W)

        # 勾配の計算
        dR = prm.S_bar .- sum(S_H .* n, dims=2)[:] .- S_F .* m  # forall i
        dW = sum((prm.E .- L_H) .* n, dims=1)[:] .- L_F .* m  # forall j

        # 勾配を結合して返す
        dRW = vcat(dR, dW)

        return dRW
    end

    function short_solve(
        prm,
        algprm, 
        RW_ini::Float64,
        m_fixed::Vector{Float64},
        n_fixed::Matrix{Float64},
        err_short::Float64,
        short_itr::Int
    )
        # 短期均衡を解く関数
        RW_before = RW_ini * prm.one_2K
        p_bar_before = prm.one_2K
        L_before = algprm.L
        t_before = 1.0
        max_value = 0
        obj_hist = Dict{Int, Float64}(0 => 1000.0)
        obj_rel_list = []
        RW_box = []
        max_value_box = []

        R_ini = RW_before[1:prm.K]
        W_ini = RW_before[(prm.K + 1):end]
    
        g = 1

        if minimum(v_cal(prm, R_ini, W_ini)) < 0
            throw(ArgumentError("v must be non-negative"))
        end

        for k in 1:short_itr
        
            # Step 1: Backtracking
            L = backtracking(prm, algprm, p_bar_before, L_before, m_fixed, n_fixed)
            dZ_SD = short_dual_df(prm, p_bar_before, m_fixed, n_fixed)

            # Step 2: 解の更新
            RW = max.(prm.RW_proj, p_bar_before .- (dZ_SD ./ L))

            # Step 3: 収束判定
            if maximum(abs.((RW .- RW_before) ./ RW_before)) < err_short
                RW_box = vcat(RW_box, RW)
                break
            end

            # Step 4: Adaptive restart
            if dot(short_dual_df(prm, RW_before, m_fixed, n_fixed), RW .- RW_before) > 0
                t_before = 1.0
            end

            # Step 5: momentum項の計算
            t = (1.0 + sqrt(1.0 + 4.0 * t_before^2)) / 2.0
            p_bar = max.(prm.RW_proj, RW .+ ((t_before - 1.0) / t) * (RW .- RW_before))

            RW_before = RW
            p_bar_before = p_bar
            L_before = L
            t_before = t
            g += 1
        end

        R_hist = RW_box[1:prm.K]
        W_hist = RW_box[(prm.K + 1):end]
    
        R_hist = convert(Vector{Float64}, R_hist)
        W_hist = convert(Vector{Float64}, W_hist)

        return R_hist, W_hist, g, obj_rel_list
    end

    # function Q(
    #     model::ShortModel,
    #     p::Vector{Float64},
    #     p_bar::Vector{Float64},
    #     L_bar::Float64,
    #     m_fixed::Vector{Float64},
    #     n_fixed::Matrix{Float64}
    # )
    #     Q = Z_SD(p_bar, m_fixed, n_fixed) +
    #         dot(p .- p_bar, calculate_short_dual_df(model, p_bar, m_fixed, n_fixed)) +
    #         0.5 * L_bar * dot(p .- p_bar, p .- p_bar)
    #     return Q
    # end


    function backtracking(prm, algprm, p_bar::Vector{Float64}, L::Float64, m_fixed::Vector{Float64}, n_fixed::Matrix{Float64})

        i = 0
        L_bar = L
        dRdW = short_dual_df(prm, p_bar, m_fixed, n_fixed)
        Z_SD_p_bar = Z_SD(prm, p_bar, m_fixed, n_fixed)

        for k in 1:10^6
            p = max.(algprm.p_proj, p_bar .- dRdW ./ L_bar)

            if Z_SD(prm, p, m_fixed, n_fixed) - (Z_SD_p_bar + dot(p .- p_bar, dRdW) + 0.5 * L_bar * dot(p .- p_bar, p .- p_bar)) <= 0.0
                break
            end

            L_bar *= algprm.eta
            i += 1
        end

        return L_bar
    end

#     function equilibrium(
#         R::Vector{Float64},
#         W::Vector{Float64},
#         m_fixed::Vector{Float64},
#         n_fixed::Matrix{Float64}
#     )
#         RW = vcat(R, W)
#         dRW = short_dual_df(RW, m_fixed, n_fixed)

#         dR = dRW[1:short.prm.K]
#         dW = dRW[(short.prm.K + 1):end]

#         equ_R = maximum(R .* dR)
#         equ_W = maximum(W .* dW)

#         return equ_R, equ_W
#     end

end #module end
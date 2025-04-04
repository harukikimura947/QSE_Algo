module par_Short

    export demand_sparse, cal_pi, pi_noex, v_cal, Z_SD, short_dual_df, short_solve,
           backtracking, equilibrium

    using LinearAlgebra
    
    function demand_sparse(prm, R::Vector{Float64}, W::Vector{Float64})

        R_S_H = prm.alpha_S_R * (R .^ -prm.power_S_H_R)
        W_S_H = prm.alpha_S_W * (W .^ -prm.power_S_H_W)

        RW_S_H = R_S_H * W_S_H'

        R_L_H = prm.alpha_L_R * (R .^ -prm.power_L_H_R)
        W_L_H = prm.alpha_L_W * (W .^ -prm.power_L_H_W)

        RW_L_H = R_L_H * W_L_H'

        R_S_F = prm.beta_S_R * (R .^ -prm.power_S_F_R)
        W_S_F = prm.beta_S_W * (W .^ -prm.power_S_F_W)

        R_L_F = prm.beta_L_R * (R .^ -prm.power_L_F_R)
        W_L_F = prm.beta_L_W * (W .^ -prm.power_L_F_W)
    
        S_H = prm.T_power .* RW_S_H
        L_H = prm.T_power .* RW_L_H

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

        R_pi = prm.beta_L_R * (R .^ -prm.power_L_F_R)
        W_pi = prm.beta_S_W * (W .^ -prm.power_S_F_W)

        pi = prm.pi_beta * (R_pi .* W_pi)

        return pi
    end

    function v_cal(prm, R::Vector{Float64}, W::Vector{Float64})
        # 間接効用関数を計算する関数

        R_v = prm.alpha_L_R * (R .^ -prm.power_L_H_R)
        W_v = prm.alpha_S_W * (W .^ -prm.power_S_H_W)
        RW_v = R_v * W_v'

        v = prm.E * W .+
            prm.T_power_v .*
            RW_v

        return v
    end

    function Z_SD(prm, RW::Vector{Float64}, m::Vector{Float64}, n::Matrix{Float64})
    # 目的関数を計算する関数

        R = RW[1:prm.K]
        W = RW[(prm.K + 1):end]
    
        R_v = prm.alpha_L_R * (R .^ -prm.power_L_H_R)
        W_v = prm.alpha_S_W * (W .^ -prm.power_S_H_W)

        # v_value = v_cal(prm, R, W)
    
        vn = 0
        for i in 1:prm.K
            R_v_i = R_v[i]  # ループ内でアクセスを減らすために変数に格納
            for j in 1:prm.K
                vn += (prm.E * W[j] + prm.T_power_v[i, j] * R_v_i * W_v[j]) * n[i, j]
            end
        end

        # F_value = sum(v_value .* n) +
        #           dot(pi_noex(prm, R, W, m), m) +
        #           prm.S_bar * sum(R)
    
        F_value = vn +
                  dot(pi_noex(prm, R, W, m), m) +
                  prm.S_bar * sum(R)

        return F_value
    end

    function short_dual_df(prm, dR::Vector{Float64}, dW::Vector{Float64}, RW::Vector{Float64}, m::Vector{Float64}, n::Matrix{Float64})

        R = @view RW[1:prm.K]
        W = @view RW[(prm.K + 1):end]

        S_H_R = (prm.alpha_S_R * (R .^ -prm.power_S_H_R))
        S_H_W = W .^ -prm.power_S_H_W

        L_H_R = R .^ -prm.power_L_H_R
        L_H_W = (prm.alpha_L_W * (W .^ -prm.power_L_H_W))

        S_F = (prm.beta_S_R * (R .^ -prm.power_S_F_R)) .* (prm.beta_S_W * (W .^ -prm.power_S_F_W))
        L_F = (prm.beta_L_R * (R .^ -prm.power_L_F_R)) .* (prm.beta_L_W * (W .^ -prm.power_L_F_W))

        for i in 1:prm.K
            Sn = sum(@views prm.T_power[i, :] .* S_H_W .* n[i, :])  # ベクトル演算
            dR[i] = prm.S_bar - S_H_R[i] * prm.alpha_S_W * Sn - S_F[i] * m[i]
        end

        for j in 1:prm.K
            Ln = sum(@views prm.T_power[:, j] .* L_H_R .* n[:, j])  # ベクトル演算
            sum_n = sum(@views n[:, j])
            dW[j] = prm.E * sum_n - L_H_W[j] * prm.alpha_L_R * Ln - L_F[j] * m[j]
        end

        return vcat(dR, dW)

    end

#     function short_dual_df(prm, RW::Vector{Float64}, m::Vector{Float64}, n::Matrix{Float64})
#         # 目的関数の勾配を計算する関数

#         R = RW[1:prm.K]
#         W = RW[(prm.K + 1):end]

#         # 需要関数の計算
#         S_H, L_H, S_F, L_F = demand_sparse(prm, R, W)
    
#         E_L_H = prm.E .- L_H
#         S_H_n = cu(S_H) .* cu(n)
#         L_H_n = cu(E_L_H) .* cu(n)

#         # 勾配の計算
#         # dR = prm.S_bar .- sum(S_H .* n, dims=2)[:] .- S_F .* m  # forall i
#         # dW = sum((prm.E .- L_H) .* n, dims=1)[:] .- L_F .* m  # forall j
    
#         dR = prm.S_bar .- sum(S_H_n, dims=2)[:] .- S_F .* m  # forall i
#         dW = sum(L_H_n .* n, dims=1)[:] .- L_F .* m  # forall j

#         # 勾配を結合して返す
#         dRW = vcat(dR, dW)

#         return dRW
#     end

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
        RW_before = ones(2 * prm.K)
        RW = zeros(2 * prm.K)
        dZ_SD = zeros(2 * prm.K)
        p_bar_before = ones(2 * prm.K)
        dRdW_dot = ones(2 * prm.K)
        RW_diff = ones(2 * prm.K)
        L_before = algprm.L
        t_before = 1.0
        max_value = 0
        obj_rel_list = []
        RW_box = []

        R_ini = RW_before[1:prm.K]
        W_ini = RW_before[(prm.K + 1):end]
    
        g = 1

        if minimum(v_cal(prm, R_ini, W_ini)) < 0
            throw(ArgumentError("v must be non-negative"))
        end

        dR = zeros(prm.K)
        dW = zeros(prm.K)

        for k in 1:short_itr
        
            println("=====================================================")
            println("iteration:", k)
        
            # println("p_bar_before:", p_bar_before)

            # Step 1: Backtracking
            L = backtracking(prm, algprm, dR, dW, p_bar_before, L_before, m_fixed, n_fixed)
        
            println("L:", L)
        
            # println("Before allocation: ", Base.gc_bytes(), " bytes")
        
            dZ_SD .= short_dual_df(prm, dR, dW, p_bar_before, m_fixed, n_fixed)
        
            # println("dRdW:", dZ_SD)
        
#             println("After allocation: ", Base.gc_bytes(), " bytes")

#             # メモリを解放
#             # dZ_SD = nothing
#             GC.gc()

#             println("After GC: ", Base.gc_bytes(), " bytes")
        
            # if k % 10 == 0
            #     GC.gc()  # 定期的にGCを呼び出す
            #     println("Iteration $k: Memory used: ", Base.gc_bytes(), " bytes")
            # end
        
            # Step 2: 解の更新
            RW .= max.(prm.RW_proj, p_bar_before .- (dZ_SD ./ L))
        
            # println("RW_before:", RW_before)
            # println("RW:", RW)

            # Step 3: 収束判定
            if maximum(abs.((RW .- RW_before) ./ RW_before)) < err_short
                RW_box = vcat(RW_box, RW)
                break
            end
        
            # Step 4: Adaptive restart
            # dRdW_dot = short_dual_df(prm, dR, dW, RW_before, m_fixed, n_fixed);
            # RW_diff = RW .- RW_before
            # println("dRdW_dot:", dRdW_dot)
            # println("RW_diff:", RW_diff)
            # println("RW_dot:", dot(dRdW_dot, RW_diff))
        
            if dot(short_dual_df(prm, dR, dW, RW_before, m_fixed, n_fixed), RW .- RW_before) > 0
                t_before = 1.0
                println("Adaptive_on")
            end
        
#             df = short_dual_df(prm, dR, dW, RW_before, m_fixed, n_fixed)
        
#             println("maximum(abs.(df)):", maximum(abs.(df)))
        
#             # Step 3: 収束判定
#             if maximum(abs.(df)) < err_short
#                 RW_box = vcat(RW_box, RW)
#                 break
#             end
        
#             # Step 4: Adaptive restart
#             if dot(df, RW .- RW_before) > 0
#                 t_before = 1.0
#             end

            # Step 5: momentum項の計算
            t = (1.0 + sqrt(1.0 + 4.0 * t_before^2)) / 2.0
            p_bar = max.(prm.RW_proj, RW .+ ((t_before - 1.0) / t) * (RW .- RW_before))

            RW_before .= RW
            p_bar_before .= p_bar
            L_before = L
            t_before = t
            g += 1
        
            # println("RW_before:", RW_before)
            # println("t_before:", t_before)
            # println("L_before:", L_before)
        
            if g == short_itr + 1
                RW_box = vcat(RW_box, RW)
            end
        
            # println("T_object:", objectid(prm.T))
        
        end
    
        println("g:", g)

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


    function backtracking(prm, algprm, dR::Vector{Float64}, dW::Vector{Float64}, p_bar::Vector{Float64}, L::Float64, m_fixed::Vector{Float64}, n_fixed::Matrix{Float64})

        i = 0
        L_bar = L
        dRdW = short_dual_df(prm, dR, dW, p_bar, m_fixed, n_fixed)
        Z_SD_p_bar = Z_SD(prm, p_bar, m_fixed, n_fixed)
    
        # println("Z_SD_p_bar:", Z_SD_p_bar)
        # println("dRdW:", dRdW)
    
        p = similar(p_bar)  # 配列を事前に確保

        for k in 1:10^6
            # p = max.(algprm.p_proj, p_bar .- dRdW ./ L_bar)
            p .= max.(algprm.p_proj, p_bar - dRdW / L_bar)  # 新しい配列の割り当てを防ぐ
            Z = Z_SD(prm, p, m_fixed, n_fixed)
        
            # println("p:", p)
            # println("Z:", Z)

            if Z - (Z_SD_p_bar + dot(p .- p_bar, dRdW) + 0.5 * L_bar * dot(p .- p_bar, p .- p_bar)) <= 0.0
                println("btra_itr:", i)
                break
            end

            L_bar *= algprm.eta
            i += 1
        end

        return L_bar
    end

    function equilibrium(
        prm, 
        dR::Vector{Float64},
        dW::Vector{Float64},
        R::Vector{Float64},
        W::Vector{Float64},
        m_fixed::Vector{Float64},
        n_fixed::Matrix{Float64}
    )

        RW = vcat(R, W)
        dRW = short_dual_df(prm, dR, dW, RW, m_fixed, n_fixed)

        dR = dRW[1:prm.K]
        dW = dRW[(prm.K + 1):end]

        equ_R = maximum(R .* dR)
        equ_W = maximum(W .* dW)

        return equ_R, equ_W
    end

end #module end
module par_Long

    export par_bond, par_breakingdown, par_Z_LP, par_long_df, par_long_solve, par_logit, par_long_armijo, par_Lyapunov

    include("par_exp_Short_FISTA.jl")

    using SparseArrays
    using LinearAlgebra
    using .par_Short

    function par_bond(prm, m::Vector{Float64}, n::Matrix{Float64})

        n_list = collect(reshape(n', :))  # 行列を1次元のベクトルに変換
        mn = vcat(m, n_list)
        return mn
    end


    function par_breakingdown(prm, mn::Vector{Float64})

        K = prm.K
        m = mn[1:prm.K]
        n_list = mn[prm.K+1:end]
        n = reshape(n_list, K, K)'
        n = convert(Matrix{Float64}, n)
        return m, n
    
    end

    function par_Z_LP(prm, mn::Vector{Float64}, RW::Vector{Float64})

        m, n = par_breakingdown(prm, mn)
        Short = Z_SD(prm, RW, m, n)

        F_value = - Short - 0.5 * m' * prm.D * m + (1 / prm.theta_firm) * (dot(m, log.(m ./ prm.M))) + (1 / prm.theta_house) * sum(n .* log.(n ./ prm.N))

        return F_value
    end

#     function long_df(mn::Vector{Float64}, RW::Vector{Float64})
#         """目的関数の勾配を計算する関数

#         Parameters
#         ----------
#         mn: Vector{Float64} (K + K * K,)
#             企業・家計分布の結合リスト
#         RW: Vector{Float64} (2 * K,)
#             価格変数の結合リスト

#         Returns
#         ----------
#         dF: Vector{Float64} (K + K * K,)
#             目的関数の勾配
#         """
#         m, n = breakingdown(mn)
#         dF_m = -pi(self.short, RW[1:self.prm.K], RW[self.prm.K+1:end], m) + (1 / self.prm.theta_firm) * (log.(m ./ self.prm.M) .+ 1)
#         dF_n = -v(self.short, RW[1:lself.prm.K], RW[self.prm.K+1:end]) .+ (1 / self.prm.theta_house) * (log.(n ./ self.prm.N) .+ 1)

#         dF = bond(dF_m, dF_n)
#         return dF
#     end

    function par_long_solve(prm, algprm, m0::Vector{Float64}, n0::Matrix{Float64}, RW_ini::Float64, err_short::Float64, err_long::Float64, obj_corr::Float64, long_itr::Int)

        mn0 = par_bond(prm, m0, n0)
        mn_before = mn0
        obj_before = 0.0

        mn_box = []
        RW_box = []
        g = 1

        for k in 1:long_itr
            m_before, n_before = par_breakingdown(prm, mn_before)

            R, W, iteration, short_obj_rel = short_solve(prm, algprm, RW_ini, m_before, n_before, err_short, 100000)
        
            RW = vcat(R, W)

            mn_d = par_logit(prm, R, W, m_before)

    #         obj = Z_LP(mn_before, RW)
    #         push!(long_iteration, k)
    #         push!(obj_list, obj)

    #         obj_rel = abs(obj - obj_corr) / abs(obj_corr)
    #         push!(obj_rel_list, obj_rel)

            alpha = 0.05
            mn = (1 - alpha) * mn_before .+ alpha * mn_d
        
            if maximum(abs.((mn .- mn_before) ./ mn_before)) < err_long
                mn_box = vcat(mn_box, mn_before)
                RW_box = vcat(RW_box, RW)
                break
            end

            mn_before = mn
            # obj_before = obj
            g += 1
        end

        mn_box = convert(Vector{Float64}, mn_box)
        RW_box = convert(Vector{Float64}, RW_box)
        m_true, n_true = par_breakingdown(prm, mn_box)
        print("Long_max_value:", g)
        println("Lyapunov:", par_Lyapunov(prm, m_true, n_true, RW_box))
    
        return m_true, n_true, RW_box, g
    end

    function par_logit(prm, R::Vector{Float64}, W::Vector{Float64}, m::Vector{Float64})

        pi = cal_pi(prm, R, W, m)
        v = v_cal(prm, R, W)

        m_d = prm.M .* exp.(prm.theta_firm * pi) ./ sum(exp.(prm.theta_firm * pi))

        n_d = prm.N .* exp.(prm.theta_house * v) ./ sum(exp.(prm.theta_house * v))

        mn_d = par_bond(prm, m_d, n_d)

        return mn_d
    end

#     function long_armijo(Z_LP, dZ, mn::Vector{Float64}, mn_d::Vector{Float64}, RW::Vector{Float64}; c_1=0.5, beta=0.6)

#         t = 1.0

#         while true
#             println("Armijo")
#             if Z_LP(mn .+ t .* mn_d, RW) < Z_LP(mn, RW) + c_1 * t * (dZ(mn, RW)' * mn_d)
#                 break
#             end
#             t *= beta
#         end

#         return t
#     end

    function par_Lyapunov(prm, m::Vector{Float64}, n::Matrix{Float64}, RW::Vector{Float64})

        R = RW[1:prm.K]
        W = RW[prm.K+1:end]

        pi = cal_pi(prm, R, W, m)
        v = v_cal(prm, R, W)

        G = sum(pi .* m) + sum(v .* n) - (prm.M * log(sum(exp.(prm.theta_firm * pi))) / prm.theta_firm) - (prm.N * log(sum(exp.(prm.theta_house * v))) / prm.theta_house) - (sum(m .* log.(m ./ prm.M)) / prm.theta_firm) - (sum(n .* log.(n ./ prm.N)) / prm.theta_house)

        return G
    end

end #module end
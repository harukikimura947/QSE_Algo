module Parameter

    export st_Parameter, st_Algo_Parameter

    struct st_Parameter
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
        power_S_H_R::Float64
        power_S_H_W::Float64
        power_L_H_R::Float64
        power_L_H_W::Float64
        power_S_F_R::Float64
        power_S_F_W::Float64
        power_L_F_R::Float64
        power_L_F_W::Float64
        T_power::Matrix{Float64}
        alpha_S_R::Float64
        alpha_S_W::Float64
        alpha_L_R::Float64
        alpha_L_W::Float64
        beta_S_R::Float64
        beta_S_W::Float64
        beta_L_R::Float64
        beta_L_W::Float64
        T_power_v::Matrix{Float64}
        pi_beta::Float64
        one_K::Vector{Float64}
        one_2K::Vector{Float64}
        RW_S_H::Matrix{Float64}
        RW_L_H::Matrix{Float64}
    end

    # アルゴリズム中のパラメータを保持するクラス
    struct st_Algo_Parameter
        L::Float64
        eta::Float64
        p_proj::Float64

    end 

end #module end
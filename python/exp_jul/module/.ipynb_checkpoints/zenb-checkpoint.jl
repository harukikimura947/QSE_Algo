using LinearAlgebra
using Random
using Printf

function create_column(Col::Int, tr_row_id::Int)
    link_dic = Dict{String, Dict}()
    node_dic = Dict{String, Dict}()
    for i in 0:Col-1
        node_id = string(tr_row_id * Col + i)
        node_dic[node_id] = Dict("x_pos" => tr_row_id, "y_pos" => i)

        next_node_id = string(parse(Int, node_id) + 1)
        if i < Col - 1
            link1_id = node_id * "_" * next_node_id
            link2_id = next_node_id * "_" * node_id

            link_dic[link1_id] = Dict("from_node_id" => node_id, "to_node_id" => next_node_id)
            link_dic[link2_id] = Dict("from_node_id" => next_node_id, "to_node_id" => node_id)
        end
    end
    return link_dic, node_dic
end

function add_column(Col::Int, network_dic::Dict, tr_row_id::Int, new_link_dic::Dict, new_node_dic::Dict)
    # 新しいリンクとノードを追加
    merge!(network_dic["link_dic"], new_link_dic)
    merge!(network_dic["node_dic"], new_node_dic)

    # ノードを接続
    for node_id in keys(new_node_dic)
        left_node_id = string(parse(Int, node_id) - Col)

        link1_id = node_id * "_" * left_node_id
        link2_id = left_node_id * "_" * node_id

        network_dic["link_dic"][link1_id] = Dict("from_node_id" => node_id, "to_node_id" => left_node_id)
        network_dic["link_dic"][link2_id] = Dict("from_node_id" => left_node_id, "to_node_id" => node_id)
    end
    return network_dic
end

function make_lattice(Col::Int)
    network_dic = Dict{String, Dict}()
    network_dic["link_dic"], network_dic["node_dic"] = create_column(Col, 1)

    for row_id in 0:Col-1
        new_link_dic, new_node_dic = create_column(Col, row_id)
        network_dic = add_column(Col, network_dic, row_id, new_link_dic, new_node_dic)
    end
    return network_dic
end

function demand_sparse(R::Vector{Float64}, W::Vector{Float64})

    R_S_H = alpha_S_R * (R .^ -power_S_H_R)
    W_S_H = alpha_S_W * (W .^ -power_S_H_W)

    RW_S_H = R_S_H * W_S_H'

    R_L_H = alpha_L_R * (R .^ -power_L_H_R)
    W_L_H = alpha_L_W * (W .^ -power_L_H_W)

    RW_L_H = R_L_H * W_L_H'

    R_S_F = beta_S_R * (R .^ -power_S_F_R)
    W_S_F = beta_S_W * (W .^ -power_S_F_W)

    R_L_F = beta_L_R * (R .^ -power_L_F_R)
    W_L_F = beta_L_W * (W .^ -power_L_F_W)

    S_H = T_power .* RW_S_H
    L_H = T_power .* RW_L_H

    S_F = R_S_F .* W_S_F
    L_F = R_L_F .* W_L_F

    return S_H, L_H, S_F, L_F

end

function cal_pi(R::Vector{Float64}, W::Vector{Float64}, m::Vector{Float64})
    # 利潤関数を計算する関数

    R_pi = beta_L_R * (R .^ -power_L_F_R)
    W_pi = beta_S_W * (W .^ -power_S_F_W)

    pi = D * m +
         pi_beta * (R_pi .* W_pi)

    return pi
end

function pi_noex(R::Vector{Float64}, W::Vector{Float64}, m::Vector{Float64})
    # 需要関数を設定する関数

    R_pi = beta_L_R * (R .^ -power_L_F_R)
    W_pi = beta_S_W * (W .^ -power_S_F_W)

    pi = pi_beta * (R_pi .* W_pi)

    return pi
end

function v_cal(R::Vector{Float64}, W::Vector{Float64})
    # 間接効用関数を計算する関数

    R_v = alpha_L_R * (R .^ -power_L_H_R)
    W_v = alpha_S_W * (W .^ -power_S_H_W)
    RW_v = R_v * W_v'

    v = E * W .+
        T_power_v .*
        RW_v

    return v
end

function Z_SD(RW::Vector{Float64}, m::Vector{Float64}, n::Matrix{Float64})
# 目的関数を計算する関数

    R = RW[1:K]
    W = RW[(K + 1):end]

    R_v = alpha_L_R * (R .^ -power_L_H_R)
    W_v = alpha_S_W * (W .^ -power_S_H_W)

    # v_value = v_cal(R, W)

    vn = 0
    for i in 1:K
        R_v_i = R_v[i]  # ループ内でアクセスを減らすために変数に格納
        for j in 1:K
            vn += (E * W[j] + T_power_v[i, j] * R_v_i * W_v[j]) * n[i, j]
        end
    end

    F_value = vn +
              dot(pi_noex(R, W, m), m) +
              S_bar * sum(R)

    return F_value
end

function short_dual_df(dR::Vector{Float64}, dW::Vector{Float64}, RW::Vector{Float64}, m::Vector{Float64}, n::Matrix{Float64})

    R = @view RW[1:K]
    W = @view RW[(K + 1):end]

    S_F = (beta_S_R * (R .^ -power_S_F_R)) .* (beta_S_W * (W .^ -power_S_F_W))
    L_F = (beta_L_R * (R .^ -power_L_F_R)) .* (beta_L_W * (W .^ -power_L_F_W))
    
    for i in 1:K
        Sn = sum(@views T_power[i, :] .* W .* n[i, :])  # ベクトル演算
        dR[i] = S_bar - R[i] * Sn - S_F[i] * m[i]
    end

    for j in 1:K
        Ln = sum(@views T_power[:, j] .* R .* n[:, j])  # ベクトル演算
        sum_n = sum(@views n[:, j])
        dW[j] = E * sum_n - W[j] * Ln - L_F[j] * m[j]
    end

    return vcat(dR, dW)

end

function short_solve(
    RW_ini::Float64,
    m_fixed::Vector{Float64},
    n_fixed::Matrix{Float64},
    err_short::Float64,
    short_itr::Int
)

    # 短期均衡を解く関数
    RW_before = ones(2 * K)
    RW = zeros(2 * K)
    dZ_SD = zeros(2 * K)
    p_bar_before = ones(2 * K)
    L_before = par_L
    t_before = 1.0
    max_value = 0
    obj_rel_list = []
    RW_box = []

    R_ini = RW_before[1:K]
    W_ini = RW_before[(K + 1):end]

    g = 1

    if minimum(v_cal(R_ini, W_ini)) < 0
        throw(ArgumentError("v must be non-negative"))
    end

    dR = zeros(K)
    dW = zeros(K)

    for k in 1:short_itr

        # Step 1: Backtracking
        L = backtracking(dR, dW, p_bar_before, L_before, m_fixed, n_fixed)

        dZ_SD .= short_dual_df(dR, dW, p_bar_before, m_fixed, n_fixed)
        
        # Step 2: 解の更新
        RW .= max.(RW_proj, p_bar_before .- (dZ_SD ./ L))

        # Step 3: 収束判定
        if maximum(abs.((RW .- RW_before) ./ RW_before)) < err_short
            RW_box = vcat(RW_box, RW)
            break
        end

        # Step 4: Adaptive restart
        if dot(short_dual_df(dR, dW, RW_before, m_fixed, n_fixed), RW .- RW_before) > 0
            t_before = 1.0
        end

        # Step 5: momentum項の計算
        t = (1.0 + sqrt(1.0 + 4.0 * t_before^2)) / 2.0
        p_bar = max.(RW_proj, RW .+ ((t_before - 1.0) / t) * (RW .- RW_before))

        RW_before .= RW
        p_bar_before .= p_bar
        L_before = L
        t_before = t
        g += 1

        if g == short_itr + 1
            RW_box = vcat(RW_box, RW)
        end

    end

    println("g:", g)

    R_hist = RW_box[1:K]
    W_hist = RW_box[(K + 1):end]

    R_hist = convert(Vector{Float64}, R_hist)
    W_hist = convert(Vector{Float64}, W_hist)

    return R_hist, W_hist, g, obj_rel_list
end

function backtracking(dR, dW, p_bar::Vector{Float64}, L::Float64, m_fixed::Vector{Float64}, n_fixed::Matrix{Float64})

    i = 0
    L_bar = L
    dRdW = short_dual_df(dR, dW, p_bar, m_fixed, n_fixed)
    Z_SD_p_bar = Z_SD(p_bar, m_fixed, n_fixed)

    p = similar(p_bar)  # 配列を事前に確保

    for k in 1:10^6
        # p = max.(p_proj, p_bar .- dRdW ./ L_bar)
        p .= max.(p_proj, p_bar - dRdW / L_bar)  # 新しい配列の割り当てを防ぐ

        if Z_SD(p, m_fixed, n_fixed) - (Z_SD_p_bar + dot(p .- p_bar, dRdW) + 0.5 * L_bar * dot(p .- p_bar, p .- p_bar)) <= 0.0
            break
        end

        L_bar *= eta
        i += 1
    end

    return L_bar
end

function equilibrium(
    
    dR::Vector{Float64},
    dW::Vector{Float64},
    R::Vector{Float64},
    W::Vector{Float64},
    m_fixed::Vector{Float64},
    n_fixed::Matrix{Float64}
)

    RW = vcat(R, W)
    dRW = short_dual_df(dR, dW, RW, m_fixed, n_fixed)

    dR = dRW[1:K]
    dW = dRW[(K + 1):end]

    equ_R = maximum(R .* dR)
    equ_W = maximum(W .* dW)

    return equ_R, equ_W
end

# 初期パラメータ設定
E = 5
Col = 30
K = Col * Col
M = 1.0
N = 1.0
Scaling = 10.0 / Col
alter_T_num = 0.5

# 格子ネットワークを生成
city_network = make_lattice(Col)

# 通勤費用パラメータ
t = 0.1
# 距離抵抗パラメータ
tau = 0.5
# 総土地供給量と平均床面積
S_total = 100
S_bar = S_total / K

# 座標データの生成
Coordinate_Data = [ 
    (city_network["node_dic"][string(i)]["x_pos"] * Scaling,
     city_network["node_dic"][string(i)]["y_pos"] * Scaling) 
    for i in 0:(K - 1)
]

# 距離行列の作成
distance_matrix = [sqrt(sum((Coordinate_Data[i] .- Coordinate_Data[j]).^2))
                   for i in 1:K, j in 1:K]

# モデルパラメータ設定
theta_firm = 1.0
theta_house = 1.0
alpha_1 = 0.4
alpha_2 = 0.4
beta_1 = 0.4
beta_2 = 0.4
# L = 1.0
# eta = 5.0
par_L = 0.2
eta = 1.2
p_proj = 1e-3
RW_ini = 1.0
RW_proj = 1e-3
err_short = 1e-5
err_long = 1e-3
obj_corr = 1.0
m_per = M / K
m0 = fill(m_per, K)
n0 = fill(N / (K * K), (K, K))
method = "FISTA"
dic = "Scaling=30"

K = Col * Col
T = max.(Scaling * t * alter_T_num, t * distance_matrix)
D = exp.(-tau * distance_matrix)

power_S_H_R = (1 - alpha_2) / (1 - alpha_1 - alpha_2)
power_S_H_W = alpha_2 / (1 - alpha_1 - alpha_2)
power_L_H_R = alpha_1 / (1 - alpha_1 - alpha_2)
power_L_H_W = (1 - alpha_1) / (1 - alpha_1 - alpha_2)

power_S_F_R = (1 - beta_2) / (1 - beta_1 - beta_2)
power_S_F_W = beta_2 / (1 - beta_1 -beta_2)
power_L_F_R = beta_1 / (1 - beta_1 - beta_2)
power_L_F_W = (1 - beta_1) / (1 - beta_1 - beta_2)

T_power = (1 ./ exp.(T)) .^ (1 / (1 - alpha_1 - alpha_2))
alpha_S_R = alpha_1 ^ power_S_H_R
alpha_S_W = alpha_2 ^ power_S_H_W
alpha_L_R = alpha_1 ^ power_L_H_R
alpha_L_W = alpha_2 ^ power_L_H_W

beta_S_R = beta_1 ^ power_S_F_R
beta_S_W = beta_2 ^ power_S_F_W
beta_L_R = beta_1 ^ power_L_F_R
beta_L_W = beta_2 ^ power_L_F_W

RW_S_H = zeros(K, K)
RW_L_H = zeros(K, K)

T_power_v = (1 - alpha_1 - alpha_2) * T_power
pi_beta = 1 - beta_1 - beta_2

one_K = ones(K)
one_2K = ones(2 * K)

R_ini_list = ones(K)
W_ini_list = ones(K)
RW_ini_list = ones(2 * K)

short_itr = 10
long_itr = 1000

short_solve(RW_ini, m0, n0, err_short, short_itr)
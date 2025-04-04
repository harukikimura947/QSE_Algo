# MK_0705_Logit_S-P_ver_1_0

# 0708 logit式をUniversal Gravity型にアップデートしたい

# 短期均衡（主問題）の解をlogitダイナミクスでアップデート

# 数量変数だけで計算
# 0705 長期・短期同時均衡モデルとの解の一致を確認

println("【 短期均衡（主問題）モデルのlogitダイナミクスアップデート 】")
println("【 間接効用関数を使わずに、行う（双対問題解けなくても/作れなくても大丈夫か確認） 】")
# モデル構築・解法
module RCGE
using JuMP, Ipopt
# ==========
# パラメタ用構造体
# ==========
struct Param
    I :: Int # 地域数
    N :: Int # コモディティ数
    M :: Int # 財・サービス（産業）数
    H :: Float64 # 家計数
    G :: Vector{Float64} # 産業m(=1,⋯,M)の企業総数
    μ :: Matrix{Float64} # μ[i,n]: 地域i(=1,⋯,I)家計のコモディティn(=1,⋯,N)に対する効用の重み
    ω :: Matrix{Float64} # ω[i,n]: 地域i(=1,⋯,I)の家計のコモディティn(=1,⋯,N)の保有量
    ω_0 :: Vector{Float64} #ω_0[i]: 地域i(=1,⋯,I)の家計のnumerarire財の保有量(=初期所得)
    ϕ_0 :: Vector{Float64} # ϕ_0[i]: 地域i(=1,⋯,I)固有の固有
    β :: Float64 # 所得の限界効用
    γ :: Array{Float64, 3} # γ[m,n,i]: 地域i(=1,⋯,I)の産業m(=1,⋯,M)の生産関数に対するコモディティn(=1,⋯,N)の投入係数
    Ȳ :: Matrix{Float64} # Ȳ[i,m]: 地域i(=1,⋯,I)企業m(=1,⋯,M)の基準生産量
    κ :: Matrix{Float64} # κ[i,m]: 地域i(=1,⋯,I)企業m(=1,⋯,M)の生産の収穫低減係数
    W_0 :: Matrix{Float64} # W_0[i,n]: 地域i(=1,⋯,I)の不在所有者のコモディティn(=1,⋯,N)の保有量
    T :: Array{Float64, 3} # T[i,j,m]: 財mの地域i(=1,⋯,I)から地域j(=1,⋯,I)への財m(=1,⋯,M)の輸送費用
    θ :: Float64 # 家計の異質性パラメタ
    η :: Vector{Float64} # η[m]: 企業m(=1,⋯,M)の異質性パラメタ
    σ :: Vector{Float64} # σ[m]: 財m(=1,⋯,M)の輸送・販売業者の異質性パラメタ
    # コンストラクタ
    Param(I, N, M, H, G, μ, ω, ω_0, ϕ_0, β, γ, Ȳ, κ, W_0, T, θ, η, σ) = new(
        I, N, M, H, G, μ, ω, ω_0, ϕ_0, β, γ, Ȳ, κ, W_0, T, θ, η, σ
    )
end

# ==========
# 共用関数
# ==========
# ----------
# 家計
# ----------
# 直接効用関数
function _u(prm,c,i)
    return sum( prm.μ[i,n]*log(c[i,n]) for n in 1:prm.N )
end
# 間接効用関数
function _v(prm,p,i)
    return ( - sum( prm.μ[i,n]*log(p[i,n]) for n in 1:prm.N ) + 
            prm.β * sum( p[i,n]*prm.ω[i,n] for n in 1:prm.N ) + prm.ϕ_0[i] + 
            sum(prm.μ[i,n]*log.(prm.μ[i,n]) for n in 1:prm.N) + prm.β*prm.ω_0[i]- log(prm.β) - 1)
end
# 需要関数
function _c_opt(prm,p,i,n)
    return ( 1/prm.β*prm.μ[i,n]/p[i,n] )
end

# ----------
# 企業
# ----------
# 生産関数
function _y(prm,b,i,m)
    return ( prm.Ȳ[i,m]*prod( b[i,m,n]^prm.γ[m,n,i] for n in 1:prm.N )^prm.κ[i,m] )
end
# 投入要素の価格指数
function _P_index(prm,p,i,m)
    return prod( ( p[i,n]/prm.γ[m,n,i] )^prm.γ[m,n,i] for n in 1:prm.N )
end
# 最適生産量関数
function _y_opt(prm,p,i,m)
    return ( prm.Ȳ[i,m]^( 1/(1-prm.κ[i,m]) ) *
            ( prm.κ[i,m]*p[i,m]/_P_index(prm,p,i,m) )^( prm.κ[i,m]/(1-prm.κ[i,m]) ) )
end
# 最適投入量関数
function _b_opt(prm,p,i,m,n)
    return prm.κ[i,m]*p[i,m]*_y_opt(prm,p,i,m)*prm.γ[m,n,i]/p[i,n]
end
# 利潤関数
function _Π(prm,p,i,m)
    return (1-prm.κ[i,m])*p[i,m]*_y_opt(prm,p,i,m)
end
# ----------
# 輸送販売業者
# ----------
# 最適輸送量関数
function _z_opt(prm,p,i,j,m)
    return exp(prm.σ[m]*( p[j,m]-p[i,m]-prm.T[i,j,m] ))
end
# 利潤関数
function _Π_t(prm,p,m)
    return 1/prm.σ[m]*sum( _z_opt(prm,p,i,j,m) for i in 1:prm.I for j in 1:prm.I )
end


# ==========
# モデル構築関数
# ==========
# ----------
# 短期モデル: 双対問題
# ----------
function short_dual(prm, h, g)
    model = JuMP.Model(Ipopt.Optimizer)
    # ---
    # 未知変数（財価格）
    @variable(model, p[i=1:prm.I,n=1:prm.N] >= 0)
    # ---
    # 各種関数
    @expressions(model, begin
            # ■ 家計部門
            # 間接効用関数
            v[i=1:prm.I], _v(prm,p,i)
            # 最適需要関数
            c[i=1:prm.I,n=1:prm.N], _c_opt(prm,p,i,n)
            # ■ 企業部門
            # 価格指数
            P_index[i=1:prm.I,m=1:prm.M], _P_index(prm,p,i,m)
            # 最適生産量関数
            y[i=1:prm.I,m=1:prm.M], _y_opt(prm,p,i,m)
            # 最適投入量関数
            b[i=1:prm.I,m=1:prm.M,n=1:prm.N], _b_opt(prm,p,i,m,n)
            # 利潤関数
            Π[i=1:prm.I,m=1:prm.M], _Π(prm,p,i,m)
            # ■ 輸送販売部門
            # 最適輸送量
            z[i=1:prm.I,j=1:prm.I,m=1:prm.M], _z_opt(prm,p,i,j,m)
            # 利潤相当
            Π_t[m=1:prm.M], _Π_t(prm,p,m)
            # ■ 需給ギャップ
            # 家計部門の総超過供給
            W_h[i=1:prm.I,n=1:prm.N], h[i]*(prm.ω[i,n] - c[i,n])
            # 企業mの超過供給
            w_f[i=1:prm.I,m=1:prm.M,n=1:prm.N], ifelse(m==n, y[i,m], 0) - b[i,m,n]
            # 企業部門の総超過供給
            W_f[i=1:prm.I,n=1:prm.N], sum( g[i,m]*w_f[i,m,n] for m in 1:prm.M )
            # 輸送部門の総超過供給
            W_t[i=1:prm.I,n=1:prm.N], if n <= prm.M sum( z[j,i,n]-z[i,j,n] for j in 1:prm.I if i!=j) else 0 end
            # 経済全体の総超過供給
            W[i=1:prm.I,n=1:prm.N],  W_h[i,n]+W_f[i,n]+W_t[i,n]+prm.W_0[i,n]
        end)
    # ---
    # 目的関数
    @expressions(model, begin
            # 家計部門
            Z_h, 1/prm.β * sum( h[i] * v[i] for i in 1:prm.I )
            # 企業部門
            Z_f, sum( g[i,m]*Π[i,m] for i in 1:prm.I for m in 1:prm.M )
            # 輸送販売部門
            Z_t, sum( Π_t[m] for m in 1:prm.M )
            # 不在保有者部門
            Z_a, sum( prm.W_0[i,n]*p[i,n] for i in 1:prm.I for n in 1:prm.N )
        end)
    # 生産関数が規模に関して一次同次な場合も取り扱えるように
    if all(prm.κ .< 1.0)
        @objective(model, Min, Z_h + Z_f + Z_a + Z_t)    
    else
        @objective(model, Min, Z_h + Z_a + Z_t) 
        @constraint(model, zero_profit[i=1:prm.I,m=1:prm.M],
            p[i,m] - P_index[i,m]/prm.Ȳ[i,m] <= 0
            )
    end
    return model
end

# ----------
# 長期・短期同時モデル: 双対問題
# ----------
function dual(prm)
    model = JuMP.Model(Ipopt.Optimizer)
    # ---
    # 未知変数（財価格）
    @variable(model, p[i=1:prm.I,n=1:prm.N] >= 0)
    # ---
    # 各種関数
    @expressions(model, begin
            # ■ 家計部門
            # 間接効用関数
            v[i=1:prm.I], _v(prm,p,i)
            # 最大間接効用関数
            V, 1/prm.θ*log( sum( exp(prm.θ*v[i]) for i in 1:prm.I) )
            # 最適需要関数
            c[i=1:prm.I,n=1:prm.N], _c_opt(prm,p,i,n)
            # 家計数
            h[i=1:prm.I], prm.H*exp(prm.θ*v[i])/exp(prm.θ*V)
            # ■ 企業部門
            # 価格指数
            P_index[i=1:prm.I,m=1:prm.M], _P_index(prm,p,i,m)
            # 最適生産量関数
            y[i=1:prm.I,m=1:prm.M], _y_opt(prm,p,i,m)
            # 最適投入量関数
            b[i=1:prm.I,m=1:prm.M,n=1:prm.N], _b_opt(prm,p,i,m,n)
            # 利潤関数
            Π[i=1:prm.I,m=1:prm.M], _Π(prm,p,i,m)
            # 最大利潤関数
            Ψ[m=1:prm.M], 1/prm.η[m]*log( sum(exp(prm.η[m]*Π[i,m]) for i in 1:prm.I) )
            # 企業数
            g[i=1:prm.I,m=1:prm.M], prm.G[m]*exp(prm.η[m]*Π[i,m])/exp(prm.η[m]*Ψ[m])
            # ■ 輸送販売部門
            # 最適輸送量
            z[i=1:prm.I,j=1:prm.I,m=1:prm.M], _z_opt(prm,p,i,j,m)
            # 利潤相当
            Π_t[m=1:prm.M], _Π_t(prm,p,m)
            # ■ 需給ギャップ
            # 家計部門の総超過供給
            W_h[i=1:prm.I,n=1:prm.N], h[i]*(prm.ω[i,n] - c[i,n])
            # 企業mの超過供給
            w_f[i=1:prm.I,m=1:prm.M,n=1:prm.N], ifelse(m==n, y[i,m], 0) - b[i,m,n]
            # 企業部門の総超過供給
            W_f[i=1:prm.I,n=1:prm.N], sum( g[i,m]*w_f[i,m,n] for m in 1:prm.M )
            # 輸送部門の総超過供給
            W_t[i=1:prm.I,n=1:prm.N], if n <= prm.M sum( z[j,i,n]-z[i,j,n] for j in 1:prm.I if i!=j) else 0 end
            # 経済全体の総超過供給
            W[i=1:prm.I,n=1:prm.N],  W_h[i,n]+W_f[i,n]+W_t[i,n]+prm.W_0[i,n]
        end)
    # ---
    # 目的関数
    @expressions(model, begin
            # 家計部門
            Z_h, 1/prm.β * prm.H * V
            # 企業部門
            Z_f, sum( prm.G[m]*Ψ[m] for m in 1:prm.M )
            # 輸送販売部門
            Z_t, sum( Π_t[m] for m in 1:prm.M )
            # 不在保有者部門
            Z_a, sum( prm.W_0[i,n]*p[i,n] for i in 1:prm.I for n in 1:prm.N)
        end)
    # ---
    # 生産関数が規模に関して一次同次な場合も取り扱えるように
    if all(prm.κ .< 1.0)
        @objective(model, Min, Z_h + Z_f + Z_a + Z_t)    
    else
        @objective(model, Min, Z_h + Z_a + Z_t) 
        @constraint(model, zero_profit[i=1:prm.I,m=1:prm.M],
            p[i,m] - P_index[i,m]/prm.Ȳ[i,m] <= 0
            )
    end         
    return model
end

# ----------
# 短期モデル: 主問題
# ----------
function short_primal(prm, h, g)
    model = JuMP.Model(Ipopt.Optimizer)
    # ---
    # 未知変数（数量）
    @variables(model, begin
            c[i=1:prm.I,n=1:prm.N] >= 0
            b[i=1:prm.I,m=1:prm.M,n=1:prm.N] >= 0
            z[i=1:prm.I,j=1:prm.I,m=1:prm.M] >= 0
        end)
    # ---
    # 各種関数
    @expressions(model, begin
            # ■ 家計部門
            # 直接効用関数
            u[i=1:prm.I], _u(prm,c,i)
            # ■ 企業部門
            # 生産関数
            y[i=1:prm.I,m=1:prm.M], _y(prm,b,i,m)
            # ■ 輸送販売部門（特になし）
            # ■ 需給ギャップ
            # 家計部門の総超過供給
            W_h[i=1:prm.I,n=1:prm.N], h[i]*(prm.ω[i,n] - c[i,n])
            # 企業mの超過供給
            w_f[i=1:prm.I,m=1:prm.M,n=1:prm.N], ifelse(m==n, y[i,m], 0) - b[i,m,n]
            # 企業部門の総超過供給
            W_f[i=1:prm.I,n=1:prm.N], sum( g[i,m]*w_f[i,m,n] for m in 1:prm.M )
            # 輸送部門の総超過供給
            W_t[i=1:prm.I,n=1:prm.N], if n <= prm.M sum( z[j,i,n]-z[i,j,n] for j in 1:prm.I if i!=j) else 0 end
            # 経済全体の総超過供給
            W[i=1:prm.I,n=1:prm.N],  W_h[i,n]+W_f[i,n]+W_t[i,n]+prm.W_0[i,n]
        end)
    # ---
    # 制約条件
    @constraints(model, begin
            # 同一地域内では輸送は行なわれない
            zero_z[i=1:prm.I,m=1:prm.M], z[i,i,m] == 0
            # 超過供給は非負
            excess_supply[i=1:prm.I,n=1:prm.N], W[i,n] >= 0
        end)
        
    # ---
    # 目的関数
    @expressions(model, begin
            # 家計部門
            Z_h, 1/prm.β*sum( h[i] * u[i] for i in 1:prm.I )
            # 輸送販売部門
            Z_t, - sum( z[i,j,m] * prm.T[i,j,m]
                for i in 1:prm.I, j in 1:prm.I if i!=j for m in 1:prm.M )
            # 輸送部門エントロピー
            H_t, - sum( ifelse(z[i,j,m]>0, 1/prm.σ[m]*z[i,j,m]*(log(z[i,j,m])-1), 0) 
                for i in 1:prm.I, j in 1:prm.I if i!=j for m in 1:prm.M )
        end)
    @objective(model, Max, Z_h+Z_t+H_t)
    return model
end


# ----------
# 長期・短期同時モデル: 主問題
# ----------
function primal(prm)
    model = JuMP.Model(Ipopt.Optimizer)
    # ---
    # 未知変数（数量変数）
    @variables(model, begin
            h[i=1:prm.I] >= 0
            g[i=1:prm.I,m=1:prm.M] >=0
            c[i=1:prm.I,n=1:prm.N] >= 0
            b[i=1:prm.I,m=1:prm.M,n=1:prm.N] >= 0
            z[i=1:prm.I,j=1:prm.I,m=1:prm.M] >= 0
        end)
    # ---
    # 各種関数
    @expressions(model, begin
            # ■ 家計部門
            # 直接効用関数
            u[i=1:prm.I], _u(prm,c,i)
            # ■ 企業部門
            # 生産関数
            y[i=1:prm.I,m=1:prm.M], _y(prm,b,i,m)
            # ■ 輸送販売部門（特になし）
            # ■ 需給ギャップ
            # 家計部門の総超過供給
            W_h[i=1:prm.I,n=1:prm.N], h[i]*(prm.ω[i,n] - c[i,n])
            # 企業mの超過供給
            w_f[i=1:prm.I,m=1:prm.M,n=1:prm.N], ifelse(m==n, y[i,m], 0) - b[i,m,n]
            # 企業部門の総超過供給
            W_f[i=1:prm.I,n=1:prm.N], sum( g[i,m]*w_f[i,m,n] for m in 1:prm.M )
            # 輸送部門の総超過供給
            W_t[i=1:prm.I,n=1:prm.N], if n <= prm.M sum( z[j,i,n]-z[i,j,n] for j in 1:prm.I if i!=j) else 0 end
            # 経済全体の総超過供給
            W[i=1:prm.I,n=1:prm.N],  W_h[i,n]+W_f[i,n]+W_t[i,n]+prm.W_0[i,n]
        end)
    # ---
    # 制約条件
    @constraints(model, begin
            # 同一地域内では輸送は行なわれない
            zero_z[i=1:prm.I,m=1:prm.M], z[i,i,m] == 0
            # 超過供給は非負
            excess_supply[i=1:prm.I,n=1:prm.N], W[i,n] >= 0
            # 家計数の総和は一定
            household_reservation, sum(h[i] for i in 1:prm.I) == prm.H
            # 企業数の総和は一定
            firm_reservation[m=1:prm.M], sum(g[i,m] for i in 1:prm.I) == prm.G[m]
        end)
   # ---
    # 目的関数
    @expressions(model, begin
            # 家計部門
            Z_h, 1/prm.β*sum( h[i] * u[i] for i in 1:prm.I )
            # 輸送販売部門
            Z_t, - sum( z[i,j,m] * prm.T[i,j,m]
                for i in 1:prm.I, j in 1:prm.I if i!=j for m in 1:prm.M )
            # 輸送部門エントロピー
            H_t, - sum( ifelse(z[i,j,m]>0, 1/prm.σ[m]*z[i,j,m]*(log(z[i,j,m])-1), -Inf) 
                for i in 1:prm.I, j in 1:prm.I if i!=j for m in 1:prm.M )
            # 家計の立地選択エントロピー
            H_h, - sum( 1/prm.θ*h[i]*(log(h[i]/prm.H)-1) for i in 1:prm.I )
            # 企業の立地選択エントロピー
            H_f, - sum( 1/prm.η[m]*g[i,m]*(log(g[i,m]/prm.G[m])-1) for i in 1:prm.I for m in 1:prm.M )
        end)
    @objective(model, Max, Z_h+Z_t+H_h+H_f+H_t)
    return model
end

end #このエンドはなんだ？


# 3地域・2産業
# パラメタ設定
import .RCGE
using JuMP
I = 3 # 地域数
N = 4 # コモディティ数(n=N-1:土地, n=N:労働)
M = 2 # 財・サービスの数

# Number of households and firms
H = 1
G = [.6, .4]

# Share parameter of utility function
μ = vcat(fill(
    [.35 .25 .25 .15]
    ,I)...)
# 家計の初期生産要素保有量（家計はn=Nの労働を1単位保有)
ω = zeros(I,N)
ω[:,N] .= 1
ω_0 = ones(I)*0
# 地域固有の効用
ϕ_0 = zeros(I)
# 所得の限界効用
β = 1
# 企業の投入係数
γ = cat(fill([
            .10 .05 .05 .80; # 産業1への投入係数
            .05 .30 .05 .60; # 産業2への投入係数
], I)..., dims=3)
# 企業の基準生産力
Ȳ = [
    2.0 1.0; # 地域1は産業1のポテンシャルが高い
    1.0 2.0; # 地域2は産業2のポテンシャルが高い
    1.2 1.2; # 地域3は産業3のポテンシャルが高い
] * 1 
# 企業の生産収穫低減係数
κ = vcat(fill(
        [.9 .5]
        ,I)...)

# 不在所有者の生産要素保有量（不在所有者はm=3の土地のみを供給）
# 不在所有者の生産要素保有量（不在所有者はn=N-1の土地のみを供給）
S_L = [4 ,5, 1]
W_0 = zeros(I,N)
W_0[:,N-1] = S_L/sum(S_L)
# 地域間輸送費用
T_0 = [
    Inf 10 12;
    10 Inf 8;
    12 8 Inf
] 
T = cat([T_0, T_0./2]..., dims=3)

# 効用・利潤のバラツキの大きさ
θ = 10
η = ones(M)*20
σ = ones(M)*1

# パラメタのインスタンスを作成
prm = RCGE.Param(I, N, M, H, G, 
    μ, ω, ω_0, ϕ_0, β, γ, Ȳ, κ, W_0, T, θ, η, σ)

# 短期均衡の求解およびlogitダイナミクス適用
# ステップ幅を設定
δ = 0.1
# 収束判定の閾値
ε = 1e-8

# 初期立地パターンを与える
h_0 = ones(I)*H/I
g_0 = vcat(fill(G',I)...)/I

# 初期値の設定
h_n = h_0 # n回目の更新における家計の立地パターン
g_n = g_0 # n回目の更新における企業の立地パターン

# iter = 0
# while true
#     iter += 1
    # モデルを解く
    sp = RCGE.short_primal(prm, h_n, g_n)
    
    # 問題を解く
    set_attribute(sp, "print_level", 0)
    set_attribute(sp, "compl_inf_tol", 1e-8)
    optimize!(sp)
    println("Iteration $iter: ", termination_status(sp))

    # 家計の直接効用の計算 OK = 間接効用関数に双対変数（コモディティ価格）を代入したもの
    c_0 = sum(dual.(sp[:excess_supply]) .* (prm.ω-value.(sp[:c])), dims=2) + prm.ω_0
    u_d = value.(sp[:u]) + prm.β*c_0

    # 企業の利潤関数の計算（生産関数ベース）OK
    y_d = [dual.(sp[:excess_supply])[i,m] .* value.(sp[:y])[i,m] - sum(dual.(sp[:excess_supply])[i,n] .* value.(sp[:b])[i,m,n] for n in 1:prm.N) for i in 1:prm.I, m in 1:prm.M]


    # 家計が地域iを選択する確率(Universal Gravity型)
    P_h = [(u_d[i]^θ)./sum(u_d[j]^θ for j in 1:prm.I) for i in 1:prm.I]
    # 家計の立地需要関数
    F_h = [P_h[i]*H for i in 1:prm.I]

    # 企業が地域iを選択する確率(Universal Gravity型)
    P_f = [(y_d[i,m]^η[m])/sum(y_d[i,m]^η[m] for i in 1:prm.I) for i in 1:prm.I, m in 1:prm.M]
    # 企業の立地需要関数
    F_f = [P_f[i,m].*G[m] for i in 1:prm.I, m in 1:prm.M]

    # Universal Gravity型を使って立地を更新 0708
    h_n1 = F_h
    g_n1 = F_f

    # りゃぷのふ関数の値を表示してみる

    # # 家計の期待最大効用の計算
    # V = (1/θ)*log(sum(exp(θ*(u_d[i])) for i in 1:prm.I))
    # # 企業の期待最大利潤の計算
    # Ψ = [1/η[m].*log(sum(exp(η[m].*(y_d[i,m])) for i in 1:prm.I)) for m in 1:prm.M]

    # # logit 式を用いて立地を更新 0702修正
    # h_n1 = [h_n[i] + δ * ((H * exp(θ * (u_d[i])) / exp(θ * V)) - h_n[i]) for i in 1:prm.I]
    # g_n1 = [g_n[i,m] + δ * ((G[m] * exp(η[m] * (y_d[i,m])) / exp(η[m] * Ψ[m])) - g_n[i,m]) for i in 1:prm.I, m in 1:prm.M]
    

    # 絶対誤差のチェック
    max_diff_h = maximum(abs.((h_n1 ./ H) .- (h_n ./ H)))

    # g_n1 と g_n を G で除した値の差の最大値を計算
    g_n1_normalized = g_n1 ./ reshape(G, 1, prm.M)  # 「G」を正しい形にブロードキャスト
    g_n_normalized = g_n ./ reshape(G, 1, prm.M)
    max_diff_g = maximum(abs.(g_n1_normalized .- g_n_normalized))

    if max_diff_h < ε && max_diff_g < ε
        println("収束条件を満たしました。繰り返しを終了します。")
        break
    end

    # 更新
    h_n = h_n1
    g_n = g_n1

    # 結果表示
    println("===== 繰り返し $iter 回目 =====")
    println("") # 空白の行を追加
    println("■ 家計の立地パターン h (繰り返し $iter 回目)")
    display(h_n1)
    println("")
    println("■ 企業の立地パターン g (繰り返し $iter 回目)")
    display(g_n1)
    println("") 
    
end#ここのエンドもなんだ？


# 0705
# 少数点以下7桁（家計），6桁（企業）の一致を確認 

# Iteration 83: LOCALLY_SOLVED
# ===== 繰り返し 83 回目 =====

# ■ 家計の立地パターン h (繰り返し 83 回目)
# 3-element Vector{Float64}:
#  0.43883121458301166
#  0.43442145041341945
#  0.12674733500356952

# ■ 企業の立地パターン g (繰り返し 83 回目)
# 3×2 Matrix{Float64}:
#  0.237986  0.148759
#  0.237458  0.213038
#  0.124556  0.0382032

# Iteration 84: LOCALLY_SOLVED
# 収束条件を満たしました。繰り返しを終了します。

# Iteration 85: LOCALLY_SOLVED
# 収束条件を満たしました。繰り返しを終了します。

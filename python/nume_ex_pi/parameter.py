import lattice

# Number of colmns
Num_Cols = 10

# Scaling * tが格子の最小距離となる。
Scaling = 0.30

city_network = lattice.make_lattice(Num_Cols)
K = Num_Cols * Num_Cols

# 通勤費用パラメータ
t = 0.1

# 距離抵抗パラメータ
tau = 0.01

# 総土地供給量
S_total = 100
S_bar = S_total / K

Coordinate_Data = \
np.array([(city_network['node_dic'][str(i)]['x_pos']*Scaling,
           city_network['node_dic'][str(i)]['y_pos']*Scaling) for i in range(K)])
distance_matrix = distance.squareform(distance.pdist(Coordinate_Data))

theta_N = 1.0
theta_M = 1.0

E = 5

alpha_1 = 0.4
alpha_2 = 0.4
beta_1 = 0.4
beta_2 = 0.4

M = 10
N = 10
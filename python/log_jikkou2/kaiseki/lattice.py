"""Module summary.

This is lattice.py.
"""
import numpy as np

#境界あり
def create_column(n, tr_row_id):
    link_dic = {}
    node_dic = {}
    for i in range(n):
        node_id = str(tr_row_id*n + i)
        node_dic[node_id] = {}
        node_dic[node_id]['x_pos'] = tr_row_id
        node_dic[node_id]['y_pos'] = i
        
        next_node_id = str(int(node_id)+1)
        if i < n-1:
            
            link1_id = node_id + '_' + next_node_id
            link2_id = next_node_id + '_' + node_id
            
            link_dic[link1_id] = {}
            link_dic[link1_id]['from_node_id'] = node_id
            link_dic[link1_id]['to_node_id'] = next_node_id

            link_dic[link2_id] = {}
            link_dic[link2_id]['from_node_id'] = next_node_id
            link_dic[link2_id]['to_node_id'] = node_id
            
    return link_dic, node_dic

# add column function
def add_column(n, network_dic, tr_row_id, new_link_dic, new_node_dic):
    
    #add new link_dic, node_dic
    network_dic['link_dic'].update(new_link_dic)
    network_dic['node_dic'].update(new_node_dic)
    
    # connect
    for node_id in new_node_dic.keys():
        left_node_id = str(int(node_id) - n)
        
        link1_id = node_id + '_' + left_node_id
        link2_id = left_node_id + '_' + node_id
            
        network_dic['link_dic'][link1_id] = {}
        network_dic['link_dic'][link1_id]['from_node_id'] = node_id
        network_dic['link_dic'][link1_id]['to_node_id'] = left_node_id

        network_dic['link_dic'][link2_id] = {}
        network_dic['link_dic'][link2_id]['from_node_id'] = left_node_id
        network_dic['link_dic'][link2_id]['to_node_id'] = node_id
        
    return network_dic

# make lattice network function
def make_lattice(n):
    network_dic = {}
    network_dic['link_dic'], network_dic['node_dic'] = create_column(n, 1)
    
    for row_id in range(n):
        new_link_dic, new_node_dic = create_column(n, row_id)
        network_dic = add_column(n, network_dic, row_id, new_link_dic, new_node_dic)
            
    return network_dic

#境界なし
class LatticeCalculator:
    def __init__(self, Num_Cols, Scaling, n):
        self.Num_Cols = Num_Cols
        self.Scaling = Scaling
        self.K = Num_Cols * Num_Cols
        self.from_node = n
        self.distance = np.zeros(Num_Cols * Num_Cols)
        self.node_index = np.zeros(Num_Cols * Num_Cols, dtype=int)

    def calc_distance(self, i, j):
        
        #左上から一行ずつ番号付けをしていく．pythonではノード番号が0から始まることに注意．
        dx = abs(i % self.Num_Cols - j % self.Num_Cols) #余りを表す
        dy = abs(i // self.Num_Cols - j // self.Num_Cols) #商を表す
        
        #都市全体の半分より離れていたら，境界を経由した方が早い
        if dx > (self.Num_Cols-1)/2:
            dx = (self.Num_Cols-1) - dx
        else:
            dx = dx
            
        #都市全体の半分より離れていたら，境界を経由した方が早い  
        if dy > (self.Num_Cols-1)/2:
            dy = (self.Num_Cols-1) - dy
        else:
            dy = dy
        
        dx = dx * self.Scaling
        dy = dy * self.Scaling

        return np.sqrt((dx*dx + dy*dy))

    def set(self):
        for i in range(self.K):
            self.node_index[i] = i
            self.distance[i] = self.calc_distance(self.from_node, i)
        # print("node_index", self.node_index)
        # print("distance", self.distance)
    
class MakePeriodicDistanceMatrix:
    def __init__(self, Num_Cols, Scaling):
        self.Num_Cols = Num_Cols
        self.Scaling = Scaling
        self.K = Num_Cols * Num_Cols
        self.distance_matrix = np.zeros((Num_Cols * Num_Cols, Num_Cols * Num_Cols))

    def add_distance_row(self):
        for n in range(self.K):
            lattice_calculator = LatticeCalculator(self.Num_Cols, self.Scaling, n)
            lattice_calculator.set()
            self.distance_matrix[n, :] = lattice_calculator.distance

# make_periodic_distance_matrix = MakePeriodicDistanceMatrix(Num_Cols)
# make_periodic_distance_matrix.add_distance_row()
# print("distance_matrix:")
# print(make_periodic_distance_matrix.distance_matrix)
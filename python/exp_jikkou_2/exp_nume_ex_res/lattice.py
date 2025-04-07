"""Module summary.

This is lattice.py.
"""

# create column function
def create_column(Col, tr_row_id):
    link_dic = {}
    node_dic = {}
    for i in range(Col):
        node_id = str(tr_row_id * Col + i)
        node_dic[node_id] = {}
        node_dic[node_id]['x_pos'] = tr_row_id
        node_dic[node_id]['y_pos'] = i
        
        next_node_id = str(int(node_id)+1)
        if i < Col - 1:
            
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
def add_column(Col, network_dic, tr_row_id, new_link_dic, new_node_dic):
    
    #add new link_dic, node_dic
    network_dic['link_dic'].update(new_link_dic)
    network_dic['node_dic'].update(new_node_dic)
    
    # connect
    for node_id in new_node_dic.keys():
        left_node_id = str(int(node_id) - Col)
        
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
def make_lattice(Col):
    network_dic = {}
    network_dic['link_dic'], network_dic['node_dic'] = create_column(Col, 1)
    
    for row_id in range(Col):
        new_link_dic, new_node_dic = create_column(Col, row_id)
        network_dic = add_column(Col, network_dic, row_id, new_link_dic, new_node_dic)
            
    return network_dic
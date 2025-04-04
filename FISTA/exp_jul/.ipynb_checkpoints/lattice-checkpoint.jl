module LatticeNetwork

export create_column, add_column, make_lattice

using Random
using Printf

# 関数 create_column
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

# 関数 add_column
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

# 関数 make_lattice
function make_lattice(Col::Int)
    network_dic = Dict{String, Dict}()
    network_dic["link_dic"], network_dic["node_dic"] = create_column(Col, 1)

    for row_id in 0:Col-1
        new_link_dic, new_node_dic = create_column(Col, row_id)
        network_dic = add_column(Col, network_dic, row_id, new_link_dic, new_node_dic)
    end
    return network_dic
end

end # module LatticeNetwork
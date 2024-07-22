#Module 3 Reaction networks analysis
## 3.2 反应网络分析
### 3.2.1 物种状态转变检索

import os
import numpy as np
import pandas as pd
from functools import lru_cache
import networkx as nx
import matplotlib.pyplot as plt
from config import work_path, bondfile_name, N_Frame, type_dic


@lru_cache(maxsize=None)
def create_reaction_matrixs():
    if os.path.exists(work_path+'state_trace_data_filter.csv'):
        
        # Load data
        data = pd.read_csv(work_path+'state_trace_data_filter.csv')
        unique_states = np.unique(data.astype(str).values)
        state_index = {state: i for i, state in enumerate(unique_states)}

        def create_transition_prob_matrix(data,unique_states,state_index): #记录某种状态转变为另一状态的次数，在某种状态转变成其它所有状态次数的百分比的矩阵

            # Initialize transition matrix
            num_states = len(unique_states)
            transition_matrix = np.zeros((num_states, num_states))

            # Calculate transition counts without counting transitions to the same state  不统计转变成自身相同转变的转变类型
            for row in data.values:
                for i in range(len(row) - 1):
                    if pd.isna(row[i]) or pd.isna(row[i + 1]):
                        continue
                    if row[i] == row[i + 1]:  # Skip transitions to the same state
                        continue
                    start_state = state_index[str(row[i])]
                    end_state = state_index[str(row[i + 1])]
                    transition_matrix[start_state, end_state] += 1  #transition_matrix内容解读，第5行第7列的值，代表索引为5的分子向索引为7的分子的转变的次数

            # Convert counts to probabilities
            row_sums = transition_matrix.sum(axis=1)
            for i, (row, total) in enumerate(zip(transition_matrix, row_sums)):
                if total != 0:
                    transition_matrix[i] = row / total

            return transition_matrix, unique_states

        
        def create_transition_count_matrix(data,unique_states,state_index): #记录某种状态转变在总共发生过多少次的矩阵
            """
            Create a matrix that counts the transitions between states in the data.
            
            Parameters:
            - data: DataFrame with rows representing entities and columns representing sequential states.

            Returns:
            - transition_count_matrix: Matrix with counts of transitions between states.
            - unique_states: List of unique states found in the data.
            """

            # Initialize transition count matrix
            num_states = len(unique_states)
            transition_count_matrix = np.zeros((num_states, num_states))

            # Calculate transition counts without counting transitions to the same state
            for row in data.values:
                for i in range(len(row) - 1):
                    if pd.isna(row[i]) or pd.isna(row[i + 1]):
                        continue
                    if row[i] == row[i + 1]:  # Skip transitions to the same state
                        continue
                    start_state = state_index[str(row[i])]
                    end_state = state_index[str(row[i + 1])]
                    transition_count_matrix[start_state, end_state] += 1

            return transition_count_matrix

        
        def create_transition_path_count_matrix(data,unique_states,state_index): #记录某种状态转变在多少条目标原子路径中发生过的矩阵
            """
            Create a matrix that counts the number of unique paths between states in the data.
            
            Parameters:
            - data: DataFrame with rows representing entities and columns representing sequential states.
            - unique_states: List of unique states found in the data.

            Returns:
            - transition_path_count_matrix: Matrix with counts of unique paths between states.
            """

            # Initialize transition path count matrix
            num_states = len(unique_states)
            transition_path_count_matrix = np.zeros((num_states, num_states))

            for row in data.values:
                path_recorded = set()  # To ensure we count each path only once for a specific transition
                for i in range(len(row) - 1):
                    if row[i] == row[i + 1]:  # Skip transitions to the same state
                        continue
                    start_state_idx = state_index[str(row[i])]
                    end_state_idx = state_index[str(row[i + 1])]
                    
                    # Check if the path has been counted for this transition
                    if (start_state_idx, end_state_idx) not in path_recorded:
                        transition_path_count_matrix[start_state_idx, end_state_idx] += 1
                        path_recorded.add((start_state_idx, end_state_idx))

            return transition_path_count_matrix
        
        def compute_net_accumulation(data,unique_states): #用于计算每个物种的净积累量，即转变成该物种的次数减去该物种转变成其它物种的次数，结果以字典的形式存在
            # Initialize dictionaries to store the counts
            N1_counts = {state: 0 for state in unique_states}
            N2_counts = {state: 0 for state in unique_states}
            
            # Count the number of transitions to each state (N1)
            for row in data.values:
                for i in range(1, len(row)):
                    N1_counts[row[i]] += 1
                    
            # Count the number of transitions from each state (N2)
            for row in data.values:
                for i in range(len(row) - 1):
                    N2_counts[row[i]] += 1
                    
            # Calculate the net transitions for each state
            net_transitions = {state: N1_counts[state] - N2_counts[state] for state in unique_states} #该字典中存放了state_trace_data_clean.csv中涉及到的所有分子的净积累量
            
            return net_transitions

        transition_prob_matrix, unique_states = create_transition_prob_matrix(data,unique_states,state_index)
        transition_count_matrix = create_transition_count_matrix(data,unique_states,state_index)
        transition_path_count_matrix = create_transition_path_count_matrix(data,unique_states,state_index)
        Target_species_net_accumulation_dic = compute_net_accumulation(data,unique_states)
        return unique_states,state_index, transition_prob_matrix, transition_count_matrix, transition_path_count_matrix, Target_species_net_accumulation_dic #一共5个输出

    else:
        print("Please run Module 3.1 to obatin state_trace_data_clean.csv")
        return None, None, None, None, None, None





def explore_next_states_info(input_state):
    # 步骤3: 调用缓存版本的create_reaction_matrixs
    matrices = create_reaction_matrixs()
    if matrices[0] is None:  # 检查create_reaction_matrixs是否成功返回数据
        print("Data not loaded properly.")
        return pd.DataFrame()  # 返回一个空的DataFrame

    unique_states,state_index, transition_prob_matrix, transition_count_matrix, transition_path_count_matrix, Target_species_net_accumulation_dic = create_reaction_matrixs()

    input_state_idx = state_index[input_state]
    
    # Fetch the transition probabilities for the input state
    transition_probs = transition_prob_matrix[input_state_idx]
    
    # Fetch the transition counts and path counts for the input state
    transition_counts_for_input = transition_count_matrix[input_state_idx]
    path_counts_for_input = transition_path_count_matrix[input_state_idx]
    
    # Compute net transitions (transitions in - transitions out)
    net_transitions = transition_counts_for_input - transition_count_matrix[:, input_state_idx]
    
    # Prepare data for DataFrame
    data_for_df = []
    sorted_indices = np.argsort(transition_probs)[::-1]
    for i in sorted_indices:
        if path_counts_for_input[i] > 0:
            data_for_df.append({
                "Next state name": unique_states[i],
                "probability of trans": transition_probs[i],
                "no. of net trans": net_transitions[i],
                "no. of times of reaction": transition_counts_for_input[i],
                "no. of path": path_counts_for_input[i]
            })
    df = pd.DataFrame(data_for_df)
    return df


def explore_previous_states_info(input_state):
    # 步骤3: 调用缓存版本的create_reaction_matrixs
    matrices = create_reaction_matrixs()
    if matrices[0] is None:  # 检查create_reaction_matrixs是否成功返回数据
        print("Data not loaded properly.")
        return pd.DataFrame()  # 返回一个空的DataFrame

    unique_states,state_index, transition_prob_matrix, transition_count_matrix, transition_path_count_matrix, Target_species_net_accumulation_dic = create_reaction_matrixs()
    # Convert the input state to index
    input_state_idx = state_index[input_state]
    

    # Fetch the transition probabilities to the input state
    transition_probs_to_input = transition_prob_matrix[:, input_state_idx]
    
    # Fetch the transition counts and path counts to the input state
    transition_counts_to_input = transition_count_matrix[:, input_state_idx]
    path_counts_to_input = transition_path_count_matrix[:, input_state_idx]
    
    # Compute net transitions (transitions in - transitions out)
    net_transitions = transition_counts_to_input - transition_count_matrix[input_state_idx, :]
    
    # Prepare data for DataFrame
    data_for_df = []
    sorted_indices = np.argsort(transition_probs_to_input)[::-1]
    for i in sorted_indices:
        if path_counts_to_input[i] > 0:
            data_for_df.append({
                "Previous state name": unique_states[i],
                "probability of trans": transition_probs_to_input[i],
                "no. of net trans": net_transitions[i],
                "no. of times of reaction": transition_counts_to_input[i],
                "no. of path": path_counts_to_input[i]
            })

    # Convert the list of dictionaries to a DataFrame
    df = pd.DataFrame(data_for_df)
    return df



#### 3.2.2 物种反应网络绘制
def atuo_plot_reaction_networks(top_n,plot_graph=False):
    # 步骤3: 调用缓存版本的create_reaction_matrixs
    matrices = create_reaction_matrixs()
    if matrices[0] is None:  # 检查create_reaction_matrixs是否成功返回数据
        print("Data not loaded properly.")
        return pd.DataFrame()  # 返回一个空的DataFrame

    unique_states,state_index, transition_prob_matrix, transition_count_matrix, transition_path_count_matrix, Target_species_net_accumulation_dic = create_reaction_matrixs()
#     """ 
#     1. 区分不同净反应次数的物种，根据Target_species_net_accumulation_dic，净积累量小于零，等于0和大于0三种，用三种颜色区分
#     2.考虑逆向转变问题
#     3.搜索功能分为：
#     a.自动搜索，即以净累计量小于0的分子，且其绝对值排名前5为起点，且搜索到净积累量排名前5的分子为终点，遇到则停止搜索（该功能适合于那些对反应物种没有了解的情况下使用，不建议使用）；
#     b. 手动输入搜索起点和终点.（该功能适合对物种信息有所了解的前提下使用，建议先使用前面的模块对物种信息进行检索）
#     """
#     matrices = create_reaction_matrixs()
#     if matrices[0] is None:  # 检查create_reaction_matrixs是否成功返回数据
#         print("Data not loaded properly.")
#         return pd.DataFrame()  # 返回一个空的DataFrame
#     unique_states, transition_prob_matrix, transition_count_matrix, transition_path_count_matrix, Target_species_net_accumulation_dic = create_reaction_matrixs()



    def build_state_transition_dict(top_n, transition_prob_matrix, transition_count_matrix, unique_states,state_index,Target_species_net_accumulation_dic):
        # 内部函数：创建净转变次数矩阵
        def create_net_transition_matrix(transition_count_matrix):
            num_states = transition_count_matrix.shape[0]
            net_transition_matrix = np.zeros((num_states, num_states))

            for i in range(num_states):
                for j in range(num_states):
                    net_transition_matrix[i, j] = transition_count_matrix[i, j] - transition_count_matrix[j, i]

            return net_transition_matrix

        net_transition_matrix = create_net_transition_matrix(transition_count_matrix)

        Target_species = list(Target_species_net_accumulation_dic.keys())

        state_transition_dict = {}
        for state in Target_species:
            current_index = state_index[state]

            # 根据转移概率矩阵寻找下一状态
            next_indices = np.argsort(-transition_prob_matrix[current_index])[:top_n]
            
            # 为每个下一状态收集转变概率和净转变次数
            next_states_info = []
            for i in next_indices:
                if net_transition_matrix[current_index, i] > 0: # 过滤掉净转变次数为负数的转变
                    transition_prob = transition_prob_matrix[current_index, i]
                    net_transition_count = net_transition_matrix[current_index, i]
                    next_states_info.append((unique_states[i], transition_prob, net_transition_count))

            state_transition_dict[state] = next_states_info

        return state_transition_dict

    state_transition_dict = build_state_transition_dict(top_n, transition_prob_matrix, transition_count_matrix, unique_states,state_index,Target_species_net_accumulation_dic)

    def is_on_single_path(G, node):
        """
        Check if the node is on a single path, which means:
        - It has exactly one predecessor, and
        - Its predecessor has exactly one successor (the current node itself).
        """
        predecessors = list(G.predecessors(node))
        if len(predecessors) != 1:
            return False
        if len(list(G.successors(predecessors[0]))) != 1:
            return False
        return True

    def filter_paths(G, state_index, Target_species_net_accumulation_dic, unique_states):
        while True:
            # 标记是否进行了节点移除
            removed_node = False

            # 第一步：过滤游离节点和小集群
            wcc = list(nx.connected_components(G.to_undirected()))
            largest_wcc = max(wcc, key=len)
            nodes_to_remove = {node for wcc in wcc for node in wcc if wcc != largest_wcc}
            if nodes_to_remove:
                G.remove_nodes_from(nodes_to_remove)
                removed_node = True

            # 第二步：过滤掉以净积累量小于等于0的物种为终端的分支路径
            end_nodes = [node for node in G.nodes() if G.out_degree(node) == 0]
            for node in end_nodes:
                state = unique_states[node]
                net_acc = Target_species_net_accumulation_dic[state]
                if net_acc <= 0:
                    current_node = node
                    while is_on_single_path(G, current_node):
                        predecessor = next(G.predecessors(current_node), None)
                        if predecessor is None or G.in_degree(current_node) > 1:
                            break
                        G.remove_node(current_node)
                        removed_node = True
                        current_node = predecessor

            # 第三步：过滤掉以净积累量大于等于0的物种为起始端的分支路径
            start_nodes = [node for node in G.nodes() if G.in_degree(node) == 0]
            for node in start_nodes:
                state = unique_states[node]
                net_acc = Target_species_net_accumulation_dic[state]
                if net_acc >= 0:
                    G.remove_node(node)
                    removed_node = True
                    continue
                    current_node = node
                    successors = list(G.successors(current_node)) if node in G else []
                    while successors:
                        if len(successors) != 1 or G.in_degree(successors[0]) != 1:
                            break
                        next_node = successors[0]
                        G.remove_node(next_node)
                        if next_node not in G:
                            break
                        current_node = next_node
                        successors = list(G.successors(current_node)) if current_node in G else []

            # 如果这一轮没有移除任何节点，则结束过滤
            if not removed_node:
                break


        

    def read_species_time_centroid(filename):
        """从Excel文件读取物种的时间权重"""
        df = pd.read_excel(filename)
        # 假设文件中有'Species_name'和'Time Centroid'列
        time_centroids = df.set_index('Species_name')['Time Centroid'].to_dict()
        return time_centroids
                
    time_centroids_dic = read_species_time_centroid(work_path+"Species_statistical_info_filter.xlsx"  )             
                
    def visualize_updated_graph(state_transition_dict, unique_states, state_index, Target_species_net_accumulation_dic):
        G = nx.DiGraph()
        
        node_color_map_dic = {}
        for state in unique_states:
            net_acc = Target_species_net_accumulation_dic[state]
            node_color_map_dic[state] = 'red' if net_acc < 0 else ('green' if net_acc > 0 else 'lightblue')

        for state in unique_states:
            G.add_node(state_index[state])

        for state, next_states_info in state_transition_dict.items():
            current_index = state_index[state]
            for next_state_info in next_states_info:
                next_state, _, _ = next_state_info
                next_index = state_index[next_state]
                G.add_edge(current_index, next_index)

        # 执行过滤操作
        filter_paths(G, state_index, Target_species_net_accumulation_dic, unique_states)
        
        
        node_colors = [node_color_map_dic[unique_states[n]] for n in G.nodes()]
        
        #分子式形式的label
        labels_dic = {}
        for n in G.nodes():
            labels_dic[n] = unique_states[n]
            
        #按照Time_centord权重对对节点生成新的time_centord_label
        def read_species_time_centroid(filename):
            """从Excel文件读取物种的时间权重"""
            df = pd.read_excel(filename)
            # 假设文件中有'Species_name'和'Time Centroid'列
            time_centroids = df.set_index('Species_name')['Time Centroid'].to_dict()
            return time_centroids

        time_centroids_dic = read_species_time_centroid(work_path+"Species_statistical_info_filter.xlsx"  ) 

        G_node_unique_states_list = [unique_states[n] for n in G.nodes()]
        G_node_order_list = [n for n in G.nodes()]
        time_centord_value = [] 
        for i in G_node_unique_states_list:
            time_centord_value.append(time_centroids_dic[i])

        time_centord_label = np.argsort(np.argsort(time_centord_value)) + 1
        time_centord_label = time_centord_label.tolist()
        time_centord_label_dic = dict(zip(G_node_order_list, time_centord_label))
        
        time_centord_label_unique_states_dic = dict(zip(time_centord_label,G_node_unique_states_list))
        time_centord_label_unique_states_df = pd.DataFrame(list(time_centord_label_unique_states_dic.items()), columns=['Order', 'Species'])
        # 根据Order列排序
        time_centord_label_unique_states_df = time_centord_label_unique_states_df.sort_values(by='Order')
        if plot_graph:
            print(time_centord_label_unique_states_df.to_string(index = False))        
            plt.figure(figsize=(12, 12))
            pos = nx.kamada_kawai_layout(G)
            nx.draw(G,pos,node_color=node_colors, with_labels=True,node_size=200, edge_color='gray', labels=time_centord_label_dic  ,font_size=8)#labels=labels_dic,或者labels=time_centord_label_dic 
            plt.show()
        return G
    G = visualize_updated_graph(state_transition_dict, unique_states,state_index,Target_species_net_accumulation_dic)
    return G




def manual_plot_reaction_networks(top_n,start_states=[], end_states=[],plot_graph=True):
    G = atuo_plot_reaction_networks(top_n,plot_graph=False)
    # 步骤3: 调用缓存版本的create_reaction_matrixs
    matrices = create_reaction_matrixs()
    if matrices[0] is None:  # 检查create_reaction_matrixs是否成功返回数据
        print("Data not loaded properly.")
        return pd.DataFrame()  # 返回一个空的DataFrame
    unique_states,state_index, transition_prob_matrix, transition_count_matrix, transition_path_count_matrix, Target_species_net_accumulation_dic = create_reaction_matrixs()
    
    # 步骤3: 调用缓存版本的create_reaction_matrixs

    def remove_unwanted_starting_paths(sub_G, start_states, end_states, state_index):
        """
        修改版：删除非设定起点但具有出射没有入射的节点及其分支路径，同时保护设定的终点不被删除。
        :param sub_G: 子图
        :param start_states: 设定的起点状态列表
        :param end_states: 设定的终点状态列表
        :param state_index: 状态到节点索引的映射
        """
        unwanted_starts = [n for n in sub_G.nodes if sub_G.in_degree(n) == 0 and sub_G.out_degree(n) > 0 and unique_states[n] not in start_states]

        for start_node in unwanted_starts:
            nodes_to_remove = set()
            queue = [start_node]
            while queue:
                current_node = queue.pop(0)
                if current_node in nodes_to_remove:
                    continue
                # 如果当前节点是设定的起点或终点之一，停止删除
                if unique_states[current_node] in start_states or unique_states[current_node] in end_states:
                    break
                nodes_to_remove.add(current_node)
                queue.extend(sub_G.successors(current_node))
            
            # 从待删除节点中移除设定的终点
            end_state_indexes = [state_index[state] for state in end_states if state in state_index]
            nodes_to_remove = nodes_to_remove - set(end_state_indexes)

            sub_G.remove_nodes_from(nodes_to_remove)

    def remove_outgoing_paths_from_end_states(sub_G, end_states, state_index):
        """
        删除以end_states为终点的发出的路径。
        :param sub_G: 子图
        :param end_states: 终点状态列表
        :param state_index: 状态到节点索引的映射
        """
        for end_state in end_states:
            end_index = state_index[end_state]
            # 检查终点是否有出度
            while sub_G.out_degree(end_index) > 0:
                # 删除所有发出的边
                out_edges = list(sub_G.out_edges(end_index))
                sub_G.remove_edges_from(out_edges)
                # 如果删除发出的边后，终点成为了孤立节点，则停止处理
                if sub_G.degree(end_index) == 0:
                    break
                # 删除因删除发出的边而变成孤立的节点
                isolated_nodes = [n for n in sub_G.nodes if sub_G.degree(n) == 0 and n != end_index]
                sub_G.remove_nodes_from(isolated_nodes)

    def visualize_reaction_network(G, unique_states, state_index, Target_species_net_accumulation_dic, start_states=[], end_states=[]):
        sub_G = nx.DiGraph()

        # 如果指定了起点
        if start_states:
            start_indices = [state_index[s] for s in start_states if s in state_index]
        else:
            start_indices = list(G.nodes)

        # 如果指定了终点
        if end_states:
            end_indices = [state_index[s] for s in end_states if s in state_index]
        else:
            end_indices = list(G.nodes)

        # 寻找路径
        for start_index in start_indices:
            for end_index in end_indices:
                if start_index != end_index:
                    for path in nx.all_simple_paths(G, source=start_index, target=end_index):
                        nx.add_path(sub_G, path)

        # 删除终点的发出路径
        if end_states:
            remove_outgoing_paths_from_end_states(sub_G, end_states, state_index)

        # 删除非设定起点的节点
        if start_states:
            while True:
                original_node_count = len(sub_G.nodes())
                remove_unwanted_starting_paths(sub_G, start_states, end_states, state_index)
                if len(sub_G.nodes()) == original_node_count:
                    break

        # 根据净积累量给节点上色
        node_color_map = {}
        for state, net_acc in Target_species_net_accumulation_dic.items():
            if net_acc < 0:  # 净积累量小于0，红色
                node_color_map[state_index[state]] = "red"
            elif net_acc > 0:  # 净积累量大于0，绿色
                node_color_map[state_index[state]] = "green"
            else:  # 净积累量等于0，蓝色
                node_color_map[state_index[state]] = "lightblue"

                    
        node_colors = [node_color_map.get(node, "grey") for node in sub_G.nodes()]  # 使用灰色作为默认颜色

        #分子式形式的label
        labels_dic = {}
        for n in sub_G.nodes():
            labels_dic[n] = unique_states[n]
            
        #按照Time_centord权重对对节点生成新的time_centord_label
        def read_species_time_centroid(filename):
            """从Excel文件读取物种的时间权重"""
            df = pd.read_excel(filename)
            # 假设文件中有'Species_name'和'Time Centroid'列
            time_centroids = df.set_index('Species_name')['Time Centroid'].to_dict()
            return time_centroids

        time_centroids_dic = read_species_time_centroid(work_path+"Species_statistical_info_filter.xlsx"  ) 

        sub_G_node_unique_states_list = [unique_states[n] for n in sub_G.nodes()]
        sub_G_node_order_list = [n for n in sub_G.nodes()]
        time_centord_value = [] 
        for i in sub_G_node_unique_states_list:
            time_centord_value.append(time_centroids_dic[i])

        time_centord_label = np.argsort(np.argsort(time_centord_value)) + 1
        time_centord_label = time_centord_label.tolist()
        time_centord_label_dic = dict(zip(sub_G_node_order_list, time_centord_label))
        
        time_centord_label_unique_states_dic = dict(zip(time_centord_label,sub_G_node_unique_states_list))
        time_centord_label_unique_states_df = pd.DataFrame(list(time_centord_label_unique_states_dic.items()), columns=['Order', 'Species'])
        # 根据Order列排序
        time_centord_label_unique_states_df = time_centord_label_unique_states_df.sort_values(by='Order')

        #labels_dic = {n: unique_states[n] for n in sub_G.nodes()}
        if plot_graph:
            print(time_centord_label_unique_states_df.to_string(index = False)) 
            plt.figure(figsize=(12, 12))
            pos = nx.kamada_kawai_layout(sub_G)
            nx.draw(sub_G, pos, node_color=node_colors, with_labels=True, labels=time_centord_label_dic, font_size=12)#time_centord_label_dic
            plt.show()
        return sub_G
    sub_G = visualize_reaction_network(G, unique_states, state_index, Target_species_net_accumulation_dic, start_states, end_states)
    return sub_G



    # sub_G = visualize_reaction_network(G, unique_states, state_index, Target_species_net_accumulation_dic, start_states, end_states)
    # return sub_G


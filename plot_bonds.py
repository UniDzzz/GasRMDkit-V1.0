## 2.3 Plot bonds
###函数输入target形式形如['Mo-S','S-S'......]
from config import work_path, bondfile_name, N_Frame, type_dic
import numpy as np
import matplotlib.pyplot as plt


def bonds_time_evolution(target,list_N_Frame,atom_type_dic):

    def count_bonds_per_frame(list_N_Frame, atom_type_dic, type_dic):
        bond_counts = []
        for frame in list_N_Frame:
            bond_count = {}
            for atom_id, connected_atoms in frame.items():
                atom_type = atom_type_dic[atom_id]
                for connected_atom in connected_atoms:
                    connected_atom_type = atom_type_dic[connected_atom]
                    bond_type = '-'.join(sorted([type_dic[atom_type], type_dic[connected_atom_type]]))

                    if atom_id < connected_atom:  # 避免重复计数
                        if bond_type not in bond_count:
                            bond_count[bond_type] = 0
                        bond_count[bond_type] += 1

            bond_counts.append(bond_count)

        return bond_counts
    bond_counts_per_frame = count_bonds_per_frame(list_N_Frame, atom_type_dic, type_dic) 


    def generate_bond_counts_matrix(bond_counts_per_frame):
        all_bond_types = set()     # 确定所有时间步中出现的所有键类型
        for frame in bond_counts_per_frame:
            all_bond_types.update(frame.keys())
        bond_type_index = {bond_type: i for i, bond_type in enumerate(sorted(all_bond_types))}     # 为每种键类型创建索引
        bond_counts_matrix = np.zeros((len(all_bond_types), len(bond_counts_per_frame)))     # 初始化矩阵
        for frame_idx, frame in enumerate(bond_counts_per_frame):     # 填充矩阵
            for bond_type, count in frame.items():
                row_idx = bond_type_index[bond_type]
                bond_counts_matrix[row_idx, frame_idx] = count
                
        return bond_counts_matrix, bond_type_index

    bond_counts_matrix, bond_type_index = generate_bond_counts_matrix(bond_counts_per_frame)

    for bond_type in target:
        if bond_type in bond_type_index:
            idx = bond_type_index[bond_type]
            plt.plot(bond_counts_matrix[idx, :], label=bond_type)
    plt.xlabel('Time Step')
    plt.ylabel('Bond Count')
    plt.legend()
    plt.show()
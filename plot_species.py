from config import work_path, bondfile_name, N_Frame, type_dic
from species_statistical_info import sepcies_statistical_analysis
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import re
def speices_time_evolution(target,List_N_Molecular):
    """
    该函数预期实现的功能是 
    1.输入任意存在于Speices_name列表中的物种名称，均能绘制出对应的数量随时间变化曲线，比如输入['Species_name1','Species_name2'......]
    2. 输入任意元素名称，如('Mo')，再如('S'),为了区别与功能1，所以功能2这里用小括号；就能绘制所有含有该元素并且存在于Speices_name_filter列表中的物种数量随时间变化的曲线
    3. 输入任意元素或者元素组合，如{'Mo'},{'Mo','S'},{'Mo','O','S'},绘制只含有大括号内元素并且存在于Speices_name_filter列表中的物种数量随时间变化曲线。

    参数:
    - target: 目标物种名称列表、含有特定元素的元组或含特定元素组合的集合。
    - all_molecules_matrix: 所有物种数量随时间变化的矩阵。
    - all_molecule_indices: 物种名称到矩阵行索引的映射。
    - Species_name: 所有物种名称列表。
    - Species_name_filter: 经过筛选的物种名称列表。
    - type_dic: 元素类型的字典，键为元素的类型编号，值为元素符号。
    """
    Species_name_df, Species_name_filter_df = sepcies_statistical_analysis(List_N_Molecular)
    Species_name_list = Species_name_df['Species_name'].tolist()
    Species_name_filter_list = Species_name_filter_df['Species_name'].tolist()
    def generate_all_molecules_matrix(List_N_Molecular):

        # 识别出所有出现过的分子种类
        all_molecules = set()
        for frame in List_N_Molecular:
            all_molecules.update(frame)
        molecules_index = {molecule: i for i, molecule in enumerate(all_molecules)}    # 为每种分子创建索引
        all_molecules_matrix = np.zeros((len(all_molecules), len(List_N_Molecular)))    # 初始化矩阵
        for frame_idx, frame in enumerate(List_N_Molecular):    # 填充矩阵
            for molecule in frame:
                all_molecules_matrix[molecules_index[molecule], frame_idx] += 1

        return all_molecules_matrix, molecules_index
    all_molecules_matrix, all_molecule_indices = generate_all_molecules_matrix(List_N_Molecular)

    def species_contains_elements(species, elements):
        # 获取物种中所有元素的正则表达式模式
        pattern = r'[A-Z][a-z]*'
        species_elements = set(re.findall(pattern, species))
        return species_elements == elements

    if isinstance(target, list):  # 输入为物种名称列表
        for species in target:
            if species in Species_name_list:
                idx = all_molecule_indices[species]
                plt.plot(all_molecules_matrix[idx, :], label=species)
    elif isinstance(target, tuple):  # 输入为含特定元素的元组
        element = target[0]
        for species in Species_name_filter_list:
            if element in species:
                idx = all_molecule_indices[species]
                plt.plot(all_molecules_matrix[idx, :], label=species)
    elif isinstance(target, set):  # 输入为含特定元素组合的集合
        elements = {type_dic.get(e, e) for e in target}  # 从type_dic转换元素名称
        for species in Species_name_filter_list:
            if species_contains_elements(species, elements):
                idx = all_molecule_indices[species]
                plt.plot(all_molecules_matrix[idx, :], label=species)
    
    plt.xlabel('Time')
    plt.ylabel('Quantity')
    plt.legend()
    plt.show()
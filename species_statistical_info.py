# Module 2. Species & bond analysis
## 2.1 Species statistical analysis
###生成模拟过程中所有物种的数量随时间变化的矩阵 
import pandas as pd
import numpy as np
import re
import os
from config import work_path, bondfile_name, N_Frame, type_dic
from scipy.integrate import trapz, simps
from scipy.ndimage import gaussian_filter1d
def sepcies_statistical_analysis(List_N_Molecular):
    def generate_all_molecules_matrix(List_N_Molecular):
        """
        生成包含所有分子类型随时间变化的数量矩阵。

        参数:
        - List_N_Molecular: 每一帧中分子类型的嵌套列表。

        返回:
        - all_molecules_matrix: 每种分子随时间的数量变化矩阵。
        - molecules_index: 分子类型到矩阵行索引的映射。
        """
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

    ###对all_molecules_matrix进行平滑处理
    #使用高斯平滑，助于平滑短期波动，突出长期趋势。
    def gaussian_smoothing(matrix, sigma=2):
        smoothed_matrix = np.zeros_like(matrix)
        for row_idx, row in enumerate(matrix):
            smoothed_matrix[row_idx] = gaussian_filter1d(row, sigma=sigma)
        return smoothed_matrix

    #将小于1的值全部处理为0
    def threshold_filter(matrix, threshold):
        # 使用NumPy的where函数将低于阈值的数据点设置为0
        filtered_matrix = np.where(matrix < threshold, 0, matrix)
        return filtered_matrix


    ###定义Species_statistical_info表格
    #定义分子类别
    def assign_group(molecule):
        """
        根据分子的组成，将其归类为包含不同元素组合的类型。

        参数:
        - molecule: 分子式字符串。
        - type_dic: 元素类型的字典，键为元素的类型编号，值为元素符号。

        返回:
        - 分类后的分子类型字符串。
        """
        # 生成元素符号到类型编号的反向映射
        symbol_to_type = {v: k for k, v in type_dic.items()}

        # 识别分子中包含的所有元素类型
        elements = set(re.findall(r'[A-Z][a-z]*', molecule))

        # 生成分子类型字符串
        molecule_type = []
        for symbol in sorted(elements, key=lambda x: symbol_to_type[x]):  # 按type_dic中的顺序排序
            count = len(re.findall(symbol, molecule))
            molecule_type.append(f"{symbol}{count if count > 1 else ''}")

        return ''.join(molecule_type)
    #定义分子时间重心，对数量分布曲线对时间积分，求得面积为一半处对应的时间
    def compute_time_centroids(matrix):
        centroids = []
        for row in matrix:
            total_quantity = np.sum(row)
            half_quantity = total_quantity / 2
            accumulated_quantity = 0
            centroid_time = 0
            for idx, quantity in enumerate(row):
                accumulated_quantity += quantity
                if accumulated_quantity >= half_quantity:
                    centroid_time = idx
                    break
            centroids.append(centroid_time)
        return centroids
    #定义分子数量积分稳定性权重，即数量随时间演变曲线的面积
    def compute_quantity_integration_stability(matrix):
        quantity_integration_stability = []
        for row in matrix:
            y = [i for i in row]
            x = np.linspace(1, len(y), len(y))  
            # 使用辛普森规则进行积分
            area_simps = simps(y, x)
            quantity_integration_stability.append(area_simps)
        return quantity_integration_stability

    #定义分子数量均值稳定性权重，即数量随时间演变曲线的均值
    def compute_quantity_average_stability(matrix):
        quantity_average_stability = []
        for row in matrix:
            a = sum(row)
            b= np.count_nonzero(row)
            if b == 0:
                quantity_average_stability.append(0)
            else:
                quantity_average_stability.append(a/b)
        return quantity_average_stability

    # 创建Species_statistical_info表格
    def Species_statistical_info_df(matrix):
        df = pd.DataFrame({
            'Species_name':list(all_molecule_indices.keys()),
            'Integration_stability':compute_quantity_integration_stability(matrix),
            'Average_stability': compute_quantity_average_stability(matrix),
            'Time Centroid':compute_time_centroids(matrix),
            'Group': [assign_group(i) for i in list(all_molecule_indices.keys())]
        })

        return df
    #创建过滤Species_statistical_info的函数，选择稳定性占比为设定阈值以内的分子，该功能选择性使用。
    def filter_top_percent_within_group(df, top_n_percent):
        """
        对于每种“Group”中，该Group中的分子的Integration_stability权重从大到小排序，
        保留前top_n_percent个数的分子，向下取整。

        参数:
        - df: 包含分子统计信息的DataFrame。
        - top_n_percent: 保留每个组顶部百分比的分子，值范围为0到1。

        返回:
        - 筛选后的DataFrame。
        """
        def filter_group(group_df):
            # 计算每个组内应保留的分子数量,个数向下取整，抹零
            n_to_keep = int(len(group_df) * top_n_percent)
            # 如果n_to_keep为0但组内有数据，保留至少一个分子
            #n_to_keep = max(n_to_keep, 1)
            # 按Integration_stability降序排序并保留前n_to_keep个分子
            return group_df.sort_values(by="Integration_stability", ascending=False).head(n_to_keep)

        # 对DataFrame按Group进行分组，并应用筛选规则
        filtered_df = df.groupby('Group').apply(filter_group).reset_index(drop=True)
        return filtered_df

    # 使用函数处理List_N_Molecular
    all_molecules_matrix, all_molecule_indices = generate_all_molecules_matrix(List_N_Molecular)############################################
    gs_mo_matrix = gaussian_smoothing(all_molecules_matrix, sigma=2)#################
    all_molecules_matrix_smooth = threshold_filter(gs_mo_matrix, threshold=1)########################################################
    #生成Species_statistical_info_df,并进行排序处理,输出Species_statistical_info.xlsx文件
    species_statistical_info_df = Species_statistical_info_df(all_molecules_matrix_smooth)
    species_statistical_info_df = species_statistical_info_df.sort_values(by=["Group", "Average_stability"], ascending=[True, False]).reset_index(drop=True)
    species_statistical_info_df_filter = species_statistical_info_df[species_statistical_info_df['Average_stability'] > 2.11]
    species_statistical_info_df_filter = filter_top_percent_within_group(species_statistical_info_df_filter,1) #该df中存放了所有物种过滤后的信息，包括所有元素
    if os.path.exists(work_path+'Species_statistical_info.xlsx'):
        print("Hello, Species_statistical_info.xlsx have existed")
    else:
        species_statistical_info_df.to_excel(work_path+'Species_statistical_info.xlsx',index=False)
    if os.path.exists(work_path+'Species_statistical_info_filter.xlsx'):
        print("Hello, Species_statistical_info_filter.xlsx have existed")
    else:
        species_statistical_info_df_filter.to_excel(work_path+'Species_statistical_info_filter.xlsx',index=False)
    return species_statistical_info_df, species_statistical_info_df_filter
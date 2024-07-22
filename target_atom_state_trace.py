#Module 3 Reaction networks analysis
## 3.1 输出目标元素的状态跟踪文件：state_trace_data.csv（去除噪声分子版本）,反应网络基于该文件进行分析
###处理Dic_N_cluster,函数输入element_type = 1 ，该数值为type_dic中的元素对应的数字
import os
import numpy as np
import pandas as pd
from tqdm import tqdm  # 引入tqdm库
from config import work_path, bondfile_name, N_Frame, type_dic
from species_statistical_info import sepcies_statistical_analysis

element_type = 1 #type_dic中的元素对应的数字

def capture_target_atom_track_notclean(element_type,atom_type_dic,Dic_N_cluster):

    def capture_target_atom_id(element_type,atom_type_dic):
        Goal_atom_list = []
        for x,y in atom_type_dic.items():
            if y == element_type:
                Goal_atom_list.append(x)
        return Goal_atom_list
    target_atom_id_list = capture_target_atom_id(element_type,atom_type_dic)

    # 检查文件是否存在
    if os.path.exists(work_path+'state_trace_data.csv'):
        print("Hello, state_trace_data.csv have existed")
    else:
        total_frame_result = [] #用来存放所有目标原子的所有中间状态，该列表将包含长度为len(target_atom_id_list)个小list,每个小List包含N_Frame个分子式。
        for atom in tqdm(target_atom_id_list, desc="文件正在写入，耐心等待哦"):
        #for atom in target_atom_id_list: # 遍历每个目标原子
            result=[] # 创建一个空的list
            for i in range(len(Dic_N_cluster)): # 遍历每个帧
                molecular = Dic_N_cluster[i][int(atom)]
                result.append(molecular)
            total_frame_result.append(result)

        total_frame_result_df = pd.DataFrame(total_frame_result)
        total_frame_result_df.to_csv(work_path+'state_trace_data.csv',index=False)


def capture_target_atom_track_clean(element_type,atom_type_dic,Dic_N_cluster,List_N_Molecular):#过滤掉Spcies_name_noise中的噪声分子，即如果某状态为噪声分子，则该状态取上一个状态的值，如果第一个状态即为噪声分子，则第一个状态可以填该噪声分子。允许噪声进入表格
    a,b = sepcies_statistical_analysis(List_N_Molecular)
    Spcies_name_filter = b['Species_name'].tolist()

    def capture_target_atom_id(element_type,atom_type_dic):
        Goal_atom_list = []
        for x,y in atom_type_dic.items():
            if y == element_type:
                Goal_atom_list.append(x)
        return Goal_atom_list
    target_atom_id_list = capture_target_atom_id(element_type,atom_type_dic)
    # 检查文件是否存在
    if os.path.exists(work_path+'state_trace_data_clean.csv'):
        print("Hello, state_trace_data_clean.csv have existed")
    else:
        total_frame_result = [] 
        for atom in tqdm(target_atom_id_list, desc="文件正在写入，耐心等待哦"):
        #for atom in target_atom_id_list: # 遍历每个目标原子
            result=[] # 创建一个空的list
            previous_molecular = None
            for frame in Dic_N_cluster: #遍历每个帧
                molecular = frame[int(atom)]
                if molecular in Spcies_name_filter:
                    result.append(molecular)
                    previous_molecular = molecular
                else: ###这个写法是运行存在噪声分子的行，该行后续也可能出现噪声分子；如果不是噪声分子开头，整行就不会存在噪声分子。
                    if previous_molecular is None:
                        result.append(molecular)
                    else:
                        result.append(previous_molecular)
            total_frame_result.append(result)
        
        total_frame_result_df = pd.DataFrame(total_frame_result)
        total_frame_result_df.to_csv(work_path+'state_trace_data_clean.csv',index=False)



def capture_target_atom_track_filter(element_type,atom_type_dic,Dic_N_cluster,List_N_Molecular):#过滤掉Spcies_name_noise中的噪声分子，即如果某状态为噪声分子，则该状态取上一个状态的值，如果第一个状态即为噪声分子，则第一个状态可以填该噪声分子。允许噪声进入表格
    a,b = sepcies_statistical_analysis(List_N_Molecular)
    Spcies_name_filter = b['Species_name'].tolist()

    def capture_target_atom_id(element_type,atom_type_dic):
        Goal_atom_list = []
        for x,y in atom_type_dic.items():
            if y == element_type:
                Goal_atom_list.append(x)
        return Goal_atom_list
    target_atom_id_list = capture_target_atom_id(element_type,atom_type_dic)
    # 检查文件是否存在
    if os.path.exists(work_path+'state_trace_data_filter.csv'):
        print("Hello, state_trace_data_filter.csv have existed")
    else:
        total_frame_result = [] 
        for atom in tqdm(target_atom_id_list, desc="文件正在写入，耐心等待哦"):
            result = []  # 创建一个空的list
            previous_molecular = None
            first_frame = True  # 新增一个变量来标记是否是第一帧
            for frame in Dic_N_cluster:  # 遍历每个帧
                molecular = frame[int(atom)]
                if molecular in Spcies_name_filter:
                    result.append(molecular)
                    previous_molecular = molecular
                    first_frame = False  # 一旦处理了第一帧，就将此变量设为False
                else:
                    if first_frame:  # 如果是第一帧且为噪声分子
                        break  # 直接跳出循环，不再追踪此原子
                    else:
                        if previous_molecular is not None:
                            result.append(previous_molecular)
                        # 如果第一个状态就是噪声，这里不会执行到，因为已经在上面的if first_frame中break了
            if not first_frame:  # 如果处理了至少一帧（即这不是个完全跳过的原子）
                total_frame_result.append(result)
        total_frame_result_df = pd.DataFrame(total_frame_result)
        total_frame_result_df.to_csv(work_path+'state_trace_data_filter.csv',index=False)

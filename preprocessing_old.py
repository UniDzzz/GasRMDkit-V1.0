import os
from tqdm import tqdm  # 引入tqdm库
from config import work_path, bondfile_name, N_Frame, type_dic

### 数据读取 ###
# 读取原子数量
atom_num_path = os.path.join(work_path, bondfile_name)
with open(atom_num_path, "r") as atom_num_file:
    for _ in range(2):
        atom_num_file.readline()
    n_atom = int(atom_num_file.readline().split()[4])

# 读取原子id对应的原子类型
bonds_path = os.path.join(work_path, bondfile_name)
with open(bonds_path, "r") as bonds_file:
    atom_type_dic = {}
    for _ in range(7):  # 跳过前7行
        bonds_file.readline()
    for _ in range(n_atom):
        line = bonds_file.readline().split()
        atom_id = int(line[0])
        atom_type = int(line[1])
        atom_type_dic[atom_id] = atom_type

### 函数定义 ###
#输入原子id lists e.g.[1,3,5,6],对应输出原子类型列表[1,2,2,3]
def identify1(input_list):
    output_list = list(map(lambda x: atom_type_dic[x], input_list))
    return output_list
#输入原子类型lists e.g.[1,2,2,3],输出[Mo,O,O,S] 
def identify2(input_list):
    output_list = list(map(lambda x: type_dic[x],input_list))
    return output_list
#读取输入列表中[Mo,O,O,S],输出分子式为MoO2S,分子式中元素按照type_dic中输入的元素序号排序
def identify3(lists):
    type_lists = []
    for i,j in type_dic.items():
        type_lists.append(j)
    str_o = ''
    for i in type_lists:
        if lists.count(i) == 0:
            str_i=''
        elif lists.count(i) == 1:
            str_i = i #Mo
        elif lists.count(i) > 1:
            str_i = i + str(lists.count(i))
        str_o = str_o + str_i
    return str_o

def preprocess_data():
    ### 主要逻辑 ###
    # 构建每一帧原子间连接字典
    list_N_Frame = [] ####################供模块2调用
    Time_steps = []
    Timestep_path = os.path.join(work_path, bondfile_name)
    list_N_Frame_path =  os.path.join(work_path, bondfile_name)
    with open(list_N_Frame_path, "r") as bonds_file, open(Timestep_path, "r") as timestep_file:
        for i in tqdm(range(N_Frame), desc="正在加载中哦，耐心等待"):  # 使用tqdm添加进度条
            if i == 0:
                pre_line = 7
                timesteps = timestep_file.readline().split()[2]
                for n in range(n_atom+6):
                    timestep_file.readline() #读掉 Timestep 500 这一行下面的6行文字加上n_atom行原子信息。
            else:
                pre_line = 8
                timestep_file.readline()
                timesteps = timestep_file.readline().split()[2]
                for n in range(n_atom+6):
                    timestep_file.readline()
            Time_steps.append(timesteps)
            
            # 读取bonds文件
            for _ in range(pre_line):
                bonds_file.readline()
            dic_per_total = {}
            for _ in range(n_atom):
                line = bonds_file.readline().split()
                lists = [int(i) for i in line if '.' not in i]
                del lists[1:3]
                del lists[-1]
                dic_per_total[lists[0]] = lists[1:]
            list_N_Frame.append(dic_per_total)

    # 构建每一帧分子信息列表，包括团簇信息（每个分子所包含的所有原子id）；分子个数；分子式，并输出产物文件。
    List_N_cluster = []  # 存放所有帧的团簇
    List_N_cluster_Num = []  # 存放所有帧中的分子个数
    List_N_Molecular = []  # 存放所有帧的分子类型

    species_out_path = os.path.join(work_path, "species.out")

    # 团簇识别和处理
    with open(species_out_path, "w") as speicies_out:
        for step, frame in enumerate(list_N_Frame):
            timesteps = Time_steps[step]
            visited = set()  # 存储已处理过的原子
            List_per_clusters = []  # 存储每一帧的团簇信息
            dic_per_frame = {}  # 存储当前帧中每个原子ID对应的分子式信息
            # 深度优先搜索识别团簇
            for atom in frame:
                if atom not in visited:
                    molecule = []  # 当前团簇的原子列表
                    List_per_clusters.append(molecule)
                    stack = [atom]

                    while stack:
                        current_atom = stack.pop()
                        if current_atom not in visited:
                            molecule.append(current_atom)
                            visited.add(current_atom)
                            for neighbor in frame[current_atom]:
                                if neighbor not in visited:
                                    stack.append(neighbor)

            # 团簇类型识别
            List_per_molecular = [identify3(identify2(identify1(cluster))) for cluster in List_per_clusters]


            # 更新结果存储列表
            List_N_cluster.append(List_per_clusters)
            List_N_cluster_Num.append(len(List_per_clusters))
            List_N_Molecular.append(List_per_molecular)

            

                # 输出到species.out文件
            Species = {} 
            for product in List_per_molecular:
                Species[product] = Species.get(product,0)+1
            Species_lists = [product for product in Species.keys()]
            Species_counts = [counts for counts in Species.values()]
            Species_counts_str = list(map(lambda x:str(x), Species_counts)) ###把数值变成字符串，才能用于后面写入文件
            speicies_out.write("#"+" "+"Timestep"+"     "+"No_Specs"+"     "+"No_Cluster"+"     "+" ".join(Species_lists) + "\n")
            speicies_out.write(timesteps +"             " +" " + str(len(Species_lists)) + "             "+str(len(List_per_clusters))+"          " +"  ".join(Species_counts_str) +"\n")  

    return  atom_type_dic,list_N_Frame,List_N_Molecular
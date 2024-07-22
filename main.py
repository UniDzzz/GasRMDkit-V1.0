# main.py
from config import work_path, bondfile_name, N_Frame, type_dic
from preprocessing import *
from species_statistical_info import sepcies_statistical_analysis
from plot_species import speices_time_evolution
from plot_bonds import  bonds_time_evolution
from target_atom_state_trace import capture_target_atom_track_filter
from reaction_network_analysis import explore_previous_states_info, explore_next_states_info, atuo_plot_reaction_networks, manual_plot_reaction_networks
# 预处理，预处理只需要执行一次，应放在main函数外部
atom_type_dic, list_N_Frame, List_N_Molecular, Dic_N_cluster = preprocess_data()

def main():
    while True:
        
        print("\n模块选项：")
        print("1: 物种统计分析")
        print("2: 物种时间进化")
        print("3: 键时间进化")
        print("4: 目标原子状态跟踪")
        print("5: 下一转变状态探索")
        print("6: 上一转变状态探索")
        print("7: 自动绘制反应网络")
        print("8: 手动控制绘制反应网络")
        print("0: 退出")
        
        choice = input("请输入你想运行的模块编号（或输入0退出）: ")
        if choice == "1":
            species_statistical_info_df, species_statistical_info_df_filter = sepcies_statistical_analysis(List_N_Molecular)
            print(species_statistical_info_df_filter)
        elif choice == "2":
            target = input("输入目标分子式，如Mo3O9: ")
            speices_time_evolution([str(target)],List_N_Molecular)
        elif choice == "3":
            target = input("输入目标键连，如Mo-S: ")
            bonds_time_evolution([str(target)],list_N_Frame,atom_type_dic)
        elif choice == "4":
            element_type = int(input("根据反应实际编号，输入目标元素id，如1: "))  # 输入转换为整数
            capture_target_atom_track_filter(element_type, atom_type_dic, Dic_N_cluster, List_N_Molecular)
        elif choice == "5":
            target = input("输入目标分子式，如Mo3O9: ")
            print(explore_next_states_info(str(target)))
        elif choice == "6":
            target = input("输入目标分子式，如Mo3O9: ")
            print(explore_previous_states_info(str(target)))
        elif choice == "7":
            top_n = int(input("输入网络大小参数（建议从小开始测试，例如5）: "))
            atuo_plot_reaction_networks(top_n,plot_graph=True)
        elif choice == "8":
            # 对于手动绘制反应网络，处理输入为列表的情况
            top_n = int(input("输入网络大小参数（建议从小开始测试，例如5）: "))
            start_states = input("输入起始分子状态，例如Mo3O9(可以为空，使用','分隔多个状态): ").split(',')
            end_states = input("输入终点分子状态，例如MoS6,MoS7(可以为空，使用','分隔多个状态): ").split(',')
            if start_states == ['']: start_states = []  # 处理空输入
            if end_states == ['']: end_states = []
            manual_plot_reaction_networks(top_n, start_states, end_states, plot_graph=True)

        elif choice == "0":
            print("退出程序...")
            break
        else:
            print("无效选择，请输入0到8之间的数字。")

if __name__ == "__main__":
    main()

  


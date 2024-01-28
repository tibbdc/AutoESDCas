#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @FileName     :draw_tool.py
# @Time         :2024/01/27 14:09:39
# @Author       :YangChunhe
# @Email        :2393492851@qq.com
# @Description  :file content
import os 
import pandas as pd
import numpy as np 
import warnings          
warnings.filterwarnings('ignore')
import re  
from os.path import exists
import sys
sys.path.append("./utils/")
from hr_analysis_tool import get_uha_dha_hr_offtarget


def calculate_percentages(fail_target_list):
    total_failures = sum(fail_target_list)
    percentages = [(fail_target / total_failures) * 100 for fail_target in fail_target_list]
    formatted_percentages = [f'{percentage:.2f}%' for percentage in percentages]
    return formatted_percentages

def draw_venn( evaluate_result_path ):
   
    evaluate_result = pd.read_excel(evaluate_result_path)
    u_d_hr_offTarget_id, u_hr_offTarget_id, d_hr_offTarget_id = get_uha_dha_hr_offtarget(evaluate_result)
    u_d_offTarget_id = set(list(u_d_hr_offTarget_id) + list(u_hr_offTarget_id) + list(d_hr_offTarget_id))


    import matplotlib.pyplot as plt
    from matplotlib_venn import venn2

    # 假设您有两组数据集合A和B
    u_hr =  set(list(u_hr_offTarget_id) + list(u_d_hr_offTarget_id))   
    d_hr =   set(list(d_hr_offTarget_id) + list(u_d_hr_offTarget_id) )    
    # u_d_hr = set(list(u_hr_offTarget_id) + list(u_d_hr_offTarget_id) + list(d_hr_offTarget_id))

    # 创建韦恩图
    venn = venn2(subsets=(u_hr, d_hr), set_labels=('UHA', 'DHA'), set_colors=('y','c'))

    # 显示图表
    # plt.title('Number of HR off-target events')
    plt.show()  
    
    print(f"只上游同源臂脱靶：{len(u_hr_offTarget_id)}--------只下游同源臂脱靶：{len(d_hr_offTarget_id)}---------上、下同时脱靶：{len(u_d_hr_offTarget_id)}------>总脱靶：{len(u_d_offTarget_id)} ")

def draw_pie(x,y,x_label='sucesss target',y_label='fail target'):

    import matplotlib.pyplot as plt
    # 示例数据
    labels = [x_label, y_label]
    sizes = [ x, y]
    colors = ['lightblue', 'lightgreen']

    # 创建饼图
    plt.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=140)

    # 添加标题
    plt.title('Design target result statistics')

    # 显示饼图
    plt.show()

def draw_bar( values, percent, categories):

    import matplotlib.pyplot as plt

    # 创建柱状图，设置颜色为浅绿色，宽度为0.5
    bars = plt.bar(categories, values, color='skyblue', width=0.5)
    li_percent_bar = zip(bars, percent)


    # 添加数据标签
    for i in li_percent_bar:
        bar,percent_cell = i
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width() / 2, height + 1, f'({percent_cell})', 
                ha='center', va='bottom', fontsize=10, color='black')

    # 调整Y轴范围
    plt.ylim(0, max(values) + 10)

    # 添加标题和标签
    plt.title('Optimization Comparison')
    plt.xlabel('Categories')
    plt.ylabel('Values')

    # 调整横坐标的位置，使其居中
    plt.xticks(range(len(categories)), categories)

    # 显示柱状图
    plt.show()

def draw_bar(values, percent, categories):
    import matplotlib.pyplot as plt

    # 创建柱状图，设置颜色为浅绿色，宽度为0.5
    bars = plt.bar(categories, values, color='skyblue', width=0.5)
    li_percent_bar = zip(bars, percent)

    # 添加数据标签
    for i in li_percent_bar:
        bar, percent_cell = i
        height = bar.get_height()
        value = values[bars.index(bar)]  # 获取每个柱的值
        plt.text(bar.get_x() + bar.get_width() / 2, height + 1, f'{value} ({percent_cell})',
                 ha='center', va='bottom', fontsize=10, color='black')

    # 调整Y轴范围
    plt.ylim(0, max(values) + 10)

    # 添加标题和标签
    plt.title('Optimization Comparison')
    plt.xlabel('Categories')
    plt.ylabel('Values')

    # 调整横坐标的位置，使其居中
    plt.xticks(range(len(categories)), categories)

    # 显示柱状图
    plt.show()

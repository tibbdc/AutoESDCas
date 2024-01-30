#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @FileName     :cp_tool.py
# @Time         :2024/01/30 11:10:03
# @Author       :YangChunhe
# @Email        :2393492851@qq.com
# @Description  :file content

import os
from os.path import exists,splitext,dirname,splitext,basename,realpath,abspath

def excute_cmd(ori_file_path, cur_file_path):
    import os
    cmd = f"cp  {ori_file_path} {cur_file_path}"
    state = os.system(cmd)
    print(cmd, state)

def cp_one_tool( genome_name_list, opt_type, origin_dir, current_dir, file_name):

    for genome in genome_name_list:  
        for i in range(50,550,50):   
            for step in range(1,6,1):
                index = str(i)+"-"+str(step)

                if opt_type != "":
                    e_s_d_path = f"{genome}/{index}/opt_{opt_type}_whole_process_output/edit_sequence_desgin/one_plasmid_system_result"
                    cur_dir_e_s_d_path = os.path.join(current_dir, e_s_d_path)   
                    ori_dir_e_s_d_path = os.path.join(origin_dir, e_s_d_path)         
                    
                    if not exists( cur_dir_e_s_d_path ):
                        os.makedirs( cur_dir_e_s_d_path )
                    ori_file_path, cur_file_path = os.path.join(ori_dir_e_s_d_path, file_name ), os.path.join(cur_dir_e_s_d_path, file_name )
                    excute_cmd(ori_file_path, cur_file_path)
                else:
                    input_path = f"{genome}/{index}"
                    cur_dir_input_path = os.path.join(current_dir, input_path)   
                    ori_dir_input_path = os.path.join(origin_dir, input_path)    
                    if not exists( cur_dir_input_path ):
                        os.makedirs( cur_dir_input_path )
                        
                    if file_name == "":
                        ori_file_path, cur_file_path = os.path.join(ori_dir_input_path, f"{index}.csv" ), os.path.join(cur_dir_input_path, f"{index}.csv")
                        excute_cmd(ori_file_path, cur_file_path)    

def cp_time_tool(genome_name_list,opt_type,origin_dir,current_dir):

    for genome_name in genome_name_list:

        file_name = f"{genome_name}_{opt_type}.csv"

        cur_dir_time_path = os.path.join(current_dir, file_name)   
        ori_dir_time_path = os.path.join(origin_dir, file_name) 
        
        excute_cmd(ori_dir_time_path, cur_dir_time_path)

def cp_genome_tool(genome_name_list,origin_dir,current_dir):

     for genome_name in genome_name_list:

        file_name = f"{genome_name}.fna"

        cur_dir_time_path = os.path.join(current_dir, file_name)   
        ori_dir_time_path = os.path.join(origin_dir, file_name) 
        
        excute_cmd(ori_dir_time_path, cur_dir_time_path)

def main():
    import os
    from os.path import exists,splitext,dirname,splitext,basename,realpath,abspath

    genome_name_list =  ['Escherichia_coli','Corynebacterium_glutamicum', 'Bacillus_subtilis', 'Salmonella_typhimurium', 'Pseudomonas_aeruginosa']

    origin_dir = "/home/yanghe/program/data_preprocessing/genome_test"
    current_dir = "/home/yanghe/public_github/AutoESDCas/runtime"

    success_file_name = "one_plasmid_design_result.xlsx"
    fail_file_name = "one_plasmid_design_F_Result.xlsx" 

    #
    cp_one_tool(genome_name_list=genome_name_list, opt_type="before", origin_dir=origin_dir, current_dir=os.path.join(current_dir,"data"), file_name="one_plasmid_design_result.xlsx")
    cp_one_tool(genome_name_list=genome_name_list, opt_type="before", origin_dir=origin_dir, current_dir=os.path.join(current_dir,"data"), file_name="one_plasmid_design_F_Result.xlsx")

    cp_one_tool(genome_name_list=genome_name_list, opt_type="after", origin_dir=origin_dir, current_dir=os.path.join(current_dir,"data"), file_name="one_plasmid_design_result.xlsx")
    cp_one_tool(genome_name_list=genome_name_list, opt_type="after", origin_dir=origin_dir, current_dir=os.path.join(current_dir,"data"), file_name="one_plasmid_design_F_Result.xlsx")

    cp_one_tool(genome_name_list=genome_name_list, opt_type="", origin_dir=origin_dir, current_dir=os.path.join(current_dir,"data"), file_name="")
    #
    cp_time_tool(genome_name_list,"after",origin_dir, os.path.join(current_dir,"data"))
    cp_time_tool(genome_name_list,"before",origin_dir, os.path.join(current_dir,"data"))
    #
    cp_genome_tool(genome_name_list,origin_dir,os.path.join(current_dir,"data"))

if __name__ == '__main__':

    main()



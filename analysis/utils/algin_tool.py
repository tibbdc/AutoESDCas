#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @FileName     :algin_tool.py
# @Time         :2024/01/27 14:02:59
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

BOWTIE_PATH = "/home/yanghe/software/bowtie"

#序列比对工具
def df_to_fasta(df,fasta_filename,id,sequence_name):
    with open(fasta_filename, 'w') as fasta_file:
        for index, row in df.iterrows():
            sequence_id = row[id]
            sequence = row[sequence_name]
            fasta_file.write(f'>{sequence_id}\n{sequence}\n')

def parse_sam_file(sam_file_path = "alignment.sam"):

    # 定义用于存储比对结果的列表
    alignment_data = []  

    with open(sam_file_path, "r") as sam_file:
        for line in sam_file:
            if not line.startswith("@"):  # 跳过SAM文件头部
                fields = line.strip().split("\t")
                read_name = fields[0]  
                chain = fields[1]
                reference_name = fields[2]
                reference_start = int(fields[3])  
                sequence = fields[9]
                mismatch = fields[-2]
                matching_number = fields[-1]

                alignment_data.append([read_name, chain, reference_name, reference_start, sequence, mismatch, matching_number])

    # 创建DataFrame
    columns = ["ReadName","chain", "ReferenceName", "ReferenceStart", "Sequence","Mismatch", "MatchingNumber"]
    alignment_df = pd.DataFrame(alignment_data, columns=columns)

    return alignment_df

def bowtie_seq_genome(tmp_path, genome_path, fasta_filename):

    if not exists(tmp_path):
        os.makedirs(tmp_path)

    sam_path = os.path.join(tmp_path ,'output.sam')  

    index_prefix = os.path.join(tmp_path, "genome_index")

    if not exists(index_prefix):
        os.makedirs(index_prefix)

        cmd = f'{BOWTIE_PATH}/bowtie-build {genome_path} {index_prefix}'
        os.system(cmd)  
    # os.chmod(f'{config.BOWTIE_PATH}bowtie', 0o755)
    cmd = f'{BOWTIE_PATH}/bowtie -p 2 -v 3 --sam-nohead -k 1000000000 {index_prefix} -f {fasta_filename} -S {sam_path}  > /dev/null 2>&1'
    os.system(cmd)

    #解析  
    sam_df = parse_sam_file(sam_path)

    #删除文件

    return sam_df

def get_no_offtarget(no_off_target,name):


    no_off_target["len"] = no_off_target.Sequence.apply(lambda x: len(x))
    no_off_target["ReferenceStart"] = no_off_target["ReferenceStart"]
    no_off_target["ReferenceEnd"] = no_off_target["ReferenceStart"] + no_off_target["len"]
    no_off_target["ReferenceStart"]  = no_off_target["ReferenceStart"].astype("str")
    no_off_target["ReferenceEnd"] = no_off_target["ReferenceEnd"].astype("str")
    no_off_target[f"After optimization {name} location (start-end)"] = no_off_target["ReferenceName"] +":"+ no_off_target["ReferenceStart"]+"-"+no_off_target["ReferenceEnd"]
    return no_off_target

def algin(A_no_offtarget_B_no_result_primer_df,sequence_name, fasta_file, tmp_path, genome_path, chain):

    if not exists(tmp_path):
        os.makedirs(tmp_path)

    df_to_fasta(A_no_offtarget_B_no_result_primer_df,fasta_file,"ID", sequence_name)
    sam_df = bowtie_seq_genome(tmp_path, genome_path, fasta_file)

    # no_off_target = sam_df[(sam_df['MatchingNumber'] == 'XM:i:1') & (sam_df['Mismatch'] == 'NM:i:0') ] 
    temp = sam_df.groupby("ReadName").apply(lambda x: x[x["chain"]==chain] if len( x[x["chain"]==chain]) == 1 else None )    

    no_off_target = temp.reset_index(drop=True)
    no_off_target = get_no_offtarget(no_off_target, sequence_name)
    no_off_target = no_off_target[["ReadName", f"After optimization {sequence_name} location (start-end)"]]
    no_off_target = no_off_target.rename(columns={"ReadName":"ID"})

    return no_off_target
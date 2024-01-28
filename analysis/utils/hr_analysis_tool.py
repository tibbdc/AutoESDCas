#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @FileName     :hr_data_analysis_tool.py
# @Time         :2024/01/27 14:04:56
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


#同源臂数据统计工具
def lambda2cols(df,lambdaf,in_coln,to_colns):         #apply函数的助手

    if len(in_coln) == 2:
        df_=df.apply(lambda x: lambdaf(x[in_coln[0]], x[in_coln[1]]),
                     axis=1).apply(pd.Series)
    if len(in_coln) == 3:
        df_=df.apply(lambda x: lambdaf(x[in_coln[0]],x[in_coln[1]],x[in_coln[2]]),
                     axis=1).apply(pd.Series)
    elif len(in_coln) == 4:
        df_=df.apply(lambda x: lambdaf(x[in_coln[0]],x[in_coln[1]],x[in_coln[2]],x[in_coln[3]]),
                     axis=1).apply(pd.Series)
    elif len(in_coln) == 5:
         df_=df.apply(lambda x: lambdaf(x[in_coln[0]],x[in_coln[1]],x[in_coln[2]],x[in_coln[3]],x[in_coln[4]]),
                     axis=1).apply(pd.Series)
    df_.columns = to_colns
    df = df.join(df_)        
    return df

def get_duplicatedRow_by_someoneColumns(df,column_name = 'Name'):

    import pandas as pd

    # 假设您的DataFrame名为df，column_name为您要检查的列名
    column_name = column_name

    # 找出重复的行，将其标记为True
    duplicated_rows = df.duplicated(subset=[column_name], keep=False)

    # 通过布尔型Series筛选出重复的行
    duplicate_records = df[duplicated_rows]  

    return duplicate_records


def get_uha_dha_hr_offtarget(evaluate_result):
    
    evaluate_result = evaluate_result.groupby('query or source (gene) sequence id').apply(lambda x:  x if len(x[ (x['off target']=='high') | (x['off target']=='medium') ] ) > 0  else pd.DataFrame())
    evaluate_result = evaluate_result.reset_index(drop=True)
    evaluate_result = evaluate_result.drop_duplicates(subset = ['query or source (gene) sequence id'] )

    #去除无用列
    unique_off_target =  evaluate_result.drop_duplicates('subject or target (reference genome) sequence id')
    temp_df = unique_off_target['query or source (gene) sequence id'].str.split(';',expand = True )
    temp_df = temp_df.rename(columns={0:'Name',1:'position',2:'hr'})

    all_df = temp_df

    #上下游同源臂，同时脱靶----------Name相同，hr不同
    u_d_hr_offTarget = get_duplicatedRow_by_someoneColumns(df=all_df,column_name = 'Name')

    temp = pd.concat([all_df, u_d_hr_offTarget]).drop_duplicates(keep=False)
    #只上游同源臂脱靶
    u_hr_offTarget = temp[temp['hr'] == 'UHA']
    #只下游同源臂脱靶
    d_hr_offTarget = temp[temp['hr'] == 'DHA']

    u_d_hr_offTarget_id = set(u_d_hr_offTarget['Name']) 
    u_hr_offTarget_id = set(u_hr_offTarget['Name'])
    d_hr_offTarget_id = set(d_hr_offTarget['Name'])
    u_d_offTarget_id = set(list(u_d_hr_offTarget_id) + list(u_hr_offTarget_id) + list(d_hr_offTarget_id))

    return list(u_d_hr_offTarget_id), list(u_hr_offTarget_id), list(d_hr_offTarget_id)

def get_offtarget_genome_cor(evaluate_result):

    evaluate_result = evaluate_result.groupby('query or source (gene) sequence id').apply(lambda x:  x if len(x[ (x['off target']=='high') | (x['off target']=='medium') ] ) > 0  else pd.DataFrame())
    evaluate_result = evaluate_result.reset_index(drop=True)

    def work(cor,start,end):
        gene_id,chom_cor,hr = cor.split(";")
        chom, genome_start_end = chom_cor.split(":")
        genome_start, genome_end = genome_start_end.split("-")
        genome_start, genome_end = int(genome_start), int(genome_end)
        start = int(genome_start - 20000 + start-1)
        end =  int(genome_start - 20000 + end)

        genome_cor = chom + ":" + str(start)+ "-" + str(end)

        return gene_id,hr,genome_cor
    a = lambda2cols(evaluate_result,work,in_coln=["query or source (gene) sequence id","start of alignment in query","end of alignment in query"],to_colns=["locus_tag","hr","genome_cor"])
    a = a[["locus_tag","hr","genome_cor","off target"]]
    return a

def get_offtarget(a, hr_off_target_type,  u_d_hr_offTarget_id, u_hr_offTarget_id,  d_hr_offTarget_id ):
    hr_df = pd.DataFrame()
    data_id_list = [ list(u_d_hr_offTarget_id), list(u_hr_offTarget_id), list(d_hr_offTarget_id)  ]
    
    for data_id  in data_id_list:

        b = pd.DataFrame(columns=["locus_tag"], data=data_id)
        a_b = pd.merge(a, b)

        def work(x):
        
            locus_tag = list(x["locus_tag"])[0]
            
            uha = x[x["hr"]=="UHA"]
            uha_cor = set(list(uha["genome_cor"]))
            uha_target = uha[uha["off target"]=="target"]
            uha_target_cor = set(list(uha_target["genome_cor"]))

            dha = x[x["hr"]=="DHA"]
            dha_cor = set( list(dha["genome_cor"]) )
            dha_target =  dha[dha["off target"]=="target"]
            dha_target_cor = set(list(dha_target["genome_cor"]))

            df_dict = {
                "locus_tag":locus_tag,
                "UHA location (start-end)": ";".join( uha_target_cor ),
                "UHA off-target site": ";".join( uha_cor - uha_target_cor ),
                "DHA location (start-end)": ";".join(dha_target_cor),
                "DHA off-target site": ";".join(dha_cor - dha_target_cor),
                "type":hr_off_target_type
            }
            
            df = pd.DataFrame([df_dict])
            return df

        a_b_df = a_b.groupby("locus_tag").apply(lambda x: work(x))
        a_b_df = a_b_df.reset_index(drop=True)

        hr_df = hr_df.append(a_b_df)
        hr_df = hr_df.reset_index(drop=True)


    hr_df["locus_tag"] = hr_df.locus_tag.apply( lambda x: x.split("_")[0])

    return hr_df

def get_cor_A_no_offtarget_B_no_result(SEQUENCING_PRIMER_1_SEQUENCING_PRIMER_2, valid_primer_df ):

    for i,v in SEQUENCING_PRIMER_1_SEQUENCING_PRIMER_2.iterrows():
        id = v["ID"]
        if pd.isna(v["After optimization SEQUENCING_PRIMER_1 location (start-end)"]):
            temp = valid_primer_df[valid_primer_df.ID.str.contains(  f"valid_primer:{id}" )]
            
            cor = list(temp["off_target"])[0]
            cor = eval(cor)

            for k,value in cor.items():
                if len(value.split(";")) == 2 and "SEQUENCING_PRIMER_1" in k:
                    SEQUENCING_PRIMER_1_0 = value
                    chrom_cor_strand = SEQUENCING_PRIMER_1_0.split(";")[1]
                    chrom_cor,strand = chrom_cor_strand.split("|")
                    chrom,cor = chrom_cor.split(":")
                    SEQUENCING_PRIMER_1 = chrom_cor

            v["After optimization SEQUENCING_PRIMER_1 location (start-end)"] = SEQUENCING_PRIMER_1
        if pd.isna(v["After optimization SEQUENCING_PRIMER_2 location (start-end)"]):
            
            temp = valid_primer_df[valid_primer_df.ID.str.contains(  f"valid_primer:{id}" )]
            
            cor = list(temp["off_target"])[0]
            cor = eval(cor)

            for k,value in cor.items():
                if  len(value.split(";")) == 2 and "SEQUENCING_PRIMER_2" in k:
                    SEQUENCING_PRIMER_2_0 = value
                    chrom_cor_strand = SEQUENCING_PRIMER_2_0.split(";")[1]
                    chrom_cor,strand = chrom_cor_strand.split("|")
                    chrom,cor = chrom_cor.split(":")
                    SEQUENCING_PRIMER_2 = chrom_cor

            v["After optimization SEQUENCING_PRIMER_2 location (start-end)"] = SEQUENCING_PRIMER_2
    return SEQUENCING_PRIMER_1_SEQUENCING_PRIMER_2

def add_uha_dha_coord(hr_df_145,evaluate_result,hr_length ):
    temp = pd.DataFrame( columns=["gene"] ,data=set(evaluate_result["query or source (gene) sequence id"].str.replace("UHA", "").str.replace("DHA", "")) )
    temp = temp["gene"].str.split(";", expand=True)
    temp["locus_tag"]= temp[0].str.replace("_del", "")
    temp = temp.drop(columns=[0,2])
    temp = temp.rename(columns={1:"genome_coord"})
    temp_df = pd.merge(temp, hr_df_145)

    def work(genome_coord,uha,dha):
        chom,coord = genome_coord.split(":")
        start,end = coord.split("-")
        start,end = int(start),int(end)
        if uha == "" or pd.isna(uha):
            uha_cor = chom + ":" + str(start-hr_length) +"-"+ str(end)
            uha = uha_cor
        return uha

    temp_df["UHA location (start-end)"] = temp_df.apply(lambda x: work(x["genome_coord"], x["UHA location (start-end)"], x["DHA location (start-end)"]),axis=1)

    def work(genome_coord,uha,dha):
        chom,coord = genome_coord.split(":")
        start,end = coord.split("-")
        start,end = int(start),int(end)

        if dha == "" or pd.isna(dha):
            dha_cor = chom + ":" + str(end) +"-"+ str(end + hr_length)
            dha = dha_cor
        return dha

    temp_df["DHA location (start-end)"] = temp_df.apply(lambda x: work(x["genome_coord"], x["UHA location (start-end)"], x["DHA location (start-end)"]),axis=1)

    temp_df = temp_df.drop(columns=["genome_coord"])

    return temp_df
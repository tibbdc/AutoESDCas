#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @FileName     :test_analysis.py
# @Time         :2024/01/30 11:09:31
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
input_dir = "/home/yanghe/public_github/AutoESDCas/runtime/data"
from os.path import exists,splitext,dirname,splitext,basename,realpath,abspath


def get_duplicatedRow_by_someoneColumns(df,column_name = 'Name'):

    import pandas as pd

    # 假设您的DataFrame名为df，column_name为您要检查的列名
    column_name = column_name

    # 找出重复的行，将其标记为True
    duplicated_rows = df.duplicated(subset=[column_name], keep=False)

    # 通过布尔型Series筛选出重复的行
    duplicate_records = df[duplicated_rows]  

    return duplicate_records

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

def bowtie_seq_genome(tmp_genome_index ,tmp_sam_path, genome_path, fasta_filename):

    if not exists(tmp_genome_index):
        os.makedirs(tmp_genome_index)  
    
    if not exists(tmp_sam_path):
        os.makedirs(tmp_sam_path)
    
    sam_path = os.path.join(tmp_sam_path ,'output.sam')  
    index_prefix = os.path.join(tmp_genome_index, "genome_index")

    if not exists(index_prefix):
        os.makedirs(index_prefix)

        cmd = f'{BOWTIE_PATH}/bowtie-build {genome_path} {index_prefix}'
        os.system(cmd)  
    # os.chmod(f'{config.BOWTIE_PATH}bowtie', 0o755)
    cmd = f'{BOWTIE_PATH}/bowtie -p 2 -v 3 --sam-nohead -k 1000000000 {index_prefix} -f {fasta_filename} -S {sam_path}'
    os.system(cmd)

    #解析  
    sam_df = parse_sam_file(sam_path)

    #删除文件

    return sam_df

def get_no_offtarget(no_off_target,name):

    if len( no_off_target ) > 0:
        no_off_target["len"] = no_off_target.Sequence.apply(lambda x: len(x))
        no_off_target["ReferenceStart"] = no_off_target["ReferenceStart"]
        no_off_target["ReferenceEnd"] = no_off_target["ReferenceStart"] + no_off_target["len"]
        no_off_target["ReferenceStart"]  = no_off_target["ReferenceStart"].astype("str")
        no_off_target["ReferenceEnd"] = no_off_target["ReferenceEnd"].astype("str")
        no_off_target[f"After optimization {name} location (start-end)"] = no_off_target["ReferenceName"] +":"+ no_off_target["ReferenceStart"]+"-"+no_off_target["ReferenceEnd"]
    else:
        no_off_target = pd.DataFrame( columns = ["ReadName", f"After optimization {name} location (start-end)"]  )

    return no_off_target

def algin(A_no_offtarget_B_no_result_primer_df,sequence_name, fasta_file, tmp_genome_index, tmp_sam_path, genome_path, chain):

   

    if not exists( dirname(fasta_file) ):
        os.makedirs( dirname(fasta_file) )

    df_to_fasta(A_no_offtarget_B_no_result_primer_df,fasta_file,"ID", sequence_name) 

    # /home/yanghe/public_github/AutoESDCas/runtime/data/Escherichia_coli/xxx/
    # /home/yanghe/public_github/AutoESDCas/runtime/data/Escherichia_coli/50-1/opt_after_whole_process_output/edit_sequence_desgin/primer_bowtie/xxx/

   
    sam_df = bowtie_seq_genome(tmp_genome_index ,tmp_sam_path, genome_path, fasta_file)

    # no_off_target = sam_df[(sam_df['MatchingNumber'] == 'XM:i:1') & (sam_df['Mismatch'] == 'NM:i:0') ] 
    temp = sam_df.groupby("ReadName").apply(lambda x: x[x["chain"]==chain] if len( x[x["chain"]==chain]) == 1 else None )    

   
    no_off_target = temp.reset_index(drop=True)
    no_off_target = get_no_offtarget(no_off_target, sequence_name)
    no_off_target = no_off_target[["ReadName", f"After optimization {sequence_name} location (start-end)"]]
    no_off_target = no_off_target.rename(columns={"ReadName":"ID"})

    return no_off_target

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

def get_invalid_primer(invalid_priemr_df):
    li =[]
    for i,v in invalid_priemr_df.iterrows(): 

        if v["ID"] !=  "invalid_primer:NC_002516.2_44_substitution;NC_002516.2:3600412-3600477:SEQUENCING_PRIMER_1_AND_SEQUENCING_PRIMER_2:off_target" \
            and  v["ID"] != "invalid_primer:NC_002516.2_215_insertion;NC_002516.2:1769525-1769526:SEQUENCING_PRIMER_1_AND_SEQUENCING_PRIMER_2:off_target" \
            and v["ID"] != "invalid_primer:NC_002516.2_167_insertion;NC_002516.2:6129485-6129486:SEQUENCING_PRIMER_1_AND_SEQUENCING_PRIMER_2:off_target" \
            and v["ID"] != "invalid_primer:NC_002516.2_396_insertion;NC_002516.2:258276-258277:SEQUENCING_PRIMER_1_AND_SEQUENCING_PRIMER_2:off_target": 
        #id = v["ID"].split(":")[1]  
            state,id,cor,name,offtarget = v["ID"].split(":")  #PCR
            id = id + ":" + cor                               #PCR
            off_target = v["off_target"]
            off_target = eval(off_target)

            print(id)   #Pseudomonas_aeruginosa 500-1
        

            #目标引物  
            SEQUENCING_PRIMER_1 = ""
            SEQUENCING_PRIMER_1_offtarget = ""
            SEQUENCING_PRIMER_2 = ""
            SEQUENCING_PRIMER_2_offtarget = ""

            #脱靶引物
            coordinates_1 = []  
            coordinates_2 = []
            for i,v in off_target.items():
                
                if len(v.split(";")) == 2 and "SEQUENCING_PRIMER_1" in i:
                    SEQUENCING_PRIMER_1_0 = v    
                    chrom_cor_strand = SEQUENCING_PRIMER_1_0.split(";")[1]
                    chrom_cor,strand = chrom_cor_strand.split("|")
                    chrom,cor = chrom_cor.split(":")
                    SEQUENCING_PRIMER_1 =  chrom_cor
                
                if  len(v.split(";")) == 2 and "SEQUENCING_PRIMER_2" in i:
                    SEQUENCING_PRIMER_2_0 = v
                    chrom_cor_strand = SEQUENCING_PRIMER_2_0.split(";")[1]
                    chrom_cor,strand = chrom_cor_strand.split("|")
                    chrom,cor = chrom_cor.split(":")
                    SEQUENCING_PRIMER_2 = chrom_cor 
            
            # for i,v in off_target.items(): 
                chrom_cor_strand = v.split(";")[1]
                chrom_cor,strand = chrom_cor_strand.split("|")

                if "SEQUENCING_PRIMER_1" in i and len(v.split(";")) != 2:  
                    coordinates_1.append( chrom_cor )
                if "SEQUENCING_PRIMER_2" in i and len(v.split(";")) != 2: 
                    coordinates_2.append( chrom_cor )
            
            off_target_type = "" 
            if  SEQUENCING_PRIMER_1 !="": 
                SEQUENCING_PRIMER_1_offtarget = ";".join(coordinates_1)  
                off_target_type = "SEQUENCING_PRIMER_1 offtarget" 
            elif  SEQUENCING_PRIMER_2 !="":
                SEQUENCING_PRIMER_2_offtarget = ";".join(coordinates_2)
                off_target_type = "SEQUENCING_PRIMER_2 offtarget"
            if SEQUENCING_PRIMER_1 !="" and SEQUENCING_PRIMER_2 !="":
                SEQUENCING_PRIMER_1_offtarget = ";".join(coordinates_1)
                SEQUENCING_PRIMER_2_offtarget = ";".join(coordinates_2)
                off_target_type = "SEQUENCING_PRIMER_1 offtarget" + "&" +  "SEQUENCING_PRIMER_2 offtarget"
                
            dict_offtarget = {
                    "locus_tag" : id.split(";")[0],  
                    "SEQUENCING_PRIMER_1 location (start-end)":SEQUENCING_PRIMER_1,
                    "SEQUENCING_PRIMER_1 off-target site":SEQUENCING_PRIMER_1_offtarget,
                    "SEQUENCING_PRIMER_2 location (start-end)":SEQUENCING_PRIMER_2,
                    "SEQUENCING_PRIMER_2 off-target site":SEQUENCING_PRIMER_2_offtarget,
                    "type (B: Before optimization);(A:After optimization)":off_target_type
                }
            li.append(dict_offtarget)
        else:
            state,id,cor,name,offtarget = v["ID"].split(":")  #PCR
            id = id + ":" + cor                               #PCR
            off_target = v["off_target"]
            off_target = eval(off_target)

            print(id)   #Pseudomonas_aeruginosa 350-3
        

            #目标引物  
            SEQUENCING_PRIMER_1 = ""
            SEQUENCING_PRIMER_1_offtarget = ""
            SEQUENCING_PRIMER_2 = ""
            SEQUENCING_PRIMER_2_offtarget = ""

            #脱靶引物
            coordinates_1 = []  
            coordinates_2 = []
            for i,v in off_target.items():
                
                if len(v.split(";")) == 1 and "SEQUENCING_PRIMER_1" in i:
                    SEQUENCING_PRIMER_1_0 = v    
                    chrom_cor_strand = SEQUENCING_PRIMER_1_0.split(";")[0]
                    chrom_cor,strand = chrom_cor_strand.split("|")
                    chrom,cor = chrom_cor.split(":")
                    SEQUENCING_PRIMER_1 =  chrom_cor
                
                if  len(v.split(";")) == 1 and "SEQUENCING_PRIMER_2" in i:
                    SEQUENCING_PRIMER_2_0 = v
                    chrom_cor_strand = SEQUENCING_PRIMER_2_0.split(";")[0]
                    chrom_cor,strand = chrom_cor_strand.split("|")
                    chrom,cor = chrom_cor.split(":")
                    SEQUENCING_PRIMER_2 = chrom_cor 
            
            # for i,v in off_target.items(): 
                chrom_cor_strand = v.split(";")[0]
                chrom_cor,strand = chrom_cor_strand.split("|")
                
                if "SEQUENCING_PRIMER_1" in i and len(v.split(";")) != 1:  
                    coordinates_1.append( chrom_cor )
                if "SEQUENCING_PRIMER_2" in i and len(v.split(";")) != 1: 
                    coordinates_2.append( chrom_cor )

            off_target_type = "" 
            if  SEQUENCING_PRIMER_1 !="": 
                SEQUENCING_PRIMER_1_offtarget = ";".join(coordinates_1)  
                off_target_type = "SEQUENCING_PRIMER_1 offtarget" 
            elif  SEQUENCING_PRIMER_2 !="":
                SEQUENCING_PRIMER_2_offtarget = ";".join(coordinates_2)
                off_target_type = "SEQUENCING_PRIMER_2 offtarget"
            if SEQUENCING_PRIMER_1 !="" and SEQUENCING_PRIMER_2 !="":
                SEQUENCING_PRIMER_1_offtarget = ";".join(coordinates_1)
                SEQUENCING_PRIMER_2_offtarget = ";".join(coordinates_2)
                off_target_type = "SEQUENCING_PRIMER_1 offtarget" + "&" +  "SEQUENCING_PRIMER_2 offtarget"
                
            dict_offtarget = {
                    "locus_tag" : id.split(";")[0],  
                    "SEQUENCING_PRIMER_1 location (start-end)":SEQUENCING_PRIMER_1,
                    "SEQUENCING_PRIMER_1 off-target site":SEQUENCING_PRIMER_1_offtarget,
                    "SEQUENCING_PRIMER_2 location (start-end)":SEQUENCING_PRIMER_2,
                    "SEQUENCING_PRIMER_2 off-target site":SEQUENCING_PRIMER_2_offtarget,
                    "type (B: Before optimization);(A:After optimization)":off_target_type
                }
            li.append(dict_offtarget)
    
    
    
    
    
    
    if len(li)>0:
        df = pd.DataFrame(li)
    else:
         
        df = pd.DataFrame(columns=["locus_tag"])
        

    return df

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

def count_optimal_valid_invalid_fail_nums(one_plasmid_design_F_Result, Test_primer_G, Primer_g_offTarget ):

    #1.primer3设计出来的所有引物
    Test_primer_G_primer_id = set(Test_primer_G['ID'])

    #2.脱靶的所有引物id
    # off_target_primer_id = set(Primer_g_offTarget['ID'].str.split(':',expand=True)[1])
    off_target_primer_id = set( Primer_g_offTarget['ID'].str.split(':',expand=True)[1]+":"+Primer_g_offTarget['ID'].str.split(':',expand=True)[2] )  #PCR

    #3.设计失败的所有引物id
    if len( one_plasmid_design_F_Result ) >0:
        one_plasmid_design_F_Result_id = set(one_plasmid_design_F_Result['ID'])
    else:
        one_plasmid_design_F_Result_id = set()

    #4.所有引物id
    all_gene_id =  set(list(Test_primer_G_primer_id) + list(one_plasmid_design_F_Result_id))
    
    #5.无效引物df
    invalid_priemr_df = Primer_g_offTarget[Primer_g_offTarget.ID.str.contains('invalid_primer')]
	
	#6.无效引物id
    if len(invalid_priemr_df) >0 : 
        invalid_priemr_id = set( invalid_priemr_df['ID'].str.split(':',expand=True)[1] + ":" + invalid_priemr_df['ID'].str.split(':',expand=True)[2] ) #PCR
		# invalid_priemr_id = set(invalid_priemr_df['ID'].str.split(':',expand=True)[1])
    else:
        invalid_priemr_id = set()
    
    #6.最优引物optimal
    optimal_primer_id = Test_primer_G_primer_id - off_target_primer_id

    #7.有效引物valid
    valid_primer_id = off_target_primer_id - invalid_priemr_id


    print(f"全部引物：{len(all_gene_id)}------最优引物：{len(optimal_primer_id)}-------有效引物:{len(valid_primer_id)}-----------无效引物：{len(invalid_priemr_id)}----------设计失败引物：{len(one_plasmid_design_F_Result_id)}")


    all_offtarget =  set(invalid_priemr_df[invalid_priemr_df.ID.str.contains('SEQUENCING_PRIMER_1_AND_SEQUENCING_PRIMER_2')]['ID'])   #PCR
    one_offtarget = set(invalid_priemr_df[~ invalid_priemr_df.ID.str.contains('SEQUENCING_PRIMER_1_AND_SEQUENCING_PRIMER_2')]['ID'] )  #PCR
	
	#all_offtarget =  set(invalid_priemr_df[invalid_priemr_df.Name.str.contains('SEQUENCING_PRIMER_1_AND_SEQUENCING_PRIMER_2')]['ID'])
	#one_offtarget = set(invalid_priemr_df[~ invalid_priemr_df.Name.str.contains('SEQUENCING_PRIMER_1_AND_SEQUENCING_PRIMER_2')]['ID'] ) 

    print(f"左右引物全部脱靶导致无效引物：{len( all_offtarget )}-------左或右引物脱靶:{len( one_offtarget )}----------primer3设计不出来的引物：{len(one_plasmid_design_F_Result_id)}")


    return invalid_priemr_df, one_plasmid_design_F_Result_id

def rename_invalid_primer(invalid_primer, opt_type ):

    invalid_primer = invalid_primer.rename( columns={
                                                        "SEQUENCING_PRIMER_1 location (start-end)": f"{opt_type} optimization SEQUENCING_PRIMER_1 location (start-end)",
                                                        "SEQUENCING_PRIMER_1 off-target site":f"{opt_type} optimization SEQUENCING_PRIMER_1 off-target site",
                                                        "SEQUENCING_PRIMER_2 location (start-end)":f"{opt_type} optimization SEQUENCING_PRIMER_2 location (start-end)",
                                                        "SEQUENCING_PRIMER_2 off-target site":f"{opt_type} optimization  SEQUENCING_PRIMER_2 off-target site"
                                                    } )
    
    if  opt_type == "Before":
        invalid_primer["type (B: Before optimization);(A:After optimization)"] = invalid_primer["type (B: Before optimization);(A:After optimization)"].apply(lambda x: "B:"+x)
        invalid_primer["type (B: Before optimization);(A:After optimization)"] = invalid_primer["type (B: Before optimization);(A:After optimization)"].astype(str)
    elif opt_type == "After":
        invalid_primer["type (B: Before optimization);(A:After optimization)"] = invalid_primer["type (B: Before optimization);(A:After optimization)"].astype(str)
        invalid_primer["type (B: Before optimization);(A:After optimization)"] = invalid_primer["type (B: Before optimization);(A:After optimization)"].apply(lambda x:  "All primer templates offtarget"  if x == "" else x )
        invalid_primer["type (B: Before optimization);(A:After optimization)"] =  "A:" + invalid_primer["type (B: Before optimization);(A:After optimization)"]

    
    return invalid_primer

def extract_failture_primer(one_plasmid_design_F_Result_id):

    fail_primer_df = pd.DataFrame( columns=["locus_tag"] ,data = list(one_plasmid_design_F_Result_id) )  
    fail_primer_df["locus_tag"] = fail_primer_df.locus_tag.apply(lambda x: x.split(";")[0] )
    fail_primer_df["type (B: Before optimization);(A:After optimization)"] = "No primer3 design results"
    
    return fail_primer_df

def extract_B_offtarget_A_no_offtarget( opt_before_invalid_primer, opt_after_invalid_primer, Primer_g_offTarget ,Test_primer_G, algin_path):

    fasta_file, genome_index, tmp_sam_path, genome_path = algin_path 

    A_no_offtarget_B_no_result = set(opt_before_invalid_primer["locus_tag"]) - set(opt_after_invalid_primer["locus_tag"])

    A_no_offtarget_B_no_result_df = pd.DataFrame(columns=["ID"], data = list(A_no_offtarget_B_no_result) )

    Test_primer_G["ID"] = Test_primer_G.ID.apply(lambda x: x.split(";")[0]  )
    A_no_offtarget_B_no_result_primer_df = pd.merge( Test_primer_G, A_no_offtarget_B_no_result_df )

    A_no_offtarget_B_no_result_primer_df = A_no_offtarget_B_no_result_primer_df[["ID","SEQUENCING_PRIMER_1","SEQUENCING_PRIMER_2"]]

    print( A_no_offtarget_B_no_result_primer_df )


    SEQUENCING_PRIMER_1 = algin(A_no_offtarget_B_no_result_primer_df, "SEQUENCING_PRIMER_1",fasta_file, genome_index, tmp_sam_path, genome_path, "0")
    SEQUENCING_PRIMER_2 = algin(A_no_offtarget_B_no_result_primer_df, "SEQUENCING_PRIMER_2",fasta_file, genome_index, tmp_sam_path, genome_path, "16")  

    SEQUENCING_PRIMER_1_SEQUENCING_PRIMER_2 = pd.merge(SEQUENCING_PRIMER_1, SEQUENCING_PRIMER_2, how="outer")

    valid_primer_df = Primer_g_offTarget[ Primer_g_offTarget.ID.str.contains('valid_primer') &  (~Primer_g_offTarget.ID.str.contains('invalid_primer')) ]

    SEQUENCING_PRIMER_1_SEQUENCING_PRIMER_2 = get_cor_A_no_offtarget_B_no_result(SEQUENCING_PRIMER_1_SEQUENCING_PRIMER_2, valid_primer_df )

    SEQUENCING_PRIMER_1_SEQUENCING_PRIMER_2 = SEQUENCING_PRIMER_1_SEQUENCING_PRIMER_2.rename(columns={"ID":"locus_tag"})

    SEQUENCING_PRIMER_1_SEQUENCING_PRIMER_2["type (B: Before optimization);(A:After optimization)"] = "A:No offtarget"

    return SEQUENCING_PRIMER_1_SEQUENCING_PRIMER_2

def extract_opt_invalid_primer(Test_primer_G, invalid_primer, alin_path, opt_type, opt_type_value = "A:All primer templates offtarget"):

    fasta_file, genome_index, tmp_sam_path, genome_path = alin_path
    
    ab = Test_primer_G.copy()  
    ab["ID"] = list( ab.ID.apply(lambda x: x.split(";")[0]) )
    ab.index = ab["ID"]

    for i,v in invalid_primer.iterrows():
        id = v["locus_tag"]
        if v["type (B: Before optimization);(A:After optimization)"] != opt_type_value:
                
                if pd.isna( v[f"{opt_type} optimization SEQUENCING_PRIMER_1 location (start-end)"] ) or  v[f"{opt_type} optimization SEQUENCING_PRIMER_1 location (start-end)"] == "":
                    left_primer = ab.loc[id,"SEQUENCING_PRIMER_1"]
                    seq_df = pd.DataFrame(
                            [{
                            "ID":id,
                            "SEQUENCING_PRIMER_1":left_primer
                            }]
                        )
                    result_df = algin(seq_df,"SEQUENCING_PRIMER_1",fasta_file, genome_index, tmp_sam_path, genome_path, "0")
                    invalid_primer.loc[i,f"{opt_type} optimization SEQUENCING_PRIMER_1 location (start-end)"] = result_df.loc[0,"After optimization SEQUENCING_PRIMER_1 location (start-end)"]

                if pd.isna(  v[f"{opt_type} optimization SEQUENCING_PRIMER_2 location (start-end)"]  ) or  v[f"{opt_type} optimization SEQUENCING_PRIMER_2 location (start-end)"]=="":
                    

                    right_primer = ab.loc[id,"SEQUENCING_PRIMER_2"]
                    
                    seq_df = pd.DataFrame(
                            [{
                            "ID":id,
                            "SEQUENCING_PRIMER_2":right_primer
                            }]
                            
                        )
                    result_df = algin(seq_df,"SEQUENCING_PRIMER_2",fasta_file, genome_index, tmp_sam_path, genome_path, "16")
            
                    invalid_primer.loc[i,f"{opt_type} optimization SEQUENCING_PRIMER_2 location (start-end)"] = result_df.loc[0,"After optimization SEQUENCING_PRIMER_2 location (start-end)"]

    return invalid_primer

def extract_df(genome, index, opt_type):  

    edit_result_path = f"{genome}/{index}/opt_{opt_type}_whole_process_output/edit_sequence_desgin/one_plasmid_system_result/one_plasmid_design_result.xlsx"
    edit_result_path = os.path.join(input_dir, edit_result_path)
    one_plasmid_design_F_Result_path = f"{genome}/{index}/opt_{opt_type}_whole_process_output/edit_sequence_desgin/one_plasmid_system_result/one_plasmid_design_F_Result.xlsx"
    one_plasmid_design_F_Result_path = os.path.join(input_dir, one_plasmid_design_F_Result_path)


    one_plasmid_design_F_Result = pd.read_excel(one_plasmid_design_F_Result_path, sheet_name='Test_primer_G')
    Test_primer_G = pd.read_excel(edit_result_path,sheet_name='Test_primer_G')
    Primer_g_offTarget =  pd.read_excel(edit_result_path,sheet_name='Primer_g_offTarget') 

    return Test_primer_G, Primer_g_offTarget, one_plasmid_design_F_Result 

def extract_alin_path(genome, index, opt_type):


    bowtie_workdir = f"{genome}/{index}/opt_{opt_type}_whole_process_output/edit_sequence_desgin/primer_bowtie"



    fasta_file_dir = os.path.join(input_dir,bowtie_workdir)
    fasta_file = os.path.join( fasta_file_dir, "sequence.fasta")

    genome_index = os.path.join(input_dir, f"{genome}/genome_index")
    genome_path = os.path.join(input_dir, f"{genome}.fna")
    tmp_sam_path =os.path.join(input_dir, ) 

    alin_path = fasta_file, genome_index, tmp_sam_path, genome_path

    return alin_path

def merge_opt_primer(opt_before_invalid_primer, opt_after_invalid_primer):

    df = pd.merge( opt_before_invalid_primer, opt_after_invalid_primer,on="locus_tag",how="outer")
    df = df.fillna("")
    if len(df) > 0: 
        df["type (B: Before optimization);(A:After optimization)_x"] = df["type (B: Before optimization);(A:After optimization)_x"].apply(lambda x:  "B:No offtarget" if x == "" else x )
        df["type (B: Before optimization);(A:After optimization)_y"] = df["type (B: Before optimization);(A:After optimization)_y"].apply(lambda x:  "A:No offtarget" if x == "" else x )

        df["type (B: Before optimization);(A:After optimization)"] = df["type (B: Before optimization);(A:After optimization)_x"] +";"+ df["type (B: Before optimization);(A:After optimization)_y"] 

        df = df.drop(columns=["type (B: Before optimization);(A:After optimization)_x", "type (B: Before optimization);(A:After optimization)_y"])
    else:

        columns = ['locus_tag',
                'Before optimization SEQUENCING_PRIMER_1 location (start-end)',
                'Before optimization SEQUENCING_PRIMER_1 off-target site',
                'Before optimization SEQUENCING_PRIMER_2 location (start-end)',
                'Before optimization  SEQUENCING_PRIMER_2 off-target site',
                'After optimization SEQUENCING_PRIMER_1 location (start-end)',
                'After optimization SEQUENCING_PRIMER_2 location (start-end)',
                'type (B: Before optimization);(A:After optimization)']

        df = pd.DataFrame(columns= columns)

    return df

def decorate_time( opt_before_all_df, opt_type="A" ):

    opt_before_all_df.set_index( ['strain', 'design target'], inplace=True )
    opt_before_all_df_mean = opt_before_all_df.groupby(["strain","design target"]).mean()
    opt_before_all_df_mean = opt_before_all_df_mean.rename(columns={"all_time":"mean time"})
    opt_before_all_df_mean_all  = pd.merge( opt_before_all_df ,opt_before_all_df_mean, left_index=True, right_index=True )

    opt_before_all_df_mean_all = opt_before_all_df_mean_all.reset_index()
    opt_before_all_df_mean_all["design target"]= opt_before_all_df_mean_all["design target"].astype("str")
    opt_before_all_df_mean_all["design target index"]= opt_before_all_df_mean_all["design target index"].astype("str")
    opt_before_all_df_mean_all["design target"] =  opt_before_all_df_mean_all["design target"] + "-" + opt_before_all_df_mean_all["design target index"]
    opt_before_all_df_mean_all.drop(columns=["design target index"], inplace=True)
    opt_before_all_df_mean_all.rename(columns={"all_time":f"{opt_type}:time"},inplace=True)
    opt_before_all_df_mean_all.rename(columns={"mean time": f"{opt_type}:mean time"},inplace=True)

    return opt_before_all_df_mean_all

def count_time(genome_name_list, columns_index, opt_type="before" ):

    opt_before_all_df = pd.DataFrame()
    opt_before_li_mean = []
    # opt_after_li_mean = []

    for genome_name in genome_name_list:
        # opt_type = "before"
        genome_path = f"{input_dir}/{genome_name}_{opt_type}.csv"
        df = pd.read_csv(genome_path)
        
        df = df.loc[:,columns_index]
        temp = df["design target"].str.split("-",expand=True)
        df["design target"] = temp[0]
        df["design target index"] = temp[1]
        df["design target"] = df["design target"].astype(int)
        df = df.sort_values("design target")

        opt_before_all_df = opt_before_all_df.append(df)

        data = df["all_time"].values.reshape(10,5)
        data_mean = np.mean(data, axis=1)
        opt_before_li_mean.append( data_mean )
        
    opt_before_data = np.vstack( opt_before_li_mean )
    opt_before_data = np.transpose( opt_before_data ) 

    return opt_before_all_df, opt_before_data

def compute_all_time( genome_name_list, columns_index, opt_before_after_time_path, draw_path):

    # genome_name_list  = ['Escherichia_coli','Corynebacterium_glutamicum', 'Pseudomonas_aeruginosa', 'Bacillus_subtilis', 'Salmonella_typhimurium']
    # columns_index = ['strain', 'design target', 'time']

    opt_before_all_df, opt_before_data = count_time(genome_name_list, columns_index, opt_type="before" )
    opt_before_all_df_mean_all = decorate_time( opt_before_all_df, opt_type="B" )

    opt_after_all_df, opt_after_data = count_time(genome_name_list, columns_index, opt_type="after" )
    opt_after_all_df_mean_all = decorate_time( opt_after_all_df, opt_type="A" )

    opt_before_after_all_df_mean_all = pd.merge( opt_before_all_df_mean_all, opt_after_all_df_mean_all )
    opt_before_after_all_df_mean_all.to_csv(opt_before_after_time_path)


    labels = [i for i in range(50, 550, 50)]
    draw_bar_gram(labels ,genome_name_list , opt_before_data, opt_after_data,  draw_path )

    return opt_before_data,opt_after_data,opt_before_after_all_df_mean_all
    
def draw_bar_gram(labels ,categories , opt_before_data, opt_after_data, path ):

    import matplotlib.pyplot as plt
    import numpy as np

    # 生成数据
    # labels = [i for i in range(50, 550, 50)]
    # categories = ['Escherichia_coli', 'Corynebacterium_glutamicum', 'Pseudomonas_aeruginosa', 'Bacillus_subtilis', 'Salmonella_typhimurium']

    # 设置图形大小和子图
    fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(15, 6))

    # 绘制条形统计图 - opt-before
    bar_width = 0.15
    for i, category in enumerate(categories):
        x = np.arange(len(labels)) + i * bar_width
        axes[0].bar(x, opt_before_data[:, i], width=bar_width, label=category)

    # 设置横坐标和图例
    axes[0].set_xticks(np.arange(len(labels)) + 2 * bar_width)
    axes[0].set_xticklabels(labels)
    axes[0].legend(loc='upper left')
    axes[0].set_title('opt-before')
    axes[0].set_xlabel('target')
    axes[0].set_ylabel('time')

    # 设置y轴范围
    axes[0].set_ylim(0, 5000)

    # 绘制条形统计图 - opt-after
    for i, category in enumerate(categories):
        x = np.arange(len(labels)) + i * bar_width
        axes[1].bar(x, opt_after_data[:, i], width=bar_width, label=category)

    # 设置横坐标和图例
    axes[1].set_xticks(np.arange(len(labels)) + 2 * bar_width)
    axes[1].set_xticklabels(labels)
    axes[1].legend(loc='upper left')
    axes[1].set_title('opt-after')
    axes[1].set_xlabel('target')
    axes[1].set_ylabel('time')

    # 设置y轴范围
    axes[1].set_ylim(0, 5000)

    # 调整子图之间的间距
    plt.tight_layout()

    # 保存为PDF
    
    plt.savefig(path, format='pdf')

    # 显示图形
    plt.show() 
    plt.close()

def drawing_line_chart():

    import matplotlib.pyplot as plt
    import random
    import matplotlib.font_manager as fm
    import pandas as pd 
    df = pd.read_csv('time.csv')

    def drawing_line_chart(x,y,color,name):
        x = x
        y = y    
        plt.plot(x, y,'o-',color=color,label=name)
        # plt.ylim()
        plt.xticks(x)

    names = ['Escherichia_coli', 'Corynebacterium_glutamicum', 'Pseudomonas_aeruginosa', 'Bacillus_subtilis', 'Salmonella_typhimurium']
    colors = ['red', 'orange', 'yellow','green','blue']

    values = []

    for name,color in zip(names,colors):
        strain_df = df[ df['strain'] == name]
        design_target = strain_df['design target']
        time = strain_df['all_time']
        values.append((name, color, design_target, time))

    # 设置斜体字体
    font_properties = fm.FontProperties(style='italic')

    #画图
    plt.figure(figsize=(8, 6))
    for value in values:
        name, color, design_target, time = value
        name = ' '.join( name.split('_') )
        drawing_line_chart(design_target,time,color, name) 

    plt.legend(loc='best')
    plt.xlabel('Design number',fontsize=14)
    plt.ylabel('Time (s)',fontsize=14)
    plt.yticks([200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000])
    plt.savefig('opt-after-target-time.pdf')

    plt.show()
    plt.close()

def extand_mean_to_df(BB_count_no_offTarget_type_means):

    df = pd.DataFrame()
    for i,v in BB_count_no_offTarget_type_means.iterrows():
        dict_temp_li = []
        for i in range(1,6):
            dict_temp = {
                "strain": v["strain"],
                "design target":  v["design target"]+"-"+str(i),
                "B mean:off target numbers": v["B:off target numbers"],
            }
            dict_temp_li.append(dict_temp)
    
        temp = pd.DataFrame( dict_temp_li)
        df = df.append(temp)

    return df

def compute_mean(  B_A_strain_index_df, opt_type="B" ):
    
    B = B_A_strain_index_df[[opt_type, "strain", "design target"]]
    B_count = B.groupby(["strain", "design target"]).count()
    B_count.rename(columns={ opt_type:'B:off target numbers'},inplace=True)
    B.set_index(["strain", "design target"], inplace=True)
    BB_count = pd.merge( B, B_count, left_index=True, right_index=True )
    BB_count = BB_count.reset_index()


    BB_count= BB_count[[ "strain", "design target", "B:off target numbers" ]].drop_duplicates()
    BB_count_no_offTarget_type = BB_count.copy()
    temp = BB_count_no_offTarget_type["design target"].str.split("-", expand=True)
    BB_count_no_offTarget_type["design target"] = temp[0]
    BB_count_no_offTarget_type["design target index"] = temp[1]
    BB_count_no_offTarget_type_means = BB_count_no_offTarget_type.groupby(["strain", "design target"]).mean(f"{opt_type}:off target numbers")
    BB_count_no_offTarget_type_means = BB_count_no_offTarget_type_means.reset_index()

    return BB_count, BB_count_no_offTarget_type_means

def compute_mean_number(path):
    
    a = pd.read_csv(path)
    a.rename(columns={"index":"design target"}, inplace=True)
    B_A = a[["type (B: Before optimization);(A:After optimization)", "strain", "design target"]]["type (B: Before optimization);(A:After optimization)"].str.split(";", expand=True).rename(columns={0:"B", 1:"A"})
    B_A_strain_index = a[[ "strain", "design target"]]
    B_A_strain_index_df = pd.merge(B_A, B_A_strain_index, left_index=True, right_index=True )

    BB_count, BB_count_no_offTarget_type_means = compute_mean(  B_A_strain_index_df, opt_type="B" )
    BB_count_no_offTarget_type_means = extand_mean_to_df(BB_count_no_offTarget_type_means)
    B_offtarget_num_df =  pd.merge( BB_count, BB_count_no_offTarget_type_means,how="outer")
    B_offtarget_num_df = B_offtarget_num_df.fillna(0)

    return B_offtarget_num_df

def correlation_analysis(x,y,draw_path):

    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt_1

    # 使用 numpy 中的 corrcoef 计算相关性矩阵
    correlation_matrix = np.corrcoef(x,y)
    # 提取相关性系数
    correlation_coefficient = correlation_matrix[0, 1]
    print(f"相关性系数：{correlation_coefficient}")
    
    # 绘制散点图
    plt_1.scatter(x,y)
    plt_1.title('numbers of off target vs runtime of off target')
    plt_1.xlabel('off target numbers')
    plt_1.ylabel('runtime(s)')

     # 保存为PDF
    plt_1.savefig(draw_path, format='pdf')

    plt_1.show()
    plt_1.close()

def compute_mean_numbber_and_time(B_offtarget_num_df, opt_before_after_all_df_mean_all, runtime_offtarget_numbers_path, draw_path):
 
    #合并
    opt_before_after_time_and_num = pd.merge( opt_before_after_all_df_mean_all, B_offtarget_num_df )
    opt_before_after_time_and_num["A:time - B:time"] = opt_before_after_time_and_num["A:time"] -  opt_before_after_time_and_num["B:time"] 
    opt_before_after_time_and_num.to_csv(runtime_offtarget_numbers_path)

    #相关性分析
    x,y = opt_before_after_time_and_num['B:off target numbers'], opt_before_after_time_and_num["A:time - B:time"]
    correlation_analysis(x,y,draw_path)

def compute_all_off_target( genome_name_list, output_all_offtarget_path):
        # genome_name_list = [ "Salmonella_typhimurium" ]
    all_df = pd.DataFrame()

    for genome in genome_name_list:  
        for i in range(50,550,50):   
                for step in range(1,6,1):
                    
                        index = str(i)+"-"+str(step)
                        print(genome,index)
                    # genome  = "Salmonella_typhimurium"  

                    # if genome == "Salmonella_typhimurium" and index ==str(i)+"-"+str(step) == "100-3":
                    #     print("lkfjdskl")
                    #     pass
                        # if index == "150-1":
                        #     print("fkjsdhks")
                        # print(genome,index)
                        opt_type = "before"
                        Test_primer_G, Primer_g_offTarget, one_plasmid_design_F_Result = extract_df(genome, index, opt_type)
                        alin_path = extract_alin_path(genome, index, opt_type)

                        invalid_priemr_df, one_plasmid_design_F_Result_id = count_optimal_valid_invalid_fail_nums(one_plasmid_design_F_Result, Test_primer_G, Primer_g_offTarget )
                        opt_befor_invalid_primer = get_invalid_primer(invalid_priemr_df)
                        fail_primer_df = extract_failture_primer(one_plasmid_design_F_Result_id)
                        opt_before_invalid_primer =  pd.merge( opt_befor_invalid_primer, fail_primer_df, how="outer" )
                        opt_type = "Before"
                        opt_before_invalid_primer = rename_invalid_primer(opt_before_invalid_primer, opt_type )
                        opt_before_invalid_primer = extract_opt_invalid_primer(Test_primer_G, opt_before_invalid_primer, alin_path, opt_type,"B:No primer3 design results")

                        opt_type = "after"
                        Test_primer_G, Primer_g_offTarget, one_plasmid_design_F_Result = extract_df(genome, index, opt_type)
                        alin_path = extract_alin_path(genome, index, opt_type)

                        opt_type = "After"
                        invalid_priemr_df, one_plasmid_design_F_Result_id = count_optimal_valid_invalid_fail_nums(one_plasmid_design_F_Result, Test_primer_G, Primer_g_offTarget )
                        opt_after_invalid_primer = get_invalid_primer(invalid_priemr_df)
                        fail_primer_df = extract_failture_primer(one_plasmid_design_F_Result_id)
                        if len( opt_after_invalid_primer ) > 0:
                            opt_after_invalid_primer = rename_invalid_primer(opt_after_invalid_primer, opt_type )
                            opt_after_invalid_primer = extract_opt_invalid_primer(Test_primer_G, opt_after_invalid_primer, alin_path, opt_type, "A:All primer templates offtarget")


                        #################################################################################################################################################
                        if len( opt_after_invalid_primer ) < len( opt_before_invalid_primer ):
                            B_offtarget_A_no_offtarget_df = extract_B_offtarget_A_no_offtarget( opt_before_invalid_primer, opt_after_invalid_primer, Primer_g_offTarget,Test_primer_G, alin_path)
                            opt_after_invalid_primer = opt_after_invalid_primer.append( B_offtarget_A_no_offtarget_df )

                        opt_after_invalid_primer = opt_after_invalid_primer.reset_index(drop=True)

                        df = merge_opt_primer(opt_before_invalid_primer, opt_after_invalid_primer)

                        path = f"{input_dir}/{genome}/{index}"
                        
                        df["strain"] = genome  
                        df["index"] = index 

                        df.to_csv( os.path.join( path,  f"{genome}_{index}.csv") ,index=False)
                  

                   
                        all_df = all_df.append(df)

        print(all_df)
        all_df.to_csv(output_all_offtarget_path, index=False)

def save_input_file(xlsx_file="data_preprocessing/genome_test/all_input.xlsx"):

    genome_name_list  =  ['Escherichia_coli','Corynebacterium_glutamicum', 'Pseudomonas_aeruginosa', 'Bacillus_subtilis', 'Salmonella_typhimurium']
    #
    genome_dict = {}
    for genome in genome_name_list: 
        genome_input_df = pd.DataFrame() 
        for i in range(50,550,50):   
            for step in range(1,6,1):
                index = str(i)+"-"+str(step)
                genome_index_path = f"{input_dir}/{genome}/{index}/{index}.csv"
                temp =  pd.read_csv(genome_index_path)
                temp["index"] = index
                temp["strain"] = genome
                genome_input_df = genome_input_df.append(  temp )
        genome_dict.update( {genome:genome_input_df})

    with pd.ExcelWriter(xlsx_file) as writer:
        for k,v in genome_dict.items():
            v.to_excel(writer,sheet_name = k,index_label='No.')

def save_sample(xlsx_file):

    genome_name_list  =  ['Escherichia_coli','Corynebacterium_glutamicum', 'Pseudomonas_aeruginosa', 'Bacillus_subtilis', 'Salmonella_typhimurium']
    #   
    genome_dict = {}

    for genome in genome_name_list: 
            genome_input_df = pd.DataFrame()
            index = str(500)+"-"+str(1)
            genome_index_path = f"{input_dir}/{genome}/{index}/{index}.csv"
            temp =  pd.read_csv(genome_index_path)
            temp = temp.drop(columns=["position"])
            temp["Sequence upstream of the manipulation site (>100bp)"] = temp["Sequence upstream of the manipulation site (>100bp)"].apply(lambda x: x[-150:])
            genome_dict.update({genome: temp})

    xlsx_file=f"{input_dir}/all_input.xlsx"
    with pd.ExcelWriter(xlsx_file) as writer:
            for k,v in genome_dict.items():
                v.to_excel(writer,sheet_name = k,index_label='No.')

def main():  

    #分析结果的路径
    analysis_path =  "/home/yanghe/public_github/AutoESDCas/runtime/output/"
    if not exists (analysis_path):
        os.makedirs(analysis_path)
    
    #统计所用脱靶数量
    genome_name_list  =  ['Escherichia_coli','Corynebacterium_glutamicum', 'Bacillus_subtilis', 'Salmonella_typhimurium', 'Pseudomonas_aeruginosa']

    output_all_offtarget_path = os.path.join(analysis_path, "opt_before_and_opt_after_sequence_primer.csv")
    compute_all_off_target( genome_name_list, output_all_offtarget_path)   
    # #统计计算优化前，检测出的脱靶的数量
    B_offtarget_num_df = compute_mean_number(output_all_offtarget_path) 

    #统计计算优化前，优化后的时间
    opt_before_after_time_path = os.path.join(analysis_path, "opt_before_after_time.csv")
    columns_index = ['strain', 'design target', 'all_time']
    draw_path = os.path.join(analysis_path, 'whole_process_design.pdf')
    opt_before_data, opt_after_data, opt_before_after_all_df_mean_all = compute_all_time( genome_name_list, columns_index, opt_before_after_time_path, draw_path)

    #合并脱靶数量与时间
    runtime_offtarget_numbers_path =  os.path.join(analysis_path, "runtime_offtarget_numbers.csv")
    draw_path = os.path.join(analysis_path, "runtime_offtargets.pdf")
    compute_mean_numbber_and_time(B_offtarget_num_df, opt_before_after_all_df_mean_all, runtime_offtarget_numbers_path, draw_path)

if __name__ == '__main__':

    main()    



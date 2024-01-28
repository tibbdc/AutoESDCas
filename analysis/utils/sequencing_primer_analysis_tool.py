
#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @FileName     :sequencing_primer_analysis_tool.py
# @Time         :2024/01/27 14:07:19
# @Author       :YangChunhe
# @Email        :2393492851@qq.com
# @Description  :file content
import os 
import pandas as pd
import numpy as np 
import warnings          
warnings.filterwarnings('ignore')
import re,sys  
from os.path import exists
sys.path.append("./utils/")
from algin_tool import algin
   
#测序引物数据统计工具
def get_invalid_primer(invalid_priemr_df):
    li =[]
    for i,v in invalid_priemr_df.iterrows():   

        id = v["ID"].split(":")[1]
     
        off_target = v["off_target"]
        off_target = eval(off_target)
        
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
    if len(li)>0:
        df = pd.DataFrame(li)
    else:
        df = pd.DataFrame(columns=["locus_tag"])
    return df

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

def extract_opt_invalid_primer(Test_primer_G, invalid_primer, alin_path, opt_type, opt_type_value = "A:All primer templates offtarget"):

    fasta_file, tmp_path, genome_path = alin_path
    
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
                    result_df = algin(seq_df,"SEQUENCING_PRIMER_1",fasta_file, tmp_path, genome_path, "0")
                    invalid_primer.loc[i,f"{opt_type} optimization SEQUENCING_PRIMER_1 location (start-end)"] = result_df.loc[0,"After optimization SEQUENCING_PRIMER_1 location (start-end)"]

                if pd.isna(  v[f"{opt_type} optimization SEQUENCING_PRIMER_2 location (start-end)"]  ) or  v[f"{opt_type} optimization SEQUENCING_PRIMER_2 location (start-end)"]=="":
                    

                    right_primer = ab.loc[id,"SEQUENCING_PRIMER_2"]
                    
                    seq_df = pd.DataFrame(
                            [{
                            "ID":id,
                            "SEQUENCING_PRIMER_2":right_primer
                            }]
                            
                        )
                    result_df = algin(seq_df,"SEQUENCING_PRIMER_2",fasta_file, tmp_path, genome_path, "16")
            
                    invalid_primer.loc[i,f"{opt_type} optimization SEQUENCING_PRIMER_2 location (start-end)"] = result_df.loc[0,"After optimization SEQUENCING_PRIMER_2 location (start-end)"]

    return invalid_primer

def count_optimal_valid_invalid_fail_nums(one_plasmid_design_F_Result, Test_primer_G, Primer_g_offTarget ):

    #1.primer3设计出来的所有引物
    Test_primer_G_primer_id = set(Test_primer_G['ID'])

    #2.脱靶的所有引物id
    off_target_primer_id = set(Primer_g_offTarget['ID'].str.split(':',expand=True)[1])
    # off_target_primer_id = set( Primer_g_offTarget['ID'].str.split(':',expand=True)[1]+":"+Primer_g_offTarget['ID'].str.split(':',expand=True)[2] )  #PCR

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
        invalid_priemr_id = set(invalid_priemr_df['ID'].str.split(':',expand=True)[1])
        # invalid_priemr_id = set( invalid_priemr_df['ID'].str.split(':',expand=True)[1] + ":" + invalid_priemr_df['ID'].str.split(':',expand=True)[2] ) #PCR
		
    else:
        invalid_priemr_id = set()
    
    #6.最优引物optimal
    optimal_primer_id = Test_primer_G_primer_id - off_target_primer_id

    #7.有效引物valid
    valid_primer_id = off_target_primer_id - invalid_priemr_id


    print(f"全部引物：{len(all_gene_id)}------最优引物：{len(optimal_primer_id)}-------有效引物:{len(valid_primer_id)}-----------无效引物：{len(invalid_priemr_id)}----------设计失败引物：{len(one_plasmid_design_F_Result_id)}")

    # all_offtarget =  set(invalid_priemr_df[invalid_priemr_df.ID.str.contains('SEQUENCING_PRIMER_1_AND_SEQUENCING_PRIMER_2')]['ID'])   #PCR
    # one_offtarget = set(invalid_priemr_df[~ invalid_priemr_df.ID.str.contains('SEQUENCING_PRIMER_1_AND_SEQUENCING_PRIMER_2')]['ID'] )  #PCR

    one_offtarget = set(invalid_priemr_df[~ invalid_priemr_df.Name.str.contains('SEQUENCING_PRIMER_1_AND_SEQUENCING_PRIMER_2')]['ID']) 
    all_offtarget = set(invalid_priemr_df[invalid_priemr_df.Name.str.contains('SEQUENCING_PRIMER_1_AND_SEQUENCING_PRIMER_2')]['ID'])
	
    success_target = len(optimal_primer_id) + len(valid_primer_id)

    fail_target = len(invalid_priemr_id) + len(one_plasmid_design_F_Result_id)

    print(f"左右引物全部脱靶导致无效引物：{len( all_offtarget )}-------左或右引物脱靶:{len( one_offtarget )}----------primer3设计不出来的引物：{len(one_plasmid_design_F_Result_id)}")


    return invalid_priemr_df, one_plasmid_design_F_Result_id, success_target, fail_target

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

def extract_B_offtarget_A_no_offtarget( opt_before_invalid_primer, opt_after_invalid_primer, Primer_g_offTarget ,Test_primer_G, algin_path):

    fasta_file, tmp_path, genome_path = algin_path 

    A_no_offtarget_B_no_result = set(opt_before_invalid_primer["locus_tag"]) - set(opt_after_invalid_primer["locus_tag"])

    A_no_offtarget_B_no_result_df = pd.DataFrame(columns=["ID"], data = list(A_no_offtarget_B_no_result) )

    Test_primer_G["ID"] = Test_primer_G.ID.apply(lambda x: x.split(";")[0]  )
    A_no_offtarget_B_no_result_primer_df = pd.merge( Test_primer_G, A_no_offtarget_B_no_result_df )

    A_no_offtarget_B_no_result_primer_df = A_no_offtarget_B_no_result_primer_df[["ID","SEQUENCING_PRIMER_1","SEQUENCING_PRIMER_2"]]


    SEQUENCING_PRIMER_1 = algin(A_no_offtarget_B_no_result_primer_df, "SEQUENCING_PRIMER_1",fasta_file, tmp_path, genome_path, "0")
    SEQUENCING_PRIMER_2 = algin(A_no_offtarget_B_no_result_primer_df, "SEQUENCING_PRIMER_2",fasta_file, tmp_path, genome_path, "16")  

    SEQUENCING_PRIMER_1_SEQUENCING_PRIMER_2 = pd.merge(SEQUENCING_PRIMER_1, SEQUENCING_PRIMER_2, how="outer")

    valid_primer_df = Primer_g_offTarget[ Primer_g_offTarget.ID.str.contains('valid_primer') &  (~Primer_g_offTarget.ID.str.contains('invalid_primer')) ]

    SEQUENCING_PRIMER_1_SEQUENCING_PRIMER_2 = get_cor_A_no_offtarget_B_no_result(SEQUENCING_PRIMER_1_SEQUENCING_PRIMER_2, valid_primer_df )

    SEQUENCING_PRIMER_1_SEQUENCING_PRIMER_2 = SEQUENCING_PRIMER_1_SEQUENCING_PRIMER_2.rename(columns={"ID":"locus_tag"})

    SEQUENCING_PRIMER_1_SEQUENCING_PRIMER_2["type (B: Before optimization);(A:After optimization)"] = "A:No offtarget"

    return SEQUENCING_PRIMER_1_SEQUENCING_PRIMER_2
  
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
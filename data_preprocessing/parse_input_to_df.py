#!/usr/bin/env python
# -*- coding:utf-8 -*-
# @FileName     :parse_input_to_df.py
# @Time         :2023/10/27 11:19:24
# @Author       :YangChunhe
# @Email        :2393492851@qq.com
# @Description  :file content  

import pandas as pd
from Bio import SeqIO 
import os,sys
import warnings   
warnings.filterwarnings('ignore')
import configparser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse   
import json
from os.path import exists,splitext,dirname,splitext,basename,realpath,abspath
from Bio import SeqIO


def get_geneID_from_gb( gb_file, feature_type):
    
    """_summary_

    Returns:
        _type_: _description_
    """

    records = SeqIO.parse(gb_file, "genbank")

    gene_id_list = []
    

    for record in records:
        for feature in record.features:

            if feature.type == feature_type:
                gene_id = feature.qualifiers.get("locus_tag", [])
                gene_id = ','.join(gene_id)
                gene_id_list.append(gene_id)

    df = pd.DataFrame()
    df['Gene id'] = gene_id_list
    df.insert(0,'Name', '_del')
    df['Name'] = df['Gene id'] + df['Name']
    df.insert(2,'Up region', 0)
    df.insert(3, 'Down region', 0)
    df.insert(4,'Inserted sequence', '-')
    df.insert(5, 'Manipulation type', 'deletion')


    return df 

def del_Unnamed(df):
    """
    Deletes all the unnamed columns

    :param df: pandas dataframe
    """
    cols_del=[c for c in df.columns if 'Unnamed' in c]
    return df.drop(cols_del,axis=1)

def gb_2_fna(gb_file,fna_file):

    # 读取 GenBank 文件
    gb_records = SeqIO.parse(gb_file, "genbank")

    # 将 DNA 序列写入 FASTA 文件
    
    with open(fna_file, "w") as output_handle:
        SeqIO.write(gb_records, output_handle, "fasta")

def get_seq_by_geneid(gb_file, gene_id, up_region, down_region):
    
    # 提取基因序列和其位置
    gene_id = gene_id  # 将YOUR_GENE_ID替换为你要提取的基因ID
    print('取序列的地址：' ,gb_file)
    records = SeqIO.parse(gb_file, "genbank")
    
    ref, mutation_pos_index, chrom, strand = '', '', '', ''

    for record in records:

        for feature in record.features:
            if feature.type == "gene" and gene_id in feature.qualifiers.get("locus_tag", []):
                #取序列
                gene_seq = feature.extract(record.seq)
                if feature.strand == -1:  # 如果是反义链，则需要取反向互补序列
                    gene_seq = gene_seq.reverse_complement()
                
                #取出染色体上的基因的位置
                gene_start = feature.location.start.position

                gene_end = feature.location.end.position
                start = gene_start + up_region
                end = gene_end - down_region

                if start == end:
                    ref = '-'

                elif int(up_region) !=0 and int(down_region) == 0:
                    ref = '-'

                elif end > start:

                    print(gene_seq)
                    gene_seq = str(gene_seq)
                    region = f'{0+up_region} : {len(gene_seq)-down_region}' 
                    ref = gene_seq[0+up_region : len(gene_seq)-down_region ]
                    print(start,end)

                else:
                    ValueError('There is an issue with the editing location you provided')

                chrom = record.id
                strand = 'plus'
                mutation_pos_index = start      

                # #取目标序列的前100bp
                # start_pos = max(gene_start - 100, 0)
                # before_gene_seq = record.seq[start_pos:gene_start]
                
                # print(f">{gene_id} {record.id}:{gene_start}-{gene_end}\n{ref}")
                break

        if ref != '':
            break

    if ref == '':  
        error_message = f'{gene_id}:{up_region},{down_region}:The region you will edit is not on the target genome'
        raise ValueError(error_message)

    print("取到的东西：",ref, mutation_pos_index, chrom, strand)
    return ref, mutation_pos_index, chrom, strand

def revComp(seq):
    complementSeq=seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))
    revcompSeq = complementSeq[::-1]
    return revcompSeq

def conf_read(filename): 
    config = configparser.ConfigParser()
    config.read(filename)
    res = dict(config._sections["point"])
    return res

def input_to_primer_template(input_file_path, genome, workdir, scene, gb_file=''):
    """
    Arguments:
        input_file_path[str]:  input information
        genome[str]:  edit reference
        config[dic]: points information   
    Return [dict]
        {   
            "key1":{  
                "seq_uha_max_whole":"",
                "seq_dha_max_whole":"",
                "seq_altered":"",
                "type":"",   # [substitution,deletion,insertion]
                "ref":"",
                "uha_upstream": seq_uha_max_whole  up 100bp  sequence,
                "dha_downstream":seq_dha_max_whole  down 100bp sequence,
            },
            "key2":{
                "seq_uha_max_whole":"",
                "seq_dha_max_whole":"",
                "seq_altered":"",  
                "type":"",
                "ref":"",
                "uha_upstream": seq_uha_max_whole  up 100bp  sequence,
                "dha_downstream":seq_dha_max_whole  down 100bp sequence,
            }
        }
    """
    input_format=input_file_path.split('.')[-1]
    with open(input_file_path,'r') as f:
        input_header = f.readlines()[0]
        # print(input_header)

    primer_template = {}

    # blast error file
    with open(os.path.join(workdir,'error.txt'),'w') as blast_error_handler:
        blast_error_handler.write('ID\tERROR\n')
        if input_format != 'csv':
            error_message = "The input file format needs to be 'csv'"
            print(error_message)
            return error_message
        elif 'Sequence upstream' in input_header:
            # input type 1: upstream
            print('processing upstream input file ...')
            record_dict = SeqIO.to_dict(SeqIO.parse(genome, "fasta"))
            blast_search_dict = blast_search(input_file_path,genome,workdir)
            df=pd.read_csv(input_file_path)

            num_lines_df=df.shape[0]
            for i in range(num_lines_df):
                data=df.loc[i].values
                #print(data)
                mun_id=str(data[0])
                if len(data) < 2:
                    error_message = "Some necessary input information is missing in the line "+mun_id+" of the input file"
                    raise ValueError(error_message)
                    print(error_message)
                    return  error_message
                else:
                    upstream = data[1].strip().upper()
                    ref = data[2]
                    name = mun_id
                    if scene == 'both_sgRNA_primer' or scene == 'only_primer':
                        alt = data[3]
                        mutation_type = data[4].strip().lower()
                       
                    if mun_id not in blast_search_dict:
                        # no blast
                        error_message = "The upstream sequence of " + mun_id + " can not be mapped to the target genome. please check whether the target sequence is located on the target genome."
                        blast_error_handler.write(mun_id + '\t' +  error_message+ '\n')
                        raise ValueError(error_message)
                        return  error_message
                    elif blast_search_dict[mun_id]["unique_mapped"] > 1:
                        # Compare 100 times
                        error_message = "The upstream sequence of " + mun_id + "  can be mapped to multiple loci in the target genome, %s, Please provide a longer upstream seqeunce." % blast_search_dict[mun_id]["description"]
                        blast_error_handler.write(mun_id+'\t'+ error_message+'\n')
                        raise ValueError(error_message)
                        return  error_message
                    elif blast_search_dict[mun_id]["unique_mapped"] == 0:
                        # No 100 comparison
                        error_message = "The upstream sequence of " + mun_id + " can not be uniquely mapped to the target genome. Please check whether the target sequence is located on the target genome."
                        blast_error_handler.write(mun_id+'\t'+error_message+'\n')
                        raise ValueError(error_message)
                        return  error_message
                    elif blast_search_dict[mun_id]["unique_mapped"] == 1:
                        # Index of the base starting to mutate on gene
                        if blast_search_dict[mun_id]["reverse_mapped"]:
                            record = revComp(str(record_dict[blast_search_dict[mun_id]["chrom"]].seq))
                            upstream_start_index = len(record) - int(blast_search_dict[mun_id]["start"])
                            strand = "-"  
                            error_message = "The upstream sequence of " + mun_id + "can be mapped to the antisense strand.  Please rightly prepare input file for target manipulation as the example of 2,3-BD" 
                            raise ValueError(error_message)
                            return  error_message
                        else:
                            record = str(record_dict[blast_search_dict[mun_id]["chrom"]].seq)   
                            upstream_start_index = int(blast_search_dict[mun_id]["start"])-1
                            strand = "+"
                        chrom = blast_search_dict[mun_id]["chrom"]
                        mutation_pos_index = upstream_start_index + len(upstream)

                        # get mutation info dict 

                        if scene == 'both_sgRNA_primer' or scene == 'only_primer':
                            res = create_mutation_info(mutation_pos_index,strand,chrom,name,ref,mutation_type,alt,record,mun_id)
                        elif scene == 'only_sgRNA':
                            if ref != '-':
                                genome_ref = record[mutation_pos_index:mutation_pos_index+len(ref)]  
                                if genome_ref.upper() == ref.upper():
                                    res =  {
                                        "ref":ref,
                                        "strand":"plus" if strand =="+" else "minus",
                                        "mutation_pos_index":mutation_pos_index,
                                        "geneid":chrom,
                                        "name":name,
                                        "region":chrom+ ':' +  str(mutation_pos_index) +'-'+ str(int(mutation_pos_index)+len(ref))
                                    }
                                else:
                                    error_message = "The target mutation ref of " + mun_id + " can not be found in reference, please check."
                                    raise ValueError(error_message)
                                    return error_message  
                            else:
                                res =  {
                                        "ref":ref,
                                        "strand":"plus" if strand =="+" else "minus",
                                        "mutation_pos_index":mutation_pos_index,
                                        "geneid":chrom,
                                        "name":name,
                                        "region":chrom+ ':' +  str(mutation_pos_index) +'-'+ str(int(mutation_pos_index)+len(ref))
                                    }
                            
                        if isinstance(res,str):
                            blast_error_handler.write(mun_id+'\t'+res+'\n')
                        else:
                            primer_template[mun_id] = res

  
        elif 'Chr,Pos,Strand' in input_header:
            # input type 2: vcf
            print('processing vcf input file ...')
            record_dict = SeqIO.to_dict(SeqIO.parse(genome, "fasta"))
            df=pd.read_csv(input_file_path)
            num_lines_df=df.shape[0]
            for i in range(num_lines_df):
                data=df.loc[i].values
                #print(data)
                mun_id=str(data[0])
                if len(data) < 7:
                    error_message = "Some necessary input information is missing in the line "+mun_id+" of the input file"
                    print(error_message)
                    return  error_message
                else:
                    chrom = data[1].strip()
                    pos = data[2]
                    strand = data[3]
                    ref = data[4].strip()
                    alt = data[5]
                    mutation_type = data[6]
                    name = mun_id

                    if strand == '+':
                        record = str(record_dict[chrom].seq)
                        mutation_pos_index = int(pos) - 1
                    elif strand == '-':
                        record = revComp(str(record_dict[chrom].seq))
                        mutation_pos_index = len(record) - int(pos)
                    # get mutation info dict
                    res = create_mutation_info(
                        record,mutation_type,mutation_pos_index,
                        ref,alt,strand,chrom,
                        name,mun_id
                        )
                    if isinstance(res,str):
                        blast_error_handler.write(mun_id + '\t'+ res +'\n')
                    else:
                        primer_template[mun_id] = res


        elif 'Gene id,Up region,Down region' in input_header:
            #input type 3: region

            df=pd.read_csv(input_file_path)
           
            for i,v in df.iterrows():
                    
                ref, mutation_pos_index, chrom, strand = get_seq_by_geneid(gb_file, v['Gene id'], v['Up region'], v['Down region'])


                if 'Inserted sequence,Manipulation type' in input_header:
                        mu_type = v['Manipulation type']
                        seq_altered = v['Inserted sequence']
                        name = v['Name']  
                        res =  {
                                "name":name,
                                "ref":ref,
                                "strand":strand,
                                "mutation_pos_index":mutation_pos_index,
                                "geneid":chrom,
                                "region":chrom+ ':' +  str(mutation_pos_index) +'-'+ str(int(mutation_pos_index)+len(ref)),
                                "type": mu_type,
                                "seq_altered":seq_altered
                                }
                else:
                        name = v['Name']   
                        res =  {
                                "name":name,
                                "ref":ref,
                                "strand":strand,
                                "mutation_pos_index":mutation_pos_index,
                                "geneid":chrom,
                                "region":chrom+ ':' +  str(mutation_pos_index) +'-'+ str(int(mutation_pos_index)+len(ref))
                                }
                    
                primer_template[name] = res
        
        else:
            error_message = "The input file format not supported, Please rightly prepare input file for target manipulation as the example of 2,3-BD."
            raise ValueError(error_message)
            return  error_message
    return primer_template

def blast_search(input_file_path,genome,workdir):
    """
    Arguments:
        input_file_path[str]:  input information
        genome[str]:  edit reference
        workdir[str]: dir path
    Return[dict]
        {
            "key1":{
                "chrom":"NC_001133.9",
                "start":"1000",
                "unique_mapped":1,  100% comparison
                "mapped":1,  frequency
                "reverse_mapped":1, frequency
                "description":"NC_001133.9:26529-26708;NC_001133.9:26664-26843;",  
            },
        }

    """
    blast_output_file_path=os.path.join(workdir,'blast_search.txt')
    ref_lib=genome.split('/')[-1].split('.')[0]
    input_fasta = os.path.join(workdir,'blast.fasta')
    fasta_length_dict = {}
    my_records = []
    with open(input_file_path,'r') as ifile:
        index = 0 
        for line in ifile:
            if index != 0 :
                linelist = line.split(',')
                seqid = linelist[0]
                seq = linelist[1].strip()
                rec = SeqRecord(Seq(seq),id=seqid)

                fasta_length_dict[seqid] = len(seq)
                my_records.append(rec)
            index += 1
    # input fasta
    SeqIO.write(my_records, input_fasta, "fasta")
    # run blast
    os.system("makeblastdb -in "+genome+" -dbtype nucl -parse_seqids -out "+ref_lib)
    os.system("blastn -query "+input_fasta+" -db "+ref_lib+" -outfmt 6 -task blastn -out "+blast_output_file_path+" -evalue 1e-30 ")
    os.system("rm %s.n*" % ref_lib)

    # return
    dictall = {}
    with open(blast_output_file_path,"r") as f:
        for i in f:
            linelist = i.split('\t')
            key = linelist[0]
            chrom = linelist[1]
            identity = linelist[2]
            allength = linelist[3]
            start = linelist[8]
            end = linelist[9]
            
            if key not in dictall:
                dictall[key] = {
                    "chrom":chrom,
                    "start":start,
                    "unique_mapped":1 if (int(float(identity)) == 100 and int(float(allength))==fasta_length_dict[key]) else 0,
                    "mapped":1,
                    "reverse_mapped":1 if (int(start) > int(end) and int(float(identity)) == 100 and int(float(allength))==fasta_length_dict[key]) else 0,
                    "description":'%s:%s-%s;' %(chrom,start,end),
                }
            else:
                dictall[key]["mapped"] += 1
                if int(float(identity)) == 100 and int(float(allength))==fasta_length_dict[key] :
                    dictall[key]["unique_mapped"] += 1
                    dictall[key]["chrom"] = chrom
                    dictall[key]["start"] = start
                    if int(start) > int(end):
                        dictall[key]["reverse_mapped"] += 1
                    dictall[key]["description"] += '%s:%s-%s;' %(chrom,start,end)
    return dictall

def create_mutation_info(mutation_pos_index,strand,chrom,name,ref,mutation_type,alt,record,mun_id):
    """  
    Arguments:
        record[str]: mutation site located genome
        mutation_type[str]: deletion insertion substitution
        mutation_pos_index[int]:  mutation_pos_index information
        ref[str]:  mutation_pos_index ref
        alt[str]: alter sequence
        strand[str]:  + plus   - minus
        chrom[str]:  mutation site located genome name
        name[str]: mutation name
        mun_id[str]: mutation id
        config[dict]:  information
    Return [dict]
        {
            "seq_uha_max_whole":"",
            "seq_dha_max_whole":"",
            "seq_altered":"",
            "type":"",   # [substitution,deletion,insertion]
            "ref":"",
            "uha_upstream": seq_uha_max_whole  up 100bp  sequence,
            "dha_downstream":seq_dha_max_whole  down 100bp sequence,
        }   
    
    """   
    if mutation_type == "insertion":
        info_dict = {
            "seq_altered":alt,
            "type":mutation_type,
            "ref":ref,
            "strand":"plus" if strand =="+" else "minus",
            "mutation_pos_index":mutation_pos_index,
            "geneid":chrom,
            "name":name,
            "region":chrom+ ':' +  str(mutation_pos_index) +'-'+ str(int(mutation_pos_index)+len(ref))
        }
        return info_dict      
#   return info_dict
    if mutation_type in ["deletion","substitution"]:
        genome_ref = record[mutation_pos_index:mutation_pos_index+len(ref)]
        
        # print(genome_ref)
        if genome_ref.upper() == ref.upper():
            info_dict = {
                "seq_altered":"-" if alt=='-' else alt,
                "type":mutation_type,
                "ref":ref,
                "strand":"plus" if strand == "+" else "minus",
                "mutation_pos_index":mutation_pos_index,
                "geneid":chrom,
                "name":name,   
                "region": chrom+ ':' +  str(mutation_pos_index) +'-'+ str(int(mutation_pos_index)+len(ref))
            }
            return info_dict
        else:
            error_message = "The target mutation ref of " + mun_id + " can not be found in reference, please check."
            raise ValueError(error_message)
            return error_message
    else:
        error_message = "The target manipulation type of " + mun_id + " must be equal to 'insertion,substitution or deletion', Please rightly prepare input file for target manipulation as the example of 2,3-BD."
        raise ValueError(error_message)
        return  error_message
    
def dict_to_df(dict_input_seq):
    info_input_df = pd.DataFrame()
    for item in dict_input_seq:
        df = pd.DataFrame([dict_input_seq[item]])
        # df.insert(loc=0,value=item,column='id')
        info_input_df = info_input_df.append(df)
    info_input_df = info_input_df.reset_index(drop=True)
    return info_input_df

def execute_input_2_chopchop_input(input_file_path,  genome_path, convert_input_file_chopchopInput_workdir, chopchop_input, scene, gb_file=''):

    before_info_input_df = pd.read_csv(input_file_path)
    before_info_input_df.columns = [i.lower() for i in before_info_input_df.columns]

    print("基因组:",genome_path)  
    dict_input_seq = input_to_primer_template(input_file_path, genome_path, convert_input_file_chopchopInput_workdir, scene, gb_file)
      
    info_input_df = dict_to_df(dict_input_seq)
    if scene == 'only_primer':
        info_input_df = pd.merge(before_info_input_df[['name','crrna']], info_input_df, on='name', how='inner')
    info_input_df.to_csv(chopchop_input,index=False)
   
def main(data): 

    #1.读取参数路径
    genome_path = data['ref_genome']
    path,name = splitext(genome_path)

    convert_input_file_chopchopInput_workdir = data['data_preprocessing_workdir']
    input_file_path = data['input_file_path']
    scene = data['scene']

    #2.生成workdir
    if not os.path.exists(convert_input_file_chopchopInput_workdir):
        os.makedirs(convert_input_file_chopchopInput_workdir)

    #3.生成成下游任务的输入文件路径
    chopchop_input = os.path.join(  
        data['data_preprocessing_workdir'],
        'info_input.csv'
    )

    #4.若上传的是基因组gb文件，生成fna文件
    if 'gb' in name:
        genome_fna = 'xxx' + '.fna'
        gb_file = genome_path
        genome_path = os.path.join(data['data_preprocessing_workdir'], genome_fna)
        #解析gb文件生成fna文件
        gb_2_fna(gb_file, genome_path)

        #5.根据不同场景，解析输入的编辑信息文件生成下游任务的标准输入
        print('场景：',scene)
        execute_input_2_chopchop_input(input_file_path, genome_path, convert_input_file_chopchopInput_workdir, chopchop_input, scene, gb_file)
        return [chopchop_input, genome_path]
    else:
        #5.根据不同场景，解析输入的编辑信息文件生成下游任务的标准输入
        print('场景：',scene)
        execute_input_2_chopchop_input(input_file_path, genome_path, convert_input_file_chopchopInput_workdir, chopchop_input, scene)
        
        genome_name = 'xxx'+ '.fna'
        genome_path1 = os.path.join(data['data_preprocessing_workdir'], genome_name)

        os.system(f"cp {genome_path} {genome_path1}")

        return chopchop_input

if __name__ == '__main__':
    call_method = 1  

    if  call_method == 1:
        data1 = {
                    "input_file_path":"/home/yanghe/program/data_preprocessing/input/editor_info.csv",
                    "ref_genome":"/home/yanghe/program/data_preprocessing/input/GCA_000011325.1_ASM1132v1_genomic.fna",
                    "data_preprocessing_workdir":"/home/yanghe/tmp/data_preprocessing/output/",
                    "scene":"only_sgRNA",  
                }
        data2 = {
                    "input_file_path":"/home/yanghe/program/data_preprocessing/input/editor_info123.csv",
                    "ref_genome":"/home/yanghe/program/data_preprocessing/input/GCA_000011325.1_ASM1132v1_genomic.fna",
                    "data_preprocessing_workdir":"/home/yanghe/tmp/data_preprocessing/output/",
                    "scene":"both_sgRNA_primer",
                }
        data3 = {
                    "input_file_path":"./input/sgRNA_editing_info.csv",
                    "ref_genome":"./input/GCA_000011325.1_ASM1132v1_genomic.fna",
                    "data_preprocessing_workdir":"/home/yanghe/tmp/data_preprocessing/output/",
                    "scene":"only_primer",  
                }
        
        data4 = {
                    "input_file_path":"./input/4-21-input.csv",
                    "ref_genome":"./input/eco.gb",
                    "data_preprocessing_workdir":"/home/yanghe/tmp/data_preprocessing/output/",
                    "scene":"only_sgRNA",
                }  
        data5 = {
                    "input_file_path":"/home/yanghe/program/data_preprocessing/input/05-10-input.csv",   
                    "ref_genome":"/home/yanghe/program/data_preprocessing/input/GCF_000005845.2_ASM584v2_genomic.gbff",
                    "data_preprocessing_workdir":"/home/yanghe/tmp/data_preprocessing/output/",
                    "scene":"both_sgRNA_primer",
                }
        data6 = {
                "input_file_path":"./input/4-23-input.csv",
                "ref_genome":"./input/GCF_000005845.2_ASM584v2_genomic.gbff",    
                "data_preprocessing_workdir":"/home/yanghe/tmp/data_preprocessing/output/",
                "scene":"only_primer",
            }

        #大肠全敲
        data7 = {
                    "input_file_path":"/home/yanghe/program/data_preprocessing/input/eco_all_knock_out.csv",
                    "ref_genome":"/home/yanghe/program/data_preprocessing/input/eco.gb",
                    "data_preprocessing_workdir":"/home/yanghe/tmp/data_preprocessing/output/",
                    "scene":"both_sgRNA_primer",
                }
        #大肠rna敲除
        data8 = {
                    "input_file_path":"/home/yanghe/program/data_preprocessing/input/10-13-input.csv",
                    "ref_genome":"/home/yanghe/program/data_preprocessing/input/GCF_000005845.2_ASM584v2_genomic.gbff",
                    "data_preprocessing_workdir":"/home/yanghe/tmp/data_preprocessing/output/",
                    "scene":"both_sgRNA_primer",
                }
          
        data = data8

    elif call_method == 2:
        parser = argparse.ArgumentParser()
        parser.add_argument('--input', '-i', help='input config file', required=True)   
        arguments = parser.parse_args()
        input_file_path = arguments.input

        with open(input_file_path,'r',encoding='utf8') as fp:
            data = json.load(fp)
 
    import time
    time1=time.time()
    a = main(data2)
    print(a)

    time2=time.time()
    print('耗时：'+str(time2-time1)+'s')     


  
 
    



      
    
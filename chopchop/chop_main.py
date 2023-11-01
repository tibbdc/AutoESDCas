# -*- coding:utf-8 -*-
# @FileName     :chop_main.py
# @Time         :2023/10/27 15:21:18
# @Author       :YangChunhe
# @Email        :2393492851@qq.com
# @Description  :file content

import os,sys
from os import makedirs  
import pandas as pd
from os.path import exists,splitext,dirname,splitext,basename,realpath,abspath
os.getcwd()
import json  
import argparse
import math   
import subprocess
from Bio import SeqIO
import multiprocessing as mp



import subprocess

def conda_env_list():
    command = ['conda', 'env', 'list']
    result = subprocess.run(command, capture_output=True, text=True)
    output_lines = result.stdout.strip().split('\n')[2:]  # 忽略前两行

    env_dict = {}
    for line in output_lines:
        parts = line.split()
        if len(parts) >= 2:
            name = parts[0]
            path = parts[1]
            env_dict[name] = path

    return env_dict


def excecute_one_chopchop(env,chopchop_params,parent_output):
    
    env = env
    base_path = base_path = os.path.abspath(os.path.dirname(__file__)) + '/'
    genome_name = chopchop_params['genome_name']
    output = chopchop_params['output']
    info_list_one = chopchop_params['info_list_one']
    PAM = chopchop_params['PAM']
    scoringMethod = chopchop_params['scoringMethod']
    maxMismatches = chopchop_params['maxMismatches']
    guideSize = chopchop_params['guideSize']

    if int(guideSize) > 40 or int(guideSize) < 5:
        raise ValueError('Please enter a guideSize between 5 and 40')
    if int(maxMismatches) < 1  or int(maxMismatches) > 4:
        raise ValueError('Please enter a maxMismatches between 1 and 4')

    cmd = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s" % \
                                    (env + '/bin/python',
                                        base_path+'/'+'chopchop.py',
                                        '-G',
                                        genome_name,
                                        '-o', 
                                        output,
                                        '-Target',
                                        info_list_one,   
                                        '-scoringMethod',     # ["XU_2015", "DOENCH_2014", "DOENCH_2016", "MORENO_MATEOS_2015", "CHARI_2015", "G_20", "KIM_2018", "ALKAN_2018", "ZHANG_2019", "ALL"]
                                        scoringMethod,
                                        '-v',    #[0, 1, 2, 3]
                                        maxMismatches,
                                        '-M',    #The PAM motif
                                        PAM,
                                        '-g',
                                        guideSize,
                                        '-config',
                                        parent_output+'/'+'config.json'
                                    )   
    # logger.info(f"cmd = {cmd}")  
    
    print('执行的命令:',cmd)   


    runbashcmd(cmd)  
    try:
        temp_df = pd.read_csv(output+'result.csv')
    except Exception as e:
        # logger.error(f'Get files list failed. | Exception: {e}')
        print(e)
        temp_df = pd.DataFrame()
    
    temp_df['Region'] = info_list_one
    return temp_df


def call_chopchop(output,chopchop_params,info_list_one, env):   
  
    temp_output = output +'/temp/'+info_list_one + '/'
    if not exists(temp_output):
        makedirs(temp_output)

    chopchop_params.update({'output':temp_output,"info_list_one":info_list_one})

    #统计计算时间
    import time
    start_time = time.time()   
   
    temp_df = excecute_one_chopchop(env,chopchop_params,output)

    end_time = time.time()
    arr = info_list_one.split(":")[1].split("-")
    temp_df['edit_region_length'] = int(arr[1]) - int(arr[0])
    temp_df['time(s)'] =  end_time - start_time

    return temp_df



def check_create_path(config):

    genome_name = config['chopchop_params']['genome_name']
    genome_path = config['genome_path']

    #2bit
    bit_datat = config['bit_datat']
    genome_2bit ='{}{}{}'.format(bit_datat,genome_name,'.2bit')
    #ebwt
    ebwt_datat = config['ebwt_datat']
    genome_ebwt_datat = '{}{}'.format(ebwt_datat,genome_name)

      #judge
    if not exists(bit_datat):  
        makedirs(bit_datat)
    if not exists(genome_2bit):
        cmd = '{}\t{}\t{}'.format('faToTwoBit',genome_path,genome_2bit)
        runbashcmd(cmd)
    if not exists(genome_ebwt_datat):
        makedirs(genome_ebwt_datat)

        print(genome_ebwt_datat)
        genome_ebwt_datat = genome_ebwt_datat + '/'+ genome_name
        cmd = '{}\t{}\t{}'.format('bowtie-build',genome_path,genome_ebwt_datat)
        print('build:-----------------------------------------------',cmd)
        runbashcmd(cmd)


def write_config(config,input_path,ref_genome,chopchop_params,output, bowtie='bowtie', twoBitToFa='twoBitToFa'):

    chopchop_workdir = output

    #-------------------------------------------
    config['genome_path'] = ref_genome
    config['chopchop_params'] = chopchop_params
    config['info_path'] = input_path 
    config['output_path'] = output


    #-------------------------------------------

    #-----------------------------------------
    if chopchop_workdir in config['PATH']['BOWTIE_INDEX_DIR']: 
        pass
    else:
        config['PATH']['BOWTIE'] = bowtie
        config['PATH']['TWOBITTOFA'] = twoBitToFa
        config['PATH']['TWOBIT_INDEX_DIR'] = chopchop_workdir +'/'+ config['PATH']['TWOBIT_INDEX_DIR']
        config['PATH']['PRIMER3'] = chopchop_workdir +'/'+ config['PATH']['PRIMER3']
        config['PATH']['BOWTIE_INDEX_DIR'] = chopchop_workdir +'/'+ config['PATH']['BOWTIE_INDEX_DIR'] +'/'+ config['chopchop_params']['genome_name']
    #------------------------------------------

    #--------------------------------------------
    config['ebwt_datat'] = chopchop_workdir + '/' + config['ebwt_datat']  
    config['bit_datat'] = chopchop_workdir + '/' + config['bit_datat']
    #--------------------------------------------


    check_create_path(config)
    


    #write new_config.json
    json_str = json.dumps(config)
   
    with open(chopchop_workdir+"/"+'config.json', 'w') as json_file:
        json_file.write(json_str)
    

    return chopchop_workdir+"/"+'config.json'


def runbashcmd(cmd,test=False,logf=None):

    # from super_beditor.lib.global_vars import dirs2ps 
    print('before:',cmd)
    # cmd = cmd.replace("$BIN", dirs2ps['binDir'])
    # cmd = cmd.replace("$PYTHON", dirs2ps['pyp'])
    # cmd = cmd.replace("$SCRIPT", dirs2ps['scriptDir'])
    print('after:',cmd)
    if test:
        print(cmd)
    
    err=subprocess.call(cmd,shell=True,stdout=logf,stderr=subprocess.STDOUT)
    print(err)
    if err!=0:
        print('bash command error: {}\n{}\n'.format(err,cmd))
        sys.exit(1)


def extract_seq_from_genome(genome,gene_id,start=0,end=0):
    # 读取基因组文件
    records = SeqIO.parse(genome,'fasta')
    
    # 遍历基因组文件中的所有记录
    for record in records:
        # 如果当前记录的ID与所需的ID匹配
        if record.id == gene_id:
            # 提取坐标范围内的序列
            if start ==0 and end ==0:
                return str(record.seq)
            else:
                sequence = str(record.seq[start:end])
                if sequence != '':
                    return  sequence.upper()
                else:
                    return sequence

  
def filter(df,genome):  
    df = df[df['region'] != '']
    df['chromosome_seq_len'] = df.region.apply(lambda x: len(extract_seq_from_genome(genome=genome, gene_id=x.split(":")[0])))
    def work(region, chromosome_seq_len):
            target_coords=region.split(":")[1]
            target_coords = int(target_coords.split('-')[0]),int(target_coords.split('-')[1])
            start = target_coords[0]
            end = target_coords[1]
            if start < 80 or chromosome_seq_len - end < 80 or start >= end:
                return 'no'
            else:
                return 'yes'
    df['tag'] = df.apply(lambda x: work(x['region'], x['chromosome_seq_len']), axis=1)

    df = df[df['tag'] == 'yes']
    return df


# 定义一个函数，用于在多进程中调用
def process_region(params):
    output,chopchop_params, region = params
    # 调用call_chopchop函数
    print(output,chopchop_params, region)   
    result = call_chopchop(output, chopchop_params, region)
    return result



#生成前端需要的json
def sgRNAdf_to_jsion(sgRNA):
    li = []
    def work(x):
    #     print(x['Region'])
        name = list(x['Name'])[0]
       
        cc = list(x.T.to_dict().values())
        temp_dict = {}
        temp_dict.update({'Name':name,'Detail':cc})
        li.append(temp_dict)
    sgRNA.groupby(by='Region').apply(lambda x: work(x))
    return li

  
def main(event):

    #
    base_path = os.path.abspath(os.path.dirname(__file__)) + "/chopchop/"
 
    print(base_path)
    print(os.listdir(base_path))  

    #env
    if os.getenv('CONDA_PREFIX') == None:
        env = event.get('env')
    else:
        env = os.environ['CONDA_PREFIX']
        if os.path.basename(  env ) != 'chopchop':
            env = conda_env_list()['chopchop']
            print(env)
        
    #parse event
    input_path = event["input_file_path"]
    ref_genome = event["ref_genome"] 
    chopchop_params = event["chopchop_config"]
    output = event["chopchop_workdir"] 
    chopchop_params['genome_name'] = splitext(basename(ref_genome))[0]

    #read config
    chopchop_local_config = base_path + '/data/input/' + 'config.json'

    with open(chopchop_local_config, "r") as f:
        data = json.load(f)

    #write config to temp
    tem_config = write_config(data, input_path, ref_genome, chopchop_params, output)
    print('temp中的config位置',tem_config)

    #
    df = pd.read_csv(input_path)
    #过滤无用信息
    df = filter(df,ref_genome)   

    print(df.columns) 
  
    # info_list_one = 'CM001534.1:231516-232905'      
    # sgRNA = call_chopchop(output,base_path,chopchop_params,info_list_one)
    # ------------------------------------------------------
    import time
    start_time = time.time()
    #parallel processing  
    
    temp = df.region.apply(lambda x: call_chopchop(output,chopchop_params,x, env))


    # pool = mp.Pool()
    #   # 使用map函数在多个进程中调用process_region函数
    # results = pool.map(process_region, [(output,chopchop_params, region) for region in df.region])

    # # 关闭进程池
    # pool.close()

    sgRNA_df = pd.concat([row for i,row in temp.items()])
    sgRNA_df = pd.merge(df[['name','region']],sgRNA_df,left_on=['region'],right_on=['Region'],how='inner')
    sgRNA_df.drop(columns=['region'],inplace=True)
    sgRNA_df = sgRNA_df.rename(columns={'name':'Name'})
    # sgRNA_df.rename(columns={'Name':'ID'},inplace=True)

    sgRNA_df = sgRNA_df[['Name','Rank',"Region",'Genomic location','Strand','Target sequence',"GC content (%)","MM0","MM1","MM2","MM3","Efficiency"]]  

    #输出json
    li_sgRNA = sgRNAdf_to_jsion(sgRNA_df)  
    with open(output+"/"+"sgRNA.json", 'w') as input_json:
        json.dump(li_sgRNA, input_json, indent=4, ensure_ascii=False)

    sgRNA_df.to_csv(output+'/'+'sgRNA.csv',index=False)   
    end_time = time.time()
    print("Execution time:", end_time - start_time, "seconds")       
    # --------------------------------------------------------  


    import shutil
    # remove .sam files as they take up wayyy to much space
    print(output)
    if exists(output +'/'+ 'temp'):
        shutil.rmtree(output +'/'+ 'temp')
    return output + '/' + 'sgRNA.csv'  


if __name__ == '__main__':
    # parser = argparse.ArgumentParser()
    # parser.add_argument('--input', '-i', help='input params file', required=True) 
    # args = parser.parse_args()
    # input_path =  args.input  

    "XU_2015", "DOENCH_2014", "DOENCH_2016", "MORENO_MATEOS_2015", "CHARI_2015", "G_20"


    # "env": "/home/yanghe/anaconda3/envs/crispr_hr_editor/"
  
    call_method = 1
    if  call_method == 1:  
        event1 = {
            "input_file_path":"/home/yanghe/tmp/data_preprocessing/output/info_input.csv",
            "ref_genome":"/home/yanghe/program/data_preprocessing/input/GCA_000011325.1_ASM1132v1_genomic.fna",
            "chopchop_workdir":"/home/yanghe/tmp/chopchop/output/", 
            "chopchop_config":{
                "PAM": "NGG", 
                "guideSize": 20,
                "maxMismatches": 3,
                "scoringMethod": "XU_2015"
            }
        }
        event2 = {
            "input_file_path":"/home/yanghe/tmp/data_preprocessing/output/info_input.csv",
            "ref_genome":"/home/yanghe/tmp/data_preprocessing/output/xxx.fna",
            "chopchop_workdir":"/home/yanghe/tmp/chopchop/output/", 
            "chopchop_config":{
                "PAM": "NGG", 
                "guideSize": 20,
                "maxMismatches": 3,
                "scoringMethod": "DOENCH_2014"
            }
        }
        event = event2
    
    elif call_method == 2:
        parser = argparse.ArgumentParser()
        parser.add_argument('--input', '-i', help='input config file', required=True)   
        arguments = parser.parse_args()
        input_file_path = arguments.input

        with open(input_file_path,'r',encoding='utf8') as fp:
            event = json.load(fp)
    


    a=main(event)
    print(a)

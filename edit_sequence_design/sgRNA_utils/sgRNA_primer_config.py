
# -*- coding:utf-8 -*-
# @FileName     :sgRNA_primer_config.py
# @Time         :2023/11/01 12:46:44
# @Author       :YangChunhe
# @Email        :2393492851@qq.com
# @Description  :file content

import os

class config: 

    TMP = '' 

    #单双点设计引物的全局参数  
    GLOBAL_ARGS = {
                'PRIMER_PICK_ANYWAY':1,
                'PRIMER_PRODUCT_SIZE_RANGE': 0,
                'PRIMER_NUM_RETURN':1000
        }
   
    S_GLOBAL_ARGS = {  
            'PRIMER_OPT_SIZE': 20,   
            'PRIMER_MIN_SIZE': 18,
            'PRIMER_MAX_SIZE': 25,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 35.0,
            'PRIMER_MAX_TM': 85.0,
            'PRIMER_MIN_GC': 30.0,
            'PRIMER_OPT_GC':55,
            'PRIMER_MAX_GC': 70.0,
    }
    UHA_ARGS = {
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 55.0,
        'PRIMER_MAX_TM': 75.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_OPT_GC':55,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_MIN_SIZE':18,
        'PRIMER_MAX_SIZE':25,
        'PRIMER_OPT_SIZE': 20
    }
    SEQ_ALTERED_ARGS = {
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 55.0,
        'PRIMER_MAX_TM': 75.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_OPT_GC':55,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_MIN_SIZE':18,
        'PRIMER_MAX_SIZE':25,
        'PRIMER_OPT_SIZE': 20
    }
    DHA_ARGS = {
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 55.0,
        'PRIMER_MAX_TM': 75.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_OPT_GC':55,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_MIN_SIZE':18,
        'PRIMER_MAX_SIZE':25,
        'PRIMER_OPT_SIZE': 20
    }
    UP_SGRNA_ARGS = {
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 55.0,
        'PRIMER_MAX_TM': 75.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_OPT_GC':55,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_MIN_SIZE':18,
        'PRIMER_MAX_SIZE':25,
        'PRIMER_OPT_SIZE': 20
    }
    DOWN_SGRNA_ARGS = {
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 55.0,
        'PRIMER_MAX_TM': 75.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_OPT_GC':55,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_MIN_SIZE':18,
        'PRIMER_MAX_SIZE':25,
        'PRIMER_OPT_SIZE': 20
    }

    #测序引物设计的全局参数
    Q_ARGS = {
            'PRIMER_OPT_SIZE': 20,
            'PRIMER_MIN_SIZE': 18,  
            'PRIMER_MAX_SIZE': 25,
            'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 55.0,  
            'PRIMER_MAX_TM': 75.0,    
            'PRIMER_MIN_GC': 20.0,
            'PRIMER_OPT_GC':55,
            'PRIMER_MAX_GC': 80.0        
    }

    PLASMID_Q_ARGS = {
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 55.0,  
        'PRIMER_MAX_TM': 75.0,    
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_OPT_GC':55,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_MIN_SIZE': 18,  
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_SIZE': 20
    }
    GENOME_Q_ARGS = {
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 55.0,  
        'PRIMER_MAX_TM': 75.0,    
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_OPT_GC':55,
        'PRIMER_MAX_GC': 50.0,
        'PRIMER_MIN_SIZE': 18,  
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_SIZE': 20
    }

    UHA_DHA_CONFIG = {
        "max_right_arm_seq_length": 0,  
        "max_left_arm_seq_length": 150,   
        "min_left_arm_seq_length": 1,   
        "min_right_arm_seq_length": 150   
    }
    UHA_DHA_LENGTH = {
        'uha':145,
        'dha':145  
    }

    PRIMER_TEMPLATE_EXTEND = 2000


    workdir = './'
    tmp_path = './'    


    PLASMID_LABEL = {
        'ccdb_label':'',
        'promoter_terminator_label':'',
        'n_20_label':'',
        'promoter_label':''
    }



    
    Q_GLOBAL_ARGS = {   
                'PRIMER_PICK_ANYWAY':1,    
                'PRIMER_TASK': 'pick_primer_list', 
        }  

    #提取target_gene_seq和marker_seq需要的常量
    AMPLICONIC_GENE_TARGET_SEQ_LENGTH = 0
    
    #提取marker时需要的常量
    AMPLICONIC_MARKER_SEQ_LENGTH = 0

    INPUT_FILE_PATH = ''   
    
    PLASMID_FILE_NAME = 'xxx.gb' 

    INPUT_FILE_NAME = ''  
    NEW_OUTPUT_FILE_PATH=''
    OUTPUT_ORDERS_NAME = 'order.xlsx'  
    INPUT_IMG_PATH = '/input/'
    IMG_ORDERS_NAME = 'tib.png'

    SEQUENCING_PRIMER_SUCCESS = 'sequencing_primer_success.xlsx'
    SEQUENCING_PRIMER_FAILTRUE = 'sequencing_primer_failtrue.xlsx'

    OUTPUT_FILE_NAME_PLASMID_MUTATION = "plasmid_mutation.gb" 
    DATA_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..")) + '/data'

    OUTPUT_FILE_PATH = DATA_ROOT + '/output/'

    BOWTIE_PATH = ""




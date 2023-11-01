import re
import pandas as pd
import primer3
from sgRNA_utils.sgRNA_primer_config import config      
from Bio import SeqIO
from os.path import exists,splitext,dirname,splitext,basename,realpath,abspath 


#检查DNA序列中非规范碱基的位置
def find_non_agct_positions(sequence):
    non_standard_bases = [m.start() for m in re.finditer('[^AGCTagct]', sequence)]
    return non_standard_bases

def del_Unnamed(df):
    """
    Deletes all the unnamed columns

    :param df: pandas dataframe
    """
    cols_del=[c for c in df.columns if 'Unnamed' in c]
    return df.drop(cols_del,axis=1)

def dfswapcols(df,cols):                #交换列
    df[f"_{cols[0]}"]=df[cols[0]].copy()
    df[cols[0]]=df[cols[1]].copy()
    df[cols[1]]=df[f"_{cols[0]}"].copy()
    df=df.drop([f"_{cols[0]}"],axis=1)
    return df

#取互补链
def revComp(seq):
    complementSeq=seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))
    revcompSeq = complementSeq[::-1]
    return revcompSeq

#发现字符串的位置
def find_coordinate_by_pattern(pattern,seq):
    sub_seq_cor = {}
    i = 0
    for m in re.finditer(pattern, seq):
        sub_seq_cor.update({f'{i}':(m.start(), m.end())})
        i = i + 1

    for m in re.finditer(pattern, revComp(seq)):
        sub_seq_cor.update({f'{i}':(m.start(), m.end())})
        i = i + 1
    return sub_seq_cor,i

#使id位于第一列
def make_id_in_first_columns(df,id,columns):
    assert id in columns
    first_columns_list=df[id].tolist()
    df.drop(columns = [id],inplace=True)
    df.insert(loc=0,column =id ,value = first_columns_list)
    return df

def read_excel(p,sheet_name=None,):
    xl = pd.ExcelFile(p)
    xl.sheet_names  # see all sheet names
    if sheet_name is None:
        sheet_name=input(', '.join(xl.sheet_names))
    return xl.parse(sheet_name) 

def to_excel(sheetname2df,datap,):
    writer = pd.ExcelWriter(datap)
    for sn in sheetname2df:
        sheetname2df[sn].to_excel(writer,sn)
    writer.save() 
   
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

#将一个dataframe拆分列，变成堆加行
def columns_2_row_from_one_df(sgRNA_df,in_col,to_col):
    in_col1 = [in_col[0],in_col[1],in_col[2]]
    in_col2 = [in_col[0],in_col[1],in_col[3]]
    uha_primer_template = sgRNA_df[in_col1].rename(columns={in_col[2]:'primer_template'})
    dha_prime_template = sgRNA_df[in_col2].rename(columns={in_col[3]:'primer_template'})
    uha_primer_template['type'] = 'uha'
    dha_prime_template['type'] = 'dha'
    primer_template = pd.concat([uha_primer_template,dha_prime_template])
    return primer_template

#创建uha_dha
def columns_2_row_by_groupby(uha_dha_primer_df,in_col,ou_col,type='u_d'):
    print(uha_dha_primer_df.columns)
    
    def work(x):
        uha = x[x['Type']=='uha'].reset_index(drop=True)
        dha = x[x['Type']=='dha'].reset_index(drop=True)

        if type == 'u_d':
            if len(uha)>0 and len(dha)>0:
                data = [[   uha.loc[0,in_col[0]],
                            uha.loc[0,in_col[1]],
                            uha.loc[0,in_col[2]],
                            uha.loc[0,in_col[3]],
                            dha.loc[0,in_col[2]], 
                            dha.loc[0,in_col[3]]
                        ]]
                U_D_df = pd.DataFrame(columns=ou_col,data=data)
            else:
                U_D_df = pd.DataFrame()  
            return U_D_df
            
        elif type == 'primer':
            data = [[    uha.loc[0,in_col[0]], 
                         uha.loc[0,in_col[1]],
                         uha.loc[0,in_col[2]],
                         uha.loc[0,in_col[3]],
                         uha.loc[0,in_col[4]],
                         uha.loc[0,in_col[5]],

                         dha.loc[0,in_col[2]],
                         dha.loc[0,in_col[3]],
                         dha.loc[0,in_col[4]],
                         dha.loc[0,in_col[5]]
                    ]]
            u_primer_df = pd.DataFrame(columns=ou_col,data=data)
            return u_primer_df
    UHA_DHA_df = uha_dha_primer_df.groupby([in_col[0],in_col[1]]).apply(lambda x: work(x))   
    UHA_DHA_df.reset_index(drop=True,inplace=True)
    return UHA_DHA_df   


#读取质粒的特征坐标     
def get_feature_coordinate(target_gene_label,gb_path, selected_feature_type='misc_feature'):
    gb = SeqIO.read(gb_path, "genbank") 

    for ele in gb.features:
        if ele.type == selected_feature_type:
            if target_gene_label in ele.qualifiers.get('label'):  
                target_gene_start = ele.location.start.position
                target_gene_end = ele.location.end.position
                print( target_gene_start,target_gene_end )
                break
            else:
                target_gene_start = -1
                target_gene_end = -1
                print( target_gene_start,target_gene_end )
    return int(target_gene_start),int(target_gene_end)

#换名字
def replace_primer3Name_to_peopleReadName(df,type=''):
    names = df.columns.to_list()
    df =df.rename(columns={
                            names[2]:f"primer_{type}_f_seq_(5'-3')",
                            names[3]:f"primer_{type}_r_seq_(5'-3')",
                            names[4]:f"primer_{type}_f_Tm",
                            names[5]:f"primer_{type}_r_Tm",
                            names[1]:f"{type}_product_sequence_length",
                            names[0]:f"{type}_product_sequence"
                        } )
    return df  

#设计引物   primer_design(seqId = 'first_last', seqTemplate = plasmid_backbone, primer_type='plasmid',stype = 'left_right')
def primer_design(seqId,
                  seqTemplate,
                  stype,
                  mute_type='single',
                  primer_type = '',
                  global_args=config.GLOBAL_ARGS,
                  ):  

    #检查引物设计模板序列
    non_agct_positions = find_non_agct_positions(seqTemplate)
    if len( non_agct_positions ) != 0:
        raise ValueError(f'Template is error, position is {non_agct_positions}')

    if mute_type == 'single':
        config.GLOBAL_ARGS.update(config.S_GLOBAL_ARGS)
        global_args = config.GLOBAL_ARGS
        if primer_type == 'uha':
            global_args.update(config.UHA_ARGS)
            config.UHA_DHA_CONFIG['max_left_arm_seq_length'] == config.UHA_DHA_CONFIG['min_left_arm_seq_length']
            stype = 'left_right'
            print(global_args)
        elif primer_type == 'dha':
            global_args.update(config.DHA_ARGS)
            config.UHA_DHA_CONFIG['max_left_arm_seq_length'] == config.UHA_DHA_CONFIG['min_left_arm_seq_length']
            stype = 'left_right'

        elif primer_type == 'plasmid':
            global_args.update(config.DOWN_SGRNA_ARGS)  
        elif primer_type == 'sgRNA':
            global_args.update(config.UP_SGRNA_ARGS)
        elif primer_type == 'seq_altered':
            global_args.update(config.SEQ_ALTERED_ARGS)

    elif mute_type == 'sequencing': 
        config.Q_GLOBAL_ARGS.update(config.Q_ARGS)
        global_args =config.Q_GLOBAL_ARGS 
        if primer_type == 'plasmid_seq':
            global_args.update(config.PLASMID_Q_ARGS)
        elif primer_type == 'genome_seq':
            global_args.update(config.GENOME_Q_ARGS)
 
    #序列参数
    seqlength = len(seqTemplate)   
    seq_args = {
                'SEQUENCE_ID': seqId,
                'SEQUENCE_TEMPLATE': seqTemplate,
                'SEQUENCE_FORCE_LEFT_START':0,
                'SEQUENCE_FORCE_RIGHT_START': seqlength-1
        }  
     
    #选择类型，设定序列参数
    if mute_type == 'single': #单点突变
        if stype == "left":   #target上游引物设计
            seq_args['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = [[0,50,-1,-1]]
            seq_args['SEQUENCE_FORCE_LEFT_START'] = -1
            size_range = [seqlength-50, seqlength] 
        elif stype == "right":
            seq_args['SEQUENCE_FORCE_RIGHT_START'] = -1  
            seq_args['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = [[-1,-1, seqlength-50,50]]
            size_range = [seqlength-50, seqlength]                  

        elif stype == "left_right":   #target下游引物设计  
            seq_args['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = [[-1,-1,-1,-1]]  
            seq_args['SEQUENCE_FORCE_LEFT_START'] = 0  
            seq_args['SEQUENCE_FORCE_RIGHT_START'] = seqlength-1       
            size_range = [seqlength,seqlength]

        elif stype == "left_right_change":
            seq_args['SEQUENCE_PRIMER_PAIR_OK_REGION_LIST'] = [[0,40,seqlength-40,40]]
            seq_args['SEQUENCE_FORCE_LEFT_START'] = -1
            seq_args['SEQUENCE_FORCE_RIGHT_START'] = -1
            size_range = [seqlength-80,seqlength]  
            # primer_num = 5 

       
    elif mute_type == 'sequencing':  #测序引物
        seq_args['SEQUENCE_ID'] = seqId
        seq_args['SEQUENCE_TEMPLATE'] = seqTemplate  
        size_range = [25,seqlength]
        # primer_num = 5

    #设定全局参数   
    global_args['PRIMER_PRODUCT_SIZE_RANGE']=size_range   
    if primer_type == 'plasmid_seq':
        global_args['PRIMER_NUM_RETURN'] = 10    
    elif primer_type == 'genome_seq':
        global_args['PRIMER_NUM_RETURN'] = 1000
    # print(global_args['PRIMER_NUM_RETURN'])
    #调用工具
    print(seq_args,global_args)
    primer3_result = primer3.bindings.designPrimers(seq_args, global_args) 
         
    return primer3_result  

#输出引物设计成功的
def result_output_success_df(plasmid_name,primers_dict,type=''):
    all_df = pd.DataFrame()
    
    for key in primers_dict:
        df =pd.DataFrame(primers_dict[key])    
        df['id'] = key   
        all_df=all_df.append(df)  
    #换名子
    all_df = replace_primer3Name_to_peopleReadName(all_df,type)  
    columns_name = all_df.columns.to_list()
    # #列的顺序
    all_df = all_df[[columns_name[-1],columns_name[2], columns_name[3],columns_name[4],columns_name[5],columns_name[0],columns_name[1]]]
    all_df['template']=plasmid_name
    return all_df  

#输出设计引物失败  
def result_output_failtrue_df(plasmid_name,primers_dict_failture):

    df = pd.DataFrame()
    for k,v in primers_dict_failture.items():
        temp = pd.DataFrame([v])
        temp.insert(0,'id',k)
        df = df.append(temp)
    df['template']=plasmid_name
    return df  

#构建酶库
def create_enzyme(path):
    enzyme_df = pd.DataFrame(data=[['BsaI','GGTCTC',1,4,'CCA'],
                                   ['BbsI','GAAGAC',2,4,'CCA'],
                                   ['SapI','GCTCTTC',1,3,'CCA'],
                                   ['BsmBI','CGTCTC',1,4,'CCA'],
                                   ['BbsI','GAAGAC',2,4,'CCA'],
                                  ], columns=['enzyme','recognition_seq','gap_len','cut_seq_len','protective_base'])
    enzyme_df.to_csv(path,index=False)
  
#标准化重命名,uha_primer_df,dha_primer_df,sgRNA_primer_df,  plasmid_sequencing_primer_df,genome_sequencing_primer_df  
def normal_primer_rename(uha_primer_df,dha_primer_df,sgRNA_primer_df,plasmid_sequencing_primer_df,genome_sequencing_primer_df,plasmid_backbone_primer_df,seq_altered_primer_df):
    print(plasmid_backbone_primer_df.columns)   

    uha_primer_df = uha_primer_df.rename(columns={'Region':'ID',
                                              "primer_f_seq_(5'-3')":"PRIMER_LEFT_WHOLE_SEQUENCE",
                                              "primer_r_seq_(5'-3')":"PRIMER_RIGHT_WHOLE_SEQUENCE",
                                              "UHA":"PRODUCT_SEQUENCE",
                                              "UHA_size":"PRODUCT_WHOLE_LENGTH",
                             })

    dha_primer_df = dha_primer_df.rename(columns={'Region':'ID',
                                                "primer_f_seq_(5'-3')":"PRIMER_LEFT_WHOLE_SEQUENCE",
                                                "primer_r_seq_(5'-3')":"PRIMER_RIGHT_WHOLE_SEQUENCE",
                                                "DHA":"PRODUCT_SEQUENCE",
                                                "DHA_size":"PRODUCT_WHOLE_LENGTH",
                                })
    
    plasmid_backbone_primer_df = plasmid_backbone_primer_df.rename(columns={
                                              "primer_f_seq_(5'-3')_joint":"PRIMER_LEFT_WHOLE_SEQUENCE",
                                            "primer_r_seq_(5'-3')_joint":"PRIMER_RIGHT_WHOLE_SEQUENCE",
                                            "product_value_joint":"PRODUCT_SEQUENCE",
                                            "product_size_joint":"PRODUCT_WHOLE_LENGTH"   
                                })
   
    sgRNA_primer_df=sgRNA_primer_df[['Region',"primer_f_seq_(5'-3')_joint","primer_r_seq_(5'-3')_joint","product_value_joint","product_size_joint"]]
    sgRNA_primer_df = sgRNA_primer_df.rename(columns={'Region':'ID',
                                                "primer_f_seq_(5'-3')_joint":"PRIMER_LEFT_WHOLE_SEQUENCE",
                                                "primer_r_seq_(5'-3')_joint":"PRIMER_RIGHT_WHOLE_SEQUENCE",
                                                "product_value_joint":"PRODUCT_SEQUENCE",
                                                "product_size_joint":"PRODUCT_WHOLE_LENGTH",
                                })
    
    seq_altered_primer_df = seq_altered_primer_df[['Region',"primer_f_seq_(5'-3')_joint","primer_r_seq_(5'-3')_joint","product_value_joint","product_size_joint"]]
    seq_altered_primer_df = seq_altered_primer_df.rename(columns={'Region':'ID',
                                                "primer_f_seq_(5'-3')_joint":"PRIMER_LEFT_WHOLE_SEQUENCE",
                                                "primer_r_seq_(5'-3')_joint":"PRIMER_RIGHT_WHOLE_SEQUENCE",
                                                "product_value_joint":"PRODUCT_SEQUENCE",
                                                "product_size_joint":"PRODUCT_WHOLE_LENGTH",
                                })

    plasmid_sequencing_primer_df = plasmid_sequencing_primer_df.rename(columns={'Region':'ID'})
    plasmid_sequencing_primer_df.reset_index(drop=True, inplace=True)
    genome_sequencing_primer_df = genome_sequencing_primer_df.rename(columns={'Region':'ID'})
    genome_sequencing_primer_df.reset_index(drop=True, inplace=True)  
       
    return uha_primer_df,dha_primer_df,sgRNA_primer_df,plasmid_sequencing_primer_df,genome_sequencing_primer_df,plasmid_backbone_primer_df,seq_altered_primer_df

#先填充a的内容，再合并两个df 
def merge_fill_two_df(temp_sgRNA_df,a):
    #填充
    a_df = pd.DataFrame()
    for i in range(len(temp_sgRNA_df)):
        a_df = a_df.append(a)
    a_df.reset_index(drop=True,inplace=True)
    #合并
    sgRNA_plasmid_primer_joint_df = pd.concat([temp_sgRNA_df,a_df],axis=1)
    return sgRNA_plasmid_primer_joint_df

# columns
def df_to_df(columns,df,index):

    r_li = []
    f_li = []
    for item in columns:
        if 'r_r' in item:
            r_li.append(item)
        elif 'r_f' in item:
            f_li.append(item)
    r_li = sorted(r_li)
    f_li = sorted(f_li)

    primer_li = []
    i = 0
    for f,r in zip(f_li,r_li):
        i = i + 1
        if 'Region' in df.columns:
            primer_dict = {
                            "Region":df.loc[index,'Region'],    
                            "primer_f_seq_(5'-3')_joint": df.loc[index,f],
                           "primer_r_seq_(5'-3')_joint": df.loc[index,r]
                          }
        else:
            primer_dict = {
                            "Region": i,    
                            "primer_f_seq_(5'-3')_joint": df.loc[index,f],
                           "primer_r_seq_(5'-3')_joint": df.loc[index,r]
                          }
        primer_li.append(primer_dict)
    df = pd.DataFrame(primer_li)
    return df

#重命名
def rename_u_d_primer_df(uha_primer_df,dha_primer_df):
    uha_primer_df = uha_primer_df.rename(columns={'Region':'ID',
                                              "primer_f_seq_(5'-3')":"PRIMER_LEFT_WHOLE_SEQUENCE",
                                              "primer_r_seq_(5'-3')":"PRIMER_RIGHT_WHOLE_SEQUENCE",
                                              "UHA":"PRODUCT_SEQUENCE",
                                              "UHA_size":"PRODUCT_WHOLE_LENGTH"
                             })

    dha_primer_df = dha_primer_df.rename(columns={'Region':'ID',
                                                "primer_f_seq_(5'-3')":"PRIMER_LEFT_WHOLE_SEQUENCE",
                                                "primer_r_seq_(5'-3')":"PRIMER_RIGHT_WHOLE_SEQUENCE",
                                                "DHA":"PRODUCT_SEQUENCE",
                                                "DHA_size":"PRODUCT_WHOLE_LENGTH"
                                })
    return uha_primer_df, dha_primer_df

def rename_sequencing_primer_df(*df_list):
    df_list_renamed = []
    for df in df_list:
        df = df.rename(columns={'Region':'ID'})
        df.reset_index(drop=True, inplace=True)
        df_list_renamed.append(df)
    return df_list_renamed

def result_primer_list_to_df(primer_result_list):
    df = pd.DataFrame()
    for item_dict in primer_result_list:
        temp = pd.DataFrame([item_dict])
        df = df.append(temp)
    return df

def rename_common_primer_df(*df_list):
    df_list_renamed = []
    for df in df_list:
        if 'index' not in df.columns:  
            df = df.rename(columns={'Region':'ID',
                                "primer_f_seq_(5'-3')_joint":"PRIMER_LEFT_WHOLE_SEQUENCE",
                                "primer_r_seq_(5'-3')_joint":"PRIMER_RIGHT_WHOLE_SEQUENCE",
                                "product_value_joint":"PRODUCT_SEQUENCE",
                                "product_size_joint":"PRODUCT_WHOLE_LENGTH",
                    })
        elif 'index' in df.columns:
             df = df.rename(columns={
                                'index':"Index",
                                'Region':'ID',
                                "primer_f_seq_(5'-3')_joint":"PRIMER_LEFT_WHOLE_SEQUENCE",
                                "primer_r_seq_(5'-3')_joint":"PRIMER_RIGHT_WHOLE_SEQUENCE",
                                "product_value_joint":"PRODUCT_SEQUENCE",
                                "product_size_joint":"PRODUCT_WHOLE_LENGTH",
                    })
        df_list_renamed.append(df)
    return df_list_renamed

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

#只要不存在off-target，可义将一切序列转化成坐标
def convert_seq_cor(gb_path,region_seq_json,strand='+',seq=''):
    if seq != '':
        gb_seq = seq.upper()
        print('给定序列，坐标转成给定序列上的坐标')
    else:
        gb = SeqIO.read(gb_path, "genbank")
        gb_seq = str(gb.seq).upper()
    
    if strand == '-': 
        gb_seq = revComp(gb_seq)
    
    region_cor_json = {}
    for k,v in region_seq_json.items():
        v = v.upper()
        start = gb_seq.find(v)
        end = start + len(v)     
        if start == -1:
            i = 1
            start_seq = ''
            end_seq = ''
            while True:
                new_start = gb_seq.find(v[:i])
                if new_start == -1:
                    break
                else:
                    if gb_seq.find(v[i:]) == 0:
                        start_seq = v[:i]
                        end_seq = v[i:]
                        break
                    else:        
                        i = i + 1
                        if i == len(v):
                            break
            if start_seq == '' and end_seq == '':
                start = -1
                end = -1
            else:
                start = gb_seq.find(start_seq)
                end = gb_seq.find(end_seq)  
                end = end + len(v)
                
        region_cor_json.update({k:f'{start},{end}'})
        
    return region_cor_json  

def check_seq_in_gb(gb_path, seq_json):
    
    cor_json_plus = convert_seq_cor(gb_path,seq_json,strand='+')
    cor_json_min = convert_seq_cor(gb_path,seq_json,strand='-')
    
    for plus,minus in zip(cor_json_plus,cor_json_min):
        if cor_json_plus[plus] =='-1,-1' and  cor_json_min[minus] == '-1,-1':
            #序列不存在GB文件
            seq_json[plus] = 'The sequence you provided is not on the plasmid file'
        elif cor_json_plus[plus]!='-1,-1' and cor_json_min[minus] == '-1,-1':
            #序列在正义链
            pass
        elif cor_json_plus[plus]!='-1,-1' and cor_json_min[minus] == '-1,-1':
            #序列在负义链
            seq_json[plus] = revComp(seq_json[plus])
        else:
            pass
    return seq_json  

def get_joint_by_enzyme(enzyme_df,enzyme_name):
    sgRNA_enzyme_df = enzyme_df[enzyme_df['name']==enzyme_name]
    sgRNA_enzyme_df = sgRNA_enzyme_df.reset_index(drop=True)
    protective_base = sgRNA_enzyme_df.loc[0,'protective_base']
    recognition_seq = sgRNA_enzyme_df.loc[0,'recognition_seq']
    cut_seq_len = sgRNA_enzyme_df.loc[0,'cut_seq_len']
    gap_len = sgRNA_enzyme_df.loc[0,'gap_len']
    joint_len = len(protective_base) + len(recognition_seq) + gap_len + cut_seq_len
    return joint_len, cut_seq_len

def write_gb(record,output_path,gb_name,gb_type='genbank'):
    SeqIO.write(record, output_path+'/'+gb_name+'.gb' ,gb_type)

def columns_to_row(n20down_primer_p_df):
    n20down_primer_p_df = n20down_primer_p_df.sort_values(by='ID')
    temp = n20down_primer_p_df[['PRIMER_LEFT_WHOLE_SEQUENCE','PRIMER_RIGHT_WHOLE_SEQUENCE']]
    index = [str(i[0])+'_'+str(i[1]) for i in temp.stack().to_frame().index]
    temp = temp.stack().to_frame()
    temp.index = index
    temp = temp.T
    return temp

def groupby_columns_to_row(df):
    def work(x):
        x = x.reset_index(drop=True)
        id = list(x['ID'])[0]
        temp = x[['PRIMER_LEFT_WHOLE_SEQUENCE','PRIMER_RIGHT_WHOLE_SEQUENCE']]    
        temp = temp.stack().to_frame()
        temp = temp.T
        temp.columns = [str(i[0]) + '_' + i[1] for i in temp.columns] 
        temp.insert(0,'ID',id)
        return temp
    sg = df.groupby(by='ID').apply(lambda x: work(x))
    sg = sg.reset_index(drop=True)
    return sg  

def create_primerCor_in_plasmid(plasmid_seq,primer):
    plasmid_seq = plasmid_seq.upper()
    primer = primer.upper()
    start = plasmid_seq.find(primer)
    end = start + len(primer)
    return start, end
   
def create_plasmid_primer_featrue_df(sequencing_primer_template,
                                     uha_primer_df,
                                     seq_altered_p_df,
                                     dha_primer_df,
                                     n20up_primer_p_df=pd.DataFrame()):
    
   
    if 'promoter_N20_terminator' in sequencing_primer_template.columns:
        sequencing_plasmid = sequencing_primer_template[['Region','plasmid','UHA','DHA','promoter_N20_terminator']].rename(columns={'Region':'ID','plasmid':'PLASMID','promoter_N20_terminator':'PROMOTER_N20_TERMINATOR'})
    else:
        if len(seq_altered_p_df)>0:
            sequencing_plasmid = sequencing_primer_template[['Region','plasmid','UHA','DHA','seq_altered']].rename(columns={'Region':'ID','plasmid':'PLASMID','seq_altered':'SEQ_ALTERED'})  
        else:
            sequencing_plasmid = sequencing_primer_template[['Region','plasmid','UHA','DHA']].rename(columns={'Region':'ID','plasmid':'PLASMID'})     
       
    uha_primer_df = uha_primer_df[['ID','PRIMER_LEFT_WHOLE_SEQUENCE','PRIMER_RIGHT_WHOLE_SEQUENCE']]
    uha_primer_df.columns = ['ID','UHA_PRIMER_LEFT_WHOLE_SEQUENCE','UHA_PRIMER_RIGHT_WHOLE_SEQUENCE']  
    
    if len(seq_altered_p_df) > 0 :
        seq_altered_p_df = seq_altered_p_df[['ID','PRIMER_LEFT_WHOLE_SEQUENCE','PRIMER_RIGHT_WHOLE_SEQUENCE']]
        seq_altered_p_df.columns = ['ID','SEQ_PRIMER_LEFT_WHOLE_SEQUENCE','SEQ_PRIMER_RIGHT_WHOLE_SEQUENCE']
    
    dha_primer_df = dha_primer_df[['ID','PRIMER_LEFT_WHOLE_SEQUENCE','PRIMER_RIGHT_WHOLE_SEQUENCE']]
    dha_primer_df.columns = ['ID','DHA_PRIMER_LEFT_WHOLE_SEQUENCE','DHA_PRIMER_RIGHT_WHOLE_SEQUENCE']
    
    if len(n20up_primer_p_df) != 0:
        n20up_primer_p_df = n20up_primer_p_df[['ID','PRIMER_LEFT_WHOLE_SEQUENCE','PRIMER_RIGHT_WHOLE_SEQUENCE']]
        n20up_primer_p_df.columns = ['ID','N20UP_PRIMER_LEFT_WHOLE_SEQUENCE','N20UP_PRIMER_RIGHT_WHOLE_SEQUENCE']
        df = pd.merge(sequencing_plasmid,uha_primer_df).merge(dha_primer_df).merge(n20up_primer_p_df)

        if len(seq_altered_p_df) > 0 :
            df = pd.merge(seq_altered_p_df,df,how='outer')
    else:
        df = pd.merge(sequencing_plasmid,uha_primer_df).merge(dha_primer_df)
        if len(seq_altered_p_df) > 0 :
            df = pd.merge(seq_altered_p_df,df,how='outer')    
    return df

#根据质粒特征label取序列
def get_sequence_by_feature_label(filename, feature_label):
    for record in SeqIO.parse(filename, "genbank"):
        for feature in record.features:
            if 'label' in  feature.qualifiers and feature.qualifiers['label'][0] == feature_label :
                print(feature.location.extract(record).seq)
                return feature.location.extract(record).seq
    return None

import re
def is_dna(sequence):
    # 将序列中的所有字母转换为大写字母
    sequence = sequence.upper()
    # 使用正则表达式匹配是否只包含A、T、G、C四种字符
    pattern = re.compile(r'^[ATGC]+$')
    if pattern.match(sequence):
        return True
    else:
        return False
    
# #打包文件为zip
import zipfile 
import os
def zip_ya(startdir, file_news,num=1):

    z = zipfile.ZipFile(file_news, 'a', zipfile.ZIP_DEFLATED)
    for dirpath, dirnames, filenames in os.walk(startdir):
        fpath = dirpath.replace(startdir, '')  # 这一句很重要，不replace的话，就从根目录开始复制
        fpath = fpath and fpath + os.sep or ''
        for filename in filenames:  
            if num == 1:
                if filename =='success.xlsx' or filename == 'failtrue.xlsx':
                    z.write(os.path.join(dirpath, filename), fpath + filename,zipfile.ZIP_DEFLATED)
            else:
                fpath=startdir.split('/')[-2]+'/'
                z.write(os.path.join(dirpath, filename), fpath + filename,zipfile.ZIP_DEFLATED)

#
def convert_twoColumns_to_oneColumns(df,id,region,name1,name2,name):
    
    temp1 = df[[id,region,name1]].rename(columns={name1:name})
    temp1[id] = temp1[id] +';'+ temp1[region]+ f';{name1}'
    temp2 = df[[id,region,name2]].rename(columns={name2:name})
    temp2[id] = temp2[id] +';'+ temp1[region] + f';{name2}'
    temp = temp1.append(temp2)  
    temp.drop(columns=['Region'], inplace=True)
    
    return temp   

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def df_to_fasta(df,fasta_filename):
    with open(fasta_filename, 'w') as fasta_file:
        for index, row in df.iterrows():
            sequence_id = row['id']
            sequence = row['sequence']
            fasta_file.write(f'>{sequence_id}\n{sequence}\n')




def bowtie_seq_genome(tmp_path, path, index_prefix, fasta_filename):

    if not exists(tmp_path):
        os.makedirs(tmp_path)

    sam_path = os.path.join(tmp_path ,'output.sam')  

    if not exists(index_prefix):
        os.makedirs(index_prefix)
        # index = os.path.join(index_prefix, 'genome_index')
        # os.chmod(f'{config.BOWTIE_PATH}bowtie-build', 0o755)
        cmd = f'{config.BOWTIE_PATH}/bowtie-build {path} {index_prefix}'
        os.system(cmd)  
    # os.chmod(f'{config.BOWTIE_PATH}bowtie', 0o755)
    cmd = f'{config.BOWTIE_PATH}/bowtie -p 2 -v 3 --sam-nohead -k 1000000000 {index_prefix} -f {fasta_filename} -S {sam_path}'
    os.system(cmd)

    #解析  
    sam_df = parse_sam_file(sam_path)

    #删除文件

    return sam_df




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


def gb_2_fna(gb_file,fna_file):

    # 读取 GenBank 文件
    gb_records = SeqIO.parse(gb_file, "genbank")

    # 将 DNA 序列写入 FASTA 文件
    with open(fna_file, "w") as output_handle:
        SeqIO.write(gb_records, output_handle, "fasta")




def convert_df_to_fastaFile(genome_path,df,id,name,input_fasta,lib_fasta):

    record_dict = SeqIO.to_dict(SeqIO.parse(genome_path, "fasta"))
    
     
    my_records = []
    fasta_length_dict = {}

    lib_records = []
    lib_length_dict = {}

    #根据gene的序列找到前后各20k

    for i,row in df.iterrows():

        seqid = row[id]
        seq = row[name]
        chrom = seqid.split(';')[1].split(':')[0]
        start,end = seqid.split(';')[1].split(':')[1].split('-')
        gene_start, gene_end, target_lib_seq = search_sequence_in_genome_seq(chrom, genome_path, start, end,  start_distance=20000, end_distance=20000)
        if target_lib_seq != None:
            rec = SeqRecord(Seq(target_lib_seq),id=seqid)
            lib_length_dict[seqid] = len(target_lib_seq)
            lib_records.append(rec)

        input_rec = SeqRecord(Seq(seq),id=seqid)
        fasta_length_dict[seqid] = len(seq)
        my_records.append(input_rec)
           

    SeqIO.write(my_records, input_fasta, "fasta")
    SeqIO.write(lib_records, lib_fasta, "fasta")

    return fasta_length_dict, lib_length_dict

def search_sequence_in_genome_fasta(fasta_file, search_sequence):
    sequence_length = len(search_sequence)
    start_index = None
    end_index = None

    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                if start_index is not None and end_index is not None:
                    return genome_sequence[start_index:end_index]

                genome_sequence = ''
                continue

            genome_sequence += line.strip()

            if start_index is None and search_sequence in genome_sequence:
                index = genome_sequence.index(search_sequence)
                start_index = max(0, index - 20000)
                end_index = min(len(genome_sequence), index + sequence_length + 20000)

    if start_index is not None and end_index is not None:
        return genome_sequence[start_index:end_index]

    return "未找到匹配的序列"

def fasta_to_dataframe(fasta_file):
    sequences = []
    headers = []
    
    # 使用Biopython的SeqIO模块读取FASTA文件
    for record in SeqIO.parse(fasta_file, "fasta"):
        headers.append(record.id)
        sequences.append(str(record.seq))
    
    # 创建DataFrame
    df = pd.DataFrame({"Header": headers, "Sequence": sequences})
    return df

def search_sequence_in_genome_seq(chrom, genome_path, start, end  , start_distance=20000, end_distance=20000):

    record_dict = SeqIO.to_dict(SeqIO.parse(genome_path, "fasta"))

    genome_sequence=str(record_dict[chrom].seq)

    # search_sequence_length = len(search_sequence)
    start_index =  max(0, int(start) - start_distance)
    end_index = min( len(genome_sequence), int(end) + end_distance )

    #基因在更改后基因组的位置 
    if start_index > 0:
        gene_start = start_distance + 1
        gene_end = gene_start + (int(end) - int(start) ) - 1
    elif start_index == 0:
        gene_start = int(start) + 1
        gene_end = gene_start + (int(end) - int(start)) - 1
        
    return gene_start , gene_end, genome_sequence[start_index: end_index]

    # result, num = find_coordinate_by_pattern(search_sequence, genome_sequence)

#Map the target sequence to the reference genome by Blast
def blast_ha(ref_genome, blast_input_file_path, blast_output_file_path):  

    ref_lib=ref_genome.split('/')[-1].split('.')[0]
    seq_length=0
    with open(blast_input_file_path,'r') as ifile:
        for line in ifile:
            if not line[0]=='>':
                seq_length += len(line)-1
                break
    if seq_length > 550:
            evalue='300'
    else:
            evalue=str(int((seq_length*0.5521-7.5856)*0.8))

    #建库比对多次
    os.system("makeblastdb -in "+ref_genome+" -dbtype nucl -parse_seqids -out "+ref_lib)
    os.system("blastn -query "+blast_input_file_path+" -db "+ref_lib+" -outfmt 6 -out "+blast_output_file_path+" -evalue 1e-"+evalue+" -max_target_seqs 5 -num_threads 4")
    os.system("rm %s.n*" % ref_lib)

#Evaluate the feasibility of design with the mapping of the homologous arm to the reference genome
def blast_output_evaluate(genome_path , workdir, blast_input, blast_output):   
    evaluate_output_file_path=workdir+'/Evaluation_result.txt'
    evaluate_input_df = pd.DataFrame()
    middle_df = pd.DataFrame()
    high_df = pd.DataFrame()
    low_df = pd.DataFrame()
    aim_all_df = pd.DataFrame()

    with open(evaluate_output_file_path,'a') as evaluate_output:
        dict_evaluate_output={}
    
        with open(blast_output,'r') as evaluate_input:
            dict_result_id = {}
            evaluate_result = {}
            for line_result in evaluate_input:
                qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, ssend, evalue, bitscore = line_result.split('\t')
                evaluate_result.update({
                    'qseqid':qseqid,
                    'sseqid':sseqid,
                    'pident':float(pident),
                    'length':int(length),
                    'mismatch':int(mismatch),
                    'gapopen':gapopen,
                    'qstart':qstart,
                    'qend':qend,
                    'sstart':sstart,
                    'ssend':ssend,
                    'evalue':evalue,
                    'bitscore':bitscore
                })
                result_id=qseqid
                
                dict_result_id[result_id] = dict_result_id.get(result_id,0) + 1
                evaluate_input_df = evaluate_input_df.append(pd.DataFrame([evaluate_result]))

            
            list_result_id_unmap=[]
            with open(blast_input,'r') as blast_input:
                for lines in blast_input:
                    if lines[0] =='>':
                        blast_input_id=lines[1:-1]
                        blast_input_id = blast_input_id.split(' ')[0]
                        if blast_input_id in dict_result_id:
                            if dict_result_id[blast_input_id] > 1:                               
                                temp = evaluate_input_df[evaluate_input_df['qseqid']==blast_input_id]
                               
                                # aaa = aaa.append(temp) 
                                temp = temp[['qseqid','sseqid','pident','length','mismatch','sstart','ssend','qstart','qend']]

                                if config.UHA_DHA_CONFIG['max_right_arm_seq_length']  == config.UHA_DHA_CONFIG['min_right_arm_seq_length']:
                                    uha_dha_length =  config.UHA_DHA_CONFIG['max_right_arm_seq_length']

                                temp['coverage'] = temp['length'] / uha_dha_length

                                id, position, hr =  blast_input_id.split(';')
                                chrom, start_end = position.split(':')
                                qstart, qend = start_end.split('-')

                                # 明天研究
                                gene_start, gene_end, target_lib_seq = search_sequence_in_genome_seq(chrom, genome_path, start = qstart, end = qend, start_distance=20000, end_distance=20000)
                                print( gene_start, gene_end )

                                if hr == 'UHA':
                                    aim_df =  temp[ ( temp['qstart'] == str(gene_start - uha_dha_length)  ) & ( temp['qend'] == str(gene_start -1) ) ]
                                elif hr == 'DHA':
                                    aim_df =  temp[ ( temp['qstart'] == str(gene_end + 1)   ) & ( temp['qend'] == str(gene_end + uha_dha_length) ) ]

                                aim_df['off target'] = 'target'
                                temp = temp[ ~ ( ( temp['pident'] == 100 ) & ( temp['coverage'] == 1) ) ]
                                
                                high_off_target_df = temp[ ( temp['pident'] > 90 ) & ( temp['coverage'] > 0.9) ]
                                high_off_target_df['off target'] = 'high'
                                
                                surplus_df = temp[~ ( (temp['pident'] > 90) & (temp['coverage'] > 0.9) ) ] 

                                middle_off_target_df =  surplus_df[ ( (surplus_df['pident'] > 90) & (surplus_df['coverage'] > 0.7) ) ]
                                middle_off_target_df['off target'] = 'medium'

                                temp_surplus_df = surplus_df[ ~ ( (surplus_df['pident'] > 90) & (surplus_df['coverage'] > 0.7) )  ]

                                low_off_target_df = temp_surplus_df[ ( (temp_surplus_df['pident'] > 90) & (temp_surplus_df['length'] >100) )  ]
                                low_off_target_df['off target'] = 'low'

                                

                                #检查序列
                                if len(high_off_target_df) >= 1 or len(middle_off_target_df) >= 1 or len(low_off_target_df) >= 1 :
                                    aim_all_df = aim_all_df.append(aim_df)

                                if len(high_off_target_df) >= 1:
                                    evaluate_output.write(blast_input_id+'\t'+'The target sequence can map to multiple positions in the reference genome. The genome editing may be mislocated. High risk of miss target. '+'\n')
                                    high_df = high_df.append( high_off_target_df )
                                    

                                if len(middle_off_target_df) >= 1:
                                    evaluate_output.write(blast_input_id+'\t'+'The target sequence can map to multiple positions in the reference genome. The genome editing may be mislocated. Medium risk of miss target. '+'\n')
                                    middle_df = middle_df.append( middle_off_target_df )

                                if len(low_off_target_df) >= 1:
                                    evaluate_output.write(blast_input_id+'\t'+'The target sequence can map to multiple positions in the reference genome. The genome editing may be mislocated. Low risk of miss target. '+'\n')
                                    low_df = low_df.append( low_off_target_df )

                        else:
                            if blast_input_id in list_result_id_unmap:
                                continue
                            else:
                                evaluate_output.write(blast_input_id+'\t'+'The target sequence can not map to the reference genome. Please check them.'+'\n')
                                list_result_id_unmap.append(blast_input_id)  
    


    return middle_df, high_df, low_df, aim_all_df

from Bio import SeqIO
def gb_2_fna(gb_file,fna_file):

    # 读取 GenBank 文件
    gb_records = SeqIO.parse(gb_file, "genbank")

    # 将 DNA 序列写入 FASTA 文件
    
    with open(fna_file, "w") as output_handle:
        SeqIO.write(gb_records, output_handle, "fasta")

def search_sequence_in_genome_fasta(fasta_file, search_sequence):
    sequence_length = len(search_sequence)
    start_index = None
    end_index = None

    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                if start_index is not None and end_index is not None:
                    return genome_sequence[start_index:end_index]

                genome_sequence = ''
                continue

            genome_sequence += line.strip()

            if start_index is None and search_sequence in genome_sequence:
                index = genome_sequence.index(search_sequence)
                start_index = max(0, index - 20000)
                end_index = min(len(genome_sequence), index + sequence_length + 20000)

    if start_index is not None and end_index is not None:
        return genome_sequence[start_index:end_index]

    return "未找到匹配的序列"

def str2num(x):
    """
    This extracts numbers from strings. eg. 114 from M114R.

    :param x: string
    """
    return int(''.join(ele for ele in x if ele.isdigit()))

def reset_index_df(*li):

    for i in li:
            if len(i)>0:
                #重置索引
                i.index = range(1,len(i)+1)
                

                #重置id
                a = i.ID.str.split(';',expand=True)
                a['ID'] = a[0].apply(lambda x:  x.split('|')[1] if len( x.split('|')  ) > 1 else x.split('|')[0] )

                i.rename(columns={'ID':'Name'}, inplace=True)
                i = i.insert(0,'ID',a['ID'])
            
            
    return li

#
def read_chromSeq_from_genome(path,chrom):

    record_dict = SeqIO.to_dict(SeqIO.parse(path, "fasta"))   
        # chrom = seq_id.split(';')[1].split(':')[0]
    seq = str(record_dict[chrom].seq)
   
    return seq  

    
def off_target_analysis_from_genome(primer_seq='12', seq_id='12' , genome_path='234', primer_type='LEFT'):

    primer_suc = {}

    chrom = seq_id.split(';')[1].split(':')[0]
    genome_seq = read_chromSeq_from_genome(genome_path, chrom)

    result = find_coordinate_by_pattern(pattern=primer_seq, seq=genome_seq)

    #若模板脱靶，直接结束
    if result[1] > 1:
        return primer_suc
    else:
    #若模板不脱靶，循环找出15bp不脱靶的序列
        if primer_type == 'RIGHT':
            primer_seq = revComp(primer_seq)

        #对每一条引物进行脱靶分析，种子序列为15bp
        j = 0
        for m in range(15,len(primer_seq)+1):                     #控制滑动次数

            pattern = primer_seq[j:m]
            result = find_coordinate_by_pattern(pattern=pattern, seq=genome_seq)

            if result[1] > 1:
                if j == len(primer_seq) - 15 + 1 :
                    break
                else:
                    j = j + 1
                    continue
            else:
                coordinate = result[0]['0']

                if primer_type == 'LEFT':
                    primer_suc.update(
                        {
                            'PRIMER_LEFT_0_SEQUENCE':pattern ,
                            'PRIMER_LEFT_0': (coordinate[0] , (coordinate[1] - coordinate[0] + 1) )
                        }
                    )
                elif primer_type == 'RIGHT':
                    primer_suc.update(
                        {
                            'PRIMER_RIGHT_0_SEQUENCE':pattern,  
                            'PRIMER_RIGHT_0': ( len(primer_seq) -1 - coordinate[0] , (coordinate[1] - coordinate[0] + 1) )
                        }
                    )
                break
        
        return primer_suc  
  

def parse_primer3_result(dict_res, primer_num):
        primer_df = pd.DataFrame()
        for u in range(primer_num) :
            dict_primer = {}
            left_primer = dict_res.get(f'PRIMER_LEFT_{u}_SEQUENCE')
            left_primer_tm = dict_res.get(f'PRIMER_LEFT_{u}_TM')
            right_primer = dict_res.get(f'PRIMER_RIGHT_{u}_SEQUENCE')
            right_primer_tm = dict_res.get(f'PRIMER_RIGHT_{u}_TM')

            dict_primer = {
                'id':u,
                'left_primer':left_primer,
                'left_primer_tm':left_primer_tm,
                'right_primer':right_primer,
                'right_primer_tm':right_primer_tm
            }
            primer_df = primer_df.append([dict_primer])
            u = u +1
        return primer_df







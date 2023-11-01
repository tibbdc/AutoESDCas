import pandas as pd
from Bio import SeqIO 
import os,sys
# sys.path.append('..')
import sgRNA_utils.sgRNA_primer_util as su
import module.sequencing_primer as sr
import module.from_gbFile_to_seq as fq
import module.parser_input_to_df as pf
import warnings 
warnings.filterwarnings('ignore')
import configparser
from Bio import SeqIO
from Bio.Seq import Seq   
from Bio.SeqRecord import SeqRecord
import math  
from sgRNA_utils.sgRNA_primer_config import config 


#默认推荐最优sgRNA   
def extract_sgRNA_from_chopchop(sgRNA_result_path, selected_sgRNA_result):

    sgRNA = su.del_Unnamed(pd.read_csv(sgRNA_result_path,index_col=False))
   

    sgRNA.index = list(range(len(sgRNA)))
    sgRNA = sgRNA.reindex(columns=['Name','Region','Rank','Target sequence','Genomic location','Strand','GC content (%)','Self-complementarity','MM0','MM1','MM2','MM3','Efficiency'])
   
    print(selected_sgRNA_result == {})
    if selected_sgRNA_result != {}:
        df = pd.DataFrame()
        for k,v in selected_sgRNA_result.items():
            print(k,v)
            temp_df = sgRNA[sgRNA['Name']==k]
            
            best_sgRNA = temp_df[temp_df['Rank']==int(v)]
            df = df.append(best_sgRNA)
    else:
        def work(x):
            x = x[x['Rank']==1]
            return x
        sgRNA = sgRNA.groupby('Region').apply(lambda x: work(x))
        sgRNA = sgRNA.reset_index(drop=True)
        df = sgRNA
       
    df['Rev Target sequence'] = df['Target sequence'].apply(lambda x: su.revComp(x))  


    return df

def design_primer(primer_template,id_name,template_name,stype):

    result_list = []
    failture_primer_df = pd.DataFrame() 
    for i,row in primer_template.iterrows():
        result_dict = {}
        failture_dict= {}
        if stype == 'u_d':
            if row['type'] == 'uha':
                primer_result = su.primer_design(row[id_name],row[template_name],'left',primer_type='uha')
            elif row['type'] == 'dha':
                primer_result = su.primer_design(row[id_name],row[template_name],'right',primer_type='dha')    
        elif stype == 'plasmid':
            primer_result = su.primer_design(row[id_name],row[template_name],'left_right',primer_type='plasmid')
        elif stype == 'sgRNA':
            primer_result = su.primer_design(row[id_name],row[template_name],'left_right',primer_type='sgRNA')
        elif stype == 'seq_altered':
            primer_result = su.primer_design(row[id_name],row[template_name],'left_right',primer_type='seq_altered')

        print(result_dict)      

        if len(primer_result.keys())<10 :
            #失败
            failture_dict['ID'] = row[id_name]
            if stype == 'u_d':
                failture_dict['type'] = row['type']
            failture_dict.update(primer_result)
            
            failture_primer_df = failture_primer_df.append([failture_dict])
        else:
            #成功
            result_dict['Region'] = row[id_name]
            result_dict["primer_f_seq_(5'-3')"] = primer_result['PRIMER_LEFT_0_SEQUENCE']
            result_dict["primer_f_Tm"] = primer_result['PRIMER_LEFT_0_TM']

            result_dict["primer_r_seq_(5'-3')"] = primer_result['PRIMER_RIGHT_0_SEQUENCE']
            result_dict["primer_r_Tm"] = primer_result['PRIMER_RIGHT_0_TM']

            result_dict['product_size'] = primer_result['PRIMER_PAIR_0_PRODUCT_SIZE']   
            result_dict['product_value'] = row[template_name][primer_result['PRIMER_LEFT_0'][0]:primer_result['PRIMER_RIGHT_0'][0]+1]
            
            result_list.append(result_dict)
                
    primer_df = pd.DataFrame(result_list)
    
    return primer_df, failture_primer_df

def create_primer_template(info_input_df,sgRNA):
    sgRNA_df = pd.merge(info_input_df,sgRNA,left_on=['name','region'],right_on=['Name','Region'],how='inner')
    primer_template = su.columns_2_row_from_one_df(sgRNA_df,in_col=['Name','Region','seq_uha_max_whole','seq_dha_max_whole'],to_col=['Name','Region','primer_template'])
    primer_template['Name_Region'] = primer_template['Name']+';'+primer_template['Region']+';'+primer_template['type']
    primer_template.drop(columns=['Name','Region'],inplace=True)
    return primer_template

def create_uha_dha_df(uha_dha_primer_df):
    in_col = ['Name','Region','product_value','product_size']
    ou_col = ['Name','Region','UHA','UHA_size','DHA','DHA_size']
    UHA_DHA_df = su.columns_2_row_by_groupby(uha_dha_primer_df,in_col,ou_col,type='u_d')
    return UHA_DHA_df

def create_uha_dha_primer_df(uha_dha_primer_df,in_col,ou_col):
    # in_col = ['Name','Region',"primer_f_seq_(5'-3')","primer_r_seq_(5'-3')",'product_value','product_size']
    # ou_col = ['Name','Region',"u_primer_f_seq_(5'-3')","u_rimer_r_seq_(5'-3')",'UHA','UHA_size',"d_primer_f_seq_(5'-3')","d_rimer_r_seq_(5'-3')",'DHA','DHA_size']
    UHA_DHA_df = su.columns_2_row_by_groupby(uha_dha_primer_df,in_col,ou_col,type='primer')
    UHA_DHA_df['Region'] = UHA_DHA_df['Name'] +';'+ UHA_DHA_df['Region']
    
    uha_primer_df = UHA_DHA_df[['Region',
        "u_primer_f_seq_(5'-3')",
                                    "u_primer_r_seq_(5'-3')",
                                    'UHA',
                                    'UHA_size']].rename(columns={"u_primer_f_seq_(5'-3')":"primer_f_seq_(5'-3')",
                                    "u_primer_r_seq_(5'-3')":"primer_r_seq_(5'-3')"})
     

    dha_primer_df = UHA_DHA_df[['Region',"d_primer_f_seq_(5'-3')","d_primer_r_seq_(5'-3')",
        'DHA',
        'DHA_size']].rename(columns={"d_primer_f_seq_(5'-3')":"primer_f_seq_(5'-3')",
        "d_primer_r_seq_(5'-3')":"primer_r_seq_(5'-3')"})     
    
    return uha_primer_df, dha_primer_df

def lambda2cols(df,lambdaf,in_coln,to_colns):         #apply函数的助手
    if len(in_coln) == 3:
        df_=df.apply(lambda x: lambdaf(x[in_coln[0]],x[in_coln[1]],x[in_coln[2]]),
                     axis=1).apply(pd.Series)
    elif len(in_coln) == 4:
        df_=df.apply(lambda x: lambdaf(x[in_coln[0]],x[in_coln[1]],x[in_coln[2]],x[in_coln[3]]),
                     axis=1).apply(pd.Series)
    elif len(in_coln) == 5:
         df_=df.apply(lambda x: lambdaf(x[in_coln[0]],x[in_coln[1]],x[in_coln[2]],x[in_coln[3]],x[in_coln[4]]),
                     axis=1).apply(pd.Series)
    df_.columns=to_colns
    df=df.join(df_)        
    return df

def plasmid_region_division_by_labels(gb_path,ccdb_label='ccdB', promoter_terminator_label='gRNA', n_20_label='N20'):
    gb = SeqIO.read(gb_path, "genbank")
    gb_seq = str(gb.seq)
    #get coordinate
    ccdb_coordinate = su.get_feature_coordinate(ccdb_label,gb_path)
    promoter_terminator_coordinate = su.get_feature_coordinate(promoter_terminator_label,gb_path)  
    #N20
    n20_coordinate = su.get_feature_coordinate(n_20_label,gb_path)
    if ccdb_coordinate[0] !=-1 and promoter_terminator_coordinate[0] != -1:
        pass
    elif ccdb_coordinate[0] == -1 and promoter_terminator_coordinate[0] != -1: 
        #无ccdb
        before_processed_seq_dict, after_processed_seq_dict = fq.get_data_from_genebank(gb_path,marker=ccdb_label,target_gene=promoter_terminator_label)
        promoter_terminator = before_processed_seq_dict['target_gene_seq']
        promoter_terminator_up_seq = before_processed_seq_dict['target_gene_up_seq']
        promoter_terminator_down_seq = before_processed_seq_dict['target_gene_down_seq']
        
        #根据坐标取出n20序列
        n20_coordinate_seq = gb_seq[n20_coordinate[0]:n20_coordinate[1]]
        start = promoter_terminator.find(n20_coordinate_seq)
        terminator_seq = promoter_terminator[start+len(n20_coordinate_seq):]
        promoter_seq = promoter_terminator[:start]
        
        #质粒被分割成几部分
        plasmid_dict_region_seq = {
            'promoter_terminator_up_seq':promoter_terminator_up_seq,
            'promoter_seq':promoter_seq,
            'n20_coordinate_seq':n20_coordinate_seq,
            'terminator_seq':terminator_seq,
            'promoter_terminator_down_seq':promoter_terminator_down_seq
        }
    elif ccdb_coordinate[0] != -1 and promoter_terminator_coordinate[0]  == -1:
        #双质粒系统：无sgRNA
        before_processed_seq_dict, after_processed_seq_dict = fq.get_data_from_genebank(gb_path,marker=promoter_terminator_label,target_gene=ccdb_label)
        ccdb = before_processed_seq_dict['target_gene_seq']
        ccdb_up_seq = before_processed_seq_dict['target_gene_up_seq']
        ccdb_down_seq = before_processed_seq_dict['target_gene_down_seq']
        
        #质粒被分割成几部分
        plasmid_dict_region_seq = {
            'ccdb_up_seq': ccdb_up_seq,
            'ccdb': ccdb,
            'ccdb_down_seq': ccdb_down_seq
        }
    return gb_seq,plasmid_dict_region_seq

def plasmid_region_division_by_coordinate(coordinate_json,region_seq_dict,plasmid_seq):

    if len(coordinate_json) == 0:
        pass
    else:
        for k,v in coordinate_json.items():
            start, end = v.split(',')
            plasmid_seq[start:end]

def design_primer_by_region_in_plasmid(first_primer_start_position, plasmid_seq, distance_dict, temp_variable=0):

    sgRNA_plasmid_seq_len = len(plasmid_seq)
    primer_result_list = []
    failture_result_list = []
    i=1
    relative_distance = ''
    primer_template_start = first_primer_start_position
    

    #针对区域设计引物
    for k,v in distance_dict.items():
        min_distance,max_distance = v[0],v[1]
        region = max_distance - min_distance
        if i == 1:
            primer_template_end = primer_template_start + max_distance
        else:
            primer_template_end = primer_template_start + max_distance + relative_distance
            
        if primer_template_end > sgRNA_plasmid_seq_len and primer_template_end-sgRNA_plasmid_seq_len < first_primer_start_position:
            template_seq_1 = plasmid_seq[primer_template_start:]
            temp_len =  primer_template_end-sgRNA_plasmid_seq_len
            template_seq_2 = plasmid_seq[:temp_len]
            template_seq = template_seq_1 + template_seq_2
            #这就可以设计引物
            primer_result = su.primer_design(k,template_seq, 'right',primer_type='plasmid')
            if primer_result.get('PRIMER_LEFT_0_SEQUENCE') == None or primer_result.get('PRIMER_RIGHT_0_SEQUENCE') == None:
                failture_dict = {}
                failture_dict['ID'] = f"region:{i}"
                failture_dict.update(primer_result)
                failture_result_list.append(failture_dict)
                return 'failture',failture_result_list
            else: 
                dict_primer_result={
                   'Region':i,
                        f"primer_f_seq_(5'-3')" : primer_result['PRIMER_LEFT_0_SEQUENCE'],
                        f"primer_r_seq_(5'-3')" : primer_result['PRIMER_RIGHT_0_SEQUENCE'],
                        f"primer_f_Tm":primer_result['PRIMER_LEFT_0_TM'],
                        f"primer_r_Tm":primer_result['PRIMER_RIGHT_0_TM'],
                        f'product_size' : primer_result['PRIMER_PAIR_0_PRODUCT_SIZE'],
                        f'product_value' : template_seq[:primer_result['PRIMER_PAIR_0_PRODUCT_SIZE']],
                }
                primer_result_list.append(dict_primer_result)
                right_primer_end = (primer_template_start + primer_result['PRIMER_PAIR_0_PRODUCT_SIZE'])-sgRNA_plasmid_seq_len
                primer_template_start = right_primer_end
                relative_distance = len(template_seq) - primer_result['PRIMER_PAIR_0_PRODUCT_SIZE']
            i = i + 1
        else:
            template_seq = plasmid_seq[primer_template_start : primer_template_end]
            #设计引物
            primer_result = su.primer_design(k,template_seq, 'right',primer_type='plasmid')
            if primer_result.get('PRIMER_LEFT_0_SEQUENCE') == None or primer_result.get('PRIMER_RIGHT_0_SEQUENCE') == None:
                failture_dict = {}
                failture_dict['ID'] = f"region:{i}"
                failture_dict.update(primer_result)
                failture_result_list.append(failture_dict)
                return 'failture',failture_result_list
            else:
                dict_primer_result={
                     'Region':i,
                        f"primer_f_seq_(5'-3')" : primer_result['PRIMER_LEFT_0_SEQUENCE'],
                        f"primer_r_seq_(5'-3')" : primer_result['PRIMER_RIGHT_0_SEQUENCE'],
                        f"primer_f_Tm":primer_result['PRIMER_LEFT_0_TM'],
                        f"primer_r_Tm":primer_result['PRIMER_RIGHT_0_TM'],
                        f'product_size' : primer_result['PRIMER_PAIR_0_PRODUCT_SIZE'],
                        f'product_value' : template_seq[:primer_result['PRIMER_PAIR_0_PRODUCT_SIZE']],
                }
                primer_result_list.append(dict_primer_result)
                right_primer_end = primer_template_start + primer_result['PRIMER_PAIR_0_PRODUCT_SIZE']
                primer_template_start = right_primer_end
                relative_distance = len(template_seq) - primer_result['PRIMER_PAIR_0_PRODUCT_SIZE']
            i = i + 1

    #处理设计最后一对引物


    #获取已经取到的最后一对引物
    last_primer_end = plasmid_seq.find(su.revComp(primer_result_list[-1]["primer_r_seq_(5'-3')"])) + len( su.revComp(primer_result_list[-1]["primer_r_seq_(5'-3')"]) )
   
    print(last_primer_end , first_primer_start_position, "hdsfjkgfhjgjghfj")
    if last_primer_end > first_primer_start_position:

        primer_template = plasmid_seq[last_primer_end:] + plasmid_seq[:first_primer_start_position-temp_variable]
        primer_result = su.primer_design('last',primer_template, 'left_right', primer_type='plasmid')
        
    else:
        # temp_len = last_product_value_seq_end - sgRNA_plasmid_seq_len

        primer_template = plasmid_seq[last_primer_end:first_primer_start_position-temp_variable]
        primer_result = su.primer_design('last',primer_template, 'left_right', primer_type='plasmid')
    
    last_num = i
    if primer_result.get('PRIMER_LEFT_0_SEQUENCE') == None or primer_result.get('PRIMER_RIGHT_0_SEQUENCE') == None:
        failture_dict = {}
        failture_dict['ID'] = f"region:{last_num}"
        failture_dict.update(primer_result)
        failture_result_list.append(failture_dict)
        return 'failture',failture_result_list
    else:
        dict_primer_result={
                            'Region':last_num,
                            f"primer_f_seq_(5'-3')" : primer_result['PRIMER_LEFT_0_SEQUENCE'],
                            f"primer_r_seq_(5'-3')" : primer_result['PRIMER_RIGHT_0_SEQUENCE'],
                            f"primer_f_Tm":primer_result['PRIMER_LEFT_0_TM'],
                            f"primer_r_Tm":primer_result['PRIMER_RIGHT_0_TM'],
                            f'product_size' : primer_result['PRIMER_PAIR_0_PRODUCT_SIZE'],
                            f'product_value' : template_seq[:primer_result['PRIMER_PAIR_0_PRODUCT_SIZE']],
                        }
        primer_result_list.append(dict_primer_result)

    return 'success',primer_result_list


def get_plasmid_backbone_by_labels(gb_path, ccdb_label='ccdB', promoter_terminator_label='gRNA', n_20_label='N20'):

    gb = SeqIO.read(gb_path, "genbank")
    #get coordinate
    ccdb_coordinate = su.get_feature_coordinate(ccdb_label,gb_path)
    promoter_terminator_coordinate = su.get_feature_coordinate(promoter_terminator_label,gb_path)  
  

    if ccdb_coordinate[0] !=-1 and promoter_terminator_coordinate[0] != -1:
        before_processed_seq_dict, after_processed_seq_dict = fq.get_data_from_genebank(gb_path,marker=ccdb_label, target_gene=n_20_label)
       
        # ccdb = before_processed_seq_dict['marker_seq']
        n_20 = after_processed_seq_dict['target_gene_seq']
        n20_up_seq = after_processed_seq_dict['target_gene_up_seq']    
        n20_down_seq = after_processed_seq_dict['target_gene_down_seq']
             
        plasmid_backbone = n20_down_seq
    
    elif ccdb_coordinate[0] == -1 and promoter_terminator_coordinate[0] != -1:

         #双质粒系统:无ccdb
        before_processed_seq_dict, after_processed_seq_dict = fq.get_data_from_genebank(gb_path,marker=ccdb_label, target_gene=n_20_label)

        # ccdb = before_processed_seq_dict['marker_seq']
        n_20 = before_processed_seq_dict['target_gene_seq']
        n20_up_seq = before_processed_seq_dict['target_gene_up_seq']
        n20_down_seq = before_processed_seq_dict['target_gene_down_seq']
        
        plasmid_backbone = n20_up_seq + n_20 + n20_down_seq


    elif ccdb_coordinate[0] != -1 and promoter_terminator_coordinate[0]  == -1:

         #双质粒系统：无sgRNA
        before_processed_seq_dict, after_processed_seq_dict = fq.get_data_from_genebank(gb_path,marker=n_20_label,target_gene=ccdb_label)
        ccdb = before_processed_seq_dict['target_gene_seq']
        ccdb_up_seq = before_processed_seq_dict['target_gene_up_seq']
        ccdb_down_seq = before_processed_seq_dict['target_gene_down_seq']
       
        plasmid_backbone = ccdb_down_seq + ccdb_up_seq

    return  plasmid_backbone


def create_new_plasmid(gb_path, sgRNA_df):
    # ccdb_label='ccdB', promoter_terminator_label='gRNA', n_20_label='N20', promoter_label='promoter'
    
    ccdb_label = config.PLASMID_LABEL['ccdb_label']
    promoter_terminator_label = config.PLASMID_LABEL['promoter_terminator_label']
    n_20_label = config.PLASMID_LABEL['n_20_label']
    promoter_label = config.PLASMID_LABEL['promoter_label']


    gb = SeqIO.read(gb_path, "genbank")
    gb_seq = str(gb.seq)
    #get coordinate
    ccdb_coordinate = su.get_feature_coordinate(ccdb_label,gb_path)
    promoter_terminator_coordinate = su.get_feature_coordinate(promoter_terminator_label,gb_path)  
    #N20
    n20_coordinate = su.get_feature_coordinate(n_20_label,gb_path)


     
    #判断单、双质粒系统
    #type_kind=1:sgRNA启动子序列上游长度>600bp，type_kind=2:sgRNA启动子序列上游长度<600bp，type_kind=3:双质粒系统无ccdb，type_kind=4:双质粒系统无sgRNA
    if ccdb_coordinate[0] !=-1 and promoter_terminator_coordinate[0] != -1:
        #单质粒系统
            #对质粒划分
        before_processed_seq_dict, after_processed_seq_dict = fq.get_data_from_genebank(gb_path,marker=ccdb_label,target_gene=promoter_terminator_label)

        ccdb = after_processed_seq_dict['marker_seq']
        promoter_terminator = after_processed_seq_dict['target_gene_seq']
        promoter_terminator_up_seq = after_processed_seq_dict['target_gene_up_seq']
        promoter_terminator_down_seq = after_processed_seq_dict['target_gene_down_seq']

        print("启动子-n20-终止子上游序列：",promoter_terminator_up_seq)

        print("启动子—n20-终止子下游序列：",promoter_terminator_down_seq)

        n20_coordinate_seq = gb_seq[n20_coordinate[0]:n20_coordinate[1]]
        start = promoter_terminator.find(n20_coordinate_seq)
        terminator_seq = promoter_terminator[start+len(n20_coordinate_seq):]
        promoter_seq = promoter_terminator[:start]
        
        joint_need_seq = ''


        # if len(promoter_terminator_up_seq) >= 600: 
        # type_kind = 1                       
        plasmid_backbone = terminator_seq + promoter_terminator_down_seq

        real_promoter = su.get_sequence_by_feature_label(gb_path, feature_label = promoter_label)
        if real_promoter == None:
            raise ValueError('There is a problem with the label you provided')

        def work(strand, target_seq, uha, dha, seq_altered):
                
                if str( real_promoter.upper() ) != promoter_terminator[:start].upper():
                    target_seq = su.revComp(target_seq)  

                new_promoter_terminator = promoter_terminator[:start] + target_seq + terminator_seq
                promoter_up_promoter = promoter_terminator_up_seq +  promoter_seq  
                terminator_terminator_down = terminator_seq + promoter_terminator_down_seq
                if seq_altered == '-':
                    plasmid = uha + dha + promoter_terminator_up_seq +  new_promoter_terminator + promoter_terminator_down_seq
                else:
                    plasmid = uha + seq_altered + dha + promoter_terminator_up_seq +  new_promoter_terminator + promoter_terminator_down_seq

                return plasmid, promoter_up_promoter, terminator_terminator_down, new_promoter_terminator, promoter_terminator_up_seq, promoter_terminator_down_seq       

        sgRNA_df = su.lambda2cols(df=sgRNA_df, lambdaf=work, in_coln=["Strand",'Target sequence','UHA','DHA','seq_altered'], to_colns=['plasmid','n20_up_template','n20_down_template','promoter_N20_terminator','promoter_N20_terminator_up','promoter_N20_terminator_down'])
        
        promoter_terminator_up_promoter_seq = promoter_terminator_up_seq + promoter_seq
        promoter_terminator_down_terminator_seq = terminator_seq + promoter_terminator_down_seq

        type_kind = -1
        if len(promoter_terminator_up_seq) > 600 and len(promoter_terminator_down_seq) > 600:
            type_kind = 1
        elif  len(promoter_terminator_up_seq) <= 600 :
            type_kind = 2
        elif len(promoter_terminator_down_seq) <= 600 :
            type_kind = 3

        return sgRNA_df, promoter_seq, promoter_terminator_up_promoter_seq, promoter_terminator_down_terminator_seq, type_kind

    elif ccdb_coordinate[0] == -1 and promoter_terminator_coordinate[0] != -1:
        #双质粒系统:无ccdb
        before_processed_seq_dict, after_processed_seq_dict = fq.get_data_from_genebank(gb_path,marker=ccdb_label,target_gene=promoter_terminator_label)
        promoter_terminator = before_processed_seq_dict['target_gene_seq']
        promoter_terminator_up_seq = before_processed_seq_dict['target_gene_up_seq']
        promoter_terminator_down_seq = before_processed_seq_dict['target_gene_down_seq']
        
        n20_coordinate_seq = gb_seq[n20_coordinate[0]:n20_coordinate[1]]
        start = promoter_terminator.find(n20_coordinate_seq)
        terminator_seq = promoter_terminator[start+len(n20_coordinate_seq):]
        promoter_seq = promoter_terminator[:start]

        
        
        joint_need_seq = promoter_terminator
        plasmid_backbone = promoter_terminator_up_seq + promoter_terminator + promoter_terminator_down_seq
        
        real_promoter = su.get_sequence_by_feature_label(gb_path, feature_label = promoter_label)
        if real_promoter == None:
            raise ValueError('There is a problem with the label you provided')
        def work(target_seq, uha,dha, seq_altered, type):

            if real_promoter.upper() != promoter_terminator[:start].upper():  
                    target_seq = su.revComp(target_seq)

            new_promoter_terminator = promoter_terminator[:start] + target_seq + terminator_seq
            plasmid = promoter_terminator_up_seq + new_promoter_terminator + promoter_terminator_down_seq
            sgRNA_template = terminator_seq + promoter_seq
           
            return plasmid, sgRNA_template, new_promoter_terminator, promoter_terminator_up_seq, promoter_terminator_down_seq

        sgRNA_df = su.lambda2cols(df=sgRNA_df, lambdaf=work, in_coln=['Target sequence','UHA','DHA','seq_altered','type'], to_colns=['plasmid','sgRNA_template','promoter_N20_terminator','promoter_N20_terminator_up','promoter_N20_terminator_down'])        
        return sgRNA_df, promoter_seq, plasmid_backbone, promoter_seq, terminator_seq, joint_need_seq

    elif ccdb_coordinate[0] != -1 and promoter_terminator_coordinate[0]  == -1:
        #双质粒系统：无sgRNA
        before_processed_seq_dict, after_processed_seq_dict = fq.get_data_from_genebank(gb_path,marker=promoter_terminator_label,target_gene=ccdb_label)
        ccdb = before_processed_seq_dict['target_gene_seq']
        ccdb_up_seq = before_processed_seq_dict['target_gene_up_seq']
        ccdb_down_seq = before_processed_seq_dict['target_gene_down_seq']

       
        plasmid_backbone = ccdb_down_seq + ccdb_up_seq
        joint_need_seq = ''
        def work(target_seq, uha,dha, seq_altered, type):
            
            if type == 'insertion' or type == 'substitution':
                if seq_altered == '-':
                    plasmid = uha + dha + ccdb_down_seq + ccdb_up_seq  
                else:
                    plasmid = uha + seq_altered + dha + ccdb_down_seq + ccdb_up_seq
            elif type == 'deletion':
                plasmid = uha + dha + ccdb_down_seq + ccdb_up_seq
            
            return plasmid, '', ccdb, ccdb_up_seq, ccdb_down_seq

        sgRNA_df = su.lambda2cols(df=sgRNA_df, lambdaf=work, in_coln=['Target sequence','UHA','DHA','seq_altered','type'], to_colns=['plasmid','sgRNA_template','ccdb','ccdb_up','ccdb_down'])
        return sgRNA_df, plasmid_backbone, joint_need_seq


def add_joint_sgRNA_primer(sgRNA_primer_df,enzyme_df,enzyme_name,promoter_terminator_down_terminator_seq='',promoter_terminator_up_promoter_seq='',stype='sgRNA_joint'):
    
    sgRNA_enzyme_df = enzyme_df[enzyme_df['name']==enzyme_name]
    sgRNA_enzyme_df = sgRNA_enzyme_df.reset_index(drop=True)
    protective_base = sgRNA_enzyme_df.loc[0,'protective_base']
    recognition_seq = sgRNA_enzyme_df.loc[0,'recognition_seq']
    cut_seq_len = sgRNA_enzyme_df.loc[0,'cut_seq_len']
    gap_len = sgRNA_enzyme_df.loc[0,'gap_len']
    gap_seq = 'AGACTAGACTAGACTAGACTAGACTAGACTAGACTAGACTAGACT'
    gap_seq = sgRNA_enzyme_df.loc[0,'gap_seq']   


    def work(*x):       
        if len(x) == 4:  
            left_primer,right_primer,product_value,Name = x[0],x[1],x[2],x[3]
        if len(x) == 5:
            left_primer,right_primer,product_value,rev_target_seq,Name = x[0],x[1],x[2],x[3],x[4]   

        if stype == 'n20up_primer_joint'  or stype == 'sgRNA_plasmid_primer_joint' :
          
                left_temp_seq = protective_base + recognition_seq + gap_seq[:gap_len]
                right_temp_seq = protective_base + recognition_seq + gap_seq[:gap_len]

                left_primer = left_temp_seq + left_primer  

                right_primer = su.revComp(right_temp_seq)[::-1] + rev_target_seq + right_primer

                product_value_joint = left_temp_seq + product_value + su.revComp(rev_target_seq)[::-1] + right_temp_seq[::-1]
                product_value_size_joint = len(product_value_joint)

                print(stype,f'n20up_left_{Name}:', protective_base , recognition_seq , gap_seq[:gap_len], left_primer)
                print(stype,f'n20up_right_{Name}:', su.revComp(right_temp_seq)[::-1] , rev_target_seq , right_primer ,';',  protective_base , recognition_seq , gap_seq[:gap_len])

        elif stype == 'n20down_primer_joint':
            left_temp_seq = protective_base + recognition_seq + gap_seq[:gap_len]  + su.revComp(rev_target_seq)[-cut_seq_len:]
            right_temp_seq = protective_base + recognition_seq + gap_seq[:gap_len]  

            left_primer = left_temp_seq + left_primer   
            right_primer = su.revComp(right_temp_seq)[::-1]+ right_primer

            product_value_joint = left_temp_seq + product_value + right_temp_seq[::-1]   
            product_value_size_joint = len(product_value_joint)
            
            print(stype,f'n20down_left_{Name}:', protective_base , recognition_seq , gap_seq[:gap_len], su.revComp(rev_target_seq)[-cut_seq_len:], left_primer)
            print(stype,f'n20down_right_{Name}:', su.revComp(right_temp_seq)[::-1] , right_primer ,';',  protective_base , recognition_seq , gap_seq[:gap_len])

        elif stype == 'plasmid_backbone_primer_joint':
        # elif stype == 'n20down_primer_joint':
            left_temp_seq = protective_base + recognition_seq + gap_seq[:gap_len]
            right_temp_seq = protective_base + recognition_seq + gap_seq[:gap_len]  

            left_primer = left_temp_seq + left_primer   
            right_primer = su.revComp(right_temp_seq)[::-1]+ right_primer

            product_value_joint = left_temp_seq + product_value + right_temp_seq[::-1]   
            product_value_size_joint = len(product_value_joint)

            print(f'plasmid_left_{Name}:', protective_base , recognition_seq , gap_seq[:gap_len], left_primer)  
            print(f'plasmid_right_{Name}:', su.revComp(right_temp_seq)[::-1] , right_primer ,';', protective_base , recognition_seq , gap_seq[:gap_len])

        elif stype == 'seq_altered_primer_joint'or stype == 'ccdb_plasmid_primer_joint':
            left_temp_seq = protective_base + recognition_seq + gap_seq[:gap_len]  
            right_temp_seq = protective_base + recognition_seq + gap_seq[:gap_len]   

            left_primer = left_temp_seq + left_primer
            right_primer = su.revComp(right_temp_seq)[::-1] + right_primer
            
            product_value_joint = left_temp_seq + product_value + right_temp_seq[::-1]
            product_value_size_joint = len(product_value_joint)

            print(f'seq_altered_left_{Name}:', protective_base , recognition_seq , gap_seq[:gap_len], left_primer)  
            print(f'seq_altered_right_{Name}:', su.revComp(right_temp_seq)[::-1] , right_primer ,';', protective_base , recognition_seq , gap_seq[:gap_len])  

        return left_primer, right_primer, product_value_joint, product_value_size_joint

    #设计uha，dha的上下游引物接头 
    def u_d_primer_joint(x): 
        x = x.reset_index(drop=True)  
        mute_type = x.loc[0,'type']
        seq_altered	 = x.loc[0,'seq_altered']
        seq_altered_len = len(seq_altered)  
        add_len = math.floor( seq_altered_len/2 ) 

        df = pd.DataFrame(columns=['Name',
                                'Region',
                                "primer_f_seq_(5'-3')_joint",
                                "primer_r_seq_(5'-3')_joint",
                                "product_value_joint",
                                "product_size_joint",
                                'Type'])
        u_d_primer_joint_dict={}
        uha_df = x[x['Type']=='uha'].reset_index(drop=True)
        Name, Region, uha_left_primer, uha_right_primer, uha_product_value, uha_type = uha_df.loc[0,"Name"],\
                                                                        uha_df.loc[0,"Region"],\
                                                                        uha_df.loc[0,"primer_f_seq_(5'-3')"],\
                                                                        uha_df.loc[0,"primer_r_seq_(5'-3')"],\
                                                                        uha_df.loc[0,"product_value"],\
                                                                        uha_df.loc[0,"Type"]
        dha_df = x[x['Type']=='dha'].reset_index(drop=True)
        dha_left_primer, dha_right_primer, dha_product_value, dha_type = dha_df.loc[0,"primer_f_seq_(5'-3')"],\
                                                                        dha_df.loc[0,"primer_r_seq_(5'-3')"],\
                                                                        dha_df.loc[0,"product_value"],\
                                                                        dha_df.loc[0,"Type"]
    
      

       
        if uha_type == 'uha':
            #uha的左引物
            uha_left_temp_seq = protective_base + recognition_seq + gap_seq[:gap_len] + promoter_terminator_down_terminator_seq[-cut_seq_len:]
            uha_left_primer = uha_left_temp_seq + uha_left_primer  

            #uha的右引物
            if mute_type == 'deletion':
                uha_right_temp_seq = protective_base + recognition_seq + gap_seq[:gap_len]
                uha_right_primer =  su.revComp(uha_right_temp_seq)[::-1]  + uha_right_primer

            elif mute_type == 'substitution' or mute_type == 'insertion':
                if seq_altered_len <= cut_seq_len :
                        uha_right_temp_seq = protective_base + recognition_seq + gap_seq[:gap_len] + seq_altered
                        uha_right_primer = su.revComp(uha_right_temp_seq[:-seq_altered_len])[::-1] + su.revComp(seq_altered)  + uha_right_primer
                        print(stype,f'uha_right_{Name}{Region}:',su.revComp(uha_right_temp_seq[:-seq_altered_len])[::-1], su.revComp(seq_altered), uha_right_primer, ';' , protective_base, recognition_seq, gap_seq[:gap_len])

                elif seq_altered_len > cut_seq_len and  seq_altered_len <= 120:
                        uha_right_temp_seq = protective_base + recognition_seq + gap_seq[:gap_len] + seq_altered[:add_len]
                        uha_right_primer = su.revComp(uha_right_temp_seq[:-add_len])[::-1] + su.revComp(seq_altered[:add_len]) + uha_right_primer
                        print(stype,f'uha_right_{Name}{Region}:', su.revComp(uha_right_temp_seq[:-add_len])[::-1], su.revComp(seq_altered[:add_len]), uha_right_primer, ';' , protective_base, recognition_seq, gap_seq[:gap_len])

                elif seq_altered_len > 120:
                    uha_right_temp_seq = protective_base + recognition_seq + gap_seq[:gap_len] + seq_altered[:cut_seq_len]
                    uha_right_primer = su.revComp(uha_right_temp_seq[:-cut_seq_len])[::-1] + su.revComp(seq_altered[:cut_seq_len]) + uha_right_primer
                    print(stype,f'uha_right_{Name}{Region}:',  su.revComp(uha_right_temp_seq[:-cut_seq_len])[::-1], su.revComp(seq_altered[:cut_seq_len]), uha_right_primer, ';' , protective_base, recognition_seq, gap_seq[:gap_len])

            uha_product_value_joint = uha_left_temp_seq + uha_product_value + uha_right_temp_seq[::-1] 
            uha_product_value_size_joint = len(uha_product_value_joint)

            u_d_primer_joint_dict['Name']=Name
            u_d_primer_joint_dict['Region']=Region
            u_d_primer_joint_dict["primer_f_seq_(5'-3')_joint"] = uha_left_primer   
            u_d_primer_joint_dict["primer_r_seq_(5'-3')_joint"] = uha_right_primer
            u_d_primer_joint_dict["product_value_joint"] = uha_product_value_joint
            u_d_primer_joint_dict["product_size_joint"] = uha_product_value_size_joint
            u_d_primer_joint_dict['Type'] = uha_type
            print(stype,f'uha_left_{Name}{Region}:',protective_base, recognition_seq, gap_seq[:gap_len], promoter_terminator_down_terminator_seq[-cut_seq_len:],uha_left_primer)
           
            df = df.append(pd.DataFrame([u_d_primer_joint_dict]))        

        if dha_type == 'dha':
           
            #dha右引物
            dha_right_temp_seq = protective_base + recognition_seq + gap_seq[:gap_len]
            dha_right_primer = su.revComp(dha_right_temp_seq)[::-1] + su.revComp(promoter_terminator_up_promoter_seq[:cut_seq_len]) + dha_right_primer

            #dha左引物
            uha_right_temp_seq_len=len(uha_right_temp_seq)
            
            if mute_type == 'deletion':  
                dha_left_temp_seq = protective_base + recognition_seq + gap_seq[:gap_len] + su.revComp(uha_right_primer[uha_right_temp_seq_len : uha_right_temp_seq_len + cut_seq_len])[::-1]
            elif mute_type == 'substitution' or mute_type == 'insertion':
               
                if seq_altered_len <= cut_seq_len :
                    temp_len = cut_seq_len - seq_altered_len
                    temp_seq=su.revComp(uha_right_primer[uha_right_temp_seq_len : uha_right_temp_seq_len + temp_len])[::-1]
                    dha_left_temp_seq = protective_base + recognition_seq + gap_seq[:gap_len]  + temp_seq + seq_altered           #对    
                    print(stype,f'dha_left_{Name}{Region}:',protective_base, recognition_seq, gap_seq[:gap_len],temp_seq, seq_altered,dha_left_primer)

                elif seq_altered_len > cut_seq_len and  seq_altered_len <= 120: 
                    
                    if add_len >= cut_seq_len:
                        if seq_altered_len % 2 == 0:
                            dha_left_temp_seq = protective_base + recognition_seq + gap_seq[:gap_len] + seq_altered[add_len - cut_seq_len:]
                            print(stype,f'dha_left_{Name}{Region}:',protective_base, recognition_seq, gap_seq[:gap_len], seq_altered[add_len - cut_seq_len:],dha_left_primer) 
                        else:
                            dha_left_temp_seq = protective_base + recognition_seq + gap_seq[:gap_len] + seq_altered[add_len - cut_seq_len:]
                            print(stype,f'dha_left_{Name}{Region}:',protective_base, recognition_seq, gap_seq[:gap_len], seq_altered[add_len - cut_seq_len:],dha_left_primer) 
                    else:
                        temp_len = cut_seq_len - add_len
                        if seq_altered_len % 2 == 0:
                            #保护碱基 + 识别序列 +（向uha右引物借的序列+重叠序列）+ seq_altered[add_len:]
                            dha_left_temp_seq = protective_base + recognition_seq + gap_seq[:gap_len]+ (su.revComp(uha_right_primer[uha_right_temp_seq_len:temp_len][::-1]) ) + seq_altered
                            print(stype,f'dha_left_{Name}{Region}:',protective_base, recognition_seq, gap_seq[:gap_len], (su.revComp(uha_right_primer[uha_right_temp_seq_len:temp_len][::-1]) ), seq_altered) 
                        else:
                            dha_left_temp_seq = protective_base + recognition_seq + gap_seq[:gap_len]+ (su.revComp(uha_right_primer[uha_right_temp_seq_len:temp_len][::-1]) ) + seq_altered
                            print(stype,f'dha_left_{Name}{Region}:',protective_base, recognition_seq, gap_seq[:gap_len], (su.revComp(uha_right_primer[uha_right_temp_seq_len:temp_len][::-1]) ), seq_altered)
                elif seq_altered_len > 120:
                    dha_left_temp_seq = protective_base + recognition_seq + gap_seq[:gap_len] + seq_altered[-cut_seq_len:]    

            dha_left_primer = dha_left_temp_seq + dha_left_primer   
            #产物   
            dha_product_value_joint = dha_left_temp_seq + dha_product_value + dha_right_temp_seq[::-1]   
            dha_product_value_size_joint = len(dha_product_value_joint)  
            
            #赋值  
            u_d_primer_joint_dict["primer_f_seq_(5'-3')_joint"] = dha_left_primer
            u_d_primer_joint_dict["primer_r_seq_(5'-3')_joint"] = dha_right_primer
            u_d_primer_joint_dict["product_value_joint"] = dha_product_value_joint
            u_d_primer_joint_dict["product_size_joint"] = dha_product_value_size_joint 
            u_d_primer_joint_dict['Type'] = dha_type 
            print(stype,f'dha_right_{Name}{Region}:',su.revComp(dha_right_temp_seq)[::-1], su.revComp(promoter_terminator_up_promoter_seq[:cut_seq_len]), dha_right_primer, ';' ,protective_base ,recognition_seq ,gap_seq[:gap_len])
            df = df.append(pd.DataFrame([u_d_primer_joint_dict]))  
        return df
    
    if stype == 'n20up_primer_joint' or stype == 'sgRNA_plasmid_primer_joint' or stype =='n20down_primer_joint':
        sgRNA_primer_df = su.lambda2cols(sgRNA_primer_df,work,["primer_f_seq_(5'-3')","primer_r_seq_(5'-3')","product_value","Rev Target sequence","Region"],["primer_f_seq_(5'-3')_joint","primer_r_seq_(5'-3')_joint","product_value_joint","product_size_joint"])

        sgRNA_primer_df = sgRNA_primer_df.reindex(columns=[
                                                'Region',
                                                "primer_f_seq_(5'-3')",
                                                "primer_r_seq_(5'-3')",
                                                "primer_f_Tm",  
                                                "primer_r_Tm",
                                                "primer_f_seq_(5'-3')_joint",
                                                "primer_r_seq_(5'-3')_joint",
                                                'product_value',
                                                'product_size',
                                                "product_value_joint",
                                                "product_size_joint"
                                            ])

    elif stype == 'plasmid_backbone_primer_joint':
        sgRNA_primer_df = su.lambda2cols(sgRNA_primer_df,work,["primer_f_seq_(5'-3')","primer_r_seq_(5'-3')","product_value","Region"],["primer_f_seq_(5'-3')_joint","primer_r_seq_(5'-3')_joint","product_value_joint","product_size_joint"])

        sgRNA_primer_df = sgRNA_primer_df.reindex(columns=[
                                                'Region',
                                                "primer_f_seq_(5'-3')",
                                                "primer_r_seq_(5'-3')",
                                                "primer_f_Tm",  
                                                "primer_r_Tm",
                                                "primer_f_seq_(5'-3')_joint",
                                                "primer_r_seq_(5'-3')_joint",
                                                'product_value',
                                                'product_size',
                                                "product_value_joint",
                                                "product_size_joint"
                                            ])

    elif stype == 'u_d_primer_joint':
        # sgRNA_primer_df = sgRNA_primer_df.groupby(by=['Name','Region'])
        # sgRNA_primer_df = su.lambda2cols(sgRNA_primer_df,work,["primer_f_seq_(5'-3')","primer_r_seq_(5'-3')","product_value","Type","product_value"],["primer_f_seq_(5'-3')_joint","primer_r_seq_(5'-3')_joint","product_value_joint","product_size_joint"])
        
        aa = sgRNA_primer_df.groupby(by=['Name','Region']).apply(lambda x:u_d_primer_joint(x))
        aa.index = sgRNA_primer_df.index
        sgRNA_primer_df = pd.merge(aa,sgRNA_primer_df)
        
        sgRNA_primer_df = sgRNA_primer_df.reindex(columns=[
                                                'Name',
                                                'Region',
                                                "primer_f_seq_(5'-3')",
                                                "primer_r_seq_(5'-3')",
                                                "primer_f_Tm",  
                                                "primer_r_Tm",
                                                "primer_f_seq_(5'-3')_joint",
                                                "primer_r_seq_(5'-3')_joint",
                                                'product_value',
                                                'product_size',
                                                "product_value_joint",
                                                "product_size_joint",  
                                                "Type"                 
                                            ])

    elif stype == 'seq_altered_primer_joint' or stype == 'ccdb_plasmid_primer_joint':
        sgRNA_primer_df = su.lambda2cols(sgRNA_primer_df, work, ["primer_f_seq_(5'-3')","primer_r_seq_(5'-3')","product_value","Region"],["primer_f_seq_(5'-3')_joint","primer_r_seq_(5'-3')_joint","product_value_joint","product_size_joint"])
        seq_altered_primer_df = sgRNA_primer_df.reindex(columns=[
                                                'Region',
                                                "primer_f_seq_(5'-3')",
                                                "primer_r_seq_(5'-3')",
                                                "primer_f_Tm",  
                                                "primer_r_Tm",
                                                "primer_f_seq_(5'-3')_joint",
                                                "primer_r_seq_(5'-3')_joint",
                                                'product_value',
                                                'product_size',
                                                "product_value_joint",
                                                "product_size_joint"
                                            ])
     
    return sgRNA_primer_df       

#构建重组载体的测序引物
def create_sequencing_primer(gb_path,sgRNA_df,sequencing_primer,sequencing_template='plasmid',sequencing_region='plasmid_sequencing_region',seq_type=''):
    
    sgRNA_plasmid_df = pd.DataFrame()
    failture_sgRNA_plasmid_df = pd.DataFrame()
    off_target_df = pd.DataFrame()


    for i,v in sgRNA_df.iterrows():
        sgRNA_plasmid_seq=v[sequencing_template]
        # sgRNA_plasmid_seq_len = len(sgRNA_plasmid_seq)
        plasmid_sequencing_region = v[sequencing_region]
        
        #定义测序区域
        start = sgRNA_plasmid_seq.find(plasmid_sequencing_region) 
        end = start + len(plasmid_sequencing_region)
        print(start, end)
        
        target_gene_seq = plasmid_sequencing_region
        target_gene_up_seq = sgRNA_plasmid_seq[0:start]
        target_gene_down_seq =sgRNA_plasmid_seq[end:len(sgRNA_plasmid_seq)]  

        #定义用于测序的数据结构  
        dict_plasmid_seq={}
        # if seq_type == 'genome_seq':
            
        dict_plasmid_seq['target_gene_down_seq']=target_gene_down_seq
        dict_plasmid_seq['target_gene_up_seq']=target_gene_up_seq
        dict_plasmid_seq['mute_after_target_gene_seq']=target_gene_seq     #用户指定测序序列的区域序列
   
        #测序质粒的id
        if seq_type == 'genome_seq':
            gene_change_type = v['type']
            ref = v['ref']
            seq_altered = v['seq_altered']
            sucesse,failture = sequencing_primer.design_sequencing_primers(sgRNA_plasmid_seq,
                                                                            gb_path,
                                                                            v['Region'],
                                                                            dict_plasmid_seq, 
                                                                            seq_type=seq_type, 
                                                                            gene_change_type=gene_change_type, 
                                                                            ref_seq=ref, 
                                                                            seq_altered=seq_altered )
        else:
            sucesse,failture = sequencing_primer.design_sequencing_primers(sgRNA_plasmid_seq,
                                                                            gb_path,
                                                                            v['Region'],
                                                                            dict_plasmid_seq, 
                                                                            seq_type=seq_type)

        if sucesse != {}:
            result = {'Region':v['Region']}
            result.update(sucesse)
            sgRNA_plasmid_df = sgRNA_plasmid_df.append(pd.DataFrame([result]))

        if failture != {}:
          
            for k,v in failture.items():

                id = k

                if id.split(':')[-1] =='off_target':
                    #处理脱靶引物
                    temp = {'ID':id}
                    value={'off_target':v}
                    temp.update(value)   
                    off_target_df = off_target_df.append(pd.DataFrame([temp]))
                else:
                    failture_result={'ID':id}
                    value = v
                    failture_result.update(value)  

                    # failture_result.update(failture)
                    failture_sgRNA_plasmid_df = failture_sgRNA_plasmid_df.append(pd.DataFrame([failture_result]))

        
    if len(sgRNA_plasmid_df) > 0: 
        sgRNA_plasmid_df = su.make_id_in_first_columns(sgRNA_plasmid_df,id='Region',columns=sgRNA_plasmid_df.columns.tolist())


    return sgRNA_plasmid_df, failture_sgRNA_plasmid_df, off_target_df

#合并两个测序结果
def merge_sequencing_result(plasmid_sequencing_primer_df1,plasmid_sequencing_primer_df2):

     #获取uha_dha测序引物的最后一个引物数字
    last_num = (plasmid_sequencing_primer_df1.columns[-1])[-4:-3]
    new_columns = ['Region']
    for i in plasmid_sequencing_primer_df2.columns[1::2]:
        last_num = int(last_num) + 1
        new_columns.append(i[:-1] + str(last_num))  
        new_columns.append(i[:-1] + str(last_num) + '_TM')
    plasmid_sequencing_primer_df2.columns = new_columns
    plasmid_sequencing_primer_df = pd.merge(plasmid_sequencing_primer_df1,plasmid_sequencing_primer_df2,on='Region')
    return plasmid_sequencing_primer_df

def sgRNA_sgRNAprimer_merge(uha_dha_sgRNA_df,
                            sgRNA_plasmid_primer_df,
                            sgRNA_columns=['Name','Region','Rev Target sequence'],
                            sgRNA_primer_columns=['Region',"primer_f_seq_(5'-3')","primer_r_seq_(5'-3')","product_value"]):
    
    temp_sgRNA_df = uha_dha_sgRNA_df[sgRNA_columns]
    temp_sgRNA_df['Region'] = uha_dha_sgRNA_df['Name'] +';'+ uha_dha_sgRNA_df['Region']
    temp_sgRNA_df.drop(columns = 'Name', inplace=True)
    
    sgRNA_plasmid_primer_df = sgRNA_plasmid_primer_df[sgRNA_primer_columns]
    
    
    #纵向堆积----->横向堆积
    index = list(sgRNA_plasmid_primer_df['Region'])
    sgRNA_plasmid_primer_df.drop(columns='Region',inplace=True)
    sgRNA_plasmid_primer_df.index = index
    a = sgRNA_plasmid_primer_df.unstack().to_frame().T
    columns = [i[0] +'_'+ str(i[1]) for i in a.columns]
    a.columns = columns
    
    sgRNA_plasmid_primer_joint_df = su.merge_fill_two_df(temp_sgRNA_df,a)
    return sgRNA_plasmid_primer_joint_df

def add_joint_plasmid_primer(enzyme_df,enzyme_name,sgRNA_plasmid_primer_joint_df,primers_sum,primer_type='sgRNA'):
    #取出去除第一个引物和最后一个引物的剩余引物
    columns = sgRNA_plasmid_primer_joint_df.columns
    primer_columns=[i for i in columns if 'primer' in i]
   
    print(f"primer_r_seq_(5'-3')_{primers_sum}")
    primer_columns = [i for i in primer_columns if i != "primer_f_seq_(5'-3')_1" and i != f"primer_r_seq_(5'-3')_{primers_sum}" ]
    primer_columns.sort(key=lambda x:x[-1])
    primer_columns_1 = [primer_columns[i] for i in range(len(primer_columns)) if i%2==0]
    primer_columns_2 = [primer_columns[i] for i in range(len(primer_columns)) if i%2!=0]
    
    #添加接头
    # enzyme_name = 'BsaI'

    sgRNA_enzyme_df = enzyme_df[enzyme_df['name']==enzyme_name]
    sgRNA_enzyme_df = sgRNA_enzyme_df.reset_index(drop=True)
    
    protective_base = sgRNA_enzyme_df.loc[0,'protective_base']
    recognition_seq = sgRNA_enzyme_df.loc[0,'recognition_seq']
    cut_seq_len = sgRNA_enzyme_df.loc[0,'cut_seq_len']
    gap_len = sgRNA_enzyme_df.loc[0,'gap_len']
    gap_seq = 'AGACTAGACTAGACTAGACTAGACTAGACTAGACTAGACTAGACT'
    gap_seq = sgRNA_enzyme_df.loc[0,'gap_seq']

    left_temp_seq = protective_base + recognition_seq + gap_seq[:gap_len]
    right_temp_seq = su.revComp(protective_base + recognition_seq + gap_seq[:gap_len])[::-1]
    
    for i,v in sgRNA_plasmid_primer_joint_df.iterrows():
        if primer_type == 'sgRNA':
            v["primer_f_seq_(5'-3')_1"] = left_temp_seq + su.revComp(v['Rev Target sequence'])[-cut_seq_len:] + v["primer_f_seq_(5'-3')_1"]
        elif primer_type == 'ccdb':
            v["primer_f_seq_(5'-3')_1"] = left_temp_seq + v["primer_f_seq_(5'-3')_1"]

        if 'r_r' in primer_columns_1[0]:
            for r,f in zip(primer_columns_1, primer_columns_2):
                v[r] = right_temp_seq + su.revComp(v[f][:cut_seq_len]) + v[r]
                print( v[f] ,su.revComp(v[f][:cut_seq_len]))
                v[f] = left_temp_seq + v[f]
        else:
            for r,f in zip(primer_columns_2, primer_columns_1):
                v[r] = right_temp_seq + su.revComp(v[f][:cut_seq_len]) + v[r]  
                print(v[f] ,su.revComp(v[f][:cut_seq_len]))   
                v[f] = left_temp_seq + v[f]

        if primer_type == 'sgRNA':
            v[f"primer_r_seq_(5'-3')_{primers_sum}"] = right_temp_seq + v['Rev Target sequence'] + v[f"primer_r_seq_(5'-3')_{primers_sum}"]
        elif primer_type == 'ccdb':
            v[f"primer_r_seq_(5'-3')_{primers_sum}"] = right_temp_seq + v[f"primer_r_seq_(5'-3')_{primers_sum}"]
    return sgRNA_plasmid_primer_joint_df


def check_locate_primer(sgRNA_plasmid_backbone, primer_json):
    sgRNA_plasmid_backbone = sgRNA_plasmid_backbone.upper()
    
    primer_position_json={}
    failture_primer={}
    for k,v in primer_json.items():
        #引物特异性、存在性验证
        if sgRNA_plasmid_backbone.count(v.upper()) != 1:
            failture_primer.update({k:v})
        else:
        #确定引物的位置
            start = sgRNA_plasmid_backbone.find(v.upper())
            primer_position_json.update({k:start})
    return primer_position_json, failture_primer


def first_left_last_right_primer_design(gb_path, ccdb_label, promoter_terminator_label, n_20_label, last_primer_num, plasmid_backbone=''):

    primer_dict = {}
    failture_primer_dict={}

    if plasmid_backbone != '':
        primer_result = su.primer_design(seqId = 'first_last', seqTemplate = plasmid_backbone, primer_type='plasmid',stype = 'left_right')
        if primer_result.get('PRIMER_LEFT_0_SEQUENCE')==None or primer_result.get('PRIMER_RIGHT_0_SEQUENCE')==None: 
            failture_primer_dict['ID'] = 'plasmid_first_last'
            failture_primer_dict.update(primer_result)
        else:
            first_left_primer = primer_result["PRIMER_LEFT_0_SEQUENCE"]
            last_right_primer = primer_result["PRIMER_RIGHT_0_SEQUENCE"]
            primer_dict.update({f"primer_f_seq_(5'-3')_1":first_left_primer})
            primer_dict.update({f"primer_r_seq_(5'-3')_{last_primer_num}":last_right_primer})
            print('一个质粒系统应用场景')
        return  primer_dict, failture_primer_dict  
    
    else:
        before_processed_seq_dict, after_processed_seq_dict = fq.get_data_from_genebank(gb_path,marker=ccdb_label,target_gene=promoter_terminator_label)
        promoter_terminator = before_processed_seq_dict['target_gene_seq']
        promoter_terminator_up_seq = before_processed_seq_dict['target_gene_up_seq']
        promoter_terminator_down_seq = before_processed_seq_dict['target_gene_down_seq']
        
        #N20
        gb = SeqIO.read(gb_path, "genbank")
        gb_seq = str(gb.seq)
        n20_coordinate = su.get_feature_coordinate(n_20_label, gb_path)

        if n20_coordinate[0] != -1:
            #无ccdb，只有sgRNA的情况
            n20_coordinate_seq = gb_seq[n20_coordinate[0]:n20_coordinate[1]]
            promoter_seq = promoter_terminator[:promoter_terminator.find(n20_coordinate_seq)]
        
            #设计第一引物左引物
            terminator_seq = promoter_terminator[promoter_terminator.find(promoter_seq) + len(promoter_seq) + 20 :]
            first_primer_template = terminator_seq + promoter_terminator_down_seq
            primer_result = su.primer_design(seqId = 'first', seqTemplate = first_primer_template, stype = 'right',primer_type='plasmid')

            if primer_result.get('PRIMER_LEFT_0_SEQUENCE')==None or primer_result.get('PRIMER_RIGHT_0_SEQUENCE')==None: 
                failture_primer_dict['ID'] = 'plasmid_first_left'
                failture_primer_dict.update(primer_result)
            else:
                first_left_primer = primer_result["PRIMER_LEFT_0_SEQUENCE"]
                primer_dict.update({f"primer_f_seq_(5'-3')_1":first_left_primer})

            #设计最一引物右引物
            last_primer_template = promoter_terminator_up_seq + promoter_seq
            primer_result = su.primer_design(seqId = 'last', seqTemplate = last_primer_template, stype = 'left',primer_type='plasmid')
            if primer_result.get('PRIMER_LEFT_0_SEQUENCE')==None or primer_result.get('PRIMER_RIGHT_0_SEQUENCE')==None: 
                failture_primer_dict['ID'] = 'plasmid_last_right'
                failture_primer_dict.update(primer_result)
            else:
                last_right_primer = primer_result["PRIMER_RIGHT_0_SEQUENCE"]
                primer_dict.update({f"primer_r_seq_(5'-3')_{last_primer_num}":last_right_primer})
        else:
            #有ccdb的情况
            ccdb_plasmid_backbone = promoter_terminator_down_seq + promoter_terminator_up_seq
            primer_result = su.primer_design(seqId = 'fist_last', seqTemplate = ccdb_plasmid_backbone, stype = 'left_right',primer_type='plasmid')

            if primer_result.get('PRIMER_LEFT_0_SEQUENCE')==None or primer_result.get('PRIMER_RIGHT_0_SEQUENCE')==None: 
                failture_primer_dict['ID'] = 'plasmid_first_last'
                failture_primer_dict.update(primer_result)
            else:
                first_left_primer = primer_result["PRIMER_LEFT_0_SEQUENCE"]
                last_right_primer = primer_result["PRIMER_RIGHT_0_SEQUENCE"]
                primer_dict.update({f"primer_f_seq_(5'-3')_1":first_left_primer})
                primer_dict.update({f"primer_r_seq_(5'-3')_{last_primer_num}":last_right_primer})
        return primer_dict, failture_primer_dict

#对引物进行排序：
    #1.确定sgRNA_promoter_terminator_start的位置坐标
    #2.优先排列引物所在位置的起始坐标 》sgRNA_promoter_terminator_start的位置坐标
    #3.其次排列引物所在位置的起始坐标 《 sgRNA_promoter_terminator_start的位置坐标
    #4.最终确定用户的引物排列顺序
    
def sort_compose_primer(sgRNA_promoter_terminator_start,
                        primer_json,
                        primer_position_json,
                        gb_path,
                        sgRNA_plasmid_backbone,
                        n_20_label='N20',
                        ccdb_label='ccdB',
                        promoter_terminator_label = 'gRNA',
                        plasmid_type = 'two'
                        ):
    
    
    position = list(primer_position_json.values())
    first_list =  sorted([i for i in position if i > sgRNA_promoter_terminator_start])
    second_list =  sorted([i for i in position if i < sgRNA_promoter_terminator_start])
    first_list.extend(second_list)
    position_primer_json = {value: key for key, value in primer_position_json.items()}
    
    #引物数量
    primers_sum = len(first_list)*2 + 2
    primer_dict = {}
    i = 1


    failture_primer_dict_df = pd.DataFrame()  

    if plasmid_type == 'two':
        first_last_primer_dict, failture_first_last_primer_dict = first_left_last_right_primer_design(
                                            gb_path,
                                            ccdb_label,
                                            promoter_terminator_label,
                                            n_20_label,
                                            last_primer_num = int(primers_sum/2))
    else:
        first_last_primer_dict, failture_first_last_primer_dict = first_left_last_right_primer_design(
                                            gb_path,
                                            ccdb_label,
                                            promoter_terminator_label,
                                            n_20_label,
                                            last_primer_num = int(primers_sum/2),
                                            plasmid_backbone = sgRNA_plasmid_backbone
                                            )

      
    if failture_first_last_primer_dict != {}:
        failture_primer_dict_df =  pd.DataFrame([failture_first_last_primer_dict])
        return 'failture', failture_primer_dict_df,0

    primer_dict.update(first_last_primer_dict)
    #在质粒上取引物,同时设计引物
    for item in first_list:
        left_primer = primer_json[position_primer_json[item]]
        if item < 40:
            template = sgRNA_plasmid_backbone[-100:]+sgRNA_plasmid_backbone[:item]
        else:
            template = sgRNA_plasmid_backbone[:item]
        primer_result =  su.primer_design(seqId=item,seqTemplate=template, primer_type='plasmid',stype='right')
        if primer_result.get('PRIMER_RIGHT_0_SEQUENCE')==None:
            failture_dict = {}
            failture_dict['ID'] = position_primer_json[item]
            failture_dict.update(primer_result)
            failture_primer_dict_df = pd.DataFrame([failture_dict]) 
            return 'failture', failture_primer_dict_df,0
        else:
            right_primer = primer_result['PRIMER_RIGHT_0_SEQUENCE']
            primer_dict.update({f"primer_r_seq_(5'-3')_{i}":right_primer})
            i = i + 1
            primer_dict.update({f"primer_f_seq_(5'-3')_{i}":left_primer})
      
    primer_dict_df = pd.DataFrame([primer_dict])
    
    return 'success',primer_dict_df, primers_sum

def plasmid_primer(sgRNA_plasmid_primer_joint_df):
    sgRNA_columns = [i for i in sgRNA_plasmid_primer_joint_df.columns if 'primer' in i]
    df = pd.DataFrame()
    for i,row in sgRNA_plasmid_primer_joint_df.iterrows():
        df = df.append(su.df_to_df(sgRNA_columns, sgRNA_plasmid_primer_joint_df, i))
#     df = df.drop_duplicates(subset = ["primer_f_seq_(5'-3')","primer_r_seq_(5'-3')"])
    return df

def add_product_and_size(gb_path,primer_df,enzyme_df,enzyme_name='BsaI',seq=''):
       
    sgRNA_enzyme_df = enzyme_df[enzyme_df['name']==enzyme_name]
    sgRNA_enzyme_df = sgRNA_enzyme_df.reset_index(drop=True)
    protective_base = sgRNA_enzyme_df.loc[0,'protective_base']
    recognition_seq = sgRNA_enzyme_df.loc[0,'recognition_seq']
    cut_seq_len = sgRNA_enzyme_df.loc[0,'cut_seq_len']
    gap_len = sgRNA_enzyme_df.loc[0,'gap_len']
    gap_seq = 'AGACTAGACTAGACTAGACTAGACTAGACTAGACTAGACTAGACT'
    gap_seq = sgRNA_enzyme_df.loc[0,'gap_seq']

    joint_len  = len(protective_base) + len(recognition_seq) + gap_len + cut_seq_len
    print('接头长度：',joint_len)  
    if seq != '':
        gb_seq = seq.upper()
    else:
        gb = SeqIO.read(gb_path, "genbank")
        gb_seq = str(gb.seq)
        gb_seq = gb_seq.upper()

       # no_ccdb_plasmid, no_sgRNA_plasmid
    def work(f_primer,r_primer):    
        left_joint_seq = f_primer[:joint_len]
        f_primer = f_primer[joint_len:]
        right_joint_seq = r_primer[:joint_len]       
        right_joint_seq = su.revComp(right_joint_seq)
    
        n20 = ''
        if len(r_primer) >= 40:
            n20 = r_primer[joint_len:joint_len+20]
            r_primer = r_primer[joint_len+20:]
        else:
            r_primer = r_primer[joint_len:]
            
        r_primer = su.revComp(r_primer)    
        start = gb_seq.find(f_primer)
        end = gb_seq.find(r_primer)+len(r_primer)
        
        
        if end < start:
            product = gb_seq[start:] + gb_seq[:end]
        else:
            product = gb_seq[start:end]
        
        #加接头
        if n20 != '':
            product = left_joint_seq + product + su.revComp(n20) + right_joint_seq
        else:
            product = left_joint_seq + product + right_joint_seq
        
        
        return product, len(product)

    df = su.lambda2cols(primer_df, work, in_coln=["primer_f_seq_(5'-3')_joint", "primer_r_seq_(5'-3')_joint"], to_colns=["product_value_joint","product_size_joint"])
    return df
 
def create_enzymeCutSeq_and_N20(temp_sgRNA_df,promoter_seq,enzyme_cut_len=4):
    def work(sgRNA_n20, promoter_N20_terminator):
        terminator_seq = promoter_N20_terminator[len(promoter_seq)+20:]
        return (promoter_seq[-enzyme_cut_len:] + sgRNA_n20 + terminator_seq[:enzyme_cut_len])
    df = su.lambda2cols(df=temp_sgRNA_df, lambdaf=work, in_coln=['Target sequence','promoter_N20_terminator'], to_colns=['enzymeCutSeq_and_N20'])
    return df

from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
def create_gb_file(plasmid,primer_cor_dict):
    plasmid_len = len(plasmid)
    li_feature = []
    
     # 创建一个SeqRecord对象
    record = SeqRecord(Seq(plasmid),
                   id="plamid",
                   name="plasmid",
                   description="primer seq in plasmid")
    
    print('3123uha检查',primer_cor_dict)

    for k,v in primer_cor_dict.items():
        print(k,v)
        cor = v.split(';')[0]
        seq = v.split(';')[1]
        if cor != '()' or seq != '':
            cor = cor.replace('(','').replace(')',"")
            arr = cor.split(',')
            print(arr[0],arr[1])
            start, end = int(arr[0]), int(arr[1])
            temp_feature = ''
            
            if 'LEFT' in k:
                if start < 0:
                    start = plasmid_len - abs(start)
                    print('oooooooooooooooooooooo',start,end)
                    locations = [FeatureLocation(start, plasmid_len,+1), FeatureLocation(0, end,+1)]
                    temp_feature =  SeqFeature(CompoundLocation(locations), type='primer_bind',strand=+1,qualifiers={"note":k,"primer": f"{seq}"})
                else:
                    temp_feature = SeqFeature(FeatureLocation(start, end), type="primer_bind",strand=+1,qualifiers={"note":k,"primer": f"{seq}"})
                record.features.append(temp_feature)
            elif 'RIGHT' in k:
                if end > plasmid_len:
                    end = end - plasmid_len
                    print('ppppppppppppp',start,end)
                    locations = [FeatureLocation(0, end,-1), FeatureLocation(start, plasmid_len,-1)]
                    temp_feature = SeqFeature(CompoundLocation(locations), type='primer_bind',strand=-1,qualifiers={"note":k,"primer": f"{seq}"})
                else:
                    temp_feature = SeqFeature(FeatureLocation(start, end), type="primer_bind",strand=-1,qualifiers={"note":k,"primer": f"{seq}"})
            
                li_feature.append(temp_feature)
                record.features.append(temp_feature)
            elif 'SEQUENCING' in k:
                # if start < 0:
                #     start = plasmid_len - abs(start)  
                if 'SEQUENCING_PRIMER_2' == k or 'right' in k:
                    if 'right' in k:
                        k, right = k.split(';')
                    temp_feature = SeqFeature(FeatureLocation(start, end), type="primer_bind",strand=-1,qualifiers={"note":k,"primer": f"{seq}"})
                else:
                    temp_feature = SeqFeature(FeatureLocation(start, end), type="primer_bind",strand=+1,qualifiers={"note":k,"primer": f"{seq}"})
                li_feature.append(temp_feature)
                record.features.append(temp_feature)



            elif 'UHA' in k:
                temp_feature = SeqFeature(FeatureLocation(start, end), type="cds",strand=+1,qualifiers={"note":k})
                li_feature.append(temp_feature)
                record.features.append(temp_feature)
            elif 'DHA' in k:
                temp_feature = SeqFeature(FeatureLocation(start, end), type="cds",strand=+1,qualifiers={"note":k})
                li_feature.append(temp_feature)
                record.features.append(temp_feature)
            elif 'SEQ_ALTERED' in k:
                temp_feature = SeqFeature(FeatureLocation(start, end), type="cds",strand=+1,qualifiers={"note":k})  
                li_feature.append(temp_feature)
                record.features.append(temp_feature)
            elif 'PROMOTER_N20_TERMINATOR' in k:
                temp_feature = SeqFeature(FeatureLocation(start, end), type="cds",strand=+1,qualifiers={"note":k})
                li_feature.append(temp_feature)
                record.features.append(temp_feature)
            elif 'N20' in k:
                temp_feature = SeqFeature(FeatureLocation(start, end), type="cds",strand=+1,qualifiers={"note":k})
                li_feature.append(temp_feature)
                record.features.append(temp_feature)
                


    # 添加molecule_type信息到annotations
    record.annotations['molecule_type'] = 'genomic DNA'            
#     SeqIO.write(record, "my_seq.gb", "genbank")
    return record

   
def add_sequencing_primer_to_gb(plasmid_sequencing_primer_df,temp_dict):
    plasmid = temp_dict.get('PLASMID')
    id = temp_dict.get('ID')
    uha = temp_dict.get('UHA')
    dha = temp_dict.get('DHA')
    # seq_altered = temp_dict.get('SEQ_ALTERED')
    pro_n20_ter = temp_dict.get('PROMOTER_N20_TERMINATOR')

    sequencing_primer_dict={}
    sequencing_primer_temp = plasmid_sequencing_primer_df[plasmid_sequencing_primer_df['ID'] == id]
    sequencing_primer_temp.reset_index(drop=True,inplace=True)
    sequencing_primer_temp = sequencing_primer_temp.fillna('')

    
    # if seq_altered != '':
    #     seq_altered_cor = su.create_primerCor_in_plasmid(plasmid, seq_altered)
    #     sequencing_primer_dict.update({'SEQ_ALTERED': f'{seq_altered_cor};{seq_altered}'})
    if uha != '' and uha != None:
        uha_cor = su.create_primerCor_in_plasmid(plasmid, uha)
        sequencing_primer_dict.update({'UHA': f'{uha_cor};{uha}'})
    if dha != '' and dha != None:
        dha_cor = su.create_primerCor_in_plasmid(plasmid, dha)
        sequencing_primer_dict.update({'DHA':f'{dha_cor};{dha}'})
    if pro_n20_ter !='' and pro_n20_ter != None: 
        pro_n20_ter_cor = su.create_primerCor_in_plasmid(plasmid, pro_n20_ter)
        sequencing_primer_dict.update({'PROMOTER_N20_TERMINATOR':f'{pro_n20_ter_cor};{pro_n20_ter}'})
    
   
    k = 1
    for i in sequencing_primer_temp.columns:
        if 'PRIMER' in i:
            if i =='ID' or 'TM' in i:
                continue 
            else:
                if len(i.split('_')) == 3:           
                    arr1 = i.split('_')[0]
                    arr2 = i.split('_')[1]
                    arr3 = i.split('_')[2]
                    primer_seq = sequencing_primer_temp.loc[0,i]
                    if primer_seq == '':
                        continue
                    key = '_'.join([arr1, arr2, str(k)])
                    #唯一一条负义链的引物
                    if arr3 == '2':
                        primer_cor = su.create_primerCor_in_plasmid(plasmid, su.revComp(primer_seq))
                        sequencing_primer_dict.update({key:f'{primer_cor};{primer_seq}'})     
                    else:
                        primer_cor = su.create_primerCor_in_plasmid(plasmid, primer_seq)

                        if primer_cor[0] == -1:  
                            primer_cor = su.create_primerCor_in_plasmid(plasmid, su.revComp(primer_seq))
                            sequencing_primer_dict.update({key+';right':f'{primer_cor};{primer_seq}'})
                        else:
                            sequencing_primer_dict.update({key:f'{primer_cor};{primer_seq}'})

                    k = k + 1
    return  sequencing_primer_dict  


def create_gb_for_region(plasmid_primer_featrue_df, n20down_primer_p_df, joint_len, cut_seq_len, output, type='sgRNA_ccdb'):

   
    #处理特殊的质粒片段引物：不一定几对引
    
    df = pd.DataFrame()
   
    if type == 'sgRNA':
        
        # plasmid_primer_df = pd.merge(plasmid_primer_featrue_df,primer_df)
        for i,v in plasmid_primer_featrue_df.iterrows():
            plasmid  = v['PLASMID'].upper()
            ID = v['ID']
            #在质粒上找到位置坐标
            primer_df = su.groupby_columns_to_row(n20down_primer_p_df)

            primer_cor_dict={}   
            primer_dict = {}

            promoter_N20_terminator = v['PROMOTER_N20_TERMINATOR']
            promoter_N20_terminator_cor = su.create_primerCor_in_plasmid(plasmid, promoter_N20_terminator)
            primer_cor_dict.update({'PROMOTER_N20_TERMINATOR': f"{promoter_N20_terminator_cor};{promoter_N20_terminator}"})  

            for k in primer_df.columns:
                if k == 'ID':
                    continue
                value = primer_df.loc[i,k]
                temp_value = value[joint_len-cut_seq_len:]

                if 'RIGHT' in k:
                    primer_cor = su.create_primerCor_in_plasmid(plasmid, su.revComp(temp_value))
                    right_primer_cor = primer_cor[0], primer_cor[1] + (joint_len-cut_seq_len)
                    primer_cor_dict.update({k: f"{right_primer_cor};{value}"})
                    primer_dict.update({k:value})    
                else:
                    primer_cor = su.create_primerCor_in_plasmid(plasmid,temp_value)
                    left_primer_cor = primer_cor[0] - (joint_len-cut_seq_len), primer_cor[1]
                    primer_cor_dict.update({k:f'{left_primer_cor};{value}'})
                    primer_dict.update({k:value})


            gb_record = create_gb_file(plasmid,primer_cor_dict)
            name=v['ID'].split(';')[0]
          
            su.write_gb(gb_record, output_path = output, gb_name=type+'_'+name, gb_type='genbank')

            temp = pd.DataFrame(columns=['name',type+'_gb'],data=[[name, type+'_'+name+ '.gb']])
            df = df.append(temp)

    elif type == 'enzyme_cut':
        for i,v in plasmid_primer_featrue_df.iterrows():
            plasmid  = v['PLASMID'].upper()
            ID = v['ID']

            primer_cor_dict={}   
          
            promoter_N20_terminator = v['PROMOTER_N20_TERMINATOR']
            promoter_N20_terminator_cor = su.create_primerCor_in_plasmid(plasmid, promoter_N20_terminator)
            primer_cor_dict.update({'PROMOTER_N20_TERMINATOR': f"{promoter_N20_terminator_cor};{promoter_N20_terminator}"}) 

            enzymeCutSeq_and_N20 = n20down_primer_p_df
            temp = enzymeCutSeq_and_N20[enzymeCutSeq_and_N20['ID']==ID].reset_index(drop=True)
            n20 = temp.loc[0,'Target sequence']
            n20_cor = su.create_primerCor_in_plasmid(plasmid, n20)
            primer_cor_dict.update({'N20':f'{n20_cor};{n20}'})

            gb_record = create_gb_file(plasmid,primer_cor_dict)
            name=v['ID'].split(';')[0]

            type = 'sgRNA' 
            su.write_gb(gb_record, output_path = output, gb_name=type+'_'+name, gb_type='genbank')

            temp = pd.DataFrame(columns=['name',type+'_gb'],data=[[name, type+'_'+name+ '.gb']])
            df = df.append(temp)


    elif type == 'sgRNA_ccdb' or type == 'ccdb': 
        
        print(n20down_primer_p_df)
        n20down_primer_temp_df = su.columns_to_row(n20down_primer_p_df)
        print(n20down_primer_temp_df)
        
        for i,v in plasmid_primer_featrue_df.iterrows():

            ID = v['ID']
            plasmid  = v['PLASMID'].upper()

            uha_left_primer = v['UHA_PRIMER_LEFT_WHOLE_SEQUENCE']
            uha_left_nojoint_primer = uha_left_primer[joint_len:].upper()                                                                                                                                                          
            uha_left_primer_cor = su.create_primerCor_in_plasmid(plasmid,uha_left_nojoint_primer)
            uha_left_primer_cor = uha_left_primer_cor[0] - joint_len, uha_left_primer_cor[1]

            uha = v['UHA']

            uha_cor = su.create_primerCor_in_plasmid(plasmid, uha)

            uha_right_primer = v['UHA_PRIMER_RIGHT_WHOLE_SEQUENCE']
            uha_right_nojoint_primer = uha_right_primer[joint_len:].upper()
            uha_right_primer_cor = su.create_primerCor_in_plasmid(plasmid, su.revComp(uha_right_nojoint_primer))
            uha_right_primer_cor = uha_right_primer_cor[0], uha_right_primer_cor[1] + joint_len
           

            dha_left_primer = v['DHA_PRIMER_LEFT_WHOLE_SEQUENCE']
            dha_left_nojoint_primer = dha_left_primer[joint_len:].upper()
            dha_left_primer_cor = su.create_primerCor_in_plasmid(plasmid, dha_left_nojoint_primer)
            dha_left_primer_cor = dha_left_primer_cor[0] - joint_len, dha_left_primer_cor[1]  

            dha = v['DHA']
            dha_cor = su.create_primerCor_in_plasmid(plasmid, dha)  
            dha_right_primer = v['DHA_PRIMER_RIGHT_WHOLE_SEQUENCE']
            dha_right_nojoint_primer = dha_right_primer[joint_len:].upper()
            dha_right_primer_cor = su.create_primerCor_in_plasmid(plasmid, su.revComp(dha_right_nojoint_primer))
            dha_right_primer_cor = dha_right_primer_cor[0] - joint_len, dha_right_primer_cor[1]
            
            if 'SEQ_PRIMER_LEFT_WHOLE_SEQUENCE' in plasmid_primer_featrue_df.columns:

                seq_altered_left_primer = v['SEQ_PRIMER_LEFT_WHOLE_SEQUENCE']
                seq_altered_left_primer_cor = ()
                if seq_altered_left_primer != '':
                    seq_altered_left_nojoint_primer = seq_altered_left_primer[joint_len-cut_seq_len:].upper()
                    seq_altered_left_primer_cor = su.create_primerCor_in_plasmid(plasmid, seq_altered_left_nojoint_primer)
                    seq_altered_left_primer_cor = seq_altered_left_primer_cor[0] - (joint_len-cut_seq_len), seq_altered_left_primer_cor[1]

                seq_altered_right_primer =  v['SEQ_PRIMER_RIGHT_WHOLE_SEQUENCE']   
                seq_altered_right_primer_cor = ()
                if seq_altered_right_primer != '':   
                    seq_altered_right_nojoint_primer = seq_altered_right_primer[joint_len-cut_seq_len:].upper()
                    seq_altered_right_primer_cor = su.create_primerCor_in_plasmid(plasmid, su.revComp(seq_altered_right_nojoint_primer))
                    seq_altered_right_primer_cor = seq_altered_right_primer_cor[0], seq_altered_right_primer_cor[1] + (joint_len-cut_seq_len)
            else:
                seq_altered_left_primer= ""
                seq_altered_left_primer_cor=""
                seq_altered_right_primer=""
                seq_altered_right_primer_cor=""

            if type == 'sgRNA_ccdb':
                n20up_left_primer = v['N20UP_PRIMER_LEFT_WHOLE_SEQUENCE']
                n20up_left_nojoint_primer = n20up_left_primer[joint_len-cut_seq_len:].upper()
                n20up_left_primer_cor = su.create_primerCor_in_plasmid(plasmid, n20up_left_nojoint_primer)
                n20up_left_primer_cor = n20up_left_primer_cor[0]-(joint_len-cut_seq_len), n20up_left_primer_cor[1]

                n20up_right_primer = v['N20UP_PRIMER_RIGHT_WHOLE_SEQUENCE']
                n20up_right_nojoint_primer = n20up_right_primer[joint_len-cut_seq_len+20:].upper()
                n20up_right_primer_cor = su.create_primerCor_in_plasmid(plasmid, su.revComp(n20up_right_nojoint_primer))
                n20up_right_primer_cor = n20up_right_primer_cor[0], n20up_right_primer_cor[1]+(joint_len-cut_seq_len+20)

            if 'PROMOTER_N20_TERMINATOR' in plasmid_primer_featrue_df.columns:
                promoter_N20_terminator = v['PROMOTER_N20_TERMINATOR']
                promoter_N20_terminator_cor = su.create_primerCor_in_plasmid(plasmid, promoter_N20_terminator)
  
            #在质粒上找到位置坐标
            n20down_primer_temp_df = su.columns_to_row(n20down_primer_p_df)
            n20down_primer_cor_dict={}
            n20down_primer_dict = {}
            for i in n20down_primer_temp_df.columns:
                value = n20down_primer_temp_df.loc[0,i]
                temp_value = value[joint_len-cut_seq_len:]

                if 'RIGHT' in i:
                    primer_cor = su.create_primerCor_in_plasmid(plasmid, su.revComp(temp_value))
                    n20down_right_primer_cor = primer_cor[0], primer_cor[1] + (joint_len-cut_seq_len)
                    n20down_primer_cor_dict.update({i: f"{n20down_right_primer_cor};{value}"})
                    n20down_primer_dict.update({i:value})
                    
                elif 'LEFT'in i:
                    primer_cor = su.create_primerCor_in_plasmid(plasmid, temp_value)
                    if primer_cor[0] == -1:
                        primer_cor = su.create_primerCor_in_plasmid(plasmid,temp_value[cut_seq_len:])
                        n20down_left_primer_cor = primer_cor[0] - (joint_len), primer_cor[1] 
                    else:
                        n20down_left_primer_cor = primer_cor[0] - (joint_len-cut_seq_len), primer_cor[1]

                    n20down_primer_cor_dict.update({i:f'{n20down_left_primer_cor};{value}'})
                    n20down_primer_dict.update({i:value})


            if type == 'sgRNA_ccdb':
                primer_cor_dict = {
                    'UHA_PRIMER_LEFT_WHOLE_SEQUENCE':f'{uha_left_primer_cor};{uha_left_primer}',
                    'UHA':f'{uha_cor};{uha}',
                    'UHA_PRIMER_RIGHT_WHOLE_SEQUENCE':f'{uha_right_primer_cor};{uha_right_primer}',
                    'DHA_PRIMER_LEFT_WHOLE_SEQUENCE':f'{dha_left_primer_cor};{dha_left_primer}',
                    'DHA':f'{dha_cor};{dha}',
                    'DHA_PRIMER_RIGHT_WHOLE_SEQUENCE':f'{dha_right_primer_cor};{dha_right_primer}',
                    'N20UP_PRIMER_LEFT_WHOLE_SEQUENCE':f'{n20up_left_primer_cor};{n20up_left_primer}',
                    'N20UP_PRIMER_RIGHT_WHOLE_SEQUENCE':f"{n20up_right_primer_cor};{n20up_right_primer}",
                    'PROMOTER_N20_TERMINATOR':f'{promoter_N20_terminator_cor};{promoter_N20_terminator}'
                }

            elif type == 'ccdb':   
                primer_cor_dict = {
                    'UHA_PRIMER_LEFT_WHOLE_SEQUENCE':f'{uha_left_primer_cor};{uha_left_primer}',
                    'UHA':f'{uha_cor};{uha}',
                    'UHA_PRIMER_RIGHT_WHOLE_SEQUENCE':f'{uha_right_primer_cor};{uha_right_primer}',
                    'DHA_PRIMER_LEFT_WHOLE_SEQUENCE':f'{dha_left_primer_cor};{dha_left_primer}',
                    'DHA':f'{dha_cor};{dha}',
                    'DHA_PRIMER_RIGHT_WHOLE_SEQUENCE':f'{dha_right_primer_cor};{dha_right_primer}'
                }
                # primer_dict = {
                #     'UHA_PRIMER_LEFT_WHOLE_SEQUENCE':f'{uha_left_primer}',
                #     'UHA_PRIMER_RIGHT_WHOLE_SEQUENCE':f'{uha_right_primer}',
                #     'SEQ_PRIMER_LEFT_WHOLE_SEQUENCE':f'{seq_altered_left_primer}',
                #     'SEQ_PRIMER_RIGHT_WHOLE_SEQUENCE':f'{seq_altered_right_primer}',
                #     'DHA_PRIMER_LEFT_WHOLE_SEQUENCE':f'{dha_left_primer}',
                #     'DHA_PRIMER_RIGHT_WHOLE_SEQUENCE':f'{dha_right_primer}'
                # }
            
            if  seq_altered_left_primer != "":
                primer_cor_dict.update(
                    {
                        'SEQ_PRIMER_LEFT_WHOLE_SEQUENCE':f'{seq_altered_left_primer_cor};{seq_altered_left_primer}',
                        'SEQ_PRIMER_RIGHT_WHOLE_SEQUENCE':f'{seq_altered_right_primer_cor};{seq_altered_right_primer}',
                        'SEQ_PRIMER_LEFT_WHOLE_SEQUENCE':f'{seq_altered_left_primer_cor};{seq_altered_left_primer}',
                        'SEQ_PRIMER_RIGHT_WHOLE_SEQUENCE':f'{seq_altered_right_primer_cor};{seq_altered_right_primer}'  
                    }

                )

            #添加质粒片段引物
            print('---------------------------------------------------------------',n20down_primer_cor_dict)
            primer_cor_dict.update(n20down_primer_cor_dict) 

            print(primer_cor_dict,'\n')
            gb_record = create_gb_file(plasmid,primer_cor_dict)

            name=v['ID'].split(';')[0]

            su.write_gb(gb_record, output_path = output, gb_name=type+'_'+name, gb_type='genbank')
            #生成gb文件的csv文件
            temp = pd.DataFrame(columns=['name',type+'_gb'],data=[[name, type+'_'+name+ '.gb']])
            df = df.append(temp)
            
    return df

   
def create_gb_for_sequencing_region(plasmid_primer_featrue_df, plasmid_sequencing_primer_df, output, type='plasmid_sequencing'):


    df = pd.DataFrame()
    columns = list(plasmid_primer_featrue_df.columns)

    for i,v in plasmid_primer_featrue_df.iterrows():
        plasmid = v['PLASMID'].upper()
        temp_dict = {   'ID':v['ID'],
                        'PLASMID': plasmid}
        
        if type == 'genome_sequencing':
            temp_dict.update({
                'UHA':v['UHA'],
                'DHA':v['DHA']
            })
        elif 'sequencing' in type:
            if 'UHA' in columns and 'PROMOTER_N20_TERMINATOR' in columns:
                temp_dict.update({
                    'UHA':v['UHA'],
                    'DHA':v['DHA'],
                    'PROMOTER_N20_TERMINATOR':v['PROMOTER_N20_TERMINATOR']
                })
            elif 'UHA' not in columns and 'PROMOTER_N20_TERMINATOR' in columns:
                temp_dict.update({
                    'PROMOTER_N20_TERMINATOR':v['PROMOTER_N20_TERMINATOR']
                })
            elif 'UHA' in columns and 'PROMOTER_N20_TERMINATOR' not in columns:
                temp_dict.update({
                    'UHA':v['UHA'],
                    'DHA':v['DHA']
                })

        sequencing_primer_cor_dict = add_sequencing_primer_to_gb(plasmid_sequencing_primer_df, temp_dict)

        gb_record = create_gb_file(plasmid,sequencing_primer_cor_dict)
        name=v['ID'].split(';')[0]

        su.write_gb(gb_record, output_path = output, gb_name=type+'_'+name, gb_type='genbank')

        temp = pd.DataFrame(columns=['name',type+'_gb'],data=[[name, type+'_'+name+ '.gb']])
        df = df.append(temp)

    return df
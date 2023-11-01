import sys,os      
from sgRNA_utils.sgRNA_primer_config import config   
import sgRNA_utils.sgRNA_primer_util as util
import pandas as pd
from os.path import exists,splitext,dirname,splitext,basename,realpath,abspath  
from Bio import SeqIO  
import numpy as np 


# config.tmp_path = os.path.join(config.workdir,  'sequencing_primer_alignment')   

def search_primer_region( len_target, path, seq_id, seq_temp, primer_type = 'LEFT', site_target_temp=200, left_primer_region=40, right_primer_region=40, protect_region=80, max_region = 120, step= 40 ):    
     
    sequencinng_region_len = ( left_primer_region - 40 + protect_region + len_target + protect_region + right_primer_region - 40 )

    if primer_type =='LEFT': 
        primer_region = left_primer_region
    elif primer_type =='RIGHT':
        primer_region = right_primer_region

    if max_region == 120:
        region_len = primer_region
    elif max_region >120:
        region_len = sequencinng_region_len

    while region_len <= max_region:

        if primer_type=='LEFT':
            start = site_target_temp - (primer_region + protect_region)
            end = site_target_temp - protect_region  
            primer_region_seq = seq_temp[start : end]
            primer_region_seq = seq_temp[start : start + step]
        elif primer_type=='RIGHT':
            start = site_target_temp + len_target + protect_region
            end = site_target_temp + len_target + (protect_region + primer_region)
            primer_region_seq = seq_temp[start: end]
            primer_region_seq = seq_temp[end - step : end]

        #对引物模板脱靶分析

        #sam_df文件
        index_prefix_path = os.path.join(config.tmp_path, 'genome_index')
        fasta_filename = os.path.join( config.tmp_path, 'sequences.fasta')
        data={
                'id':[seq_id],  
                'sequence':[primer_region_seq]
            }
        df = pd.DataFrame(data)
        util.df_to_fasta(df,fasta_filename)  
    
        sam_df = util.bowtie_seq_genome(config.tmp_path, path, index_prefix_path, fasta_filename)
        if len(sam_df)==1:  
            #不脱靶
            back_up_primers = {'0','khfvdshkh'}
        elif len(sam_df) > 1:
            #脱靶
            back_up_primers = {}   
        if back_up_primers != {}:
            break
        else:
            primer_region = primer_region + step
            if  max_region >120 :
                if  primer_type=='LEFT':
                    sequencinng_region_len = ( primer_region - 40 + protect_region + len_target + protect_region + right_primer_region - 40 )
                elif primer_type == 'RIGHT': 
                    sequencinng_region_len = ( left_primer_region - 40 + protect_region + len_target + protect_region + primer_region - 40 ) 
                region_len = sequencinng_region_len  
            else:  
                region_len = primer_region

    return primer_region, back_up_primers

#生成测序引物模板
def design_sequencing_primers(plasmid_seq, path, dict_plasmid_id, dict_plasmid_seq, mute_position=0, seq_type='', gene_change_type='', ref_seq='', seq_altered=''):
    '''
        生成测序引物模板 (在target_gene_seq 前后加200bp) ,设计测序引物
        params:     
            dict_plasmid_id  
            dict_plasmid_seq
        returns:
    '''

    #生成引物模板
    target_gene_down_seq = dict_plasmid_seq['target_gene_down_seq']
    target_gene_up_seq = dict_plasmid_seq['target_gene_up_seq']
    target_gene_seq = dict_plasmid_seq['mute_after_target_gene_seq']

    threshold = 200 + 1800   

    # target_gene_seq 前后加200bp    
    if len(target_gene_down_seq) >= threshold:
        sequencing_peimers_template = target_gene_seq + target_gene_down_seq[ : threshold]  
    else:
        sequencing_peimers_template = target_gene_seq + target_gene_down_seq
        temp_len = threshold - len(target_gene_down_seq) 
        sequencing_peimers_template = sequencing_peimers_template + target_gene_up_seq[:temp_len]  
    if len(target_gene_up_seq) >= threshold:
        sequencing_peimers_template = target_gene_up_seq[-threshold:] + sequencing_peimers_template
    else:
        sequencing_peimers_template = target_gene_up_seq + sequencing_peimers_template  
        temp_len = threshold - len(target_gene_up_seq)
        sequencing_peimers_template = target_gene_down_seq[-temp_len:] + sequencing_peimers_template

    #设计引物 
    if mute_position ==0:
        result, failture_result = design_seq_primer(gene_change_type, ref_seq, seq_altered, plasmid_seq, path, seq_id=dict_plasmid_id, seq_temp=sequencing_peimers_template, seq_target=target_gene_seq, seq_type=seq_type)
        print(result, failture_result)
    else:
         result, failture_result = design_seq_primer(plasmid_seq, path,seq_id=dict_plasmid_id, seq_temp=sequencing_peimers_template, seq_target=target_gene_seq,mute_position=mute_position,seq_type=seq_type)

    return result, failture_result 

#循环调用primer3
def try_pcr_primer_design(len_target ,left_primer_region, right_primer_region, max_primer_region, protect_region, site_target_temp, seq_id, seq_type, seq_temp ):

    while left_primer_region <= max_primer_region and  right_primer_region <= max_primer_region:

        start = (site_target_temp - (left_primer_region + protect_region) )
        end = start + (left_primer_region + protect_region ) + len_target + (right_primer_region + protect_region)
        pcr_template = seq_temp[ start : end ]
 
        #2.设计
        dict_res_p = util.primer_design(seqId=seq_id, 
                                            seqTemplate=pcr_template,                                           #开始使用
                                            stype='left_right_change',
                                            mute_type='single',
                                            primer_type=seq_type,
                                            global_args=config.Q_GLOBAL_ARGS  
                                        )
        #测序引物设计结果，若设计不出来，扩展引物设计区域
        if dict_res_p.get('PRIMER_LEFT_NUM_RETURNED') == 0 and 'ok 0' in dict_res_p.get('PRIMER_LEFT_EXPLAIN'):
            left_primer_region = left_primer_region + 10
        if dict_res_p.get('PRIMER_RIGHT_NUM_RETURNED') == 0 and 'ok 0' in dict_res_p.get('PRIMER_RIGHT_EXPLAIN'):
            right_primer_region = right_primer_region + 10
        if dict_res_p.get('PRIMER_LEFT_NUM_RETURNED') > 1 and dict_res_p.get('PRIMER_RIGHT_NUM_RETURNED') > 1:
            break
    
    if left_primer_region > max_primer_region or right_primer_region > max_primer_region:
        print("引物区域扩张到120bp,仍然有PCR引物没设计出来,引物id为:",seq_id)

    return dict_res_p, pcr_template  

# def try_sequencing_primer_design( seq_id, seq_temp, site_target_temp, left_primer_region, protect_region, seq_type, max_primer_region):     
def try_sequencing_primer_design( seq_id, sequencing_template, left_primer_region, max_primer_region, protect_region, seq_type):
    while left_primer_region <=  max_primer_region:  
        temp_pn =  sequencing_template[max_primer_region - left_primer_region : max_primer_region]
        dict_res_pn = util.primer_design(seqId=seq_id,
                                                    seqTemplate=temp_pn,  
                                                    stype='none',
                                                    mute_type='sequencing',  
                                                    primer_type=seq_type,
                                                    global_args=config.Q_GLOBAL_ARGS
                                                )

        if dict_res_pn.get('PRIMER_LEFT_NUM_RETURNED') > 0:
            break
        else:
            left_primer_region = left_primer_region + 10
    
    if left_primer_region > max_primer_region :  
        print("引物区域扩张到120bp,仍然有测序引物没设计出来,引物id为:",seq_id)
        errorMessage = f'引物区域扩张到120bp,仍然有测序引物没设计出来,引物id为:{seq_id}'
        raise ValueError(errorMessage)   

    return dict_res_pn

def try_pcr_offtarget_design(path, len_target , seq_target, expend_traget_len , left_primer_region, right_primer_region, max_primer_region, protect_region, site_target_temp, seq_id, seq_type, seq_temp):

        i = 0
        origin_left_primer_region = left_primer_region
        origin_right_primer_region = right_primer_region
        while (left_primer_region - origin_left_primer_region + protect_region + len_target + protect_region + right_primer_region - origin_right_primer_region ) <= expend_traget_len:

            # 0.左、右引物区间是否大于120bp 
            step=40
            if left_primer_region < 120 or right_primer_region < 120:
                max_primer_region = 120
            else:
                max_primer_region = expend_traget_len

            #1.在40-120区间内探测引物模板是否脱靶
            primer_region = left_primer_region  
            left_primer_region, left_back_up_primers = search_primer_region( len_target,
                                                                                path,
                                                                                seq_id,
                                                                                seq_temp,
                                                                                primer_type = 'LEFT',
                                                                                site_target_temp = site_target_temp, 
                                                                                left_primer_region=left_primer_region, 
                                                                                right_primer_region=right_primer_region,
                                                                                protect_region=80, 
                                                                                max_region = max_primer_region, 
                                                                                step=step)
            if left_back_up_primers == {} and i == 0:
                    left_primer_region = primer_region
            elif left_back_up_primers == {} and i != 0:
                    left_primer_region = primer_region

            primer_region = right_primer_region
            right_primer_region, right_back_up_primers = search_primer_region(  len_target, 
                                                                                    path, 
                                                                                    seq_id, 
                                                                                    seq_temp, 
                                                                                    primer_type = 'RIGHT', 
                                                                                    site_target_temp = site_target_temp, 
                                                                                    left_primer_region=left_primer_region,
                                                                                    right_primer_region=right_primer_region, 
                                                                                    protect_region=80, 
                                                                                    max_region = max_primer_region,
                                                                                    step=step)
            if  right_back_up_primers == {} and i ==0:
                    right_primer_region = primer_region
            elif right_back_up_primers == {} and i !=0:
                    right_primer_region = primer_region


            #2.若左、右引物模板脱靶，延长引物模板，改变测序目标片段长度  
            if left_back_up_primers == {} and right_back_up_primers == {} and i==0: 
                step=40
                # max_primer_region = diff_len - protect_region
                primer_region = left_primer_region
                # left_primer_region = max_primer_region
                left_primer_region, left_back_up_primers = search_primer_region( len_target, 
                                                                                path,
                                                                                seq_id,   
                                                                                seq_temp, 
                                                                                primer_type = 'LEFT', 
                                                                                site_target_temp=site_target_temp, 
                                                                                left_primer_region = left_primer_region,
                                                                                right_primer_region = right_primer_region, 
                                                                                protect_region=80, 
                                                                                max_region= expend_traget_len, 
                                                                                step=step)   
                if left_back_up_primers == {}: 
                    #左引物必定脱靶，尝试右引物是否脱靶
                    left_primer_region= primer_region
                    # right_primer_region = max_primer_region
                    right_primer_region, right_back_up_primers = search_primer_region( len_target, 
                                                                                        path, 
                                                                                        seq_id, 
                                                                                        seq_temp, 
                                                                                        primer_type = 'RIGHT', 
                                                                                        site_target_temp=site_target_temp, 
                                                                                        left_primer_region = left_primer_region,
                                                                                        right_primer_region = right_primer_region, 
                                                                                        protect_region=80, 
                                                                                        max_region = expend_traget_len,
                                                                                        step=step)
                    if right_back_up_primers == {}:
                        #右引物脱靶，且左引物也脱靶  
                        right_primer_region = primer_region
                    else:
                        #右引物不脱靶，而左引物脱靶 
                        pass

                else:  
                    #左引物不脱靶, 确定右引物扩增的长度 
                    primer_region = right_primer_region
                    right_primer_region, right_back_up_primers = search_primer_region( len_target, 
                                                                                        path, 
                                                                                        seq_id, 
                                                                                        seq_temp, 
                                                                                        primer_type = 'RIGHT', 
                                                                                        site_target_temp=site_target_temp, 
                                                                                        left_primer_region = left_primer_region,
                                                                                        right_primer_region = right_primer_region, 
                                                                                        protect_region=80, 
                                                                                        max_region = expend_traget_len,
                                                                                        step=step)
                    if right_back_up_primers=={}:
                        #右引物脱靶脱靶，而左引物不脱靶
                        right_primer_region = primer_region
                    else:
                        #右引物不脱靶，且左引物也不脱靶
                        # max_primer_region = max_primer_region - (right_primer_region - 40)
                        pass

                # 检查是否仍然脱靶：
                if left_back_up_primers == {} and right_back_up_primers == {}:
                    # 引物模板脱靶
                    break
            
            #经过尝试脱靶检测后，仍然脱靶
            if left_back_up_primers == {} and right_back_up_primers == {}:
                break
            
            
            
            #生成引物模板
            start = (site_target_temp - (left_primer_region + protect_region) )
            end = start + (protect_region + left_primer_region) + len_target + (right_primer_region + protect_region)
            pcr_template = seq_temp[ start : end ] 

     
            #设计
            dict_res_p = util.primer_design(seqId=seq_id, 
                                                            seqTemplate=pcr_template,                                           #开始使用
                                                            stype='left_right_change',
                                                            mute_type='single',
                                                            primer_type=seq_type,
                                                            global_args=config.Q_GLOBAL_ARGS  
                                                        )
              
            if ('ok 0' in dict_res_p.get('PRIMER_LEFT_EXPLAIN') ) and  ('ok 0' not in dict_res_p.get('PRIMER_RIGHT_EXPLAIN') ) > 0:
                left_primer_region = left_primer_region + 10
                # right_primer_offtarget_judge_tag = 'no'
            
            elif ( 'ok 0' in dict_res_p.get('PRIMER_RIGHT_EXPLAIN') ) and  ('ok 0' not in dict_res_p.get('PRIMER_LEFT_EXPLAIN') ):
                right_primer_region = right_primer_region + 10  
                # left_primer_offtarget_judge_tag = 'no'
            
            elif ( 'ok 0' in dict_res_p.get('PRIMER_RIGHT_EXPLAIN')  and 'ok 0' in dict_res_p.get('PRIMER_LEFT_EXPLAIN')  ):
                left_primer_region = left_primer_region + 10
                right_primer_region = right_primer_region + 10 

            elif dict_res_p.get('PRIMER_LEFT_NUM_RETURNED') > 0 and dict_res_p.get('PRIMER_RIGHT_NUM_RETURNED') > 0:
                break  
            i = i + 1
        

        #左右引物模板都尝试后，仍然脱靶
        if left_back_up_primers == {} and right_back_up_primers == {}:  
            dict_res_p = {}
            pcr_template = {}  
            
        else:
            #重新确定测序目标片段区域与长度
            if left_primer_region <=120 and right_primer_region <=120:
                #测序目标片段不变
                pass
            elif (left_primer_region > 120 and right_primer_region <=120) or (right_primer_region > 120 and left_primer_region <= 120):
                if (i * 10) <= 120:  
                    if (left_primer_region > 120 and right_primer_region <=120) :
                        len_target = (left_primer_region - step - protect_region) + len_target
                        seq_target = seq_temp[site_target_temp - (left_primer_region - step - protect_region) : site_target_temp - (left_primer_region - step - protect_region) + len_target  ]
                    elif (right_primer_region > 120 and left_primer_region <= 120):

                        len_target = right_primer_region - step - protect_region  + len_target
                        seq_target = seq_temp[site_target_temp: site_target_temp + len_target  ] 
                else:
                    errorMessage = f'120bp的引物区间都没设计出引物出来,说明序列有问题'
                    raise ValueError(errorMessage) 
            elif left_primer_region > 120 and right_primer_region > 120:
                if (i * 10) <= 120:
                    len_target =  (left_primer_region - step - protect_region) + len_target +  (right_primer_region - step - protect_region)
                    seq_target = seq_temp[ site_target_temp - (left_primer_region - step - protect_region) :  site_target_temp - (left_primer_region - step - protect_region) + len_target ]
                else:
                    errorMessage = f'120bp的引物区间都没设计出引物出来,说明序列有问题'
                    raise ValueError(errorMessage)     


        return dict_res_p, seq_target, len_target, pcr_template, left_primer_region, right_primer_region


def parser_pcr_primer(seq_temp, len_target ,pcr_template, sequencing_region, dict_res_p, seq_type, plasmid_seq, path, seq_id, dict_seq_primer, dict_seq_primer_failtrue, gene_change_type,ref_seq,seq_altered):

        #需要确认引物是否涉及成功，且何种类型 
        if len(dict_res_p) > 10 and seq_type == 'genome_seq': 
            #取第一条左引物                             
            u_left = judge_primer_is_or_not( plasmid_seq,
                                    path,
                                    seq_id,
                                    dict_res_p,
                                    dict_seq_primer,
                                    dict_seq_primer_failtrue,
                                    primer_name="SEQUENCING_PRIMER_1",
                                    type='PCR_LEFT')

            #计算出引物模板在基因组上的位置
            gene_start, gene_end = seq_id.split(':')[1].split('-')
            gene_start, gene_end = int(gene_start), int(gene_end)
            template_start_in_seq_temp_position  = seq_temp.find(pcr_template)
            seq_temp_in_genome_position = gene_start - config.UHA_DHA_LENGTH['uha'] - config.PRIMER_TEMPLATE_EXTEND
            template_in_genome_position = seq_temp_in_genome_position + template_start_in_seq_temp_position


            if u_left != -1:
                #左引物不脱靶----左引物设计成功
                primer_name="SEQUENCING_PRIMER_1"
                type='LEFT'  
                dict_seq_primer[primer_name] = dict_res_p[f'PRIMER_{type}_{u_left}_SEQUENCE']
                dict_seq_primer[f'{primer_name}_TM'] = dict_res_p[f'PRIMER_{type}_{u_left}_TM']

                u_right = u_left
                primer_name="SEQUENCING_PRIMER_2"
                type='RIGHT'    
                dict_seq_primer[primer_name] = dict_res_p[f'PRIMER_{type}_{u_right}_SEQUENCE']
                dict_seq_primer[f'{primer_name}_TM'] = dict_res_p[f'PRIMER_{type}_{u_right}_TM']

                #取出此条右引物，且判断是否脱靶    
                right_primer = dict_res_p[f'PRIMER_{type}_{u_right}_SEQUENCE']
                right_primer_position = dict_res_p[f'PRIMER_RIGHT_{u_right}'][0]               #右引物最后1bp在模板上的位置 
                left_primer_position = dict_res_p[f'PRIMER_LEFT_{u_left}'][0]
                product_seq_len = dict_res_p[f'PRIMER_PAIR_{u_right}_PRODUCT_SIZE']            #取出产物的长度   

                if gene_change_type == 'deletion':
                    gene_len = len(ref_seq)
                    right_primer_in_genome_position = template_in_genome_position + gene_len + right_primer_position               #右引物在基因组上最后1bp的位置
                    left_primer_in_genome_position = template_in_genome_position  + left_primer_position
                elif gene_change_type == 'insertion':
                    gene_len = len(seq_altered)
                    right_primer_in_genome_position = template_in_genome_position - gene_len + right_primer_position               #右引物在基因组上最后1bp的位置
                    left_primer_in_genome_position = template_in_genome_position + left_primer_position
                elif gene_change_type == 'substitution':
                    right_primer_in_genome_position = template_in_genome_position  + right_primer_position 
                    left_primer_in_genome_position = template_in_genome_position  + left_primer_position

                primer_type = 'right'

                tag = judge_primer_offtarget(path, seq_id, primer_name, dict_seq_primer_failtrue, right_primer, right_primer_in_genome_position, primer_type, left_primer_in_genome_position)
            
            else:
                #左引物脱靶, 尝试右引物是否全部脱靶 
                u_right = judge_primer_is_or_not( plasmid_seq,  
                                        path,
                                        seq_id,
                                        dict_res_p,
                                        dict_seq_primer,
                                        dict_seq_primer_failtrue,
                                        primer_name="SEQUENCING_PRIMER_2",
                                        type='PCR_RIGHT')
                if u_right !=-1:
                    #左引物脱靶，右引物不脱靶
                    primer_name="SEQUENCING_PRIMER_2"
                    type='RIGHT'
                    dict_seq_primer[primer_name] = dict_res_p[f'PRIMER_{type}_{u_right}_SEQUENCE']
                    dict_seq_primer[f'{primer_name}_TM'] = dict_res_p[f'PRIMER_{type}_{u_right}_TM']

                    u_left=u_right
                    primer_name="SEQUENCING_PRIMER_1"
                    type='LEFT'  
                    dict_seq_primer[primer_name] = dict_res_p[f'PRIMER_{type}_{u_left}_SEQUENCE']
                    dict_seq_primer[f'{primer_name}_TM'] = dict_res_p[f'PRIMER_{type}_{u_left}_TM']
                    dict_seq_primer_failtrue

                    #取出此条右引物，且判断是否脱靶
                    left_primer = dict_res_p[f'PRIMER_{type}_{u_left}_SEQUENCE']
                    #右引物相对于基因组的位置
                    left_primer_position = dict_res_p[f'PRIMER_LEFT_{u_left}'][0]                                                   #左引物在基因组上最前1b的位置
                    right_primer_position = dict_res_p[f'PRIMER_RIGHT_{u_right}'][0] 
                    product_seq_len = dict_res_p[f'PRIMER_PAIR_{u_left}_PRODUCT_SIZE']            #取出产物的长度
                    #计算出左引物在在基因组上的位置  
                    if gene_change_type == 'deletion':
                        gene_len = len(ref_seq)
                        left_primer_in_genome_position = template_in_genome_position  + left_primer_position               #右引物在基因组上最后1bp的位置
                        right_primer_in_genome_position = template_in_genome_position + gene_len + right_primer_position
                    elif gene_change_type == 'insertion':
                        gene_len = len(seq_altered)
                        left_primer_in_genome_position = template_in_genome_position  + left_primer_position               #右引物在基因组上最后1bp的位置
                        right_primer_in_genome_position = template_in_genome_position - gene_len + right_primer_position
                    elif gene_change_type == 'substitution':
                        left_primer_in_genome_position = template_in_genome_position  + left_primer_position 
                        right_primer_in_genome_position = template_in_genome_position + right_primer_position 
                    primer_type = 'left'
                    tag = judge_primer_offtarget(path, seq_id, primer_name, dict_seq_primer_failtrue, left_primer, left_primer_in_genome_position, primer_type, right_primer_in_genome_position)

                else:
                    #引物模板不脱靶，引物脱靶
                    primer_name = 'SEQUENCING_PRIMER_1_AND_SEQUENCING_PRIMER_2'
                    dict_seq_primer_failtrue[f'invalid_primer:{seq_id}:{primer_name}:off_target'] = {}
                    dict_seq_primer['SEQUENCING_TARGET'] = np.nan
                    dict_seq_primer['SEQUENCING_TARGET_LENGTH'] = np.nan

            #取出产物、长度
            if  dict_res_p.get(f'PRIMER_LEFT_{u_left}') != None and dict_res_p.get(f'PRIMER_RIGHT_{u_right}') != None:
                start = dict_res_p[f'PRIMER_LEFT_{u_left}'][0]
                dict_seq_primer['SEQUENCING_TARGET'] = pcr_template[ start : dict_res_p[f'PRIMER_RIGHT_{u_right}'][0] + 1 ]
                dict_seq_primer['SEQUENCING_TARGET_LENGTH'] = dict_res_p[f'PRIMER_RIGHT_{u_right}'][0] - start + 1
            else:
                dict_seq_primer['SEQUENCING_TARGET'] = np.nan
                dict_seq_primer['SEQUENCING_TARGET_LENGTH'] = np.nan
    
        elif  len(dict_res_p) > 10 and seq_type == 'plasmid_seq' and len_target > sequencing_region: 
            u = 0
            primer_name= 'SEQUENCING_PRIMER_1'
            dict_seq_primer[primer_name] = dict_res_p[f'PRIMER_LEFT_0_SEQUENCE']
            dict_seq_primer[f'{primer_name}_TM'] = dict_res_p[f'PRIMER_LEFT_0_TM']
            primer_name = 'SEQUENCING_PRIMER_2'
            dict_seq_primer[primer_name] = dict_res_p[f'PRIMER_RIGHT_0_SEQUENCE']
            dict_seq_primer[f'{primer_name}_TM'] = dict_res_p[f'PRIMER_RIGHT_0_TM']

        elif len(dict_res_p) > 10 and seq_type == 'plasmid_seq' and len_target <= sequencing_region: 
            u = 0
            primer_name = 'SEQUENCING_PRIMER_1'  
            dict_seq_primer[primer_name] = dict_res_p[f'PRIMER_LEFT_0_SEQUENCE']
            dict_seq_primer[f'{primer_name}_TM'] = dict_res_p[f'PRIMER_LEFT_0_TM']

        elif dict_res_p == {}:
            #引物模板脱靶，导致引物脱靶
            primer_name = 'TEMPLATE_AND_SEQUENCING_PRIMER_1_AND_SEQUENCING_PRIMER_2'
            dict_seq_primer_failtrue[f'invalid_primer:{seq_id}:{primer_name}:off_target'] = {}
            dict_seq_primer['SEQUENCING_TARGET'] = np.nan  
            dict_seq_primer['SEQUENCING_TARGET_LENGTH'] = np.nan 

        return dict_seq_primer, dict_seq_primer_failtrue, primer_name

#设计测序引物
def design_seq_primer(gene_change_type, ref_seq, seq_altered, plasmid_seq, path, seq_id, seq_temp, seq_target, mute_position=0,seq_type=''):

    if seq_type == 'genome_seq':
        if gene_change_type == 'deletion':
            wt_gene_target = ref_seq
            wt_gene_target_len = len(ref_seq)
        elif gene_change_type == 'insertion':
            wt_gene_target = seq_altered
            wt_gene_target_len = len(seq_altered)
        elif gene_change_type == 'substitution':
            wt_gene_target = seq_altered
            wt_gene_target_len = len(seq_altered)    

    wt_gene_target = ''  
    wt_gene_target_len = ''   

    dict_seq_primer={}
    dict_seq_primer_failtrue={}
    len_target=len(seq_target)     
    site_target_temp=seq_temp.find(seq_target)
    
    left_primer_region = 40   #默认40-120---引物模板（600+80-----600+200），40-920---引物模板（600+80------- 680 + 120 + 800 ）
    right_primer_region = 40  #默认40-120
    protect_region = 80        
    sequencing_region = 600
    max_primer_region = 120   #120
    # max_region = 120   
   
    #首先计算总共需要多少引物
    primer_nums =int(len_target / sequencing_region) + 1 

    #扩增之后测序目标的长度
    if len_target < sequencing_region:  
        expend_traget_len = (primer_nums + 2) * sequencing_region + 2 * protect_region 
    else:
        expend_traget_len = (primer_nums + 1) * sequencing_region + 2 * protect_region 
    
    #差值
    diff_len = expend_traget_len - len_target


    analysis = False
    #
    if analysis == False:  
        #  
        if seq_type == 'genome_seq': 

            while len_target + 2 * protect_region <= expend_traget_len: 
                origin_seq_target =  seq_target
                origin_len_target =  len_target 
                dict_res_p, seq_target, len_target, pcr_template, left_primer_region, right_primer_region = try_pcr_offtarget_design( 
                                                                                                                path,
                                                                                                                len_target,
                                                                                                                seq_target,
                                                                                                                expend_traget_len,
                                                                                                                left_primer_region, 
                                                                                                                right_primer_region, 
                                                                                                                max_primer_region, 
                                                                                                                protect_region, 
                                                                                                                site_target_temp, 
                                                                                                                seq_id, 
                                                                                                                seq_type, 
                                                                                                                seq_temp)
                dict_seq_primer, dict_seq_primer_failtrue, primer_name =  parser_pcr_primer(seq_temp,
                                                                                            len_target,
                                                                                            pcr_template,
                                                                                            sequencing_region, 
                                                                                            dict_res_p, 
                                                                                            seq_type, 
                                                                                            plasmid_seq, 
                                                                                            path, 
                                                                                            seq_id, 
                                                                                            dict_seq_primer, 
                                                                                            dict_seq_primer_failtrue, 
                                                                                            gene_change_type,
                                                                                            ref_seq,
                                                                                            seq_altered)
                if primer_name =='SEQUENCING_PRIMER_1_AND_SEQUENCING_PRIMER_2': 
                    if origin_len_target == len_target: 
                        left_primer_region = left_primer_region + 10
                        right_primer_region = right_primer_region + 10  
                    else:
                        left_primer_region = 40 + 10
                        right_primer_region = 40 + 10
                    continue
                else:  
                    break
                
        if seq_type == 'plasmid_seq':  
            # 尝试在最大的引物设计区域中设计出引物 600
            dict_res_p, pcr_template = try_pcr_primer_design(len_target , left_primer_region, right_primer_region, max_primer_region, protect_region, site_target_temp, seq_id, seq_type, seq_temp)
            dict_seq_primer, dict_seq_primer_failtrue, primer_name =  parser_pcr_primer(seq_temp,
                                                                                            len_target,
                                                                                            pcr_template,
                                                                                            sequencing_region, 
                                                                                            dict_res_p, 
                                                                                            seq_type, 
                                                                                            plasmid_seq, 
                                                                                            path, 
                                                                                            seq_id, 
                                                                                            dict_seq_primer, 
                                                                                            dict_seq_primer_failtrue, 
                                                                                            gene_change_type,
                                                                                            ref_seq,
                                                                                            seq_altered)  
        

    else:
        start = (site_target_temp - (left_primer_region + protect_region) )
        end = start + (left_primer_region + protect_region ) + len_target + (right_primer_region + protect_region)
        pcr_template = seq_temp[ start : end ]

         #2.设计
        dict_res_p = util.primer_design(seqId=seq_id, 
                                            seqTemplate=pcr_template,                                           
                                            stype='left_right_change',
                                            mute_type='single',
                                            primer_type=seq_type,
                                            global_args=config.Q_GLOBAL_ARGS  
                                        )
        if dict_res_p.get('PRIMER_LEFT_NUM_RETURNED') == 0:
            #设计失败
            primer_name = 'SEQUENCING_PRIMER_1_AND_SEQUENCING_PRIMER_2:failture'
            dict_seq_primer_failtrue[f'{seq_id}:{primer_name}'] = dict_res_p 
        else:  
            #检查第一对PCR引物设计情况  
            dict_seq_primer, dict_seq_primer_failtrue, primer_name =  parser_pcr_primer(seq_temp,
                                                                                            len_target,
                                                                                            pcr_template,
                                                                                            sequencing_region, 
                                                                                            dict_res_p, 
                                                                                            seq_type, 
                                                                                            plasmid_seq, 
                                                                                            path, 
                                                                                            seq_id, 
                                                                                            dict_seq_primer, 
                                                                                            dict_seq_primer_failtrue, 
                                                                                            gene_change_type,
                                                                                            ref_seq,
                                                                                            seq_altered)

    
    #继续测序引物设计
    if 1200 < len_target:    
        i=3
        #更新测序目标的长度,去除第一条PCR引物可测的sequencing_region区域，得到更新后的seq_target，用于测序引物设计              
        seq_target = seq_target[sequencing_region:]
        site_target_temp=seq_temp.find(seq_target)      
        len_target=len(seq_target)

        if len_target >= sequencing_region:
            while len_target >= sequencing_region:
                
                #抽取seq_target的前200bp
                pre_seq = seq_temp[seq_temp.find(seq_target) - protect_region - max_primer_region : seq_temp.find(seq_target)]
                sequencing_template = pre_seq + seq_target
                left_primer_region = 40

                dict_res_pn = try_sequencing_primer_design( seq_id, sequencing_template, left_primer_region, max_primer_region, protect_region, seq_type)

                # dict_res_pn = try_sequencing_primer_design( seq_id, seq_temp, site_target_temp, left_primer_region, protect_region, seq_type, max_primer_region)
                #判断引物是否成功
                primer_name = f"SEQUENCING_PRIMER_{i}"
                dict_seq_primer[primer_name] = dict_res_pn[f'PRIMER_LEFT_0_SEQUENCE']  
                dict_seq_primer[f'{primer_name}_TM'] = dict_res_pn[f'PRIMER_LEFT_0_TM']
                #更新测序区域序列
                seq_target = seq_target[sequencing_region:]
                len_target = len(seq_target)

                i = i + 1
                
        elif 200 <= len_target < sequencing_region:
            
            pre_seq = seq_temp[seq_temp.find(seq_target) - protect_region - max_primer_region : seq_temp.find(seq_target)]
            sequencing_template = pre_seq + seq_target
            left_primer_region = 40

            dict_res_pn = try_sequencing_primer_design( seq_id, sequencing_template, left_primer_region, max_primer_region, protect_region, seq_type)
            primer_name = f"SEQUENCING_PRIMER_{i}"
            dict_seq_primer[primer_name] = dict_res_pn[f'PRIMER_LEFT_0_SEQUENCE']
            dict_seq_primer[f'{primer_name}_TM'] = dict_res_pn[f'PRIMER_LEFT_0_TM']

    return dict_seq_primer, dict_seq_primer_failtrue  
   

def judge_primer_offtarget(path, seq_id, primer_name, primer_fail, primer, primer_in_genome_position, primer_type, another_primer_in_genome_position):

    tag = -1

    if path.split('.')[-1] != 'gb':
        chrom = seq_id.split(';')[1].split(':')[0]
        seq = util.read_chromSeq_from_genome(path, chrom)
        index_prefix_path = os.path.join(config.tmp_path, 'genome_index')
    elif path.split('.')[-1] == 'gb':
        #gb文件转换成fasta文件
        gb_path = path
        fasta_path = splitext(gb_path)[0]+'.fasta'
        util.gb_2_fna(gb_path, fasta_path)
        path = fasta_path
        index_prefix_path = os.path.join(config.tmp_path, 'gb_index')

    fasta_filename = os.path.join(config.tmp_path,'sequence.fasta')
    data={
                'id':[seq_id],  
                'sequence':[primer]
            }
    df = pd.DataFrame(data)
    util.df_to_fasta(df,fasta_filename)  
    sam_df = util.bowtie_seq_genome(config.tmp_path, path, index_prefix_path, fasta_filename)


    if primer_type == 'right':
        primer_start_in_genome_position = ( primer_in_genome_position - len(primer) + 1) 
    elif primer_type == 'left':
        primer_start_in_genome_position = primer_in_genome_position



    if len(sam_df) >1:
        #脱靶
        off_target_df = sam_df
        off_target_dict = {}
        for i,v in off_target_df.iterrows():
            
            ReadName = v['ReadName']
            Sequence = v['Sequence']
            ReferenceStart = v['ReferenceStart']
            Mismatch = v['Mismatch']
            MatchingNumber = v['MatchingNumber']
            
            if (ReferenceStart - 1) == primer_start_in_genome_position :
                primer_seq = ';'.join([primer,str(ReferenceStart)])  
            else:
                primer_seq = ';'.join([Sequence,str(ReferenceStart),Mismatch,MatchingNumber])
            off_target_dict.update({f'{primer_name}_{i}':primer_seq})

        #判断此对引物为：有效引物对，无效引物对
        if primer_type == 'right':
            temp = off_target_df[(off_target_df['chain']=='16') & (off_target_df['ReferenceStart'] != primer_start_in_genome_position + 1) ]                          #负义链比对结果
            if len(temp) > 0:
                temp['distance'] = (temp['ReferenceStart']-1) - another_primer_in_genome_position
                temp = temp[temp['distance'] > 0]
                len1 = len(temp)
                temp = temp[ temp['distance'] >  (primer_in_genome_position - another_primer_in_genome_position) ]
                len2 = len(temp) 
                if len1 == len2:
                    #有效引物
                    id = f'valid_primer:{seq_id}:{primer_name}:off_target'
                else:
                    #无效引物
                    id = f'invalid_primer:{seq_id}:{primer_name}:off_target'
            else:
                #有效引物
                id = f'valid_primer:{seq_id}:{primer_name}:off_target'
        elif primer_type == 'left':
            temp = off_target_df[(off_target_df['chain']=='0') & (off_target_df['ReferenceStart'] != primer_start_in_genome_position + 1) ]                          #负义链比对结果
            if len(temp) > 0:
                temp['distance'] = another_primer_in_genome_position - (temp['ReferenceStart']-1)
                temp = temp[temp['distance'] > 0]
                len1 = len(temp)
                temp = temp[ temp['distance'] >  (another_primer_in_genome_position -  primer_in_genome_position) ]
                len2 = len(temp) 
                if len1 == len2:
                    #有效引物
                    id = f'valid_primer:{seq_id}:{primer_name}:off_target'
                else:
                    #无效引物
                    id = f'invalid_primer:{seq_id}:{primer_name}:off_target'
            else:
                #有效引物
                id = f'valid_primer:{seq_id}:{primer_name}:off_target'

        primer_fail[id] = off_target_dict
    else:
        tag = 0
    
    return tag

#判断引物是与否 
def judge_primer_is_or_not(plasmid_seq, path, seq_id, dict_res, primer_suc, primer_fail, primer_name, type='LEFT'):


    judge_type = ''
    #检索引物的个数：
    u = -1
    if path.split('.')[-1] != 'gb':
        chrom = seq_id.split(';')[1].split(':')[0]
        seq = util.read_chromSeq_from_genome(path, chrom)
        index_prefix_path = os.path.join(config.tmp_path, 'genome_index')
    elif path.split('.')[-1] == 'gb':
        #gb文件转换成fasta文件
        gb_path = path
        fasta_path = splitext(gb_path)[0]+'.fasta'
        util.gb_2_fna(gb_path, fasta_path)
        path = fasta_path
        index_prefix_path = os.path.join(config.tmp_path, 'gb_index')
    else: 
        seq = plasmid_seq

    fasta_file = os.path.join(config.tmp_path,'sequence.fasta')
    if type =='PCR_LEFT':
        type = 'LEFT'
        if dict_res.get('PRIMER_LEFT_NUM_RETURNED') == None or dict_res.get('PRIMER_LEFT_NUM_RETURNED') == 0:
            primer_fail[f'{seq_id}:{primer_name}'] = dict_res
        else:
            primer_num = dict_res.get('PRIMER_LEFT_NUM_RETURNED')
            left_primer_df = util.parse_primer3_result(dict_res, primer_num)
            left_primer_df = left_primer_df[['id', 'left_primer']].rename(columns={'left_primer': 'sequence'})
            left_primer_df = left_primer_df.dropna(how='all')
            left_primer_df['id']=left_primer_df['id'].astype('str')
            left_primer_df['id']=left_primer_df['id']+'_left'
            primer_df = left_primer_df
    elif type =='PCR_RIGHT':  
        type ='RIGHT'
        if dict_res.get('PRIMER_RIGHT_NUM_RETURNED') == None or dict_res.get('PRIMER_RIGHT_NUM_RETURNED') == 0:
            primer_fail[f'{seq_id}:{primer_name}'] = dict_res
        else:
            primer_num = dict_res.get('PRIMER_RIGHT_NUM_RETURNED')  
            right_primer_df = util.parse_primer3_result(dict_res, primer_num) 
            right_primer_df = right_primer_df[['id', 'right_primer']].rename(columns={'right_primer': 'sequence'}) 
            right_primer_df = right_primer_df.dropna(how='all')
            right_primer_df['id'] = right_primer_df['id'].astype('str')
            right_primer_df['id'] = right_primer_df['id']+'_right' 
            primer_df = right_primer_df
    elif type == 'LEFT' or type =='RIGHT':
        if dict_res.get(f'PRIMER_{type}_NUM_RETURNED') == None or dict_res.get(f'PRIMER_{type}_NUM_RETURNED') == 0:
            primer_fail[f'{seq_id}:{primer_name}'] = dict_res
        else:
            u = 0
            primer_suc[primer_name] = dict_res.get(f'PRIMER_{type}_{u}_SEQUENCE')
            primer_suc[f'{primer_name}_TM'] = dict_res.get(f'PRIMER_{type}_{u}_TM')
        return u

       
    #引物设计成功： 
    if ( dict_res.get('PRIMER_LEFT_NUM_RETURNED') != 0) or ( dict_res.get('PRIMER_RIGHT_NUM_RETURNED') != 0 ) : 
        primer_df = primer_df.drop_duplicates('sequence')
        util.df_to_fasta(primer_df,fasta_file)
        #对引物脱靶检测
        sam_df = util.bowtie_seq_genome(config.tmp_path, path, index_prefix_path, fasta_file)  
        #XM:i:1  NM:i:0  
        no_off_target = sam_df[(sam_df['MatchingNumber'] == 'XM:i:1') & (sam_df['Mismatch'] == 'NM:i:0') ]
        if len(no_off_target) > 0:
                primer = no_off_target[:1]
                primer_id = list(primer.loc[:,'ReadName'])[0] 
                primer_id, primer_type = primer_id.split('_')
                u = primer_id
                primer_type = primer_type.upper()
                primer_suc[primer_name] = dict_res.get(f'PRIMER_{primer_type}_{primer_id}_SEQUENCE')
                primer_suc[f'{primer_name}_TM'] = dict_res.get(f'PRIMER_{primer_type}_{primer_id}_TM')

    return u

#有点问题   
def blast_primer_in_genome(dict_res, type, path, ):

    df = pd.DataFrame()
    u = 0 
    for i,k in dict_res.items():
        id = f'PRIMER_{type}_{u}_SEQUENCE'
        seq=dict_res.get(f'PRIMER_{type}_{u}_SEQUENCE')
        if seq != None:
            df = df.append(pd.DataFrame([{'id':id,'seq':seq}]))
            u = u +1

    primer_workdir = config.workdir+'/primer'
    if not exists(primer_workdir):
        os.makedirs(primer_workdir)
        
    blast_output_file_path= os.path.join(primer_workdir,'blast_output.txt')
    blast_input_file_path = os.path.join(primer_workdir,'blast.fasta')
        
    #create fasta
    util.convert_df_to_fastaFile(df,id='id',name='seq',input_fasta = blast_input_file_path)

    #run blast
    if path.split('.')[-1] == 'gb':
        fna_file = config.workdir+'/primer/xxx.fna'
            #convert fasta
        util.gb_2_fna(path, fna_file)
        genome = fna_file
    else:
        genome = path

    util.blast_ha(genome, blast_input_file_path, blast_output_file_path)

    result_df = util.blast_output_evaluate(workdir = config.workdir + '/primer/' , blast_input = blast_input_file_path, blast_output= blast_output_file_path)

    return result_df 






















    # elif len(dict_res_p) < 10:

    #     primer_fail[f'{seq_id}:{primer_name}:off_target'] = off_target_dict



        # #左脱, 右不脱, 基因组测序引物
    # if seq_type == 'genome_seq':
    #     if left_back_up_primers == {} and right_back_up_primers !={}:
    #         u_left = 0
    #         primer_name="SEQUENCING_PRIMER_1"
    #         type='LEFT'  
    #         dict_seq_primer[primer_name] = dict_res_p[f'PRIMER_{type}_{u_left}_SEQUENCE']
    #         dict_seq_primer[f'{primer_name}_TM'] = dict_res_p[f'PRIMER_{type}_{u_left}_TM']

    #         primer_name="SEQUENCING_PRIMER_2"
    #         type='RIGHT'
    #         dict_seq_primer[primer_name] = dict_res_p[f'PRIMER_{type}_{u_right}_SEQUENCE']
    #         dict_seq_primer[f'{primer_name}_TM'] = dict_res_p[f'PRIMER_{type}_{u_right}_TM']



    
    # if seq_type == 'genome_seq' or (seq_type == 'plasmid_seq' and len_target > sequencing_region  ):
    #     #取出第二条引物
    #     u_right = judge_primer_is_or_not( plasmid_seq,  
    #                             path,
    #                             seq_id,
    #                             dict_res_p,
    #                             dict_seq_primer,
    #                             dict_seq_primer_failtrue,
    #                             primer_name="SEQUENCING_PRIMER_2",
    #                             type='PCR_RIGHT')  
                    
    #     #右脱，左不脱, 基因组测序引物    
    #     if  seq_type == 'genome_seq':
    #         if right_back_up_primers == {} and  left_back_up_primers !={}:
    #             u_right = 0                                                #有问题，需要匹配脱靶序列
    #             primer_name="SEQUENCING_PRIMER_2"
    #             type='RIGHT'
    #             dict_seq_primer[primer_name] = dict_res_p[f'PRIMER_{type}_{u_right}_SEQUENCE']
    #             dict_seq_primer[f'{primer_name}_TM'] = dict_res_p[f'PRIMER_{type}_{u_right}_TM']

    #     #取出产物、长度
    #     if  dict_res_p.get(f'PRIMER_LEFT_{u_left}') != None and dict_res_p.get(f'PRIMER_RIGHT_{u_right}') != None:
    #         start = dict_res_p[f'PRIMER_LEFT_{u_left}'][0]
    #         dict_seq_primer['SEQUENCING_TARGET'] = pcr_template[ start : dict_res_p[f'PRIMER_RIGHT_{u_right}'][0] + 1 ]
    #         dict_seq_primer['SEQUENCING_TARGET_LENGTH'] = dict_res_p[f'PRIMER_RIGHT_{u_right}'][0] - start + 1
    #     else:
    #         dict_seq_primer['SEQUENCING_TARGET'] = np.nan
    #         dict_seq_primer['SEQUENCING_TARGET_LENGTH'] = np.nan
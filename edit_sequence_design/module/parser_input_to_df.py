import pandas as pd
from Bio import SeqIO 
import os,sys
import warnings   
warnings.filterwarnings('ignore')
import configparser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sgRNA_utils.sgRNA_primer_util as su

def conf_read(filename): 
    config = configparser.ConfigParser()
    config.read(filename)
    res = dict(config._sections["point"])
    return res

def input_to_primer_template(input_file_path, genome,config,workdir):
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
                if len(data) < 5:
                    error_message = "Some necessary input information is missing in the line "+mun_id+" of the input file"
                    print(error_message)
                    return  error_message
                else:
                    upstream = data[1].strip().upper()
                    ref = data[2]
                    alt = data[3]
                    mutation_type = data[4].strip().lower()
                    name = mun_id
                    
                    if mun_id not in blast_search_dict:
                        # 没有blast上
                        error_message = "The upstream sequence of " + mun_id + " can not be mapped to the target genome. please check whether the target sequence is located on the target genome."
                        blast_error_handler.write(f'{mun_id}\t{error_message}\n')
                    elif blast_search_dict[mun_id]["unique_mapped"] > 1:
                        # 多次比对上 100
                        error_message = "The upstream sequence of " + mun_id + "  can be mapped to multiple loci in the target genome, %s, Please provide a longer upstream seqeunce." % blast_search_dict[mun_id]["description"]
                        blast_error_handler.write(f'{mun_id}\t{error_message}\n')
                    elif blast_search_dict[mun_id]["unique_mapped"] == 0:
                        # 无 100 比对
                        error_message = "The upstream sequence of " + mun_id + " can not be uniquely mapped to the target genome. Please check whether the target sequence is located on the target genome."
                        blast_error_handler.write(f'{mun_id}\t{error_message}\n')
                    elif blast_search_dict[mun_id]["unique_mapped"] == 1:
                        # 开始突变的碱基在genome上的索引
                        if blast_search_dict[mun_id]["reverse_mapped"]:
                            record = su.revComp(str(record_dict[blast_search_dict[mun_id]["chrom"]].seq))
                            upstream_start_index = len(record) - int(blast_search_dict[mun_id]["start"])
                            strand = "-"  
                        else:
                            record = str(record_dict[blast_search_dict[mun_id]["chrom"]].seq)
                            upstream_start_index = int(blast_search_dict[mun_id]["start"])-1
                            strand = "+"
                        chrom = blast_search_dict[mun_id]["chrom"]
                        mutation_pos_index = upstream_start_index + len(upstream)
                        # get mutation info dict
                        res = create_mutation_info(
                            record,mutation_type,mutation_pos_index,
                            ref,alt,strand,chrom,
                            name,mun_id,config
                            )
                        if isinstance(res,str):
                            blast_error_handler.write(f'{mun_id}\t{res}\n')
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
                        record = su.revComp(str(record_dict[chrom].seq))
                        mutation_pos_index = len(record) - int(pos)
                    # get mutation info dict
                    res = create_mutation_info(
                        record,mutation_type,mutation_pos_index,
                        ref,alt,strand,chrom,
                        name,mun_id,config
                        )
                    if isinstance(res,str):
                        blast_error_handler.write(f'{mun_id}\t{res}\n')
                    else:
                        primer_template[mun_id] = res
        else:
            error_message = "The input file format not supported, Please rightly prepare input file for target manipulation as the example of 2,3-BD."
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
                "unique_mapped":1,  次数 100%比对
                "mapped":1,  次数
                "reverse_mapped":1, 次数
                "description":"NC_001133.9:26529-26708;NC_001133.9:26664-26843;",  
            },
        }

    """
    blast_output_file_path=os.path.join(workdir,'blast_search.txt')
    ref_lib = genome.split('/')[-1].split('.')[0]
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

def create_mutation_info(record,mutation_type,mutation_pos_index,ref,alt,strand,chrom,name,mun_id,config):
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
    max_left_arm_seq_length=int(config['max_left_arm_seq_length'])
    max_right_arm_seq_length=int(config['max_right_arm_seq_length'])

    # length 
    if mutation_pos_index - max_left_arm_seq_length < 0:
        error_message = "The length of upstream sequence of manipulation site of " + mun_id + " must be larger than sum of 'Max Length of UHA' and 'Max Length of UIS'."
        return error_message
    
    if mutation_type ==  "insertion":
        info_dict = {
            "seq_uha_max_whole":str(record[
                mutation_pos_index - max_left_arm_seq_length : mutation_pos_index
            ]),
            "seq_dha_max_whole":str(record[
                mutation_pos_index : mutation_pos_index + max_right_arm_seq_length
            ]),
            "seq_altered":alt,
            "type":mutation_type,
            "ref":ref,
            "strand":"plus" if strand =="+" else "minus",
            "mutation_pos_index":mutation_pos_index,
            "geneid":chrom,
            "name":name,
            "region":f'{chrom}:{mutation_pos_index}-{mutation_pos_index+len(ref)}',
            "uha_upstream":str(  
                record[
                    mutation_pos_index - max_left_arm_seq_length - 100 : mutation_pos_index - max_left_arm_seq_length
                ]
            ),
            "dha_downstream":str(
                record[
                    mutation_pos_index + max_right_arm_seq_length  : mutation_pos_index + max_right_arm_seq_length  + 100
                ]
            ),

        }
        return info_dict      

#         return info_dict
    if mutation_type in ["deletion","substitution"]:
        genome_ref = record[mutation_pos_index:mutation_pos_index+len(ref)]
        # print(genome_ref)
        if genome_ref.upper() == ref.upper():
            info_dict = {
                "seq_uha_max_whole":str(record[
                    mutation_pos_index 
                    - max_left_arm_seq_length :mutation_pos_index
                    ]),
                "seq_dha_max_whole":str(record[
                    mutation_pos_index + len(ref)
                    : mutation_pos_index + len(ref) + max_right_arm_seq_length
                    ]),
                "seq_altered":"" if alt=='-' else alt,
                "type":mutation_type,
                "ref":ref,
                "strand":"plus" if strand == "+" else "minus",
                "mutation_pos_index":mutation_pos_index,
                "geneid":chrom,
                "name":name,   
                "region":f'{chrom}:{mutation_pos_index}-{mutation_pos_index+len(ref)}',
                "uha_upstream":str(  
                    record[
                        mutation_pos_index 
                    - max_left_arm_seq_length - 100
                    : mutation_pos_index 
                    - max_left_arm_seq_length

                    ]),
                "dha_downstream":str(
                    record[
                        mutation_pos_index + len(ref) + max_right_arm_seq_length
                        : mutation_pos_index + len(ref) + max_right_arm_seq_length + 100
                    ]),
            }
            return info_dict
        else:
            error_message = "The target mutation ref of " + mun_id + " can not be found in reference, please check."
            return error_message
    else:
        error_message = "The target manipulation type of " + mun_id + " must be equal to 'insertion,substitution or deletion', Please rightly prepare input file for target manipulation as the example of 2,3-BD."
        return  error_message

def dict_to_df(dict_input_seq):
    info_input_df = pd.DataFrame()
    for item in dict_input_seq:
        df = pd.DataFrame([dict_input_seq[item]])
        df.insert(loc=0,value=item,column='id')
        info_input_df = info_input_df.append(df)
    info_input_df = info_input_df.reset_index(drop=True)
    return info_input_df

def main():
    input_file_path = '/home/yanghe/githubcode/crispr_hr_editor/edit_sequence_design/data/input/editor_info.csv'
    ref_genome = '/home/yanghe/githubcode/crispr_hr_editor/edit_sequence_design/data/temp/GCA_000011325.1_ASM1132v1_genomic.fna'
    conf = '/home/yanghe/githubcode/crispr_hr_editor/edit_sequence_design/data/input/config.txt'
    workdir = '/home/yanghe/githubcode/crispr_hr_editor/edit_sequence_design/data/output/'
    chopchop_input = '/home/yanghe/githubcode/crispr_hr_editor/chopchop/data/input/genome/'
    config = conf_read(conf)  
    dict_input_seq = input_to_primer_template(input_file_path,ref_genome,config,workdir)
    info_input_df = dict_to_df(dict_input_seq)
    info_input_df.to_csv(chopchop_input+'info_input.csv')

if __name__ == 'main':
    main()  
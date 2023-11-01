# editing_sequence_design

## Project Introduction

**editing_sequence_design** serves as the primer design module of AutoESDCas, primarily focusing on automating and optimizing the primer design process for genome editing. It leverages Crispr/HR technology and offers various features, including homologous arm design, PCR primer design, sequencing primer design, and the visualization of recombinant plasmid maps. This module is designed to support multiple experimental conditions, catering to the diverse needs of biological researchers. The supported scenarios include:

1. **Same Plasmid Template - Synthesis Recombinant Plasmid:**
   - The sgRNA and homologous arm exist on the same plasmid template.
   - Recombinant plasmid synthesis is employed.

2. **Different Plasmid Templates - Synthesis Recombinant Plasmid:**
   - The sgRNA and homologous arms are on separate plasmid templates.
   - Recombinant plasmid synthesis is used.

3. **Same Plasmid Template - PCR-Based Recombinant Plasmid:**
   - The sgRNA and homologous arm are on the same plasmid template.
   - The sgRNA sequence fragment, homologous arm sequence fragment, and recombinant plasmid sequence fragment are obtained via PCR.
   - Primer design for the recombinant plasmid sequence fragment supports three scenarios:
     - i. System automatic design
     - ii. Manual provision of partial primers
     - iii. Manual provision of primer design range.

4. **Different Plasmid Templates - PCR-Based Recombinant Plasmid:**
   - The sgRNA and homologous arms are located on different plasmid templates.
   - The sgRNA sequence fragment, homologous arm sequence fragment, and two different recombinant plasmid sequence fragments are obtained through PCR.
   - Primer design for the recombinant plasmid sequence fragment supports three scenarios:
     - i. System automatic design
     - ii. Manual provision of partial primers
     - iii. Manual provision of primer design range.

5. **Different Plasmid Templates - PCR-Based Homologous Arm Sequence Fragment:**
   - The sgRNA and homologous arm are located in different plasmid templates.
   - The homologous arm sequence fragment and the recombinant plasmid sequence fragment where the homologous arm is located are obtained through PCR.
   - Primer design for the recombinant plasmid sequence fragment where the source arm is located supports three primer design scenarios:
     - i. System automatic design
     - ii. Manual provision of partial primers
     - iii. Manual provision of primer design range.
   - Additionally, sgRNA sequence fragments are obtained through primer annealing, and the recombinant plasmid fragments where the sgRNA is located are obtained through enzymatic cleavage.



## Installation

It is recommended to use Python 3.8 for **Editing Sequence Design**. To install the required dependencies, execute the following command:

```shell
pip install -r requirements.txt
```

## Usage & Example
       


### 1. Same Plasmid Template - Synthesis Recombinant Plasmid

**Input:**
- **Step 1:** Upload the genome (fna) file, data preprocessed (CSV) file, sgRNA (CSV) file, plasmid (gb) file.
- **Step 2:** Configure relevant parameters for primer desgin

    - Example configuration (data):

     ```json
        {   
            "chopchop_input": "/home/XXX/tmp/data_preprocessing/output/info_input.csv",   
            "sgRNA_result_path": "/home/XXX/tmp/chopchop/output/sgRNA.csv",
            "edit_sequence_design_workdir":"/home/XXX/tmp/edit_sequence_design/output", 
            "ref_genome": "/home/XXX/tmp/edit_sequence_design/GCA_000011325.1_ASM1132v1_genomic.fna",
            "one_plasmid_file_path":"./input/pMB1-sgRNA-wacJ1.gb",     
            "bowtie_path":"/home/XXXS/software/bowtie",
            "scene":"both_sgRNA_primer",
            "plasmid_metod":'1',

            "plasmid_label":{
                "ccdb_label":"HR arm",  
                "promoter_terminator_label":"gRNA ORF",
                "n_20_label":"wacJ",
                "promoter_label":"promoter" 
            },

            "uha_dha_config": {     
                "max_right_arm_seq_length": 145,     
                "max_left_arm_seq_length": 145,   
                "min_left_arm_seq_length": 145,        
                "min_right_arm_seq_length": 145     
            },

            "PLASMID_Q_ARGS":{
                "PRIMER_OPT_TM": 65,
                "PRIMER_MIN_TM": 55,  
                "PRIMER_MAX_TM": 75,    
                "PRIMER_MIN_GC": 20,
                "PRIMER_OPT_GC":65,
                "PRIMER_MAX_GC": 80,
                "PRIMER_MIN_SIZE":15,
                "PRIMER_MAX_SIZE":25,
                "PRIMER_OPT_SIZE":18
            },

            "GENOME_Q_ARGS":{
                "PRIMER_OPT_TM": 65,
                "PRIMER_MIN_TM": 55,     
                "PRIMER_MAX_TM": 75,    
                "PRIMER_MIN_GC": 20,
                "PRIMER_OPT_GC":65,
                "PRIMER_MAX_GC": 80,
                "PRIMER_MIN_SIZE":15,
                "PRIMER_MAX_SIZE":25,
                "PRIMER_OPT_SIZE":18
            }     
        }
     ```

**Execute:**
```shell
python edit_main.py
```

**Output:**
- `one_plasmid_system_result.zip`

These files will be generated in the `/home/XXX/tmp/edit_sequence_design/output/` directory.



### 2. Different Plasmid Templates - Synthesis Recombinant Plasmid:
**Input:**
- **Step 1:** Upload the genome (fna) file, data preprocessed (CSV) file, sgRNA (CSV) file, plasmid (gb) file.
- **Step 2:** Configure relevant parameters for primer desgin

    - Example configuration (data):

     ```json
        {   
            "chopchop_input": "/home/XXX/tmp/data_preprocessing/output/info_input.csv",   
            "sgRNA_result_path": "/home/XXX/tmp/chopchop/output/sgRNA.csv",
            "edit_sequence_design_workdir":"/home/XXX/tmp/edit_sequence_design/output", 
            "ref_genome": "/home/XXX/tmp/edit_sequence_design/GCA_000011325.1_ASM1132v1_genomic.fna",
            "no_ccdb_plasmid":"/home/XXX/program/edit_sequence_design/input/no-ccdb-pXMJ19-Cas9A-gRNA-crtYEb-Ts - ori.gb",
            "no_sgRNA_plasmid":"/home/XXX/program/edit_sequence_design/input/no-sgRNA-pXMJ19-Cas9A-gRNA-crtYEb-Ts - ori.gb",
            "bowtie_path":"/home/XXX/software/bowtie",
            "scene":"both_sgRNA_primer",
            "plasmid_metod":'1',

            "plasmid_label":{
                "ccdb_label":"HR arm",  
                "promoter_terminator_label":"gRNA ORF",
                "n_20_label":"wacJ",
                "promoter_label":"promoter" 
            },

            "uha_dha_config": {     
                "max_right_arm_seq_length": 145,     
                "max_left_arm_seq_length": 145,   
                "min_left_arm_seq_length": 145,        
                "min_right_arm_seq_length": 145     
            },

            "PLASMID_Q_ARGS":{
                "PRIMER_OPT_TM": 65,
                "PRIMER_MIN_TM": 55,  
                "PRIMER_MAX_TM": 75,    
                "PRIMER_MIN_GC": 20,
                "PRIMER_OPT_GC":65,
                "PRIMER_MAX_GC": 80,
                "PRIMER_MIN_SIZE":15,
                "PRIMER_MAX_SIZE":25,
                "PRIMER_OPT_SIZE":18
            },

            "GENOME_Q_ARGS":{
                "PRIMER_OPT_TM": 65,
                "PRIMER_MIN_TM": 55,     
                "PRIMER_MAX_TM": 75,    
                "PRIMER_MIN_GC": 20,
                "PRIMER_OPT_GC":65,
                "PRIMER_MAX_GC": 80,
                "PRIMER_MIN_SIZE":15,
                "PRIMER_MAX_SIZE":25,
                "PRIMER_OPT_SIZE":18
            }     
        }
     ```

**Execute:**
```shell
python edit_main.py
```

**Output:**
- `one_plasmid_system_result.zip`
These files will be generated in the `/home/XXX/tmp/edit_sequence_design/output/` directory.




### 3. Same Plasmid Template - PCR-Based Recombinant Plasmid

**Input:**
- **Step 1:** Upload the genome (fna) file, data preprocessed (CSV) file, sgRNA (CSV) file, plasmid (gb) file.
- **Step 2:** Configure relevant parameters for primer desgin

    - Example configuration (data):

     ```json
        {   
            "chopchop_input": "/home/XXX/tmp/data_preprocessing/output/info_input.csv",   
            "sgRNA_result_path": "/home/XXX/tmp/chopchop/output/sgRNA.csv",
            "edit_sequence_design_workdir":"/home/XXX/tmp/edit_sequence_design/output", 
            "ref_genome": "/home/XXX/tmp/edit_sequence_design/GCA_000011325.1_ASM1132v1_genomic.fna",
            "one_plasmid_file_path":"./input/pMB1-sgRNA-wacJ1.gb",    
            "bowtie_path":"/home/XXX/software/bowtie",
            "scene":"both_sgRNA_primer",
            "plasmid_metod":'0',
            "plasmid_label":{
                "ccdb_label":"HR arm",  
                "promoter_terminator_label":"gRNA ORF",
                "n_20_label":"wacJ",
                "promoter_label":"promoter" 
            },
            "enzyme":{
                "enzyme_name":"BsaI",
                "gap_sequence":"AA",  
                "protection_sequence":"CCA"   
            },
            "primer_json":{},
            "region_label":"",      
            "uha_dha_config": {     
                "max_right_arm_seq_length": 145,     
                "max_left_arm_seq_length": 145,   
                "min_left_arm_seq_length": 145,        
                "min_right_arm_seq_length": 145     
            },
            "UHA_ARGS":{
                "PRIMER_OPT_TM": 65,
                "PRIMER_MIN_TM": 55,  
                "PRIMER_MAX_TM": 75,    
                "PRIMER_MIN_GC": 20,
                "PRIMER_OPT_GC":65,
                "PRIMER_MAX_GC": 80,
                "PRIMER_MIN_SIZE":15,
                "PRIMER_MAX_SIZE":25,
                "PRIMER_OPT_SIZE":18
            },
            "SEQ_ALTERED_ARGS":{
                "PRIMER_OPT_TM": 65,
                "PRIMER_MIN_TM": 55,
                "PRIMER_MAX_TM": 75,  
                "PRIMER_MIN_GC": 20,
                "PRIMER_OPT_GC":65,
                "PRIMER_MAX_GC": 80,
                "PRIMER_MIN_SIZE":15,
                "PRIMER_MAX_SIZE":25,
                "PRIMER_OPT_SIZE":18
            },
            "DHA_ARGS":{
                "PRIMER_OPT_TM": 65,
                "PRIMER_MIN_TM": 55,
                "PRIMER_MAX_TM": 75,
                "PRIMER_MIN_GC": 20,
                "PRIMER_OPT_GC":65,
                "PRIMER_MAX_GC": 80,
                "PRIMER_MIN_SIZE":15,
                "PRIMER_MAX_SIZE":25,
                "PRIMER_OPT_SIZE":18 
            },
            "UP_SGRNA_ARGS": {
                "PRIMER_OPT_TM": 65,
                "PRIMER_MIN_TM": 55,  
                "PRIMER_MAX_TM": 75,    
                "PRIMER_MIN_GC": 20,
                "PRIMER_OPT_GC":65,
                "PRIMER_MAX_GC": 80,
                "PRIMER_MIN_SIZE":15,
                "PRIMER_MAX_SIZE":25,
                "PRIMER_OPT_SIZE":18
            },
            "DOWN_SGRNA_ARGS": {
                "PRIMER_OPT_TM": 65,
                "PRIMER_MIN_TM": 55,  
                "PRIMER_MAX_TM": 75,    
                "PRIMER_MIN_GC": 20,
                "PRIMER_OPT_GC":65,
                "PRIMER_MAX_GC": 80,
                "PRIMER_MIN_SIZE":15,
                "PRIMER_MAX_SIZE":25,
                "PRIMER_OPT_SIZE":18
            },
            "PLASMID_Q_ARGS":{
                "PRIMER_OPT_TM": 65,
                "PRIMER_MIN_TM": 55,  
                "PRIMER_MAX_TM": 75,    
                "PRIMER_MIN_GC": 20,
                "PRIMER_OPT_GC":65,
                "PRIMER_MAX_GC": 80,
                "PRIMER_MIN_SIZE":15,
                "PRIMER_MAX_SIZE":25,
                "PRIMER_OPT_SIZE":18
            },
            "GENOME_Q_ARGS":{
                "PRIMER_OPT_TM": 65,
                "PRIMER_MIN_TM": 55,     
                "PRIMER_MAX_TM": 75,    
                "PRIMER_MIN_GC": 20,
                "PRIMER_OPT_GC":65,
                "PRIMER_MAX_GC": 80,
                "PRIMER_MIN_SIZE":15,
                "PRIMER_MAX_SIZE":25,
                "PRIMER_OPT_SIZE":18
            }     
        }
     ```

**Execute:**
```shell
python edit_main.py
```

**Output:**
- `one_plasmid_system_result.zip`
These files will be generated in the `/home/XXX/tmp/edit_sequence_design/output/` directory.


### 4. Different Plasmid Templates - PCR-Based Recombinant Plasmid
**Input:**
- **Step 1:** Upload the genome (fna) file, data preprocessed (CSV) file, sgRNA (CSV) file, plasmid (gb) file.
- **Step 2:** Configure relevant parameters for primer desgin

    - Example configuration (data):

     ```json
        {   
            "chopchop_input": "/home/XXX/tmp/data_preprocessing/output/info_input.csv",   
            "sgRNA_result_path": "/home/XXX/tmp/chopchop/output/sgRNA.csv",
            "edit_sequence_design_workdir":"/home/XXX/tmp/edit_sequence_design/output", 
            "ref_genome": "/home/XXX/tmp/edit_sequence_design/GCA_000011325.1_ASM1132v1_genomic.fna",
            "no_ccdb_plasmid":"/home/XXX/program/edit_sequence_design/input/no-ccdb-pXMJ19-Cas9A-gRNA-crtYEb-Ts - ori.gb",
            "no_sgRNA_plasmid":"/home/XXX/program/edit_sequence_design/input/no-sgRNA-pXMJ19-Cas9A-gRNA-crtYEb-Ts - ori.gb",   
            "bowtie_path":"/home/XXX/software/bowtie",
            "scene":"both_sgRNA_primer",
            "plasmid_metod":'0',
            "plasmid_label":{
                "ccdb_label":"HR arm",  
                "promoter_terminator_label":"gRNA ORF",
                "n_20_label":"wacJ",
                "promoter_label":"promoter" 
            },
            "enzyme":{
                "enzyme_name":"BsaI",
                "gap_sequence":"AA",  
                "protection_sequence":"CCA"   
            },
            "sgRNA_primer_json":{},
            "ccdb_primer_json":{},   
            "sgRNA_region_label":"",
            "ccdb_region_label":"", 
            "uha_dha_config": {     
                "max_right_arm_seq_length": 145,     
                "max_left_arm_seq_length": 145,   
                "min_left_arm_seq_length": 145,        
                "min_right_arm_seq_length": 145     
            },
            "UHA_ARGS":{
                "PRIMER_OPT_TM": 65,
                "PRIMER_MIN_TM": 55,  
                "PRIMER_MAX_TM": 75,    
                "PRIMER_MIN_GC": 20,
                "PRIMER_OPT_GC":65,
                "PRIMER_MAX_GC": 80,
                "PRIMER_MIN_SIZE":15,
                "PRIMER_MAX_SIZE":25,
                "PRIMER_OPT_SIZE":18
            },
            "SEQ_ALTERED_ARGS":{
                "PRIMER_OPT_TM": 65,
                "PRIMER_MIN_TM": 55,
                "PRIMER_MAX_TM": 75,  
                "PRIMER_MIN_GC": 20,
                "PRIMER_OPT_GC":65,
                "PRIMER_MAX_GC": 80,
                "PRIMER_MIN_SIZE":15,
                "PRIMER_MAX_SIZE":25,
                "PRIMER_OPT_SIZE":18
            },
            "DHA_ARGS":{
                "PRIMER_OPT_TM": 65,
                "PRIMER_MIN_TM": 55,
                "PRIMER_MAX_TM": 75,
                "PRIMER_MIN_GC": 20,
                "PRIMER_OPT_GC":65,
                "PRIMER_MAX_GC": 80,
                "PRIMER_MIN_SIZE":15,
                "PRIMER_MAX_SIZE":25,
                "PRIMER_OPT_SIZE":18 
            },
            "UP_SGRNA_ARGS": {
                "PRIMER_OPT_TM": 65,
                "PRIMER_MIN_TM": 55,  
                "PRIMER_MAX_TM": 75,    
                "PRIMER_MIN_GC": 20,
                "PRIMER_OPT_GC":65,
                "PRIMER_MAX_GC": 80,
                "PRIMER_MIN_SIZE":15,
                "PRIMER_MAX_SIZE":25,
                "PRIMER_OPT_SIZE":18
            },
            "DOWN_SGRNA_ARGS": {
                "PRIMER_OPT_TM": 65,
                "PRIMER_MIN_TM": 55,  
                "PRIMER_MAX_TM": 75,    
                "PRIMER_MIN_GC": 20,
                "PRIMER_OPT_GC":65,
                "PRIMER_MAX_GC": 80,
                "PRIMER_MIN_SIZE":15,
                "PRIMER_MAX_SIZE":25,
                "PRIMER_OPT_SIZE":18
            },
            "PLASMID_Q_ARGS":{
                "PRIMER_OPT_TM": 65,
                "PRIMER_MIN_TM": 55,  
                "PRIMER_MAX_TM": 75,    
                "PRIMER_MIN_GC": 20,
                "PRIMER_OPT_GC":65,
                "PRIMER_MAX_GC": 80,
                "PRIMER_MIN_SIZE":15,
                "PRIMER_MAX_SIZE":25,
                "PRIMER_OPT_SIZE":18
            },
            "GENOME_Q_ARGS":{
                "PRIMER_OPT_TM": 65,
                "PRIMER_MIN_TM": 55,     
                "PRIMER_MAX_TM": 75,    
                "PRIMER_MIN_GC": 20,
                "PRIMER_OPT_GC":65,
                "PRIMER_MAX_GC": 80,
                "PRIMER_MIN_SIZE":15,
                "PRIMER_MAX_SIZE":25,
                "PRIMER_OPT_SIZE":18
            }     
        }
     ```

**Execute:**
```shell
python edit_main.py
```

**Output:**
- `two_plasmid_system_result.zip`
These files will be generated in the `/home/XXX/tmp/edit_sequence_design/output/` directory.




### 5. Different Plasmid Templates - PCR-Based Homologous Arm Sequence Fragment
**Input:**
- **Step 1:** Upload the genome (fna) file, data preprocessed (CSV) file, sgRNA (CSV) file, plasmid (gb) file.
- **Step 2:** Configure relevant parameters for primer desgin

    - Example configuration (data):

     ```json
        {   
            "chopchop_input": "/home/XXX/tmp/data_preprocessing/output/info_input.csv",   
            "sgRNA_result_path": "/home/XXX/tmp/chopchop/output/sgRNA.csv",
            "edit_sequence_design_workdir":"/home/XXX/tmp/edit_sequence_design/output", 
            "ref_genome": "/home/XXX/tmp/edit_sequence_design/GCA_000011325.1_ASM1132v1_genomic.fna",
            "no_ccdb_plasmid":"/home/XXX/program/edit_sequence_design/input/no-ccdb-pXMJ19-Cas9A-gRNA-crtYEb-Ts - ori.gb",
            "no_sgRNA_plasmid":"/home/XXX/program/edit_sequence_design/input/no-sgRNA-pXMJ19-Cas9A-gRNA-crtYEb-Ts - ori.gb",   
            "bowtie_path":"/home/XXX/software/bowtie",
            "scene":"both_sgRNA_primer",
            "plasmid_metod":'0',
            "plasmid_label":{
                "ccdb_label":"HR arm",  
                "promoter_terminator_label":"gRNA ORF",
                "n_20_label":"wacJ",
                "promoter_label":"promoter" 
            },
            "enzyme":{
                "enzyme_name":"BsaI",
                "gap_sequence":"AA",  
                "protection_sequence":"CCA"   
            },
            "ccdb_primer_json":{},   
            "ccdb_region_label":"", 
            "uha_dha_config": {     
                "max_right_arm_seq_length": 145,     
                "max_left_arm_seq_length": 145,   
                "min_left_arm_seq_length": 145,        
                "min_right_arm_seq_length": 145     
            },
            "UHA_ARGS":{
                "PRIMER_OPT_TM": 65,
                "PRIMER_MIN_TM": 55,  
                "PRIMER_MAX_TM": 75,    
                "PRIMER_MIN_GC": 20,
                "PRIMER_OPT_GC":65,
                "PRIMER_MAX_GC": 80,
                "PRIMER_MIN_SIZE":15,
                "PRIMER_MAX_SIZE":25,
                "PRIMER_OPT_SIZE":18
            },
            "SEQ_ALTERED_ARGS":{
                "PRIMER_OPT_TM": 65,
                "PRIMER_MIN_TM": 55,
                "PRIMER_MAX_TM": 75,  
                "PRIMER_MIN_GC": 20,
                "PRIMER_OPT_GC":65,
                "PRIMER_MAX_GC": 80,
                "PRIMER_MIN_SIZE":15,
                "PRIMER_MAX_SIZE":25,
                "PRIMER_OPT_SIZE":18
            },
            "DHA_ARGS":{
                "PRIMER_OPT_TM": 65,
                "PRIMER_MIN_TM": 55,
                "PRIMER_MAX_TM": 75,
                "PRIMER_MIN_GC": 20,
                "PRIMER_OPT_GC":65,
                "PRIMER_MAX_GC": 80,
                "PRIMER_MIN_SIZE":15,
                "PRIMER_MAX_SIZE":25,
                "PRIMER_OPT_SIZE":18 
            },
            "PLASMID_Q_ARGS":{
                "PRIMER_OPT_TM": 65,
                "PRIMER_MIN_TM": 55,  
                "PRIMER_MAX_TM": 75,    
                "PRIMER_MIN_GC": 20,
                "PRIMER_OPT_GC":65,
                "PRIMER_MAX_GC": 80,
                "PRIMER_MIN_SIZE":15,
                "PRIMER_MAX_SIZE":25,
                "PRIMER_OPT_SIZE":18
            },
            "GENOME_Q_ARGS":{
                "PRIMER_OPT_TM": 65,
                "PRIMER_MIN_TM": 55,     
                "PRIMER_MAX_TM": 75,    
                "PRIMER_MIN_GC": 20,
                "PRIMER_OPT_GC":65,
                "PRIMER_MAX_GC": 80,
                "PRIMER_MIN_SIZE":15,
                "PRIMER_MAX_SIZE":25,
                "PRIMER_OPT_SIZE":18
            }     
        }
     ```

**Execute:**
```shell
python edit_main.py
```

**Output:**
- `two_plasmid_system_result.zip`
These files will be generated in the `/home/XXX/tmp/edit_sequence_design/output/` directory.

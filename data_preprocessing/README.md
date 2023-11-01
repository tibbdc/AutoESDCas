  
# data_preprocessing

## Project Introduction  
data_preprocessing, as the data preprocessing module of AutoESDCas, its main function is to convert user input information into standard and standardized input information for AutoESDCas.Among them, for user input information, there are mainly two types supported: 1. the user provides the upstream sequence of the target to be edited on the genome, and 2. the user provides the coordinates of the target to be edited on the genome.At the same time, according to the main functions that AutoESDCas serves users, including: 1. designing only sgRNA, 2. designing only primers, 3. designing both sgRNA and primers, both types of input information also require relevant configurations.


## Installation


### python packages
We suggest using Python 3.8 for data_preprocessing.

```shell
pip install -r requirements.txt

```

### blast+
```shell
wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.13.0+-x64-linux.tar.gz -O ~/ncbi-blast-2.13.0+-x64-linux.tar.gz

tar -zxvf ncbi-blast-2.13.0+-x64-linux.tar.gz

export PATH=~/ncbi-blast-2.13.0+/bin:$PATH

```


## Usage & Example

### 1. User Provides Upstream Sequence of the Target for Genome Editing

**Input:**

- **Step 1:** Upload the genome (fna) file and the target information (CSV) file to be edited.

- **Step 2:** Select from different task types and provide the necessary configuration information.

   1. Designing only sgRNA:
      - Example configuration (data1):
        ```json
        {
            "input_file_path": "./input/designing_only_sgRNA_1.csv",
            "ref_genome": "./input/GCA_000011325.1_ASM1132v1_genomic.fna",
            "data_preprocessing_workdir": "/home/XXX/tmp/data_preprocessing/output/",
            "scene": "only_sgRNA"
        }
        ```

   2. Designing only primers:
      - Example configuration (data2):
        ```json
        {
            "input_file_path": "./input/designing_only_primers_1.csv",
            "ref_genome": "./input/GCA_000011325.1_ASM1132v1_genomic.fna",
            "data_preprocessing_workdir": "/home/XXX/tmp/data_preprocessing/output/",
            "scene": "only_primer"
        }
        ```

   3. Designing both sgRNA and primers:
      - Example configuration (data3):
        ```json
        {
            "input_file_path": "./input/designing_both_sgRNA_and_primers_1.csv",
            "ref_genome": "./input/GCA_000011325.1_ASM1132v1_genomic.fna",
            "data_preprocessing_workdir": "/home/XXX/tmp/data_preprocessing/output/",
            "scene": "both_sgRNA_primer"
        }
        ```

**Execute:**

```shell
python parse_input_to_df.py
```
**Output:**

- `info_input.csv` 
- `xxx.fna` 

These files will be generated in the `/home/XXX/tmp/data_preprocessing/output/` directory.






### 2.the user provides the coordinates of the target to be edited on the genome.


**Input:**

- **Step 1:** Upload the genome (gb) file and the target information (CSV) file to be edited.

- **Step 2:** Select from different task types and provide the necessary configuration information.

   1. Designing only sgRNA:
      - Example configuration (data1):
        ```json
        {
            "input_file_path": "./input/designing_only_sgRNA_2.csv",
            "ref_genome": "./input/eco.fna",
            "data_preprocessing_workdir": "/home/XXX/tmp/data_preprocessing/output/",
            "scene": "only_sgRNA"
        }
        ```

   2. Designing only primers:
      - Example configuration (data2):
        ```json
        {
            "input_file_path": "./input/designing_only_primers_2.csv",
            "ref_genome": "./input/eco.fna",
            "data_preprocessing_workdir": "/home/XXX/tmp/data_preprocessing/output/",
            "scene": "only_primer"
        }
        ```

   3. Designing both sgRNA and primers:
      - Example configuration (data3):
        ```json
        {
            "input_file_path": "./input/designing_both_sgRNA_and_primers_2.csv",
            "ref_genome": "./input/eco.fna",
            "data_preprocessing_workdir": "/home/XXX/tmp/data_preprocessing/output/",
            "scene": "both_sgRNA_primer"
        }
        ```

**Execute:**

```shell
python parse_input_to_df.py
```
**Output:**

- `info_input.csv` 
- `xxx.fna` 

These files will be generated in the `/home/XXX/tmp/data_preprocessing/output/` directory.










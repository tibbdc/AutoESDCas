# chopchop
## Project Introduction 
As the sgRNA design module of AutoESDCas, chopshop's main function is to high-throughput complete sgRNA design by calling the core algorithm of the original chopshop sgRNA design.


## Installation
The installation of the original chopshop requires the establishment of a virtual environment for Python 2.7, as well as the installation of tools such as Bowtie and twoBitToFa. On this basis, the virtual environment of python 3.8 and the installation of its related packages are established, and the original chop shop of the python 2.7 environment is finally called in the python 3.8 environment.

### python packages
We suggest using Python 3.8 and  Python 2.7 for chopchop.

**Python 2.7:**
```shell
pip install -r chopchop_requirements.txt
```

### Bowtie   
```
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.3.1/bowtie-1.3.1-linux-x86_64.zip -O /opt/bowtie.zip && \
    unzip /opt/bowtie.zip -d /opt/ && \
    mv /opt/bowtie-1.3.1-linux-x86_64  /opt/bowtie && \
    rm /opt/bowtie.zip
export PATH=~/software/bowtie/:$PATH
```

### twoBitToFa
```
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa -O /opt/twoBitToFa
export PATH=~/software/twoBitToFa/:$PATH
```

**Python 3.8:**
```shell
pip install -r requirements.txt
```

## Usage & Example

**Input:**
- **Step 1:** Upload the genome (fna) file and data preprocessed (CSV) file.

- **Step 2:** Configure relevant parameters for sgRNA design

    - Example configuration (data):
    ```json
    {
        "input_file_path":"/home/XXX/tmp/data_preprocessing/output/info_input.csv",
        "ref_genome":"/home/XXX/tmp/data_preprocessing/output/xxx.fna",
        "chopchop_workdir":"/home/XXX/tmp/chopchop/output/", 
        "chopchop_config":{
            "PAM": "NNNNGMTT", 
            "guideSize": 20,
            "maxMismatches": 3,
            "scoringMethod": "XU_2015"
        }
    }
    ```


**Execute:**
```shell
python chopchop_main.py
```

**Output:**
- `sgRNA.json` 
- `sgRNA.csv` 

These files will be generated in the `/home/XXX/tmp/chopchop/output/` directory.





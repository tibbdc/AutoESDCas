# AutoESDCas
AutoESDCas, an online tool for the whole-workflow editing sequence design for microbial genome editing based on CRISPR/Cas system. This tool facilitates all types of genetic manipulation and different CRISPR/Cas-mediated genome editing technique variants, enables biologists to quickly and efficiently obtain all editing sequences needed for the entire genome editing process, and empowers high-throughput strain modification.
   
## Project Introduction

**AutoESDCas** an online tool for the whole-workflow editing sequence design for microbial genome editing based on CRISPR/Cas system. This tool facilitates all types of genetic manipulation and different CRISPR/Cas-mediated genome editing technique variants, enables biologists to quickly and efficiently obtain all editing sequences needed for the entire genome editing process, and empowers high-throughput strain modification.

## Modules
AutoESDCas consists of the following modules:
1. [Data Preprocessing](#data_preprocessing)
2. [ChopChop](#chopchop)
3. [Editing Sequence Design](#editing_sequence_design)


## Data Preprocessing

**Data Preprocessing** serves as the data preprocessing module of AutoESDCas. Its primary function is to convert user input information into standardized and structured input data for AutoESDCas. It supports two main types of user input information:

- The user provides the upstream sequence of the target to be edited on the genome.
- The user provides the coordinates of the target to be edited on the genome.

[Learn more about Data Preprocessing](#link-to-data_preprocessing)


## chopchop

**chopchop** serves as the sgRNA design module of AutoESDCas. Its primary function is to perform high-throughput sgRNA design by utilizing the core algorithm of the original ChopChop sgRNA design. This process is based on the results obtained from the Data Preprocessing module.

[Learn more about ChopChop](#link-to-chopchop)

## editing_sequence_design

**editing_sequence_design** serves as the primer design module of AutoESDCas, with a primary focus on automating and optimizing the primer design process for genome editing. This module builds upon the results obtained from the Data Preprocessing and ChopChop modules. It offers a range of features, including the design of homologous arms, PCR primers, sequencing primers, and the visualization of recombinant plasmid maps.

[Learn more about Editing Sequence Design](#link-to-editing_sequence_design)

## Installation

AutoESDCas requires various Python packages and additional tools. Detailed installation instructions can be found in the respective module documentation.

## Usage & Examples

Each module within AutoESDCas comes with its own usage instructions and example configurations. Refer to the documentation of each module for more details.

Feel free to explore each module to learn more about its capabilities and how it can assist you in your genome editing projects.




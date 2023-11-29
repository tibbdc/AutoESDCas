# Tutorial

## User Manual of AutoESDCas WEB Server

### 1.Introduction

**AutoESDCas** an online tool for the whole-workflow editing sequence design for microbial genome editing based on CRISPR/Cas system. This tool facilitates all types of genetic manipulation and different CRISPR/Cas-mediated genome editing technique variants, enables biologists to quickly and efficiently obtain all editing sequences needed for the entire genome editing process, and empowers high-throughput strain modification. **AutoESDCas** offer two essntial modules:  

- **sgRNA Design** serves as the sgRNA design module of AutoESDCas. Its primary function is to perform high-throughput sgRNA design by utilizing the core algorithm of the original chopshop sgRNA design.
- **Genome Editing Design** serves as the primer design module of AutoESDCas, with a primary focus on automating and optimizing the primer design process for genome editing. It offers a range of features, including the design of homologous arms, PCR primers, sequencing primers, and the visualization of recombinant plasmid maps. This module is designed to accommodate various experimental conditions, catering to the diverse requirements of biological researchers.
  Supported scenarios include:  
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

Additionally, **AutoESDCas** offers multiple flexible usage options to cater to diverse user needs:

1. **Exclusive focus on sgRNA design**:
   - In this mode, users can swiftly and accurately design sgRNA sequences tailored to their experimental requirements.
   - Opt for this approach if sgRNA design is your primary concern.
2. **Sole emphasis on genome editing**:  
   - Ideal for users exclusively interested in genome editing, this option allows direct execution of genome editing and design without dealing with sgRNA-related steps.Choose this mode for a streamlined genome editing experience.
3. **Comprehensive consideration of sgRNA and genome editing**:
   - The most comprehensive choice for users needing a holistic approach. Initiate sgRNA design first and then delve into a more in-depth genome editing design based on these initial designs.
   - Opt for this approach if you require a thorough integration of sgRNA and genome editing.

### 2.PRIMARY FEATURES OFFERED BY THE WEBSITE

#### 2.1 Exclusive focus on sgRNA design

In this mode, AutoESDCas provides high-throughput sgRNA design functionality and supports two kinds of input, allowing end users to easily and quickly obtain sgRNA.

- The user provides the upstream sequence of the target to be edited on the genome.  
![Watch the video](https://autoesdcas.s3.amazonaws.com/uploads/only_sgRNA_1.mp4)
<video src="https://autoesdcas.s3.amazonaws.com/uploads/only_sgRNA_1.mp4" data-canonical-src="https://autoesdcas.s3.amazonaws.com/uploads/only_sgRNA_1.mp4" controls="controls" muted="muted" class="d-block rounded-bottom-2 width-fit" style="max-height:640px;">
- The user provides the coordinates of the target to be edited on the genome.  

#### 2.2 Sole emphasis on genome editing

In this mode,AutoESDCas provides support for five different genome editing experimental scenarios based on CRISPR/Cas technology, achieving functions such as primer design and visualization of recombinant plasmids in relevant scenarios, and still supporting  two kinds of input, in each scenario.In addition, AutoESDCas still provides three scenarios for primer design of recombinant plasmid fragments and sgRNA fragments in the scenario of constructing recombinant plasmids based on PCR method, to meet the needs of different users.

- Same Plasmid Template - Synthesis Recombinant Plasmid
  - The user provides the upstream sequence and sgRNA of the target to be edited on the genome.
  - The user provides the coordinates and sgRNA of the target to be edited on the genome.
- Different Plasmid Templates - Synthesis Recombinant Plasmid
  - The user provides the upstream sequence and sgRNA of the target to be edited on the genome.
  - The user provides the coordinates and sgRNA of the target to be edited on the genome.
- Same Plasmid Template - PCR-Based Recombinant Plasmid
  - The user provides the upstream sequence and sgRNA of the target to be edited on the genome.
    - i. System automatic design
      - ii. Manual provision of partial primers
      - iii. Manual provision of primer design range.
- Different Plasmid Templates - PCR-Based Recombinant Plasmid
  - The user provides the upstream sequence and sgRNA of the target to be edited on the genome.
    - i. System automatic design
      - ii. Manual provision of partial primers
      - iii. Manual provision of primer design range.
- Different Plasmid Templates - PCR-Based Homologous Arm Sequence Fragment
  - The user provides the upstream sequence and sgRNA of the target to be edited on the genome.
    - i. System automatic design
      - ii. Manual provision of partial primers
      - iii. Manual provision of primer design range.
  - Additionally, sgRNA sequence fragments are obtained through primer annealing, and the recombinant plasmid fragments where the sgRNA is located are obtained through enzymatic cleavage.

#### 2.3 Comprehensive consideration of sgRNA and genome editing

In this mode,AutoESDCas provides automated design functions for editing sequences required for the entire genome editing process, including sgRNA design, related primer design, etc., and supports two kinds of input.

- The user provides the upstream sequence of the target to be edited on the genome.
- The user provides the coordinates of the target to be edited on the genome.

### 3. Input and output files

#### 3.1 Exclusive focus on sgRNA design

##### 3.1.1 Input files

AutoESDCas's sgRNA design task supports two types of input files. The first group consists of the GB file of the target genome, containing information on all genes, and a CSV file with the gene ID and relative distances for the edited targets. The CSV file includes a Name column (ensuring information uniqueness), a Gene id column, an Up region column (relative distance from the start coordinate of the gene ID in the genome), and a Down region column (relative distance from the end coordinate of the gene ID). These last three columns are used to locate the editing target on the genome.

The second group of input files includes the FASTA file of the target genome, containing all sequences, and a CSV file with the upstream sequence of the edited target. The CSV file includes a Name column (ensuring information uniqueness), a Sequence upstream of the manipulation site (>100bp) column, and a Reference sequence column (wild-type sequence). The latter two columns are used to locate the editing target on the genome.

![1](https://github.com/tibbdc/AutoESDCas/blob/main/doc/img/1.png)

**Figure 1  Input file format. A1.** Example in a GB input file for the targeted genome. **B1.** The CSV input file for the standard sequence manipulation sgRNA design task includes a parts list with user-provided genome target coordinates for the desired editing location. **A2** Example in a Fasta input file for the targeted genome. **B2** The CSV input file for the standard sequence manipulation sgRNA design task includes a parts list, with users providing information about the upstream sequence of the target to be edited on the genome.

##### 3.1.2 Output files

The output of AutoESDCas's sgRNA design task is an XLSX file containing all sgRNAs. The file includes columns such as Name and Rank, where the sgRNA ranking is determined using the ChopChop sgRNA design algorithm.

![1701157308621](C:\Users\TIBD_L~1\AppData\Local\Temp\1701157308621.png)

**Figure 2  Output file of format.** Example of design results, list in XLSX output file.

#### 3.2 Sole emphasis on genome editing

##### 3.2.1 Input files

AutoESDCas's Genome Editing Design Task supports two types of input files. In the first group, a GB file of the target genome, containing information on all genes, is accompanied by a CSV file specifying the gene ID and relative distance for the edited targets. Additionally, a GB file for the plasmid template is required. The CSV file, an extension of the structure in Figure 1 B1, now includes columns for "Inserted Sequence," "Manipulation type," and "crRNA." The "Manipulation type" column can have values such as "substitution," "insertion," or "deletion."

The second group of input files includes the FASTA file of the target genome, containing all sequences, along with a CSV file detailing the upstream sequence of the edited target. Similar to the first group, a GB file for the plasmid template is required. The CSV file, an extension of the structure in Figure 1 B2, includes columns for "Inserted Sequence," "Manipulation type," and "crRNA."

![1701157538029](C:\Users\TIBD_L~1\AppData\Local\Temp\1701157538029.png)

**Figure 3  Input file format. A1.** Example in a GB input file for the targeted genome. **B1.** The CSV input file for the standard sequence manipulation Genome design task includes a parts list with user-provided genome target coordinates for the desired editing location.**C1** Example in a GB Input File for the Targeted Plasmid. **A2** Example in a Fasta input file for the targeted genome. **B2** The CSV input file for the standard sequence manipulation Genome design task includes a parts list, with users providing information about the upstream sequence of the target to be edited on the genome.**C2** Example in a GB Input File for the Targeted Plasmid.

###### 3.2.2 Output files

The output of AutoESDCas's Genome Editing Design Task includes the following files in xlsx format, providing comprehensive information for different scenarios:

- Design_result:Contains all primer information for the designed task.

  In the Synthesis Recombinant Plasmid scenario, it includes sequencing primers for the recombinant plasmid (Test_primer_P1), genome sequencing primers (Test_primer_G), off-target sequencing primers for the recombinant plasmid (Primer_p_offTarget), off-target sequencing primers for the genome (Primer_g_offTarget), and primer order information (Primer_order).

  In the PCR-Based Recombinant Plasmid scenario, it additionally includes upstream homologous arm PCR primers (Primer_UHA), downstream homologous arm PCR primers (Primer_DHA), PCR primers for the inserted fragment (Primer_inserted_fragment), sgRNA fragment primers (Primer_sgRNA_fragment), and plasmid backbone primers (Primer_plasmid_backbone).

- Evaluation:Provides results of off-target analysis for homologous arms.

- Plasmid_order:Includes ordering information for synthesizing the recombinant plasmid.

- plasmid directory:Contains all GenBank files for the synthesized recombinant plasmids.

The output for the "Same Plasmid Template - PCR-Based Recombinant Plasmid" application scenario is exemplified in various sections of Figure 4:In Figure A, an example is presented that provides details regarding PCR primers, specifically those associated with the upstream homologous arm. The information encompasses the ID of the homologous arm primer, the PCR primer pair designed for the homologous arm, the resulting product of the PCR amplification, and the length of the generated product.In Figure B, an example illustrates the output of the evaluation of off-target sequencing primers for the genome. This includes details such as the primer ID, information about occurrences of off-target effects, site ID, sequence, distance from the target site, mismatch details (NM), and the frequency of occurrences on the genome.Figure C illustrates an example providing details about sequencing primers designed specifically for the recombinant plasmid. The information includes specifics about the primer itself, such as its sequence, and its associated melting temperature (TM) value.In Figure E, information related to the off-target analysis for homologous arms is presented.Figure F provides an example of a GenBank (gb) file for the recombinant plasmid. This file encompasses characteristic details about the PCR primers used, the homologous arms, and features related to sgRNA.

![1701173503274](C:\Users\TIBD_L~1\AppData\Local\Temp\1701173503274.png)

**Figure 4  Examples of the output files . ** **A** Example in Design_result for upstream homologous arm (Primer_UHA) PCR primers.**B** Example in Design_result for Off-Target sequencing primers for the genome (Primer_g_offTarget).**C **Example in Design_result for recombinant plasmid sequencing primers（Test_primer_P1）.**D** Example in Design_result for primer ordering. **E** Example in Evaluation for Off-Target analysis report for homologous arms. **F** Example in gb files for recombinant plasmid.



#### 3.3 Comprehensive consideration of sgRNA and genome editing

AutoESDCas's sgRNA design and Genome Editing Design task supports two types of input files. In this mode (referring to the context of Figure 4, columns B1 and B2), the input format differs from the "Sole emphasis on genome editing" mode by excluding the `crRNA` columns found in Figure 2 (B1, B2).And,this mode's output can be referenced using Figure 2 and Figure 4.

![1701176981783](C:\Users\TIBD_L~1\AppData\Local\Temp\1701176981783.png)

**Figure 5  Input file format. A1.** Example in a GB input file for the targeted genome. **B1.** The CSV input file for the standard sequence manipulation Genome design task includes a parts list with user-provided genome target coordinates for the desired editing location.**C1** Example in a GB Input File for the Targeted Plasmid. **A2** Example in a Fasta input file for the targeted genome. **B2** The CSV input file for the standard sequence manipulation Genome design task includes a parts list, with users providing information about the upstream sequence of the target to be edited on the genome.**C2** Example in a GB Input File for the Targeted Plasmid.

### 4.Browser compatibility

| OS            | Version        | Version	Chrome | Firefox       | Microsoft Edge | Safari |
| ------------- | -------------- | ----------------- | ------------- | -------------- | ------ |
| Linux         | Ubuntu18.04    | 96.0.4664.93      | 92.0          | n/a            | n/a    |
| MacOS         | Ventura 13.0.1 | 108.0.5359.124    | 107.0.1       | 108.0.1462.54  | 16.1   |
| Windows	10 | 108.0.5359.125 | 108.0.1           | 108.0.1462.54 | n/a            |        |

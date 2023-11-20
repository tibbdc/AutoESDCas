# Tutorial
## User Manual of AutoESDCas WEB Server

### 1.Introduction
**AutoESDCas** an online tool for the whole-workflow editing sequence design for microbial genome editing based on CRISPR/Cas system. This tool facilitates all types of genetic manipulation and different CRISPR/Cas-mediated genome editing technique variants, enables biologists to quickly and efficiently obtain all editing sequences needed for the entire genome editing process, and empowers high-throughput strain modification. **AutoESDCas** offer two essntial modules:
**sgRNA Design** serves as the sgRNA design module of AutoESDCas. Its primary function is to perform high-throughput sgRNA design by utilizing the core algorithm of the original chopshop sgRNA design.
**Genome Editing Design** serves as the primer design module of AutoESDCas, with a primary focus on automating and optimizing the primer design process for genome editing. It offers a range of features, including the design of homologous arms, PCR primers, sequencing primers, and the visualization of recombinant plasmid maps. This module is designed to accommodate various experimental conditions, catering to the diverse requirements of biological researchers.
Supported scenarios include:
1. **Same Plasmid Template - Synthesis Recombinant Plasmid:**
   - The sgRNA and homologous arm exist on the same plasmid template.
   - Recombinant plasmid synthesis is employed.

1. **Different Plasmid Templates - Synthesis Recombinant Plasmid:**
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
**AutoESDCas** offers multiple flexible usage options to cater to diverse user needs:
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
In this mode, AutoESDCas provides high-throughput sgRNA design functionality and supports two types of input methods, allowing end users to easily and quickly obtain sgRNA.
1. The user provides the upstream sequence of the target to be edited on the genome.
<video width="640" height="360" controls>
        <source src="./video/only_sgRNA_1.mp4" type="video/mp4">
        Fig.1 sgRNA design.
</video>
    2. The user provides the coordinates of the target to be edited on the genome.
<video width="640" height="360" controls>
        <source src="./video/only_sgRNA_1.mp4" type="video/mp4">
        Fig.2 sgRNA design.
</video>  
#### 2.2 Sole emphasis on genome editing


#### 2.3 Comprehensive consideration of sgRNA and genome editing

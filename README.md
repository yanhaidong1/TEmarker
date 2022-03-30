

TEmarker
===
TEmarker developes a marker system based on transposable elements (Transposons; TEs). This tool aims to generate a genotype file using TE markers that could be applied in population genetic analysis.


# Introduction
In this study, we developed a TE marker system (TEmarker) on the basis of whole-genome shot gun sequencing data. This system contains three steps: 1) construct pan-TE insertions for all individuals in a population, and combine the close insertions together to generate candidate TE markers; 2) determine marker genotypes based on a proportion of split reads mapping on the border region of the TE insertion sites; 3) generate a ‘VCF’ output file used for the downstream analysis such as population structure and GWAS analyses. This tool can combine results from different (TE insertion polymorphisms) TIP calling tools, remove potential false positive TE insertions, and conduct high-throughput genotyping.

# Dependence and requirements
TEmarker is developed in python with modules and external tools.

Before running this pipeline, a dependency check should be performed first to make sure every dependency is correctly installed.

For information about installing the dependencies, please see below. The version numbers listed below represents the version this pipeline is developed with, and using the newest version is recommended.

## Requirements
**Python** (v3.0 or more; Developed and tested with version 3.7.1)  
Modules can be installed using pip: pip install -r requirements.txt or pip install [module_name]  

**Module version**  
biopython (Developed and tested with version 1.72)  

**hisat2** (Developed and tested with version 2.1.0)  

**samtools** (Developed and tested with version 0.1.19-44428cd)  

**McClintock**  
https://github.com/bergmanlab/mcclintock  


## Optional requirements

**CNVnator** (Tested with version 0.3.2)  

**bwa** (Tested with version 0.7.10-r789)  
Available at https://github.com/lh3/bwa.  

**Any other polymorphic TE detection tools**

# Quick start
Four steps should be running to finish the whole process:  
1.	Step 0: Prepare input data
2.	Step 1: Create pan-TEs
3.	Step 2: Genotyping
4.	Step 3: Generate output  

Please check **Manual.pdf** to learn to condcut **Step 0**.
Here we provide testing data for all the steps users can download following the command lines.

## Download testing data  
### Option 1: Users can download data from the following google drive link:  
Example_dir_step0.tar.gz  
https://drive.google.com/file/d/1TwaDBGYQ2OoHIENUgGMfQF8eiO3T1aZV/view?usp=sharing  

Example_dir_step123.zip  
https://drive.google.com/file/d/1N9hb1WFWdzrP0kQ4v_uxaH8VfeOjZS2n/view?usp=sharing  

### Option 2: Users can download data by using command lines:
**Commands (Step 1,2,3 data):**  
* > mkdir output_dir  
* > Download_test_data.py -o output_dir -s123 yes

**Outputs (Step 1,2,3 data):**
1.	Example_dir_step123 (test_genome.fasta; TE.lib; bed_dir; bwa_bam_dir; hisat2_bam_dir; material_name.txt)  
**Note:** Six files or directories shown as follows need to be prepared for Step 1, 2, 3, and their test data can be downloaded by ‘Download_test_data.py':  
1). A test genome file (**test_genome.fasta**; 5 Mb) derived from a part of rice genome (Os-Nipponbare-Reference-IRGSP-1.0).  
2). A test TE library file (**TE.lib**).  
3). A test directory including bed files that contain TE insertion information for each sample (**bed_dir**).   
4). A test directory including bam files derived from bwa tool for each sample (**bwa_bam_dir**).  
5). A test directory including bam files derived from hisat2 tool for each sample (**hisat2_bam_dir**).  
6). A test material name file (**material_name.txt**) that is used to change original name (the first column) to a new name (the second column).  

## **Step 1,2,3: Create pan-TEs; Genotyping; Generate output**
We provide an example to use the **TEmarker_pipeline.py** to run **Step 1, 2, and 3** after finishing the **Step 0**. This pipeline uses the default parameter setting, if users want to change parameters, please run these three steps seperately shown in **Manual.pdf**.

**Inputs:**
1.	test_genome.fasta
2.	TE.lib
3.	bed_dir
4.	bwa_bam_dir
5.	hisat2_bam_dir
6.	material_name.txt

**Commands:**  
* > mkdir working_dir  
mkdir output_dir  
* > TEmarker_pipeline.py \\  
-d working_dir -o output_dir \\  
-fas_f test_genome.fasta \\  
-lib_f TE.lib \\   
-b_d bed_dir \\  
-bam_bwa_d bwa_bam_dir \\  
-bam_hisat2_d hisat2_bam_dir \\  
-m_f material_name.txt \\  
-TEmarker_p /Path/to/TEmarker \\  
-pros_n 20

**Outputs:**  
1.	vcf_output_dir  
**opt.vcf** (no filtration)  
**opt_bi.vcf** (transfer ‘0/1’ to ‘1/1’)  
**opt_fltmissing_fltmaf.vcf** (remove loci based on missing sample and maf thresholds)  
3.	genos_output_dir  
**opt_genos.txt**  
**opt_genos_fltmissing_fltmaf.txt**  
4.	summary of TE families  
**opt_summary_te_family.txt**  


**Note:**
The tested outputs can be downloaded in the '**Test_output_dir**'.



## Please follow the 'Manual.pdf' to learn more information


## Appendix and FAQ

:::info
**Find this document incomplete?** Leave a comment!
:::



#!/usr/bin/env python

##updating 082921 add an option to close the bam modification
##updation 122420
##this script will be used to run step 1,2,3

##This wrapper will helps users to generate the vcf file for the population genomic analysis
##BUILT-IN MODULES
import argparse
import sys
import os
import subprocess
from distutils.spawn import find_executable

def get_parsed_args():
    parser = argparse.ArgumentParser(description="TEmarker pipeline for step 1, 2, and 3")

    ##############################
    ##parse the required arguments
    parser.add_argument("-d", dest='working_dir', help="Working directory"
                                                        "Default: ./ ")

    parser.add_argument("-o", dest='output_dir',  help="Output directory"
                                                        "Default: ./ ")

    parser.add_argument('-fas_f',dest='genome_fas',help="Provide the genome fasta file.")

    parser.add_argument('-lib_f', dest='te_library', help="Provide the TE library.")

    parser.add_argument('-b_d', dest='bed_dir',help="Users provide the bed dir that will be used to generate the pan TE insertions")

    parser.add_argument('-bam_bwa_d', dest='bam_bwa_dir',help="Users provide the bam dir generated from the TEScape_marker_bam.")

    parser.add_argument('-bam_hisat2_d', dest='bam_hisat2_dir',help="Users provide the bam dir generated from the TEScape_marker_bam.")

    parser.add_argument('-TEmarker_p', dest='TEmarker_path', help='Users provide the TEmarker directory.')

    parser.add_argument('-m_f', dest='material_file', help="If users provide the name_file, and pipeline will substitue the fastq name "
                                                       "to the name of the second column in the name_file.")

    parser.add_argument('--samtools', dest='samtools', default="/usr/bin/samtools",help="Users should provide the path of samtools.")


    ##optional
    parser.add_argument('-modify', dest='m_yes',help="Users activate the modification on the genotyping file based on the CNV detection dir")

    parser.add_argument('-cnv_d', dest='cnv_dir',help="After activation of modification, users must provide the cnv_dir generated from TEScape_marker_genos_cnv")

    parser.add_argument('-pros_n', dest='process_num',help="Users provide process number that helps to speedup the pipeline."
                                                           "Default: 10")

    parser.add_argument('-clean', dest='clean_temp_dir',help='If user set the yes, all the temp dir in the working_dir will be deletec to save the place')

    parser.add_argument('-closebam', dest='closebam_yes', help="Users close to modify bam files")

    ##parse of parameters
    args = parser.parse_args()
    return args


def main(argv=None):
    if argv is None:
        argv = sys.argv
    args = get_parsed_args()

    #######################################
    ##check the required software and files
    ##for the input files
    bed_dir = args.bed_dir
    if args.bed_dir is None:
        print('Cannot find input te bed dir, please provide that')
        return  ##import to close the script if not find the te lib

    bam_bwa_dir = args.bam_bwa_dir
    if args.bam_bwa_dir is None:
        print('Cannot find input te bam bwa dir, please provide that')
        return  ##import to close the script if not find the te lib

    bam_hisat2_dir = args.bam_hisat2_dir
    if args.bam_hisat2_dir is None:
        print('Cannot find input te bam hisat dir, please provide that')
        return  ##import to close the script if not find the te lib

    TE_lib_fl = args.te_library
    if args.te_library is None:
        print ('Cannot find input TE library file, please provide that')
        return ##import to close the script if not find the te lib
    else:
        try:
            file = open(args.te_library, 'r')  ##check if the file is not the right file
        except IOError:
            print('There was an error opening the file!')
            return

    genome_file = args.genome_fas
    if args.genome_fas is None:
        print('Cannot find input genome file, please provide that')
        return  ##import to close the script if not find the te lib
    else:
        try:
            file = open(args.genome_fas, 'r')  ##check if the file is not the right file
        except IOError:
            print('There was an error opening the file!')
            return

    cnv_dir = ''
    ##if args.m_yes has been provided
    if args.m_yes is not None:
        if args.cnv_dir is None:
            print('Please provide the cnv_dir')
            return
        else:
            cnv_dir = args.cnv_dir

    material_file = args.material_file
    if args.material_file is None:
        print('Cannot find input material file, please provide that')
        return  ##import to close the script if not find the te lib
    else:
        try:
            file = open(args.material_file, 'r')  ##check if the file is not the right file
        except IOError:
            print('There was an error opening the file!')
            return

    TEmarker_path = args.TEmarker_path
    if args.TEmarker_path is None:
        print('Cannot find input TEmarker directory path, please provide that')
        return  ##import to close the script if not find the te lib

    if args.process_num is not None:
        process_num = args.process_num
    else:
        process_num = 10

    if find_executable(args.samtools) is not None:
        samtools_exe = args.samtools
        print ('samtools executable can be found')
    else:
        print("Cannot find samtools executable, please check if it has been installed.")
        return

    output_dir = args.output_dir
    if not output_dir.endswith('/'):
        output_dir = output_dir + '/'
    else:
        output_dir = output_dir

    working_dir = args.working_dir
    if not working_dir.endswith('/'):
        working_dir = working_dir + '/'
    else:
        working_dir = working_dir

    ####################
    ##start the pipeline
    ########################
    ##Step 1: Create pan-TEs
    print('Step 1: Create pan-TEs')
    ##Step 1.1: Modify bam files
    print('Step 1.1: Modify bam files')

    ##create working and output directories
    D01_1_modify_bam = working_dir + '/01_1_modify_bam'
    if not os.path.exists(D01_1_modify_bam):
        os.makedirs(D01_1_modify_bam)

    D01_1_modify_bam_od = D01_1_modify_bam + '/output_dir'
    if not os.path.exists(D01_1_modify_bam_od):
        os.makedirs(D01_1_modify_bam_od)

    D01_1_modify_bam_wd = D01_1_modify_bam + '/working_dir'
    if not os.path.exists(D01_1_modify_bam_wd):
        os.makedirs(D01_1_modify_bam_wd)

    ##updating 082921
    if args.closebam_yes == 'yes':

        print ('Users choose to close modifying bam files')

    else:
        ##run the TEmarker_bam.py
        cmd = 'python ' + TEmarker_path + '/TEmarker_bam.py -d ' + D01_1_modify_bam_wd + ' -o ' + D01_1_modify_bam_od + \
              ' -bam_bwa_d ' + bam_bwa_dir + \
              ' -bam_hisat2_d ' + bam_hisat2_dir + \
              ' --samtools ' + samtools_exe + \
              ' -pros_n_samtls ' + str(process_num)
        print(cmd)
        subprocess.call(cmd,shell=True)

    if args.clean_temp_dir is not None:
        if args.clean_temp_dir == 'yes':
            ##check if bed_dir and two bam files exit
            cmd = 'rm -rf ' + D01_1_modify_bam_wd
            print(cmd)
            subprocess.call(cmd,shell=True)
        else:
            print("please type 'yes' behind clean argument.")
            return

    ##Step 1.2: Create pan-TEs
    print('Step 1.2: Create pan-TEs')

    ##create working_dir and output_dir
    D01_2_create_panTEs = working_dir + '/01_2_create_panTEs'
    if not os.path.exists(D01_2_create_panTEs):
        os.makedirs(D01_2_create_panTEs)

    D01_2_create_panTEs_wd = D01_2_create_panTEs + '/working_dir'
    if not os.path.exists(D01_2_create_panTEs_wd):
        os.makedirs(D01_2_create_panTEs_wd)

    D01_2_create_panTEs_od = D01_2_create_panTEs + '/output_dir'
    if not os.path.exists(D01_2_create_panTEs_od):
        os.makedirs(D01_2_create_panTEs_od)

    ##run the TEmarker_panTEs.py
    cmd = 'python ' + TEmarker_path + '/TEmarker_panTEs.py -d ' + D01_2_create_panTEs_wd + ' -o ' + D01_2_create_panTEs_od + \
          ' -b_d ' + bed_dir
    print(cmd)
    subprocess.call(cmd,shell=True)

    ####################
    ##Step 2: Genotyping
    print('Step 2: Genotyping')
    ##create working and output directories
    D02_genotyping = working_dir + '/02_genotyping'
    if not os.path.exists(D02_genotyping):
        os.makedirs(D02_genotyping)

    D02_genotyping_wd = D02_genotyping + '/working_dir'
    if not os.path.exists(D02_genotyping_wd):
        os.makedirs(D02_genotyping_wd)

    D02_genotyping_od = D02_genotyping + '/output_dir'
    if not os.path.exists(D02_genotyping_od):
        os.makedirs(D02_genotyping_od)


    ##updating 082921
    if args.closebam_yes == 'yes':
        ipt_bam_bwa_dir = bam_bwa_dir
        ipt_bam_hisat2_dir = bam_hisat2_dir
    else:
        ipt_bam_bwa_dir = D01_1_modify_bam_od + '/rm_pcr_bwa_bam_dir'
        ipt_bam_hisat2_dir = D01_1_modify_bam_od + '/rm_pcr_hisat2_bam_dir'


    ##run TEmarker_genotyping.p
    ##if users do not provide the cnv directories
    if args.m_yes is None:
        cmd = 'python ' + TEmarker_path + '/TEmarker_genotyping.py -d ' + D02_genotyping_wd + ' -o ' + D02_genotyping_od + \
              ' -panTEs_f ' + D01_2_create_panTEs_od + '/opt_pan_TEs.txt' + \
              ' -bam_bwa_d ' + ipt_bam_bwa_dir + \
              ' -bam_hisat2_d ' + ipt_bam_hisat2_dir + \
              ' -lib_f ' + TE_lib_fl + \
              ' -pros_n ' + str(process_num)
        print(cmd)
        subprocess.call(cmd,shell=True)

    ##if users provide the cnv directories
    else:
        cmd = 'python ' + TEmarker_path + '/TEmarker_genotyping.py -d ' + D02_genotyping_wd + ' -o ' + D02_genotyping_od + \
              ' -canTE_f ' + D01_2_create_panTEs_od + '/opt_pan_TEs.txt' + \
              ' -bam_bwa_d ' + ipt_bam_bwa_dir + \
              ' -bam_hisat2_d ' + ipt_bam_hisat2_dir + \
              ' -lib_f ' + TE_lib_fl + \
              ' -pros_n ' + str(process_num) + \
              ' -modify yes' + \
              ' -cnv_d ' + cnv_dir
        print(cmd)
        subprocess.call(cmd,shell=True)

    #########################
    ##Step 3: Generate output
    print('Step 3: Generate output')
    ##create working and output directories
    D03_generate_output = working_dir + '/03_generate_output'
    if not os.path.exists(D03_generate_output):
        os.makedirs(D03_generate_output)

    D03_generate_output_wd = D03_generate_output + '/working_dir'
    if not os.path.exists(D03_generate_output_wd):
        os.makedirs(D03_generate_output_wd)

    D03_generate_output_od = D03_generate_output + '/output_dir'
    if not os.path.exists(D03_generate_output_od):
        os.makedirs(D03_generate_output_od)

    ##run TEmarker_output.py
    if args.m_yes is None:
        cmd = 'python ' + TEmarker_path + '/TEmarker_output.py -d ' + D03_generate_output_wd + ' -o ' + D03_generate_output_od + \
              ' -genos_f ' + D02_genotyping_od + '/opt_te_genotype.txt' + \
              ' -m_f ' + material_file + \
              ' -fas_f ' + genome_file
        print(cmd)
        subprocess.call(cmd,shell=True)
    else:
        cmd = 'python ' + TEmarker_path + '/TEmarker_output.py -d ' + D03_generate_output_wd + ' -o ' + D03_generate_output_od + \
              ' -genos_f ' + D02_genotyping_od + '/opt_modified_te_genotype.txt' + \
              ' -m_f ' + material_file + \
              ' -fas_f ' + genome_file
        print(cmd)
        subprocess.call(cmd,shell=True)

    ##pipeline finish
    print('All the steps are complete')
    print('copy output from Step 3 to the output_dir')
    cmd = 'cp -r ' + D03_generate_output_od + '/* ' + output_dir
    print(cmd)
    subprocess.call(cmd,shell=True)

if __name__ == "__main__":
    main()

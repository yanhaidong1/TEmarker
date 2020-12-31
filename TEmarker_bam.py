#!/usr/bin/env python

##this script helps to generate bam bai and remove PCR duplicate
##updation 042620 change the output name and remove the working_dir that is not used in this script
##updaiton 043020 add working_dir for storing sort bam

import argparse
import sys
from distutils.spawn import find_executable

import os
import subprocess
import glob
import re

##SCRIPTS
def get_parsed_args():
    parser = argparse.ArgumentParser(description="Genotype TE markers")

    ##############################
    ##parse the required arguments
    parser.add_argument("-o", dest='output_dir', default="./", help="Output directory"
                                                                    "Default: ./ ")

    parser.add_argument("-d", dest='working_dir', default="./", help="Working directory"
                                                                     "Default: ./")

    parser.add_argument('-bam_bwa_d', dest='bam_bwa_dir', help="Users provide the bam dir generated from the TEScape_marker_bed."
                                                               "This dir contains bam files generated from MicClinctock tool."
                                                               "Or uses can generate the bam dir based on bwa tools.")

    parser.add_argument('-bam_hisat2_d', dest='bam_hisat_dir', help="Users provide the bam dir generated from the TEScape_marker_bed."
                                                                   "This dir contains bam files generated from hisat2 tool")

    parser.add_argument('--samtools', dest='samtools', default="/usr/bin/samtools",help="Users should provide the path of samtools.")

    ##optional
    parser.add_argument('-pros_n_samtls', dest='process_num_samtools',help="Users provide process number that helps to divide the candidate TE insertion file"
                                                                           "into multiple files, which can increase the speed of processing."
                                                                           "Default: 1")

    parser.add_argument('-clean', dest='clean_temp_dir',help='If user set the yes, all the temp dir in the working_dir will be deletec to save the place')

    ##parse of parameters
    args = parser.parse_args()
    return args

def main(argv=None):
    if argv is None:
        argv = sys.argv
    args = get_parsed_args()

    #######################################
    ##check the required software and files
    bam_bwa_dir = args.bam_bwa_dir
    if args.bam_bwa_dir is None:
        print('Cannot find input te bam bwa dir, please provide that')
        return  ##import to close the script if not find the te lib

    bam_hisat_dir = args.bam_hisat_dir
    if args.bam_hisat_dir is None:
        print('Cannot find input te bam hisat dir, please provide that')
        return  ##import to close the script if not find the te lib

    ##################
    ##for the software
    samtools_exe = ''
    if find_executable(args.samtools) is not None:
        samtools_exe = args.samtools
        print ('samtools executable can be found')
    else:
        print("Cannot find samtools executable, please check if it has been installed.")
        return

    if args.process_num_samtools is not None:
        process_num_samtools = args.process_num_samtools
    else:
        process_num_samtools = 1

    ###########################################
    ##create the output directories
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

    #################
    ##run the process
    #################
    ##for the bwa bam
    mcc_bam_file_list = glob.glob(bam_bwa_dir + '/*.bam')
    for eachbamfl in mcc_bam_file_list:
        ##change the name of bam file if there is a _2 in the dir
        ##get the fl name  ##attention for this name information
        if '_2' in eachbamfl:
            mt = re.match('.+/(.+)_2\.bam', eachbamfl)
            fl_nm = mt.group(1)
        else:
            mt = re.match('.+/(.+)\.bam', eachbamfl)
            fl_nm = mt.group(1)

        cmd = 'mv ' + eachbamfl + ' ' + bam_bwa_dir + '/' + fl_nm + '.bam'
        print(cmd)
        subprocess.call(cmd, shell=True)


    #################
    ##updation 043020
    ##for the bwa bam
    bwa_bam_wd = working_dir + '/bwa_bam_wd'
    if not os.path.exists(bwa_bam_wd):
        os.makedirs(bwa_bam_wd)

    bwa_bam_sort_dir = bwa_bam_wd + '/bwa_bam_sort_dir'
    if not os.path.exists(bwa_bam_sort_dir):
        os.makedirs(bwa_bam_sort_dir)

    ##sort the bwa bam
    bwa_bam_file_list = glob.glob(bam_bwa_dir + '/*.bam')
    for eachbamfl in bwa_bam_file_list:

        if '_2' in eachbamfl:
            mt = re.match('.+/(.+)_2\.bam', eachbamfl)
            fl_nm = mt.group(1)
        else:
            mt = re.match('.+/(.+)\.bam', eachbamfl)
            fl_nm = mt.group(1)

        ##sort the bamfile  ##should use -o
        cmd = samtools_exe + ' sort ' + eachbamfl + ' -@ ' + str(process_num_samtools) + ' -o ' + bwa_bam_sort_dir + '/' + fl_nm + '.bam > ' \
              + bwa_bam_sort_dir + '/' + fl_nm + '.bam'
        print(cmd)
        subprocess.call(cmd, shell=True)

    ##rm the pcr bwa bam
    rm_pcr_bwa_bam_dir = bwa_bam_wd + '/rm_pcr_bwa_bam_dir'
    if not os.path.exists(rm_pcr_bwa_bam_dir):
        os.makedirs(rm_pcr_bwa_bam_dir)

    bwa_sort_bam_file_list = glob.glob(bwa_bam_sort_dir + '/*.bam')
    for eachbamfl in bwa_sort_bam_file_list:
        mt = re.match('.+/(.+)\.bam', eachbamfl)
        fl_nm = mt.group(1)
        ##remove the pcr duplicates
        cmd = samtools_exe + ' rmdup ' + eachbamfl + ' ' + rm_pcr_bwa_bam_dir + '/' + fl_nm + '.bam'
        print(cmd)
        subprocess.call(cmd, shell=True)

    rm_pcr_bwa_bam_list = glob.glob(rm_pcr_bwa_bam_dir + '/*.bam')
    for eachbamfl in rm_pcr_bwa_bam_list:
        ##use samtools view to generate the index
        cmd = samtools_exe + ' index ' + eachbamfl + ' ' + eachbamfl + '.bai'
        print(cmd)
        subprocess.call(cmd, shell=True)

    ################
    ##for the hisat2
    hisat2_bam_wd = working_dir + '/hisat2_bam_wd'
    if not os.path.exists(hisat2_bam_wd):
        os.makedirs(hisat2_bam_wd)

    hisat2_bam_sort_dir = hisat2_bam_wd + '/hisat2_bam_sort_dir'
    if not os.path.exists(hisat2_bam_sort_dir):
        os.makedirs(hisat2_bam_sort_dir)

    hisat2_bam_file_list = glob.glob(bam_hisat_dir + '/*.bam')
    for eachbamfl in hisat2_bam_file_list:
        ##get the fl name  ##attention for this name information
        # fl_nm = ''
        if '_2' in eachbamfl:
            mt = re.match('.+/(.+)_2\.bam', eachbamfl)
            fl_nm = mt.group(1)
        else:
            mt = re.match('.+/(.+)\.bam', eachbamfl)
            fl_nm = mt.group(1)

        ##sor the bamfile  ##should use -o
        cmd = samtools_exe + ' sort ' + eachbamfl + ' -@ ' + str(process_num_samtools) + ' -o ' + hisat2_bam_sort_dir + '/' + fl_nm + '.bam > ' \
              + hisat2_bam_sort_dir + '/' + fl_nm + '.bam'
        print(cmd)
        subprocess.call(cmd, shell=True)

    ##updation8.27 add pcr rm analysis
    rm_pcr_hisat2_bam_dir = hisat2_bam_wd + '/rm_pcr_hisat2_bam_dir'
    if not os.path.exists(rm_pcr_hisat2_bam_dir):
        os.makedirs(rm_pcr_hisat2_bam_dir)

    hisat2_sort_bam_file_list = glob.glob(hisat2_bam_sort_dir + '/*.bam')
    for eachbamfl in hisat2_sort_bam_file_list:
        mt = re.match('.+/(.+)\.bam', eachbamfl)
        fl_nm = mt.group(1)
        ##remove the pcr duplicates
        cmd = samtools_exe + ' rmdup ' + eachbamfl + ' ' + rm_pcr_hisat2_bam_dir + '/' + fl_nm + '.bam'
        print(cmd)
        subprocess.call(cmd, shell=True)

    ##updation8.27 get the new index of the bam file
    ##updation4.19: get the index of the bam file
    hisat2_rm_pcr_bam_file_list = glob.glob(rm_pcr_hisat2_bam_dir + '/*.bam')
    # index the bam file
    for eachbamfl in hisat2_rm_pcr_bam_file_list:
        cmd = samtools_exe + ' index ' + eachbamfl + ' ' + eachbamfl + '.bai'
        print(cmd)
        subprocess.call(cmd, shell=True)

    ##move the rm_pcr_hisat2_bam_dir to the output dir
    ##mv mcc dir
    cmd = 'mv ' + rm_pcr_bwa_bam_dir + ' ' + output_dir
    print(cmd)
    subprocess.call(cmd, shell=True)

    ##mv hisat2 dir
    cmd = 'mv ' + rm_pcr_hisat2_bam_dir + ' ' + output_dir
    print(cmd)
    subprocess.call(cmd, shell=True)

    ##clean the temp file
    ##now the temp dir is the sorted hisat2 bam dir
    if args.clean_temp_dir is not None:
        if args.clean_temp_dir == 'yes':

            ##check if bed_dir and two bam files exit
            cmd = 'rm -rf ' + bwa_bam_wd + ' ' + hisat2_bam_wd
            print(cmd)
            subprocess.call(cmd,shell=True)

        else:
            print("please type 'yes' behind clean argument.")
            return

if __name__ == "__main__":
    main()
#!/usr/bin/env python

##updating 082221 add a function to do the single reads mapping
##this script will wrap bwa hisat2 and cnvnator as a pipeline to help users to develop the relative files using this script

import argparse
import sys
from distutils.spawn import find_executable

import os
import subprocess
import glob
import re
##SCRIPTS

def get_parsed_args():
    parser = argparse.ArgumentParser(description="Prepare input files for Step 1, 2, 3")

    ##############################
    ##parse the required arguments
    parser.add_argument("-d", dest='working_dir', help="Working directory"
                                                        "Default: ./ ")

    parser.add_argument("-o", dest='output_dir',  help="Output directory"
                                                        "Default: ./ ")

    parser.add_argument('-fas_f',dest='genome_fas',help="Provide the genome fasta file.")

    ########################
    ##optional for each tool
    ################
    ##for the hisat2
    parser.add_argument('-hst', dest='hisat2_yes', help="If users set hisat2 to yes, users need to provide the required hisat2 dependences."
                                                        "Please provide samtools path using --samtools and fastq_dir using -fastq_d")

    parser.add_argument('--hisat2_t', dest='hisat2', default="/usr/bin/hisat2",help="Users should provide the path of hisat2.")

    parser.add_argument('--hisat2build', dest='hisat2build', default="/usr/bin/hisat2-build",help="Users should provide the path of hisat2-build.")

    parser.add_argument('-pros_n_hisat2', dest='process_num_hisat2',help="Users provide process number that helps to increase speed in hisat2."
                                                                           "Default: 1")

    parser.add_argument('--samtools', dest='samtools', default="/usr/bin/samtools",help="Users should provide the path of samtools.")

    parser.add_argument('-fastq_d', dest='fastq_dir', help="Provide a dir to contain all the fastq file."
                                                           "Attention: please provide the absolute path of fastq_dir.")

    #############
    ##for the bwa
    parser.add_argument('-bwa', dest='bwa_yes', help= 'If user set bwa to yes, users need to provide the required bwa dependences.'
                                                      'Please provide samtools path using --samtools and fastq_dir using -fastq_d')

    parser.add_argument('--bwa_t', dest='bwa', default='/usr/bin/bwa',help='Users should provide the path of bwa.')

    parser.add_argument('-pros_n_bwa', dest='process_num_bwa', help= 'Users provide process number that helps to increase speed in bwa.'
                                                                     'Default: 1')

    ##updating 082221
    parser.add_argument('-bwa_s', dest='bwa_s_yes', help='If user set bwa_s to yes, users need to provide the fastq_d that contains single read.'
                                                      'Please provide samtools path using --samtools and fastq_dir using -fastq_d')

    ##bwa also needs samtools and fastq_d

    ##################
    ##for the cnvnator
    parser.add_argument('-cnv', dest='cnvnator_yes', help='If users set cnvnator to yes, users need to provide the required cnvnator dependences')

    parser.add_argument('-cnvnator_t', dest='cnvnator_tool', default='/usr/bin/cnvnator',help="Users provide the cnvnator tool")

    parser.add_argument('-bam_d', dest='bam_dir',help="Users provide bam dir generated from the TEmarker_bed.py,"
                                                      "Or Users provide bam generated from bwa tool using TEmarker_tool.py.")

    ##############################
    ##optional for all three tools
    parser.add_argument('-m_f', dest='name_file', help="If users provide the name_file, and pipeline will substitue the fastq name "
                                                       "to the name of the second column in the name_file.")

    parser.add_argument('-clean', dest='clean_temp_dir',help='If user set the yes, all the temp dir will be deleted to save the place')

    ##parse of parameters
    args = parser.parse_args()
    return args

def main(argv=None):
    if argv is None:
        argv = sys.argv
    args = get_parsed_args()

    ###########################################
    ##create the working and output directories
    working_dir = args.working_dir
    if not working_dir.endswith('/'):
        working_dir = working_dir + '/'
    else:
        working_dir = working_dir

    output_dir = args.output_dir
    if not output_dir.endswith('/'):
        output_dir = output_dir + '/'
    else:
        output_dir = output_dir

    #####################
    ##if hisat2 is called
    if args.hisat2_yes is not None:

        if args.hisat2_yes != 'yes':
            print("please type 'yes' behind -hst to call the hisat2 pipeline.")
            return
        else:
            #################################################
            ##the following is to conduct the hisat2 pipeline
            #######################
            ##check depedence tools
            if find_executable(args.samtools) is not None:
                samtools_exe = args.samtools
                print ('samtools executable can be found')
            else:
                print("Cannot find samtools executable, please check if it has been installed.")
                return

            if args.process_num_hisat2 is not None:
                process_num_hisat2 = args.process_num_hisat2
            else:
                process_num_hisat2 = 1

            if find_executable(args.hisat2) is not None:
                hisat2_exe = args.hisat2
                print ('hisat2 executable can be found')
            else:
                print("Cannot find hisat2 executable, please check if it has been installed.")
                return

            if find_executable(args.hisat2build) is not None:
                hisat2build_exe = args.hisat2build
                print ('hisat2-build executable can be found')
            else:
                print("Cannot find hisat2-build executable, please check if it has been installed.")
                return

            #######################
            ##check the input files
            if args.genome_fas is None:
                print ('Cannot find input genome fasta file, please provide that')
                return
            else:
                try:
                    file = open(args.genome_fas, 'r')
                except IOError:
                    print('There was an error opening the genome fasta file!')
                    return

            if args.fastq_dir is None:
                print ('Cannot find fastq dir, please provide it')
                return

            ############################################
            ##generate hisat2 working_dir and output_dir
            hisat2_working_dir = working_dir + '/hisat2_working_dir'
            if not os.path.exists(hisat2_working_dir):
                os.makedirs(hisat2_working_dir)

            hisat2_output_dir = output_dir + '/hisat2_output_dir'
            if not os.path.exists(hisat2_output_dir):
                os.makedirs(hisat2_output_dir)

            ###################################################################
            ##use hisat2 to generate another bam file to detect deletion events
            ###################################################################
            ##create a hisat2 dir
            hisat2_dir = hisat2_working_dir + '/hisat2_dir'
            if not os.path.exists(hisat2_dir):
                os.makedirs(hisat2_dir)

            ##create a sam dir to store all the sample output from the hisat2
            hisat2_sam_dir = hisat2_dir + '/hisat2_sam_dir'
            if not os.path.exists(hisat2_sam_dir):
                os.makedirs(hisat2_sam_dir)

            ##create a bam dir to store the bam output from the sam output
            hisat2_bam_dir = hisat2_output_dir + '/hisat2_bam_dir'
            if not os.path.exists(hisat2_bam_dir):
                os.makedirs(hisat2_bam_dir)

            ##prepare genome
            genome_file = args.genome_fas

            if '/' in genome_file:
                mt = re.match('.+/(.+)', genome_file)
                genome_fl_nm = mt.group(1)
                mt = re.match('(.+)\.(.+)', genome_fl_nm)
                genome_nm = mt.group(1)
                genome_suffix = mt.group(2)
            else:
                genome_fl_nm = genome_file
                mt = re.match('(.+)\.(.+)', genome_fl_nm)
                genome_nm = mt.group(1)
                genome_suffix = mt.group(2)

            ##cp the fasta file to the hisat2 dir
            cmd = 'cp ' + genome_file + ' ' + hisat2_dir
            print(cmd)
            subprocess.call(cmd, shell=True)

            ##prepare fastq dir
            sample_nm_dic = {}
            if args.name_file is not None:
                sample_nm_file = args.name_file
                ##initiate a dic to store the name into a dic
                ##key is the orginial name and value is the new name
                with open (sample_nm_file,'r') as ipt:
                    for eachline in ipt:
                        eachline = eachline.strip('\n')
                        col = eachline.strip().split()
                        sample_nm_dic[col[0]] = col[1]
            else:
                ##search name information in the fastq_dir
                fastq_dir = args.fastq_dir
                fastq_list = glob.glob(fastq_dir + '/*')
                for eachfastq in fastq_list:
                    ##extract the name from the eachfastq
                    mt = re.match('(.+)/(.+)_(\d)(\..+)', eachfastq)
                    fq_nm = mt.group(2)
                    sample_nm_dic[fq_nm] = fq_nm
            ##prepare fastq files
            fastq_dir = args.fastq_dir
            fastq_list = glob.glob(fastq_dir + '/*')

            ##create the indext for the genome
            cmd = hisat2build_exe + ' ' + hisat2_dir + '/' + genome_nm + '.' + genome_suffix + ' ' + hisat2_dir + '/' + genome_nm + '.' + genome_suffix
            print(cmd)
            subprocess.call(cmd, shell=True)

            for eachnm in sample_nm_dic:
                new_nm = sample_nm_dic[eachnm]

                mt_pair_list = []
                ##updation 11.23 revise the input fastq
                for eachfastq in fastq_list:
                    ##extract the name from the eachfastq
                    mt = re.match('(.+)/(.+)_(\d)(\..+)', eachfastq)
                    path = mt.group(1)
                    fq_nm = mt.group(2)
                    pair_num = mt.group(3)
                    type_nm = mt.group(4)

                    ##updation 11.23
                    ##new_nm should be changed to fq_nm
                    if eachnm == fq_nm:
                        mt_pair_list.append(path + '/' + fq_nm + '_' + pair_num + type_nm)

                ##run the hisat2
                cmd = hisat2_exe + ' -p ' + str(process_num_hisat2) + ' -x ' +  hisat2_dir + '/' + genome_nm + '.' + genome_suffix + ' -1 ' + \
                      mt_pair_list[0] + ' -2 ' + mt_pair_list[1] + ' -S ' + hisat2_sam_dir + '/' + new_nm + '.sam'
                print(cmd)
                subprocess.call(cmd, shell=True)

                ##use the samtools to change the sam to the bam
                cmd = samtools_exe + ' view -bS ' +  hisat2_sam_dir + '/' + new_nm + '.sam' + ' > ' + hisat2_bam_dir + '/' + \
                      new_nm + '.bam'
                print(cmd)
                subprocess.call(cmd, shell=True)

                if args.clean_temp_dir is not None:
                    if args.clean_temp_dir == 'yes':
                        ##delete the sam file to save the storage
                        cmd = 'rm ' + hisat2_sam_dir + '/' + new_nm + '.sam'
                        print(cmd)
                        subprocess.call(cmd, shell=True)
                    else:
                        print("please type 'yes' behind clean argument.")
                        return

    ##################
    ##if bwa is called
    if args.bwa_yes is not None:

        if args.bwa_yes != 'yes':
            print("please type 'yes' behind -bwa to call the bwa pipeline.")
            return
        else:
            #######################
            ##check depedence tools
            if find_executable(args.samtools) is not None:
                samtools_exe = args.samtools
                print('samtools executable can be found')
            else:
                print("Cannot find samtools executable, please check if it has been installed.")
                return

            if args.process_num_bwa is not None:
                process_num_bwa = args.process_num_bwa
            else:
                process_num_bwa = 1

            if find_executable(args.bwa) is not None:
                bwa_exe = args.bwa
                print('bwa executable can be found')
            else:
                print("Cannot find bwa executable, please check if it has been installed.")
                return

            #######################
            ##check the input files
            if args.genome_fas is None:
                print('Cannot find input genome fasta file, please provide that')
                return
            else:
                try:
                    file = open(args.genome_fas, 'r')
                except IOError:
                    print('There was an error opening the genome fasta file!')
                    return

            if args.fastq_dir is None:
                print('Cannot find fastq dir, please provide it')
                return

            ############################################
            ##generate hisat2 working_dir and output_dir
            bwa_working_dir = working_dir + '/bwa_working_dir'
            if not os.path.exists(bwa_working_dir):
                os.makedirs(bwa_working_dir)

            bwa_output_dir = output_dir + '/bwa_output_dir'
            if not os.path.exists(bwa_output_dir):
                os.makedirs(bwa_output_dir)

            ########################
            ##start the bwa pipeline
            ##create a bwa dir
            bwa_dir = bwa_working_dir + '/bwa_dir'
            if not os.path.exists(bwa_dir):
                os.makedirs(bwa_dir)

            ##create a sam dir to store all the sample output from the bwa
            bwa_sam_dir = bwa_dir + '/bwa_sam_dir'
            if not os.path.exists(bwa_sam_dir):
                os.makedirs(bwa_sam_dir)

            ##create a bam dir to store the bam output from the sam output
            bwa_bam_dir = bwa_output_dir + '/bwa_bam_dir'
            if not os.path.exists(bwa_bam_dir):
                os.makedirs(bwa_bam_dir)

            ##prepare genome
            genome_file = args.genome_fas

            if '/' in genome_file:
                mt = re.match('.+/(.+)', genome_file)
                genome_fl_nm = mt.group(1)
                mt = re.match('(.+)\.(.+)', genome_fl_nm)
                genome_nm = mt.group(1)
                genome_suffix = mt.group(2)
            else:
                genome_fl_nm = genome_file
                mt = re.match('(.+)\.(.+)', genome_fl_nm)
                genome_nm = mt.group(1)
                genome_suffix = mt.group(2)

            ##cp the fasta file to the hisat2 dir
            cmd = 'cp ' + genome_file + ' ' + bwa_dir
            print(cmd)
            subprocess.call(cmd, shell=True)


            ##updating 082221
            if args.bwa_s_yes != 'yes':

                ##prepare fastq dir
                sample_nm_dic = {}
                if args.name_file is not None:
                    sample_nm_file = args.name_file
                    ##initiate a dic to store the name into a dic
                    ##key is the orginial name and value is the new name
                    with open(sample_nm_file, 'r') as ipt:
                        for eachline in ipt:
                            eachline = eachline.strip('\n')
                            col = eachline.strip().split()
                            sample_nm_dic[col[0]] = col[1]
                else:
                    ##search name information in the fastq_dir
                    fastq_dir = args.fastq_dir
                    fastq_list = glob.glob(fastq_dir + '/*')
                    for eachfastq in fastq_list:
                        ##extract the name from the eachfastq
                        mt = re.match('(.+)/(.+)_(\d)(\..+)', eachfastq)
                        fq_nm = mt.group(2)
                        sample_nm_dic[fq_nm] = fq_nm

                ##prepare fastq files
                fastq_dir = args.fastq_dir
                fastq_list = glob.glob(fastq_dir + '/*')

                ##create index of genome
                ##bwa will directly use genome_nm as the index name
                cmd = bwa_exe + ' index ' + bwa_dir + '/' + genome_nm + '.' + genome_suffix
                print(cmd)
                subprocess.call(cmd, shell=True)

                for eachnm in sample_nm_dic:
                    new_nm = sample_nm_dic[eachnm]

                    mt_pair_list = []
                    ##updation 11.23 revise the input fastq
                    for eachfastq in fastq_list:
                        ##extract the name from the eachfastq
                        mt = re.match('(.+)/(.+)_(\d)(\..+)', eachfastq)
                        path = mt.group(1)
                        fq_nm = mt.group(2)
                        pair_num = mt.group(3)
                        type_nm = mt.group(4)

                        ##updation 11.23
                        ##new_nm should be changed to fq_nm
                        if eachnm == fq_nm:
                            mt_pair_list.append(path + '/' + fq_nm + '_' + pair_num + type_nm)

                    ##run the bwa
                    cmd = bwa_exe + ' mem -t ' + str(process_num_bwa) + ' -v 0 -R \'@RG\\tID:\'' + new_nm + '\'\\tSM:\'' + \
                          new_nm + ' ' + bwa_dir + '/' + genome_nm + '.' + genome_suffix + ' ' + mt_pair_list[0] + ' ' + mt_pair_list[1] + ' > ' + \
                          bwa_sam_dir + '/' + new_nm + '.sam'
                    print(cmd)
                    subprocess.call(cmd, shell=True)

                    ##use the samtools to change the sam to the bam
                    cmd = samtools_exe + ' view -bS ' +  bwa_sam_dir + '/' + new_nm + '.sam' + ' > ' + bwa_bam_dir + '/' + \
                          new_nm + '.bam'
                    print(cmd)
                    subprocess.call(cmd, shell=True)

                    if args.clean_temp_dir is not None:
                        if args.clean_temp_dir == 'yes':
                            ##delete the sam file to save the storage
                            cmd = 'rm ' + bwa_sam_dir + '/' + new_nm + '.sam'
                            print(cmd)
                            subprocess.call(cmd, shell=True)
                        else:
                            print("please type 'yes' behind clean argument.")
                            return

            else:

                ##prepare fastq dir
                sample_nm_dic = {}
                if args.name_file is not None:
                    sample_nm_file = args.name_file
                    ##initiate a dic to store the name into a dic
                    ##key is the orginial name and value is the new name
                    with open(sample_nm_file, 'r') as ipt:
                        for eachline in ipt:
                            eachline = eachline.strip('\n')
                            col = eachline.strip().split()
                            sample_nm_dic[col[0]] = col[1]
                else:
                    ##search name information in the fastq_dir
                    fastq_dir = args.fastq_dir
                    fastq_list = glob.glob(fastq_dir + '/*')
                    for eachfastq in fastq_list:
                        ##extract the name from the eachfastq
                        mt = re.match('.+/(.+)',eachfastq)
                        flnm = mt.group(1)

                        fq_nm = ''
                        if 'fastq' in flnm:
                            mt = re.match('(.+)\.fast.+',flnm)
                            fq_nm = mt.group(1)
                        if 'fq' in flnm:
                            if 'gz' in flnm:
                                mt = re.match('(.+)\.fq.gz',flnm)
                                fq_nm = mt.group(1)
                            else:
                                mt = re.match('(.+)\.fq',flnm)
                                fq_nm = mt.group(1)

                        sample_nm_dic[fq_nm] = fq_nm

                    ##prepare fastq files
                    fastq_dir = args.fastq_dir
                    fastq_list = glob.glob(fastq_dir + '/*')

                    ##create index of genome
                    ##bwa will directly use genome_nm as the index name
                    cmd = bwa_exe + ' index ' + bwa_dir + '/' + genome_nm + '.' + genome_suffix
                    print(cmd)
                    subprocess.call(cmd, shell=True)

                    for eachnm in sample_nm_dic:
                        new_nm = sample_nm_dic[eachnm]

                        mt_pair_list = []
                        ##updation 11.23 revise the input fastq
                        for eachfastq in fastq_list:
                            ##extract the name from the eachfastq
                            #mt = re.match('(.+)/(.+)_(\d)(\..+)', eachfastq)
                            #path = mt.group(1)
                            #fq_nm = mt.group(2)
                            #pair_num = mt.group(3)
                            #type_nm = mt.group(4)


                            mt = re.match('(.+)/(.+)',eachfastq)
                            path = mt.group(1)
                            flnm = mt.group(2)

                            type_nm = ''
                            fq_nm = ''
                            if 'fastq' in flnm:
                                mt = re.match('(.+)\.(fast.+)', flnm)
                                fq_nm = mt.group(1)
                                type_nm = mt.group(2)
                            if 'fq' in flnm:
                                if 'gz' in flnm:
                                    mt = re.match('(.+)\.(fq.gz)', flnm)
                                    fq_nm = mt.group(1)
                                    type_nm = mt.group(2)
                                else:
                                    mt = re.match('(.+)\.(fq)', flnm)
                                    fq_nm = mt.group(1)
                                    type_nm = mt.group(2)

                            ##updation 11.23
                            ##new_nm should be changed to fq_nm
                            if eachnm == fq_nm:
                                mt_pair_list.append(path + '/' + fq_nm + '.' + type_nm)

                        ##run the bwa
                        cmd = bwa_exe + ' mem -t ' + str(
                            process_num_bwa) + ' -v 0 -R \'@RG\\tID:\'' + new_nm + '\'\\tSM:\'' + \
                              new_nm + ' ' + bwa_dir + '/' + genome_nm + '.' + genome_suffix + ' ' + mt_pair_list[
                                  0] + ' > ' + \
                              bwa_sam_dir + '/' + new_nm + '.sam'
                        print(cmd)
                        subprocess.call(cmd, shell=True)

                        ##use the samtools to change the sam to the bam
                        cmd = samtools_exe + ' view -bS ' + bwa_sam_dir + '/' + new_nm + '.sam' + ' > ' + bwa_bam_dir + '/' + \
                              new_nm + '.bam'
                        print(cmd)
                        subprocess.call(cmd, shell=True)

                        if args.clean_temp_dir is not None:
                            if args.clean_temp_dir == 'yes':
                                ##delete the sam file to save the storage
                                cmd = 'rm ' + bwa_sam_dir + '/' + new_nm + '.sam'
                                print(cmd)
                                subprocess.call(cmd, shell=True)
                            else:
                                print("please type 'yes' behind clean argument.")
                                return

    #######################
    ##if cnvnator is called
    if args.cnvnator_yes is not None:

        if args.cnvnator_yes != 'yes':
            print("please type 'yes' behind -cnv to call the cnvnator pipeline.")
            return
        else:

            ##check files
            genome_fasta_file = args.genome_fas
            if args.genome_fas is None:
                print('Cannot find input genome fasta file, please provide that')
                return  ##import to close the script if not find the te lib
            else:
                try:
                    file = open(args.genome_fasta_file, 'r')  ##check if the file is not the right file
                except IOError:
                    print('There was an error opening the genome fasta file!')
                    return

            bam_dir = args.bam_dir
            if args.bam_dir is None:
                print('Cannot find input te bam dir, please provide that')
                return  ##import to close the script if not find the te lib

            ##checking tools
            if find_executable(args.cnvnator_tool) is not None:
                cnvnator_tool = args.cnvnator_tool
                print('cnvnator_tool executable can be found')
            else:
                print("Cannot find cnvnator_tool executable, please check if it has been installed.")
                return

            ##create dir for the cnvnator
            cnvnator_working_dir = working_dir + '/cnvnator_working_dir'
            if not os.path.exists(cnvnator_working_dir):
                os.makedirs(cnvnator_working_dir)

            cnvnator_output_dir = output_dir + '/cnvnator_output_dir'
            if not os.path.exists(cnvnator_output_dir):
                os.makedirs(cnvnator_output_dir)

            ##run the processes
            bam_file_list = glob.glob(bam_dir + '/*.bam')

            for eachfl in bam_file_list:

                ##updation11.31
                if '_2' in eachfl:
                    mt = re.match('.+/(.+)_2\.bam', eachfl)
                    bam_nm = mt.group(1)
                else:
                    mt = re.match('.+/(.+)\.bam', eachfl)
                    bam_nm = mt.group(1)

                ##create dir in the input_working_dir to store the temp output file from cnvnator tool
                sample_bam_dir = cnvnator_working_dir + '/' + bam_nm
                if not os.path.exists(sample_bam_dir):
                    os.makedirs(sample_bam_dir)

                ##run the cnvnator
                ##tree
                cmd = cnvnator_tool + ' -root ' + sample_bam_dir + '/' + bam_nm + '.root ' + \
                      '-tree ' + eachfl
                print(cmd)
                subprocess.call(cmd, shell=True)

                ##his
                cmd = cnvnator_tool + ' -root ' + sample_bam_dir + '/' + bam_nm + '.root ' + \
                      '-his 100 -d ' + genome_fasta_file
                print(cmd)
                subprocess.call(cmd, shell=True)

                ##stat
                cmd = cnvnator_tool + ' -root ' + sample_bam_dir + '/' + bam_nm + '.root ' + \
                      '-stat 100'
                print(cmd)
                subprocess.call(cmd, shell=True)

                ##partition
                cmd = cnvnator_tool + ' -root ' + sample_bam_dir + '/' + bam_nm + '.root ' + \
                      '-partition 100'
                print(cmd)
                subprocess.call(cmd, shell=True)

                ##call cnv
                cmd = cnvnator_tool + ' -root ' + sample_bam_dir + '/' + bam_nm + '.root ' + \
                      '-call 100 > ' + sample_bam_dir + '/cnv_temp.txt'
                print(cmd)
                subprocess.call(cmd, shell=True)

                ##reanalysis his for 1000 bin to avoid warning
                cmd = cnvnator_tool + ' -root ' + sample_bam_dir + '/' + bam_nm + '.root ' + \
                      '-his 1000 -d ' + genome_fasta_file
                print(cmd)
                subprocess.call(cmd, shell=True)

                ##generate genotype file
                cmd = 'awk \'{ print $2 } END { print \"exit\" }\' ' + sample_bam_dir + '/cnv_temp.txt' + ' | ' + \
                      cnvnator_tool + ' -root ' + sample_bam_dir + '/' + bam_nm + '.root ' + \
                      '-genotype 100 > ' + sample_bam_dir + '/opt_genotype.txt'
                print(cmd)
                subprocess.call(cmd, shell=True)

                ##cp the genotype file to the output_dir
                cmd = 'cp ' + sample_bam_dir + '/opt_genotype.txt ' + cnvnator_output_dir + '/' + bam_nm + '_genotype.txt'
                print(cmd)
                subprocess.call(cmd, shell=True)


    if args.hisat2_yes is None and args.cnvnator_yes is None and args.bwa_yes is None:
        print('Please at least set one tool that will be used')
        return


if __name__ == "__main__":
    main()



















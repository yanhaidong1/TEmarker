#!/usr/bin/env python

##updating 052422 make annotations
##updating 052122 check directory of genome file
##updating 051922 close the repeatmasker
##updation 122420 modify a little for the argument description
##updation 100320 do not import others script since python2 cannot recognize these scripts
##conda install -c conda-forge biopython
##updation 100320 directly add another arguments for jitterbug and TEFLoN
##updation 100320 the mcclintock has not added the --mem and --make_annotations arguments we will directly use it
##updation later we can add them into the scripts
##updation 100220 use the new mcclintock pipeline and add other tools remove the histat2 and hisat2build
##updation 100220 since the jitterbug and TEFloN needs the python2 we need to generate another pipline

##updation 042120 mv the dir we needed to the output_dir
##updation 042120 change the fasta genome to fa or fasta
##updation 11.23 remove input_fastq_dir since it doest not work and it allows no bam file generated from hisat2
##updation 11.23 use : in the tool selection argument
##updation 11.15 set the file name as optional argument
##updation 11.15 set a tool selection argument
##default is RetroSeq;TEMP;TE-locate;ngs_te_mapper
##updation 11.15 cp the results to the output dir
##updation 11.15 set a argument to delete the temp dir to save space.
##updation 10.28 allow gzip fastq file as input file
##updation 9.30 delete PoPoolationTE
##updation 9.23 add another two tools in the pipeline PoPoolationTE and ngs_te_mapper

##this wrapper should have a argument to notify whether intermediate file will be stored or not
##updation 4.3 do not use new name in the script
##if there is gz, it should unzip first, since the hisat2 only accept the unzip fastq
##the input fastq should be _1.fq or _1.fastq do not _R1.fq or _R1.fastq
##the input fastq should be uncompressed or compressed (gz)

##updation 4.1 wrap the retroseq to the script
##updation 3.31
##change the reference file in the mcclintock that allows to run the RetroSeq

##This wrapper will helps users to generate the vcf file for the population genomic analysis
##BUILT-IN MODULES
import argparse
import sys
import os
import subprocess
import re
import glob
from Bio import SeqIO
from distutils.spawn import find_executable

##SCRIPTS
def get_parsed_args():
    parser = argparse.ArgumentParser(description="Prepare input files for Step 1, 2, 3")

    ##############################
    ##parse the required arguments
    parser.add_argument("-d", dest='working_dir', default="./", help="Working directory to store the output from the micclintock,"
                                                                     "and from the Hisat2 by the pipline recommanded software. "
                                                                     "Attention: please provide the absolute path of working_dir")

    parser.add_argument("-o", dest='output_dir', default="./", help="output directory to store the final insertion TE bed from software. "
                                                                    "Attention: please provide the absolute path of output_dir")

    parser.add_argument('--mcclintock', dest='mcclintock', help="Before using the mcclintock, please make sure mcclintock software can work."
                                                                "Attention: please provide the absolute path of mcclintock.")

    parser.add_argument('-fastq_d',dest='fastq_dir',help="Provide a dir to contain all the fastq file."
                                                          "Attention: please provide the absolute path of fastq_dir.")

    parser.add_argument('-fas_f',dest='genome_fas',help="Provide the genome fasta file."
                                                      "Attention: please provide the absolute path of genome_fas.")

    parser.add_argument('-lib_f', dest='te_library', help="Provide the TE library."
                                                        "Attention: please provide the absolute path of te_library.")

    ##updation 100320
    parser.add_argument('--TEFLoN_d', dest='TEFLoN', help="Specify a TEFLoN folder and use python2 conda environment to run")
    parser.add_argument('--jitterbug_d', dest='jitterbug', help="Provide a jitterbug folder and Use python2 conda environment to run")
    parser.add_argument('--repeatmasker', dest='repeatmasker',default="/usr/bin/RepeatMasker", help="jitterbug and TEFLoN need TE location files generated from the repeatmasker")
    parser.add_argument('-bam_bwa_d',dest='bam_bwa_dir',help='TEmarker_tool.py can help to generate bam files used for the jitterbug')
    parser.add_argument('--bwa', dest='bwa',default="bwa",help="TEFLoN needs to create the bam. bwa is provided to generated the bam file")
    parser.add_argument('--samtools',dest='samtools',default="samtools",help='TEFLoN needs the samtools to transfer the sam to the bam file')



    ##optional
    parser.add_argument('-m_f', dest="name_file",help="If users provide the name_file, and pipeline will substitue the fastq name to "
                                                          "the name of the second column in the name_file.")

    ##updation 11.15
    parser.add_argument('-tool', dest='tool_string',help="If users want to choose the tool they want to use in the Mcclintock"
                                                         "Default: all the tools")
    ##updation 11.15
    parser.add_argument('-clean', dest='clean_files',help='If user set the clean_yes, the large size files in the working_dir will be deleted to save the place')

    parser.add_argument('-pros_n', dest='processor',help="Provide the processor,and default will be 12")

    #parser.add_argument('-M', dest='mem',help='Provide the RAM, and default will be 16')

    ##parse of parameters
    args = parser.parse_args()
    return args

def change_TEnm_lib (input_TE_lib_fl):

    store_seq_dic = {}
    for seq_record in SeqIO.parse(input_TE_lib_fl,'fasta'):
        seq_nm = seq_record.id

        fam_col = seq_nm.split('_')
        class_nm = fam_col[0]

        new_name = seq_nm + '#' + class_nm
        store_seq_dic[new_name] = str(seq_record.seq)

    return (store_seq_dic)


def change_gff_to_gff3 (input_gff_from_rpmsk_fl):

    store_final_line_list = []
    with open (input_gff_from_rpmsk_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            if not eachline.startswith('#'):

                col = eachline.strip().split()

                mt = re.match('\"Motif:(.+)\"',col[9])
                TE_nm = mt.group(1)

                mt = re.match('(.+)_TE\d+$',TE_nm)
                TE_fam_nm = mt.group(1)

                annot_line = 'ID=' + TE_nm + ';Name=' + TE_nm + ';Family=' + TE_fam_nm
                final_line = col[0] + '\t' + col[1] + '\t' + col[2] + '\t' + col[3] + '\t' + col[4] + '\t' + \
                             col[5] + '\t' + col[6] + '\t' + col[7] + '\t' + annot_line
                store_final_line_list.append(final_line)

    return (store_final_line_list)

def change_format_jitterbug (jitterbug_gff_fl):

    ##for the jitterbug
    jitterbug_fl_path = jitterbug_gff_fl

    store_final_line_list = []
    with open (jitterbug_fl_path,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split('\t')

            annot_str = col[8]

            annot_col = annot_str.split('; ')

            store_annot_dic = {}
            for eachitem in annot_col:
                if '=' in eachitem:
                    mt = re.match('(.+)=(.+)',eachitem)
                    store_annot_dic[mt.group(1)] = mt.group(2)

            insert_str = store_annot_dic['softclipped_pos']
            mt = re.match('\((.+), (.+)\)',insert_str)
            insert_st = mt.group(1)
            insert_ed = mt.group(2)

            zygosity_value = store_annot_dic['zygosity']
            family_nm = store_annot_dic['predicted_superfam']

            final_line = col[0] + '\t' + insert_st + '\t' + insert_ed + '\t' + family_nm + '\t' + zygosity_value
            store_final_line_list.append(final_line)

    return (store_final_line_list)


def change_format_TEFLoN (TEFLoN_fl):

    ##for the TEFLoN
    TEFLoN_fl_path = TEFLoN_fl

    store_final_line_list = []
    with open(TEFLoN_fl_path, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()

            if col[6] == '-':


                loc_st = col[1]
                loc_ed = col[2]

                if loc_st == '-':
                    new_loc_st = loc_ed
                else:
                    new_loc_st = loc_st

                if loc_ed == '-':
                    new_loc_ed = loc_st
                else:
                    new_loc_ed = loc_ed

                genotype_type = col[12]
                family_nm = col[3]

                final_line = col[0] + '\t' + new_loc_st + '\t' + new_loc_ed + '\t' + family_nm + '\t' + genotype_type
                store_final_line_list.append(final_line)

    return (store_final_line_list)


def main(argv=None):
    if argv is None:
        argv = sys.argv
    args = get_parsed_args()

    #######################################
    ##check the required software and files
    ##for the input files
    if args.te_library is None:
        print ('Cannot find input te lib, please provide that')
        return ##import to close the script if not find the te lib
    else:
        try:
            file = open(args.te_library, 'r')  ##check if the file is not the right file
        except IOError:
            print('There was an error opening the file!')
            return

    if args.genome_fas is None:
        print ('Cannot find input genome fasta file, please provide that')
        return
    else:
        try:
            file = open(args.genome_fas, 'r')
        except IOError:
            print('There was an error opening the file!')
            return

    if args.fastq_dir is None:
        print ('Cannot find fastq dir, please provide it')
        return

    ##################
    ##for the software
    if args.mcclintock is not None:
        if args.jitterbug is not None or args.TEFLoN is not None:
            print('mcclintock cannot run with the jitterbug and TEFLoN!')
            return

    if args.jitterbug is not None or args.TEFLoN is not None:
        if args.mcclintock is not None:
            print('mcclintock cannot run with the jitterbug and TEFLoN!')
            return
        #if args.repeatmasker is None:
        #    print('provide repeatmasker when jitterbug is activated')
        #    return

    if args.mcclintock is None and args.jitterbug is None and args.TEFLoN is None:
        print('Please provide the tools')
        return

    if args.jitterbug is not None:
        if args.bam_bwa_dir is None:
            print('provide bwa bam dir')
            return

    if find_executable(args.samtools) is not None:
        print ('samtools executable can be found')
    else:
        print("Cannot find samtools executable, please check if it has been installed.")
        return

    if find_executable(args.bwa) is not None:
        print ('bwa executable can be found')
    else:
        print("Cannot find bwa executable, please check if it has been installed.")
        return


    #if args.TEFLoN is not None:
    #    if args.bwa is None:
    #        print('provide bwa tool')
    #        return
    #    if args.samtools is None:
    #        print('provide samtools tool')
    #        return

    #######################
    ##for the optional file
    ##check for the name_file
    if args.name_file is not None:
        try:
            file = open(args.name_file, 'r')  ##check if the file is not the right file
        except IOError:
            print('There was an error opening the name_file!')
            return

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

    #################################
    ##create a dir in the working dir
    sample_bed_dir = output_dir + '/bed_dir'
    if not os.path.exists(sample_bed_dir):
        os.makedirs(sample_bed_dir)

    #############################################
    ##prepare the stuffs before running the tools
    ##extract the name information
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
        fastq_fl_list = glob.glob(fastq_dir + '/*')
        for eachfastq_fl in fastq_fl_list:
            mt = re.match('.+/(.+)',eachfastq_fl)
            fl_nm = mt.group(1)
            ##extract the name from the eachfastq
            mt = re.match('(.+)_(\d)(\..+)', fl_nm)
            fq_nm = mt.group(1)
            sample_nm_dic[fq_nm] = fq_nm

    ##read the fastq_dir
    fastq_dir = args.fastq_dir
    fastq_list = glob.glob(fastq_dir + '/*')
    genome_file = args.genome_fas
    library_file = args.te_library

    ##create the pro number
    if args.processor is not None:
        pro_num = args.processor
    else:
        pro_num = '1'

    #################
    ##mcclintock tool
    #################
    if args.mcclintock is not None:

        print('step')

        if args.tool_string is not None:
            ##extract tool information
            tool_str = args.tool_string
        else:
            tool_str = 'all'

        ##create a dir in the working_dir
        bwa_bam_dir = output_dir + '/bwa_bam_dir'
        if not os.path.exists(bwa_bam_dir):
            os.makedirs(bwa_bam_dir)

        #########################
        ##1. run the first sample
        #########################
        ##run the first sample to get the gff and hierachy file
        ##generate the gff and hie file
        ##run the first sample to remove the tes that are not available in the RetroSeq software
        ##generate another gff and hie file for the RetroSeq

        ##get the first name in the sample_nm_dic
        #fst_nm_ori = list(sample_nm_dic.keys())[0]

        ##updation 10.28 do not cp the raw data to the temp fastq dir to save space
        ##do not unzip fastq if it is wrapped
        ##initiate a list to store the first pair match
        #first_mt_pair_list = []
        ##find the fastq pair for the first sample
        #for eachfastq_fl in fastq_list:
            ##extract the name from the eachfastq
        #    mt = re.match('(.+)/(.+)', eachfastq_fl)
        #    path = mt.group(1)
        #    fl_nm = mt.group(2)
        #    mt = re.match('(.+)_(\d)(\..+)',fl_nm)
        #    fq_nm = mt.group(1)
        #    pair_num = mt.group(2)
        #    type_nm = mt.group(3)

            ##updation 11.21 change fst_nm_new to the fq_nm
        #    if fst_nm_ori == fq_nm:
                ##store the fist match pair
        #        first_mt_pair_list.append(path + '/' +  fq_nm + '_' + pair_num + type_nm)

        ########################################################
        ##now, we have first pair samples in the input_fastq_dir
        ##run the mcclintock
        ##get the ram and processor information
        ##check for the p and M
        ##create the dir to store the preliminary results
        output_mcclintock_dir = working_dir + '/output_mcclintock'
        if not os.path.exists(output_mcclintock_dir):
            os.makedirs(output_mcclintock_dir)

        ##create the dir for the new first name dir
        ##use another name to that, since all the files will be re-analyzed
        first_nm_dir = output_mcclintock_dir + '/annotation'
        if not os.path.exists(first_nm_dir):
            os.makedirs(first_nm_dir)

        ##updating 052422 only make the annotation file
        mcclintock_exe = args.mcclintock
        cmd = 'python ' + mcclintock_exe + \
              ' -r ' + genome_file + \
              ' -c ' + library_file + \
              ' -p ' + str(pro_num) + \
              ' -o ' + first_nm_dir + \
              ' --make_annotations'
        subprocess.call(cmd, shell=True)



        #############################################
        ##2. extract the gff and hierachy information
        #############################################
        ##udpation 042120
        ##get the genome file without .fasta

        if '/' in genome_file:
            mt = re.match('.+/(.+)', genome_file)
            genome_fl_nm = mt.group(1)
        else:
            genome_fl_nm = genome_file

        mt = re.match('(.+)\.(.+)', genome_fl_nm)
        genome_nm = mt.group(1)
        ##udpation get the genome suffix
        #genome_suffix = mt.group(2)

        gff_path = first_nm_dir + '/' + genome_nm + '/reference_te_locations/unaugmented_inrefTEs.gff'
        hie_path = first_nm_dir + '/' + genome_nm + '/te_taxonomy/unaugmented_taxonomy.tsv'


        ######################
        ##3. run all the files
        ######################
        samples_output_dir = output_mcclintock_dir + '/samples_output_dir'
        if not os.path.exists(samples_output_dir):
            os.makedirs(samples_output_dir)

        ##For the te-locate and temp
        ##find the fastq pair
        for eachnm in sample_nm_dic:

            ##get the new name of fastq
            new_nm = sample_nm_dic[eachnm]

            ##create the directories
            opt_mcc_temp_dir = samples_output_dir + '/' + new_nm
            if not os.path.exists(opt_mcc_temp_dir):
                os.makedirs(opt_mcc_temp_dir)

            ##initiate a list to store the pair list
            mt_pair_list = []
            ##cp the raw fastq file to a new fastq file with new name
            for eachfastq_fl in fastq_list:
                ##extract the name from the eachfastq
                mt = re.match('(.+)/(.+)', eachfastq_fl)
                path = mt.group(1)
                fl_nm = mt.group(2)
                mt = re.match('(.+)_(\d)(\..+)', fl_nm)
                fq_nm = mt.group(1)
                pair_num = mt.group(2)
                type_nm = mt.group(3)

                ##updation 11.23 change new_nm to fq_nm in the append
                if eachnm == fq_nm:

                    mt_pair_list.append(path + '/' + fq_nm + '_' + pair_num + type_nm)

            if tool_str == 'all':
                cmd = 'python ' + mcclintock_exe + \
                      ' -r ' + genome_file + \
                      ' -c ' + library_file + \
                      ' -1 ' + mt_pair_list[0] + \
                      ' -2 ' + mt_pair_list[1] + \
                      ' -p ' + str(pro_num) + \
                      ' -o ' + opt_mcc_temp_dir + \
                      ' -g ' + gff_path + \
                      ' -t ' + hie_path
                subprocess.call(cmd, shell=True)
            else:
                cmd = 'python ' + mcclintock_exe + \
                      ' -r ' + genome_file + \
                      ' -c ' + library_file + \
                      ' -1 ' + mt_pair_list[0] + \
                      ' -2 ' + mt_pair_list[1] + \
                      ' -p ' + str(pro_num) + \
                      ' -m ' + tool_str + \
                      ' -o ' + opt_mcc_temp_dir + \
                      ' -g ' + gff_path + \
                      ' -t ' + hie_path
                subprocess.call(cmd, shell=True)

        ##############################################################
        ##4. extract all the bed and bam results to one output bed dir
        ##store the method and species name that fail to generate the bed files
        store_failed_line_list = []

        ##extract all the bed information
        samples_results_dir_list = glob.glob(output_mcclintock_dir + '/samples_output_dir/*')
        for eachsample_dir in samples_results_dir_list:
            mt = re.match('.+/(.+)', eachsample_dir)
            sample_nm = mt.group(1)

            ##updating 052422 sometimes the sample nm is SRR800842_1
            subsample_nm = ''
            allfl_list = glob.glob(eachsample_dir + '/*')
            for eachfl in allfl_list:
                mt = re.match('.+/(.+)',eachfl)
                flnm = mt.group(1)
                if sample_nm in flnm:
                    subsample_nm = flnm


            method_res_dir_list = glob.glob(eachsample_dir + '/' + subsample_nm + '/results/*')
            for eachmethod_dir in method_res_dir_list:
                mt = re.match('.+/(.+)',eachmethod_dir)
                eachmethod_nm = mt.group(1)

                if eachmethod_nm != 'coverage' and eachmethod_nm != 'summary':

                    method_res_fl_list = glob.glob(eachmethod_dir + '/*')
                    right_bed_file_count = 0
                    for eachmethod_res_fl in method_res_fl_list:
                        mt = re.match('.+/(.+)', eachmethod_res_fl)
                        med_fl_nm = mt.group(1)
                        if '_nonredundant.bed' in med_fl_nm:
                            bed_file_path = eachmethod_res_fl
                            cmd = 'cp ' + bed_file_path + ' ' + sample_bed_dir
                            subprocess.call(cmd,shell=True)
                            right_bed_file_count += 1

                        ##clean the files
                        if args.clean_files is not None:
                            if args.clean_files == 'yes':
                                if 'unfiltered' in med_fl_nm:
                                    cmd = 'rm -rf ' + med_fl_nm
                                    subprocess.call(cmd, shell=True)
                            else:
                                print("please type 'yes' behind clean argument.")
                                return

                    if right_bed_file_count == 0:
                        failed_line = sample_nm + '\t' + eachmethod_nm
                        store_failed_line_list.append(failed_line)

                if eachmethod_nm == 'coverage':
                    ##clean the files
                    if args.clean_files is not None:
                        if args.clean_files == 'yes':
                            cmd = 'rm ' + eachmethod_dir + '/input/*.sam'
                            subprocess.call(cmd, shell=True)
                        else:
                            print("please type 'yes' behind clean argument.")
                            return

            bam_res_dir_list = glob.glob(eachsample_dir + '/' + subsample_nm + '/intermediate/mapped_reads/*')
            sorted_bam_count = 0
            for eachfl in bam_res_dir_list:
                mt = re.match('.+/(.+)', eachfl)
                fl_nm = mt.group(1)
                if 'sorted.bam' in fl_nm:
                    bam_fl_path = eachfl
                    cmd = 'cp ' + bam_fl_path + ' ' + bwa_bam_dir
                    subprocess.call(cmd,shell=True)
                    sorted_bam_count += 1

                ##clean the files
                if args.clean_files is not None:
                    if args.clean_files == 'yes':
                        if '.sam' in fl_nm:
                            sam_fl_path = eachfl
                            cmd = 'rm ' + sam_fl_path
                            subprocess.call(cmd,shell=True)
                    else:
                        print("please type 'yes' behind clean argument.")
                        return

            ##store the wrong line
            if sorted_bam_count == 0:
                failed_line = sample_nm + '\t' + 'miss_bamfile'
                store_failed_line_list.append(failed_line)

            ##clean the the fastq
            if args.clean_files is not None:
                if args.clean_files == 'yes':
                    cmd = 'rm -rf ' + eachsample_dir + '/' + sample_nm + '/intermediate/fastq'
                    subprocess.call(cmd,shell=True)
                else:
                    print("please type 'yes' behind clean argument.")
                    return

    ################
    ##jitterbug tool
    ################
    if args.jitterbug is not None:

        if find_executable(args.repeatmasker) is not None:
            print('repeatmasker executable can be found')
        else:
            print("Cannot find repeatmasker executable, please check if it has been installed.")
            return

        ##create the working_dir
        output_jitterbug_dir = working_dir + '/output_jitterbug'
        if not os.path.exists(output_jitterbug_dir):
            os.makedirs(output_jitterbug_dir)

        repeatmasker_dir = output_jitterbug_dir + '/repeatmasker_dir'
        if not os.path.exists(repeatmasker_dir):
            os.makedirs(repeatmasker_dir)

        temp_gff3_dir = output_jitterbug_dir + '/temp_gff3_dir'
        if not os.path.exists(temp_gff3_dir):
            os.makedirs(temp_gff3_dir)

        samples_output_dir = output_jitterbug_dir + '/samples_output_dir'
        if not os.path.exists(samples_output_dir):
            os.makedirs(samples_output_dir)

        jitterbug_dir = args.jitterbug

        ##prepare the files
        repeatmasker_tool = args.repeatmasker
        ##create gff from the repeatmasker
        cmd = repeatmasker_tool + \
              ' ' + genome_file + \
              ' -lib ' + library_file + \
              ' -gff' + \
              ' -dir ' + repeatmasker_dir + \
              ' -par ' + str(pro_num) + \
              ' -nolow'
        subprocess.call(cmd,shell=True)

        ##transfer the gff file
        store_final_line_gff3_list = change_gff_to_gff3(repeatmasker_dir + '/opt_rm_te_genome.fa.out.gff')
        with open(temp_gff3_dir + '/opt_te.gff3', 'w+') as opt:
            for eachline in store_final_line_gff3_list:
                opt.write(eachline + '\n')

        bam_fl_dir = args.bam_bwa_dir

        ##find the fastq pair
        for eachnm in sample_nm_dic:

            ##get the new name of fastq
            new_nm = sample_nm_dic[eachnm]

            ##create the directories
            opt_jitterbug_temp_dir = samples_output_dir + '/' + new_nm
            if not os.path.exists(opt_jitterbug_temp_dir):
                os.makedirs(opt_jitterbug_temp_dir)

            ##initiate a list to store the pair list
            mt_pair_list = []
            ##cp the raw fastq file to a new fastq file with new name
            for eachfastq_fl in fastq_list:
                ##extract the name from the eachfastq
                mt = re.match('(.+)/(.+)', eachfastq_fl)
                path = mt.group(1)
                fl_nm = mt.group(2)
                mt = re.match('(.+)_(\d)(\..+)', fl_nm)
                fq_nm = mt.group(1)
                pair_num = mt.group(2)
                type_nm = mt.group(3)

                ##updation 11.23 change new_nm to fq_nm in the append
                if eachnm == fq_nm:
                    mt_pair_list.append(path + '/' + fq_nm + '_' + pair_num + type_nm)

            ##prepare the bam files
            ipt_bam_fl_path = bam_fl_dir + '/' + new_nm + '.bam'

            ##run each sample
            cmd = 'chmod -R 777 ' + opt_jitterbug_temp_dir
            subprocess.call(cmd,shell=True)

            cmd = 'python ' + jitterbug_dir + '/jitterbug.py' + \
                  ' --numCPUs ' + str(pro_num) + \
                  ' --output_prefix ' + opt_jitterbug_temp_dir + '/' + new_nm + \
                  ' ' + ipt_bam_fl_path + \
                  ' ' + temp_gff3_dir + '/opt_te.gff3'
            subprocess.call(cmd,shell=True)

            ##filter the results
            cmd = 'python ' + jitterbug_dir + '/tools/jitterbug_filter_results_func.py' + \
                  ' -g ' + opt_jitterbug_temp_dir + '/' + new_nm + '.TE_insertions_paired_clusters.gff3' + \
                  ' -c ' + opt_jitterbug_temp_dir + '/' + new_nm + '.filter_config.txt' + \
                  ' -o ' + opt_jitterbug_temp_dir + '/' + new_nm + '.TE_insertions_paired_clusters.filtered.gff3'
            subprocess.call(cmd,shell=True)

            ##write the results to the final output dir
            ##modify the output and directly write to the output dir
            jitterbug_change_format_final_list = change_format_jitterbug(opt_jitterbug_temp_dir + '/' + new_nm + '.TE_insertions_paired_clusters.filtered.gff3')
            with open(sample_bed_dir + '/' + new_nm + '_jitterbug.bed', 'w+') as opt:
                for eachline in jitterbug_change_format_final_list:
                    opt.write(eachline + '\n')

    #############
    ##TEFLoN tool
    #############
    if args.TEFLoN is not None:

        if find_executable(args.repeatmasker) is not None:
            print('repeatmasker executable can be found')
        else:
            print("Cannot find repeatmasker executable, please check if it has been installed.")
            return

        ##create the working_dir
        output_TEFLoN_dir = working_dir + '/output_TEFLoN'
        if not os.path.exists(output_TEFLoN_dir):
            os.makedirs(output_TEFLoN_dir)

        samples_output_dir = output_TEFLoN_dir + '/samples_output_dir'
        if not os.path.exists(samples_output_dir):
            os.makedirs(samples_output_dir)

        TEFLoN_dir = args.TEFLoN
        repeatmasker_tool = args.repeatmasker

        ##modify the lib
        store_seq_dic = change_TEnm_lib(library_file)

        with open(output_TEFLoN_dir + '/opt_modified_TElib.fa', 'w+') as opt:
            for eachnm in store_seq_dic:
                opt.write('>' + eachnm + '\n' + store_seq_dic[eachnm] + '\n')

        bwa_tool = args.bwa
        samtools_tool = args.samtools

        ##find the fastq pair
        for eachnm in sample_nm_dic:
            ##get the new name of fastq
            new_nm = sample_nm_dic[eachnm]

            ##initiate a list to store the pair list
            mt_pair_list = []
            ##cp the raw fastq file to a new fastq file with new name
            for eachfastq_fl in fastq_list:
                ##extract the name from the eachfastq
                mt = re.match('(.+)/(.+)', eachfastq_fl)
                path = mt.group(1)
                fl_nm = mt.group(2)
                mt = re.match('(.+)_(\d)(\..+)', fl_nm)
                fq_nm = mt.group(1)
                pair_num = mt.group(2)
                type_nm = mt.group(3)
                ##updation 11.23 change new_nm to fq_nm in the append
                if eachnm == fq_nm:
                    mt_pair_list.append(path + '/' + fq_nm + '_' + pair_num + type_nm)
            ##create the directories
            opt_TEFLoN_temp_dir = samples_output_dir + '/' + new_nm
            if not os.path.exists(opt_TEFLoN_temp_dir):
                os.makedirs(opt_TEFLoN_temp_dir)

            ##step 1
            cmd = 'python ' + TEFLoN_dir + '/teflon_prep_custom.py' + \
                  ' -wd ' + opt_TEFLoN_temp_dir + \
                  ' -e ' + repeatmasker_tool + \
                  ' -g ' + genome_file + \
                  ' -l ' + output_TEFLoN_dir + '/opt_modified_TElib.fa' + \
                  ' -p ' + new_nm + \
                  ' -t ' + str(pro_num)
            subprocess.call(cmd,shell=True)

            ##step 2
            cmd = bwa_tool + ' index ' + opt_TEFLoN_temp_dir + '/' + new_nm + '.prep_MP/' + new_nm + '.mappingRef.fa'
            subprocess.call(cmd, shell=True)

            cmd = bwa_tool + ' mem -t ' + str(pro_num) + \
                  ' ' + opt_TEFLoN_temp_dir + '/' + new_nm + '.prep_MP/' + new_nm + '.mappingRef.fa' + \
                  ' ' + mt_pair_list[0] + \
                  ' ' + mt_pair_list[1] + \
                  ' > ' + opt_TEFLoN_temp_dir + '/alignment.sam'
            subprocess.call(cmd, shell=True)

            ##step 3
            cmd = samtools_tool + ' view -bS -h ' + \
                  opt_TEFLoN_temp_dir + '/alignment.sam -@ ' + str(pro_num) + ' > ' + opt_TEFLoN_temp_dir + '/alignment.bam'
            subprocess.call(cmd, shell=True)

            ##step 4
            cmd = samtools_tool + ' sort ' + opt_TEFLoN_temp_dir + '/alignment.bam' + \
                  ' -o ' + opt_TEFLoN_temp_dir + '/alignment_sorted.bam' + \
                  ' -@ ' + str(pro_num) + \
                  ' > ' + opt_TEFLoN_temp_dir + '/alignment_sorted.bam'
            subprocess.call(cmd, shell=True)

            ##step 5
            ##samtools index alignment_sorted.bam
            cmd = samtools_tool + ' index ' + opt_TEFLoN_temp_dir + '/alignment_sorted.bam'
            subprocess.call(cmd, shell=True)

            ##step 6
            with open(opt_TEFLoN_temp_dir + '/samples.txt', 'w+') as opt:
                opt.write(opt_TEFLoN_temp_dir + '/alignment_sorted.bam' + '\t' + new_nm)

            cmd = 'python ' + TEFLoN_dir + '/teflon.v0.4.py' + \
                  ' -wd ' + opt_TEFLoN_temp_dir + \
                  ' -d ' + opt_TEFLoN_temp_dir + '/' + new_nm + '.prep_TF' + \
                  ' -s ' + opt_TEFLoN_temp_dir + '/samples.txt' + \
                  ' -i ' + new_nm + \
                  ' -eb ' + bwa_tool + \
                  ' -es ' + samtools_tool + \
                  ' -l1 family -l2 family' + \
                  ' -q 10' + \
                  ' -t ' + str(pro_num)
            subprocess.call(cmd, shell=True)

            ##step 7
            cmd = 'python ' + TEFLoN_dir + '/teflon_collapse.py' + \
                  ' -wd ' + opt_TEFLoN_temp_dir + \
                  ' -d ' + opt_TEFLoN_temp_dir + '/' + new_nm + '.prep_TF' + \
                  ' -s ' + opt_TEFLoN_temp_dir + '/samples.txt' + \
                  ' -es ' + samtools_tool + \
                  ' -n1 1 -n2 1 -q 10 -t ' + str(pro_num)
            subprocess.call(cmd, shell=True)

            ##step 8
            cmd = 'python ' + TEFLoN_dir + '/teflon_count.py' + \
                  ' -wd ' + opt_TEFLoN_temp_dir + \
                  ' -d ' + opt_TEFLoN_temp_dir + '/' + new_nm + '.prep_TF' + \
                  ' -s ' + opt_TEFLoN_temp_dir + '/samples.txt' + \
                  ' -i ' + new_nm + \
                  ' -eb ' + bwa_tool + \
                  ' -es ' + samtools_tool + \
                  ' -l2 family' + \
                  ' -q 10 ' + \
                  ' -t ' + str(pro_num)
            subprocess.call(cmd, shell=True)

            ##step 9
            cmd = 'python ' + TEFLoN_dir + '/teflon_genotype.py' + \
                  ' -wd ' + opt_TEFLoN_temp_dir + \
                  ' -d ' + opt_TEFLoN_temp_dir + '/' + new_nm + '.prep_TF' + \
                  ' -s ' + opt_TEFLoN_temp_dir + '/samples.txt' + \
                  ' -dt pooled'
            subprocess.call(cmd, shell=True)

            TEFLoN_change_format_final_list = change_format_TEFLoN(opt_TEFLoN_temp_dir + '/genotypes/' + new_nm + '.genotypes.txt')
            with open(sample_bed_dir + '/' + new_nm + '_TEFLoN.bed', 'w+') as opt:
                for eachline in TEFLoN_change_format_final_list:
                    opt.write(eachline + '\n')


if __name__ == "__main__":
    main()
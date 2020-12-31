#!/usr/bin/env python

##updation 122420 modify the clean function for the annotation step
##updation 122120 correct the dividing files
##udpation 122020 change annotation blast database construction
##updation 121920 add another threshold to increase the range of combine case
##updation 121920 add a threshold to define the range of single case we need to analyze
##updation 110220 modify combine ref close loci to remove the c or s covered with the o type
##updation 103120 add the panTE combine and a threshold
##updation 102820 add the gap version
##updation 102520 use a new genotyping to have an analysis this step now only contains script to treat the combined location
##not for the reference and single

##updation 101720 increase default threads add the total cover threshold
##updation 060820 equally split the file
##updaiton 12.31 rm the bai construction in this script
##updation 12.03 add the threads number choosing for the samtools
##updation 10.28 genos modification add output path
##this script pipeline is to genotype the TE insertions

##BUILT-IN MODULES
import argparse
import sys
from multiprocessing import Pool
import os
import subprocess
import glob
import re


from TEmarker_utils import TM_genos_construct_parallele as genos_construct
#from TEmarker_utils import TM_genos_modify as genos_modify
from TEmarker_utils import TM_genos_combine_ref_close_loci as genos_combine_loci


##SCRIPTS
def get_parsed_args():
    parser = argparse.ArgumentParser(description="Genotype TE markers")

    ##############################
    ##parse the required arguments
    parser.add_argument("-d", dest='working_dir', default="./", help="Working directory"
                                                                     "Default: ./ ")

    parser.add_argument("-o", dest='output_dir', default="./", help="Output directory"
                                                                    "Default: ./ ")

    parser.add_argument('-panTEs_f', dest='candidate_TE_fl', help="Users provide the panTEs file generated from TEmarker_panTEs.py")

    parser.add_argument('-bam_bwa_d', dest='bam_bwa_dir', help="Users provide the bam dir generated from the TEmarker_bam.py.")

    parser.add_argument('-bam_hisat2_d', dest='bam_hisat2_dir', help="Users provide the bam dir generated from the TEmarker_bam.py.")

    parser.add_argument('-lib_f', dest='TE_lib_fl', help="Users provide the TE lib file to decide the TE family of unknown insertion")

    ##option
    parser.add_argument('-range_bp', dest='support_thr', help="Users provide a range (bp) to decide if reads support TE insertion."
                                                         "Default: 15")

    parser.add_argument('-pros_n', dest='process_num', help="Users provide process number that helps to divide the candidate TE insertion file"
                                                            "into multiple files, which can increase the speed of processing."
                                                            "Default: 10")

    parser.add_argument('-acc_thr', dest='accuracy_thr', help="Users set a threshold to decide if reads support the deletion of a TE."
                                                             "For example, accuracy is 0.9. If we set acc_th as 0.8, suggesting this read support"
                                                             "the deletion event of a TE."
                                                             "Default: 0.8")

    parser.add_argument('-heter_thr', dest='genotype_thr', help="Users set a threshold to decide the homozygous or heterozygous of the candidate TE markers."
                                                               "For example, we will genotype a TE marker diplaying insertions in samples comparing with "
                                                               "the reference TEs"
                                                               "Proportion of insertion reads is calculated first."
                                                               "Compare genos_th to this proportion."
                                                               "If we set genos_th as 0.7"
                                                               "0/0: 0 <= proportion < 0.3"
                                                               "1/1: 0.7 < proportion <= 1"
                                                               "0/1: 0.3 < proportion <= 0.7."
                                                               "Default: 0.7")

    parser.add_argument('-t_cover_thr', dest='total_cover_thr', help="Users set a threshold to remove candidate insertions covering less than the total_cover_thr."
                                                                     "Default: 3")

    parser.add_argument('-clrd_thr', dest='clipped_read_thr', help="Users set a threshold to decide average of clipped reads for each sample."
                                                                   "Default: 2")

    parser.add_argument('-TSD_thr', dest='TSD_thr', help="Users set a threshold to decide the TSD length. If the TSD length is larger than the threshold, "
                                                         "the two locations will be regarded as two independent insertions."
                                                          "Default: 15")

    parser.add_argument('-comb_thr', dest='comb_close_thr', help="Users define a combined close threshold to determine the close genotyping location will be combined."
                                                                 "Default: 3")

    parser.add_argument('-antgap_thr',dest='annot_gap_thr', help="The nearby TE insertions will be combined within a gap during annotation."
                                                                 "Default: 15")

    parser.add_argument('-miss_thr', dest='thr_miss', help="Users provide a threshold to generate a sample number threshold."
                                                           "For example, if thr_miss is equal to 0.7, and the sample number threshold is equal to 0.7*total_sample_number. "
                                                           "If sample number with genotype information is over than the 0.7*total_sample_number, this location will be added to compare to choose the best location in the annotation modifcation step."
                                                           "Default: 0.7")

    ##udpation 121920
    parser.add_argument('-s_rg_thr', dest='s_rg_thr', help="Users provide a range to search if there is other candidate insertion around the single insertion or combine insertion case labeled with s in the panTEs.txt."
                                                           "Default: 50")

    ##remove the c_rg_thr
    #parser.add_argument('-c_rg_thr', dest='c_rg_thr', help="Users provide a range to search if there is other candidate insertion around the combine insertion case labeled with c in the panTEs.txt."
    #                                                       "Default: 50")
    ##updating 122220
    parser.add_argument('-antclgap_thr',dest='annot_close_gap_thr', help="The nearby TE with very close distance will be combined within a gap during annotation even the family is not the same."
                                                                         "Default:5")


    parser.add_argument('-modify', dest='m_yes', help="Users activate the modification on the genotyping file based on the CNV detection dir")

    parser.add_argument('-cnv_d', dest='cnv_dir', help="After activation of modification, users must provide the cnv_dir generated from TEScape_marker_genos_cnv")

    ##add choice for the genotyping and annotation
    parser.add_argument('-disable_g',dest='disable_g', help='Users do not conduct the genotyping.')

    parser.add_argument('-disable_a',dest='disable_a', help='Users do not conduct the annotation.')


    ##parse of parameters
    args = parser.parse_args()
    return args


def multi_run_wrapper_genos_construct(args):
    return genos_construct.genos_construct(*args)

def multi_run_wrapper_genos_annotation(args):
    return genos_construct.annotation(*args)


def main(argv=None):
    if argv is None:
        argv = sys.argv
    args = get_parsed_args()

    #######################################
    ##check the required software and files
    ##for the input files

    candidate_TE_fl = args.candidate_TE_fl
    if args.candidate_TE_fl is None:
        print ('Cannot find input candidate TE file, please provide that')
        return ##import to close the script if not find the te lib
    else:
        try:
            file = open(args.candidate_TE_fl, 'r')  ##check if the file is not the right file
        except IOError:
            print('There was an error opening the file!')
            return

    bam_bwa_dir = args.bam_bwa_dir
    if args.bam_bwa_dir is None:
        print('Cannot find input te bam bwa dir, please provide that')
        return  ##import to close the script if not find the te lib

    bam_hisat2_dir = args.bam_hisat2_dir
    if args.bam_hisat2_dir is None:
        print('Cannot find input te bam hisat dir, please provide that')
        return  ##import to close the script if not find the te lib

    TE_lib_fl = args.TE_lib_fl
    if args.TE_lib_fl is None:
        print ('Cannot find input TE library file, please provide that')
        return ##import to close the script if not find the te lib
    else:
        try:
            file = open(args.TE_lib_fl, 'r')  ##check if the file is not the right file
        except IOError:
            print('There was an error opening the file!')
            return

    cnv_dir = ''
    ##if args.m_yes has been provided
    if args.m_yes is not None:
        if args.m_yes != '1':
            print('Please write 1 to initate the modification')
            return
        else:
            if args.cnv_dir is None:
                print('Please provide the cnv_dir')
                return
            else:
                cnv_dir = args.cnv_dir

    if args.support_thr is not None:
        support_thr = args.support_thr
    else:
        support_thr = 15

    ##updation 101720
    if args.total_cover_thr is not None:
        total_cover_thr = args.total_cover_thr
    else:
        total_cover_thr = 3

    if args.process_num is not None:
        process_num = args.process_num
    else:
        process_num = 10

    if args.accuracy_thr is not None:
        accuracy_thr = args.accuracy_thr
    else:
        accuracy_thr = 0.8

    if args.genotype_thr is not None:
        genotype_thr = args.genotype_thr
    else:
        genotype_thr = 0.7

    if args.thr_miss is not None:
        thr_miss = args.thr_miss
    else:
        thr_miss = 0.7

    if args.clipped_read_thr is not None:
        clrd_thr = args.clipped_read_thr
    else:
        clrd_thr = '2'

    if args.TSD_thr is not None:
        TSD_thr = args.TSD_thr
    else:
        TSD_thr = '15'

    if args.comb_close_thr is not None:
        comb_close_thr = args.comb_close_thr
    else:
        comb_close_thr = '3'

    if args.annot_gap_thr is not None:
        annot_gap_thr = args.annot_gap_thr
    else:
        annot_gap_thr = '15'

    ##updation 121920
    if args.s_rg_thr is not None:
        s_rg_thr = args.s_rg_thr
    else:
        s_rg_thr = '50'

    ##udpation 122120
    #if args.c_rg_thr is not None:
    #    c_rg_thr = args.c_rg_thr
    #else:
    #    c_rg_thr = '50'

    ##updating 122220
    if args.annot_close_gap_thr is not None:
        annot_close_gap_thr = args.annot_close_gap_thr
    else:
        annot_close_gap_thr = '15'




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

    ##run the scripts
    ##step 1: create the genotype file
    rm_pcr_mcc_bam_dir = bam_bwa_dir
    rm_pcr_hisat2_bam_dir = bam_hisat2_dir

    sample_number = genos_construct.generate_sample_dic(rm_pcr_mcc_bam_dir)
    #sample_number = len(list(store_sample_dic.key()))

    D01_construct_genos_dir = working_dir + '/D01_construct_genos_dir'
    if not os.path.exists(D01_construct_genos_dir):
        os.makedirs(D01_construct_genos_dir)

    genos_construct.create_working_genos_construct_dir(D01_construct_genos_dir)
    #genos_construct.generate_whole_te_library(TE_lib_fl, D01_construct_genos_dir + '/temp_store_whole_library_dir')

    ##updation 102520
    ##generate D02 for the annotation
    D02_genos_annotation_dir = working_dir + '/D02_genos_annotation_dir'
    if not os.path.exists(D02_genos_annotation_dir):
        os.makedirs(D02_genos_annotation_dir)

    ##input the multi processor num
    pool = Pool(int(process_num))


    ##updation 103120 make some modification for the panTEs file
    store_modifed_panTEs_line_list = genos_construct.combine_panTEs_location (candidate_TE_fl,D01_construct_genos_dir)
    with open (D01_construct_genos_dir + '/opt_temp_panTEs.txt','w+') as opt:
        for eachline in store_modifed_panTEs_line_list:
            opt.write(eachline + '\n')

    ##decide how many files will be splited
    ##cp the loc into the temp split dir
    cmd = 'cp ' + D01_construct_genos_dir + '/opt_temp_panTEs.txt' + ' ' + D01_construct_genos_dir + '/temp_split_dir'
    print(cmd)
    subprocess.call(cmd, shell=True)

    ##udpation 122120 change the candidate TE fl to the opt_temp_panTEs.txt
    ##updation 060820
    count_can_TEs = 0
    #with open (candidate_TE_fl,'r') as ipt:
    #    for eachline in ipt:
    #        count_can_TEs += 1
    with open (D01_construct_genos_dir + '/temp_split_dir/opt_temp_panTEs.txt','r') as ipt:
        for eachline in ipt:
            count_can_TEs += 1

    div_num = int(count_can_TEs/int(process_num)) ##even 5.5 it is 5
    if div_num != count_can_TEs / int(process_num):
        #filter_te_loc_fl_list = glob.glob(D01_construct_genos_dir + '/temp_split_dir' + '/*')
        cmd = 'split -l ' + str(div_num) + ' ' + D01_construct_genos_dir + '/temp_split_dir/opt_temp_panTEs.txt' + ' ' + \
              D01_construct_genos_dir + '/temp_split_dir' + '/temp_split_'
        print(cmd)
        subprocess.call(cmd,shell=True)
        ##combine the last two files
        ##updation 122120
        split_te_loc_fl_list = sorted(glob.glob(D01_construct_genos_dir + '/temp_split_dir' + '/temp_split_*'))

        ##updation 122120
        total_num_split_files = len(split_te_loc_fl_list)
        ##caculate how many files needs to be combine
        add_num_split_files = total_num_split_files - int(process_num)
        ##we need to combine the added files and the final file
        ##if we have two more files, the last_but_one_fl_nm will be split_te_loc_fl_list[-3]
        mt = re.match('.+/(.+)', split_te_loc_fl_list[-(add_num_split_files + 1)])
        last_but_one_fl_nm = mt.group(1)
        ##we need to combine the -1 -2 -3 together
        combine_fl_str = ''
        for i in range(-(add_num_split_files + 1), 0):
            print(i)
            combine_fl_str = combine_fl_str + ' ' + split_te_loc_fl_list[i]
        cmd = 'cat ' + combine_fl_str + ' > ' + D01_construct_genos_dir + '/temp_split_dir/temp_split_lasttwo'
        print(cmd)
        subprocess.call(cmd, shell=True)

        cmd = 'rm ' + combine_fl_str
        print(cmd)
        subprocess.call(cmd, shell=True)

        cmd = 'mv ' + D01_construct_genos_dir + '/temp_split_dir/temp_split_lasttwo ' + D01_construct_genos_dir + '/temp_split_dir/' + last_but_one_fl_nm
        print(cmd)
        subprocess.call(cmd,shell=True)

        ##extract the name of the last but one
        #mt = re.match('.+/(.+)',split_te_loc_fl_list[-2])
        #last_but_one_fl_nm = mt.group(1)
        ##combine
        #cmd = 'cat ' + split_te_loc_fl_list[-1] + ' ' + split_te_loc_fl_list[-2] + ' > ' + D01_construct_genos_dir + '/temp_split_dir/temp_split_lasttwo'
        #print(cmd)
        #subprocess.call(cmd,shell=True)
        ##rm the last two
        #cmd = 'rm ' + split_te_loc_fl_list[-1] + ' ' + split_te_loc_fl_list[-2]
        #print(cmd)
        #subprocess.call(cmd,shell=True)
        ##change the name
        #cmd = 'mv ' + D01_construct_genos_dir + '/temp_split_dir/temp_split_lasttwo ' + D01_construct_genos_dir + '/temp_split_dir/' + last_but_one_fl_nm
        #print(cmd)
        #subprocess.call(cmd,shell=True)
    else:
        ##do not need to combine the last two
        #filter_te_loc_fl_list = glob.glob(D01_construct_genos_dir + '/temp_split_dir' + '/*')
        cmd = 'split -l ' + str(div_num) + ' ' + D01_construct_genos_dir + '/temp_split_dir/opt_temp_panTEs.txt' + ' ' + \
              D01_construct_genos_dir + '/temp_split_dir' + '/temp_split_'
        print(cmd)
        subprocess.call(cmd, shell=True)

    ##updation 122120
    ##sort the temp_split files
    temp_split_file_unsorted_list = glob.glob(D01_construct_genos_dir + '/temp_split_dir' + '/temp_split*')
    for eachfl in temp_split_file_unsorted_list:
        mt = re.match('.+/(.+)',eachfl)
        fl_nm = mt.group(1)
        cmd = 'sort -k1,1V -k2,2n ' + eachfl + ' > ' + D01_construct_genos_dir + '/temp_sorted_split_dir/' + fl_nm
        subprocess.call(cmd, shell=True)


    ##get the temp split location files
    #temp_split_file_list = glob.glob(D01_construct_genos_dir + '/temp_split_dir' + '/temp_split*')
    ##updation 122120
    temp_split_file_list = glob.glob(D01_construct_genos_dir + '/temp_sorted_split_dir' + '/temp_split*')

    ##updation 8.22
    ##generate the temp blast dir for each run process
    ##like the temp sam file dir
    ##get the temp blast location files
    # temp_blast_dir_list = glob.glob(temp_store_blast_index_dir + '/*')

    ##updation 8.23
    ##generate the changing list dir
    #for x in range(0, int(process_num)):
    #    dir_code = x + 1
    #    temp_changingte_dir = D01_construct_genos_dir + '/temp_store_changing_te_dir' + '/temp_chaning_0' + str(dir_code) + 'dir'
    #    if not os.path.exists(temp_changingte_dir):
    #        os.makedirs(temp_changingte_dir)


    ##updation 102520
    ##temp_loc_line_dir
    for x in range(0, int(process_num)):
        dir_code = x + 1
        temp_loc_line_dir = D01_construct_genos_dir + '/temp_loc_line_dir' + '/temp_loc_line_0' + str(dir_code) + 'dir'
        if not os.path.exists(temp_loc_line_dir):
            os.makedirs(temp_loc_line_dir)

    ##updation 102520
    ##generate a temp store output
    for x in range(0, int(process_num)):
        dir_code = x + 1
        temp_store_opt_dir = D01_construct_genos_dir + '/temp_split_opt_dir' + '/temp_store_opt_0' + str(dir_code) + 'dir'
        if not os.path.exists(temp_store_opt_dir):
            os.makedirs(temp_store_opt_dir)



    temp_store_opt_dir_list = glob.glob(D01_construct_genos_dir + '/temp_split_opt_dir' + '/*')



    if args.disable_g is not None:
        if args.disable_g != '1':
            print ('Please write 1 to disable the genotyping')
            return
        else:
            print('Users stop the genotyping')
    else:
        print('Conduct the genotyping')

        ##generate a list to store the define_hetero function
        run_list = []
        ##set a loop to run the splited file
        ##updation 102520
        # store_final_line_list = [] #do not use the final line list we will directly write out the lines
        #opt_te_genotype_fl = open(output_dir + '/opt_te_genotype.txt', 'w')  # w when u wanna write sth on the file
        ##(input_panTEs_file, input_bam_fl_dir, input_output_dir,opt_te_genotype_fl,sample_num, clrd_thr, TSD_thr, t_cover_thr, Her_thd,pa_num)

        ##updation 121920 add the s_rg_thr
        ##updation 122120 remove the c_rg_thr
        for x in range(0, int(process_num)):
            each_func_argument = (temp_split_file_list[x],
                                  rm_pcr_mcc_bam_dir,
                                  rm_pcr_hisat2_bam_dir,
                                  D01_construct_genos_dir,
                                  temp_store_opt_dir_list[x],
                                  sample_number,
                                  clrd_thr,
                                  TSD_thr,
                                  total_cover_thr,
                                  genotype_thr,
                                  accuracy_thr,
                                  comb_close_thr,
                                  s_rg_thr,
                                  x)
            run_list.append(each_func_argument)
        pool.map(multi_run_wrapper_genos_construct, run_list)
        #opt_te_genotype_fl.close()

        #for x in range(0, int(process_num)):
        #    each_func_argument = (rm_pcr_mcc_bam_dir, temp_split_file_list[x], rm_pcr_hisat2_bam_dir,
        #                          D01_construct_genos_dir,store_sample_dic,accuracy_thr, genotype_thr, 3, x,support_thr,total_cover_thr)
        #    run_list.append(each_func_argument)
        #res_dic_list = pool.map(multi_run_wrapper, run_list)

        #count = 0
        #for res_dic in res_dic_list:
        #    count += 1
        #    with open(D01_construct_genos_dir + '/temp_split_opt_dir' + '/opt_classify_hetero' + '_split0' + str(count) + '.txt','w+') as opt:
        #        for eachline in res_dic:
        #            opt.write(eachline + '\n')
        ##combine all the opt_classify_hetero_split files into one file
        ##set a cat string
        file_string = ''
        for x in range(0, int(process_num)):
            dir_code = x + 1
            classify_hetero_file_list = glob.glob(D01_construct_genos_dir + '/temp_split_opt_dir' + '/temp_store_opt_0' + str(dir_code) + 'dir/*')
            for eachfl in classify_hetero_file_list:
                file_string = file_string + ' ' + eachfl

        ##run the cat
        cmd = 'cat ' + file_string + ' > ' + D01_construct_genos_dir + '/temp_te_genotype.txt'
        subprocess.call(cmd, shell=True)

    ##updaton 102520
    ############
    ##annotation
    ############
    if args.disable_a is not None:
        if args.disable_a != '1':
            print ('Please write 1 to disable the annotation')
            return
        else:
            print('Users stop the annotation')
    else:
        print('Conduct the annotation')

        genos_construct.create_working_genos_annotation_dir(D02_genos_annotation_dir)
        genos_construct.rc_te_lib(TE_lib_fl, D02_genos_annotation_dir)

        ##updation 102520
        cmd = 'cp ' + D01_construct_genos_dir + '/temp_te_genotype.txt' + ' ' + D02_genos_annotation_dir + '/temp_split_dir'
        print(cmd)
        subprocess.call(cmd, shell=True)

        ##if we use the fitration step we need to use the opt_te_genotype_filtration.txt
        count_can_TEs = 0
        with open(D01_construct_genos_dir + '/temp_te_genotype.txt', 'r') as ipt:
            for eachline in ipt:
                count_can_TEs += 1

        div_num = int(count_can_TEs / int(process_num))  ##even 5.5 it is 5
        if div_num != count_can_TEs / int(process_num):
            #filter_te_loc_fl_list = glob.glob(D02_genos_annotation_dir + '/temp_split_dir' + '/*')
            cmd = 'split -l ' + str(div_num) + ' ' + D02_genos_annotation_dir + '/temp_split_dir' + '/temp_te_genotype.txt' + ' ' + \
                  D02_genos_annotation_dir + '/temp_split_dir' + '/temp_split_'
            print(cmd)
            subprocess.call(cmd, shell=True)
            ##combine the last two files
            split_te_loc_fl_list = sorted(glob.glob(D02_genos_annotation_dir + '/temp_split_dir' + '/temp_split_*'))


            ##updation 122120
            total_num_split_files = len(split_te_loc_fl_list)
            ##caculate how many files needs to be combine
            add_num_split_files = total_num_split_files - int(process_num)
            ##we need to combine the added files and the final file
            ##if we have two more files, the last_but_one_fl_nm will be split_te_loc_fl_list[-3]
            mt = re.match('.+/(.+)', split_te_loc_fl_list[-(add_num_split_files + 1)])
            last_but_one_fl_nm = mt.group(1)
            ##we need to combine the -1 -2 -3 together
            combine_fl_str = ''
            for i in range(-(add_num_split_files + 1), 0):
                print(i)
                combine_fl_str = combine_fl_str + ' ' + split_te_loc_fl_list[i]
            cmd = 'cat ' + combine_fl_str + ' > ' + D02_genos_annotation_dir + '/temp_split_dir/temp_split_lasttwo'
            print(cmd)
            subprocess.call(cmd, shell=True)

            cmd = 'rm ' + combine_fl_str
            print(cmd)
            subprocess.call(cmd, shell=True)

            cmd = 'mv ' + D02_genos_annotation_dir + '/temp_split_dir/temp_split_lasttwo ' + D02_genos_annotation_dir + '/temp_split_dir/' + last_but_one_fl_nm
            print(cmd)
            subprocess.call(cmd, shell=True)



            ##extract the name of the last but one
            #mt = re.match('.+/(.+)', split_te_loc_fl_list[-2])
            #last_but_one_fl_nm = mt.group(1)
            ##combine
            #cmd = 'cat ' + split_te_loc_fl_list[-1] + ' ' + split_te_loc_fl_list[
            #    -2] + ' > ' + D02_genos_annotation_dir + '/temp_split_dir/temp_split_lasttwo'
            #print(cmd)
            #subprocess.call(cmd, shell=True)
            ##rm the last two
            #cmd = 'rm ' + split_te_loc_fl_list[-1] + ' ' + split_te_loc_fl_list[-2]
            #print(cmd)
            #subprocess.call(cmd, shell=True)
            ##change the name
            #cmd = 'mv ' + D02_genos_annotation_dir + '/temp_split_dir/temp_split_lasttwo ' + D02_genos_annotation_dir + '/temp_split_dir/' + last_but_one_fl_nm
            #print(cmd)
            #subprocess.call(cmd, shell=True)
        else:
            ##do not need to combine the last two
            #filter_te_loc_fl_list = glob.glob(D02_genos_annotation_dir + '/temp_split_dir' + '/*')
            cmd = 'split -l ' + str(div_num) + ' ' + D02_genos_annotation_dir + '/temp_split_dir' + '/temp_te_genotype.txt' + ' ' + \
                  D02_genos_annotation_dir + '/temp_split_dir' + '/temp_split_'
            print(cmd)
            subprocess.call(cmd, shell=True)

        ##get the temp split location files
        temp_split_file_list = glob.glob(D02_genos_annotation_dir + '/temp_split_dir' + '/temp_split*')

        for x in range(0, int(process_num)):
            dir_code = x + 1
            temp_blast_dir = D02_genos_annotation_dir + '/temp_store_blast_index_dir' + '/temp_blast_0' + str(dir_code) + 'dir'
            if not os.path.exists(temp_blast_dir):
                os.makedirs(temp_blast_dir)

        ##mkdir for the temp sam dir under the temp_multi_sam_save_dir
        for x in range(0, int(process_num)):
            dir_code = x + 1
            temp_sam_dir = D02_genos_annotation_dir + '/temp_multi_sam_save_dir' + '/temp_sam_0' + str(dir_code) + 'dir'
            if not os.path.exists(temp_sam_dir):
                os.makedirs(temp_sam_dir)

        ##updation 102520
        ##generate a temp store output
        for x in range(0, int(process_num)):
            dir_code = x + 1
            temp_store_opt_dir = D02_genos_annotation_dir + '/temp_split_opt_dir' + '/temp_store_opt_0' + str(dir_code) + 'dir'
            if not os.path.exists(temp_store_opt_dir):
                os.makedirs(temp_store_opt_dir)

        temp_store_opt_dir_list = glob.glob(D02_genos_annotation_dir + '/temp_split_opt_dir' + '/*')

        ##updation 122020
        ##mkdir for the temp preTE_db_dir
        for x in range(0, int(process_num)):
            dir_code = x + 1
            temp_db_dir = D02_genos_annotation_dir + '/temp_preTE_db_dir/temp_preTEdb_0' + str(dir_code) + 'dir'
            if not os.path.exists(temp_db_dir):
                os.makedirs(temp_db_dir)


        # annotation(opt_te_genotype_fl, opt_for_res_te_fl, input_bam_fl_dir, input_output_dir)
        ##generate a list to store the define_hetero function
        run_list = []
        ##set a loop to run the splited file
        ##updation 102520
        # store_final_line_list = [] #do not use the final line list we will directly write out the lines
        #opt_te_genotype_annot_fl = open(output_dir + '/opt_te_genotype_annot.txt', 'w')  # w when u wanna write sth on the file
        ##annotation (opt_te_genotype_annot_fl,input_output_dir + '/opt_te_genotype.txt', input_output_dir + '/opt_for_res_te.lib',input_bam_fl_dir,input_output_dir,pa_num)
        for x in range(0, int(process_num)):
            each_func_argument = (temp_store_opt_dir_list[x],
                                  temp_split_file_list[x],
                                  D02_genos_annotation_dir + '/temp_lib_db_dir',
                                  rm_pcr_mcc_bam_dir,
                                  D02_genos_annotation_dir,
                                  TE_lib_fl,
                                  x)
            run_list.append(each_func_argument)
        pool.map(multi_run_wrapper_genos_annotation, run_list)

        file_string = ''
        for x in range(0, int(process_num)):
            dir_code = x + 1
            classify_hetero_file_list = glob.glob(D02_genos_annotation_dir + '/temp_split_opt_dir' + '/temp_store_opt_0' + str(dir_code) + 'dir/*')
            for eachfl in classify_hetero_file_list:
                file_string = file_string + ' ' + eachfl

        ##run the cat
        cmd = 'cat ' + file_string + ' > ' + D02_genos_annotation_dir + '/temp_te_genotype_annot.txt'
        subprocess.call(cmd, shell=True)

        ##udpation 122020 modify output of the temp_te_genotype_annot.txt
        ##since it contains additional string of TEname
        store_rm_addcol_te_genotype_annot_line_list = []
        with open(D02_genos_annotation_dir + '/temp_te_genotype_annot.txt', 'r') as ipt:
            for eachline in ipt:
                eachline = eachline.strip('\n')
                col = eachline.strip().split()

                ##updation 122420 col[5]
                if col[5].startswith('Class'):
                    ##it means it has double name
                    other_str = ''
                    for i in range(6,len(col)):
                        other_str = other_str + '\t' + col[i]
                    final_line = col[0] + '\t' + col[1] + '\t' + col[2] + '\t' + col[3] + '\t' + col[4] + other_str
                    store_rm_addcol_te_genotype_annot_line_list.append(final_line)

                #if len(col) == 7:
                    ##inidcate there is an additional col
                #    final_line = col[0] + '\t' + col[1] + '\t' + col[2] + '\t' + col[3] + '\t' + col[4] + '\t' + col[6]
                #    store_rm_addcol_te_genotype_annot_line_list.append(final_line)
                else:
                    store_rm_addcol_te_genotype_annot_line_list.append(eachline)

        with open(D02_genos_annotation_dir + '/temp_te_genotype_annot_clean.txt', 'w+') as opt:
            for eachline in store_rm_addcol_te_genotype_annot_line_list:
                opt.write(eachline + '\n')


        ##updation 103120
        ##add the modification of the genotype
        ##updation 110520 generate sample number above the missing threshold
        sample_num_above_missing = int(int(sample_number)*float(thr_miss))

        store_final_annot_line_list = genos_construct.modify_annotation_file (D02_genos_annotation_dir + '/temp_te_genotype_annot_clean.txt',D02_genos_annotation_dir,annot_gap_thr,sample_num_above_missing,annot_close_gap_thr)
        with open (D02_genos_annotation_dir + '/temp_te_genotype_annot_modified.txt','w+') as opt:
            for eachline in store_final_annot_line_list:
                opt.write(eachline + '\n')

    ##updation 072620
    ##sorted the modified genotype to make sure the order is right
    cmd = 'sort -k1,1V -k2,2n ' + D02_genos_annotation_dir + '/temp_te_genotype_annot_modified.txt'  + ' > ' + D02_genos_annotation_dir + '/temp_te_genotype_annot_modified_sorted.txt'
    subprocess.call(cmd, shell=True)

    ##the combine close loci already take all as many as samples to be consider, eg. if a has 0/1 b do not have 0/1 we will select the a 0/1, so this location has a genotype to be 0/1
    store_final_line_list,store_remove_c_s_line_list = genos_combine_loci.combine_close_loci(D02_genos_annotation_dir + '/temp_te_genotype_annot_modified_sorted.txt', D02_genos_annotation_dir)
    with open(output_dir + '/opt_te_genotype.txt', 'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')
    with open(D02_genos_annotation_dir + '/temp_remove_te_genotype_cs_cover_o.txt', 'w+') as opt:
        for eachline in store_remove_c_s_line_list:
            opt.write(eachline + '\n')

    ########
    ##step 2: modify the genotyping file
    ##if users notify this situation
    #if args.m_yes is not None:
    #    final_line_list = genos_modify.modification(output_dir + '/opt_te_genotype.txt', cnv_dir, genotype_thr)
    #    with open(output_dir + '/opt_modified_te_genotype.txt', 'w+') as opt:
    #        for eachline in final_line_list:
    #            opt.write(eachline + '\n')


if __name__ == "__main__":
    main()
#!/usr/bin/env python

##updation 122220 for modifying the annotation, if we identify same location, we need to filter one recover 558 finish
##updation 122220 recover 558 774 finish
##updation 122220 modify the genotyping calling since there are several cases that are not right
##updation 122120 we need to recover 536 576 and 752  all recover
##updation 122120 change the searching range for the c and s to allow them not to be overlapped and only allows one rg parameter
##updation 122020 change annotation firstly based on the mcclintock results and then the full annotation by blasting
##udpation 121920
##updation 110520 we need also compare the missing sample information for the modify_annotation_file since some samples may be filtered out because of the missing we will choose the second good one
##updation 103120 add comb_close_thr (default is 3)
##updation 102820 add a script to combine the close temp clipped location or we will choose the first one to keep
##updation 102820 add the gap version
##updation 102620 change some temp file
##this version is the parallele version wrapped in the TEmarker_genotyping.py
##updation 102520 we will use a new way to conduct the genotyping
##first we need to decide a insertion location based on the bed file generated from the bam by combining all the individual samples
##second annotation, we extract all the clipped reads and map to the one location and extract the match rate to be 100% to do the further comparision

##import modules
import re
import glob
import subprocess
import os
import os.path
from os import path
#from scipy.stats import binom_test
from Bio import SeqIO

#input_panTEs_file = sys.argv[1]
#input_bam_fl_dir = sys.argv[2]
#input_te_lib_fl = sys.argv[3]
#input_output_dir = sys.argv[4]


#MAF_filtration = '1'
#TSD_thr = '15'
#clrd_thr = '2' ##the average read for the inserted samples ##so the total clipped reads threshold = 2*(176*2*0.05)
#sample_num = '176'
#t_cover_thr = '3'
#Her_thd = '0.7'
#s_rg_thr = 10

##step 0: generate sample dic to calculate the sample number
def generate_sample_dic (input_bam_fl_dir):

    ##updation 8.11
    ##generate a sample list that will be used in the next function
    store_sample_dic = {}

    ##this step is to generate the bam index for all the bam files
    mcc_bam_file_list = glob.glob(input_bam_fl_dir + '/*.bam')
    for eachbamfl in mcc_bam_file_list:
        ##change the name of bam file if there is a _2 in the dir
        ##get the fl name  ##attention for this name information
        if '_2' in eachbamfl:
            mt = re.match('.+/(.+)_2\.bam', eachbamfl)
            fl_nm = mt.group(1)
        else:
            mt = re.match('.+/(.+)\.bam', eachbamfl)
            fl_nm = mt.group(1)

        cmd = 'mv ' + eachbamfl + ' ' + input_bam_fl_dir + '/' + fl_nm + '.bam'
        #print(cmd)
        subprocess.call(cmd, shell=True)
        store_sample_dic[fl_nm] = 1

    sample_number = len(list(store_sample_dic.keys()))

    return (sample_number)

##updation 8.24 create temp dir under the working_dir
def create_working_genos_construct_dir (input_working_dir):

    ##updation8.24 create dir under this function
    ##makeblastdb for eachid

    temp_split_dir = input_working_dir + '/temp_split_dir'
    if not os.path.exists(temp_split_dir):
        os.makedirs(temp_split_dir)

    temp_split_opt_dir = input_working_dir + '/temp_split_opt_dir'
    if not os.path.exists(temp_split_opt_dir):
        os.makedirs(temp_split_opt_dir)

    ##updation 102520
    temp_loc_line_dir = input_working_dir + '/temp_loc_line_dir'
    if not os.path.exists(temp_loc_line_dir):
        os.makedirs(temp_loc_line_dir)

    ##updation 122120
    ##sort the panTEs
    temp_sorted_split_dir = input_working_dir + '/temp_sorted_split_dir'
    if not os.path.exists(temp_sorted_split_dir):
        os.makedirs(temp_sorted_split_dir)


def create_working_genos_annotation_dir (input_working_dir):

    ##updation8.24 create dir under this function
    ##makeblastdb for eachid
    temp_multi_sam_save_dir = input_working_dir + '/temp_multi_sam_save_dir'
    if not os.path.exists(temp_multi_sam_save_dir):
        os.makedirs(temp_multi_sam_save_dir)

    temp_split_dir = input_working_dir + '/temp_split_dir'
    if not os.path.exists(temp_split_dir):
        os.makedirs(temp_split_dir)

    temp_store_blast_index_dir = input_working_dir + '/temp_store_blast_index_dir'
    if not os.path.exists(temp_store_blast_index_dir):
        os.makedirs(temp_store_blast_index_dir)

    temp_split_opt_dir = input_working_dir + '/temp_split_opt_dir'
    if not os.path.exists(temp_split_opt_dir):
        os.makedirs(temp_split_opt_dir)

    ##updation 122020
    ##create dir for storing the temp te library database
    temp_preTE_db_dir = input_working_dir + '/temp_preTE_db_dir'
    if not os.path.exists(temp_preTE_db_dir):
        os.makedirs(temp_preTE_db_dir)




##updation 103120 combine the panTEs location
def combine_panTEs_location (input_panTEs_file,input_working_dir):

    ##the output file will be stored in the working_dir
    ##step 1 sort the file
    ##sort the bed fl
    cmd = 'sort -k1,1V -k2,2n ' + input_panTEs_file  + ' > ' + input_working_dir + '/temp_panTEs_sorted.txt'
    subprocess.call(cmd, shell=True)

    store_final_line_list = []
    store_c_s_type_line_list = []
    with open(input_working_dir + '/temp_panTEs_sorted.txt', 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            te_type = col[6]  ##c or s or o
            if te_type == 'c' or te_type == 's':
                store_c_s_type_line_list.append(eachline)
            else:
                store_final_line_list.append(eachline)

    line_count = 0
    store_all_line_list = []
    line_list = []
    final_line_count = len(store_c_s_type_line_list)
    store_all_line_dic = {} ##key is the location information and value is line from col[3] to end
    ##this dic will be used for the further line extraction
    for eachline in store_c_s_type_line_list:

        line_count += 1
        col = eachline.strip().split()
        current_chr = col[0]
        current_start = col[1]
        current_end = col[2]
        current_id = current_chr + '_' + current_start + '_' + current_end

        other_line_str = ''
        for i in range(3,len(col)):
            other_line_str = other_line_str + '\t' + col[i]
        store_all_line_dic[current_id] = other_line_str

        store_line_dic = {'chr': current_chr, 'st': current_start, 'ed': current_end}

        if line_count == 1:
            line_list.append(store_line_dic)
        else:
            pre_line_dic = line_list[-1]
            pre_chr = pre_line_dic['chr']
            pre_start = pre_line_dic['st']
            pre_ed = pre_line_dic['ed']

            if pre_chr == current_chr:
                if int(current_start) <= int(pre_ed) and int(current_end) >= int(pre_start):
                    line_list.append(store_line_dic)
                    if line_count == final_line_count:
                        store_all_line_list.append(line_list)

                else:
                    store_all_line_list.append(line_list)
                    line_list = []
                    line_list.append(store_line_dic)
                    if line_count == final_line_count:
                        store_all_line_list.append(line_list)
            else:
                store_all_line_list.append(line_list)
                line_list = []
                line_list.append(store_line_dic)
                if line_count == final_line_count:
                    store_all_line_list.append(line_list)

    ##analyze on each list to find the start and end
    for eachlist in store_all_line_list:
        first_dic = eachlist[0]
        last_dic = eachlist[-1]
        final_st = first_dic['st']
        final_ed = last_dic['ed']
        final_chr = first_dic['chr']
        final_loc =  final_chr + '_' + final_st + '_' + final_ed

        if final_loc in store_all_line_dic:
            final_line = final_chr + '\t' + final_st + '\t' + final_ed + store_all_line_dic[final_loc]
            store_final_line_list.append(final_line)
        else:
            final_line = final_chr + '\t' + final_st + '\t' + final_ed + '\t' + 'combine' + '\t' + 'combine' + '\t' + 'combine' + '\t' + 'c' + '\t' + 'Sample_TE_infor;combine'
            store_final_line_list.append(final_line)

    return (store_final_line_list)



##first step: for each location bam to bed and combine all individuals
##we add the final output dir that directly write something on the output
##udpation 122120 remove the c_rg_thr
def genos_construct (input_panTEs_file,
                     input_bam_fl_dir,
                     rm_pcr_hisat2_bam_dir,
                     input_output_dir,
                     temp_split_opt_dir,
                     sample_num,clrd_thr,TSD_thr,t_cover_thr,Her_thd,Acc_thd,comb_close_thr,s_rg_thr,pa_num):

    print('Genotyping')

    ##generate list of dir for the parallele analysis
    temp_loc_line_dir_list = glob.glob(input_output_dir + '/temp_loc_line_dir/*')

    #temp_changing_dir_list = glob.glob(input_output_dir + '/temp_store_changing_te_dir/*')
    #temp_blast_dir_list = glob.glob(input_output_dir + '/temp_store_blast_index_dir/*')


    ##store the location information
    ##this step will be changed when we analyze in a pipeline
    store_panTE_loc_dic = {}
    loc_id = 0
    with open (input_panTEs_file,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            loc_id += 1
            store_panTE_loc_dic[str(loc_id)] = {'chr':col[0],'st':col[1],'ed':col[2],'type':col[6],'sample_str':col[7]}
            ##col[6] is s o or c

    ##updation 122120
    ##store a new version of panTE with only containing c or s and the id is consective
    store_panTE_cs_loc_dic = {}
    cs_id = 0
    with open (input_panTEs_file,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            if col[6] == 'c' or col[6] == 's':
                cs_id += 1
                ##chr1_100_101_s
                cs_loc = col[0] + '_' + col[1] + '_' + col[2] + '_' + col[6]
                store_panTE_cs_loc_dic[str(cs_id)] = {'cs_loc':cs_loc}
    last_cs_id = cs_id


    loc_count = 0
    ##open a file that will write lines to it
    temp_te_genotype_fl = open(temp_split_opt_dir + '/temp_te_genotype.txt', 'w')


    #################
    ##udpation 122120
    ##we need to store the all the potential connected s and c with different less than s_rg_thr
    ##store all the information in the dic_list
    ##next we need to assign the loc information for each location that will be directly used for the next loop
    cs_id = 0
    dic_te = {}
    dic_list = []

    for eachloc_id in store_panTE_loc_dic:
        chr = store_panTE_loc_dic[eachloc_id]['chr']
        bg = store_panTE_loc_dic[eachloc_id]['st']
        ed = store_panTE_loc_dic[eachloc_id]['ed']
        loc_type = store_panTE_loc_dic[eachloc_id]['type']
        comb_type = store_panTE_loc_dic[eachloc_id]['type'] ##add comb_type since the original key has this information
        line = 'none' ##we have no information for the line so we use the line to represent

        if loc_type == 'c' or loc_type == 's':
            cs_id += 1
            if cs_id == 1:  ##if the id is 1, it will directly store in the dic_te
                dic_te[str(cs_id)] = {'chr': chr, 'begin': bg, 'end': ed, 'comb_type': comb_type, 'line': line}
                if cs_id == last_cs_id:  ##should be out of the previous loop
                    dic_list.append(dic_te)

            else:  ##if the id is over 1
                ##if the chr is the same as the previous one
                if dic_te[str(int(cs_id) - 1)]['chr'] == chr:


                    #if dic_te[str(int(cs_id) - 1)]['comb_type'] == comb_type:  ##if the comb_type is the same
                    #pre_st = dic_te[str(int(cs_id) - 1)]['begin']
                    pre_ed = dic_te[str(int(cs_id) - 1)]['end']

                    diff_cs = int(bg) - int(pre_ed)
                    ##we add safe base pair 2
                    if diff_cs <= (int(s_rg_thr) * 2 + 2):
                        dic_te[str(cs_id)] = {'chr': chr, 'begin': bg, 'end': ed, 'comb_type': comb_type, 'line': line}
                    else:
                        ##directly store the previous dic and assign a new dic id
                        dic_list.append(dic_te)
                        dic_te = {}
                        dic_te[str(cs_id)] = {'chr': chr, 'begin': bg, 'end': ed, 'comb_type': comb_type, 'line': line}

                    if cs_id == last_cs_id:  ##should be out of the previous loop
                        dic_list.append(dic_te)

                else:  ##if the chr is not the same as the previous one
                    ##store the dic_te which has been stored the te informations
                    dic_list.append(dic_te)
                    dic_te = {}
                    dic_te[str(cs_id)] = {'chr': chr, 'begin': bg, 'end': ed, 'comb_type': comb_type, 'line': line}

                    if cs_id == last_cs_id:  ##should be out of the previous loop
                        dic_list.append(dic_te)

    ##debug
    ##there is no problem for the dic_list
    #print(dic_list)


    ##now we will generate the loc_line information for each location:  loc_chr + '_' + loc_st + '_' + loc_ed + '_' + loc_type
    ##no problem for the following scripts
    store_cs_loc_locline_dic = {}
    for eachdic_te in dic_list:

        ##sort the id in the eachdic_te
        #wrong eachdic_te_sort = dict(sorted(eachdic_te.items()))
        eachdic_te_sort = {int(k): v for k, v in eachdic_te.items()}


        #print(eachdic_te_sort)

        ##it means there is cs TEs are connected
        if len(list(eachdic_te_sort.keys())) != 1:

            first_cs_id = list(eachdic_te_sort.keys())[0]
            last_cs_id = list(eachdic_te_sort.keys())[-1]

            #print('first id is ' + first_cs_id)
            #print('last is is ' + last_cs_id)


            for eachcs_id in eachdic_te_sort:

                loc_chr = eachdic_te_sort[eachcs_id]['chr']
                loc_st = eachdic_te_sort[eachcs_id]['begin']
                loc_ed = eachdic_te_sort[eachcs_id]['end']
                loc_type = eachdic_te_sort[eachcs_id]['comb_type']
                loc_nm = loc_chr + '_' + loc_st + '_' + loc_ed + '_' + loc_type

                if str(eachcs_id) == str(first_cs_id):

                    #print('analyze first id')

                    next_cs_id = list(eachdic_te_sort.keys())[1]
                    next_st = eachdic_te_sort[next_cs_id]['begin']

                    diff_cs = int(next_st) - int(loc_ed)
                    s_rg_down_thr = int(diff_cs / 2) - 1
                    loc_line = loc_chr + ':' + str(int(loc_st) - int(s_rg_thr)) + '-' + str(int(loc_ed) + int(s_rg_down_thr))

                    store_cs_loc_locline_dic[loc_nm] = loc_line

                    #print(store_cs_loc_locline_dic)

                else:
                    ##if the id is not the first one and not the last one
                    if eachcs_id != last_cs_id:
                        curr_index = list(eachdic_te_sort.keys()).index(eachcs_id)
                        pre_cs_id = list(eachdic_te_sort.keys())[int(curr_index) - 1]
                        next_cs_id = list(eachdic_te_sort.keys())[int(curr_index) + 1]

                        pre_ed = eachdic_te_sort[pre_cs_id]['end']
                        next_st = eachdic_te_sort[next_cs_id]['begin']

                        up_diff_cs = int(loc_st) - int(pre_ed)
                        s_rg_up_thr = int(up_diff_cs / 2) - 1

                        down_diff_cs = int(next_st) - int(loc_ed)
                        s_rg_down_thr = int(down_diff_cs / 2) - 1

                        loc_line = loc_chr + ':' + str(int(loc_st) - int(s_rg_up_thr)) + '-' + str(
                            int(loc_ed) + int(s_rg_down_thr))

                        store_cs_loc_locline_dic[loc_nm] = loc_line

                        #print(store_cs_loc_locline_dic)

                    else:
                        #print('analyze last id')

                        ##if the id is the last one
                        pre_cs_id = list(eachdic_te_sort.keys())[-2]
                        pre_ed = eachdic_te_sort[pre_cs_id]['end']

                        up_diff_cs = int(loc_st) - int(pre_ed)
                        s_rg_up_thr = int(up_diff_cs / 2) - 1

                        loc_line = loc_chr + ':' + str(int(loc_st) - int(s_rg_up_thr)) + '-' + str(int(loc_ed) + int(s_rg_thr))
                        store_cs_loc_locline_dic[loc_nm] = loc_line

                        #print(store_cs_loc_locline_dic)

        else:
            ##if only single TE
            eachcs_id = list(eachdic_te_sort.keys())[0]
            loc_chr = eachdic_te_sort[eachcs_id]['chr']
            loc_st = eachdic_te_sort[eachcs_id]['begin']
            loc_ed = eachdic_te_sort[eachcs_id]['end']
            loc_type = eachdic_te_sort[eachcs_id]['comb_type']
            loc_nm = loc_chr + '_' + loc_st + '_' + loc_ed + '_' + loc_type

            loc_line = loc_chr + ':' + str(int(loc_st) - int(s_rg_thr)) + '-' + str(int(loc_ed) + int(s_rg_thr))
            store_cs_loc_locline_dic[loc_nm] = loc_line

    #print('total_dic is ')
    #print(store_cs_loc_locline_dic)



    ###########################
    ##now we conduct genotyping
    for eachloc_id in store_panTE_loc_dic:

        loc_chr = store_panTE_loc_dic[eachloc_id]['chr']
        loc_st = store_panTE_loc_dic[eachloc_id]['st']
        loc_ed = store_panTE_loc_dic[eachloc_id]['ed']
        loc_type = store_panTE_loc_dic[eachloc_id]['type']

        ##this is the temp dir we will delete it later
        opt_loc_line_dir = temp_loc_line_dir_list[int(pa_num)]
        # opt_loc_line_dir = input_output_dir + '/temp_' + loc_line
        # if not os.path.exists(opt_loc_line_dir):
        #    os.makedirs(opt_loc_line_dir)

        ##generate a temp bam for the target locus
        opt_loc_line_sam_dir = opt_loc_line_dir + '/loc_sam_dir'
        if not os.path.exists(opt_loc_line_sam_dir):
            os.makedirs(opt_loc_line_sam_dir)

        ##generate a temp bed dir to store all the output bed
        opt_loc_line_bam2bed_dir = opt_loc_line_dir + '/bam2bed_dir'
        if not os.path.exists(opt_loc_line_bam2bed_dir):
            os.makedirs(opt_loc_line_bam2bed_dir)

        if loc_type == 'c' or loc_type == 's':

            ##updation 122020
            sample_te_string = store_panTE_loc_dic[eachloc_id]['sample_str']
            ##if the loc_type is c or s
            store_te_name_dic = {}
            sample_te_list = sample_te_string.split(';')
            for nm_count in range(1,len(sample_te_list)):
                mt = re.match('(.+)\:(.+)', sample_te_list[nm_count])
                te_nm = mt.group(2)
                store_te_name_dic[te_nm] = 1
            ##generate tename str
            te_nm_list = []
            for eachtenm in store_te_name_dic:
                te_nm_list.append(eachtenm)
            te_nm_str = ','.join(te_nm_list)


            ##udpation 122120
            loc_nm = loc_chr + '_' + loc_st + '_' + loc_ed + '_' + loc_type
            loc_line = store_cs_loc_locline_dic[loc_nm]
            ##wait for the search range
            ##updation 121920
            #if loc_type == 'c':
            #    loc_line = loc_chr + ':' + str(int(loc_st) - int(c_rg_thr)) + '-' + str(int(loc_ed) + int(c_rg_thr))
            #else:
            #    loc_line = loc_chr + ':' + str(int(loc_st) - int(s_rg_thr)) + '-' + str(int(loc_ed) + int(s_rg_thr))
            #print(loc_line)


            ##analyze each bam file
            bed_fl_str = ''
            bam_fl_list = glob.glob(input_bam_fl_dir + '/*.bam')
            for each_bam_fl in bam_fl_list:
                ##get the fl name
                if '_2' in each_bam_fl:
                    mt = re.match('.+/(.+)_2\.bam', each_bam_fl)
                    fl_nm = mt.group(1)
                else:
                    mt = re.match('.+/(.+)\.bam', each_bam_fl)
                    fl_nm = mt.group(1)

                ##we need to extract pan loc first
                cmd = 'samtools view -h -bS ' + each_bam_fl + ' ' + loc_line + ' > ' + \
                      opt_loc_line_sam_dir + '/temp.bam'   ##generate a new dir to store the temp.sam file
                subprocess.call(cmd, shell=True)

                ##transfer the target bam to bed file
                cmd = 'bedtools bamtobed -i ' + opt_loc_line_sam_dir + '/temp.bam' + ' -cigar > ' + opt_loc_line_bam2bed_dir + '/' + fl_nm + '.bed'
                subprocess.call(cmd,shell=True)

                bed_fl_str = bed_fl_str + ' ' + opt_loc_line_bam2bed_dir + '/' + fl_nm + '.bed'


            ##combine all the bed together
            cmd = 'cat ' + bed_fl_str + ' > ' + opt_loc_line_dir + '/temp_combined_bam2bed_fl.bed'
            subprocess.call(cmd,shell=True)

            ##sort the bed fl
            cmd = 'sort -k1,1V -k2,2n ' + opt_loc_line_dir + '/temp_combined_bam2bed_fl.bed > ' +  opt_loc_line_dir + '/temp_combined_bam2bed_fl_sorted.bed'
            subprocess.call(cmd,shell=True)

            ##calculate the number of st clipped reads
            store_clipped_st_reads_dic = {} ##key is the 4443026 and value is the clipped count
            store_clipped_ed_reads_dic = {}
            with open (opt_loc_line_dir + '/temp_combined_bam2bed_fl_sorted.bed','r') as ipt:
                for eachline in ipt:
                    eachline = eachline.strip('\n')
                    col = eachline.strip().split()
                    st_loc = col[1]
                    ed_loc = col[2]
                    cigar = col[6]

                    if 'S' in cigar:


                        ##updation 122220
                        ##we need to consider whether the clipped position is the right position
                        if not re.match('^\d+S.+\d+S$',cigar):
                            if re.match('.+\d+S$',cigar):
                                ##it means the clipped read is in the end
                                ##for the ed ver
                                if ed_loc in store_clipped_ed_reads_dic:
                                    store_clipped_ed_reads_dic[ed_loc] += 1
                                else:
                                    store_clipped_ed_reads_dic[ed_loc] = 1

                            if re.match('^\d+S.+',cigar):
                                ##it means the clipped read is in the start
                                if st_loc in store_clipped_st_reads_dic:
                                    store_clipped_st_reads_dic[st_loc] += 1
                                else:
                                    store_clipped_st_reads_dic[st_loc] = 1

                            ##for the st ver
                            #if st_loc in store_clipped_st_reads_dic:
                            #    store_clipped_st_reads_dic[st_loc] += 1
                            #else:
                            #    store_clipped_st_reads_dic[st_loc] = 1
                            ##for the ed ver
                            #if ed_loc in store_clipped_ed_reads_dic:
                            #    store_clipped_ed_reads_dic[ed_loc] += 1
                            #else:
                            #    store_clipped_ed_reads_dic[ed_loc] = 1

            ##generate a new file for the clipped reads number information
            clrd_total_thr = int(int(clrd_thr) * int(sample_num) * 2 * 0.05)
            #clrd_total_thr = 3
            with open (opt_loc_line_dir + '/temp_clipped_reads_flt_num.txt','w+') as opt:
                for eachloc in store_clipped_st_reads_dic:
                    if store_clipped_st_reads_dic[eachloc] >= clrd_total_thr:
                        opt.write('loc' + '\t' + eachloc + '\t' + str(store_clipped_st_reads_dic[eachloc]) + '\t' + 'st' + '\n')
                for eachloc in store_clipped_ed_reads_dic:
                    if store_clipped_ed_reads_dic[eachloc] >= clrd_total_thr:
                        opt.write('loc' + '\t' + eachloc + '\t' + str(store_clipped_ed_reads_dic[eachloc]) + '\t' + 'ed' + '\n')

            ##sort the file
            cmd = 'sort -k1,1V -k2,2n ' + opt_loc_line_dir + '/temp_clipped_reads_flt_num.txt > ' +  opt_loc_line_dir + '/temp_clipped_reads_flt_num_sorted.txt'
            subprocess.call(cmd,shell=True)

            ##decide the location
            store_final_loc_line_list = []

            total_line_count = 0
            with open (opt_loc_line_dir + '/temp_clipped_reads_flt_num_sorted.txt','r') as ipt:
                for eachline in ipt:
                    total_line_count += 1

            store_st_ed_dic_list = []
            store_st_ed_dic = {}
            line_count = 0
            with open (opt_loc_line_dir + '/temp_clipped_reads_flt_num_sorted.txt','r') as ipt:
                for eachline in ipt:
                    eachline = eachline.strip('\n')
                    col = eachline.strip().split()
                    current_loc_type = col[3] ##st or ed
                    current_loc = col[1] ##location information eg. 305
                    line_count += 1

                    store_st_ed_dic[current_loc] = current_loc_type

                    if line_count == 1:
                        store_st_ed_dic_list.append(current_loc)

                        ##updation 122120
                        if line_count == total_line_count:
                            final_loc_line = loc_chr + '\t' + current_loc
                            store_final_loc_line_list.append(final_loc_line)

                    else:

                        last_loc = store_st_ed_dic_list[-1]
                        last_loc_type = store_st_ed_dic[last_loc]

                        if last_loc_type == 'st':
                            ##we will compare the last loc and the current loc
                            if current_loc_type == 'ed':
                                ##check if they are overlapped
                                if (int(current_loc) - int(last_loc)) > 0 and (int(current_loc) - int(last_loc)) <= int(TSD_thr):
                                    ##store the final loc
                                    final_loc_line = loc_chr + '\t' + last_loc
                                    store_final_loc_line_list.append(final_loc_line)
                                    ##we do not append this current loc to the store_st_ed_dic_list
                                    ##since this current_loc is ed and we choose the last st loc to go to the final output
                                else:
                                    final_loc_line = loc_chr + '\t' + last_loc
                                    store_final_loc_line_list.append(final_loc_line)

                                    if line_count == total_line_count:
                                        final_loc_line = loc_chr + '\t' + current_loc
                                        store_final_loc_line_list.append(final_loc_line)
                                    else:
                                        ##we need to stoer the current_loc that is 'ed'
                                        store_st_ed_dic_list.append(current_loc)

                            else:
                                final_loc_line = loc_chr + '\t' + last_loc
                                store_final_loc_line_list.append(final_loc_line)

                                if line_count == total_line_count:
                                    final_loc_line = loc_chr + '\t' + current_loc
                                    store_final_loc_line_list.append(final_loc_line)
                                else:
                                    ##we need to stoer the current_loc since it is 'st'
                                    store_st_ed_dic_list.append(current_loc)

                        else:
                            final_loc_line = loc_chr + '\t' + last_loc
                            store_final_loc_line_list.append(final_loc_line)
                            ##we need to stoer the current_loc that is 'st' or 'ed'
                            if line_count == total_line_count:
                                final_loc_line = loc_chr + '\t' + current_loc
                                store_final_loc_line_list.append(final_loc_line)
                            else:
                                store_st_ed_dic_list.append(current_loc)

            with open (opt_loc_line_dir + '/opt_decide_loc.txt','w+') as opt:
                for eachline in store_final_loc_line_list:
                    opt.write(eachline + '\n')

            ##we need to set a unique for the opt_decide_loc_fl since the above function will generate duplicate location information
            ##it does not matter
            cmd = 'sort -u ' + opt_loc_line_dir + '/opt_decide_loc.txt > ' +  opt_loc_line_dir + '/opt_decide_loc_uniq.txt'
            subprocess.call(cmd,shell=True)

            ##we need to modify the opt_decide file that allows to not be connected
            store_loc_infor_list = []
            store_final_line_list = []
            count = 0
            with open(opt_loc_line_dir + '/opt_decide_loc_uniq.txt', 'r') as ipt:
                for eachline in ipt:
                    eachline = eachline.strip('\n')
                    col = eachline.strip().split()
                    count += 1
                    current_start = col[1]
                    if count == 1:
                        store_loc_infor_list.append(current_start)
                        store_final_line_list.append(eachline)
                    else:
                        pre_start = store_loc_infor_list[-1]
                        if int(current_start) - int(pre_start) > int(comb_close_thr):
                            store_final_line_list.append(eachline)
                            store_loc_infor_list.append(current_start)
                        else:
                            ##it means these two locations are too close, we need to select the first one
                            store_loc_infor_list.append(current_start)

            with open (opt_loc_line_dir + '/opt_decide_loc_uniq_flt.txt','w+') as opt:
                for eachline in store_final_line_list:
                    opt.write(eachline + '\n')

            #############################
            ##step 2 check the genotyping
            ##use the loc in the opt_decide_loc fl to check the genotype
            ##analyze each bam file
            ##generate a temp bed dir to store all the output bed
            temp_check_geno_sam_dir = opt_loc_line_dir + '/temp_check_geno_sam_dir'
            if not os.path.exists(temp_check_geno_sam_dir):
                os.makedirs(temp_check_geno_sam_dir)


            with open (opt_loc_line_dir + '/opt_decide_loc_uniq_flt.txt','r') as ipt:
                for eachloc_line in ipt:

                    loc_count += 1
                    print('genotype te count ' + str(loc_count))

                    eachloc_line = eachloc_line.strip('\n')
                    eachloc_line_col = eachloc_line.strip().split()
                    chr = eachloc_line_col[0]
                    start = eachloc_line_col[1]

                    ##we do not use start - 1 since there is one case that shows 100M at the start -1 location
                    loc_line = chr + ':' + str(int(start)) + '-' + str(int(start) + 1)

                    fl_nm_string = ''
                    bam_fl_list = glob.glob(input_bam_fl_dir + '/*.bam')
                    for each_bam_fl in bam_fl_list:

                        ##get the fl name
                        if '_2' in each_bam_fl:
                            mt = re.match('.+/(.+)_2\.bam', each_bam_fl)
                            fl_nm = mt.group(1)
                        else:
                            mt = re.match('.+/(.+)\.bam', each_bam_fl)
                            fl_nm = mt.group(1)

                        ##use samtools view to extract loc_line information that save as temp.sam file
                        ##updation 8.24
                        cmd = 'samtools view ' + each_bam_fl + ' ' + loc_line + ' > ' + \
                              temp_check_geno_sam_dir + '/temp.sam'  ##generate a new dir to store the temp.sam file
                        subprocess.call(cmd, shell=True)

                        ##detect the H information: H means clipped read in the sam file
                        total_num = 0
                        H_num = 0
                        ##updation 8.23
                        with open(temp_check_geno_sam_dir + '/temp.sam' , 'r') as ipt:
                            for each_cover_line in ipt:
                                total_num += 1
                                each_cover_line = each_cover_line.strip('\n')
                                col_tmp_sam = each_cover_line.strip().split()
                                CIGAR = col_tmp_sam[5]
                                if 'S' in CIGAR:
                                    H_num += 1

                        if int(total_num) != 0:

                            ##updation 101720 add the total number threshold
                            if int(total_num) >= int(t_cover_thr):

                                H_pro = float(H_num / total_num)

                                ##this binom_test will be atlernative choice for users
                                ##we will determine if the genotyping is biased towards the clipped reads
                                ##H0: π <= 0.5 (the genotyping is not biased towards the clipped reads) p value >= 0.05 it will be hetero or homo
                                ##Ha π > 0.5 (the genotyping is biased towards the clipped reads) p value < 0.05 it will be homo
                                ##if n is larger than 16, so we need to use less
                                ##to have an alternative test.
                                ##if we accept the alternative test the p should below 0.05 others we need to think no differences as 0.5 and it would be hetero
                                #p = binom_test(15, n=int(total_num), p=0.5, alternative='less')
                                if (1 - float(Her_thd)) <= float(H_pro) and float(H_pro) <= float(Her_thd):
                                    fl_nm_string = fl_nm_string + '\t' + fl_nm + ':' + str(H_pro) + ';' + \
                                                   str(int(H_num)) + '/' + str(total_num) + ';' + '0/1'

                                ##the same as the consensuse TEs, but different as reference this will be the 1/1
                                if float(Her_thd) < float(H_pro) and float(H_pro) <= 1.0:
                                    fl_nm_string = fl_nm_string + '\t' + fl_nm + ':' + str(H_pro) + ';' + \
                                                   str(int(H_num)) + '/' + str(total_num) + ';' + '1/1'
                                ##different as the consensus TEs, but same as the reference this will be the 0/0
                                if 0 <= float(H_pro) < (1 - float(Her_thd)):
                                    fl_nm_string = fl_nm_string + '\t' + fl_nm + ':' + str(H_pro) + ';' + \
                                                   str(int(H_num)) + '/' + str(total_num) + ';' + '0/0'

                    if fl_nm_string != '':
                        ##updation 122020 add the te_nm_str that only contains the TEname eg. ClassI_LTR_ltr_Others_TE967,ClassI_LTR_Gypsy_TE1672
                        final_line = chr + '\t' + start + '\t' + str(int(start) + 1) + '\t' + loc_type + '\t' + te_nm_str + fl_nm_string
                        temp_te_genotype_fl.write(final_line + '\n')
                        #store_final_line_list.append(final_line)

            ##debug
            ##remove the opt_loc_line_dir after all the analysis
            cmd = 'rm -rf ' + opt_loc_line_dir
            subprocess.call(cmd,shell=True)

        else:

            sample_te_string = store_panTE_loc_dic[eachloc_id]['sample_str']
            ##if the loc_type is o
            sample_te_list = sample_te_string.split(';')
            first_sp_te = sample_te_list[1]
            mt = re.match('(.+)\:(.+)', first_sp_te)
            te_nm = mt.group(2)

            ##This will divide into two cases
            ##Case A
            ##the read has information for the N
            ##extract the bam file with start
            start = loc_st
            end = loc_ed

            loc_line = loc_chr + ':' + str(int(start) - 1) + '-' + str(int(start) + 1)

            ##generate a string of fl name information that will be added to the input_filter_te_loc file
            ##example: sample1:XX   sample2:XX ....   XX means depth
            fl_nm_string = ''

            ##analyze each bam file from the hisat2 software output

            ##NOTE:the sample name should be paid attention

            bam_fl_list = glob.glob(rm_pcr_hisat2_bam_dir + '/*.bam')

            for each_bam_fl in bam_fl_list:

                ##get the fl name  ##attention for this name information
                fl_nm = ''
                if '_2' in each_bam_fl:
                    mt = re.match('.+/(.+)_2\.bam', each_bam_fl)
                    fl_nm = mt.group(1)
                else:
                    mt = re.match('.+/(.+)\.bam', each_bam_fl)
                    fl_nm = mt.group(1)

                ##use samtools view to extract loc_line information that save as temp.sam file
                ##generate a new dir to store the temp.sam file
                ##updation 8.23
                cmd = 'samtools view ' + each_bam_fl + ' ' + loc_line + ' > ' + opt_loc_line_sam_dir + '/temp.sam'
                #print(cmd)
                subprocess.call(cmd, shell=True)

                ##store each read information. key is the read nm and value is the marker
                store_read_dic = {}

                ##initiate a number to store all the valid read with N
                Total_TE_minus_cover_num = 0
                ##initiate a number to store the clipped read without N
                # clip_cover_num_noN = 0
                ##initiate a number to store final right read
                Total_TE_plus_cover_num = 0

                ##initiate a string to store the string of total_cover_read and right_cover_read
                total_string = ''
                cover_string = ''

                ##updation4.16
                ##initiate a dic to store number big gap situation
                ##this dic use the st_read as key and the value is a list of dic_list and output of this loop, we will
                ##calcualte the how many same N in this value part
                big_gap_dic = {}

                ##updation 4.16
                ##initiate a list to store the dic_list that will be used for the big gap analysis
                # dic_list_list = []

                ##updation 8.23
                with open(opt_loc_line_sam_dir + '/temp.sam', 'r') as ipt_sam:

                    for eachline_sam in ipt_sam:

                        if not eachline_sam.startswith('@'):

                            eachline_sam = eachline_sam.strip('\n')
                            col_sam = eachline_sam.strip().split()
                            st_read = col_sam[3]

                            ##add the match reads number to the start of read, and if the results are close the insert position +-1
                            CIGAR = col_sam[5]

                            if CIGAR != '*':

                                ##extract the number of D
                                list_CIGAR = re.findall('\d+|\D+', CIGAR)
                                n = int(len(list_CIGAR) / 2)

                                dic_list = []
                                for i in range(0, int(n)):
                                    dic = {list_CIGAR[(2 * i + 1)]: int(list_CIGAR[(2 * i)])}
                                    dic_list.append(dic)

                                ##do not consider the S and N
                                ##do not consider the longest N
                                ##For the TE+, we will calculate all the ranges of N, and detect their accuracy
                                ##For the TE-, we will detect all the M regions, and detect if these regions cover the start region

                                ##divide two situations
                                ##first: N is in the CIGAR
                                ##second: N is not in the CIGAR

                                ###################
                                ##first situations:
                                if 'N' in CIGAR:

                                    ##identify the range of M and N,
                                    ##for the M, detect if M cover the start
                                    ##for the N, detect accuracy of the N

                                    ##initiate a list to store all the range information (dic) for each read
                                    store_range_list = []

                                    ##initate a start point for the letters
                                    letter_start = st_read

                                    for eachdic in dic_list:

                                        letter = list(eachdic.keys())[0]
                                        letter_len = eachdic[list(eachdic.keys())[0]]

                                        ##do not consider the I P and S
                                        if letter != 'I' and letter != 'P' and letter != 'S':
                                            ##initiate a dic to store the range information
                                            ##the key is the letter and the value is the range
                                            store_letter_dic = {}

                                            store_letter_dic[letter] = {'bg': int(letter_start), 'ed': 0}

                                            letter_start = int(letter_start) + int(letter_len)

                                            store_letter_dic[letter]['ed'] = int(letter_start)

                                            ##store the dic to the store_range_list
                                            store_range_list.append(store_letter_dic)

                                    ################################
                                    ##detect if N is accurate or not
                                    N_TE_plus_count = 0
                                    N_TE_minus_count = 0
                                    for eachdic in store_range_list:

                                        if list(eachdic.keys())[0] == 'N':

                                            letter_bg = eachdic['N']['bg']
                                            letter_ed = eachdic['N']['ed']

                                            ##updation 110620 recover this case
                                            ##updation 4.16
                                            ##consider whether the reads are partially covered or all covered.
                                            ##if N covers all the reference read
                                            ##in this situation, we will detect the accuracy and if not fit the accuracy,
                                            ##we will consider how many reads showing similar bp in the N gap,
                                            ##set a parameter to control the read number
                                            if int(letter_bg) < int(start) and int(letter_ed) > int(end):

                                                ##whether fit the accuracy
                                                union_left = int(letter_bg)
                                                inter_left = int(start)

                                                union_right = int(letter_ed)
                                                inter_right = int(end)

                                                union_len = union_right - union_left
                                                inter_len = inter_right - inter_left + 1  ##forbid to be zero

                                                BIAS = union_len / inter_len
                                                Accuracy = 1 / BIAS

                                                ##get the right cover
                                                if float(Accuracy) >= float(Acc_thd):
                                                    # cover_string = cover_string + ';' + st_read
                                                    N_TE_plus_count += 1

                                                ##if N does not fit the cover
                                                else:
                                                    ##consider how many reads show the similar N bp gap
                                                    ##store the big N gap for a dic, and we will calculate the number of
                                                    ##similar gap
                                                    if st_read in big_gap_dic:

                                                        big_gap_dic[st_read]['read_count'] += 1

                                                        big_gap_dic[st_read]['dic_list_list'].append(dic_list)

                                                    else:
                                                        dic_list_list = []

                                                        dic_list_list.append(dic_list)

                                                        big_gap_dic[st_read] = {'read_count': 1,
                                                                                'dic_list_list': dic_list_list}


                                            ##if N do not cover all the reference read,
                                            ## it will detect TE+ or TE- showing as follows
                                            else:

                                                ##for the union_left and inter_left
                                                if int(start) >= int(letter_bg):
                                                    union_left = int(letter_bg)
                                                    inter_left = int(start)
                                                else:
                                                    union_left = int(start)
                                                    inter_left = int(letter_bg)

                                                ##for the union_right and inter_right
                                                if int(end) >= int(letter_ed):
                                                    union_right = int(end)
                                                    inter_right = int(letter_ed)
                                                else:
                                                    union_right = int(letter_ed)
                                                    inter_right = int(end)

                                                union_len = union_right - union_left
                                                inter_len = inter_right - inter_left

                                                ##forbid to be zero
                                                if inter_len == 0:

                                                    inter_len_new = inter_len + 1

                                                    BIAS = union_len / inter_len_new
                                                    Accuracy = 1 / BIAS

                                                    ##get the right cover
                                                    if float(Accuracy) >= float(Acc_thd):
                                                        # cover_string = cover_string + ';' + st_read
                                                        N_TE_plus_count += 1
                                                else:
                                                    BIAS = union_len / inter_len
                                                    Accuracy = 1 / BIAS

                                                    ##get the right cover
                                                    if float(Accuracy) >= float(Acc_thd):
                                                        # cover_string = cover_string + ';' + st_read
                                                        N_TE_plus_count += 1

                                    ##if the N do not fit the length of TE,
                                    ##the second step is to detect if it belongs to the TE-
                                    if N_TE_plus_count == 0:

                                        for eachdic in store_range_list:
                                            if list(eachdic.keys())[0] == 'M':
                                                letter_bg = eachdic['M']['bg']
                                                letter_ed = eachdic['M']['ed']
                                                if letter_bg < int(start) and letter_ed > int(start):
                                                    N_TE_minus_count += 1

                                    ##detect if the TE is TE+ or TE-
                                    if N_TE_plus_count != 0:
                                        Total_TE_plus_cover_num += 1
                                    if N_TE_minus_count != 0:
                                        Total_TE_minus_cover_num += 1

                                ##the second situation
                                if 'N' not in CIGAR:
                                    total_order_num = 0
                                    # sed_count = 0
                                    for eachdic in dic_list:
                                        # sed_count += 1
                                        # if sed_count == len(dic_list):
                                        # if list(dic_list[-1].keys())[0] == 'S':  ##detect the last one as S
                                        #    break
                                        # else:
                                        if list(eachdic.keys())[0] != 'I' and list(eachdic.keys())[0] != 'P' \
                                                and list(eachdic.keys())[0] != 'S':
                                            total_order_num = total_order_num + \
                                                              int(eachdic[list(eachdic.keys())[0]])

                                    ##get the stop location for the read
                                    stop_loc = int(st_read) + total_order_num

                                    if stop_loc > int(
                                            start):  ##if the stop_loc is over the start, this read means not deletion

                                        total_string = total_string + ';' + st_read

                                        Total_TE_minus_cover_num += 1

                ##updation110620 recover this step
                ##updation4.16
                ##detect the situations in the big_gap_dic
                ##if the key number of the big_gap_dic is over or equal to the gap_cover_th
                ##big_gap_dic stores the key of eachst_nm

                for eachst_nm in big_gap_dic:

                    ##extract the number of start name
                    st_nm_num = big_gap_dic[eachst_nm]['read_count']

                    ##gap_cover_thd = 3 that is the same as the t_cover_thr
                    if st_nm_num >= int(t_cover_thr):

                        ##initiate a dic to store number of N showing same number
                        num_N_dic = {}
                        for eachdic_list in big_gap_dic[eachst_nm]['dic_list_list']:

                            for gap_eachdic in eachdic_list:

                                if str(list(gap_eachdic.keys())[0]) == 'N':

                                    N_num = str(gap_eachdic['N'])

                                    if N_num in num_N_dic:
                                        num_N_dic[N_num] += 1
                                    else:
                                        num_N_dic[N_num] = 1

                        ##calculate the number of key number in the num_N_dic
                        for each_N_num in num_N_dic:

                            if int(num_N_dic[each_N_num]) >= st_nm_num:  ##the number N may be higher than the st_nm_num since other read has similar N.

                                Total_TE_plus_cover_num += int(st_nm_num)

                total_cover_num = Total_TE_minus_cover_num + Total_TE_plus_cover_num

                if int(total_cover_num) != 0:

                    ##use the Her_thd to decide the TE is hetero or homo

                    ##marker the TEs
                    junc_rt = float(Total_TE_plus_cover_num / total_cover_num)
                    if (1 - float(Her_thd)) <= float(junc_rt) and float(junc_rt) < float(Her_thd):
                        fl_nm_string = fl_nm_string + '\t' + fl_nm + ':' + str(junc_rt) + ';' + \
                                       str(Total_TE_plus_cover_num) + '/' + str(
                            total_cover_num) + ';' + '0/1'
                    ##this situation maybe rarely happen
                    ##this should be 1/1 because 0/0 indicates the same as the reference the reference has insertion.
                    if float(junc_rt) >= float(Her_thd) and float(junc_rt) <= 1.0:
                        fl_nm_string = fl_nm_string + '\t' + fl_nm + ':' + str(junc_rt) + ';' + \
                                       str(Total_TE_plus_cover_num) + '/' + str(
                            total_cover_num) + ';' + '1/1'
                    ##If this situation is the same, so we will use 0/0
                    if 0 <= float(junc_rt) and float(junc_rt) < (1 - float(Her_thd)):  ##this should be 0/0
                        fl_nm_string = fl_nm_string + '\t' + fl_nm + ':' + str(junc_rt) + ';' + \
                                       str(Total_TE_plus_cover_num) + '/' + str(
                            total_cover_num) + ';' + '0/0'

            ##add the fl_nm_string to the eachline end to generate the final line
            if fl_nm_string != '':
                final_line = loc_chr + '\t' + loc_st + '\t' + loc_ed + '\t' + loc_type + '\t' + te_nm + fl_nm_string
                temp_te_genotype_fl.write(final_line + '\n')

    temp_te_genotype_fl.close()

##this step will filter MAF < 0.05 location to save time for the location annotation
##in a pipeline we need to check whether the MAF_filtration is equal to 1
##if it is we will conduct the filtration
##if not, we do not initiate this function
#def filtration (opt_te_genotype_fl,input_output_dir):



def getRC(inputseq):
    out = inputseq.reverse_complement()    #saves the reverse compliment to out
    return out

def rc_te_lib (input_te_lib_fl,input_output_dir):

    ##we need to create the + and - TE library sequence
    ##inverse the lib sequence first
    store_rc_seq_lib_dic = {}
    for seq_record in SeqIO.parse(input_te_lib_fl,'fasta'):
        rc_seq = getRC(seq_record.seq)
        store_rc_seq_lib_dic['rc_' + seq_record.id] = str(rc_seq)

    with open (input_output_dir + '/opt_rc_te.lib','w+') as opt:
        for eachname in store_rc_seq_lib_dic:
            opt.write('>' + eachname + '\n' + store_rc_seq_lib_dic[eachname] + '\n')

    ##generate a combined lib
    cmd = 'cat ' + input_te_lib_fl + ' ' + input_output_dir + '/opt_rc_te.lib' + ' > ' + input_output_dir + '/opt_for_res_te.lib'
    subprocess.call(cmd,shell=True)

    ##updation 122020
    ##this is for the all te name library that is used for the unknown TEs after several round annotation
    ##generate a databased for the te lib
    temp_allte_lib_db_dir = input_output_dir + '/temp_allte_lib_db_dir'
    if not os.path.exists(temp_allte_lib_db_dir):
        os.makedirs(temp_allte_lib_db_dir)

    cmd = 'makeblastdb -in ' + input_output_dir + '/opt_for_res_te.lib' + ' -dbtype nucl -out ' + temp_allte_lib_db_dir + '/te_db'
    subprocess.call(cmd,shell=True)


def annotation (temp_split_opt_dir,
                opt_te_genotype_fl,
                temp_lib_db_dir,
                input_bam_fl_dir,
                input_output_dir,
                TE_lib_fl,
                pa_num):

    print('annotation')
    ##store the te genotype information
    temp_blast_dir_list = glob.glob(input_output_dir + '/temp_store_blast_index_dir/*')
    temp_sam_dir_list = glob.glob(input_output_dir + '/temp_multi_sam_save_dir/*')

    ##generate a temp bed dir to store all the output bed
    #temp_sam_dir = input_output_dir + '/temp_sam_dir'
    #if not os.path.exists(temp_sam_dir):
    #    os.makedirs(temp_sam_dir)

    temp_sam_dir = temp_sam_dir_list[int(pa_num)]
    temp_blast_dir = temp_blast_dir_list[int(pa_num)]
    #temp_blast_dir = input_output_dir + '/temp_blast_dir'
    #if not os.path.exists(temp_blast_dir):
    #    os.makedirs(temp_blast_dir)

    #temp_lib_db = input_output_dir + '/temp_lib_db'
    #if not os.path.exists(temp_lib_db):
    #    os.makedirs(temp_lib_db)

    ##generate a databased for the te lib
    #cmd = 'makeblastdb -in ' + opt_for_res_te_fl + ' -dbtype nucl -out ' + temp_lib_db + '/te_db'
    #subprocess.call(cmd,shell=True)

    ##updation 122020
    ##generate a dir to store the temp library that contains pre TEs name information
    temp_preTE_db_dir_list = glob.glob(input_output_dir + '/temp_preTE_db_dir/*')
    temp_preTE_db_dir = temp_preTE_db_dir_list[int(pa_num)]


    temp_te_genotype_annot_fl = open(temp_split_opt_dir + '/temp_te_genotype_annot.txt', 'w')

    #store_final_annot_te_line_list = []
    store_match_te_report_for_eachloc_line_list = []
    te_count = 0
    with open (opt_te_genotype_fl,'r') as ipt:
        for eachline in ipt:
            te_count += 1
            print('annot te count is ' + str(te_count))

            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            chr = col[0]
            start = col[1]
            end = col[2]
            loc_type = col[3]

            #if start == '28888041':
            #    print('28888041 can be found')

            if loc_type == 'c' or loc_type == 's':

                ##updation 122020
                ##check the TEname for each genotyping line
                te_nm_str = col[4]
                ##check if the te_nm_str has multiple name
                ##if it has contained only one te_nm we will directly regard it as the te_nm
                ##otherwise, we will check the blast results based on the library constructed by these te
                ##if it is still unknown, we will check the blast results based on all the TEs from the library
                te_nm_str_list = te_nm_str.split(',')
                if len(te_nm_str_list) != 1:

                    ##construct the db in the temp_preTE_db_dir
                    ##extract the seq inforamtion
                    store_preTEs_seq_dic = {}
                    for seq_record in SeqIO.parse(TE_lib_fl,'fasta'):
                        if seq_record.id in te_nm_str_list:
                            store_preTEs_seq_dic[seq_record.id] = str(seq_record.seq)

                    with open (temp_preTE_db_dir + '/temp_preTEs.fa','w+') as opt:
                        for eachid in store_preTEs_seq_dic:
                            opt.write('>' + eachid + '\n' + store_preTEs_seq_dic[eachid] + '\n')

                    ##construct the db
                    cmd = 'makeblastdb -in ' + temp_preTE_db_dir + '/temp_preTEs.fa' + ' -dbtype nucl -out ' + temp_preTE_db_dir + '/preTE_db'
                    subprocess.call(cmd, shell=True)

                    loc_line = chr + ':' + start + '-' + end

                    store_cl_rd_dic = {} ##key is the seq id and value is the seq
                    bam_fl_list = glob.glob(input_bam_fl_dir + '/*.bam')
                    seq_id = 0
                    for each_bam_fl in bam_fl_list:
                        cmd = 'samtools view -h ' + each_bam_fl + ' ' + loc_line + ' > ' + \
                              temp_sam_dir + '/temp.sam'   ##generate a new dir to store the temp.sam file
                        subprocess.call(cmd, shell=True)

                        with open(temp_sam_dir + '/temp.sam', 'r') as ipt:
                            for each_cover_line in ipt:

                                if not each_cover_line.startswith('@'):
                                    each_cover_line = each_cover_line.strip('\n')
                                    col_tmp_sam = each_cover_line.strip().split()

                                    if len(col_tmp_sam) > 3:

                                        CIGAR = col_tmp_sam[5]

                                        if 'S' in CIGAR:

                                            ##collect the clipped read sequence information
                                            list_CIGAR = re.findall('\d+|\D+', CIGAR)
                                            n = int(len(list_CIGAR) / 2)
                                            dic_list = []
                                            for i in range(0, int(n)):
                                                dic = {list_CIGAR[(2 * i + 1)]: list_CIGAR[(2 * i)]}
                                                dic_list.append(dic)

                                            ##the first item of the dic_list should be 'S':'number'
                                            if re.match('.+(H|S)$', CIGAR):
                                                S_num = dic_list[-1]['S']
                                            else:
                                                S_num = dic_list[0]['S']

                                            read_string = col_tmp_sam[9]
                                            ##extract the clipped reads
                                            cl_seq = read_string[:int(S_num)]
                                            ##use the cl_seq to blast each TE in the library
                                            ##write out a sequence to temp_blast_dir
                                            ##updation8.24
                                            seq_id += 1
                                            store_cl_rd_dic[str(seq_id)] = cl_seq

                    with open(temp_blast_dir + '/temp_cl_seq.fasta', 'w+') as opt:
                        for eachid in store_cl_rd_dic:
                            opt.write('>' + eachid + '\n' + store_cl_rd_dic[eachid] + '\n')

                    ##updation 122020
                    ##blast to the preTEs_db
                    cmd = 'blastn -task blastn-short -db ' + temp_preTE_db_dir + '/preTE_db' + \
                          ' -query ' + temp_blast_dir + '/temp_cl_seq.fasta -outfmt 6 > ' + \
                          temp_blast_dir + '/final_blast.out'
                    subprocess.call(cmd, shell=True)

                    ##make a blast
                    ##run the blastn -task blastn_short
                    #cmd = 'blastn -task blastn-short -db ' + temp_lib_db_dir + '/te_db' + \
                    #      ' -query ' + temp_blast_dir + '/temp_cl_seq.fasta -outfmt 6 > ' + \
                    #      temp_blast_dir + '/final_blast.out'
                    #subprocess.call(cmd, shell=True)

                    store_final_blast_line_list = []
                    if str(path.isfile(temp_blast_dir + '/final_blast.out')) == 'True':
                        #print(temp_blast_dir + '/final_blast.out' + ' is exist')
                        # if os.path.getsize(temp_blast_dir_list[input_x] + '/final_blast.opt') > 0:
                        #blast_line_count = 0
                        with open(temp_blast_dir + '/final_blast.out') as ipt:
                            for eachline_blast in ipt:
                                eachline_blast = eachline_blast.strip('\n')
                                col_blast = eachline_blast.strip().split('\t')

                                ##only collect the match rate that is equal to 100
                                if int(float(col_blast[2])) == 100:
                                    final_te_score = col_blast[-1]
                                    final_blast_line = col_blast[1] + '\t' + final_te_score
                                    store_final_blast_line_list.append(final_blast_line)
                    else:
                        print('this blast file is not existing')

                    #print(store_final_blast_line_list)

                    if store_final_blast_line_list == []:

                        #final_te_name = 'Unknown'
                        ##updation 122020 the TEs will further blast to the allTE.lib

                        ##run the blastn -task blastn_short
                        cmd = 'blastn -task blastn-short -db ' + temp_lib_db_dir + '/te_db' + \
                              ' -query ' + temp_blast_dir + '/temp_cl_seq.fasta -outfmt 6 > ' + \
                              temp_blast_dir + '/final_blast.out'
                        subprocess.call(cmd, shell=True)

                        ##and further check if the final_blast.out has TEs
                        if str(path.isfile(temp_blast_dir + '/final_blast.out')) == 'True':
                            # print(temp_blast_dir + '/final_blast.out' + ' is exist')
                            # if os.path.getsize(temp_blast_dir_list[input_x] + '/final_blast.opt') > 0:
                            # blast_line_count = 0
                            with open(temp_blast_dir + '/final_blast.out') as ipt:
                                for eachline_blast in ipt:
                                    eachline_blast = eachline_blast.strip('\n')
                                    col_blast = eachline_blast.strip().split('\t')

                                    ##only collect the match rate that is equal to 100
                                    if int(float(col_blast[2])) == 100:
                                        final_te_score = col_blast[-1]
                                        final_blast_line = col_blast[1] + '\t' + final_te_score
                                        store_final_blast_line_list.append(final_blast_line)
                        else:
                            print('this blast file is not existing')

                        ##if the store_final_blast_line_list is still []
                        if store_final_blast_line_list == []:
                            final_te_name = 'Unknown'

                        else:
                            ##select the largest one
                            store_te_name_count_dic = {}
                            for eachblast_line in store_final_blast_line_list:
                                col_blast = eachblast_line.split()

                                ##check the start of blast otuput
                                te_name = col_blast[0]
                                if 'rc_' in te_name:
                                    mt = re.match('rc_(.+)', te_name)
                                    true_te_nm = mt.group(1)
                                else:
                                    true_te_nm = te_name

                                if true_te_nm in store_te_name_count_dic:
                                    store_te_name_count_dic[true_te_nm] += 1
                                else:
                                    store_te_name_count_dic[true_te_nm] = 1

                            for eachte_nm in store_te_name_count_dic:
                                final_te_nm_match_line = loc_line + '\t' + eachte_nm + '\t' + str(
                                    store_te_name_count_dic[eachte_nm])
                                store_match_te_report_for_eachloc_line_list.append(final_te_nm_match_line)
                            z = [0]
                            while store_te_name_count_dic:
                                key, value = store_te_name_count_dic.popitem()
                                if value > z[0]:
                                    z = [value, [key]]
                                elif value == z[0]:
                                    z[1].append(key)
                            ##[120, ['a', 'b']]
                            ##it means the number of the key a and b are the same
                            ##we need to compare their average score
                            if len(z[1]) > 1:
                                store_te_nm_score_dic = {}  ##key is the target te name value is the sum of the score
                                for eachblast_line in store_final_blast_line_list:
                                    col_blast = eachblast_line.split()
                                    ##it allows to focus on the target te
                                    ##check the start of blast otuput
                                    te_name = col_blast[0]
                                    if 'rc_' in te_name:
                                        mt = re.match('rc_(.+)', te_name)
                                        true_te_nm = mt.group(1)
                                    else:
                                        true_te_nm = te_name
                                    if true_te_nm in z[1]:
                                        te_nm_score = float(col_blast[1])
                                        if true_te_nm in store_te_nm_score_dic:
                                            store_te_nm_score_dic[true_te_nm] += te_nm_score
                                        else:
                                            store_te_nm_score_dic[true_te_nm] = te_nm_score
                                ##select the max TE name and maybe they have same score so we will randomly select one
                                final_te_name = max(store_te_nm_score_dic, key=store_te_nm_score_dic.get)
                            else:
                                final_te_name = z[1][0]

                    else:
                        ##select the largest one
                        store_te_name_count_dic = {}
                        for eachblast_line in store_final_blast_line_list:
                            col_blast = eachblast_line.split()

                            ##check the start of blast otuput
                            te_name = col_blast[0]

                            if 'rc_' in te_name:
                                mt = re.match('rc_(.+)',te_name)
                                true_te_nm = mt.group(1)
                            else:
                                true_te_nm = te_name

                            if true_te_nm in store_te_name_count_dic:
                                store_te_name_count_dic[true_te_nm] += 1
                            else:
                                store_te_name_count_dic[true_te_nm] = 1

                        for eachte_nm in store_te_name_count_dic:
                            final_te_nm_match_line = loc_line + '\t' + eachte_nm + '\t' + str(store_te_name_count_dic[eachte_nm])
                            store_match_te_report_for_eachloc_line_list.append(final_te_nm_match_line)

                        z = [0]
                        while store_te_name_count_dic:
                            key, value = store_te_name_count_dic.popitem()
                            if value > z[0]:
                                z = [value, [key]]
                            elif value == z[0]:
                                z[1].append(key)
                        ##[120, ['a', 'b']]

                        ##it means the number of the key a and b are the same
                        ##we need to compare their average score
                        if len(z[1]) > 1:
                            store_te_nm_score_dic = {} ##key is the target te name value is the sum of the score
                            for eachblast_line in store_final_blast_line_list:
                                col_blast = eachblast_line.split()
                                ##it allows to focus on the target te

                                ##check the start of blast otuput
                                te_name = col_blast[0]

                                if 'rc_' in te_name:
                                    mt = re.match('rc_(.+)', te_name)
                                    true_te_nm = mt.group(1)
                                else:
                                    true_te_nm = te_name

                                if true_te_nm in z[1]:
                                    te_nm_score = float(col_blast[1])
                                    if true_te_nm in store_te_nm_score_dic:
                                        store_te_nm_score_dic[true_te_nm] += te_nm_score
                                    else:
                                        store_te_nm_score_dic[true_te_nm] = te_nm_score

                            ##select the max TE name and maybe they have same score so we will randomly select one
                            final_te_name = max(store_te_nm_score_dic, key=store_te_nm_score_dic.get)
                        else:
                            final_te_name = z[1][0]

                    other_line = ''
                    for i in range(4,len(col)):
                        other_line = other_line + '\t' + col[i]

                    final_line = chr + '\t' + start + '\t' + end + '\t' + loc_type + '\t' + final_te_name + other_line
                    #store_final_annot_te_line_list.append(final_line)
                    temp_te_genotype_annot_fl.write(final_line + '\n')


                else:
                    ##updation 122020
                    ##since the te_nm_str only has one name we will use this name
                    #if start == '28888041':
                    #    print ('28888041 can be generated output')
                    temp_te_genotype_annot_fl.write(eachline + '\n')

            else:
                temp_te_genotype_annot_fl.write(eachline + '\n')

    with open (input_output_dir + '/opt_genotype_te_match_count.txt','w+') as opt:
        for eachline in store_match_te_report_for_eachloc_line_list:
            opt.write(eachline + '\n')

    temp_te_genotype_annot_fl.close()







##updation 110520 add missing information default thr is 0.7 (missing_thr)
##updation 103120
##modify the annotation file
##1) filter out duplicated lines since the panTEs step already modifed the overlapped TEs, so no overlapped found in this step
##2) if location with the same family are close within 15 bp (default), we compare these locations and select one with the most number of 0/1 or 1/1
##this is for the c or s cases not for the reference
def modify_annotation_file (input_annotation_fl,input_working_dir,annot_gap_thr,sample_num_above_missing,annot_close_gap_thr):

    ##updation 110520
    ##sample_num_above_missing means if the samples with genotype above this sample number threshold, this location will be chosen to compare
    ##if all the insertions show the sample number below the missing threshold number, we will choose the largest one

    ##annot_gap_thr is 15 bp for default

    ##annotation fl is temp_te_genotype_annot.txt
    ##sort the annotation file
    cmd = 'sort -k1,1V -k2,2n ' + input_annotation_fl  + ' > ' + input_working_dir + '/temp_te_genotype_annot_sorted.txt'
    subprocess.call(cmd, shell=True)

    store_final_line_list = []
    store_c_s_type_line_list = []
    with open(input_working_dir + '/temp_te_genotype_annot_sorted.txt', 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            te_type = col[3]  ##c or s or o
            if te_type == 'c' or te_type == 's':
                store_c_s_type_line_list.append(eachline)
            else:
                store_final_line_list.append(eachline)

    store_all_line_list = []
    line_list = []
    final_line_count = len(store_c_s_type_line_list)
    store_all_line_dic = {} ##key is the location information and value is line from col[3] to end
    ##this dic will be used for the further line extraction
    line_count = 0
    for eachline in store_c_s_type_line_list:

        line_count += 1
        col = eachline.strip().split()
        current_chr = col[0]
        current_start = col[1]
        current_end = col[2]
        current_id = current_chr + '_' + current_start + '_' + current_end
        current_te_nm = col[4]

        other_line_str = ''
        for i in range(3,len(col)):
            other_line_str = other_line_str + '\t' + col[i]
        store_all_line_dic[current_id] = other_line_str

        ##updation 110520 add the sample information for the item from 5 to len(col)
        sample_num = len(col) - 5

        store_line_dic = {'chr': current_chr, 'st': current_start, 'ed': current_end, 'te':current_te_nm, 'otherline': other_line_str,'sample_num':sample_num}

        if line_count == 1:
            line_list.append(store_line_dic)
            if line_count == final_line_count:
                store_all_line_list.append(line_list)
        else:
            pre_line_dic = line_list[-1]
            pre_chr = pre_line_dic['chr']
            #pre_start = pre_line_dic['st']
            pre_ed = pre_line_dic['ed']
            pre_te = pre_line_dic['te']

            if pre_chr == current_chr:

                ##updating 122220
                ##add a case that two TEs are the same location or have less than 5 bp gap
                ##if the TEs are not the same.
                ##if the TEs are the same, it will consider in the annot_gap_thr

                if (int(current_start) - int(pre_ed)) <= int(annot_gap_thr):
                    if current_te_nm == pre_te:
                        line_list.append(store_line_dic)
                        if line_count == final_line_count:
                            store_all_line_list.append(line_list)

                    else:
                        if (int(current_start) - int(pre_ed)) <= int(annot_close_gap_thr):
                            line_list.append(store_line_dic)
                            if line_count == final_line_count:
                                store_all_line_list.append(line_list)
                        else:
                            store_all_line_list.append(line_list)
                            line_list = []
                            line_list.append(store_line_dic)
                            if line_count == final_line_count:
                                store_all_line_list.append(line_list)

                else:
                    store_all_line_list.append(line_list)
                    line_list = []
                    line_list.append(store_line_dic)
                    if line_count == final_line_count:
                        store_all_line_list.append(line_list)

                #if (int(current_start) - int(pre_ed)) <= int(annot_gap_thr) and current_te_nm == pre_te:
                #    line_list.append(store_line_dic)
                #    if line_count == final_line_count:
                #        store_all_line_list.append(line_list)

                #else:
                #    store_all_line_list.append(line_list)
                #    line_list = []
                #    line_list.append(store_line_dic)
                #    if line_count == final_line_count:
                #        store_all_line_list.append(line_list)
            else:
                store_all_line_list.append(line_list)
                line_list = []
                line_list.append(store_line_dic)
                if line_count == final_line_count:
                    store_all_line_list.append(line_list)

    ##compare the list from the store_all_line_list to select one TE
    for eachlist in store_all_line_list:
        if len(eachlist) == 1:
            line_dic = eachlist[0]
            final_line = line_dic['chr'] + '\t' + line_dic['st'] + '\t' + line_dic['ed'] + line_dic['otherline']
            store_final_line_list.append(final_line)
        else:

            dic_id = 0
            store_dic_dic = {}
            store_compared_dic = {}
            for eachdic in eachlist:
                dic_id += 1
                store_dic_dic[str(dic_id)] = eachdic

                other_line = eachdic['otherline']
                other_line_col = other_line.split()

                ##updation 110520
                sample_num = eachdic['sample_num']

                genos_score = 0
                for i in range(2,len(other_line_col)):
                    sp_str_col = other_line_col[i].split(';')
                    genos = sp_str_col[-1]

                    if genos == '0/1':
                        genos_score += 1
                    if genos == '1/1':
                        genos_score += 2

                ##updation 110520
                if sample_num >= int(sample_num_above_missing):
                    store_compared_dic[str(dic_id)] = genos_score

            ##updation 110520
            ##if all locations in the list do not meet the threshold
            ##we need to reanalyze this list and choose the largest one
            if store_compared_dic == {}:
                for eachdic in eachlist:
                    dic_id += 1
                    store_dic_dic[str(dic_id)] = eachdic
                    other_line = eachdic['otherline']
                    other_line_col = other_line.split()
                    genos_score = 0
                    for i in range(2, len(other_line_col)):
                        sp_str_col = other_line_col[i].split(';')
                        genos = sp_str_col[-1]
                        if genos == '0/1':
                            genos_score += 1
                        if genos == '1/1':
                            genos_score += 2
                    store_compared_dic[str(dic_id)] = genos_score


            ##select the max id
            max_dic_id = max(store_compared_dic, key=store_compared_dic.get)
            max_dic = store_dic_dic[max_dic_id]

            final_line = max_dic['chr'] + '\t' + max_dic['st'] + '\t' + max_dic['ed'] + max_dic['otherline']
            store_final_line_list.append(final_line)

    return (store_final_line_list)













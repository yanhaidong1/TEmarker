#!/usr/bin/env python

##updation 110520 add the MAF filtration
##updation 110120 add an opt to generate vcf file

##import modules
from Bio import SeqIO
import re
import sys
import subprocess

def generate_sample_list (input_all_sample_list_file):
    all_name_list = []
    store_change_sp_nm_dic = {}
    with open (input_all_sample_list_file,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            all_name_list.append(col[0])
            store_change_sp_nm_dic[col[0]] = col[1]
    return (all_name_list,store_change_sp_nm_dic)

##def a function to store the name of chr and its relative code number
def store_chr (input_genome_file):
    ##initiate a dic to store the chr and its relative number
    chr_num_dic = {}
    ##chr will be the key and number will be the value
    count = 0
    for seq_record in SeqIO.parse(input_genome_file,'fasta'):
        count += 1
        id = seq_record.id
        chr_num_dic[id] = str(count)

    return (chr_num_dic)

def generate_vcf (input_opt_classify,chr_num_dic,all_name_list,store_change_sp_nm_dic):
    ##initiate a final dic to store information
    store_final_line_list = []
    #final_line_dic = {}

    ##generate sample name string
    name_string = ''
    for eachnm in all_name_list:
        name_string = name_string + '\t' + store_change_sp_nm_dic[eachnm]

    ##collect the sa_nm_string to generate the title information
    name_list_dic = {}
    ##initiate a title string
    basic_title_string = '#CHROM' + '\t' + 'POS' + '\t' + 'ID' + '\t' + 'REF' + '\t' + 'ALT' + '\t' + 'QUAL' + '\t' + \
                         'FILTER' + '\t' + 'INFO' + '\t' + 'FORMAT' + name_string
    store_final_line_list.append(basic_title_string)
    ##updation 8.19
    #calculate_fam_dic = {}

    with open (input_opt_classify,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()

            ##get the chr information
            ##the vcf file must need chr ordered by position and the chr name should be number information
            chr_num_id = chr_num_dic[col[0]]

            ##get the position information
            ##the position information will extract the first nucleotide for the opt_classify file
            posi = col[1]

            ##get the ID information
            ##vcf id is the .
            ##updation 12.03
            ID = chr_num_id + '_' + posi

            ##get the REF
            ##test using A and T and later will be revised to another symbol to indicate
            ##set A as the consensus reference sequence
            REF = 'R'

            ##set B as the ALT
            ALT = 'A'

            ##get the QUAL
            ##there is no information for that, so we use 100 to indicate the QUAL
            QUAL = '.'

            ##get the FILTER
            ##there is no information for that, so we use . to indicate the FILTER
            FILTER = '.'

            ##do not generate the INFO information
            tetype = col[3] #c s or o
            te_nm = col[4]
            if tetype == 'o':
                infor_str = 'REF:' + te_nm + ';ALT:na'
            else:
                infor_str = 'REF:na;ALT:' + te_nm

            #INFO = 'na'

            ##get the FORMAT
            FORMAT = 'GT'

            basic_infor_line = chr_num_id + '\t' + posi + '\t' + ID + '\t' + REF + '\t' + ALT + '\t' + QUAL + '\t' \
                               + FILTER + '\t' + infor_str + '\t' + FORMAT

            ##the annotated information starts from col[7]
            total_col_num = len(col)

            ##store the genotype into a name dic
            name_genotype_dic = {}  ##key is the sample name and genotype is the value
            for i in range(5,total_col_num):

                ##get the sample name information
                sample_col = col[i].split(':')
                sample_nm = sample_col[0]

                ##some name contains the '-' that will be removed
                ##updation8.10 delete this filteration
                #if '_' not in sample_nm:
                #print(sample_nm)
                name_list_dic[sample_nm] = 1
                ##some information is not right and it is about the whole name of the sample
                ##should check it out
                ##get the sample information
                ##For the GT
                ############################################################
                ##becareful: change the 0/0 to 1/1 and change the 1/1 to 0/0
                #GT = ''
                ##only use one time
                #if genotype == '1/1':
                #    GT = '0/0'
                #if genotype == '0/0':
                #    GT = '1/1'
                #if genotype == '0/1':
                #    GT = '0/1'
                ##############################################################
                #GT = genotype
                #SA = GT

                ##extract the genotype information for each sample

                annot_col = col[i].split(';')
                genotype = annot_col[2]
                name_genotype_dic[sample_nm] = genotype

            ##for each name in the all_sample_list
            ##if there is no information in the classify_geno file, we will add genotype 0/0 for this sample
            ##generate a string to contain sample information
            samples_string = ''

            for eachnm in all_name_list:
                if eachnm in name_genotype_dic.keys():
                    samples_string = samples_string + '\t' + name_genotype_dic[eachnm]
                else:
                    ##updation 8.16 change the missing data to ./.
                    samples_string = samples_string + '\t' + './.'

            final_line = basic_infor_line + samples_string

            ############
            ##filter the samples_string that only allows to number samples
            number_sample = len(all_name_list)

            ##filter all the samples which are the same genotypes
            final_line_col_list = final_line.split('\t')

            geno_num_dic = {}
            for eachcol in final_line_col_list:
                if eachcol in geno_num_dic:
                    geno_num_dic[eachcol] += 1
                else:
                    geno_num_dic[eachcol] = 1

            ##check the number of geno_num_dic
            num_homo = 0
            num_heter = 0
            num_ful_heter = 0
            ##updation
            num_missing = 0
            if '0/0' in geno_num_dic.keys():
                num_homo = geno_num_dic['0/0']
            if '0/1' in geno_num_dic.keys():
                num_heter = geno_num_dic['0/1']
            if '1/1' in geno_num_dic.keys():
                num_ful_heter = geno_num_dic['1/1']
            if './.' in geno_num_dic.keys():
                num_missing = geno_num_dic['./.']

            if num_homo != number_sample and \
                    num_heter != number_sample and \
                    num_ful_heter != number_sample and \
                    num_missing != number_sample:
                store_final_line_list.append(final_line)

    return (store_final_line_list)

def remove_duplicate (store_vcf_line_list):

    store_final_line_list = []
    store_replicate_name = {}
    row_name_dic = {}
    dup_count = 0
    for eachline in store_vcf_line_list:
        eachline = eachline.strip('\n')
        col = eachline.strip().split()

        row_nm = col[0] + '_' + col[1]

        if row_nm in row_name_dic.keys():
            row_name_dic[row_nm] += 1
        else:
            row_name_dic[row_nm] = 1

    for eachrow in row_name_dic:
        if row_name_dic[eachrow] > 1:
            store_replicate_name[eachrow] = 1

    for eachline in store_vcf_line_list:
        col = eachline.strip().split()
        rownm = col[0] + '_' + col[1]

        if rownm not in store_replicate_name:
            store_final_line_list.append(eachline)
        else:
            dup_count += 1

    print('After the duplication filtration, ' + str(dup_count) + ' TEs are removed')

    return (store_final_line_list)

def filter_out_missing (genotype_file,miss_thr):

    ##get the genotype file index: vcf, hmp or genos
    mt = re.match('.+/(.+)',genotype_file)
    fl_nm = mt.group(1)

    fl_index = ''
    if '.vcf' in fl_nm:
        fl_index = 'vcf'
    if 'genos.txt' in fl_nm:
        fl_index = 'genos'
    if 'hmp.txt' in fl_nm:
        fl_index = 'hmp'

    ##generate_final_line_list
    final_line_list = []
    ##filter out missing loc
    missing_count = 0
    with open (genotype_file,'r') as ipt:
        for eachline in ipt:
            if not eachline.startswith('#'):
                eachline = eachline.strip('\n')
                col_list = eachline.strip().split()

                if fl_index == 'vcf':
                    ##calculate sample number
                    sample_num = len(col_list[9:])
                    thr_miss_sample_num = int(float(sample_num)*(1-float(miss_thr)))
                    missing_sample_num = 0
                    for eachitem in col_list:
                        if './.' == eachitem:
                            missing_sample_num += 1
                    if missing_sample_num <= thr_miss_sample_num:
                        final_line_list.append(eachline)
                    else:
                        missing_count += 1

                if fl_index == 'hmp':
                    sample_num = len(col_list[11:])
                    thr_miss_sample_num = int(float(sample_num) * (1 - float(miss_thr)))
                    missing_sample_num = 0
                    for eachitem in col_list:
                        if 'NN' == eachitem:
                            missing_sample_num += 1
                    if missing_sample_num <= thr_miss_sample_num:
                        final_line_list.append(eachline)
                    else:
                        missing_count += 1

                if fl_index == 'genos':
                    sample_num = len(col_list[4:])
                    thr_miss_sample_num = int(float(sample_num) * (1 - float(miss_thr)))
                    missing_sample_num = 0
                    for eachitem in col_list:
                        if 'NA' == eachitem:
                            missing_sample_num += 1
                    if missing_sample_num <= thr_miss_sample_num:
                        final_line_list.append(eachline)
                    else:
                        missing_count += 1
            else:
                eachline = eachline.strip('\n')
                final_line_list.append(eachline)

    print('After the missing location filtration, ' + str(missing_count) + ' TEs are removed')

    return (final_line_list)

def MAF_filtration (opt_flt_missing_vcf_fl,thr_maf):

    store_final_line_list = []
    loc_count = 0
    low_allfreq_count = 0
    with open(opt_flt_missing_vcf_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            loc_count += 1
            if loc_count != 1:

                col = eachline.strip().split()

                ##two cases but the two cases has the same calculation for the maf
                ##they both calculate the 0/1 and 1/1 divided by the 0/0
                ##0/1 and 1/1 suggest a insertion  0/0 is the reference
                ##calculate all the allele counts but not the ./.
                all_all_c = 0 ##need to double
                heter_c = 0 ##0/1
                fl_heter_c = 0 ##1/1

                for i in range(9, len(col)):
                    if col[i] != './.':
                        all_all_c += 1
                        if col[i] == '0/1':
                            heter_c += 1
                        if col[i] == '1/1':
                            fl_heter_c += 1

                final_all_c = 2*all_all_c
                af = (2*fl_heter_c + heter_c)/final_all_c
                if af >= 0.5:
                    maf = 1 - af
                else:
                    maf = af

                if float(maf) >= float(thr_maf):
                    store_final_line_list.append(eachline)
                else:
                    low_allfreq_count += 1

            else:
                store_final_line_list.append(eachline)

    print('After the MAF filtration, ' + str(low_allfreq_count) + ' TEs are removed')

    return (store_final_line_list)

def unknown_filtration (vcf_fl):

    store_final_line_list = []
    loc_count = 0
    with open(vcf_fl, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            loc_count += 1
            if loc_count != 1:
                col = eachline.strip().split()
                if 'Unknown' not in col[7]:
                    store_final_line_list.append(eachline)

            else:
                store_final_line_list.append(eachline)

    return (store_final_line_list)

def transfer_vcf_to_genos (vcf_fl):

    store_final_line_list = []
    count = 0
    with open (vcf_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()

            count += 1
            if count != 1:
                ##get the chr information

                ##generate rs
                rs = 'loc_' + col[2]
                rs_col_nm = 'loc_' + col[2]

                ##generate chrom
                chrom = col[0]
                ##generate pos
                pos = col[1]

                basic_infor_line = rs_col_nm + '\t' + rs  + '\t' + chrom + '\t' + pos

                ##the annotated information starts from col[7]
                total_col_num = len(col)

                ##store the genotype into a name dic
                name_genotype_dic = {}  ##key is the sample name and genotype is the value
                GT_str = ''
                for i in range(9,total_col_num):
                    genotype = col[i]

                    GT = ''
                    if '0/1' in genotype:
                        GT = str(int(0))
                    if '1/1' in genotype:
                        GT = str(int(-1))
                    if '0/0' in genotype:
                        GT = str(int(1))
                    if './.' in genotype:
                        GT = 'NA'

                    GT_str = GT_str + '\t' + GT

                final_line = basic_infor_line + GT_str
                store_final_line_list.append(final_line)

            else:
                basic_first_line = 'rs_col_nm' + '\t' + 'rs' + '\t' + 'chr' + '\t' + 'bp'
                frist_line_str = ''
                for i in range(9,len(col)):
                    frist_line_str = frist_line_str + '\t' + col[i]
                final_first_line = basic_first_line + frist_line_str
                store_final_line_list.append(final_first_line)

    return (store_final_line_list)


def change_to_bialleles (opt_fltmissing_fltmaf_fl):

    store_final_line_list = []
    count = 0
    with open (opt_fltmissing_fltmaf_fl,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            count += 1
            if count != 1:
                col = eachline.strip().split()

                new_line = col[0] + '\t' + col[1] + '\t' + col[2] + '\t' + col[3] + '\t' + col[4] + '\t' + col[5] + '\t' \
                           + col[6] + '\t' + col[7] + '\t' + col[8]
                for eachitem in col[9:]:
                    new_item = ''
                    if eachitem == '0/1':
                        new_item = '1|1'
                    if eachitem == '0/0':
                        new_item = '0|0'
                    if eachitem == '1/1':
                        new_item = '1|1'
                    if eachitem == './.':
                        new_item = './.'
                    #if eachitem != '0/1' and eachitem != '1/1' and eachitem != '0/0':
                    #    print(eachitem)
                    new_line = new_line + '\t' + new_item

                store_final_line_list.append(new_line)

            else:
                store_final_line_list.append(eachline)

    return (store_final_line_list)


def modify_wrong_family (te_fam_nm):
    if '_D._' in te_fam_nm:
        new_te_fam_nm = te_fam_nm.replace("_D._", "_DNA_")
    else:
        new_te_fam_nm = te_fam_nm
    return (new_te_fam_nm)

def cal_fam_from_vcf (vcf_file):

    store_final_line_list = []
    store_insertion_dic = {}
    store_deletion_dic = {}
    store_all_dic = {}
    with open (vcf_file,'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            if not eachline.startswith('#'):
                col = eachline.strip().split()
                TE_fam_string = col[7]
                if 'REF:na' in TE_fam_string:
                    if re.match('REF:na;ALT:(.+)_TE\d+$',TE_fam_string):
                        mt = re.match('REF:na;ALT:(.+)_TE\d+$',TE_fam_string)
                        te_fam_nm = modify_wrong_family(mt.group(1))
                    else:
                        mt = re.match('REF:na;ALT:(.+)', TE_fam_string)
                        te_fam_nm = modify_wrong_family(mt.group(1))

                    if te_fam_nm in store_insertion_dic:
                        store_insertion_dic[te_fam_nm] += 1
                    else:
                        store_insertion_dic[te_fam_nm] = 1

                    if te_fam_nm in store_all_dic:
                        store_all_dic[te_fam_nm] += 1
                    else:
                        store_all_dic[te_fam_nm] = 1

                else:
                    mt = re.match('REF:(.+)_TE\d+;ALT:na', TE_fam_string)
                    te_fam_nm = modify_wrong_family(mt.group(1))

                    if te_fam_nm in store_deletion_dic:
                        store_deletion_dic[te_fam_nm] += 1
                    else:
                        store_deletion_dic[te_fam_nm] = 1

                    if te_fam_nm in store_all_dic:
                        store_all_dic[te_fam_nm] += 1
                    else:
                        store_all_dic[te_fam_nm] = 1

    ##store final line
    for eachtenm in store_insertion_dic:
        final_line = 'Insertion' + '\t' + eachtenm + '\t' + str(store_insertion_dic[eachtenm])
        store_final_line_list.append(final_line)

    for eachtenm in store_deletion_dic:
        final_line = 'Deletion' + '\t' + eachtenm + '\t' + str(store_deletion_dic[eachtenm])
        store_final_line_list.append(final_line)

    for eachtenm in store_all_dic:
        final_line = 'All' + '\t' + eachtenm + '\t' + str(store_all_dic[eachtenm])
        store_final_line_list.append(final_line)

    return (store_final_line_list)

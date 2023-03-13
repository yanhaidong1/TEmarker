#!/usr/bin/env python

##updating 030823 we do not summarize this file as the name user provide are different
##updating 041121 add an option to filter out read coverage
##udpation 122420 add an original vcf file and change output file name
##updation 110520 add maf filtration argument
##updation 110320 analyze on the new output file
##updation 110120 for a new genotyping file
##updation 101720 add the split reads information and full te location information on the opt_te_loc_fam.txt that contains all the information without filtration
##updation 042720 remove the opt_calculate_loc_num.txt
##updation clean the output and change name of output
##updation 01.22.20 add fl miss for the output with filtering no pos
##updation 12.20 add filter out unknown position in the vcf
##updation 12.16 add family information in the vcf and calculate the number of family
##updation 12.03 add a output showing the TE fam for each locus
##updation 11.11 add thr to filter out missing data and generate output
##updation 10.24 add material number for the bi vcf output
##updation 10.12 add the no fam version
##this script pipeline is to generate the output format from the TEScape_marker_genos

##BUILT-IN MODULES
import argparse
import sys
import os
import subprocess

##updation 110320
from TEmarker_utils import TM_opt_vcf as opt_vcf
#from TEmarker_utils import TM_opt_vcf_no_fam as opt_vcf_no_fam
#from TEmarker_utils import TM_opt_vcf_sort as opt_vcf_sort
#from TEmarker_utils import TM_opt_vcf_change_nm as opt_vcf_ch_nm
#from TEmarker_utils import TM_opt_vcf_bi_allele as opt_vcf_bi_al
#from TEmarker_utils import TM_opt_hmp as opt_hmp
#from TEmarker_utils import TM_opt_hmp_sort as opt_hmp_sort
#from TEmarker_utils import TM_opt_hmp_change_nm as opt_hmp_ch_nm
#from TEmarker_utils import TM_opt_genos as opt_genos
#from TEmarker_utils import TM_opt_genos_sort as opt_genos_sort
#from TEmarker_utils import TM_opt_genos_change_nm as opt_genos_ch_nm
#from TEmarker_utils import TM_opt_filter_out_miss as opt_fl_out_miss
##updation 1203
#from TEmarker_utils import TM_opt_loc_te_fam as opt_te_loc_fam
##updaiton 1216
#from TEmarker_utils import TM_opt_add_fam_cal_fam_num as opt_add_cal_fam_num
##updation 1220
#from TEmarker_utils import TM_opt_filter_out_unknown_pos as opt_fl_out_un_pos




##SCRIPTS
def get_parsed_args():
    parser = argparse.ArgumentParser(description="Generate final output")

    ##############################
    ##parse the required arguments
    parser.add_argument("-o", dest='output_dir', default="./", help="Output directory"
                                                                    "Default: ./ ")

    parser.add_argument("-d", dest='working_dir', default="./", help="Working directory"
                                                                    "Default: ./ ")

    parser.add_argument('-genos_f', dest='genotype_file', help="Users provide the genotype file from TEmarker_genotyping")

    parser.add_argument('-m_f', dest='material_file', help="Users provide the materials information."
                                                           "Two columns including orignial sample name (left) and new sample "
                                                           "name (right)."
                                                           "This file also helps to create opt file that requires material number")

    parser.add_argument('-fas_f', dest='genome_file', help="Users provide the genome file.")

    ####################
    ##optional arguments
    parser.add_argument('-miss_thr', dest='thr_miss', help="Users provide a threshold to filter out missing sample in the output."
                                                           "For example, if thr_miss is equal to 0.7, and the proportion of the missing "
                                                           "samples is over 0.3, and this locus will be filtered out"
                                                           "Default: 0.7")

    parser.add_argument('-maf_thr', dest='thr_maf', help="Users provide maf filtration threshold."
                                                         "Dfault: 0.05")


    parser.add_argument('-read_sup', dest='read_sup', help="Users provide the minimum reads that support the TE deletion or insertion."
                                                           "For example, if there are 3 reads support insertion of TEs, we will consider this location."
                                                           "Otherwise, this location will be considered as NULL or ./."
                                                           "Dfault: 3")

    ##parse of parameters
    args = parser.parse_args()
    return args


def main(argv=None):
    if argv is None:
        argv = sys.argv
    args = get_parsed_args()

    ##check the input file
    genotype_file = args.genotype_file
    if args.genotype_file is None:
        print('Cannot find input genotype file, please provide that')
        return  ##import to close the script if not find the te lib
    else:
        try:
            file = open(args.genotype_file, 'r')  ##check if the file is not the right file
        except IOError:
            print('There was an error opening the file!')
            return

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

    genome_file = args.genome_file
    if args.material_file is None:
        print('Cannot find input genome file, please provide that')
        return  ##import to close the script if not find the te lib
    else:
        try:
            file = open(args.genome_file, 'r')  ##check if the file is not the right file
        except IOError:
            print('There was an error opening the file!')
            return

    ##set working_dir and output_dir
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

    ##run the processes
    ##updation 11.11 create dir to store the output format seperately
    ##for vcf

    ##set parameter
    if args.thr_miss is not None:
        thr_miss = args.thr_miss
    else:
        thr_miss = 0.7

    if args.thr_maf is not None:
        thr_maf = args.thr_maf
    else:
        thr_maf = 0.05

    if args.read_sup is not None:
        read_sup = args.read_sup
    else:
        read_sup = 3


    ##########################
    ##generate vcf output file
    chr_num_dic = opt_vcf.store_chr(genome_file)

    ##updating 031323
    ##write a plot to show the chr ID and IDs in the vcf file
    with open (working_dir + '/opt_original_chrID_and_IDs_in_output_VCF.txt','w+') as opt:
        for eachchrID in chr_num_dic:
            outputID = str(chr_num_dic[eachchrID])
            final_line = eachchrID + '\t' + outputID
            opt.write(final_line + '\n')


    all_name_list, store_change_sp_nm_dic = opt_vcf.generate_sample_list(material_file)
    store_final_vcf_line_list,store_final_list_add_ratio_noflt_line = opt_vcf.generate_vcf(genotype_file, chr_num_dic, all_name_list, store_change_sp_nm_dic,read_sup)
    store_final_vcf_uniq_line_list = opt_vcf.remove_duplicate(store_final_vcf_line_list)
    with open(working_dir + '/opt.vcf', 'w+') as opt:
        for eachline in store_final_vcf_uniq_line_list:
            opt.write(eachline + '\n')

    ##updating 041121
    ##generate a raw vcf that contains all ratio and ratio value information
    with open(working_dir + '/temp_raw_noflt_add_ratio.vcf', 'w+') as opt:
        for eachline in store_final_list_add_ratio_noflt_line:
            opt.write(eachline + '\n')


    ##sorted vcf
    #cmd = 'cat ' + working_dir + '/opt.vcf' + ' | (sed -u 1q; sort -k1,1V -k2,2n) > ' + working_dir + '/opt_sorted.vcf'
    #subprocess.call(cmd, shell=True)
    first_line = ''
    store_other_line_list = []
    with open (working_dir + '/opt.vcf','r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            if eachline.startswith('#'):
                first_line = eachline
            else:
                store_other_line_list.append(eachline)

    with open (working_dir + '/temp_unsorted_noheader.vcf','w+') as opt:
        for eachline in store_other_line_list:
            opt.write(eachline + '\n')

    with open (working_dir + '/temp_firstline.vcf','w+') as opt:
        opt.write(first_line + '\n')

    cmd = 'sort -k1,1V -k2,2n ' + working_dir + '/temp_unsorted_noheader.vcf > ' + working_dir + '/temp_unsorted_noheader_sorted.vcf'
    subprocess.call(cmd, shell=True)

    cmd = 'cat ' + working_dir + '/temp_firstline.vcf ' +  working_dir + '/temp_unsorted_noheader_sorted.vcf > ' + working_dir + '/opt_sorted.vcf'
    subprocess.call(cmd, shell=True)

    ##updation 122420 add an original vcf file
    cmd = 'cp ' + working_dir + '/opt_sorted.vcf' + ' ' + output_dir + '/opt.vcf'
    subprocess.call(cmd, shell=True)

    ##updation 122420 add the original genos file
    store_final_genos_line_list = opt_vcf.transfer_vcf_to_genos(output_dir + '/opt.vcf')
    with open(output_dir + '/opt_genos.txt', 'w+') as opt:
        for eachline in store_final_genos_line_list:
            opt.write(eachline + '\n')


    ##updation 122420 add the opt_bi.vcf
    final_line_bi_list = opt_vcf.change_to_bialleles(output_dir + '/opt.vcf')
    with open (output_dir + '/opt_bi.vcf','w+') as opt:
        for eachline in final_line_bi_list:
            opt.write(eachline + '\n')


    ##filter out missing
    final_line_list = opt_vcf.filter_out_missing(working_dir + '/opt_sorted.vcf', thr_miss)
    with open(working_dir + '/opt_fltmissing.vcf', 'w+') as opt:
        for eachline in final_line_list:
            opt.write(eachline + '\n')

    ##filter out MAF < 0.05 (default)
    ##need a script to have filtration
    final_line_maf_list = opt_vcf.MAF_filtration (working_dir + '/opt_fltmissing.vcf',thr_maf)
    with open(working_dir + '/opt_fltmissing_fltmaf_notfltunknown.vcf', 'w+') as opt:
        for eachline in final_line_maf_list:
            opt.write(eachline + '\n')

    ##filter out unknown
    final_line_without_unknown_list = opt_vcf.unknown_filtration(working_dir + '/opt_fltmissing_fltmaf_notfltunknown.vcf')
    with open(output_dir + '/opt_fltmissing_fltmaf.vcf', 'w+') as opt:
        for eachline in final_line_without_unknown_list:
            opt.write(eachline + '\n')

    ##generate genos output
    store_final_genos_fltmissing_fltmaf_line_list = opt_vcf.transfer_vcf_to_genos(output_dir + '/opt_fltmissing_fltmaf.vcf')
    with open(output_dir + '/opt_genos_fltmissing_fltmaf.txt', 'w+') as opt:
        for eachline in store_final_genos_fltmissing_fltmaf_line_list:
            opt.write(eachline + '\n')

    ##generate vcf file that only contains 0/0 and 1/1 that could be used for the PopLDdecay
    ##change 0/0 and 1/1 to 0|0 and 1|1
    #final_line_bi_list = opt_vcf.change_to_bialleles(output_dir + '/opt_fltmissing_fltmaf_fltukn.vcf')
    #with open (output_dir + '/opt_fltmissing_fltmaf_fltukn_bialleles.vcf','w+') as opt:
    #    for eachline in final_line_bi_list:
    #        opt.write(eachline + '\n')

    ##generate summary information
    ##updating 030823 we do not summarize this file as the name user provide are different
    #store_summary_line_list = opt_vcf.cal_fam_from_vcf(output_dir + '/opt_fltmissing_fltmaf.vcf')
    #with open (output_dir + '/opt_summary_te_family.txt','w+') as opt:
    #    for eachline in store_summary_line_list:
    #        opt.write(eachline + '\n')

    ##updating 041121 add an output for the final ratio vcf
    final_line_list = opt_vcf.keep_same_to_final(output_dir + '/opt_fltmissing_fltmaf.vcf',working_dir + '/temp_raw_noflt_add_ratio.vcf')
    with open(output_dir + '/opt_fltmissing_fltmaf_addratio.vcf', 'w+') as opt:
        for eachline in final_line_list:
            opt.write(eachline + '\n')


if __name__ == "__main__":
    main()


















#!/usr/bin/env python

##updation 12.03 add the opt_temp_sort_consensus_TE.csv to the D03 working_dir
##updation 11.28 revise back to the single process to try if the new TESM_conse_construct work
##this script pipeline is to construct the consensus TE location from the bed dir

##BUILT-IN MODULES
import argparse
import sys
import os

from TEmarker_utils import TM_panTEs_construct as construct_consensus
from TEmarker_utils import TM_panTEs_modify_consensus as modify_consensus
from TEmarker_utils import TM_panTEs_filter as filter_consensus

##SCRIPTS
def get_parsed_args():
    parser = argparse.ArgumentParser(description="Generate pan TE insertions")

    ##############################
    ##parse the required arguments
    parser.add_argument("-d", dest='working_dir', default="./", help="Working directory"
                                                                     "Default: ./ ")

    parser.add_argument("-o", dest='output_dir', default="./", help="Output directory"
                                                                    "Default: ./ ")

    parser.add_argument('-b_d', dest='bed_dir',help="Users provide the bed dir that will be used to generate the pan TE insertions")


    ##optional
    parser.add_argument('-c_bp',dest='combine_bp',help="Users provide a range that connect insertion. For example, in a location,"
                                                      "sample 1 has a TE insertion with location of 300 in chr01, and sample 2 has a TE insertion with location of 305 in chr01."
                                                      "If cbp is set as 5, these two insertions will be combined."
                                                      "If users do not want to combine, they can set this value as 0."
                                                      "Default: 5.")

    parser.add_argument('-f_th',dest='filter_thr',help="Users provide a threshold (decimals; from 0 to 1) to filter out location based on sample proportion."
                                                                   "For example, fncth is set to 0.05, if sample proportion is 0.04, this location will be filtered out."
                                                                   "Default: 0.05."
                                                                   "Sample proportion: if a TE is inserted in a location, we calcualte the number (A) of samples with this insertion,"
                                                                   " and the number of samples without this insertion is B. If A is smaller than B, sample proportion = A/(A+B),"
                                                                   " otherwise, sample proportion = B/(A+B).")

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

    if args.combine_bp is not None:
        cbp_value = args.combine_bp
    else:
        cbp_value = 5

    ##updation 050220
    ##we have update the na to specific value
    ##so na is in the combine situation
    #if args.filter_combine_thr is not None:
    #    fcth_value = args.filter_combine_thr
    #else:
    #    fcth_value = 5

    if args.filter_thr is not None:
        fth_value = args.filter_thr
    else:
        fth_value = 0.05

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

    ##run the scripts
    ########
    ##step 1: construct consensus TE insertions
    ##generate working_dir for the first step
    D01_construct_consensus_dir = working_dir + '/D01_construct_consensus_dir'
    if not os.path.exists(D01_construct_consensus_dir):
        os.makedirs(D01_construct_consensus_dir)

    store_final_line_list,sample_number = construct_consensus.construct_consensus_loc(bed_dir)

    ##write out the final output
    with open(D01_construct_consensus_dir + '/opt_filter_TE_loc.txt', 'w+') as opt:
        for eachline in store_final_line_list:
            opt.write(eachline + '\n')

    ########
    ##step 2: modify the consensus TE insertions (combine the close insertions)
    D02_modify_consensus_dir = working_dir + '/D02_modify_consensus_dir'
    if not os.path.exists(D02_modify_consensus_dir):
        os.makedirs(D02_modify_consensus_dir)

    ##updation add D02_modify_consensus_dir as argument
    line_dic, store_id_list, final_line_list = modify_consensus.modification_consensus(D01_construct_consensus_dir + '/opt_filter_TE_loc.txt', cbp_value,D02_modify_consensus_dir)
    final_line_list = modify_consensus.combine_close_TE(line_dic, store_id_list, final_line_list)

    ##updation 050220
    new_final_line_list = modify_consensus.calcualte_pro_to_replace_na_in_combine_loc (final_line_list,sample_number)
    with open(D02_modify_consensus_dir + '/opt_modify_filter_TE_loc.txt', 'w+') as opt:
        for eachline in new_final_line_list:
            opt.write(eachline + '\n')

    ########
    ##step 3: filter out insertions in the comdified TE insertion file
    #D03_filter_consensus_dir = working_dir + '/D03_filter_consensus_dir'
    #if not os.path.exists(D03_filter_consensus_dir):
    #    os.makedirs(D03_filter_consensus_dir)

    store_line_list = filter_consensus.check_filter_opt(D02_modify_consensus_dir + '/opt_modify_filter_TE_loc.txt', fth_value)
    with open(output_dir + '/opt_pan_TEs.txt', 'w+') as opt:
        for eachline in store_line_list:
            opt.write(eachline + '\n')


if __name__ == "__main__":
    main()
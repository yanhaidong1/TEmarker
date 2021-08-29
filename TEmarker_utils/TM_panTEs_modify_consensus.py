#!/usr/bin/env python
##updating 032821 change the csv to the txt for opt_temp_sort_consensus_TE

##import modules
import pandas as pd
import re

##define a function to calcualte the total line number
def calcualte_all_line_number (file):
    line_count = 0
    with open (file,'r') as ipt:
        for eachline in ipt:
            line_count += 1
    return (line_count)

def modification_consensus (input_consensus_TE_file,thresd,D02_modify_consensus_dir):
    final_line_list = []
    ##sort the TE file according
    ltr_fl = pd.read_csv(input_consensus_TE_file, header=None, delimiter=r"\s+")
    ##sort chr and start region
    ##updation8.11 add the sample information
    ltr_fl.columns = ['chr', 'start', 'end','te_nm','direct','num1', 'num2','ratio','sample_string']
    ltr_sort_fl = ltr_fl.sort_values(by=['chr', 'start'])
    ltr_sort_fl.to_csv(D02_modify_consensus_dir + '/opt_temp_sort_consensus_TE.txt', sep='\t')

    ##get the total line number of the file
    total_line_num = calcualte_all_line_number (input_consensus_TE_file)

    temp_id_list = []
    line_dic = {}
    store_id_list = []

    id = 0
    line_count = 0
    with open (D02_modify_consensus_dir + '/opt_temp_sort_consensus_TE.txt','r') as ipt:
        for eachline in ipt:
            line_count += 1
            if line_count != 1:
                #if not eachline.startswith('chr'):

                eachline = eachline.strip('\n')

                col = eachline.strip().split()

                if len(col) == 10:

                    chr_nm = col[1]
                    start = col[2]
                    end = col[3]
                    te_nm = col[4]
                    direct = col[5]
                    num_1 = col[6]
                    num_2 = col[7]
                    ratio = col[8]
                    samples = col[9]

                    ##only compare the start region
                    ##only consider the non-reference situation
                    diff = int(end) - int(start)


                    if diff == 1 or diff == 2:
                        id += 1
                        ##store all the line information
                        line_dic[str(id)] = {'chr': chr_nm, 'st': start, 'ed': end,
                                             'num_1': num_1, 'num_2': num_2,'ratio':ratio,
                                             'te_nm':te_nm,'dir':direct,'samples':samples}
                        if id == 1:
                            temp_id_list.append(str(id))
                        ##if id != the final line
                        else:
                            if line_count != (total_line_num + 1):
                            #if id != total_line_num:
                                diff_start = int(start) - int(line_dic[str(temp_id_list[-1])]['st'])
                                if 0 <= diff_start <= int(thresd):
                                    if chr_nm == line_dic[str(temp_id_list[-1])]['chr']:
                                        temp_id_list.append(str(id))
                                    else:
                                        store_id_list.append(temp_id_list)
                                        temp_id_list = []
                                        temp_id_list.append(str(id))
                                else:
                                    ##transfer the temp_id_list to a final list
                                    ##this final list contains all the single id or id list.
                                    store_id_list.append(temp_id_list)
                                    temp_id_list = []
                                    temp_id_list.append(str(id))

                            else:
                                #print('id is the final one')
                                diff_start = int(start) - int(line_dic[str(temp_id_list[-1])]['st'])
                                if 0 <= diff_start and diff_start <= int(thresd):
                                    if chr_nm == line_dic[str(temp_id_list[-1])]['chr']:
                                        temp_id_list.append(str(id))
                                        store_id_list.append(temp_id_list)
                                    else:
                                        store_id_list.append(temp_id_list)
                                        last_list = [str(id)]
                                        store_id_list.append(last_list)
                                else:
                                    store_id_list.append(temp_id_list)
                                    last_list = [str(id)]
                                    store_id_list.append(last_list)
                                    #print(temp_id_list)
                    else:
                        ##updation 8.11
                        ##generate the sample information item
                        ##SRRXXX:ClassIXXX;SRRXXX:ClassIXXXX
                        ##splite the sample string into list
                        sample_te_string = 'Sample_TE_infor'
                        sample_list = samples.split(';')
                        for eachsp in sample_list:
                            single_sample_te_string = eachsp + ':' + te_nm
                            sample_te_string = sample_te_string + ';' + single_sample_te_string

                        final_line = chr_nm + '\t' + str(start) + '\t' + str(end) + '\t' + \
                             str(num_1) + '\t' + str(num_2) + '\t' + str(ratio) + '\t' + 'o' + '\t' + sample_te_string
                        final_line_list.append(final_line)


    return (line_dic,store_id_list,final_line_list)


def combine_close_TE (line_dic,store_id_list,final_line_list):

    ##do not consider the ratio for the combination
    ##the ratio will be selected before this step suing the check_opt_filter_geno.py
    ##do not need to filter, since all the location should be 1
    ##but we can add num_1 and num_2 respectively

    #final_line_list = []

    ##this function combine all close TE and generate the new final line file
    ##combined TEs will be marked with c, otherwise, it will be marked with s
    ##we will calculate new ratio for the combined TE
    for eachid_list in store_id_list:
        ##if no combined id occurred
        if len(eachid_list) == 1:

            ##updation 8.11
            ##generate the sample information
            samples = line_dic[eachid_list[0]]['samples']
            te_nm = line_dic[eachid_list[0]]['te_nm']
            sample_te_string = 'Sample_TE_infor'
            sample_list = samples.split(';')
            for eachsp in sample_list:
                single_sample_te_string = eachsp + ':' + te_nm
                sample_te_string = sample_te_string + ';' + single_sample_te_string

            final_line = line_dic[eachid_list[0]]['chr'] + '\t' + line_dic[eachid_list[0]]['st'] + '\t' + \
                         line_dic[eachid_list[0]]['ed'] + '\t' + line_dic[eachid_list[0]]['num_1'] + '\t' + \
                         line_dic[eachid_list[0]]['num_2'] + '\t' + line_dic[eachid_list[0]]['ratio'] + '\t' + 's' + '\t' + \
                         sample_te_string
            final_line_list.append(final_line)

        else:
            chr_nm = line_dic[eachid_list[0]]['chr']
            start = line_dic[eachid_list[0]]['st']
            end = line_dic[eachid_list[-1]]['ed']

            ##updation 8.11
            ##generate the sample information
            sample_te_string = 'Sample_TE_infor'
            for eachid in eachid_list:
                samples = line_dic[eachid]['samples']
                te_nm = line_dic[eachid]['te_nm']
                sample_list = samples.split(';')
                for eachsp in sample_list:
                    single_sample_te_string = eachsp + ':' + te_nm
                    sample_te_string = sample_te_string + ';' + single_sample_te_string

            ##decide which number should be added together
            ##only if the close TE have same MAF type, otherwise, we will write eachid of this list
            #id_len = len(eachid_list)
            #total_num_1_type_num = 0 ##this calculates the number of num_1_type
            #total_num_2_type_num = 0 ##this calculates the number of num_2_type

            new_num_1 = 0
            for eachid in eachid_list:
                id_num_1 = int(line_dic[eachid]['num_1'])
                new_num_1 = new_num_1 + id_num_1

            new_num_2 = 0
            for eachid in eachid_list:
                id_num_2 = int(line_dic[eachid]['num_2'])
                new_num_2 = new_num_2 + id_num_2

            #new_num_2 = int(total_sample_num) - int(new_num_1)
            #new_ratio = float(new_num_1/new_num_2)
            final_line = chr_nm + '\t' + str(start) + '\t' + str(end) + '\t' + \
                         str(new_num_1) + '\t' + str(new_num_2) + '\t' + str('na') + '\t' + 'c' + '\t' + sample_te_string
            final_line_list.append(final_line)

    return (final_line_list)



def not_combine_close_TE (line_dic,store_id_list,final_line_list):

    ##this function combine all close TE and generate the new final line file
    ##combined TEs will be marked with c, otherwise, it will be marked with s
    ##we will calculate new ratio for the combined TE
    for eachid_list in store_id_list:
        ##if no combined id occurred
        if len(eachid_list) == 1:

            ##updation 8.11
            ##generate the sample information
            samples = line_dic[eachid_list[0]]['samples']
            te_nm = line_dic[eachid_list[0]]['te_nm']
            sample_te_string = 'Sample_TE_infor'
            sample_list = samples.split(';')
            for eachsp in sample_list:
                single_sample_te_string = eachsp + ':' + te_nm
                sample_te_string = sample_te_string + ';' + single_sample_te_string

            final_line = line_dic[eachid_list[0]]['chr'] + '\t' + line_dic[eachid_list[0]]['st'] + '\t' + \
                         line_dic[eachid_list[0]]['ed'] + '\t' + line_dic[eachid_list[0]]['num_1'] + '\t' + \
                         line_dic[eachid_list[0]]['num_2'] + '\t' + line_dic[eachid_list[0]]['ratio'] + '\t' + 's' + '\t' + \
                         sample_te_string
            final_line_list.append(final_line)

        else:

            ##updation 10.22 do not combine close loci
            ##updation 8.11
            ##generate the sample information

            for eachid in eachid_list:

                #chr_nm = line_dic[eachid]['chr']
                #start = line_dic[eachid]['st']
                #end = line_dic[eachid]['ed']

                samples = line_dic[eachid]['samples']
                te_nm = line_dic[eachid]['te_nm']
                sample_te_string = 'Sample_TE_infor'
                sample_list = samples.split(';')
                for eachsp in sample_list:
                    single_sample_te_string = eachsp + ':' + te_nm
                    sample_te_string = sample_te_string + ';' + single_sample_te_string

                final_line = line_dic[eachid_list[0]]['chr'] + '\t' + line_dic[eachid_list[0]]['st'] + '\t' + \
                             line_dic[eachid_list[0]]['ed'] + '\t' + line_dic[eachid_list[0]]['num_1'] + '\t' + \
                             line_dic[eachid_list[0]]['num_2'] + '\t' + line_dic[eachid_list[0]]['ratio'] + '\t' + 's' + '\t' + \
                             sample_te_string
                final_line_list.append(final_line)

    return (final_line_list)

##updation 050220
##calculate the proportion of sample to replace the na in the opt_modify_filter_TE_loc.txt
def calcualte_pro_to_replace_na_in_combine_loc (final_line_list,sample_number):

    store_new_final_line_list = []
    for eachline in final_line_list:
        col = eachline.strip().split()
        if col[5] == 'na':
            ##check the sample_TE_infor to only identify the number of unique samples
            sp_col = col[7].split(';')

            store_sp_dic = {}
            for i in range(1,len(sp_col)):
                mt = re.match('(.+):.+',sp_col[i])
                sp_nm = mt.group(1) ##SRRXXX
                store_sp_dic[sp_nm] = 1

            sample_num = len(list(store_sp_dic.keys()))

            pro = sample_num/int(sample_number)

            new_line = col[0] + '\t' + col[1] + '\t' + col[2]  + '\t' + col[3]  + '\t' + col[4] + \
                       '\t' + str(pro) + '\t' + col[6] + '\t' + col[7]
            store_new_final_line_list.append(new_line)
        else:
            store_new_final_line_list.append(eachline)

    return (store_new_final_line_list)

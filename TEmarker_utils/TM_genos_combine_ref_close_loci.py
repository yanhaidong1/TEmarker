#!/usr/bin/env python

##updation 122120 change the c and s cover case
##updation 121920 add the case when first line is the last one
##updation 110220this version did not consider the s or c is covered by the reference we need to have a check
##we need to generate a dic to store the s or c TEs that are covered by the o type and finally, we will filter out these locations

##this script we will combine the close loci with the overlapped regions

##import modules
import re

def combine_close_loci (genotype_fl,working_dir):

    ########
    ##step 1: store the id line and store the same o to one dic

    ##generate a temp genotype file with id information
    store_id_line_list = []
    count = 0
    with open (genotype_fl,'r') as ipt:
        for eachline in ipt:
            count += 1
            eachline = eachline.strip('\n')
            new_line = str(count) + '\t' + eachline
            store_id_line_list.append(new_line)

    with open (working_dir + '/temp_te_genotype_annot_add_id.txt','w+') as opt:
        for eachline in store_id_line_list:
            opt.write(eachline + '\n')

    ##initial a list contain all the TE situations
    dic_list = []
    ##initial a dictionary to store the target lines
    dic_te = {}

    ##get the last id
    with open(working_dir + '/temp_te_genotype_annot_add_id.txt', 'r') as ipt_rmk_out:
        last_line = ipt_rmk_out.readlines()[-1]
        last_col = last_line.strip().split()
        last_id = last_col[0]


    ##updation 110220
    store_covered_c_with_o_loc_dic = {}
    ##store the te infor
    with open(working_dir + '/temp_te_genotype_annot_add_id.txt', 'r') as ipt_rmk_out:

        for line in ipt_rmk_out:
            col = line.strip().split()

            chr = col[1]
            #te_nm = col[4]
            id = col[0]
            bg = col[2]
            ed = col[3]
            #dir = col[3]
            #lib_bg = col[5]
            #lib_ed = col[6]
            #lib_left = col[7]
            comb_type = col[4]

            if id == str(1):  ##if the id is 1, it will directly store in the dic_te
                dic_te[id] = {'chr': chr, 'begin': bg, 'end': ed, 'comb_type':comb_type, 'line':line}
                if id == last_id:  ##should be out of the previous loop
                    dic_list.append(dic_te)

            else:  ##if the id is over 1
                ##if the chr is the same as the previous one
                if dic_te[str(int(id) - 1)]['chr'] == chr:

                    if dic_te[str(int(id) - 1)]['comb_type'] == comb_type:  ##if the comb_type is the same
                        if comb_type  == 'o':
                            ##detect whether they are overlapped
                            pre_st = dic_te[str(int(id) - 1)]['begin']
                            pre_ed = dic_te[str(int(id) - 1)]['end']

                            ##it means there is a overlap
                            if int(ed) >= int(pre_st) and int(bg) <= int(pre_ed):
                                dic_te[id] = {'chr': chr, 'begin': bg, 'end': ed, 'comb_type': comb_type, 'line': line}

                                if id == last_id:  ##should be out of the previous loop
                                    dic_list.append(dic_te)

                            else: ##if there is no overlap we will store the previous dic and assign a new dic id
                                dic_list.append(dic_te)
                                dic_te = {}
                                dic_te[id] = {'chr': chr, 'begin': bg, 'end': ed, 'comb_type': comb_type, 'line': line}

                                if id == last_id:  ##should be out of the previous loop
                                    dic_list.append(dic_te)

                        else:
                            ##directly store the previous dic and assign a new dic id
                            dic_list.append(dic_te)
                            dic_te = {}
                            dic_te[id] = {'chr': chr, 'begin': bg, 'end': ed, 'comb_type': comb_type, 'line': line}

                            if id == last_id:  ##should be out of the previous loop
                                dic_list.append(dic_te)

                    ##if the comb_type is not the same
                    ##we also need to store
                    else:

                        ##updation 110520
                        ##it means the pre comb_type is o and current is c or s
                        pre_comb_type = dic_te[str(int(id) - 1)]['comb_type']
                        if pre_comb_type == 'o':
                            ##the current is s or c

                            ##or the pre comb_type is c or s and current o
                            ##updation 110220 we need to check whether the c is covered with the o
                            ##detect whether they are overlapped
                            pre_st = dic_te[str(int(id) - 1)]['begin']
                            pre_ed = dic_te[str(int(id) - 1)]['end']
                            ##it means there is a overlap
                            if int(ed) >= int(pre_st) and int(bg) <= int(pre_ed):
                                loc_str = chr + '_' + bg + '_' + ed
                                store_covered_c_with_o_loc_dic[loc_str] = 1
                        else:
                            ##updation 122120
                            ##if the pre is s or c
                            ##the current could be o
                            ##or could be s or c

                            if comb_type == 'o':
                                ##if the pre is s or c
                                ##the current is o
                                pre_st = dic_te[str(int(id) - 1)]['begin']
                                pre_ed = dic_te[str(int(id) - 1)]['end']
                                if int(ed) >= int(pre_st) and int(bg) <= int(pre_ed):
                                    loc_str = chr + '_' + pre_st + '_' + pre_ed
                                    store_covered_c_with_o_loc_dic[loc_str] = 1

                            ##if current is not o so it would be c or s
                            ##the reason to case this case is because we enlarge the searching range from s and c case
                            ##and allow there are some overlapping
                            ##in this case, we need to follow the single case TE since the combined case

                            ##there is no else since we will modify the searching range that allows there is no cover for the s and c
                            #else:

                        dic_list.append(dic_te)
                        dic_te = {}
                        dic_te[id] = {'chr': chr, 'begin': bg, 'end': ed, 'comb_type': comb_type, 'line': line}

                        if id == last_id:  ##should be out of the previous loop
                            dic_list.append(dic_te)

                else:  ##if the chr is not the same as the previous one
                    ##store the dic_te which has been stored the te informations
                    dic_list.append(dic_te)
                    dic_te = {}
                    dic_te[id] = {'chr': chr, 'begin': bg, 'end': ed, 'comb_type': comb_type, 'line': line}

                    if id == last_id:  ##should be out of the previous loop
                        dic_list.append(dic_te)


    #print(dic_list)
    print(store_covered_c_with_o_loc_dic)

    ########
    ##step 2: combine the overlapped o ref loci
    store_final_line_list = []
    ##updation 110220 write a file to show whether we remove some locations for the c and s
    store_remove_c_s_line_list = []
    #comb_count = 0
    for each_te_dic in dic_list:

        if len(each_te_dic.keys()) == 1:

            final_line_no_id_str = ''
            for eachid in each_te_dic: ##key is id 1,2,3,4,5
                line = each_te_dic[eachid]['line']
                col_line = line.strip().split()
                first_item = col_line[1]
                final_line_no_id_str = first_item

            for eachid in each_te_dic: ##key is id 1,2,3,4,5
                line = each_te_dic[eachid]['line']
                col_line = line.strip().split()
                for i in range (2,len(col_line)):
                    final_line_no_id_str = final_line_no_id_str + '\t' + col_line[i]

            ##updation 110220 check location of c or s TE
            wrong_count = 0
            for eachid in each_te_dic:
                chr = each_te_dic[eachid]['chr']
                bg = each_te_dic[eachid]['begin']
                ed = each_te_dic[eachid]['end']
                loc_infor = chr + '_' + bg + '_' + ed
                if loc_infor in store_covered_c_with_o_loc_dic:
                    wrong_count += 1

            if wrong_count == 0:
                store_final_line_list.append(final_line_no_id_str)
            else:
                store_remove_c_s_line_list.append(final_line_no_id_str)


        else:
            #comb_count += 1
            #print(len(each_te_dic.keys()) )

            #print(each_te_dic)
            ##it means we need to decide a new o locus
            ##first we need to select the smallest bg
            bg_list = []
            ed_list = []
            chr = ''

            ##since we need to generate a new genotype line we need to extract the the first several col information
            ori_id = ''
            #ori_num = ''
            #ori_total = ''
            #ori_pro = ''
            ori_geno = 'o'
            #ori_sap_infor_str = ''
            ori_te_nm = ''

            for eachid in each_te_dic:
                bg = int(each_te_dic[eachid]['begin'])
                bg_list.append(bg)

                ed = int(each_te_dic[eachid]['end'])
                ed_list.append(ed)

                chr = each_te_dic[eachid]['chr']

                ori_line = each_te_dic[eachid]['line']
                ori_line_col = ori_line.split()
                ori_id = ori_line_col[0]
                #ori_num = ori_line_col[4]
                #ori_total = ori_line_col[5]
                #ori_pro =  ori_line_col[6]
                ori_te_nm = ori_line_col[5]
                #ori_sap_infor_str = ori_line_col[8]


            smallest_bg = min(bg_list)
            largest_ed = max(ed_list)

            ##second we need to store the name of the sample
            store_sp_name_dic = {}
            for eachid in each_te_dic:

                loc_line = each_te_dic[eachid]['line']
                loc_col = loc_line.split()

                for i in range(6,len(loc_col)):
                    mt = re.match('(.+):.+',loc_col[i])
                    sp_nm = mt.group(1)
                    store_sp_name_dic[sp_nm] = 1

            #print(store_sp_name_dic)


            ##third we need analyze on each sample
            store_sp_str = ''  ##DRR054229:0.0;0/15;0/0;Unknown_TE       DRR054234:0.0;0/4;0/0;Unknown_TE
            for eachsp in store_sp_name_dic:

                #store_geno_value_list = [] ##1/1: 2, 0/1:1, 0/0:0
                ##then we need to compare the value of them and select the largest one
                sp_id = 0 ##
                store_same_sp_loc_dic = {}
                ##key is the sp_id (0 or 1 or 2...) and value has two part, first is the geno line (eg. DRR054241:0.0;0/7;0/0;Unknown_TE) and second is geno_value
                for eachid in each_te_dic:
                    loc_line = each_te_dic[eachid]['line']
                    loc_col = loc_line.split()
                    for i in range(6, len(loc_col)):


                        mt = re.match('(.+):.+', loc_col[i])
                        sp_nm = mt.group(1)

                        if eachsp == sp_nm:

                            sp_id += 1
                            ##so the sp_id is the location id since if we allow the eachsp == sp_nm so there is only one value from 6 to len(loc_col)

                            ##since this is the genos file so there is no missing information
                            geno_line = loc_col[i]
                            geno_col = geno_line.split(';')
                            geno = geno_col[2]
                            geno_value = ''
                            if geno == '0/0':
                                geno_value = 0
                            if geno == '0/1':
                                geno_value = 1
                            if geno == '1/1':
                                geno_value = 2

                            #store_geno_value_list.append(geno_value)

                            store_same_sp_loc_dic[str(sp_id) + '_' + geno_line] = geno_value

                ##eg 1_DRR054229:0.0;0/15;0/0;LTR
                max_geno_id = max(store_same_sp_loc_dic, key=store_same_sp_loc_dic.get)
                #print(max_geno_id)
                mt = re.match('.+?_(.+)',max_geno_id)
                real_id = mt.group(1) ##DRR054229:0.0;0/15;0/0;LTR
                #print(real_id)
                store_sp_str = store_sp_str + '\t' + real_id

            ##generate the final line
            final_line = chr + '\t' + str(smallest_bg) + '\t' + str(largest_ed) + '\t' + 'o' + '\t' + ori_te_nm + store_sp_str
            store_final_line_list.append(final_line)

    return (store_final_line_list,store_remove_c_s_line_list)















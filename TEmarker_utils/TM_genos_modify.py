#!/usr/bin/env python

##import modules
import re
import glob

def modification (input_classify_hetero_file,input_CNV_dir,Her_thd):

    ##do not consider the 0 reads insisting the location
    ##only consider location which is 1 bp or 2 bp difference
    ##later we will change the range after combine the close location together.

    ##generate dic store the sample name and sample file location
    sample_file_loc_dic = {}
    cnv_file_list = glob.glob(input_CNV_dir + '/*')
    for eachfl in cnv_file_list:
        ##get the sample name
        mt = re.match('.+/(.+)',eachfl)
        fl_nm = mt.group(1)
        mt = re.match('(.+)_genotype.txt',fl_nm)
        sample_nm = mt.group(1)
        sample_file_loc_dic[sample_nm] = eachfl


    ##store the final_output
    final_line_list = []
    count = 0
    with open (input_classify_hetero_file,'r') as ipt:
        for eachline in ipt:

            count += 1
            print('the analzyed te is ' + str(count))

            eachline = eachline.strip('\n')
            col = eachline.strip().split()

            heter_start = col[1]
            heter_end = col[2]
            heter_chr = col[0]

            loc_type = col[6]


            new_line_string = heter_chr + '\t' + heter_start + '\t' + heter_end + \
                              '\t' + col[3] + '\t' + col[4] + '\t' + col[5]

            ##updation8.9
            if loc_type == 'c' or loc_type == 's':

            #if (int(heter_end) - int(heter_start)) == 1 or (int(heter_end) - int(heter_start)) == 2:

                sample_list = col[8:]  ##updation8.13 change 7 to 8 since we added the specific family information
                ##detect if sample has reads supporting the insertions
                for eachsp_line in sample_list:
                    ##get the supported read information
                    col_sp = eachsp_line.split(';')
                    read_item = col_sp[1]
                    te_fam = col_sp[3]
                    ##get the reads
                    mt = re.match('(.+)/(.+)',read_item)
                    read_num = mt.group(1)
                    total_read_num = mt.group(2)

                    if int(read_num) != 0:
                        ##consider if there is CNV event on the TE region
                        ##search all the files in the CNV dir to find the
                        ##get the name of the sample
                        mt = re.match('(.+)\:(.+)',eachsp_line)
                        sample_nm = mt.group(1)

                        ##detect if the CNV file has covered the insertion location
                        sample_fl_path = sample_file_loc_dic[sample_nm]

                        final_copy_num = 1
                        with open (sample_fl_path,'r') as ipt:
                            for eachline in ipt:
                                eachline = eachline.strip('\n')
                                col_cnv_fl = eachline.strip().split()
                                cnv_loc = col_cnv_fl[1]
                                mt = re.match('(.+)\:(.+)\-(.+)',cnv_loc)
                                cnv_chr = mt.group(1)
                                cnv_start = mt.group(2)
                                cnv_end = mt.group(3)
                                copy_num = col_cnv_fl[3]

                                ##check if the heter location is within the cnv location
                                if heter_chr == cnv_chr:
                                    if int(cnv_start) <= int(heter_start) and int(cnv_end) >= int(heter_end):
                                        if float(copy_num) > float(0):
                                            final_copy_num = copy_num

                        ##normalize the total read number
                        normalized_total_read_num = float(total_read_num)/float(final_copy_num)

                        ##get the new proportion of the read ratio
                        #original_ratio = float(int(read_num)/int(total_read_num))
                        normalized_ratio = float(read_num)/normalized_total_read_num

                        ##use the previous genotyping.
                        ##later we will revise it
                        if (1-float(Her_thd)) <= float(normalized_ratio) and float(normalized_ratio) <= float(Her_thd):
                            new_line_string = new_line_string + '\t' + sample_nm + ':' + str(normalized_ratio) + ';' + \
                                           str(int(read_num)) + '/' + str(normalized_total_read_num) + ';' + '0/1' + ';' + te_fam

                        ##the same as the consensuse TEs, this will be the 0/0
                        ##updation8.9 change the 0/0 to 1/1, all should be followed by reference genome
                        if float(Her_thd) < float(normalized_ratio) and float(normalized_ratio) <= 1.0:
                            new_line_string = new_line_string + '\t' + sample_nm + ':' + str(normalized_ratio) + ';' + \
                                           str(int(read_num)) + '/' + str(normalized_total_read_num) + ';' + '1/1' + ';' + te_fam
                        ##different as the consensus TEs, this will be the 1/1
                        ##updation8.9 change the 1/1 to 0/0, all should be followed by reference genome
                        if 0 <= float(normalized_ratio) and float(normalized_ratio) < (1-float(Her_thd)):
                            new_line_string = new_line_string + '\t' + sample_nm + ':' + str(normalized_ratio) + ';' + \
                                           str(int(read_num)) + '/' + str(normalized_total_read_num) + ';' + '0/0' + ';' + te_fam

                        ##if normalized_ratio > 1, we will not normalize
                        if float(normalized_ratio) > float(1):
                            new_line_string = new_line_string + '\t' + eachsp_line

                    else:
                        new_line_string = new_line_string + '\t' + eachsp_line

                final_line_list.append(new_line_string)

            else:
                final_line_list.append(eachline)

    return (final_line_list)























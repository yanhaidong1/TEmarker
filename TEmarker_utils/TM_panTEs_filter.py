#!/usr/bin/env python

def check_filter_opt (input_file,fth):

    store_line_list = []
    with open (input_file, 'r') as ipt:
        for eachline in ipt:
            eachline = eachline.strip('\n')
            col = eachline.strip().split()
            #num = col[5]
            ##updation 8.26 only consider the int(col[3]) > 5, and do not filter out col[4]
            ##since some location may have all the non-reference insertion
            ##na means the combination case
            ##updation 050220 no na, since we already transfer na to a proportion
            ##so we directly use fth instead of fcth and fncth
            #if col[5] == 'na' and int(col[3]) > int(fcth):
            #if col[5] == 'na' and int(col[3]) > 5 and int(col[4]) > 5:
            #    store_line_list.append(eachline)
            #if col[5] != 'na':
            if float(col[5]) > float(fth):
                store_line_list.append(eachline)

    return (store_line_list)





#!/usr/bin/env python

##updating 033022 we will use the cyverse to download data
##this script will download the testing data for step 1,2,3
##Also download the testing data for step 0

##BUILT-IN MODULES
import argparse
import sys
import os
import subprocess
from distutils.spawn import find_executable

##SCRIPTS
def get_parsed_args():
    parser = argparse.ArgumentParser(description="Download testing data")

    ##############################
    ##parse the required arguments
    parser.add_argument("-o", dest='output_dir', default="./", help="Output directory"
                                                                    "Default: ./ ")

    parser.add_argument('-s0', dest='s0_yes', help="If users initiate this argument, they will download testing data for step 0.")

    parser.add_argument('-s123', dest='s123_yes', help="If users initiate this argument, they will download all the testing data for step123.")

    ##parse of parameters
    args = parser.parse_args()
    return args


##define a function to download testing data
def download_data_step123 (store_download_dir,name,cyverse_path):

    cmd = 'wget ' + cyverse_path + ' -O ' + store_download_dir + '/' + name + '.zip'
    subprocess.call(cmd,shell=True)

    #cmd = "wget --load-cookies /tmp/cookies.txt \"https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate " + \
    #      "\'https://docs.google.com/uc?export=download&id=" + google_drive_path + "\' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\\1\\n/p')&id=" + google_drive_path + "\" -O " + \
    #      store_download_dir + "/" + name + ".zip " + "&& rm -rf /tmp/cookies.txt"
    #subprocess.call(cmd, shell=True)

    #cmd = 'unzip ' +  store_download_dir + '/' + name + '.zip -d ' + store_download_dir + '/'
    #subprocess.call(cmd,shell=True)

    #cmd = 'rm ' + store_download_dir + '/' + name + '.zip'
    #subprocess.call(cmd,shell=True)

    #cmd = 'rm -rf ' + store_download_dir + '/__MACOSX' + ' ' + store_download_dir + '/' + name + '/Icon'
    #subprocess.call(cmd,shell=True)

def download_data_step0 (store_download_dir,name,cyverse_path):

    #cmd = "wget --load-cookies /tmp/cookies.txt \"https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate " + \
    #      "\'https://docs.google.com/uc?export=download&id=" + google_drive_path + "\' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\\1\\n/p')&id=" + google_drive_path + "\" -O " + \
    #      store_download_dir + "/" + name + ".tar.gz " + "&& rm -rf /tmp/cookies.txt"
    #subprocess.call(cmd, shell=True)

    cmd = 'wget ' + cyverse_path + ' -O ' + store_download_dir + '/' + name + '.zip'
    subprocess.call(cmd,shell=True)



def main(argv=None):
    if argv is None:
        argv = sys.argv
    args = get_parsed_args()

    output_dir = args.output_dir
    if not output_dir.endswith('/'):
        output_dir = output_dir + '/'
    else:
        output_dir = output_dir

    ##if users want to download testing data for Step 1,2,3
    if args.s123_yes is not None:
        cyverse_path = 'https://data.cyverse.org/dav-anon/iplant/home/yanhaidong1991/Example_dir_step123.zip'
        download_data_step123(output_dir,'Example_dir_step123',cyverse_path)

    if args.s0_yes is not None:
        cyverse_path = 'https://data.cyverse.org/dav-anon/iplant/home/yanhaidong1991/Example_dir_step0.tar.gz'
        download_data_step0(output_dir, 'Example_dir_step0', cyverse_path)


if __name__ == "__main__":
    main()




'''
This python module contains functions, that process (save, call ymap upon it, 
zip the output) a String or a plain Txt input file recieved via HTML forms.
'''

from django.http import HttpResponse
from django.shortcuts import render
from django.http import Http404
from .forms import *
from django.http import HttpResponseRedirect
from shutil import make_archive
from numba.pycc.decorators import process_input_files
import logging
import re
import csv
import os.path
import sys
from .ymap import *
import re
from celery import task
import datetime

'''
######################################### Setting up logger (For debugging in command line) ###################
'''

'''
usage: logger.debug("string message to print in the command line, when execution hits this line.")
'''
# set up logger, even if imported and not running as __main__
FORMAT = '%(asctime)s %(levelname)s:%(name)s:%(message)s'
VERBOSITY = os.environ.get('YMAP_SERVER_VERBOSITY', 1)

if int(VERBOSITY) in [1, 2, 3, 4, 5]:
    logging.basicConfig(format=FORMAT, level=int(
            VERBOSITY) * 10, datefmt='%Y-%m-%d %H:%M:%S')
    logger = logging.getLogger(__name__)
else:
    logging.basicConfig(format=FORMAT, level=logging.INFO,
                        datefmt='%Y-%m-%d %H:%M:%S')
    logger = logging.getLogger(__name__)
    logger.warning("Unrecognized verbosity level " + VERBOSITY +
                ". Set verbosity from 1 as DEBUG to 5 as CRITICAL. Now, using 2 (INFO) instead.")


#to do: show results on an HTML page


@task()
def validate_file_type(title):
    '''
    Checks if the file's extension is .txt.
    Warning! This extension check may not be sufficient.
    If you want to perform more thorough file format checks, try magic, django-clamav-upload
    '''
    return bool(re.match('.+\.txt$', title))


  
@task()
def is_protein(file):
    '''
    This function returns true if the tab delimited file has 2 columns
    @param file: should be given as: path + filename + extension. file to analyse.
    '''  
    num_cols = 0
    with open(file) as f:
        reader = csv.reader(f, delimiter='\t', skipinitialspace=True)
        first_row = next(reader)
        num_cols = len(first_row)
    logger.debug("number of columns in the uploaded file: " + str(num_cols))
    return bool(num_cols is 2)


@task()  
def is_gene(file):
    '''
    This function returns true if the tab delimited file has 5 columns
    @param file: should be given as: path + filename + extension. file to analyse.
    '''
    num_cols = 0
    with open(file) as f:
        reader = csv.reader(f, delimiter='\t', skipinitialspace=True)
        first_row = next(reader)
        num_cols = len(first_row)
    logger.debug("number of columns in the uploaded file: " + str(num_cols))
    return bool(num_cols is 5)



@task()
def wipe_folder(path):
    '''
    This function erases all files and all folders (recursively) from a given path
    @param path: folder to wipe clean
    '''
    folder = path
    if os.path.isdir(folder):
        for the_file in os.listdir(folder):
            file_path = os.path.join(folder, the_file)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path): shutil.rmtree(file_path)
            except Exception as e:
                print(e)

@task()
def clean_up(input_path, output_path, archive_path):
    '''
    Emtpies the folders that are given as argument
    '''

    try:

        wipe_folder(input_path)
        wipe_folder(output_path)
        wipe_folder(archive_path)

        if os.path.isdir(input_path):
            os.rmdir(input_path)
        
        if os.path.isdir(output_path):
            os.rmdir(output_path)
        
        if os.path.isdir(archive_path):
            os.rmdir(archive_path)

    except Exception as e:
        print(e)

    logger.debug("hello from cleanup")

@task()            
def run_yproteins():
    '''
    runs ymap_proteins from ymap. (it is wrapped so it the function in ymap changes, it has to only be changed here)
    '''
    ymap_proteins()

@task()    
def run_ygenes():
    '''
    runs ymap_genes from ymap. (it is wrapped so it the function in ymap changes, it has to only be changed here)
    '''
    ymap_genes()
    
@task()
def is_too_old(file, age_min = 15):
    '''
    Checks whether a given input file os older than X minutes. Default: 15 minutes.
    @param file - absolute path to the file including filename
    @param age_min - the maximum allowed age of file in minutes (default: 15 minutes)
    @return True, if the file is older than age_min minutes, False otherwise
    '''

    if sys.version_info.minor > 5:
        time_stamp = os.path.getmtime(file)
    else:
        time_stamp = os.path.getmtime(str(file))
    file_age = datetime.datetime.fromtimestamp(time_stamp)
    accepted_minutes_ago = datetime.datetime.now()-datetime.timedelta(minutes = age_min)
    return file_age < accepted_minutes_ago

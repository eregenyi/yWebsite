'''
This python module contains functions, that process (save, call ymap upon it, 
zip the output) a String or a plain Txt input file recieved via HTML forms.
To Do:
'''

from django.http import HttpResponse
from django.shortcuts import render
from django.http import Http404
from .forms import *
from django.http import HttpResponseRedirect
from .process_input import *
from shutil import make_archive
from numba.pycc.decorators import process_input_files
import logging
import re
from _overlapped import NULL
import csv
import os.path
from .ymap import *
import re


wd = os.getcwd()
input_path = os.path.join(wd, 'ymap_webtool', 'data', 'input')
output_path = os.path.join(wd, 'ymap_webtool', 'data', 'output') 
archive_path = os.path.join(wd, 'ymap_webtool', 'data', 'archive') 
input_file_name = "mutation.txt"
output_file_name = "results.zip"
archive_name = "results"
download_name = "results.zip"
gene_level_input_name = "mutated_proteins.txt"

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

'''
Checks if the file's extension is .txt.
Warning! This extension check may not be sufficient.
If you want to perform more thorough file format checks, try magic, django-clamav-upload
'''

def validate_file_type(title):
    return bool(re.match('.+\.txt$', title))


'''
This function returns true if the tab delimited file has 2 columns
@param file: should be given as: path + filename + extension. file to analyse.
'''    

def is_protein(file):
    num_cols = 0
    with open(file) as f:
        reader = csv.reader(f, delimiter='\t', skipinitialspace=True)
        first_row = next(reader)
        num_cols = len(first_row)
    logger.debug("number of columns in the uploaded file: " + str(num_cols))
    return bool(num_cols is 2)

'''
This function returns true if the tab delimited file has 5 columns
@param file: should be given as: path + filename + extension. file to analyse.
'''
    
def is_gene(file):
    num_cols = 0
    with open(file) as f:
        reader = csv.reader(f, delimiter='\t', skipinitialspace=True)
        first_row = next(reader)
        num_cols = len(first_row)
    logger.debug("number of columns in the uploaded file: " + str(num_cols))
    return bool(num_cols is 5)


'''
This function erases all files and all folders (recursively) from a given path
@param path: folder to wipe clean
'''
def wipe_folder(path):
    folder = path
    for the_file in os.listdir(folder):
        file_path = os.path.join(folder, the_file)
        try:
            if os.path.isfile(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path): shutil.rmtree(file_path)
        except Exception as e:
            print(e)

def clean_up(input_path, output_path, archive_path):
    input_path = input_path
    output_path = output_path
    archive_path = archive_path
    #erase the input, output, but keep the archive under another name
    wipe_folder(input_path)
    wipe_folder(output_path)
    logger.debug("hello from cleanup")
    #os.rename(os.path.join(archive_path, archive_name + ".zip"), os.path.join(archive_path, "archived" + str(i) + ".zip"))

def save_string(string, path, name):
    '''
    string: query text provided by user, to save as a text file
    path: where should the file be saved on the server side
    name: the file is going to be saved on the server side under this name (with extension)
    '''
    with open(os.path.join(path, name), 'w') as f:
        f.write(string)
        
def save_file(file, path, name):
    '''
    file: file to save
    path: where should the file be saved on the server side
    name: the file is going to be saved on the server side under this name (with extension)
    '''
    f = file
    with open(os.path.join(path, name), 'wb+') as destination:
        for chunk in f.chunks():
            destination.write(chunk)
            
def run_yproteins():
    ymap_proteins()
    
def run_ygenes():
    ymap_genes()

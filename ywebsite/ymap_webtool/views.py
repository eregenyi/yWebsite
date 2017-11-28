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


wd = os.getcwd()
input_path = os.path.join(wd, 'ymap_webtool', 'data', 'input')
output_path = os.path.join(wd, 'ymap_webtool', 'data', 'output') 
input_file_name = "mutation.txt"
output_file_name = "results.zip"
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
############################################### Views that are added to URLs #########################################
'''

def index(request):
    #specify the template to the Main page here. Django by default looks for templates in the templates folder
    return render(request, 'ymap_webtool/index.html')

def manual(request):
    #specify the template to the Manual/help page here.
    return render(request, 'ymap_webtool/manual.html')

def about(request):
    #specify the template to the About page here
    return render(request, 'ymap_webtool/about.html')

def publications(request):
    #Specify the template to the Publications page
    return render(request, 'ymap_webtool/publications.html')

def finished(request):
    #Specify template to the page to display after submission
    logger.debug("time to give back the finished page with the script for downloading the results")
    return render(request, 'ymap_webtool/finished.html')

def submission_fail(request):
    #Specify template to the page to display if submission fails
    return render(request, 'ymap_webtool/submission_fail.html')


'''
Checks whether the input file contains protein or gene level mutations, and runs a ymap function accordingly to map.
zips the results
deletes the original results, keeps only the zip file ?
@return redirects to ymap_webtool/finished.
'''
def processing(request):
    logger.debug("started processing the job")
    if is_protein(os.path.join(input_path, input_file_name)) is True:
        logger.debug("it was a protein input file. running yproteins now")       
        #run_yproteins()
    elif is_gene(os.path.join(input_path, input_file_name)) is True:
        logger.debug("it was a gene input file. running ygenes now")
        #os.rename(os.path.join(input_path, input_file_name), os.path.join(input_path, gene_level_input_name))
        logger.debug("now the file " + input_file_name + " is renamed to " + gene_level_input_name)
        #run_ygenes()
        logger.debug("finished processing the job")
    else:
        logger.debug("ERROR! Neither gene nor protein file!")
        return HttpResponseRedirect('submission_fail')
    
    #make_archive(os.path.join(output_path,'results'), 'zip', root_dir = output_path)
    logger.debug("made an archive")
    # wipe original results folder
    # maybe keep the archive?
    return HttpResponseRedirect('finished')


def submitted(request):
    logger.debug("the job should be submitted and passed to processing")
    return render(request, 'ymap_webtool/submitted.html')


'''
takes a string entered through the HTML form, 
validates the form
if the for is valid, saves it into a text file and redirects to the processing page
else it assignes an empty object to the form and redirects to the submission_fail page
@param request
'''
def get_string(request):
    # if this is a POST request we need to process the form data
    if request.method == 'POST':
        # create a form instance and populate it with data from the request:
        form = StringForm(request.POST)
        # check whether it's valid:
        if form.is_valid():
            # process the data in form.cleaned_data as required
            # redirect to a new URL:
            save_string(form.cleaned_data.get('query'), input_path, input_file_name)
            #Process data here with a separate pydev class?
            return HttpResponseRedirect('submitted')

    # if a GET (or any other method) we'll create a blank form
    else:
        form = StringForm()
    return HttpResponseRedirect('submission_fail')

'''
takes a file that has been uploaded with the HTML form, 
validates its extension (only .txt is allowed!)
if the extension is fine, saves that file under the input_path specified above, and redirects to the processing page
else it redirects to the submission_fail page
@param request
'''
def get_file(request):
    if request.method == 'POST':
        form = UploadFileForm(request.POST, request.FILES)
        if validate_file_type(request.FILES['myfile'].name) is True:
            save_file(request.FILES['myfile'], input_path, input_file_name)
            logger.debug("get_file ran succcesfully")
            return HttpResponseRedirect('submitted')
        else:
            return HttpResponseRedirect('submission_fail')
    

def download_result(request):
    #download the results to the user
    resp_file = os.path.join(wd, output_path, output_file_name)
    response = HttpResponse(resp_file, content_type='application/force-download')
    response['Content-Disposition'] = 'attachment; filename="%s"' % download_name
    logger.debug("download the zip file")
    return response    

def test_page(request):
    '''
    This view is for any functionalities that temporarily need to be tested.Shall be deleted upon completion of the app
    '''
    logger.debug("Print this to the command line please")
    response = HttpResponse('Nothing to test now') 
    return response

'''
######################################### Functions ###################################
'''

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


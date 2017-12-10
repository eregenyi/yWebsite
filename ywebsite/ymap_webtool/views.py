from django.http import HttpResponse
from django.shortcuts import render
from django.http import Http404
from .forms import *
from django.http import HttpResponseRedirect
from .tasks import *
from shutil import make_archive
from numba.pycc.decorators import process_input_files
import logging
import re
from _overlapped import NULL
import csv
import os.path
from pathlib import Path
from django.core.mail import send_mail , BadHeaderError  # import to send an email @test
from django.shortcuts import render_to_response # import to send an email @test
from django.template import RequestContext, Context # import to send an email @test
from django.template.loader import get_template # import to send an email @test
from django.core.mail.message import EmailMessage 

wd = os.getcwd()
input_path = os.path.join(wd, 'ymap_webtool', 'data', 'input')
output_path = os.path.join(wd, 'ymap_webtool', 'data', 'output')
archive_path = os.path.join(wd, 'ymap_webtool', 'data', 'archive')
input_file_name = "mutation.txt"
output_file_name = "results.zip"
archive_name = "results"
download_name = "results.zip"
gene_level_input_name = "mutated_proteins.txt"
i = 0



'''
######################################### Setting up logger (For debugging in command line) ########################
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
	logger.warning("Unrecognized verbosity level " + VERBOSITY + ". Set verbosity from 1 as DEBUG to 5 as CRITICAL. Now, using 2 (INFO) instead.")


# to do: show results on an HTML page

'''
############################################### Views that are added to URLs #########################################
'''

def index(request):
<<<<<<< HEAD
    #specify the template to the Main page here. Django by default looks for templates in the templates folder
    f = Path("running_job")
    if not f.is_file():
        logger.debug("No jobs are running currently (no job_running file)")
        return render(request, 'ymap_webtool/index.html')
    elif is_too_old(f, 15):
        logger.debug("The program is stuck (too old job_running file)")
        os.remove("running_job")
        clean_up(input_path, output_path, archive_path)
        return render(request, 'ymap_webtool/index.html')
    else:
        logger.debug("Another job is being processed. Come back later.")
        return render(request, 'ymap_webtool/server_is_busy.html')

def publications(request):
    #Specify the template to the Publications page
    return render(request, 'ymap_webtool/publications.html')
=======
	# specify the template to the Main page here. Django by default looks for templates in the templates folder
	return render(request, 'ymap_webtool/index.html')

def publications(request):
	# Specify the template to the Publications page
	return render(request, 'ymap_webtool/publications.html')
>>>>>>> 5316802b43aa5a59ce172376e2f492e745b65d89

def overview (request):
	# specify the template to the about/overview page here.
	return render(request, 'ymap_webtool/overview.html')

def pipeline (request):
	# specify the template to the about/pipeline page here.
	return render(request, 'ymap_webtool/pipeline.html')

def databases (request):
	# specify the template to the about/databases page here.
	return render(request, 'ymap_webtool/databases.html')

def input_format (request):
	# specify the template to the about/input-format page here.
	return render(request, 'ymap_webtool/input_format.html')

def results (request):
	# specify the template to the about/results page here.
	return render(request, 'ymap_webtool/results.html')

def how_to_cite(request):
	# specify the template to the about/how-to-cite page here
	return render(request, 'ymap_webtool/how_to_cite.html')

def troubleshooting(request):
	# specify the template to the troubleshooting page here
	return render(request, 'ymap_webtool/troubleshooting.html')

def contact(request):
	# specify the template to the contact page here
	# render the django form into the contact template
	form_class = ContactForm 
	
	return render(request, 'ymap_webtool/contact.html', {'form': form_class,})

def email_sent(request):
	# specify the template to the email_sent page here #to send email @test
	return render(request, 'ymap_webtool/contemail_sentact.html')

def finished(request):
<<<<<<< HEAD
	#Specify template to the page to display after submission
	logger.debug("time to give back the finished page with the script for downloading the results")
	is_running = False
	return render(request, 'ymap_webtool/finished.html')

def submission_fail(request):
	#Specify template to the page to display if submission fails
	#clean_up(input_path, output_path, archive_path)
	return render(request, 'ymap_webtool/submission_fail.html')

def processing(request):
	'''
	Checks whether the input file contains protein or gene level mutations, and runs a ymap function accordingly to map.
	zips the results
	deletes the original results, keeps only the zip file ?
	@return redirects to ymap_webtool/finished.
	'''
	logger.debug("started processing the job")
	if is_protein(os.path.join(input_path, input_file_name)) is True:
		logger.debug("it was a protein input file. running yproteins now")       
		run_yproteins()
	elif is_gene(os.path.join(input_path, input_file_name)) is True:
		logger.debug("it was a gene input file. running ygenes now")
		os.rename(os.path.join(input_path, input_file_name), os.path.join(input_path, gene_level_input_name))
		logger.debug("now the file " + input_file_name + " is renamed to " + gene_level_input_name)
		run_ygenes()
		logger.debug("finished processing the job")
	else:
		logger.debug("ERROR! Neither gene nor protein file!")
		return HttpResponseRedirect('submission_fail')
    
	make_archive(os.path.join(archive_path, archive_name), 'zip', root_dir = output_path)
	logger.debug("made an archive")
	# wipe original results folder
	# maybe keep the archive?
	return HttpResponseRedirect('finished')
=======
	# Specify template to the page to display after submission
	logger.debug("time to give back the finished page with the script for downloading the results")
	return render(request, 'ymap_webtool/finished.html')

def submission_fail(request):
	#Specify template to the page to display if submission fails
	#clean_up(input_path, output_path, archive_path)
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
		run_yproteins()
	elif is_gene(os.path.join(input_path, input_file_name)) is True:
		logger.debug("it was a gene input file. running ygenes now")
		os.rename(os.path.join(input_path, input_file_name), os.path.join(input_path, gene_level_input_name))
		logger.debug("now the file " + input_file_name + " is renamed to " + gene_level_input_name)
		run_ygenes()
		logger.debug("finished processing the job")
	else:
		logger.debug("ERROR! Neither gene nor protein file!")
		return HttpResponseRedirect('submission_fail')

	make_archive(os.path.join(archive_path, archive_name), 'zip', root_dir=output_path)
	logger.debug("made an archive")
	# wipe original results folder
	# maybe keep the archive?
	return HttpResponseRedirect('finished')
>>>>>>> 5316802b43aa5a59ce172376e2f492e745b65d89

def submitted(request):
<<<<<<< HEAD
	logger.debug("the job should be submitted and passed to processing")
	is_running = True
	return render(request, 'ymap_webtool/submitted.html')


def get_string(request):
	'''
	takes a string entered through the HTML form, 
	validates the form
	if the for is valid, saves it into a text file and redirects to the processing page
	else it assignes an empty object to the form and redirects to the submission_fail page
	@param request
	'''
	if is_running is False:
		is_running = True
		# if this is a POST request we need to process the form data
		if request.method == 'POST':
			# create a form instance and populate it with data from the request:
			form = StringForm(request.POST)
			# check whether it's valid:
			if form.is_valid():
				# process the data in form.cleaned_data as required
				# redirect to a new URL:
				save_string(form.cleaned_data.get('query'), input_path, input_file_name)
				return HttpResponseRedirect('submitted')
		# if the method was not POST we'll create a blank form
		else:
			form = StringForm()
		return HttpResponseRedirect('submission_fail')
	else:
		return HttpResponseRedirect('submission_fail')
 
=======
	logger.debug("the job should be submitted and passed to processing")
	return render(request, 'ymap_webtool/submitted.html')


>>>>>>> 5316802b43aa5a59ce172376e2f492e745b65d89

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
			# Process data here with a separate pydev class?
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
<<<<<<< HEAD
	'''
	takes a file that has been uploaded with the HTML form, 
	validates its extension (only .txt is allowed!)
	if the extension is fine, saves that file under the     input_path specified above, and redirects to the processing page
	else it redirects to the submission_fail page
	@param request
	'''
	f = Path("running_job")
	if not f.is_file():
		f = open("running_job", "w+")
		f.close()
		clean_up(input_path, output_path, archive_path)
		#handle file upload
		if request.method == 'POST':
		# create a form instance of the class UploadFileForm (cf form.py) and populate it with data from the request:
			form = UploadFileForm(request.POST, request.FILES)
		# check if the extention of the file is .txt (ps: this is already checked on the html)
			if validate_file_type(request.FILES['myfile'].name) is True:
				save_file(request.FILES['myfile'], input_path, input_file_name)
				logger.debug("get_file ran succcesfully")
				return HttpResponseRedirect('submitted')
		else:
			return HttpResponseRedirect('submission_fail')
	else:
		return HttpResponseRedirect('submission_fail')

def download_result(request):
    '''
    wraps the result zip file into an HttpResponse object and returns it
    '''
    #download the results to the user
    resp_file = open(os.path.join(archive_path, archive_name + ".zip"), 'rb') #the rb flag is needed on windows!motherwise r should be sufficient
    response = HttpResponse(resp_file, content_type='application/force-download')
    response['Content-Disposition'] = 'attachment; filename="%s"' % download_name
    logger.debug("download the zip file")
    f = Path("running_job")
    if f.is_file() is True:
        os.remove("running_job")
    return response    
=======
	clean_up(input_path, output_path, archive_path)
# handle file upload
	if request.method == 'POST':
	# create a form instance of the class UploadFileForm (cf form.py) and populate it with data from the request:
		form = UploadFileForm(request.POST, request.FILES)
	# check if the extention of the file is .txt (ps: this is already checked on the html)
		if validate_file_type(request.FILES['myfile'].name) is True:
			save_file(request.FILES['myfile'], input_path, input_file_name)
			logger.debug("get_file ran succcesfully")
			return HttpResponseRedirect('submitted')
	else:
		return HttpResponseRedirect('submission_fail')


def download_result(request):
	# download the results to the user
	resp_file = open(os.path.join(archive_path, archive_name + ".zip"), 'rb')  # the rb flag is needed on windows!motherwise r should be sufficient
	response = HttpResponse(resp_file, content_type='application/force-download')
	response['Content-Disposition'] = 'attachment; filename="%s"' % download_name
	logger.debug("download the zip file")
	return response
>>>>>>> 5316802b43aa5a59ce172376e2f492e745b65d89


def test_page(request):
	'''
	This view is for any functionalities that temporarily need to be tested.Shall be deleted upon completion of the app
	'''
	logger.debug("Print this to the command line please")
	response = HttpResponse('Nothing to test now')
	return response

#to send email @test
def send_email(request):
	logger.debug(request)
	
	if request.method == 'POST':
		form = ContactForm(data=request.POST)
		logger.debug(form)
		if form.is_valid():
			logger.debug("contact form: is valid = true")
			
			#store the information from the form into variable
			name = request.POST.get('name', '')
			from_email = request.POST.get('from_email', '')
			subject = request.POST.get('subject', '')
			message = request.POST.get('message', '')
			
			logger.debug("form input in variables")
			
			#use these variable to fill the email template
			template = get_template('ymap_webtool/_email_template.txt')
			
			context = {
				'name': name,
				'from_email': from_email,
				'subject': subject,
				'message': message,
			}
			
			content = template.render(context)
			logger.debug("content ready")
		
			
			#create an email from content (=filled template)
			new_email = EmailMessage(
				"New contact form submission",
				content,
				"yMap web" +'',
				['paquayleila@hotmail.com'],
				headers = {'Reply-To': from_email }
			)
			logger.debug("email ready")
			
			#send email
			new_email.send()

	else:
		form = ContactForm()		
	return render(request, "ymap_webtool/email_sent.html", {'form': form})
	

def success(request):

	return HttpResponse('Success! Thank you for your message.')


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
	# erase the input, output, but keep the archive under another name
	wipe_folder(input_path)
	wipe_folder(output_path)
	logger.debug("hello from cleanup")
	# os.rename(os.path.join(archive_path, archive_name + ".zip"), os.path.join(archive_path, "archived" + str(i) + ".zip"))

class AjaxRedirect(object):
	def process_response(self, request, response):
		if request.is_ajax():
			if type(response) == HttpResponseRedirect:
				response.status_code = 278
		return response

   
'''
This PyDev module contains functions, that process (save, call ymap upon it, 
zip the output) a String or a plain Txt input file recieved via HTML forms.
To Do:
    -    Check the file: is it the right input format?
    -     átnevezni függvényeket: save_string ---> upload_string, save_file ---> upload_file
    -    Implement in view: if its not the right input format, then remove the file 
         and inform the user that the format wasnt right
         
         AMI MŰKÖDIK:
             save_string
'''
import os.path
from .ymap import *




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
    
    
def download(wd, path, name, download_name = "results.zip"):
    '''
    Serves the user a file as an attachment using HttpResponse
    path: the file's path  and 
    wd: the working directory
    name: the name of the file to serve (with extension)
    download_name: default: results.zip. The name of the file on the client side
    '''
    resp_file = os.path.join(wd, path, name)
    response = HttpResponse(resp_file, content_type='application/force-download')
    response['Content-Disposition'] = 'attachment; filename="%s"' % download_name
    return response

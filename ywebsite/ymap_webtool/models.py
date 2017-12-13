from django.db import models
from celery import task
import os.path

@task()
def save_string(string, path, name):
    '''
    string: query text provided by user, to save as a text file
    path: where should the file be saved on the server side
    name: the file is going to be saved on the server side under this name (with extension)
    '''
    with open(os.path.join(path, name), 'w') as f:
        f.write(string)

@task()       
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



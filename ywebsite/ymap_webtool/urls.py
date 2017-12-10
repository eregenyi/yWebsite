from django.conf.urls import url
from . import views

# define URL patterns, link them to the corresponding views, and name those URLs
#!r' mean regex expression so avoid name that are subset of an other
#ex: avoi process and processing
urlpatterns = [
    url(r'index', views.index, name='index'),
    url(r'publications', views.publications, name='publications'),
    url(r'overview', views.overview , name='overview'),
<<<<<<< HEAD
    #url(r'process', views.process , name='process'),
=======
    url(r'pipeline', views.pipeline , name='pipeline'),
>>>>>>> 919ac5b84e3308fadf0d6ade22e99a1dadb070cf
    url(r'databases', views.databases , name='databases'),
    url(r'input_format', views.input_format , name='input_format'),
    url(r'results', views.results , name='results'),
    url(r'how_to_cite', views.how_to_cite, name='how_to_cite'),
    url(r'troubleshooting', views.troubleshooting, name='troubleshooting'),
    url(r'contact', views.contact, name='contact'),
    url(r'finished', views.finished, name='finished'),
    url(r'submission_fail', views.submission_fail, name='submission_fail'),
    url(r'get_string', views.get_string, name='get_string'),
    url(r'get_file', views.get_file, name='get_file'),
    url(r'test', views.test_page, name='test_page'),
    url(r'processing', views.processing, name='processing'),
    url(r'result', views.download_result, name='download_result'),
    url(r'clean_up', views.clean_up, name='clean_up'),
    url(r'submitted', views.submitted, name='submitted'),
]


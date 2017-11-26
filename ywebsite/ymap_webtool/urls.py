from django.conf.urls import url
from . import views

# define URL patterns, link them to the corresponding views, and name those URLs
urlpatterns = [
    url(r'index', views.index, name='index'),
    url(r'manual', views.manual, name='manual'),
    url(r'about', views.about, name='about'),
    url(r'publications', views.publications, name='publications'),
    url(r'post_submission', views.post_submission, name='post_submission'),
    url(r'submission_fail', views.submission_fail, name='submission_fail'),
    url(r'get_string', views.get_string, name='get_string'),
    url(r'get_file', views.get_file, name='get_file'),
    url(r'test', views.test_page, name='test_page'),
    url(r'processing', views.processing, name='processing'),
    url(r'result', views.download_result, name='download_result'),
    url(r'job_submitted', views.job_submitted, name='job_submitted'),
]
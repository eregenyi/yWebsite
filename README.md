# yWebsite
This is a repository dedicated to developing a web application for the yMap python package.
https://github.com/CSB-KUL/yMap


---------------------------------------------------------------------------------------------------------------------------------
7dec17: SECOND TRY 
issue with github, changes lost

did cahnge in mega folder not in the github repository
need to replace the one from github tomorrow with help

in index.html - form
I have no clue if you need me to add the weird stuff from your template
you should annotated so i don't spend time trying to figure out what u already figured out

what is the django message thing?
is it linking to an other template 
OR
diplsying a message with Django message framework
?

add location of result zip file in result_page

check in css "form commun":
input[type=sumbit]:hover have two background-color
does not seems to makke something not work but we never know

-----------------------------------------------------
FIRST TRY HTML TO DJANGO 
I added my html pages to the folder template/ymap_webtool
 
I updates my html name in view.py
	++ process change to pipeline
	! index is called home. I don't know if i have to change it for some technical reason
	
I updated the urls.py file 

In the template files I updates the url links
example:
	href="how-to-cite.html"
 was changed to 
	href={% url 'howtocite' %}
	
I added what i guess is the django dynamic part after the input form from home.html

			 <!--Django dynamic part-->
			{% if uploaded_file_url %}
			<p>File uploaded at: <a href="{{ uploaded_file_url }}">{{ uploaded_file_url }}</a></p>
			{% endif %}
			

I updated the css link in every html files

	<link rel="stylesheet" type="text/css" href="ymap.css">
  was changed to 
	{% load static %}
	<link rel="stylesheet" type="text/css" href="{% static 'polls/ymap.css' %}" />
	
I created a submission fail page
I created a finished page and added the script to it
I created a submitted page and added the script to it

QUESTION:
- should i had 	{% csrf_token %}?
- shoulnot we removethe komodo project from the static folder

from django import forms


class StringForm(forms.Form):
	'''
	This class is an abstraction of string input
	'''
	query = forms.CharField(label='Query')
	#label: specifies a human-friendly label
	#max_length: maximum nr of characters (the browser will prevent you to enter a longer text)
	#the Form class has an isValid() method, that checks whether the upper constraints are enforced
	#right now, the form saves values into a dictionary. W could leverage this to specify multiple input fields (e-mail, job name, etc)
	
class UploadFileForm(forms.Form):
	'''
	This class contains input files uploaded by the UserWarning
	'''
	title = forms.CharField(max_length=50)
	myfile = forms.FileField()
  

class ContactForm(forms.Form):
	'''
	This class is for the contact form @test
	'''
	name = forms.CharField(
		max_length = 30,
		widget=forms.TextInput(
			attrs={
				'placeholder': 'Write your name here'  
			}
		)				
	)
				   
	from_email = forms.EmailField(
		max_length= 50,
		widget=forms.EmailInput(
			attrs={
				'placeholder': 'Write your email here'  
			}
		)				
	)
	
	subject = forms.CharField(
		max_length = 50,
		widget=forms.TextInput(
			attrs={
				'placeholder': 'Write the subject of your message here'  
			}
		)									  					  
	)
	
	message = forms.CharField(
		max_length=2000,
		widget=forms.Textarea(
			attrs={
				'placeholder': 'Write your message here'  
			}	
		)
	)
	
def clean(self):
	'''
	This class is for the contact form @test
	'''
	cleaned_data = super(ContactForm, self).clean()
	name = cleaned_data.get('name')
	from_email = cleaned_data.get('from_email')
	subject = cleaned_data.get('subject')
	message = cleaned_data.get('message')
	if not name and not email and not subject and not message:
		raise forms.ValidationError('You have to write something!')
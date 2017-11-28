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
    

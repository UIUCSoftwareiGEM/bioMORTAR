from django import forms

class FeedbackForm(forms.Form):
  subject = forms.CharField(max_length=100)
  message = forms.CharField(widget=forms.Textarea())
  sender = forms.CharField(initial='Anonymous', required=False)
  sender_email = forms.EmailField(initial='anonymous@xyz.edu', required=False)

from django import forms
from django.core import validators
from django.core.exceptions import ValidationError

def base_validator(value):
    value = set(value.lower())
    req_set = {'a', 'u', 'g' , 'c', 't'}
    if bool(value.issubset(req_set)) == False:
        raise ValidationError("The entered sequence contains unknown bases!")

choices = [
        ('method_1', 'Takefuji'),
        ('method_2', 'Olke___'),
        ('method_3', 'Neural_'),
]
class AddSequence(forms.Form):
    seq = forms.CharField(min_length=5, validators=[base_validator], label="Enter sequence:")
    method = forms.CharField(label='Choose the algorithm:', 
                             widget=forms.RadioSelect(choices=choices))
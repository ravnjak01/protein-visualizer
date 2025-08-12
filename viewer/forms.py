from django import forms
from .models import ProteinFile

class ProteinFileForm(forms.ModelForm):
    class Meta:
        model = ProteinFile
        fields = ['name', 'file']

from django.db import models

# Create your models here.

class ProteinFile(models.Model):
    name = models.CharField(max_length=200)
    file = models.FileField(upload_to='pdb_files/')
    uploaded_at = models.DateTimeField(auto_now_add=True)
   
    def __str__(self):
        return self.name
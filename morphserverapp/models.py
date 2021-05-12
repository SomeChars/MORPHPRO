from django.db import models
import datetime
from django.core.files.storage import FileSystemStorage

my_storage = FileSystemStorage(location='./static/file_storage/')
# Create your models here.

class User(models.Model):

    name = models.CharField(max_length=256)

    email = models.EmailField(unique=True)

    password = models.CharField(max_length=256)

    user_created_at = models.DateTimeField()
    user_updated_at = models.DateTimeField()


    @classmethod
    def create_user(cls,data_dict):
        user = cls()
        user.name=data_dict['name']
        user.email=data_dict['email']
        user.password=data_dict['password']
        current_time = datetime.datetime.now()
        user.user_created_at=current_time
        user.user_updated_at=current_time
        return user

class MorphRequest(models.Model):

    author = models.IntegerField()

    protein_a_name = models.CharField(max_length=256)
    protein_b_name = models.CharField(max_length=256)

    created_at = models.DateTimeField()


    @classmethod
    def create_request(cls, data_dict):
        morph_request = cls()
        morph_request.author = data_dict['author']
        morph_request.protein_a_name = data_dict['protein_a_name']
        morph_request.protein_b_name = data_dict['protein_b_name']
        morph_request.created_at = datetime.datetime.now()
        return morph_request



class Pdb(models.Model):
    name = models.CharField(max_length=256)
    file = models.FileField(null=True,storage=my_storage)

    @classmethod
    def create_pdb(cls,name,file):
        pdb = cls()
        pdb.name = name
        pdb.file = file
        return pdb
# Generated by Django 3.2 on 2021-04-21 18:08

import django.core.files.storage
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('morphserverapp', '0011_alter_pdb_file'),
    ]

    operations = [
        migrations.AlterField(
            model_name='pdb',
            name='file',
            field=models.FileField(null=True, storage=django.core.files.storage.FileSystemStorage(location='./morphserverapp/static/file_storage'), upload_to=''),
        ),
    ]
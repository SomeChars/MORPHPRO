# Generated by Django 3.2 on 2021-04-21 15:48

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('morphserverapp', '0006_rename_morping_count_morphrequest_morphing_count'),
    ]

    operations = [
        migrations.RenameField(
            model_name='pdb',
            old_name='pdb_file',
            new_name='file',
        ),
        migrations.RenameField(
            model_name='pdb',
            old_name='pdb_name',
            new_name='name',
        ),
    ]
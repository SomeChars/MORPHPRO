# Generated by Django 3.2 on 2021-04-18 15:03

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('morphserverapp', '0005_auto_20210418_1500'),
    ]

    operations = [
        migrations.RenameField(
            model_name='morphrequest',
            old_name='morping_count',
            new_name='morphing_count',
        ),
    ]

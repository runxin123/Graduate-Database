# Generated by Django 2.2 on 2019-04-23 07:52

import calc.models
from django.db import migrations, models
import django.db.models.deletion
import djongo.models.fields


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='BatchCluster',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, unique=True)),
                ('name', models.CharField(max_length=50, null=True)),
                ('pretty_formula', models.CharField(max_length=50, null=True)),
                ('element_num', models.IntegerField(null=True)),
                ('volume', models.FloatField(null=True)),
                ('density', models.FloatField(null=True)),
                ('activity_coefficient', models.FloatField(null=True)),
                ('lattice_plane', models.TextField(null=True)),
                ('lattice_a', models.FloatField(null=True)),
                ('lattice_b', models.FloatField(null=True)),
                ('lattice_c', models.FloatField(null=True)),
            ],
        ),
        migrations.CreateModel(
            name='reactant',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=100, null=True)),
                ('element_types', models.TextField(null=True)),
                ('element_nums', models.IntegerField(null=True)),
                ('pretty_formula', models.CharField(max_length=100, null=True)),
                ('dos', models.TextField(null=True)),
                ('ldos', models.TextField(null=True)),
                ('energy', models.FloatField(null=True)),
                ('frequency', models.TextField(null=True)),
            ],
        ),
        migrations.CreateModel(
            name='surface',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, unique=True)),
                ('name', models.CharField(max_length=50, null=True)),
                ('pretty_formula', models.CharField(max_length=50, null=True)),
                ('element_num', models.IntegerField(null=True)),
                ('element_types', models.TextField(null=True)),
                ('volume', models.FloatField(null=True)),
                ('density', models.FloatField(null=True)),
                ('activity_coefficient', models.FloatField(null=True)),
                ('lattice_plane', models.CharField(max_length=50)),
                ('lattice_a', models.FloatField(null=True)),
                ('lattice_b', models.FloatField(null=True)),
                ('lattice_c', models.FloatField(null=True)),
                ('layer_thick', models.IntegerField(null=True)),
                ('size', models.TextField(null=True)),
                ('fix_layer', models.IntegerField(null=True)),
                ('dos', models.TextField()),
                ('ldos', models.TextField(null=True)),
                ('energy', models.FloatField(null=True)),
                ('zpe', models.FloatField(null=True)),
                ('workfuction', models.FloatField(null=True)),
                ('dcenter', models.TextField(null=True)),
                ('band_structure', models.TextField(null=True)),
                ('e_fermi', models.FloatField(null=True)),
            ],
        ),
        migrations.CreateModel(
            name='reaction',
            fields=[
                ('id', models.AutoField(primary_key=True, serialize=False, unique=True)),
                ('name', models.CharField(max_length=100, null=True)),
                ('element_types', models.TextField(null=True)),
                ('surface_t', djongo.models.fields.EmbeddedModelField(model_container=calc.models.surface, null=True)),
                ('reactant_t', djongo.models.fields.EmbeddedModelField(model_container=calc.models.reactant, null=True)),
                ('adsorption_energy', models.FloatField(null=True)),
                ('desorption_energy', models.FloatField(null=True)),
                ('binding_energy', models.FloatField(null=True)),
                ('catalyst_activity', models.FloatField(null=True)),
                ('selectivity', models.FloatField(null=True)),
                ('area_ratio', models.FloatField(null=True)),
                ('stability', models.FloatField(null=True)),
                ('coverage', models.FloatField(null=True)),
                ('dcenter', models.FloatField(null=True)),
                ('frequency', models.TextField(null=True)),
                ('band_structure', models.TextField(null=True)),
                ('adsorption_site', models.TextField(null=True)),
                ('dos', models.TextField(null=True)),
                ('layer_thick', models.IntegerField(null=True)),
                ('reactant_id', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='calc.reactant')),
                ('surface_id', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='calc.surface')),
            ],
        ),
    ]

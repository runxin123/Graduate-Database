from django.db import models
import django_mysql.models
import datetime
import djongo.models
# Create your models here.
# class surface(djongo.models.Model):
#     id=djongo.models.AutoField(primary_key=True,unique=True)
#     name =djongo.models.CharField(max_length=50,null=True)
#     pretty_formula = djongo.models.CharField(max_length=50,null=True)
#     element_num=djongo.models.IntegerField(null=True)
#     element_types=django_mysql.models.JSONField(null=True)
#     # VASP INPUT OUTPUT
#
#     # SLAB slab PROPERTY
#     #catalst_id = models.IntegerField(ForeginKey='reaction')
#     volume = djongo.models.FloatField(null=True)
#     density = djongo.models.FloatField(null=True)
#     activity_coefficient = djongo.models.FloatField(null=True)
#     lattice_plane  = djongo.models.CharField(max_length=50)
#
#     layer_thick = djongo.models.IntegerField(null=True)
#     size = django_mysql.models.JSONField(null=True) # {x:1,y:1,z:1}
#     fix_layer = djongo.models.IntegerField(null=True)
#
#     # calculate the Property
#     dos = djongo.models.TextField()
#     ldos = django_mysql.models.JSONField(null=True) # {'px':{/data},'py':{/data/},'pz':{/data/}}
#     energy = djongo.models.FloatField()
#     zpe = djongo.models.FloatField(null=True)
#     workfuction = djongo.models.FloatField(null=True)
#     dcenter = djongo.models.TextField(null=True)
#     band_structure = djongo.models.TextField(null=True)
#     e_fermi = djongo.models.FloatField(null=True)
#     def __str__(self):
#         return self.name
#     def test(self):
#         return self.element_types
#     class Meta:
#         app_label='calc'
#
# class reactant(django_mysql.models.Model):
#     id = models.IntegerField(primary_key=True,null=False,unique=True)
#     name = models.CharField(max_length=100,null=True)
#
#     # ELEMENT TYPES AND NUMS JSONFIELD
#
#     element_types = django_mysql.models.JSONField(null=True)#{'H':2,'O':1}
#     element_nums = models.IntegerField(null=True)#16
#     pretty_formula = models.CharField(max_length=100,null=True)#H2O2 == HO
#
#     # VASP INPUT OUTPUT
#
#     # calculate property
#     dos = models.TextField(null=True)
#     ldos = django_mysql.models.JSONField(null=True)
#     energy = models.FloatField(null=True)
#     frequency = django_mysql.models.JSONField(null=True) # {'H-O':2131}
#
#
#     # FUNCTION FOR REACTANTS
#     def __str__(self):
#         return self.name
#     def testJSON(self):
#         return self.element_types
#     class Meta:
#         app_label = "calc"
#
# class reaction(django_mysql.models.Model):
#     id = models.AutoField(primary_key=True,null=False,unique=True)
#     name = models.CharField(max_length=100,null=True)
#     element_types = django_mysql.models.JSONField(null=True)  # {types:['surface','bulk','cluster','particle'],'base':{'Pt':12,'Pd':12},'absorbate':{'H':2,'O':1}}
#     surface_id = models.ForeignKey(surface,on_delete=models.CASCADE)
#     reactant_id = models.ForeignKey(reactant,on_delete=models.CASCADE)
#     surface_t = djongo.models.EmbeddedModelField(
#         model_container=surface
#     )
#     reactant_t = djongo.models.EmbeddedModelField(
#         model_container=reactant
#     )
#
#     # calculate property
#     adsorption_energy = models.FloatField(null=True)
#     desorption_energy = models.FloatField(null=True)
#     binding_energy =  models.FloatField(null=True)
#     catalyst_activity =  models.FloatField(null=True)
#     selectivity = models.FloatField(null=True)
#     area_ratio = models.FloatField(null=True)
#     stability = models.FloatField(null=True)
#     coverage = models.FloatField(null=True)
#     dcenter = models.FloatField(null=True)
#     #JSON FOR BOND_LENGTH AND FREQUENCY
#     frequency = django_mysql.models.JSONField(null=True) # {'H-O':2131}
#     band_structure = django_mysql.models.JSONField(null=True) # {'url':'/*/','data':{/data/}}
#     adsorption_site = django_mysql.models.JSONField(null=True) #{"atop":[1,2,3],"bridge":[1,2,3,4],"hellow":[1,2,3,4,5,6]}
#     #TEXT FEILD REACTION PROPERTY
#     dos = models.TextField(null=True)
#     layer_thick = models.IntegerField(null=True)
#     class Meta:
#         app_label= 'calc'

class vaspinout(django_mysql.models.Model):
    id = models.AutoField(unique=True,primary_key=True)
    catagory = models.IntegerField(null=True) ##1 = surface; 2 = reactant; 3 = reaction; 4 = cluster
    name  = models.CharField(max_length=100,null=True)
    compare_id = models.IntegerField(null=True) # {1..500} [1..16]
    #VASP INPUT & OUTPUT
    # input
    incar = models.TextField(null=True)
    kpoints = models.TextField(null=True)
    poscar = models.TextField(null=True)
    # output
    contcar = models.TextField(null=True)
    enginval = models.TextField(null=True)
    oszicar = models.TextField(null=True)
    doscar = models.TextField(null=True)
    vasprun = models.TextField(null=True)
    xdatcar = models.TextField(null=True)
    wavecar = models.TextField(null=True)
    pub_date = models.DateField(auto_now_add=False,default=datetime.date.today)
    mod_date = models.DateField(auto_now_add=False,default=datetime.date.today)
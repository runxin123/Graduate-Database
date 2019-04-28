from djongo import models
import djongo.models.json
# Create your models here.
import django_mysql.models
import djongo.models.fields
import datetime
from mongoengine import *
connect('catalystpro',username='test',password='test')
class BatchCluster(models.Model):
    id=models.AutoField(primary_key=True,unique=True)
    name =models.CharField(max_length=50,null=True)
    pretty_formula = models.CharField(max_length=50,null=True)
    element_num=models.IntegerField(null=True)
    volume = models.FloatField(null=True)
    density = models.FloatField(null=True)
    activity_coefficient = models.FloatField(null=True)
    lattice_plane  = models.TextField(null=True) #{'111':energy,'110':energy}
    lattice_a = models.FloatField(null=True)
    lattice_b = models.FloatField(null=True)
    lattice_c = models.FloatField(null=True)

    class Meta:
        app_label='calc'

class surface(models.Model):
    id=models.AutoField(primary_key=True,unique=True)
    name =models.CharField(max_length=50,null=True)
    pretty_formula = models.CharField(max_length=50,null=True)
    element_num=models.IntegerField(null=True)

    # JSONFIELD
    # element_types = djongo.models.json.JSONFieldBase(null=True)

    element_types = models.TextField(null=True)
    # VASP INPUT OUTPUT

    # SLAB slab PROPERTY
    # catalst_id = models.IntegerField(ForeginKey='reaction')
    volume = models.FloatField(null=True)
    density = models.FloatField(null=True)
    activity_coefficient = models.FloatField(null=True)
    lattice_plane  = models.CharField(max_length=50)
    lattice_a = models.FloatField(null=True)
    lattice_b = models.FloatField(null=True)
    lattice_c = models.FloatField(null=True)
    layer_thick = models.IntegerField(null=True)

    # size = djongo.models.json.JSONFieldBase({'x':1,'y':1,'z':1})  # {x:1,y:1,z:1}
    size = models.TextField(null=True)  # {x:1,y:1,z:1}
    # JSONFIELD


    fix_layer = models.IntegerField(null=True)

    # calculate the Property
    dos = models.TextField()
    ldos = models.TextField(null=True) # {'px':{/data},'py':{/data/},'pz':{/data/}}
    energy = models.FloatField(null=True)
    zpe = models.FloatField(null=True)
    workfuction = models.FloatField(null=True)
    dcenter = models.TextField(null=True)
    band_structure = models.TextField(null=True)
    e_fermi = models.FloatField(null=True)
    objects = models.DjongoManager()
    def __str__(self):
        return self.name
    def test(self):
        return self.element_types
    class Meta:
        app_label='calc'

class reactant(models.Model):
    # id = models.IntegerField(primary_key=True,null=False,unique=True)
    name = models.CharField(max_length=100,null=True)

    # ELEMENT TYPES AND NUMS JSONFIELD
    # JSONFIELD
    # element_types = djongo.models.json.JSONFieldBase(null=True)#{'H':2,'O':1}
    element_types = models.TextField(null=True)#{'H':2,'O':1}
    element_nums = models.IntegerField(null=True)#16
    pretty_formula = models.CharField(max_length=100,null=True)#H2O2 == HO

    # VASP INPUT OUTPUT

    # calculate property
    dos = models.TextField(null=True)


    # JSONFIELD
    ldos = models.TextField(null=True)
    energy = models.FloatField(null=True)

    #JSONFIELD
    frequency = models.TextField(null=True) # {'H-O':2131}



    # FUNCTION FOR REACTANTS
    def __str__(self):
        return self.name
    def testJSON(self):
        return self.element_types
    class Meta:
        app_label = "calc"

class reaction(models.Model):
    id = models.AutoField(primary_key=True,null=False,unique=True)
    name = models.CharField(max_length=100,null=True)

    #JSONFIELD
    # element_types = djongo.models.json.JSONFieldBase(null=True)  # {types:['surface','bulk','cluster','particle'],'base':{'Pt':12,'Pd':12},'absorbate':{'H':2,'O':1}}
    element_types = models.TextField(null=True)  # {types:['surface','bulk','cluster','particle'],'base':{'Pt':12,'Pd':12},'absorbate':{'H':2,'O':1}}
    surface_id = models.ForeignKey(surface,on_delete=models.CASCADE)
    reactant_id = models.ForeignKey(reactant,on_delete=models.CASCADE)
    surface_t = models.EmbeddedModelField(
        model_container=surface
    )
    reactant_t = models.EmbeddedModelField(
        model_container=reactant
    )

    # calculate property
    adsorption_energy = models.FloatField(null=True)
    desorption_energy = models.FloatField(null=True)
    binding_energy =  models.FloatField(null=True)
    catalyst_activity =  models.FloatField(null=True)
    selectivity = models.FloatField(null=True)
    area_ratio = models.FloatField(null=True)
    stability = models.FloatField(null=True)
    coverage = models.FloatField(null=True)
    dcenter = models.FloatField(null=True)
    #JSON FOR BOND_LENGTH AND FREQUENCY


    #JSONFIELD
    # frequency = djongo.models.json.JSONField({'H-O':2111}) # {'H-O':2131}
    frequency = models.TextField(null=True)
    band_structure = models.TextField(null=True) # {'url':'/*/','data':{/data/}}
    adsorption_site = models.TextField(null=True) #{"atop":[1,2,3],"bridge":[1,2,3,4],"hellow":[1,2,3,4,5,6]}


    #TEXT FEILD REACTION PROPERTY
    dos = models.TextField(null=True)
    layer_thick = models.IntegerField(null=True)
    class Meta:
        app_label='calc'

class someItem(Document):
    title = StringField(max_length=100,null=True)
    data_modified = DateTimeField(default=datetime.datetime.utcnow())
    bob = ListField()
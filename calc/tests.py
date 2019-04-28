
# Create your tests here.
from django.test import TestCase
from calc.models import surface, reactant, reaction
from django.db import models

class surfaceTestCase(TestCase):
    multi_db = True
    def setUp(self):
        reactant_c = reactant.objects.create(name="CO", pretty_formula="CO")
        # reaction_c = reaction.objects.create(name='Pt-111-CO')
        surface_c = surface.objects.create(name="Pt-111", pretty_formula="Pt20")
    def tearDown(self):
        surface.objects.filter(name='Pt-110').delete()
        auths = surface.objects.all().values()
        print(auths)

    def test_surface_test(self):
        auths = surface.objects.all().values()
        print(surface.energy)

    class Meta:
        app_label = 'calc'
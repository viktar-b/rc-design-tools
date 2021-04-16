import unittest
from CrackWidthCalculator import *

class TestSteelProperties(unittest.TestCase):

    f_yk = 500 

    def setUp(self):
        self.steel_properties = SteelProperties(self.f_yk)
        
    def test_properties(self):
        self.assertEqual(self.steel_properties.f_yk, self.f_yk)

        self.assertEqual(self.steel_properties.f_yd, self.f_yk / 1.15)

class TestConcreteProperties(unittest.TestCase):
    pass


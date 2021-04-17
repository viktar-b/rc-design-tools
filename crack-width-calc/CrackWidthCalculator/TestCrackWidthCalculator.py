import unittest
from unittest.mock import patch

from CrackWidthCalculator import *

class TestSteelProperties(unittest.TestCase):

    f_yk = 500 

    def setUp(self):
        self.steel_properties = SteelProperties(self.f_yk)
        
    def test_properties(self):
        self.assertEqual(self.steel_properties.f_yk, self.f_yk)
        self.assertEqual(self.steel_properties.f_yd, self.f_yk / 1.15)

    @patch('builtins.print')
    def test_print(self, mock_print):
        print(repr(self.steel_properties))
        mock_print.assert_called_with(
            f"E_s = 200000\n"\
            "gamma_s = 1.15\n" + \
            "f_yk = 500\n" + \
            "f_yd = 434.78\n"
        )

class TestConcreteProperties(unittest.TestCase):
    f_ck = 35

    def setUp(self):
        self.concrete_properties = ConcreteProperties(self.f_ck)

    def test_properties(self):
        self.assertEqual(self.concrete_properties.f_ck, self.f_ck)
        self.assertEqual(self.concrete_properties.f_cm, self.f_ck + 8)


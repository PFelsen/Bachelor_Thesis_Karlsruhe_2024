import unittest
import numpy as np
import shutil
import sys
from os.path import abspath, join
from os.path import dirname, realpath
from os.path import isfile
sys.path.append("..")
from utilities import vtu_utilities


class TestVtu(unittest.TestCase):
    @classmethod
    def setUp(cls):
        cls.sample_path = abspath(join(dirname(realpath(__file__)), "test_vtu_files"))

  
class TestPlotting(TestVtu):

    # helper functions for the tests:

    def try_add_pcolormesh(self, file_name, parallel=True):
        s = vtu_utilities.VtuPlot()
        s.set_wd(abspath(join(self.sample_path, file_name)))

        if parallel:
            s.add_pcolormesh('u.pvtu') 
        else:
            s.add_pcolormesh('u.vtu') 
      
        s.save("output.png")
        s.close()
        self.assertEqual(isfile("output.png"), True)

    def try_add_vtu(self, file_name, parallel=True):
        s = vtu_utilities.VtuPlot()
        s.set_wd(abspath(join(self.sample_path, file_name)))
        
        if parallel:
            s.add_vtu('u.pvtu') 
        else:
            s.add_vtu('u.vtu') 
       
        s.save("output.png")
        s.close()
        self.assertEqual(isfile("output.png"), True)

    def try_add_mesh(self, file_name, parallel=True):
        s = vtu_utilities.VtuPlot()
        s.set_wd(abspath(join(self.sample_path, file_name)))
        
        if parallel:
            s.add_mesh('u.pvtu') 
        else:
            s.add_mesh('u.vtu') 

        s.save("output.png")
        s.close()
        self.assertEqual(isfile("output.png"), True)

    def try_non_default_plotting(self, file_name, parallel=True):
        s = vtu_utilities.VtuPlot()
        s.set_wd(abspath(join(self.sample_path, file_name)))

        if parallel:
            s.add_imshow('u.pvtu') 
        else:
            s.add_imshow('u.vtu') 

        s.save("output.png")
        s.close()
        self.assertEqual(isfile("output.png"), True)
     
    def try_all(self, file_name, parallel=True):
        self.try_add_pcolormesh(file_name, parallel)
        self.try_add_vtu(file_name, parallel)
        self.try_add_mesh(file_name, parallel)
        self.try_non_default_plotting(file_name, parallel)

    # test functions

    def test_genplot(self):
        s = vtu_utilities.VtuPlot()
        s.set_wd(self.sample_path)
        s.save("output.png")
        self.assertEqual(isfile("output.png"), True)

    def test_parallel_ascii(self):
        self.try_all("p_ascii")
        
    def test_nonparallel_ascii(self):
        self.try_all("np_ascii", parallel=False)
        
    def test_parallel_raw(self):
        self.try_all("p_raw")

    def test_nonparallel_raw(self):
        self.try_all("np_raw", parallel=False)

    def test_parallel_base64(self):
        self.try_all("p_base64")

    # still produces error: TypeError: 'NoneType' object is not subscriptable in line 245 of vtu_utilities.py
    # def test_nonparallel_base64(self):
        # self.try_all("np_base64", parallel=False)


if __name__ == '__main__':
    unittest.main()

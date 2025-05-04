import os
import sys, subprocess
from setuptools import setup, find_packages,Extension
from setuptools.command.build_ext import build_ext

# Modified from pyscf setup.py
class FortranBuild(build_ext):
    def run (self):
        extension = self.extensions[0]
        assert extension.name == "chemPackage_fortran_library"
        self.build_make(extension)
    
    def build_make(self, extension):
        self.announce("Building mfunc", level = 3)
        src_dir = os.path.abspath(os.path.join(__file__, '..', 'src', 'mfunc'))
        cmd = ['make', '-C', src_dir]
        # if self.dry_run:
        #     self.announce(' '.join(cmd))
        # else:
        self.spawn(cmd)

        self.announce("Building f2py for chemPackage", level = 3)
        src_dir = os.path.abspath(os.path.join(__file__, '..', 'src', 'chemPackage', 'f2py'))
        cmd = ['make', '-C', src_dir]
        self.spawn(cmd)

setup(
    ext_modules = [Extension('chemPackage_fortran_library',[])],
    cmdclass={'build_ext': FortranBuild },
)

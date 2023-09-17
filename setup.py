import os
import sys, subprocess
from setuptools import setup, find_packages,Extension
from setuptools.command.build_ext import build_ext

# Modified from pyscf setup.py
class FortranBuild(build_ext):
    def run (self):
        extension = self.extensions[0]
        assert extension.name == "chemPackage_fortran_placeholder"
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

    # To remove the infix string like cpython-37m-x86_64-linux-gnu.so
    # Python ABI updates since 3.5
    # https://www.python.org/dev/peps/pep-3149/
    def get_ext_filename(self, ext_name):
        ext_path = ext_name.split('.')
        filename = build_ext.get_ext_filename(self, ext_name)
        name, ext_suffix = os.path.splitext(filename)
        return os.path.join(*ext_path) + ext_suffix


setup(
    name='chemPackage',
    version='1.1.0',
    description='Python package for interacting with data from quantum chemistry programs',
    url='https://github.com/jensengrouppsu/chemPackage',
    author='Jensen Research Group',
    author_email='jeff.becca@gmail.com',
    license='GPL v3.0',
    packages=find_packages('src'),
    package_dir={'':'src'},
    include_package_data=True,
    install_requires=['numpy',
                      'natsort',
                      ],
    ext_modules = [Extension('chemPackage_fortran_placeholder',[])],
    cmdclass={'build_ext': FortranBuild },
    classifiers=['Development Status :: 1 - Planning',
                 'Intended Audience :: Science/Research',
                 'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
                 'Natural Language :: English',
                 'Programming Language :: Python :: 3,',
                 'Topic :: Scientific/Engineering :: Chemistry',
                 'Topic :: Scientific/Engineering :: Physics',],
)

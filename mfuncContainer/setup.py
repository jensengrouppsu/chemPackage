from setuptools import setup

setup(
    name='mfunc',
    version='1.1.0',
    description='Python package for interacting with data from quantum chemistry programs',
    url='https://github.com/jensengrouppsu/chemPackage',
    author='Jensen Research Group',
    author_email='jeff.becca@gmail.com',
    license='GPL v3.0',
    packages=['mfunc'],

    install_requires=['numpy',
                      'scipy'],
    classifiers=['Development Status :: 1 - Planning',
                 'Intended Audience :: Science/Research',
                 'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
                 'Natural Language :: English',
                 'Programming Language :: Python :: 3,',
                 'Topic :: Scientific/Engineering :: Chemistry',
                 'Topic :: Scientific/Engineering :: Physics',],
)

from __future__ import print_function, division
from ..chemdata import ChemData
from ..constants import elem
from numpy import array, where, empty_like, reshape
from sys import exit

class TDSPEC(ChemData):
    
    def __init__(self, name):
        import os

        ChemData.__init__(self)
        self.program = 'TDSPEC'
        # Find extention
        ftype = os.path.splitext(name)[1]
        # Zhongwei: 'hpol' for hyperpolarizability
        #           'sfg' for sfg hyperpolarizability
        if ftype not in ('.out','.pol','.hpol','.sfg'):
            raise ValueError (ftype+' not a recognized extension')
        self.filetype = ftype[1:]
        self.filename = name

    def _collect(self, abort=False):

        # Set the abort flag
        self._abort = abort

        # Read in file
        from read_file import read_file
        f, indices = read_file(self)

        # Read input block
        from input_block import collect_input
        collect_input(self, f, indices)

        # Determine calculation type
        self.__det_calc_type()

        # Return now if an input file (not implemented yet)
#        if self.filetype != 'out': return

        # Techical properties
        #from tdspecproperties import collect_technical
        #collect_technical(self, f, indices)

        # Collect all exciting information
        from tdspecproperties import collect_tdspec
        collect_tdspec(self, f, indices)

        # Collect timing
        #from tdspecproperties import collect_timing
        #collect_timing(self, f, indices)

        # Collect printed normal modes
        if 'NORMAL MODE' in indices:
            from polarizability import collect_frequency
            collect_frequency(self, f, indices)
            self.calctype.add('FREQUENCIES')

        # Collect Excitation Energy
        if 'EXCITATION ENERGY' in indices:
            from tdspecproperties import collect_e_frequency
            collect_e_frequency(self, f, indices)

        # Collect polarizability (alpha) tensor
        if 'POLARIZABILITY' in indices:
            from polarizability import collect_alpha
            collect_alpha(self, f, indices)
            self.calctype.add('POLARIZABILITY')
            self.calctype.add('RAMAN')
            self.subkey.add('LIFETIME')

        # Collect A-tensor
        if 'A-TENSOR' in indices:
            from polarizability import collect_Atensor
            collect_Atensor(self, f, indices)
            self.subkey.add('A-TENSOR')

        # Collect As-tensor
        if 'As-TENSOR' in indices:
            from polarizability import collect_Astensor
            collect_Astensor(self, f, indices)
            self.subkey.add('As-TENSOR')

        # Collect G-tensor
        if 'G-TENSOR' in indices:
            from polarizability import collect_Gtensor
            collect_Gtensor(self, f, indices)
            self.subkey.add('G-TENSOR')
            if 'A-TENSOR' in self.subkey:
                self.calctype.add('ROA')

        # Collect Gs-tensor
        if 'Gs-TENSOR' in indices:
            from polarizability import collect_Gstensor
            collect_Gstensor(self, f, indices)
            self.subkey.add('Gs-TENSOR')

        # Collect C-tensor
        if 'C-TENSOR' in indices:
            from polarizability import collect_Ctensor
            collect_Ctensor(self, f, indices)
            self.subkey.add('C-TENSOR')

        # Collect D-tensor
        if 'D-TENSOR' in indices:
            from polarizability import collect_Dtensor
            collect_Dtensor(self, f, indices)
            self.subkey.add('D-TENSOR')

        # Collect D-tensor
        if 'Ds-TENSOR' in indices:
            from polarizability import collect_Dstensor
            collect_Dstensor(self, f, indices)
            self.subkey.add('Ds-TENSOR')

        # Zhongwei
        # Collect hyperpolarizability (beta) tensor
        if 'HYPERPOLARIZABILITY' in indices:
            from polarizability import collect_beta
            collect_beta(self, f, indices)
            self.calctype.add('HYPERPOLARIZABILITY')
            self.calctype.add('HYPERRAMAN')
            self.subkey.add('LIFETIME')


    def __det_calc_type(self):
        '''Determine the calculation type based on what is found in the keys'''
        self.calctype.add('TDSPEC')
        if 'OPA' in self.key:
            self.calctype.add('OPA')
        elif 'RRS' in self.key:
            self.calctype.add('RRS')
        elif 'TPA' in self.key:
            self.calctype.add('TPA')
        elif 'RHRS' in self.key:
            self.calctype.add('RHRS')

from __future__ import print_function
from numpy import array, zeros, arange, append

class fde_class():
    '''Class to collect and store all of the FDE-related inputs
    and outputs.'''

    def __init__(self, f, indices):

        # Electrostatic energy
        key = 'FDE ELECTROSTATIC ENERGY'
        if key in indices:
            for s in indices[key]:
                try:
                    self.electrostatic_energy.append(float(f[s].split()[3]))
                except AttributeError:
                    self.electrostatic_energy = [float(f[s].split()[3])]
            self.electrostatic_energy = array(self.electrostatic_energy)

        # Error
        key = 'FDE ERROR'
        if key in indices:
            li = -1
            for s in indices[key]:
                if s > li:
                    self.error = float(f[s].split()[4])
                    li = s

        # Interaction Energy
        key = 'FDE INTERACTION ENERGY'
        if key in indices:
            for s in indices[key]:
                try:
                    self.interaction_energy.append(float(f[s].split()[4]))
                except AttributeError:
                    self.interaction_energy = [float(f[s].split()[4])]
            self.interaction_energy = array(self.interaction_energy)

        # Subsystem A
        key = 'FDE SUBSYSTEM A'
        if key in indices:
            for s in indices[key]:
                try:
                    self.subsys_A.append(float(f[s].split()[5]))
                except AttributeError:
                    self.subsys_A = [float(f[s].split()[5])]
            self.subsys_A = array(self.subsys_A)
            if (self.subsys_A == 0.).all(): del(self.subsys_A)

        # Subsystem B
        key = 'FDE SUBSYSTEM B'
        if key in indices:
            for s in indices[key]:
                try:
                    self.subsys_B.append(float(f[s].split()[5]))
                except AttributeError:
                    self.subsys_B = [float(f[s].split()[5])]
            self.subsys_B = array(self.subsys_B)
            if (self.subsys_B == 0.).all(): del(self.subsys_B)

        # All subsystems
        key = 'FDE ALL SUBSYSTEMS'
        if key in indices:
            for s in indices[key]:
                try:
                    self.subsys_AB.append(float(f[s].split()[5]))
                except AttributeError:
                    self.subsys_AB = [float(f[s].split()[5])]
            self.subsys_AB = array(self.subsys_AB)

        # Embedding energy
        key = 'FDE EMBEDDING ENERGY'
        if key in indices:
            for s in indices[key]:
                try:
                    self.embedding_energy.append(float(f[s].split()[4]))
                except AttributeError:
                    self.embedding_energy = [float(f[s].split()[4])]
            self.embedding_energy = array(self.embedding_energy)

        # Potentials
        if 'FDE VT' in indices:
            for s in indices['FDE VT']:
                self.kinetic_potential = float(f[s].split()[4])
        if 'FDE VT LDA' in indices:
            for s in indices['FDE VT LDA']:
                self.kinetic_potential_lda = float(f[s].split()[4])
        if 'FDE ENUC' in indices:
            for s in indices['FDE ENUC']:
                self.nuclear_energy = float(f[s].split()[4])
        if 'FDE J[A,B]' in indices:
            for s in indices['FDE J[A,B]']:
                self.coulomb_energy = float(f[s].split()[5])

        # Kinetic energies
        if 'FDE T A' in indices:
            for s in indices['FDE T A']:
                self.kinetic_energy_A = float(f[s].split()[4])
        if 'FDE T B' in indices:
            for s in indices['FDE T B']:
                self.kinetic_energy_B = float(f[s].split()[5])
        if 'FDE T AB' in indices:
            for s in indices['FDE T AB']:
                self.kinetic_energy_AB = float(f[s].split()[5])
        if 'FDE T AB - B' in indices:
            for s in indices['FDE T AB - B']:
                self.kinetic_energy_ABmB = float(f[s].split()[4])
        if 'FDE T AB - A - B' in indices:
            for s in indices['FDE T AB - A - B']:
                self.kinetic_energy_ABmAmB = float(f[s].split()[4])

        # Kinetic energies (LDA)
        if 'FDE T A LDA' in indices:
            for s in indices['FDE T A LDA']:
                self.kinetic_energy_A_lda = float(f[s].split()[4])
        if 'FDE T B LDA' in indices:
            for s in indices['FDE T B LDA']:
                self.kinetic_energy_B_lda = float(f[s].split()[5])
        if 'FDE T AB LDA' in indices:
            for s in indices['FDE T AB LDA']:
                self.kinetic_energy_AB_lda = float(f[s].split()[4])
        if 'FDE T AB - B LDA' in indices:
            for s in indices['FDE T AB - B LDA']:
                self.kinetic_energy_ABmB_lda = float(f[s].split()[4])
        if 'FDE T AB - A - B LDA' in indices:
            for s in indices['FDE T AB - A - B LDA']:
                self.kinetic_energy_ABmAmB_lda = float(f[s].split()[4])

        # XC energies
        if 'FDE XC A' in indices:
            for s in indices['FDE XC A']:
                self.xc_energy_A = float(f[s].split()[4])
        if 'FDE XC B' in indices:
            for s in indices['FDE XC B']:
                self.xc_energy_B = float(f[s].split()[5])
        if 'FDE XC AB' in indices:
            for s in indices['FDE XC AB']:
                self.xc_energy_AB = float(f[s].split()[5])
        if 'FDE XC AB - B' in indices:
            for s in indices['FDE XC AB - B']:
                self.xc_energy_ABmB = float(f[s].split()[4])
        if 'FDE XC AB - A - B' in indices:
            for s in indices['FDE XC AB - A - B']:
                self.xc_energy_ABmAmB = float(f[s].split()[4])

        # XC energies (LDA)
        if 'FDE XC A LDA' in indices:
            for s in indices['FDE XC A LDA']:
                self.xc_energy_A_lda = float(f[s].split()[4])
        if 'FDE XC B LDA' in indices:
            for s in indices['FDE XC B LDA']:
                self.xc_energy_B_lda = float(f[s].split()[5])
        if 'FDE XC AB LDA' in indices:
            for s in indices['FDE XC AB LDA']:
                self.xc_energy_AB_lda = float(f[s].split()[4])
        if 'FDE XC AB - B LDA' in indices:
            for s in indices['FDE XC AB - B LDA']:
                self.xc_energy_ABmB_lda = float(f[s].split()[4])
        if 'FDE XC AB - A - B LDA' in indices:
            for s in indices['FDE XC AB - A - B LDA']:
                self.xc_energy_ABmAmB_lda = float(f[s].split()[3])

        # FDEc uncoupled states
        if 'FDEU STATES' in indices:
            self.fdeu_energies = []
            self.fdeu_transdip = []
            for s in indices['FDEU STATES']:
                for i in range(s,s+1000):
                    if len(f[i].split())==0: break
                    if len(f[i].split())==8:
                        self.fdeu_energies.append(float(f[i].split()[4]))
                        self.fdeu_transdip.append(f[i].split()[5:8])
                self.fdeu_energies = array(self.fdeu_energies, dtype=float) * 0.0367493088 # convert to hartrees from eV
                self.fdeu_transdip = array(self.fdeu_transdip, dtype=float)
                self.fdec_nexci = len(self.fdeu_energies)

        if 'FDEC ENERGIES' in indices:
            for s in indices['FDEC ENERGIES']:
                e = s + self.fdec_nexci
                self.fdec_energies = array([f[i].split()[1] for i in xrange(s,e)], dtype=float)
                self.fdec_energies *= 0.0367493088 # convert from eV to hartrees
                self.fdec_oscillator_strength = array([f[i].split()[2] for i in xrange(s,e)], dtype=float)

        if 'FDEC DIPOLE MOMENTS' in indices:
            for s in indices['FDEC DIPOLE MOMENTS']:
                e = s + self.fdec_nexci
                self.fdec_transdip = array([f[i].split()[2:5] for i in xrange(s,e)], dtype=float)

        # FDEc polarizability
        if 'FDEC POLARIZABILITY' in indices:
            s = indices['FDEC POLARIZABILITY']
            try:
                self.fdec_pol = array([f[i].split()[1:4] for i in xrange(s,s+3)], dtype=float)
            except ValueError:
                self.fdec_pol = zeros((3,3))
            except IndexError:
                self.fdec_pol = zeros((3,3))

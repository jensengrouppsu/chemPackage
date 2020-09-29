'''

Class Attributes
^^^^^^^^^^^^^^^^^

.. py:attribute:: program

    For ADF, this will be the string "ADF", for NWChem, this will be
    "NWChem", etc.

.. py:attribute:: filetype

    This is 'out' for an output file and 'in' for an input file.

.. py:attribute:: filename

  The name of the file the data was collected from.

.. py:attribute:: host

  The name of the computer the calculation ran on.

.. py:attribute:: nprocs

  The number of processors used for the calculation as an integer.

.. py:attribute:: start

  A datetime object containing the date and time of day on which the
  calculation began.

.. py:attribute:: real_time

  A timedelta object containing the real time, that it took for the
  calculation to complete.

.. py:attribute:: cpu_time

  A timedelta object containing the cpu time, that it took for the
  calculation to complete.

.. py:attribute:: routine_times

   A dict containing timing information for the individual routines.
   The keys are the routine names, and each routine is a three element
   tuple containing:
   
   .. hlist::
    :columns: 3
    
    - the real time
    - the system time
    - the number of calls

.. py:attribute:: termination

  A string contatining the temination status of the calculation.  It
  returns 'NORMAL TERMINATION' if the calculation was successful,
  'Error: message', where message describes an error that occured,
  or 'Unknown Calculation Failure' if there was some unknown failure.

.. py:attribute:: title

  The title of the calculation, if any appeared.   

.. py:attribute:: calctype

  The type of calculation performed.  It is a set all valid types.

.. py:attribute:: key

  This is a dictionary of all the keywords and their associated
  parameters found in the input block, with the keywords as the
  keys and the keyword parameters in tuple form as the values.  If a
  keyword was a pure keyword with no parameters, then the value is
  True.

.. py:attribute:: subkey

  This is a set of important subkeywords that might affect how
  data is interpreted and would not be seen from the keywords only.
 
.. py:attribute:: mixedkeys

  This is a tuple of all the mixed-type keywords known to this
  class.  This means that the keyword is followed by either a block of
  parameters terminating with the termination keyword or a line only, 
  or both.

.. py:attribute:: blockkeys

  This is a tuple of all the block-type keywords known to this
  class.  This means that the keyword is followed by a block of
  parameters terminating with the termination keyword.

.. py:attribute:: linekeys

  This is a tuple of all the line-type keywords known to this
  class.  This means that the parameters are found on the same line as
  the keyword.

.. py:attribute:: singlekeys

  This is a tuple of all the single-type keywords known to this
  class.  This means that the keyword stands alone, with no parameters.

  .. note:: 
   For ADF, the keyword **EFIELD** is purposely missing from the three
   above tuples because it can be block- or line-type.

.. py:attribute:: natoms

  This is the number of QM atoms in the calculation.

.. py:attribute:: atoms

  This is a numpy vector of the QM atoms in the calculation, given
  as 'Xx', where  Xx is the element.  It has length of :py:attr:`natoms`.

.. py:attribute:: elements

  This is a set of the QM elements in the calculation.  There are
  no duplicates.

.. py:attribute:: nelements

  This is the number of QM elements in the calculation.

.. py:attribute:: coordinates

  This is a numpy array of the coordinates of the atoms.  It has
  dimensions of :py:attr:`natoms` x 3.

.. py:attribute:: initial_geometry

  For a geometry optimization or similar calculation, this stores
  the geometry given on input.  If the coordinates are not changed
  during the calculation, this is identical to :py:attr:`coordinates`.
  Either way, it has the same dimensions of :py:attr:`coordinates`.

.. py:attribute:: symmetry

  The symmetry group of the molecule.

.. py:attribute:: pauli

  This is the bonding repulsion due to the Pauli principle.  Given
  in Hartrees.

.. py:attribute:: electrostatic

  This is the bonding energy due to electrostatic interactions.

.. py:attribute:: steric

  This is the sum of the Pauli and the electrostatic interactions.
  Given in Hartrees.

.. py:attribute:: orbital

  This is the orbital bonding energy.  It is in the form of a dictionary
  where the keys are the symmetry groups of the orbitals and 'Total'
  for the total orbital interactions.  Values are given in Hartrees.

.. py:attribute:: total

  This is the total bonding energy, in Hartrees.

.. py:attribute:: charges

  This is a numpy vector of length :py:attr:`natoms` containing the atomic
  charges partitioned using Voronoi analysis.

.. py:attribute:: qmcharge

  This is the charge of the QM system, given in atomic units.

.. py:attribute:: dipole

  This is a numpy array with the permanent dipole moment (debye) in
  vecor form.

.. py:attribute:: nmos

  The number of molecular orbitals (MOs) that were found.  If the
  calculation was unrestricted, then it will be a tuple, where spin
  1 (alpha) is first and and spin 2 (beta) is second.

.. py:attribute:: nocc

  The number of occupied molecular orbitals (MOs) that were found.  If the
  calculation was unrestricted, then it will be a tuple, where spin
  1 (alpha) is first and spin 2 (beta) is second. 

.. py:attribute:: nvir

  The number of virtual molecular orbitals (MOs) that were found.  If the
  calculation was unrestricted, then it will be a tuple, where spin
  1 (alpha) is first and spin 2 (beta) is second. 

.. py:attribute:: orbital_energies

  The energies, in eV, of each of the MOs.  It is in the form of a numpy
  vector of length :py:attr:`nmos`.

.. py:attribute:: orbital_ids

  The number and symmetry group of each MO in a numpy vector of length
  :py:attr:`nmos`.

.. py:attribute:: orbital_occ

  The occulations of each orbital as a numpy vector of length
  :py:attr:`nmos`.

.. py:attribute:: orbital_spin

  The spin of the orbitals.  0 is a restricted calculation, 1 is spin up
  or alpha, and 2 is spin down or beta.  It is a numpy vector of length
  :py:attr:`nmos`.

.. py:attribute:: atomic_orbitals

  The atomic orbitals that make up each MO.  It is given as a nested
  tuple with the major axis has length :py:attr:`nmos`.  The number of AOs
  in each MO is determined by the calculation.  Each MO consists of a
  numpy record array with 'pcent' is the percent contirbutions, 'ao_id'
  is the AO ids (what atom or fragment it is), and 'sym' is the
  symmetry groups the AOs belongs to.  Each numpy record for each MO
  is as long as the number of AOs for that MO. 

  .. note::
   In ADF, these are not technically atomic orbitals but 'fragment'
   orbitals.  Typically, these orbitals are of the atoms, so they are
   effectively atomic orbitals, but when a fragment analysis is done
   they belong to those fragments.

.. py:attribute:: deltas

  These are the dimensionless deltas needed for TDSPEC.  They are
  stored in a numpy vector of length :py:attr:`nmodes`.
   
.. py:attribute:: dgdip

  These are the derivatives of the permanent dipole moment used by 
  TDSPEC for DR-SFG calculations.  They are stored in a numpy array
  of dimensions :py:attr:`nmodes` x 3.

.. py:attribute:: dtdip

  These are the derivatives of the transition dipole moment used by
  TDSPEC for calculations involving Herzberg-Teller terms.  They are
  stored in a numpy array of dimensions :py:attr:`nexci` x `nmodes` x 3.

.. py:attribute:: HOMO

  A tuple containing the ID of the HOMO and it's energy in eV.
  If unrestricted, this has two indexes, the first being spin up,
  the second spin down.

.. py:attribute:: LUMO

  A tuple containing the ID of the LUMO and it's energy in eV.
  If unrestricted, this has two indexes, the first being spin up,
  the second spin down.

.. py:attribute:: SOMO

  A tuple containing the ID of the SOMO and it's energy in eV.
  There usually is no SOMO, and it is only valid for restrticted
  calculations.

.. py:attribute:: nexcite

  The number of excitation energies calculated.

.. py:attribute:: excitation_energies

  The excitation energies, in Hartrees, that were calculated.  It is
  in the form of a numpy vector of length :py:attr:`nexcite`.

.. py:attribute:: oscillator_strengths

  The oscillator strengths of the excitations, in a.u., as a numpy
  vector of length :py:attr:`nexcite`.

.. py:attribute:: opt_rot_strengths

  The optical rotatory strengths of the excitations, in a.u., as a 
  numpy vector of length :py:attr:`nexcite`.

.. py:attribute:: excitation_symmetries

  The symmetry groups of the excitations, as a numpy vector of length
  :py:attr:`nexcite`.

.. py:attribute:: excitation_type

  The type groups of the excitations, as a numpy vector of length
  :py:attr:`nexcite`.  Choices are 'SS' for singlet-singlet, or 'ST' for
  singlet-triplet.  Unrestricted calculations may show 'UN', or if
  it is implemented, will show 'SS' or 'ST'.

.. py:attribute:: TDM

  The transition dipole moments for each excitation, in a.u., given as
  a numpy array of dimensions :py:attr:`nexcite` x 3.

.. py:attribute:: MDM

  The magnetic transition dipole moments for each excitation, in a.u., 
  given as a numpy array of dimensions :py:attr:`nexcite` x 3.

.. py:attribute:: transitions

  This is a nested tuple of the orbital transitions for each excitation.
  It has length on the major axis of .. py:attribute:: nexcite.  Each index
  contains all the transitions listed from the calculation.  For each
  excitation, there is a numpy record array where 'occ' is the occupied
  orbitals, 'unocc' is the unoccupied orbitals, and 'pcent' is the percent
  contibution to the exctiation.  Each numpy record array is as long
  as the number of transitions for that excitation.

.. py:attribute:: STPM

  The two-photon transition moments for each excitation, in a.u., given
  as a numpy array of dimensions :py:attr:`nexcite` x 6.

.. py:attribute:: linear_tpa_strengths

  The two-photon transition strengths for linear polarized light, in a.u.,
  given as a numpy vector of length :py:attr:`nexcite`.

.. py:attribute:: linear_sigma_tpa
  
  The two-photon absorbance cross section (unbroadened) for linear
  polarized light, in Goeppert-Mayer (GM).  The conversion between 
  GM and c.g.s. units is: 1 GM = 1x10^{-50} cm^4 s / photon.  This
  is given as a numpy vector of length :py:attr:`nexcite`.

.. py:attribute:: circular_tpa_strengths

  The two-photon transition strengths for circular polarized light, in a.u.,
  given as a numpy vector of length :py:attr:`nexcite`.

.. py:attribute:: circular_sigma_tpa

  The two-photon absorbance cross section (unbroadened) for circularly
  polarized light, in Goeppert-Mayer (GM).  The conversion between
  GM and c.g.s. units is: 1 GM = 1x10^{-50} cm^4 s / photon.  This is
  given as a numpy vector of length :py:attr:`nexcite`.

.. py:attribute:: tpa_polarization_ratio

  The two-photon absorption polarization ratio, given as a numpy vector of
  length :py:attr:`nexcite`.

.. py:attribute:: T3PM

  The three-photon transition moments for each excitation, in a.u., given
  as a numpy array of dimensions :py:attr:`nexcite` x 27.

.. py:attribute:: linear_3pa_strengths

  The three-photon transition strengths for linear polarized light, in a.u.,
  given as a numpy vector of length :py:attr:`nexcite`.

.. py:attribute:: linear_sigma_3pa
  
  The three-photon absorbance cross section (unbroadened) for linear
  polarized light.  The c.g.s. units for a 3PA absorbance cross section
  are: cm^6 s^2 / photon.  This is given as a numpy vector of length 
  :py:attr:`nexcite`.

.. py:attribute:: circular_3pa_strengths

  The three-photon transition strengths for circular polarized light, in a.u.,
  given as a numpy vector of length :py:attr:`nexcite`.

.. py:attribute:: circular_sigma_tpa

  The three-photon absorbance cross section (unbroadened) for circularly
  polarized light.  The c.g.s. units for a 3PA absorbance cross section
  are: cm^6 s^2 / photon.  This is given as a numpy vector of length 
  :py:attr:`nexcite`.

.. py:attribute:: 3pa_polarization_ratio

  The three-photon absorption polarization ratio, given as a numpy vector of
  length :py:attr:`nexcite`.

.. py:attribute:: gs_gradient

  This is the ground state gradient wrt the nuclear coordinates as a
  dictionary containing numpy vectors of length :py:attr:`natoms` * 3.

.. py:attribute:: es

  For an excited state calculation, this is the excited state being
  calculated.

.. py:attribute:: es_dipole

  This is the excited state dipole moment as a 3 element numpy vector.

.. py:attribute:: es_gradient

  This is the excited state gradient wrt the optimized ground state
  geometry as a dictionary of numpy vectors of length :py:attr:`natoms` * 3.

.. py:attribute:: nmodes

  The number of normal modes that were calculated, if any.  This is also
  the degrees of freedom.

.. py:attribute:: v_frequencies

  The frequency of the normal modes (vibrational frequencies) as a numpy
  vector.  It has a length of :py:attr:`nmodes`
  and units of :math:`cm^{-1}`.

.. py:attribute:: normal_modes

  The coordinates of the normal modes, given with respect to the optimized
  geometry, i.e. the vectors for each atom indicating the movement for the
  mode.  It is a numpy array of dimensions
  :py:attr:`nmodes` x  :py:attr:`natoms` x 3.
  
  The normal modes are normalized.

.. py:attribute:: IR

   IR intensities, in :math:`\\frac{km}{mole}`, for the normal modes.
   This is given as a numpy vector of length :py:attr:`nmodes`.

.. py:attribute:: npol

  The number of polarizabilities calculated, if any.  If a Raman
  calculation was performed, then this is equal to :py:attr:`nmodes`.

.. py:attribute:: nhpol

  The number of hyperpolarizabilities calculated, if any.

.. py:attribute:: e_frequencies

  The electric frequencies of the polarizability calculations, as a
  numpy vector. It has a length of :py:attr:`npol` and units of Hartrees.

.. py:attribute:: b_e_frequencies

  The B (electric) frequencies of a frequency-dependent (second) hyperpolarizability 
  calculation, as a numpy vector. It has a length of :py:attr:`nhpol` and 
  units of Hartrees.

.. py:attribute:: c_e_frequencies

  The C (electric) frequencies of a frequency-dependent (second) hyperpolarizability 
  calculation, as a numpy vector. It has a length of :py:attr:`nhpol` and 
  units of Hartrees.

.. py:attribute:: d_e_frequencies

  The D (electric) frequencies of a frequency-dependent second hyperpolarizability 
  calculation, as a numpy vector. It has a length of :py:attr:`nhpol` and 
  units of Hartrees.

.. py:attribute:: ord

  The total optical rotation tensor for each frequency as a numpy array.
  It has dimensions of :py:attr:`npol` x 3 x 3.  If the imaginary part of
  the optical rotation was also calculated, the array will be complex.
  They are given in atomic units.

.. py:attribute:: polarizability

  The total polarizability tensor for each frequency as a numpy array.
  It has dimensions of :py:attr:`npol` x 3 x 3.  If the imaginary part of
  the polarizability was also calculated, the array will be complex.
  They are given in atomic units.

.. py:attribute:: polarizability_atomic

  The atomic polarizability tensors for each frequency as a numpy array.
  It has dimensions of :py:attr:`npol` x :py:attr:`natoms` x 3 x 3.  Atomic
  units are used, and they are partitioned with Mulliken analysis.

.. py:attribute:: hyperpolarizability

   The total hyperpolarizability tensor for the static system, given
   as a numpy float array of 3 x 3 x 3.

.. py:attribute:: dhpol

   The derivatives of the hyperpolarizability tensor for the static system, 
   given as a numpy float array of: nmodes x 3 x 3 x 3.

.. py:attribute:: hyperinvarients

   The hyperpolarizability invarients as a dictionary.

.. py:attribute:: secondhyperpolarizability

   The total secondhyperpolarizability tensor for the static system, given
   as a numpy float array of 3 x 3 x 3 x 3.

.. py:attribute:: dshpol

   The derivatives of the secondhyperpolarizability tensor for the static system, 
   given as a numpy float array of: nmodes x 3 x 3 x 3 x 3.

.. py:attribute:: ndim

  The number of DIM atoms in the calculation.

.. py:attribute:: ndim_elements

  The number of DIM elements in the calculation.

.. py:attribute:: dim_atoms

  The DIM atoms found in the calculation as a numpy vector of length
  :py:attr:`ndim`.

.. py:attribute:: dim_elements

  This is a set of the DIM elements in the calculation.  There are
  no duplicates.

.. py:attribute:: dim_coordinates

  The coordinates of the DIM atoms.  It is given as a numpy array with
  dimensions :py:attr:`ndim` x 3 and units of angstroms.

.. py:attribute:: dim_charges

  A numpy vector of length :py:attr:`ndim` that contains the induced
  charges of the DIM system.

.. py:attribute:: dim_dipoles

  A numpy array of dimensions :py:attr:`ndim` x 3 that contains the induced
  dipoles of the DIM system.

.. py:attribute:: dim_dipole_tot

  This is a numpy array with the total dipole moment of the DIM
  system (Debyes) in vecor form.

.. py:attribute:: dim_pol

  A numpy array of dimensions :py:attr:`npol` x 3 x 3 that contains the
  polarizability tensor of the DIM system alone for each frequency.

.. py:attribute:: dim_efficiencies

  A natomsx3 numpy array containing the absorption, scattering and
  extinction efficiencies for each frequency. This is a unitless quantity.

.. py:attribute:: dim_cross_sections

  A natomsx3 numpy array containing the absorption, scattering and
  extinction ross-sections in units of NM$^3$ for each freqeuncy

.. py:attribute:: dim_energy

  A float conatining the DIM system energy.

.. py:attribute:: dimqm_energy

  A float conatining the DIM/QM interaction energy.

.. py:attribute:: OrbitalPlot

  A numpy array containing the MO information from a Gaussian cube file.
  This has variable dimension.

.. py:attribute:: OrbitalDim

  A numpy vector of length 3 defining the number of points used to generate
  an MO from a Gaussian cube file.

.. py:attribute:: spin_multiplicity

  An integer containing the spin multiplicity of the calculation.

.. py:attribute:: tdspec_wavelength

  A numpy vector containing the wavelength scan from a TDSPEC calculation
  for one-photon absorption (OPA) or two-photon absorption (TPA)

.. py:attribute:: tdspec_wavenumber

  A numpy vector containing the wavenumbers (Raman shifts) from a TDSPEC
  calculation for resonance Raman scattering (RRS), resonance hyper-Raman
  scattering (RHRS), or doubly-resonant sum-frequency generation (DR-SFG).

.. py:attribute:: tdspec_opa

  A numpy vector containing OPA cross sections from TDSPEC.

.. py:attribute:: tdspec_rrs

  A numpy vector containing RRS differential cross sections from TDSPEC.

.. py:attribute:: tdspec_tpa

  A numpy vector containing TPA cross sections from TDSPEC.

.. py:attribute:: tdspec_rhrs

  A numpy vector containing RHRS differential cross sections from TDSPEC.

.. py:attribute:: tdspec_drsfg

  A numpy vector containing DR-SFG squared hyperpolarizabilities from TDSPEC.

Methods
^^^^^^^
   
'''

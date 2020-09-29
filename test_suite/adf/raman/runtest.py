from chem import collect

d=collect('raman-benzene.orig.out')

print('Number of Normal Modes:\n%d'%(d.nmodes))
print('Array Shape of Polarizability Derivatives:')
print(d.polarizability_derivatives.shape)

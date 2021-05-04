Subroutine rotate (rotmat, coords, natoms)

!  Rotate coordinates.  Note that we use the overwrite feature
!  of f2py to 1) prevent using more memory than needed and
!  2) reduce the time required to pass data from python to FORTRAN

   use Constants

   Implicit None

   Integer,     Intent(In)  :: natoms
   Real(KINDR), Intent(In)  :: rotmat(3,3)
   Real(KINDR)              :: coords(natoms,3)
!f2py Real(KINDR), Intent(In,Out,Overwrite) :: coords(natoms,3)

   Real(KINDR)              :: dummy(natoms,3)

!  Use fast BLAS matrix multiply to get it done
   call dgemm ('N', 'T', natoms, 3, 3, ONE, coords, natoms, rotmat, 3, ZERO, dummy, natoms)
!  Place new coordinates into output.
   coords = dummy

End Subroutine rotate

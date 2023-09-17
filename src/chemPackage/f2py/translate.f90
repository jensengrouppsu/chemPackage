Subroutine translate (transvec, coords, natoms)

!  Translate coordinates. Note that we use the overwrite feature
!  of f2py to 1) prevent using more memory than needed and
!  2) reduce the time required to pass data from python to FORTRAN

   use Constants

   Implicit None

   Integer,     Intent(In)  :: natoms
   Real(KINDR), Intent(In)  :: transvec(3)
   Real(KINDR)              :: coords(natoms, 3)
!f2py Real(KINDR), Intent(In,Out,Overwrite) :: coords(natoms,3)

!  For each atom, each component is shifted the same,
!  so do all X's, all Y's, and all Z's.
   coords(:,1) = coords(:,1) + transvec(1)
   coords(:,2) = coords(:,2) + transvec(2)
   coords(:,3) = coords(:,3) + transvec(3)

End Subroutine translate

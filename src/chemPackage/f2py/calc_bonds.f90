Subroutine calc_bonds (coords, radii, scalefac, natoms, indices, nbonds)

!  Calculate which atom pairs have a bond

   use Constants

   Implicit None

   Integer,     Intent(In)  :: natoms
   Real(KINDR), Intent(In)  :: coords(3,natoms)
   Real(KINDR), Intent(In)  :: radii(natoms)
   Real(KINDR), Intent(In)  :: scalefac
   Integer,     Intent(Out) :: indices(2,(natoms*(natoms-1))/2)
   Integer,     Intent(Out) :: nbonds

   Integer                  :: k, m, n
   Real(KINDR)              :: d, b

   Real(KINDR),  External   :: dnrm2

!  Loop over each of the atom pairs in one loop
   nbonds = 0
   m = 1
   n = m
   do k = 1, (natoms*(natoms-1))/2

!     Increment n
      n = n + 1

!     Find distance between points
      d = dnrm2(3, coords(:,n) - coords(:,m), 1)
!     Find the maximum bond distance
      b = scalefac * ( radii(m) + radii(n) )

!     If the distance is less than the radii between the atoms
!     (times a scaling factor) then it is a bond, and keep it
      if (d < b) then
         nbonds = nbonds + 1
         indices(1,nbonds) = m - 1
         indices(2,nbonds) = n - 1
      end if

!     If n is on the last atom, increment m and n
      if (n == natoms) then
         m = m + 1
         n = m
      end if

   end do

End Subroutine calc_bonds

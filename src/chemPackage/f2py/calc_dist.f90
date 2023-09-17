Subroutine calc_dist (coords, natoms, point, dists)

!  Calculate the distances from a point in space

   use Constants

   Implicit None

   Integer,     Intent(In)  :: natoms
   Real(KINDR), Intent(In)  :: coords(3,natoms)
   Real(KINDR), Intent(In)  :: point(3)
   Real(KINDR), Intent(Out) :: dists(natoms)

   Integer                  :: n
   Real(KINDR)              :: d

   Real(KINDR),  External   :: dnrm2

!  Loop over each atom, and find distance between each point (using norm)
   do n = 1, natoms
      dists(n) = dnrm2(3, point - coords(:,n), 1)
   end do

End Subroutine calc_dist

Subroutine minmax_pdist (coord, nc, mindist, maxdist)

   use Constants

   Implicit None

!  Calculate the minimum and maximum distance between a set of points

   Integer,     Intent(In)  :: nc
   Real(KINDR), Intent(In)  :: coord(3,nc)
   Real(KINDR), Intent(Out) :: mindist
   Real(KINDR), Intent(Out) :: maxdist

   Integer                  :: k, m, n
   Real(KINDR)              :: d

   Real(KINDR),  External   :: dnrm2

!  Loop over each of the atom pairs in one loop
   maxdist = TINY(maxdist)
   mindist = HUGE(mindist)
   m = 1
   n = m
   do k = 1, ( nc * ( nc - 1 ) ) / 2

!     Increment n
      n = n + 1

!     Find distance between points
      d = dnrm2(3, coord(:,n) - coord(:,m), 1)

!     If the distance is outside the current extrema, make it the extrema
      if (d < mindist) mindist = d
      if (d > maxdist) maxdist = d

!     If n is on the last atom of coord2, increment m and reset n
      if (n == nc) then
         m = m + 1
         n = m
      end if

   end do

End Subroutine minmax_pdist

Subroutine minmax_cdist (coord1, coord2, nc1, nc2, mindist, maxdist)

   use Constants

   Implicit None

!  Calculate the minimum and maximum distance between two sets of points

   Integer,     Intent(In)  :: nc1
   Integer,     Intent(In)  :: nc2
   Real(KINDR), Intent(In)  :: coord1(3,nc1)
   Real(KINDR), Intent(In)  :: coord2(3,nc2)
   Real(KINDR), Intent(Out) :: mindist
   Real(KINDR), Intent(Out) :: maxdist

   Integer                  :: k, m, n, span
   Real(KINDR)              :: d, b

   Real(KINDR),  External   :: dnrm2

!  Loop over each of the atom pairs in one loop
   maxdist = TINY(maxdist)
   mindist = HUGE(mindist)
   m = 1
   n = 1
   do k = 1, ( nc1 * nc2 ) / 2

!     Find distance between points
      d = dnrm2(3, coord2(:,n) - coord1(:,m), 1)

!     If the distance is outside the current extrema, make it the extrema
      if (d < mindist) mindist = d
      if (d > maxdist) maxdist = d

!     Increment n
      n = n + 1

!     If n is on the last atom of coord2, increment m and reset n
      if (n == nc2) then
         m = m + 1
         n = 1
      end if

   end do

End Subroutine minmax_cdist

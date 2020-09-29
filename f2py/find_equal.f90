Subroutine find_equal (array1, array2, arrayA, arrayB, n1, n2, indx)

!  Function to find where two data sets are equal, using two
!  arrays of comparable data for each.
!
!  Array1 and 2 correspond, and array 3 and 4 correspond. The length of
!  array 1 and 3 are the same, and 2 and 4 are the same length.  Arrays
!  1 and 3 are the same length or shorter than arrays 2 and 4, and the
!  index returned is the same length as arrays 2 and 4.

   use Constants

   Implicit None

   Integer,     Intent(In)  :: n1
   Integer,     Intent(In)  :: n2
   Real(KINDR), Intent(In)  :: array1(n1)
   Real(KINDR), Intent(In)  :: array2(n2)
   Real(KINDR), Intent(In)  :: arrayA(n1)
   Real(KINDR), Intent(In)  :: arrayB(n2)
   Integer,     Intent(Out) :: indx(n2)

   Logical                  :: skip(n2)
   Integer                  :: i, j

!  Initiallize
   skip = .false.
   indx = -1

!  Loop over each of the atom pairs in one loop
   do i = 1, n1
      do j = 1, n2

!        Skip if this one has been found
         if (skip(j)) cycle

!        See if Array1 and Array2 match
         if (approx_equal(array1(i), array2(j))) then
!           If the above match, see if arrayA and arrayB match
            if (approx_equal(arrayA(i), arrayB(j))) then
!              If both match, save this index
               indx(j) = i - 1
               skip(j) = .true.
               exit
            end if 
         end if

      end do
   end do

Contains

   Logical Function approx_equal (one, two)

      Real(KINDR), Intent(In) :: one
      Real(KINDR), Intent(In) :: two

      Real(KINDR)             :: s
      Real(KINDR)             :: scaled_one
      Real(KINDR)             :: scaled_two

      Real(KINDR), Parameter  :: TEN  = 10.0_KINDR
      Real(KINDR), Parameter  :: HALF = 0.5_KINDR

!     See if they are equal
      if (one == two) then
         approx_equal = .true.
!     See if they are close
      else
!        Normalize numbers to be in the range (-10.0, 10.0)
         s = TEN**FLOOR(LOG10(HALF * (ABS(one) + ABS(two))))
         scaled_one = one / s
         scaled_two = two / s
!        Now check if they are approx equal
         approx_equal = ABS(scaled_one - scaled_two) < TEN**(-3)
      end if

   End Function approx_equal

End Subroutine find_equal

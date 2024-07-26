Subroutine error_function (y, erf)

!  ======================================================================
!  purpose:  evaluation of the error function.
!
!  input  :  y      - the argument of the function
!  output :  erf    - the function value
!
!  This was lifted directly from ADF's library, and was only edited
!  insofar as addiong the implicit none (and therefore declaring all
!  variables) adding intent in and out, and cleaning up the parameter
!  declarations.  I take no blame for the lack of comments or the
!  extensive use of GOTOs.
!  ======================================================================

   use Constants

   Implicit None

   Real(KINDR), Intent(In)  :: y
   Real(KINDR), Intent(Out) :: erf

   Real(KINDR) :: x
   Real(KINDR) :: xx
   Real(KINDR) :: md
   Real(KINDR) :: sn
   Real(KINDR) :: sd

   Integer :: i
   Integer :: isw
   Integer :: kret

   Real(KINDR), Parameter :: C     = 1.1283791670955_KINDR
   Real(KINDR), Parameter :: R0P47 = 0.47_KINDR
   Real(KINDR), Parameter :: R5P5  = 5.5_KINDR
   Real(KINDR), Parameter :: SMALL = 1.0E-15_KINDR
   Real(KINDR), Parameter :: SQR1D2 =  0.70710678118655_KINDR

   Real(KINDR), Parameter :: a(10) = (/ONE,                      &
                                       THIRD,                    &
                                       0.1_KINDR,                &
                                      -0.02380923809238_KINDR,   &
                                       0.4629629629630E-2_KINDR, &
                                     - ONE / 1320.0_KINDR,       &
                                       1.0683760683761E-4_KINDR, &
                                      -1.3227513227513E-5_KINDR, &
                                       1.4589169000934E-6_KINDR, &
                                      -1.45038552223150E-7_KINDR/)

   Real(KINDR), Parameter :: p(8) = (/883.4789426085_KINDR,  &
                                      1549.6793124037_KINDR, &
                                      1347.1941340976_KINDR, &
                                      723.04000277753_KINDR, &
                                      255.50049469496_KINDR, &
                                      59.240010112914_KINDR, &
                                      8.3765310814197_KINDR, &
                                      0.56418955944261_KINDR/)

   Real(KINDR), Parameter :: q(9) = (/883.4789426085_KINDR,  &
                                      2546.5785458098_KINDR, &
                                      3337.2213699893_KINDR, &
                                      2606.7120152651_KINDR, &
                                      1333.5699756800_KINDR, &
                                      460.28512369160_KINDR, &
                                      105.50025439769_KINDR, &
                                      14.847012237523_KINDR, &
                                      ONE/)

   kret = 2

   x  = y
   md = 0

   if (x<ZERO) then
      isw = 1
      x   = - x
   else
      isw = 2
   end if

   if (kret/=2) goto 120
   if (x<=R5P5) goto 40

   erf = ONE
   if (isw==1) erf = - ONE
   goto 150

!  -----------------------------------------------------------
!  abs(x) less than .47 compute erf by taylor series expansion
!  -----------------------------------------------------------

 40 if (x<=R0P47) goto 90
   kret = 1

!  -----------------------------------------------------
!  abs(x) between .47 and 10. compute complemented error
!  function by a rational function of x
!  -----------------------------------------------------

 50 sn = p(8)
   do i = 7, 1, -1
      sn = sn*x + p(i)
   end do 

   sd = q(9)
   do i = 8, 1, -1
      sd = sd*x + q(i)
   end do 

   erf = (sn/sd)*exp(-x*x)
   if (kret/=1) goto 80

!  -----------------------------------
!  compute complemented error function
!  -----------------------------------

   erf = ONE - erf
   if (isw/=2) erf = - erf
   goto 150

 80 if (isw/=2) erf = ONE - erf
   if (md/=0) erf = HALF*erf
   goto 150

 90 xx = x*x

   erf = a(10)
   do i = 9, 1, -1
      erf = erf*xx + a(i)
   end do 

   erf = C*erf*x

   if (kret==2) then
      if (isw/=2) erf = - erf
      goto 150
   else
      erf = one - erf
      goto 80
   end if

 120 if (x>EIGHT) then
      erf = ZERO
      goto 80

   else if (x>R0P47) then
      kret = 2
      goto 50

   else if (x>SMALL) then
      goto 90

   else
      erf = ONE
      goto 80
   end if

 150 return

End Subroutine error_function

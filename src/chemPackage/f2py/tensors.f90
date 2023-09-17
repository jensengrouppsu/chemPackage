  subroutine search(list, nl, test, lsearch)
    use Constants
    integer,     intent(in)  :: nl
    real(kindr), intent(in)  :: list(nl)
    real(kindr), intent(in)  :: test
    logical,     intent(out) :: lsearch
    integer                  :: i
    lsearch = .false.
    do i = 1, nl
        if (list(i) <= test) then
            lsearch = .true.
            exit
        end if
    end do
  end subroutine search


  subroutine t0_unscreened (invdist, t)
    use Constants
    Real(KINDR), Intent(In)  :: invdist ! 1 / dist
    Real(KINDR), Intent(out) :: t

    t = ZERO

!   1 / dist
    t = invdist

  End subroutine t0_unscreened

  subroutine t0_ret (dist, invdist, kret, t) 
    use Constants
    real(kindr), intent(in)     :: dist
    real(kindr), intent(in)     :: invdist
    Real(kindr), intent(In)     :: kret
    complex(kindr), intent(out) :: t

    complex(kindr)          :: wave_k

    t = ZERO_C

    wave_k = exp(I_C*kret*dist)
    t = wave_k * invdist
  
  end subroutine t0_ret

  subroutine t0_screened (dist, invdist, scrn1, scrn2, t)
    use Constants
    Real(KINDR), Intent(In)  :: dist    ! Length of vector r
    Real(KINDR), Intent(In)  :: invdist ! 1 / dist
    Real(KINDR), Intent(In)  :: scrn1      ! Atom 1
    Real(KINDR), Intent(In)  :: scrn2      ! Atom 2
    Real(KINDR), Intent(out) :: t

    Real(KINDR) :: invRqq

    t = ZERO

    invRqq = ONE / SQRT(scrn1**2 + scrn2**2)
    t      = ERF(dist * invRqq) * invdist

  End subroutine t0_screened

  subroutine  t0_ret_screened (dist, invdist, scrn1, scrn2, kret, t)
    use Constants
    Real(KINDR), Intent(In) :: dist    ! Length of vector r
    Real(KINDR), Intent(In) :: invdist ! 1 / dist
    Real(KINDR), Intent(In) :: scrn1      ! Atom 1
    Real(KINDR), Intent(In) :: scrn2      ! Atom 2
    Real(KINDR), Intent(In) :: kret    ! k scalar
    Complex(KINDR)          :: t
    Real(KINDR)             :: t_0


    ! T^{(0)}_{ret+scrn} = T^{(0)}_{ret} + T^{(0)}_{scrn} -T^{(0)}_{dda}
    call t0_ret(dist, invdist, kret, t)     ! t1 without screening and with    retardation
    call t0_screened(dist, invdist, scrn1, scrn2, t_0)   ! t1 with    screening and without retardation
    t   = t + cmplx(t_0, ZERO, KINDR)
    call t0_unscreened(invdist, t_0)                 ! t1 without screening and without retardation
    t   = t - cmplx(t_0, ZERO, KINDR)

  End subroutine t0_ret_screened

  subroutine t1_unscreened (r, invdist, t) 
    use Constants
    Real(KINDR), Intent(In)  :: r(3)    ! Vector between atom 1 and 2
    Real(KINDR), Intent(In)  :: invdist ! 1 / dist
    Real(KINDR), intent(out) :: t(3)

    t = ZERO
    
!   - r_a / dist^3
    t(1:3) = - r(1:3) * invdist * invdist * invdist

  End subroutine t1_unscreened

  subroutine t1_ret (r, dist, invdist, kret, t) 
    use Constants
    Real(KINDR), Intent(In)     :: r(3)    ! Vector between atom 1 and 2
    Real(KINDR), Intent(In)     :: dist    ! Length of vector r
    Real(KINDR), Intent(In)     :: invdist ! 1 / dist
    Real(KINDR), Intent(In)     :: kret    ! k scalar
    Complex(KINDR), Intent(Out) :: t(3)

    Complex(KINDR)          :: wave_k 

    t = ZERO_C

    wave_k = ( EXP(I_C*kret*dist) )
    t(1:3) = -r(1:3) * invdist * invdist * invdist * wave_k*(1-I_C*kret*dist)

  End subroutine t1_ret

  
  subroutine t1_screened (r, dist, invdist, scrn1, scrn2, t) 
    use Constants
    Real(KINDR), Intent(In)  :: r(3)    ! Vector between atom 1 and 2
    Real(KINDR), Intent(In)  :: dist    ! Length of vector r
    Real(KINDR), Intent(In)  :: invdist ! 1 / dist
    Real(KINDR), Intent(In)  :: scrn1      ! Atom 1
    Real(KINDR), Intent(In)  :: scrn2      ! Atom 2
    Real(KINDR), Intent(out) :: t(3)

    Real(KINDR) :: invRpq
    Real(KINDR) :: rtmp

    t = ZERO

  !   Calculate R values
    invRpq = ONE / SQRT(scrn1**2 + scrn2**2)
  !   Calculate the screening value and place directly into the tensor
    rtmp   = ( ERF(dist * invRpq) - TWOINVSQRTPI &
           * ( dist * invRpq ) * EXP(- ( dist * invRpq )**2) )
    t(1:3) = - rtmp * r(1:3) * invdist * invdist * invdist

  End subroutine t1_screened

  subroutine t1_ret_screened (r, dist, invdist, scrn1, scrn2, kret, t)
    use Constants
    Real(KINDR), Intent(In)     :: r(3)    ! Vector between atom 1 and 2
    Real(KINDR), Intent(In)     :: dist    ! Length of vector r
    Real(KINDR), Intent(In)     :: invdist ! 1 / dist
    Real(KINDR), Intent(In)     :: scrn1      ! Atom 1
    Real(KINDR), Intent(In)     :: scrn2      ! Atom 2
    Real(KINDR), Intent(In)     :: kret    ! k scalar
    Complex(KINDR), Intent(out) :: t(3)
    Real(KINDR)                 :: t_1(3)

    ! T^{(1)}_{ret+scrn} = T^{(1)}_{ret} + T^{(1)}_{scrn} -T^{(1)}_{dda}
    call t1_ret(r, dist, invdist, kret, t)   ! t1 without screening and with    retardation
    call t1_screened(r, dist, invdist, scrn1, scrn2, t_1) ! t1 with    screening and without retardation
    t   = t + cmplx(t_1, ZERO, KINDR)
    call t1_unscreened(r, invdist, t_1)               ! t1 without screening and without retardation
    t   = t - cmplx(t_1, ZERO, KINDR)
  
  end subroutine  t1_ret_screened

  subroutine t2_unscreened (r, invdist, t)
    use Constants
    Real(KINDR), Intent(In)  :: r(3)    ! Vector between atom 1 and 2
    Real(KINDR), Intent(In)  :: invdist ! 1 / dist
    Real(KINDR), Intent(out) :: t(3,3)

    Real(KINDR) :: invdist3
    Real(KINDR) :: invdist5_3

    t = ZERO

!   Powers of the inverse distance
!   Include the prefactors into the distance for speedup
    invdist3   = invdist  * invdist * invdist
    invdist5_3 = invdist3 * invdist * invdist * THREE

!   Now fill in the matrix in the order it is laid out in memory.
    t(1,1) = r(1)**2 * invdist5_3 - invdist3
    t(2,1) = r(1) * r(2) * invdist5_3
    t(3,1) = r(1) * r(3) * invdist5_3
    t(1,2) = t(2,1)
    t(2,2) = r(2)**2 * invdist5_3 - invdist3
    t(3,2) = r(2) * r(3) * invdist5_3
    t(1,3) = t(3,1)
    t(2,3) = t(3,2)
    t(3,3) = r(3)**2 * invdist5_3 - invdist3

  End subroutine t2_unscreened


  subroutine t2_ret (r, dist, invdist, kret, t) 
    use Constants
    Real(KINDR), Intent(In)     :: r(3)    ! Vector between atom 1 and 2
    Real(KINDR), Intent(In)     :: dist    ! Length of vector r
    Real(KINDR), Intent(In)     :: invdist ! 1 / dist
    Real(KINDR), Intent(In)     :: kret    ! k scalar
    Complex(KINDR), Intent(out) :: t(3,3)

    Real(KINDR)      :: dist2
    Real(KINDR)      :: invdist5
    Complex(KINDR)   :: temp1
    Complex(KINDR)   :: temp2
    Complex(KINDR)   :: temp3
    Complex(KINDR)   :: temp4
    Complex(KINDR)   :: temp5
    Complex(KINDR)   :: wave_k

    t = ZERO_C

!   Powers of the inverse distance and distance
    dist2    = dist * dist
    invdist5 = invdist**5

    wave_k = EXP(I_C*kret*dist)

    temp1 = wave_k*invdist5
    temp2 = (kret**2)*dist2
    temp3 = ONE - I_C*kret*dist 
    temp4 = temp1*(-temp2+three*temp3)
    temp5 = temp1*dist2*(temp2-temp3)
        
!   Now fill in the matrix in the order it is laid out in memory.
    t(1,1) = (r(1)*r(1)*temp4 + temp5)
    t(2,1) = r(1)*r(2)*temp4
    t(3,1) = r(1)*r(3)*temp4
    t(1,2) = t(2,1)
    t(2,2) = (r(2)*r(2)*temp4 + temp5)
    t(3,2) = r(2)*r(3)*temp4
    t(1,3) = t(3,1)
    t(2,3) = t(3,2)
    t(3,3) = (r(3)*r(3)*temp4 + temp5)

  End subroutine t2_ret
  
   subroutine t2_screened (r, dist, invdist, scrn1, scrn2, t) 
    use Constants
    Real(KINDR), Intent(In)  :: r(3)    ! Vector between atom 1 and 2
    Real(KINDR), Intent(In)  :: dist    ! Length of vector r
    Real(KINDR), Intent(In)  :: invdist ! 1 / dist
    Real(KINDR), Intent(In)  :: scrn1      ! Atom 1
    Real(KINDR), Intent(In)  :: scrn2      ! Atom 2
    Real(KINDR), Intent(out) :: t(3,3)

    Real(KINDR) :: dist3
    Real(KINDR) :: invdist3
    Real(KINDR) :: invdist5
    Real(KINDR) :: temp1
    Real(KINDR) :: temp2
    Real(KINDR) :: Rpp
    Real(KINDR) :: invRpp

    t = ZERO

    !   Powers of the inverse distance and distance
    dist3    = dist * dist * dist
    invdist3 = invdist  * invdist * invdist
    invdist5 = invdist3 * invdist * invdist

    Rpp    = SQRT(scrn1**2 + scrn2**2)
    invRpp = ONE / Rpp
    
    !   Calculate the screening values
    temp1  = ERF(dist * invRpp) - TWOINVSQRTPI &
             * ( dist * invRpp ) * EXP(- ( dist * invRpp )**2)
    temp2  = FOUR / ( SQRTPI * Rpp**3 ) * EXP(- ( dist * invRpp )**2)
    !   Include these prefactors into the inverse distance terms
    invdist5 = invdist5 * ( temp1 * THREE - temp2 * dist3 )
    invdist3 = invdist3 * temp1

    !   Now fill in the matrix in the order it is laid out in memory.
    t(1,1) = r(1)**2 * invdist5 - invdist3
    t(2,1) = r(1) * r(2) * invdist5
    t(3,1) = r(1) * r(3) * invdist5
    t(1,2) = t(2,1)
    t(2,2) = r(2)**2 * invdist5 - invdist3
    t(3,2) = r(2) * r(3) * invdist5
    t(1,3) = t(3,1)
    t(2,3) = t(3,2)
    t(3,3) = r(3)**2 * invdist5 - invdist3

  End subroutine t2_screened



  subroutine t2_ret_screened (r, dist, invdist, scrn1, scrn2, kret, t) 
    use Constants
    Real(KINDR), Intent(In)     :: r(3)    ! Vector between atom 1 and 2
    Real(KINDR), Intent(In)     :: dist    ! Length of vector r
    Real(KINDR), Intent(In)     :: invdist ! 1 / dist
    Real(KINDR), Intent(In)     :: scrn1      ! Atom 1
    Real(KINDR), Intent(In)     :: scrn2      ! Atom 2
    Real(KINDR), Intent(In)     :: kret    ! k scalar
    Complex(KINDR), Intent(out) :: t(3,3)
    Real(KINDR)                 :: t_2(3,3)

    ! T^{(2)}_{ret+scrn} = T^{(2)}_{ret} + T^{(2)}_{scrn} -T^{(2)}_{dda}
    call t2_ret(r,dist,invdist,kret, t)     ! t2 without screening and with    retardation
    call t2_screened(r,dist,invdist,scrn1,scrn2, t_2)    ! t2 with    screening and without retardation
    t   = t + cmplx(t_2, ZERO, KINDR)   
    call t2_unscreened(r,invdist, t_2)               ! t2 without screening and without retardation
    t   = t - cmplx(t_2, zero, KINDR)   

  End subroutine t2_ret_screened

subroutine drawing_math(x, nx, y, ny, z, nz, coord, nc1, nc2, dipole, nd1, nd2, &
                atoms, na, efield, ne, charge, nc, zeros, dir, scrn, smear)
    use Constants
  !$ use omp_lib
    implicit none
    ! Input variables
    integer,        intent(in) :: nx, ny, nz, nc1, nc2, nd1, nd2, na, ne, nc
!f2py integer, intent(hide) :: nx, ny, nz, nc1, nc2, nd1, nd2, na, ne, nc
    real(kindr),    intent(in) :: charge(nc), x(nx), y(ny), z(nz)
    real(kindr),    intent(in) :: coord(nc1,nc2), atoms(na)
    complex(kindr), intent(in) :: dipole(nd1, nd2)
    logical,        intent(in) :: zeros
    integer,        intent(in) :: dir
!f2py integer, intent(in) :: dir
    integer,        intent(in) :: scrn 
    real(kindr),    intent(in) :: smear
    real(kindr)                :: t1(3), t2(3,3)
    complex(kindr)             :: t_1(3), t_2(3,3)
    real(kindr)                :: xyz(3), scrn1,  dist, invdist

    ! Internal Variables
    real(kindr), parameter     :: const=0.001_KINDR
    real(kindr)                :: nan          
    real(kindr)    :: rx(nc2), ry(nc2), rz(nc2), r(nc2)
    complex(kindr) :: ex, ey, ez
    integer :: i, iAtom, idir, idir2
    logical :: lsearch1, lsearch2
    ! Output variables
    complex(kindr) :: efield(ne)

    ! Zero out the output
    efield = ZERO_C    

    nan = 0.0_KINDR
   !const=0.001_KINDR

!$OMP PARALLEL private(i,rx,ry,rz,r, ex,ey,ez, iAtom, t_1, t1, t_2, t2, xyz, scrn1, dist, invdist, lsearch1, lsearch2, idir, idir2)  
!$OMP DO
    do i = 1, nx
        ! Determine distance from each atom this point in space
        rx = x(i) - coord(1,:)
        ry = y(i) - coord(2,:)
        rz = z(i) - coord(3,:)
        r  = sqrt(rx**2 + ry**2 + rz**2)

        call search(r, nc2, const, lsearch1)
        call search(r, nc2, atoms(1), lsearch2)

        if (lsearch1) then
          if (zeros) then
            efield(i) = ZERO
          else
            efield(i) = nan/nan ! 0.0/0.0
          end if
          cycle
        else if (scrn /=1 .and. scrn /= 2 .and. lsearch2 ) then 
          if (zeros) then
            efield(i) = ZERO
          else
            efield(i) = nan/nan ! 0.0/0.0
          end if
          cycle
        else 
          ex = ZERO
          ey = ZERO
          ez = ZERO
          if (scrn == 1) then
            do iAtom = 1, na
              xyz =  (/rx(iAtom),ry(iAtom),rz(iAtom) /)
              scrn1 = atoms(iAtom)
              dist = r(iAtom)
              invdist = ONE/dist
              call t2_screened(xyz, dist, invdist, scrn1, smear, t2)
              t_2 = cmplx(t2, ZERO, KINDR)

              ex = ex + t_2(1,1)*dipole(1, iAtom) & 
             &        + t_2(1,2)*dipole(2, iAtom) &
             &        + t_2(1,3)*dipole(3, iAtom)
              ey = ey + t_2(2,1)*dipole(1, iAtom) & 
             &        + t_2(2,2)*dipole(2, iAtom) &
             &        + t_2(2,3)*dipole(3, iAtom)
              ez = ez + t_2(3,1)*dipole(1, iAtom) & 
             &        + t_2(3,2)*dipole(2, iAtom) &
             &        + t_2(3,3)*dipole(3, iAtom)

              !Charge 
              call t1_screened(xyz, dist, invdist, scrn1, smear, t1)
              t_1 = cmplx(t1, ZERO, KINDR)
              ex = ex + t_1(1) * charge(iAtom) 
              ey = ey + t_1(2) * charge(iAtom) 
              ez = ez + t_1(3) * charge(iAtom) 
              
            end do ! iAtom

          else ! without screening
            do iAtom = 1, na
              xyz =  (/rx(iAtom),ry(iAtom),rz(iAtom) /)
              dist = r(iAtom)
              invdist = ONE/dist
              call t2_unscreened(xyz, invdist, t2)
              t_2 = cmplx(t2, ZERO, KINDR)

              ex = ex + t_2(1,1)*dipole(1, iAtom) & 
             &        + t_2(1,2)*dipole(2, iAtom) &
             &        + t_2(1,3)*dipole(3, iAtom)
              ey = ey + t_2(2,1)*dipole(1, iAtom) & 
             &        + t_2(2,2)*dipole(2, iAtom) &
             &        + t_2(2,3)*dipole(3, iAtom)
              ez = ez + t_2(3,1)*dipole(1, iAtom) & 
             &        + t_2(3,2)*dipole(2, iAtom) &
             &        + t_2(3,3)*dipole(3, iAtom)

              !Charge 
              call t1_unscreened(xyz, invdist, t1)
              t_1 = cmplx(t1, ZERO, KINDR)
              ex = ex + t_1(1) * charge(iAtom) 
              ey = ey + t_1(2) * charge(iAtom) 
              ez = ez + t_1(3) * charge(iAtom) 
            end do !iAtom
          end if

        ! Add in external field contribution
          if(dir .eq. 0) then
            ex = ex + 1
          else if(dir .eq. 1) then
            ey = ey + 1
          else if(dir .eq. 2) then        
            ez = ez + 1
          else if (dir .eq. 4) then
            continue
          end if
          ! Store the magnitude squared of the local field vector
          !efield(i) = (ez*conjg(ez))
          efield(i) = (ex*conjg(ex) + ey*conjg(ey) + ez*conjg(ez))

        end if
    end do
!$OMP END DO
!$OMP END PARALLEL 
    return 
end subroutine drawing_math

subroutine drawing_math_vectors(x, nx, y, ny, z, nz, coord, nc1, nc2, dipole, nd1, nd2, &
                atoms, na, efieldx, nex, efieldy, ney, efieldz, nez, &
                charge, nc, zeros, dir, scrn, smear)
    use Constants
  !$ use omp_lib
    implicit none
    ! Input variables
    integer,        intent(in) :: nx, ny, nz, nc1, nc2, nd1, nd2, na, nc, nex, ney, nez
!f2py integer, intent(hide) :: nx, ny, nz, nc1, nc2, nd1, nd2, na, ne, nc
    real(kindr),    intent(in) :: charge(nc), x(nx), y(ny), z(nz)
    real(kindr),    intent(in) :: coord(nc1,nc2), atoms(na)
    complex(kindr), intent(in) :: dipole(nd1, nd2)
    logical,        intent(in) :: zeros
    integer,        intent(in) :: dir
!f2py integer, intent(in) :: dir
    integer,        intent(in) :: scrn 
    real(kindr),    intent(in) :: smear
    real(kindr)                :: t1(3), t2(3,3)
    complex(kindr)             :: t_1(3), t_2(3,3)
    real(kindr)                :: xyz(3), scrn1,  dist, invdist

    ! Internal Variables
    real(kindr), parameter     :: const=0.001_KINDR
    real(kindr)                :: nan          
    real(kindr)    :: rx(nc2), ry(nc2), rz(nc2), r(nc2)
    complex(kindr) :: ex, ey, ez
    integer :: i, iAtom
    logical :: lsearch1, lsearch2
    
       
    ! Output variables
    complex(kindr) :: efieldx(nex), efieldy(ney), efieldz(nez)

    ! Zero ou the outputs
    efieldx = ZERO_C
    efieldy = ZERO_C
    efieldz = ZERO_C

   !const=0.001_KINDR
    nan = 0.0_KINDR

!$OMP PARALLEL private(rx,ry,rz,r, ex, ey, ez,i, iAtom, t_1, t1, t_2, t2, xyz, scrn1, dist, invdist, lsearch1, lsearch2)  
!$OMP Do
    do i = 1, nx
        ! Determine distance from each atom this point in space
        rx = x(i) - coord(1,:)
        ry = y(i) - coord(2,:)
        rz = z(i) - coord(3,:)
        r = sqrt(rx**2 + ry**2 + rz**2)

        call search(r, nc2, const, lsearch1)
        call search(r, nc2, atoms(1), lsearch2)

        if (lsearch1) then
          if (zeros) then
            efieldx(i) = ZERO
            efieldy(i) = ZERO
            efieldz(i) = ZERO
          else
            efieldx(i) = nan/nan ! 0.0/0.0
            efieldy(i) = nan/nan ! 0.0/0.0
            efieldz(i) = nan/nan ! 0.0/0.0
          end if
          cycle
        else if (scrn /=1 .and. scrn /= 2 .and. lsearch2 ) then 
          if (zeros) then
            efieldx(i) = ZERO
            efieldy(i) = ZERO
            efieldz(i) = ZERO
          else
            efieldx(i) = nan/nan ! 0.0/0.0
            efieldy(i) = nan/nan ! 0.0/0.0
            efieldz(i) = nan/nan ! 0.0/0.0
          end if
          cycle
        else 
          ex = ZERO
          ey = ZERO
          ez = ZERO
          if (scrn == 1) then
            do iAtom = 1, na
              xyz =  (/rx(iAtom),ry(iAtom),rz(iAtom) /)
              scrn1 = atoms(iAtom)
              dist = r(iAtom)
              invdist = ONE/dist
              call t2_screened(xyz, dist, invdist, scrn1, smear, t2)
              t_2 = cmplx(t2, ZERO, KINDR)

              ex = ex + t_2(1,1)*dipole(1, iAtom) & 
             &        + t_2(1,2)*dipole(2, iAtom) &
             &        + t_2(1,3)*dipole(3, iAtom)
              ey = ey + t_2(2,1)*dipole(1, iAtom) & 
             &        + t_2(2,2)*dipole(2, iAtom) &
             &        + t_2(2,3)*dipole(3, iAtom)
              ez = ez + t_2(3,1)*dipole(1, iAtom) & 
             &        + t_2(3,2)*dipole(2, iAtom) &
             &        + t_2(3,3)*dipole(3, iAtom)

              !Charge 
              call t1_screened(xyz, dist, invdist, scrn1, smear, t1)
              t_1 = cmplx(t1, ZERO, KINDR)
              ex = ex + t_1(1) * charge(iAtom) 
              ey = ey + t_1(2) * charge(iAtom) 
              ez = ez + t_1(3) * charge(iAtom) 
            end do ! iAtom

          else ! without screening
            do iAtom = 1, na
              xyz =  (/rx(iAtom),ry(iAtom),rz(iAtom) /)
              dist = r(iAtom)
              invdist = ONE/dist
              call t2_unscreened(xyz, invdist, t2)
              t_2 = cmplx(t2, ZERO, KINDR)

              ex = ex + t_2(1,1)*dipole(1, iAtom) & 
             &        + t_2(1,2)*dipole(2, iAtom) &
             &        + t_2(1,3)*dipole(3, iAtom)
              ey = ey + t_2(2,1)*dipole(1, iAtom) & 
             &        + t_2(2,2)*dipole(2, iAtom) &
             &        + t_2(2,3)*dipole(3, iAtom)
              ez = ez + t_2(3,1)*dipole(1, iAtom) & 
             &        + t_2(3,2)*dipole(2, iAtom) &
             &        + t_2(3,3)*dipole(3, iAtom)

              !Charge 
              call t1_unscreened(xyz, invdist, t1)
              t_1 = cmplx(t1, ZERO, KINDR)
              ex = ex + t_1(1) * charge(iAtom) 
              ey = ey + t_1(2) * charge(iAtom) 
              ez = ez + t_1(3) * charge(iAtom) 
            end do !iAtom
          end if

        ! Add in external field contribution
        if(dir .eq. 0) then
            ex = ex + 1
        else if(dir .eq. 1) then
            ey = ey + 1
        else if(dir .eq. 2) then        
            ez = ez + 1
        else if (dir .eq. 4) then
            continue
        end if
        ! Store the magnitude squared of the local field vector
        efieldx(i) = ex
        efieldy(i) = ey
        efieldz(i) = ez
      end if
    end do
!$OMP END DO
!$OMP END PARALLEL 
    return 
end subroutine drawing_math_vectors



subroutine drawing_math_ret(x, nx, y, ny, z, nz, coord, nc1, nc2, dipole, nd1, nd2, &
                atoms, na,efield, ne, charge, nc, zeros, dir, scrn, smear, freq, nsolv, avec)
    use Constants
  !$ use omp_lib
    implicit none
    ! Input variables
    integer,        intent(in) :: nx, ny, nz, nc1, nc2, nd1, nd2, na, ne, nc
!f2py intent(hide) :: nx, ny, nz, nc1, nc2, nd1, nd2, na, ne, nc
    real(kindr),    intent(in) :: charge(nc), x(nx), y(ny), z(nz)
    real(kindr),    intent(in) :: coord(nc1,nc2), atoms(na)
    complex(kindr), intent(in) :: dipole(nd1, nd2)
    integer,        intent(in) :: avec
    logical,        intent(in) :: zeros
    integer,        intent(in) :: dir
!f2py intent(in) :: dir
    integer,        intent(in) :: scrn 
    real(kindr),    intent(in) :: smear
    real(kindr),    intent(in) :: freq, nsolv
    complex(kindr)             :: t_1(3), t_2(3,3)
    Real(kindr)                :: t2(3,3)
    real(kindr)                :: xyz(3), scrn1,  dist, invdist

    ! Internal Variables
    real(kindr), parameter     :: const=0.001_KINDR
    real(kindr)                :: nan          
    real(kindr)    :: rx(nc2), ry(nc2), rz(nc2), r(nc2)
    complex(kindr) :: ex, ey, ez
    real(kindr)    :: kret 
    Real(kindr)    :: Xexp
    integer :: i, iAtom, iDir, idir2
    logical :: lsearch1, lsearch2
    real(KINDR) :: a_vec(3)
    
       
    ! Output variables
    complex(kindr), intent(inout) :: efield(ne)

    efield = ZERO_C

   !const=0.001_KINDR
    nan =0.0_KINDR

    kret = freq/ LIGHT_AU * nsolv

    a_vec = ZERO
    a_vec(avec+1) = ONE

!$OMP PARALLEL private(rx,ry,rz,r, ex,ey,ez,Xexp,i ,iAtom, t_1, t_2, xyz, scrn1, dist, invdist, lsearch1, lsearch2, iDir, iDir2 )  
!$OMP Do
    do i = 1,nx
      ! Determine distance from each atom this point in space
      rx = x(i) - coord(1,:)
      ry = y(i) - coord(2,:)
      rz = z(i) - coord(3,:)
      r = sqrt(rx**2 + ry**2 + rz**2)

      call search(r, nc2, const, lsearch1)
      call search(r, nc2, atoms(1), lsearch2)

      ! if grid point is too close to center of atom
      if (lsearch1) then
        if (zeros) then
          efield(i) = ZERO
        else
          efield(i) = nan/nan ! 0.0/0.0
        end if
        cycle
      ! if grid point is inside of the atom and there is no specified screening.
      else if (scrn /=1 .and. scrn /= 2 .and. lsearch2 ) then 
        if (zeros) then
          efield(i) = ZERO
        else
          efield(i) = nan/nan ! 0.0/0.0
        end if
        cycle
      else 
        ex = ZERO
        ey = ZERO
        ez = ZERO
        if (scrn == 1) then
          do iAtom = 1, na
            xyz =  (/rx(iAtom),ry(iAtom),rz(iAtom) /)
            scrn1 = atoms(iAtom)
            dist = r(iAtom)
            invdist = ONE/dist
            call t2_ret_screened(xyz, dist, invdist,scrn1,smear, kret, t_2)

            ex = ex + t_2(1,1)*dipole(1, iAtom) & 
           &        + t_2(1,2)*dipole(2, iAtom) &
           &        + t_2(1,3)*dipole(3, iAtom)
            ey = ey + t_2(2,1)*dipole(1, iAtom) & 
           &        + t_2(2,2)*dipole(2, iAtom) &
           &        + t_2(2,3)*dipole(3, iAtom)
            ez = ez + t_2(3,1)*dipole(1, iAtom) & 
           &        + t_2(3,2)*dipole(2, iAtom) &
           &        + t_2(3,3)*dipole(3, iAtom)

            !Charge. DIM does not currently support the calculation of charge with retardation
        !   call  t1_ret_screened(xyz, dist, invdist, scrn1, smear, kret, t_1)
        !   ex = ex + t_1(1) * charge(iAtom) 
        !   ey = ey + t_1(2) * charge(iAtom) 
        !   ez = ez + t_1(3) * charge(iAtom) 
          end do ! iAtom

        else ! without screeening
          do iAtom = 1, na
            xyz =  (/rx(iAtom),ry(iAtom),rz(iAtom) /)
            dist = r(iAtom)
            invdist = ONE/dist
            call t2_ret(xyz, dist, invdist, kret, t_2)

            ex = ex + t_2(1,1)*dipole(1, iAtom) & 
           &        + t_2(1,2)*dipole(2, iAtom) &
           &        + t_2(1,3)*dipole(3, iAtom)
            ey = ey + t_2(2,1)*dipole(1, iAtom) & 
           &        + t_2(2,2)*dipole(2, iAtom) &
           &        + t_2(2,3)*dipole(3, iAtom)
            ez = ez + t_2(3,1)*dipole(1, iAtom) & 
           &        + t_2(3,2)*dipole(2, iAtom) &
           &        + t_2(3,3)*dipole(3, iAtom)

            !Charge. DIM does not currently support the calculation of charge with retardation
         !  call t1_ret(xyz, dist, invdist, kret, t_1)
         !  ex = ex + t_1(1) * charge(iAtom) 
         !  ey = ey + t_1(2) * charge(iAtom) 
         !  ez = ez + t_1(3) * charge(iAtom) 
          end do !iAtom
        end if

        Xexp = ZERO
        Xexp = Xexp + a_vec(1)*kret*x(i)
        Xexp = Xexp + a_vec(2)*kret*y(i)
        Xexp = Xexp + a_vec(3)*kret*z(i)
        if(dir .eq. 0) then
            ex = ex + exp(I_C * Xexp)
        else if(dir .eq. 1) then
            ey = ey + exp(I_C * Xexp)
        else if(dir .eq. 2) then        
            ez = ez + exp(I_C * Xexp)
        else if (dir .eq. 4) then
            continue
        end if
        ! Store the magnitude squared of the local field vector
        !efield(i) = (ez*conjg(ez))
        efield(i) = (ex*conjg(ex) + ey*conjg(ey) + ez*conjg(ez))
      end if
    end do
!$OMP end do 
!$OMP end PARALLEL
    return 

end subroutine drawing_math_ret

subroutine drawing_math_vectors_ret(x, nx, y, ny, z, nz, coord, nc1, nc2, dipole, nd1, nd2, &
                atoms, na, efieldx, nex, efieldy, ney, efieldz, nez, charge, nc, zeros, &
                dir, scrn, smear, freq, nsolv, avec)
    use Constants
  !$ use omp_lib
    implicit none
    ! Input variables
    integer,        intent(in) :: nx, ny, nz, nc1, nc2, nd1, nd2, na, nex, ney, nez, nc
!f2py integer, intent(hide) :: nx, ny, nz, nc1, nc2, nd1, nd2, na, nex, ney, nez, nc
    real(kindr),    intent(in) :: charge(nc), x(nx), y(ny), z(nz)
    real(kindr),    intent(in) :: coord(nc1,nc2), atoms(na)
    complex(kindr), intent(in) :: dipole(nd1, nd2)
    logical,        intent(in) :: zeros
    integer,        intent(in) :: avec
    integer, optional, intent(in)  :: dir
!f2py integer, optional, intent(in) :: dir
    integer,        intent(in) :: scrn
    real(kindr),    intent(in) :: smear
    real(kindr),    intent(in) :: freq, nsolv
    complex(kindr)             :: t_1(3), t_2(3,3)
    real(kindr)                :: xyz(3), scrn1,  dist, invdist
    
    ! Internal Variables
    real(kindr), parameter     :: const=0.001_KINDR
    real(kindr)                :: nan          

    real(kindr)    :: rx(nc2), ry(nc2), rz(nc2), r(nc2)
    complex(kindr) :: ex, ey, ez
    real(kindr)    :: kret
    integer :: i, iAtom
    Real(kindr)    ::  Xexp
    logical :: lsearch1, lsearch2
    real(KINDR):: a_vec(3)
    
    ! Output variables
    complex(kindr) :: efieldx(nex), efieldy(ney), efieldz(nez)

    ! zero out the outputs
    efieldx = ZERO_C
    efieldy = ZERO_C
    efieldz = ZERO_C

    kret = freq/ LIGHT_AU * nsolv

   !const = 0.001_KINDR
    nan = 0.0_KINDR

    a_vec = ZERO
    a_vec(avec+1) = ONE

    !Math!
!$OMP PARALLEL private(rx,ry,rz,r, ex,ey,ez,Xexp,i,iAtom, t_1, t_2, xyz, scrn1, dist, invdist, lsearch1, lsearch2 )  
!$OMP Do
    do i = 1, nx
        rx = x(i) - coord(1,:)
        ry = y(i) - coord(2,:)
        rz = z(i) - coord(3,:)
        r = sqrt(rx**2 + ry**2 + rz**2)

        call search(r, nc2, const, lsearch1)
        call search(r, nc2, atoms(1), lsearch2)

        if (lsearch1) then
          if (zeros) then
            efieldx(i) = ZERO
            efieldy(i) = ZERO
            efieldz(i) = ZERO
          else
            efieldx(i) = nan/nan ! 0.0/0.0
            efieldy(i) = nan/nan ! 0.0/0.0
            efieldz(i) = nan/nan ! 0.0/0.0
          end if
          cycle
        else if (scrn /=1 .and. scrn /= 2 .and. lsearch2 ) then 
          if (zeros) then
            efieldx(i) = ZERO
            efieldy(i) = ZERO
            efieldz(i) = ZERO
          else
            efieldx(i) = nan/nan ! 0.0/0.0
            efieldy(i) = nan/nan ! 0.0/0.0
            efieldz(i) = nan/nan ! 0.0/0.0
          end if
          cycle
        else 
          ex = ZERO
          ey = ZERO
          ez = ZERO
          if (scrn == 1) then
            do iAtom = 1, na
              xyz =  (/rx(iAtom),ry(iAtom),rz(iAtom) /)
              scrn1 = atoms(iAtom)
              dist = r(iAtom)
              invdist = ONE/dist
              call t2_ret_screened(xyz, dist, invdist, scrn1, smear, kret, t_2)

              ex = ex + t_2(1,1)*dipole(1, iAtom) & 
             &        + t_2(1,2)*dipole(2, iAtom) &
             &        + t_2(1,3)*dipole(3, iAtom)
              ey = ey + t_2(2,1)*dipole(1, iAtom) & 
             &        + t_2(2,2)*dipole(2, iAtom) &
             &        + t_2(2,3)*dipole(3, iAtom)
              ez = ez + t_2(3,1)*dipole(1, iAtom) & 
             &        + t_2(3,2)*dipole(2, iAtom) &
             &        + t_2(3,3)*dipole(3, iAtom)

            !Charge. DIM does not currently support the calculation of charge with retardation
           !  call  t1_ret_screened(xyz, dist, invdist, scrn1, smear, kret, t_1)
           !  ex = ex + t_1(1) * charge(iAtom) 
           !  ey = ey + t_1(2) * charge(iAtom) 
           !  ez = ez + t_1(3) * charge(iAtom) 
            end do ! iAtom

          else  ! without screening
            do iAtom = 1, na
              xyz =  (/rx(iAtom),ry(iAtom),rz(iAtom) /)
              dist = r(iAtom)
              invdist = ONE/dist
              call t2_ret(xyz, dist, invdist, kret, t_2)

              ex = ex + t_2(1,1)*dipole(1, iAtom) & 
             &        + t_2(1,2)*dipole(2, iAtom) &
             &        + t_2(1,3)*dipole(3, iAtom)
              ey = ey + t_2(2,1)*dipole(1, iAtom) & 
             &        + t_2(2,2)*dipole(2, iAtom) &
             &        + t_2(2,3)*dipole(3, iAtom)
              ez = ez + t_2(3,1)*dipole(1, iAtom) & 
             &        + t_2(3,2)*dipole(2, iAtom) &
             &        + t_2(3,3)*dipole(3, iAtom)

            !Charge. DIM does not currently support the calculation of charge with retardation
            ! call t1_ret(xyz, dist, invdist, kret, t_1)
            ! ex = ex + t_1(1) * charge(iAtom) 
            ! ey = ey + t_1(2) * charge(iAtom) 
            ! ez = ez + t_1(3) * charge(iAtom) 
            end do !iAtom
          end if
        if(present(dir)) then
          Xexp = ZERO
          Xexp = Xexp + a_vec(1)*kret*x(i)
          Xexp = Xexp + a_vec(2)*kret*y(i)
          Xexp = Xexp + a_vec(3)*kret*z(i)

          if(dir .eq. 0) then
              ex = ex + exp(I_C * Xexp)
          else if(dir .eq. 1) then
              ey = ey + exp(I_C * Xexp)
          else if(dir .eq. 2) then        
              ez = ez + exp(I_C * Xexp)
          else if (dir .eq. 4) then
              continue
          end if
        end if
        efieldx(i) =  ex
        efieldy(i) =  ey
        efieldz(i) =  ez
      end if
    end do
!$OMP END DO
!$OMP END PARALLEL
    return 

end subroutine drawing_math_vectors_ret




subroutine drawing_math_ret_magnetic(x, nx, y, ny, z, nz, coord, nc1, nc2, dipole, nd1, nd2, &
                atoms, na,hfield, ne, charge, nc, zeros, dir, scrn, smear, freq, nsolv, avec)
    use Constants
  !$ use omp_lib
    implicit none
    ! Input variables
    integer,        intent(in) :: nx, ny, nz, nc1, nc2, nd1, nd2, na, ne, nc
!f2py intent(hide) :: nx, ny, nz, nc1, nc2, nd1, nd2, na, ne, nc
    real(kindr),    intent(in) :: charge(nc), x(nx), y(ny), z(nz)
    real(kindr),    intent(in) :: coord(nc1,nc2), atoms(na)
    complex(kindr), intent(in) :: dipole(nd1, nd2)
    logical,        intent(in) :: zeros
    integer,        intent(in) :: avec
    integer,        intent(in) :: dir
!f2py intent(in) :: dir
    integer,        intent(in) :: scrn 
    real(kindr),    intent(in) :: smear
    real(kindr),    intent(in) :: freq, nsolv

    complex(kindr)             :: t_1(3)
    real(kindr)                :: xyz(3), scrn1,  dist, invdist

    ! Internal Variables
    real(kindr)    :: rx(nc2), ry(nc2), rz(nc2), r(nc2)
    complex(kindr) :: hx, hy, hz
    real(kindr)    :: kret
    Real(kindr)    :: a_vec(3), Xexp, avecscale
    integer :: i, iAtom
    logical :: lsearch1, lsearch2
    real(kindr), parameter     :: const=0.001_KINDR
    real(kindr)                :: nan          
    
       
    ! Output variables
    complex(kindr), intent(inout) :: hfield(ne)

    ! Zero out output
    hfield = ZERO_C

   !const=0.001_KINDR
    nan = ZERO
    kret = freq/ LIGHT_AU * nsolv

    a_vec = ZERO
    a_vec(avec+1) = ONE
    avecscale = sqrt(a_vec(1)**2+a_vec(2)**2+a_vec(3)**2)


!$OMP PARALLEL private(rx,ry,rz,r, hx,hy,hz, Xexp,iAtom, t_1, xyz, scrn1, dist, invdist, lsearch1, lsearch2)  
!$OMP Do
    do i = 1,nx 
      ! Determine distance from each atom this point in space
      rx = x(i) - coord(1,:)
      ry = y(i) - coord(2,:)
      rz = z(i) - coord(3,:)
      r = sqrt(rx**2 + ry**2 + rz**2)

      call search(r, nc2, const, lsearch1)
      call search(r, nc2, atoms(1), lsearch2)

      if (lsearch1) then
        if (zeros) then
          hfield(i) = ZERO
        else
          hfield(i) = nan/nan !0.0/0.0
        end if
        cycle
      
      else if (scrn /=1 .and. scrn /= 2 .and. lsearch2 ) then 
        if (zeros) then
          hfield(i) = ZERO
        else
          hfield(i) = nan/nan !0.0/0.0
        end if
        cycle

      else
        hx = ZERO
        hy = ZERO
        hz = ZERO
        if (scrn == 1) then
          do iAtom = 1, na
            xyz =  (/rx(iAtom),ry(iAtom),rz(iAtom) /)
            scrn1 = atoms(iAtom)
            dist = r(iAtom)
            invdist = ONE/dist
            call t1_ret_screened(xyz, dist, invdist, scrn1, smear, kret, t_1)
            
            hx = hx + t_1(2)*dipole(3, iAtom) - t_1(3)* dipole(2,iAtom)
            hy = hy + t_1(3)*dipole(1, iAtom) - t_1(1)* dipole(3,iAtom)
            hz = hz + t_1(1)*dipole(2, iAtom) - t_1(2)* dipole(1,iAtom)
          end do ! iAtom

        else ! withoug screening
          do iAtom = 1, na
            xyz =  (/rx(iAtom),ry(iAtom),rz(iAtom) /)
            scrn1 = atoms(iAtom)
            dist = r(iAtom)
            invdist = ONE/dist
            call t1_ret(xyz, dist, invdist, kret, t_1)

            hx = hx + t_1(2)*dipole(3, iAtom) - t_1(3)* dipole(2,iAtom)
            hy = hy + t_1(3)*dipole(1, iAtom) - t_1(1)* dipole(3,iAtom)
            hz = hz + t_1(1)*dipole(2, iAtom) - t_1(2)* dipole(1,iAtom)
          end do ! iAtom
        end if

        ! add factor of i*k/c to field to match the definition
        hx = I_C*kret*hx/LIGHT_AU 
        hy = I_C*kret*hy/LIGHT_AU
        hz = I_C*kret*hz/LIGHT_AU

        Xexp = ZERO
        Xexp = Xexp + a_vec(1)*kret*x(i)
        Xexp = Xexp + a_vec(2)*kret*y(i)
        Xexp = Xexp + a_vec(3)*kret*z(i)
  
        if(dir .eq. 0) then
            hy = hy + (a_vec(3)/avecscale)*exp(I_C * Xexp)/LIGHT_AU
            hz = hz - (a_vec(2)/avecscale)*exp(I_C * Xexp)/LIGHT_AU
        else if(dir .eq. 1) then
            hx = hx - (a_vec(3)/avecscale)*exp(I_C * Xexp)/LIGHT_AU
            hz = hz + (a_vec(1)/avecscale)*exp(I_C * Xexp)/LIGHT_AU
        else if(dir .eq. 2) then        
            hx = hx + (a_vec(2)/avecscale)*exp(I_C * Xexp)/LIGHT_AU
            hy = hy - (a_vec(1)/avecscale)*exp(I_C * Xexp)/LIGHT_AU
        else if (dir .eq. 4) then
            continue
        end if

        ! Store the magnitude squared of the local field vector
        !hfield(i) = (hz*conjg(hz))
        hfield(i) = (hx*conjg(hx) + hy*conjg(hy) + hz*conjg(hz))
      end if 
    end do ! i
!$OMP end do 
!$OMP end PARALLEL
    return 

end subroutine drawing_math_ret_magnetic

subroutine drawing_math_vectors_ret_magnetic(x, nx, y, ny, z, nz, coord, nc1, nc2, dipole, nd1, nd2, &
                atoms, na, hfieldx, nex, hfieldy, ney, hfieldz, nez, charge, nc, zeros, &
                dir, scrn, smear, freq, nsolv, avec)
    use Constants
  !$ use omp_lib
    implicit none
    ! Input variables
    integer,        intent(in) :: nx, ny, nz, nc1, nc2, nd1, nd2, na, nex, ney, nez, nc
!f2py integer, intent(hide) :: nx, ny, nz, nc1, nc2, nd1, nd2, na, nex, ney, nez, nc
    real(kindr),    intent(in) :: charge(nc), x(nx), y(ny), z(nz)
    real(kindr),    intent(in) :: coord(nc1,nc2), atoms(na)
    complex(kindr), intent(in) :: dipole(nd1, nd2)
    logical,        intent(in) :: zeros
    integer,        intent(in) :: avec
    integer, optional, intent(in)  :: dir
!f2py integer, optional, intent(in) :: dir
    integer,        intent(in) :: scrn
    real(kindr),    intent(in) :: smear
    real(kindr),    intent(in) :: freq, nsolv

    complex(kindr)             :: t_1(3)
    real(kindr)                :: xyz(3), scrn1,  dist, invdist
    
    ! Internal Variables
    real(kindr), parameter     :: const=0.001_KINDR
    real(kindr)                :: nan          
    real(kindr)    :: rx(nc2), ry(nc2), rz(nc2), r(nc2)
    complex(kindr) :: hx, hy, hz
    real(kindr)    :: kret
    integer :: i, iAtom
    Real(kindr)    :: a_vec(3),Xexp, avecscale
    logical :: lsearch1, lsearch2
    
    ! Output variables
    complex(kindr) :: hfieldx(nex), hfieldy(ney), hfieldz(nez)

    ! Zero out the outputs 
    hfieldx = ZERO_C
    hfieldy = ZERO_C
    hfieldz = ZERO_C

   !const = 0.001_KINDR
    nan = 0.0_KINDR
    kret = freq/ LIGHT_AU * nsolv

    a_vec = ZERO
    a_vec(avec+1) = ONE
    avecscale = sqrt(a_vec(1)**2+a_vec(2)**2+a_vec(3)**2)
    
!$OMP PARALLEL private(rx,ry,rz,r, hx,hy,hz, Xexp,iAtom, t_1, xyz, scrn1, dist, invdist, lsearch1, lsearch2)  
!$OMP Do
    do i = 1, nx
      rx = x(i) - coord(1,:)
      ry = y(i) - coord(2,:)
      rz = z(i) - coord(3,:)
      r = sqrt(rx**2 + ry**2 + rz**2)

      call search(r, nc2, const, lsearch1)
      call search(r, nc2, atoms(1), lsearch2)

      if (lsearch1) then
        if (zeros) then
          hfieldx(i) = ZERO
          hfieldy(i) = ZERO
          hfieldz(i) = ZERO
        else
          hfieldx(i) = nan/nan !0.0/0.0
          hfieldy(i) = nan/nan !0.0/0.0
          hfieldz(i) = nan/nan !0.0/0.0
        end if
        cycle
      
      else if (scrn /=1 .and. scrn /= 2 .and. lsearch2 ) then 
        if (zeros) then
          hfieldx(i) = ZERO
          hfieldy(i) = ZERO
          hfieldz(i) = ZERO
        else
          hfieldx(i) = nan/nan !0.0/0.0
          hfieldy(i) = nan/nan !0.0/0.0
          hfieldz(i) = nan/nan !0.0/0.0
        end if
        cycle
  
      else
        hx = ZERO
        hy = ZERO
        hz = ZERO
        if (scrn == 1) then
          do iAtom = 1, na
            xyz =  (/rx(iAtom),ry(iAtom),rz(iAtom) /)
            scrn1 = atoms(iAtom)
            dist = r(iAtom)
            invdist = ONE/dist
            call t1_ret_screened(xyz, dist, invdist, scrn1, smear, kret, t_1)

            hx = hx + t_1(2)*dipole(3, iAtom) - t_1(3)* dipole(2,iAtom)
            hy = hy + t_1(3)*dipole(1, iAtom) - t_1(1)* dipole(3,iAtom)
            hz = hz + t_1(1)*dipole(2, iAtom) - t_1(2)* dipole(1,iAtom)
          end do ! iAtom
  
        else  ! without screening
          do iAtom = 1, na
            xyz =  (/rx(iAtom),ry(iAtom),rz(iAtom) /)
            scrn1 = atoms(iAtom)
            dist = r(iAtom)
            invdist = ONE/dist
            call t1_ret(xyz, dist, invdist, kret, t_1)

            hx = hx + t_1(2)*dipole(3, iAtom) - t_1(3)* dipole(2,iAtom)
            hy = hy + t_1(3)*dipole(1, iAtom) - t_1(1)* dipole(3,iAtom)
            hz = hz + t_1(1)*dipole(2, iAtom) - t_1(2)* dipole(1,iAtom)
          end do ! iAtom
        end if
  
        ! add factor of i*k/c to field to match the definition
        hx = I_C*kret*hx/LIGHT_AU 
        hy = I_C*kret*hy/LIGHT_AU
        hz = I_C*kret*hz/LIGHT_AU
  
        Xexp = ZERO
        Xexp = Xexp + a_vec(1)*kret*x(i)
        Xexp = Xexp + a_vec(2)*kret*y(i)
        Xexp = Xexp + a_vec(3)*kret*z(i)
    
        if(dir .eq. 0) then
            hy = hy + (a_vec(3)/avecscale)*exp(I_C * Xexp)/LIGHT_AU
            hz = hz - (a_vec(2)/avecscale)*exp(I_C * Xexp)/LIGHT_AU
        else if(dir .eq. 1) then
            hx = hx - (a_vec(3)/avecscale)*exp(I_C * Xexp)/LIGHT_AU
            hz = hz + (a_vec(1)/avecscale)*exp(I_C * Xexp)/LIGHT_AU
        else if(dir .eq. 2) then        
            hx = hx + (a_vec(2)/avecscale)*exp(I_C * Xexp)/LIGHT_AU
            hy = hy - (a_vec(1)/avecscale)*exp(I_C * Xexp)/LIGHT_AU
        else if (dir .eq. 4) then
            continue
        end if

        hfieldx(i) =  hx
        hfieldy(i) =  hy
        hfieldz(i) =  hz
      end if
    end do
!$OMP END DO
!$OMP END PARALLEL
    return 


end subroutine drawing_math_vectors_ret_magnetic

subroutine drawing_math(x, nx, y, ny, z, nz, coord, nc1, nc2, dipole, nd1, nd2, &
                atoms, na, efield, ne, charge, nc, zeros, dir, scrn, smear)
    use Constants
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
    real            const   

    ! Internal Variables
    real(kindr)    :: rx(nc2), ry(nc2), rz(nc2), r(nc2)
    real(kindr)    :: screen(nc2), sf1(nc2), sf2(nc2)
    real(kindr)    :: Rpp(na), invRpp(na)
    complex(kindr) :: dot(nc2), ex, ey, ez
    integer :: i
    
       
    ! Output variables
    complex(kindr) :: efield(ne)
    ! Loop of all points in space
!$OMP PARALLEL
    const=0.001
    do i = 1, nx
        ! Determine distance from each atom this point in space
        rx = x(i) - coord(1,:)
        ry = y(i) - coord(2,:)
        rz = z(i) - coord(3,:)
        r = sqrt(rx**2 + ry**2 + rz**2)
        if (searchr(r, nc2, const)) then
            efield(i) = ZERO
            r = 1e8
        end if
        screen = ONE
        if(scrn == 1) then
          Rpp = sqrt(atoms**2  + smear**2)
          invRpp = ONE / Rpp
          sf1 = (TWO * r * invRpp / SQRTPI) * EXP(-((r * invRpp)**2))
          sf2 = (FOUR * invRpp**3 / SQRTPI) * EXP(-((r * invRpp)**2))
          screen = ERF(r*invRpp) - sf1
        else if(scrn == 2) then
         screen = ONE
        else
        ! If a distance is less than radius of the atom we are currently on
        ! set the distance to NaN by dividing by zero.
        ! We only do this if there's no screening.
          if (search(r, nc2, atoms(1))) then
            if (zeros) then
              efield(i) = 0.0
            else
              efield(i) = 0.0/0.0
            end if
            cycle
          end if
        end if 
           
        ! Dot product of dipole moments and distance
        dot = dipole(1,:)*rx + dipole(2,:)*ry + dipole(3,:)*rz
        ! Dipole contribution to induced field
        ex = sum(screen*((3*dot)*rx/(r**5) - dipole(1,:)/(r**3)))
        ey = sum(screen*((3*dot)*ry/(r**5) - dipole(2,:)/(r**3)))
        ez = sum(screen*((3*dot)*rz/(r**5) - dipole(3,:)/(r**3)))
        if(scrn == 1) then
          ex = ex - sum(sf2*dot*rx/(r**2))
          ey = ey - sum(sf2*dot*ry/(r**2))
          ez = ez - sum(sf2*dot*rz/(r**2))
        end if
        ! Monopole contribution to induced field
        ex = ex + sum(screen*(charge(:)*rx)/(r**3))
        ey = ey + sum(screen*(charge(:)*ry)/(r**3))
        ez = ez + sum(screen*(charge(:)*rz)/(r**3))
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
    end do
    return 
!$OMP END PARALLEL
contains

logical function search(list, nl, test)
    integer,     intent(in) :: nl
    real(kindr), intent(in) :: list(nl)
    real(kindr), intent(in) :: test
    integer                 :: i
    search = .false.
    do i = 1, nl
        if (list(i) <= test) then
            search = .true.
            exit
        end if
    end do
end function search

logical function searchr(list, nl, test)
    integer,     intent(in) :: nl
    real(kindr), intent(in) :: list(nl)
    real                        test
    integer                 :: i
    searchr = .false.
    do i = 1, nl
        if (list(i) <= test) then
            searchr = .true.
            exit
        end if
    end do
end function searchr

end subroutine drawing_math

subroutine drawing_math_vectors(x, nx, y, ny, z, nz, coord, nc1, nc2, dipole, nd1, nd2, &
                atoms, na, efieldx, nex, efieldy, ney, efieldz, nez, &
                charge, nc, zeros, dir, scrn, smear)
    use Constants
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
    real            const   

    ! Internal Variables
    real(kindr)    :: rx(nc2), ry(nc2), rz(nc2), r(nc2)
    real(kindr)    :: screen(nc2), sf1(nc2), sf2(nc2)
    real(kindr)    :: Rpp(na), invRpp(na)
    complex(kindr) :: dot(nc2), ex, ey, ez
    integer :: i
    
       
    ! Output variables
    complex(kindr) :: efieldx(nex), efieldy(ney), efieldz(nez)
    ! Loop of all points in space
!$OMP PARALLEL
    const=0.001
    do i = 1, nx
        ! Determine distance from each atom this point in space
        rx = x(i) - coord(1,:)
        ry = y(i) - coord(2,:)
        rz = z(i) - coord(3,:)
        r = sqrt(rx**2 + ry**2 + rz**2)
        if (searchr(r, nc2, const)) then
            efieldx(i) = ZERO
            efieldy(i) = ZERO
            efieldz(i) = ZERO
            r = 1e8
        end if
        screen = ONE
        if(scrn == 1) then
          Rpp = sqrt(atoms**2  + smear**2)
          invRpp = ONE / Rpp
          sf1 = (TWO * r * invRpp / SQRTPI) * EXP(-((r * invRpp)**2))
          sf2 = (FOUR * invRpp**3 / SQRTPI) * EXP(-((r * invRpp)**2))
          screen = ERF(r*invRpp) - sf1
        else if(scrn == 2) then
         screen = ONE
        else
        ! If a distance is less than radius of the atom we are currently on
        ! set the distance to NaN by dividing by zero.
        ! We only do this if there's no screening.
          if (search(r, nc2, atoms(1))) then
            if (zeros) then
              efieldx(i) = 0.0
              efieldy(i) = 0.0
              efieldz(i) = 0.0
            else
              efieldx(i) = 0.0/0.0
              efieldy(i) = 0.0/0.0
              efieldz(i) = 0.0/0.0
            end if
            cycle
          end if
        end if 
           
        ! Dot product of dipole moments and distance
        dot = dipole(1,:)*rx + dipole(2,:)*ry + dipole(3,:)*rz
        ! Dipole contribution to induced field
        ex = sum(screen*((3*dot)*rx/(r**5) - dipole(1,:)/(r**3)))
        ey = sum(screen*((3*dot)*ry/(r**5) - dipole(2,:)/(r**3)))
        ez = sum(screen*((3*dot)*rz/(r**5) - dipole(3,:)/(r**3)))
        if(scrn == 1) then
          ex = ex - sum(sf2*dot*rx/(r**2))
          ey = ey - sum(sf2*dot*ry/(r**2))
          ez = ez - sum(sf2*dot*rz/(r**2))
        end if
        ! Monopole contribution to induced field
        ex = ex + sum(screen*(charge(:)*rx)/(r**3))
        ey = ey + sum(screen*(charge(:)*ry)/(r**3))
        ez = ez + sum(screen*(charge(:)*rz)/(r**3))
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
        efieldx(i) = ex
        efieldy(i) = ey
        efieldz(i) = ez
    end do
    return 
!$OMP END PARALLEL
contains

logical function search(list, nl, test)
    integer,     intent(in) :: nl
    real(kindr), intent(in) :: list(nl)
    real(kindr), intent(in) :: test
    integer                 :: i
    search = .false.
    do i = 1, nl
        if (list(i) <= test) then
            search = .true.
            exit
        end if
    end do
end function search

logical function searchr(list, nl, test)
    integer,     intent(in) :: nl
    real(kindr), intent(in) :: list(nl)
    real                        test
    integer                 :: i
    searchr = .false.
    do i = 1, nl
        if (list(i) <= test) then
            searchr = .true.
            exit
        end if
    end do
end function searchr

end subroutine drawing_math_vectors


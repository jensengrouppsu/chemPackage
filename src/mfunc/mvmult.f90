
! LEVEL 1 BLAS : VECTOR-VECTOR

Subroutine vecvecr (x, y, z, n)

   use Constants
   
   Implicit None

   Integer,     Intent(In)  :: n
   Real(KINDR), Intent(In)  :: x(n)
   Real(KINDR), Intent(In)  :: y(n)
   Real(KINDR), Intent(Out) :: z

   Real(KINDR), External    :: ddot

   z = ddot (n, x, 1, y, 1)

End Subroutine vecvecr

Subroutine vecvecc (x, y, z, n)

   use Constants

   Implicit None

   Integer,        Intent(In)  :: n
   Complex(KINDR), Intent(In)  :: x(n)
   Complex(KINDR), Intent(In)  :: y(n)
   Complex(KINDR), Intent(Out) :: z

   Complex(KINDR), External    :: zdotc

   z = zdotc (n, x, 1, y, 1)

End Subroutine vecvecc

Subroutine normr (x, z, n)

   use Constants

   Implicit None

   Integer,      Intent(In)  :: n
   Real(KINDR),  Intent(In)  :: x(n)
   Real(KINDR),  Intent(out) :: z

   Real(KINDR),  External    :: dnrm2

   z = dnrm2(n, x, 1)

ENd Subroutine normr

Subroutine normc (x, z, n)

   use Constants

   Implicit None

   Integer,        Intent(In)  :: n
   Complex(KINDR), Intent(In)  :: x(n)
   Real(KINDR),    Intent(out) :: z

   Real(KINDR),    External    :: dznrm2

   z = dznrm2(n, x, 1)

ENd Subroutine normc

! LEVEL 2 BLAS : MATRIX-VECTOR

Subroutine vecmatr (x, a, y, nr, nc)

   use Constants

   Implicit None

   Integer,     Intent(In)  :: nr
   Integer,     Intent(In)  :: nc
   Real(KINDR), Intent(In)  :: a(nr,nc)
   Real(KINDR), Intent(In)  :: x(nr)
   Real(KINDR), Intent(Out) :: y(nc)

   call dgemv ('T', nr, nc, ONE, a, nr, x, 1, ZERO, y, 1)

End Subroutine vecmatr

Subroutine vecmatc (x, a, y, nr, nc)

   use Constants

   Implicit None

   Integer,        Intent(In)  :: nr
   Integer,        Intent(In)  :: nc
   Complex(KINDR), Intent(In)  :: a(nr,nc)
   Complex(KINDR), Intent(In)  :: x(nr)
   Complex(KINDR), Intent(Out) :: y(nc)

   call zgemv ('T', nr, nc, ONE_C, a, nr, x, 1, ZERO_C, y, 1)

End Subroutine vecmatc

Subroutine matvecr (a, x, y, nr, nc)

   use Constants

   Implicit None

   Integer,      Intent(In)  :: nr
   Integer,      Intent(In)  :: nc
   Real(KINDR),  Intent(In)  :: a(nr,nc)
   Real(KINDR),  Intent(In)  :: x(nc)
   Real(KINDR),  Intent(Out) :: y(nr)

   call dgemv ('N', nr, nc, ONE, a, nr, x, 1, ZERO, y, 1)

End Subroutine matvecr

Subroutine matvecc (a, x, y, nr, nc)

   use Constants

   Implicit None

   Integer,        Intent(In)  :: nr
   Integer,        Intent(In)  :: nc
   Complex(KINDR), Intent(In)  :: a(nr,nc)
   Complex(KINDR), Intent(In)  :: x(nc)
   Complex(KINDR), Intent(Out) :: y(nr)

   call zgemv ('N', nr, nc, ONE_C, a, nr, x, 1, ZERO_C, y, 1)

End Subroutine matvecc

! LEVEL 3 BLAS : MATRIX-MATRIX

Subroutine matmatr (a, b, c, nra, nca, nrb, ncb)

   use Constants

   Implicit None

   Integer,      Intent(In)  :: nra
   Integer,      Intent(In)  :: nca
   Integer,      Intent(In)  :: nrb
   Integer,      Intent(In)  :: ncb
   Real(KINDR),  Intent(In)  :: a(nra,nca)
   Real(KINDR),  Intent(In)  :: b(nrb,ncb)
   Real(KINDR),  Intent(Out) :: c(nra,ncb)

   call dgemm ('N', 'N', nra, ncb, nca, ONE, a, nra, b, nrb, ZERO, c, nra)

End Subroutine matmatr

Subroutine matmatc (a, b, c, nra, nca, nrb, ncb)

   use Constants

   Implicit None

   Integer,        Intent(In)  :: nra
   Integer,        Intent(In)  :: nca
   Integer,        Intent(In)  :: nrb
   Integer,        Intent(In)  :: ncb
   Complex(KINDR), Intent(In)  :: a(nra,nca)
   Complex(KINDR), Intent(In)  :: b(nrb,ncb)
   Complex(KINDR), Intent(Out) :: c(nra,ncb)

   call zgemm ('N', 'N', nra, ncb, nca, ONE_C, a, nra, b, nrb, ZERO_C, c, nra)

End Subroutine matmatc

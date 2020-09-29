Module Constants

   Implicit None

!  Precision
   Integer,        Parameter :: KINDR = kind(0d0)

   Real(KINDR),    Parameter :: PI          = 3.14159265358979323846_KINDR
   REAL(KINDR),    Parameter :: SQRTPI      = 1.7724538509055159_KINDR
   Real(KINDR),    Parameter :: HBAR        = 6.5821189916e-16_KINDR !units = eV s
   Real(KINDR),    Parameter :: nm2eV       = 1239.0_KINDR

   Real(KINDR),    Parameter :: ZERO        =  0.0e0_KINDR
   Real(KINDR),    Parameter :: ONE         =  1.0e0_KINDR
   Real(KINDR),    Parameter :: TWO         =  2.0e0_KINDR
   Real(KINDR),    Parameter :: THREE       =  3.0e0_KINDR
   Real(KINDR),    Parameter :: FOUR        =  4.0e0_KINDR
   Real(KINDR),    Parameter :: FIVE        =  5.0e0_KINDR
   Real(KINDR),    Parameter :: SIX         =  6.0e0_KINDR
   Real(KINDR),    Parameter :: SEVEN       =  7.0e0_KINDR
   Real(KINDR),    Parameter :: EIGHT       =  8.0e0_KINDR
   Real(KINDR),    Parameter :: NINE        =  9.0e0_KINDR
   Real(KINDR),    Parameter :: TEN         = 10.0e0_KINDR
   Real(KINDR),    Parameter :: ELEVEN      = 11.0e0_KINDR
   Real(KINDR),    Parameter :: TWELVE      = 12.0e0_KINDR
   Real(KINDR),    Parameter :: ONEFOUR     = 14.0e0_KINDR
   Real(KINDR),    Parameter :: ONEFIVE     = 15.0e0_KINDR
   Real(KINDR),    Parameter :: SIXTEEN     = 16.0e0_KINDR
   Real(KINDR),    Parameter :: ONESEVEN    = 17.0e0_KINDR
   Real(KINDR),    Parameter :: EIGHTEEN    = 18.0e0_KINDR
   Real(KINDR),    Parameter :: ONENINE     = 19.0e0_KINDR
   Real(KINDR),    Parameter :: TWENTY      = 20.0e0_KINDR
   Real(KINDR),    Parameter :: TWOTWO      = 22.0e0_KINDR
   Real(KINDR),    Parameter :: TWOFOUR     = 24.0e0_KINDR
   Real(KINDR),    Parameter :: THREETWO    = 32.0e0_KINDR
   Real(KINDR),    Parameter :: THREESIX    = 36.0e0_KINDR
   Real(KINDR),    Parameter :: FOUREIGHT   = 48.0e0_KINDR
   Real(KINDR),    Parameter :: SEVENTWO    = 72.0e0_KINDR
   Real(KINDR),    Parameter :: FOURTH      = 0.25e0_KINDR
   Real(KINDR),    Parameter :: HALF        =  0.5e0_KINDR
   Real(KINDR),    Parameter :: THIRD       =  1.0e0_KINDR / 3.0e0_KINDR
   Real(KINDR),    Parameter :: TWOTHIRD    =  2.0e0_KINDR / 3.0e0_KINDR
   Real(KINDR),    Parameter :: SIXTH       =  1.0e0_KINDR / 6.0e0_KINDR
   Real(KINDR),    Parameter :: THREESECOND =  3.0e0_KINDR / 2.0e0_KINDR
   Real(KINDR),    Parameter :: FOURTHIRD   =  4.0e0_KINDR / 3.0e0_KINDR

   Complex(KINDR), Parameter :: I_C     = CMPLX(ZERO,  ONE,  KINDR)
   Complex(KINDR), Parameter :: ZERO_C  = CMPLX(ZERO,  ZERO, KINDR)
   Complex(KINDR), Parameter :: ONE_C   = CMPLX(ONE,   ZERO, KINDR)
   Complex(KINDR), Parameter :: TWO_C   = CMPLX(TWO,   ZERO, KINDR)
   Complex(KINDR), Parameter :: THREE_C = CMPLX(THREE, ZERO, KINDR)
   Complex(KINDR), Parameter :: FOUR_C  = CMPLX(FOUR,  ZERO, KINDR)
   Complex(KINDR), Parameter :: FIVE_C  = CMPLX(FIVE,  ZERO, KINDR)
   Complex(KINDR), Parameter :: SIX_C   = CMPLX(SIX,   ZERO, KINDR)
   Complex(KINDR), Parameter :: HALF_C  = CMPLX(HALF,  ZERO, KINDR)

End Module Constants

!  ============================================================
!  == Constants
!  ============================================================

module constants              ! defines global constants
   implicit none              ! all supposed to be public

   integer,    parameter :: kr   = selected_real_kind(8)
   integer,    parameter :: ki   = selected_int_kind(8)
   real(kr),   parameter :: pi   = 4.0_kr*atan(1.0_kr)
   real(kr),   parameter :: fpi  = 4.0_kr*pi
   real(kr),   parameter :: zero = 0.0_kr
   real(kr),   parameter :: one  = 1.0_kr
   complex(kr),parameter :: ci   = (0.0_kr,1.0_kr)
   real(kr),   parameter :: HartreeToEV = 27.2113845_kr
 
end module constants







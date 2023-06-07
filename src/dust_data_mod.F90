module dust_data_mod

   use chem_constants,  only : kind_chem

   implicit none

   integer, parameter :: ndust = 5

   ! real(kind=kind_chem), dimension(ndust), parameter :: dust_den            = (/   2500.,  2650.,  2650.,  2650.,  2650. /)
   real(kind=kind_chem), dimension(ndust), parameter :: dust_reff          = (/ 0.73D-6, 1.4D-6, 2.4D-6, 4.5D-6, 8.0D-6 /)
   real(kind=kind_chem), dimension(ndust), parameter :: dust_frac_s        = (/     0.1,   0.25,   0.25,   0.25,   0.25 /)
   real(kind=kind_chem), dimension(ndust), parameter :: dust_lower_radius  = (/  0.1D-6, 1.0D-6, 1.8D-6, 3.0D-6, 6.0D-6 /)
   real(kind=kind_chem), dimension(ndust), parameter :: dust_upper_radius  = (/  1.0D-6, 1.8D-6, 3.0D-6, 6.0D-6,10.0D-6 /)

   real(kind=kind_chem), parameter :: dust_alpha = 0.3
   real(kind=kind_chem), parameter :: dust_gamma = 1.3
   real(kind=kind_chem), parameter :: kvhmax=2e-5
   real(kind=kind_chem), parameter :: dust_den=2650. 

end module dust_data_mod

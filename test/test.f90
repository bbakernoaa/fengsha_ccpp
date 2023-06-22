program test

   use catchem_constants, only: rk => kind_chem
   use dust_fengsha_mod, only: DustEmissionFENGSHA

   implicit none


   ! Inputs
   real(rk) :: slc, clay, sandfrac, ssm, rdrag, airdens, ustar, uthrs

   ! Result
   real(rk) :: emissions(5)

   ! Case 1: Full fractions
   slc = 1
   clay = 1
   sandfrac = 1
   ssm = 1
   rdrag = 1
   airdens = 1.2_rk
   ustar = 0.1_rk
   uthrs = 1.0_rk

   call DustEmissionFENGSHA(slc, clay, sandfrac, ssm, rdrag, airdens, ustar, uthrs, emissions)

   print *, emissions

end program test

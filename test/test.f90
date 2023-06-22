program test

   use catchem_constants, only: rk => kind_chem
   use dust_fengsha_mod, only: DustEmissionFENGSHA

   implicit none


   ! Inputs
   real(rk) :: slc, clay, sandfrac, ssm, rdrag, airdens, ustar, uthrs

   ! Result
   real(rk) :: emissions(5)

   ! Full fractions
   slc = 1
   clay = 1
   sandfrac = 1
   ssm = 1
   rdrag = 1
   airdens = 1.2_rk

   ! uthrs < ustar
   ustar = 0.1_rk
   uthrs = 0.05_rk
   call DustEmissionFENGSHA(slc, clay, sandfrac, ssm, rdrag, airdens, ustar, uthrs, emissions)
   print *, emissions
   call assert(all(emissions > 0), "Emissions exist when uthrs < ustar")

   ! No emissions when uthrs >> ustar
   uthrs = 1._rk
   call DustEmissionFENGSHA(slc, clay, sandfrac, ssm, rdrag, airdens, ustar, uthrs, emissions)
   call assert(all(abs(emissions) <= epsilon(emissions)), "No emissions when uthrs >> ustar")

contains

   subroutine assert(condition, message)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: message

      if (.not. condition) then
         print "('Failed assertion: ', a)", message
         error stop
      end if
   end subroutine assert

end program test

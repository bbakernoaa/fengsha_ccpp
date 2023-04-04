module dust_fengsha_mod

  subroutine DustEmissionFENGSHA(slc, clay, silt,  &
                                  ssm, rdrag, airdens, ustar, uthrs, alpha, gamma, &
                                  kvhmax, grav, rhop, area, upper_bin, lower_bin, reff, fecan_soil_moisture, emissions)
    
    ! !USES:
    implicit NONE
    
! !INPUT PARAMETERS:
    REAL(kind_phys), intent(in) :: slc      ! liquid water content of soil layer, volumetric fraction [1]
    REAL(kind_phys), intent(in) :: clay     ! fractional clay content [1]
    REAL(kind_phys), intent(in) :: silt     ! fractional silt content [1]
    REAL(kind_phys), intent(in) :: ssm      ! erosion map [1]
    REAL(kind_phys), intent(in) :: rdrag    ! drag partition [1/m]
    REAL(kind_phys), intent(in) :: airdens  ! air density at lowest level [kg/m^3]
    REAL(kind_phys), intent(in) :: ustar    ! friction velocity [m/sec]
    REAL(kind_phys), intent(in) :: uthrs    ! threshold velocity [m/2]
    REAL(kind_phys), intent(in) :: alpha    ! scaling factor [1]
    REAL(kind_phys), intent(in) :: gamma    ! scaling factor [1]
    REAL(kind_phys), intent(in) :: kvhmax   ! max. vertical to horizontal mass flux ratio [1]
    REAL(kind_phys), intent(in) :: grav     ! gravity [m/sec^2]
    REAL(kind_phys), intent(in) :: rhop     ! soil class density [kg/m^3]
    REAL(kind_phys), intent(in) :: area     ! area of dust emission (can be fractional area of grid cell) 
    INT(kind_phys),  intent(in) :: nbins    ! number of dust bins
    real(kind_phys), dimension(:), intent(in) :: upper_bin ! upper bin size 
    real(kind_phys), dimension(:), intent(in) :: reff ! upper bin size 
    real(kind_phys), dimension(:), intent(in) :: lower_bin ! lower bin size 
    logical        , dimension(:), intent(in) :: fecan_soil_moisture ! use fecan soil mositure or shao soil moisture
    
    ! !OUTPUT PARAMETERS:
    REAL(kind_phys), dimension(:), intent(inout) :: emissions ! binned surface emissions [kg/(m^2 sec)]
    
    ! !DESCRIPTION: Compute dust emissions using NOAA/ARL FENGSHA model
    !
    ! !REVISION HISTORY:
    !
    ! 22Feb2020 B.Baker/NOAA    - Original implementation
    ! 29Mar2021 R.Montuoro/NOAA - Refactored for process library
    ! 09Aug2022 B.Baker/NOAA    - Adapted for CCPP-Physics
    
    ! !Local Variables
    real(kind_phys)                  :: alpha_grav
    real(kind_phys)                  :: h
    real(kind_phys)                  :: kvh
    real(kind_phys)                  :: q
    real(kind_phys)                  :: rustar
    real(kind_phys)                  :: total_emissions
    real(kind_phys)                  :: u_sum, u_thresh
    
    REAL(kind_phys), parameter:: clay_thresh = 0.2
    
!EOP
!-------------------------------------------------------------------------
!  Begin

!  Initialize emissions
!  --------------------
   emissions = 0.

!  Prepare scaling factor
!  ----------------------
   alpha_grav = alpha / grav

   ! Compute vertical-to-horizontal mass flux ratio
   ! ----------------------------------------------
   if (clay > clay_thresh) then
       kvh = kvhmax
   else
       kvh = 10.0**(13.4*clay-6.0)
   end if

   ! Compute total emissions
   ! -----------------------
   total_emissions = alpha_grav * (ssm ** gamma) * airdens * kvh

   !  Compute threshold wind friction velocity using drag partition
   !  -------------------------------------------------------------
   rustar = rdrag * ustar

   !  Now compute size-dependent total emission flux
   !  ----------------------------------------------
   ! Fecan moisture correction
   ! -------------------------
   if (fecan_soil_moisture .eq. .true.) then
     h = moistureCorrectionFecan(slc, sand, clay, rhop)
   else
   ! Shao soil mositure
   !--------------------
     if (slc <= 0.03) then
          h = 
   ! Adjust threshold
   ! ----------------
   u_thresh = uthrs * h
   
   u_sum = rustar + u_thresh
   
   ! Compute Horizontal Saltation Flux according to Eq (9) in Webb et al. (2020)
   ! ---------------------------------------------------------------------------
   q = max(0., rustar - u_thresh) * u_sum * u_sum
   
   ! Distribute emissions to bins and convert to mass flux (kg s-1)
   ! --------------------------------------------------------------
   emissions = distribution(n) * total_emissions * q


 end subroutine DustEmissionFENGSHA
!-----------------------------------------------------------------
   subroutine DustAerosolDistributionKok ( radius, rLow, rUp, distribution )

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   real, dimension(:), intent(in)  :: radius      ! Dry particle bin effective radius [um]
   real, dimension(:), intent(in)  :: rLow, rUp   ! Dry particle bin edge radii [um]

! !OUTPUT PARAMETERS:
   real, dimension(:), intent(out) :: distribution    ! Normalized dust aerosol distribution [1]

! !DESCRIPTION: Computes lognormal aerosol size distribution for dust bins according to
!               J.F.Kok, PNAS, Jan 2011, 108 (3) 1016-1021; doi:10.1073/pnas.1014798108
!
! !REVISION HISTORY:
!
! 22Feb2020 B.Baker/NOAA    - Original implementation
! 01Apr2021 R.Montuoro/NOAA - Refactored for GOCART process library
!

! !Local Variables
   integer :: n, nbins
   real    :: diameter, dlam, dvol

! !CONSTANTS
   real, parameter    :: mmd    = 3.4          ! median mass diameter [um]
   real, parameter    :: stddev = 3.0          ! geometric standard deviation [1]
   real, parameter    :: lambda = 12.0         ! crack propagation length [um]
   real, parameter    :: factor = 1.e0 / (sqrt(2.e0) * log(stddev))  ! auxiliary constant

   character(len=*), parameter :: myname = 'DustAerosolDistributionKok'

!EOP
!-------------------------------------------------------------------------
!  Begin...

   distribution = 0.

!  Assume all arrays are dimensioned consistently
   nbins = size(radius)

   dvol = 0.
   do n = 1, nbins
     diameter = 2 * radius(n)
     dlam = diameter/lambda
     distribution(n) = diameter * (1. + erf(factor * log(diameter/mmd))) * exp(-dlam * dlam * dlam) * log(rUp(n)/rLow(n))
     dvol = dvol + distribution(n)
   end do

!  Normalize distribution
   do n = 1, nbins
     distribution(n) = distribution(n) / dvol
   end do

   end subroutine DustAerosolDistributionKok
!-----------------------------------------------------------------
  real function soilMoistureConvertVol2Grav(vsoil, sandfrac, rhop)

! !USES:
    implicit NONE

! !INPUT PARAMETERS:
    REAL(kind_phys), intent(in) :: vsoil       ! volumetric soil moisture fraction [1]
    REAL(kind_phys), intent(in) :: sandfrac    ! fractional sand content [1]
    REAL(kind_phys), intent(in) :: rhop        ! dry dust density [kg m-3]

! !DESCRIPTION: Convert soil moisture fraction from volumetric to gravimetric.
!
! !REVISION HISTORY:
!
!  02Apr2020, B.Baker/NOAA    - Original implementation
!  01Apr2020, R.Montuoro/NOAA - Adapted for GOCART process library

!  !Local Variables
    real :: vsat

!  !CONSTANTS:
    REAL(kind_phys), parameter :: rhow = 1000.    ! density of water [kg m-3]

!EOP
!-------------------------------------------------------------------------
!  Begin...

!  Saturated volumetric water content (sand-dependent) ! [m3 m-3]
    vsat = 0.489 - 0.00126 * ( 100. * sandfrac )

!  Gravimetric soil content
    soilMoistureConvertVol2Grav = vsoil * rhow / (rhop * (1. - vsat))

  end function soilMoistureConvertVol2Grav
!----------------------------------------------------------------
  real function moistureCorrectionFecan(slc, sand, clay, rhop)

! !USES:
    implicit NONE

! !INPUT PARAMETERS:
    REAL(kind_phys), intent(in) :: slc     ! liquid water content of top soil layer, volumetric fraction [1]
    REAL(kind_phys), intent(in) :: sand    ! fractional sand content [1]
    REAL(kind_phys), intent(in) :: clay    ! fractional clay content [1]
    REAL(kind_phys), intent(in) :: rhop    ! dry dust density [kg m-3]

! !DESCRIPTION: Compute correction factor to account for Fecal soil moisture
!
! !REVISION HISTORY:
!
!  02Apr2020, B.Baker/NOAA    - Original implementation
!  01Apr2020, R.Montuoro/NOAA - Adapted for GOCART process library

!  !Local Variables
    real :: grvsoilm
    real :: drylimit

!EOP
!---------------------------------------------------------------
!  Begin...

!  Convert soil moisture from volumetric to gravimetric
    grvsoilm = soilMoistureConvertVol2Grav(slc, sand, 2650.)

!  Compute fecan dry limit
    drylimit = clay * (14.0 * clay + 17.0)

!  Compute soil moisture correction
    moistureCorrectionFecan = sqrt(1.0 + 1.21 * max(0., grvsoilm - drylimit)**0.68)

  end function moistureCorrectionFecan
!---------------------------------------------------------------
  real function DustFluxV2HRatioMB95(clay, kvhmax)

! !USES:
    implicit NONE

! !INPUT PARAMETERS:
    REAL(kind_phys), intent(in) :: clay      ! fractional clay content [1]
    REAL(kind_phys), intent(in) :: kvhmax    ! maximum flux ratio [1]

!  !CONSTANTS:
    REAL(kind_phys), parameter :: clay_thresh = 0.2    ! clay fraction above which the maximum flux ratio is returned

! !DESCRIPTION: Computes the vertical-to-horizontal dust flux ratio according to
!               B.Marticorena, G.Bergametti, J.Geophys.Res., 100(D8), 164!               doi:10.1029/95JD00690
!
! !REVISION HISTORY:
!
! 22Feb2020 B.Baker/NOAA    - Original implementation
! 01Apr2021 R.Montuoro/NOAA - Adapted for GOCART process library
!
!EOP
!-------------------------------------------------------------------------
!  Begin...

    if (clay > clay_thresh) then
       DustFluxV2HRatioMB95 = kvhmax
    else
       DustFluxV2HRatioMB95 = 10.0**(13.4*clay-6.0)
    end if

  end function DustFluxV2HRatioMB95
  
  
  
end module dust_fengsha_mod

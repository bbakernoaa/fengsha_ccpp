module dust_fengsha_mod

   use catchem_config, only: kind_chem, DUST_OPT_FENGSHA_FECAN
   use catchem_constants
   use dust_data_mod

   implicit none

   private

   public :: DustEmissionFENGSHA


CONTAINS


   subroutine DustEmissionFENGSHA(slc, clay, ssm, rdrag, airdens, ustar, uthrs, area, emissions)

      ! !USES:
      implicit NONE

      ! !INPUT PARAMETERS:
      real(kind=kind_chem), intent(in) :: slc      ! liquid water content of soil layer, volumetric fraction [1]
      real(kind=kind_chem), intent(in) :: clay     ! fractional clay content [1] - range: [0 1]
      ! real(kind=kind_chem), intent(in) :: silt     ! fractional silt content [1] - range: [0 1]
      real(kind=kind_chem), intent(in) :: ssm      ! erosion map [1] - range: [0 1]
      real(kind=kind_chem), intent(in) :: rdrag    ! drag partition [1/m] - range: [0 1]
      real(kind=kind_chem), intent(in) :: airdens  ! air density at lowest level [kg/m^3]
      real(kind=kind_chem), intent(in) :: ustar    ! friction velocity [m/sec]
      real(kind=kind_chem), intent(in) :: uthrs    ! threshold velocity [m/2]
      ! real(kind=kind_chem), intent(in) :: kvhmax   ! max. vertical to horizontal mass flux ratio [1]
      ! real(kind=kind_chem), intent(in) :: rhop     ! soil class density [kg/m^3]
      real(kind=kind_chem), intent(in) :: area     ! area of dust emission (can be fractional area of grid cell)

      ! this will need to be read in by a configuration file
      ! right now I just hard coded it in dust_data_mod.F90
      !=====================================================
      ! real(kind=kind_chem), intent(in) :: alpha    ! scaling factor [1]
      ! real(kind=kind_chem), intent(in) :: gamma    ! scaling factor [1]
      ! integer(kind=kind_chem),  intent(in) :: nbins    ! number of dust bins
      ! real(kind=kind_chem), dimension(:), intent(in) :: upper_bin ! upper bin size
      ! real(kind=kind_chem), dimension(:), intent(in) :: reff ! upper bin size
      ! real(kind=kind_chem), dimension(:), intent(in) :: lower_bin ! lower bin size
      !=====================================================

      ! !OUTPUT PARAMETERS:
      REAL(kind=kind_chem), dimension(:), intent(inout) :: emissions ! binned surface emissions [kg/(m^2 sec)]

      ! !DESCRIPTION: Compute dust emissions using NOAA/ARL FENGSHA model
      !
      ! !REVISION HISTORY:
      !
      ! 22Feb2020 B.Baker/NOAA    - Original implementation
      ! 29Mar2021 R.Montuoro/NOAA - Refactored for process library
      ! 09Aug2022 B.Baker/NOAA    - Adapted for CCPP-Physics

      ! !Local Variables
      real(kind=kind_chem) :: alpha_grav
      real(kind=kind_chem) :: h
      real(kind=kind_chem) :: kvh
      real(kind=kind_chem) :: q
      real(kind=kind_chem) :: rustar
      real(kind=kind_chem) :: total_emissions
      real(kind=kind_chem) :: u_sum, u_thresh
      real(kind=kind_chem) :: distribution(size(emissions))
      real(kind=kind_chem) :: vsat, sandfrac, grvsoilm, vsoil, drylimit

      real(kind=kind_chem), parameter:: clay_thresh = 0.2
      real(kind=kind_chem), parameter :: rhow = 1000.

      !EOP
      !-------------------------------------------------------------------------
      !  Begin

      !  Initialize emissions
      !  --------------------
      emissions = 0.

      !  Prepare scaling factor
      !  ----------------------
      alpha_grav = dust_alpha / con_g

      ! Compute vertical-to-horizontal mass flux ratio
      ! ----------------------------------------------
      if (clay > clay_thresh) then
         kvh = kvhmax
      else
         kvh = 10.0**(13.4*clay-6.0)
      end if

      ! Compute total emissions
      ! -----------------------
      total_emissions = alpha_grav * (ssm ** dust_gamma) * airdens * kvh

      !  Compute threshold wind friction velocity using drag partition
      !  -------------------------------------------------------------
      rustar = rdrag * ustar

      !  Now compute size-dependent total emission flux
      !  ----------------------------------------------
      ! Fecan moisture correction
      ! -------------------------
      if (DUST_OPT_FENGSHA_FECAN .eqv. .true.) then
         !  Saturated volumetric water content (sand-dependent) ! [m3 m-3]
         vsat = 0.489 - 0.00126 * ( 100. * sandfrac )

         !  Gravimetric soil content
         grvsoilm = vsoil * rhow / (dust_den * (1. - vsat))

         !   Compute fecan dry limit
         drylimit = clay * (14.0 * clay + 17.0)

         !  Compute soil moisture correction
         h = sqrt(1.0 + 1.21 * max(0._kind_chem, grvsoilm - drylimit)**0.68)
      else
         ! Shao soil mositure
         !--------------------
         if (slc <= 0.03) then
            h = exp(22.7 * slc)
         else
            h = exp(93.5 * slc - 2.029)
         end if
      end if

      ! Adjust threshold
      ! ----------------
      u_thresh = uthrs * h

      u_sum = rustar + u_thresh

      ! Compute Horizontal Saltation Flux according to Eq (9) in Webb et al. (2020)
      ! ---------------------------------------------------------------------------
      q = max(0._kind_chem, rustar - u_thresh) * u_sum * u_sum

      ! get distribution 
      call DustAerosolDistributionKok(dust_reff, dust_lower_radius, dust_upper_radius, distribution)

      ! Distribute emissions to bins and convert to mass flux (kg s-1)
      ! --------------------------------------------------------------
      emissions = distribution * total_emissions * q


   end subroutine DustEmissionFENGSHA
   !-----------------------------------------------------------------
   subroutine DustAerosolDistributionKok ( radius, rLow, rUp, distribution )

      ! !USES:
      implicit NONE

      ! !INPUT PARAMETERS:
      real(kind=kind_chem), dimension(:), intent(in)  :: radius      ! Dry particle bin effective radius [um]
      real(kind=kind_chem), dimension(:), intent(in)  :: rLow, rUp   ! Dry particle bin edge radii [um]

      ! !OUTPUT PARAMETERS:
      real(kind=kind_chem), dimension(:), intent(out) :: distribution    ! Normalized dust aerosol distribution [1]

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
      real(kind=kind_chem)    :: diameter, dlam, dvol

      !   !CONSTANTS
      real(kind=kind_chem), parameter    :: mmd    = 3.4          ! median mass diameter [um]
      real(kind=kind_chem), parameter    :: stddev = 3.0          ! geometric standard deviation [1]
      real(kind=kind_chem), parameter    :: lambda = 12.0         ! crack propagation length [um]
      real(kind=kind_chem), parameter    :: factor = 1.e0 / (sqrt(2.e0) * log(stddev))  ! auxiliary constant

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

end module dust_fengsha_mod

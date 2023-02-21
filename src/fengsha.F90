module dust_fengsha_mod

  subroutine fengsha_drag(z0,R)
    implicit none

    real(kind_phys), intent(in) :: z0
    real(kind_phys), intent(out) :: R
    real(kind_phys), parameter :: z0s = 1.0e-04 !Surface roughness for ideal bare surface [m]
    ! ------------------------------------------------------------------------
    ! Function: Calculates the MacKinnon et al. 2004 Drag Partition Correction
    !
    !   R = 1.0 - log(z0 / z0s) / log( 0.7 * (12255./z0s) ** 0.8)
    !
    !--------------------------------------------------------------------------
    ! Drag partition correction. See MacKinnon et al. (2004),
    !     doi:10.1016/j.geomorph.2004.03.009
    R = 1.0 - log(z0 / z0s) / log( 0.7 * (12255./z0s) ** 0.8)

    ! Drag partition correction. See Marticorena et al. (1997),
    !     doi:10.1029/96JD02964
    !R = 1.0 - log(z0 / z0s) / log( 0.7 * (10./z0s) ** 0.8)

    return
  end subroutine fengsha_drag

  subroutine DustEmissionFENGSHA(slc, clay, sand, silt,  &
                                  ssm, rdrag, airdens, ustar, uthrs, alpha, gamma, &
                                  kvhmax, grav, rhop, emissions)
    
    ! !USES:
    implicit NONE
    
! !INPUT PARAMETERS:
    REAL(kind_phys), intent(in) :: slc      ! liquid water content of soil layer, volumetric fraction [1]
    REAL(kind_phys), intent(in) :: clay     ! fractional clay content [1]
    REAL(kind_phys), intent(in) :: sand     ! fractional sand content [1]
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
    
    ! !OUTPUT PARAMETERS:
    REAL(kind_phys), intent(inout) :: emissions ! binned surface emissions [kg/(m^2 sec)]
    
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
   kvh = DustFluxV2HRatioMB95(clay, kvhmax)

   ! Compute total emissions
   ! -----------------------
   emissions = alpha_grav * (ssm ** gamma) * airdens * kvh

   !  Compute threshold wind friction velocity using drag partition
   !  -------------------------------------------------------------
   rustar = rdrag * ustar

   !  Now compute size-dependent total emission flux
   !  ----------------------------------------------
   ! Fecan moisture correction
   ! -------------------------
   h = moistureCorrectionFecan(slc, sand, clay, rhop)
   
   ! Adjust threshold
   ! ----------------
   u_thresh = uthrs * h
   
   u_sum = rustar + u_thresh
   
   ! Compute Horizontal Saltation Flux according to Eq (9) in Webb et al. (2020)
   ! ---------------------------------------------------------------------------
   q = max(0., rustar - u_thresh) * u_sum * u_sum
   
   ! Distribute emissions to bins and convert to mass flux (kg s-1)
   ! --------------------------------------------------------------
   emissions = emissions * q


 end subroutine DustEmissionFENGSHA
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

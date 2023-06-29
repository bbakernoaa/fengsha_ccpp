!>\file catchem_chem_dust_wrapper.F90
!! This file is CATChem dust wrapper with CCPP coupling to FV3
!! Barry.Baker@noaa.gov 05/2023

module catchem_chem_dust_wrapper

   use catchem_constants, only : kind_phys => kind_chem, g => con_g, pi => con_pi
   use catchem_config
!    use dust_gocart_mod, only : gocart_dust_driver (comment out for now)
!    use dust_afwa_mod,   only : gocart_dust_afwa_driver (comment out for now)
   use dust_fengsha_mod, only : DustEmissionFENGSHA
   use dust_data_mod

   implicit none

   private

   public :: catchem_chem_dust_wrapper_init, catchem_chem_dust_wrapper_run, catchem_chem_dust_wrapper_finalize

contains

!> \brief Brief description of the subroutine
!!
   subroutine catchem_chem_dust_wrapper_init()
   end subroutine catchem_chem_dust_wrapper_init

!> \brief Brief description of the subroutine
!!
!! \section arg_table_catchem_chem_dust_wrapper_finalize Argument Table
!!
   subroutine catchem_chem_dust_wrapper_finalize()
   end subroutine catchem_chem_dust_wrapper_finalize

!> \defgroup catchem_chem_dust_group catchem Chem seas wrapper Module
!! This is the catchem chemistry
!>\defgroup catchem_chem_dust_wrapper catchem Chem seas wrapper Module
!> \ingroup catchem_chem_dust_group
!! This is the catchem Chem seas wrapper Module
!! \section arg_table_catchem_chem_dust_wrapper_run Argument Table
!! \htmlinclude catchem_chem_dust_wrapper_run.html
!!
!>\section catchem_chem_dust_wrapper catchem Chemistry Scheme General Algorithm
!> @{
   subroutine catchem_chem_dust_wrapper_run( &
      im, kte, kme, ktau, dt, garea, land, &
      ustar, &
      pr3d, ph3d,phl3d, prl3d, tk3d, us3d, vs3d, spechum, &
      nsoil, smc, vegtype, soiltyp, snow_cplchm, &
      dust_in, ntrac, ntdust1,ntdust2,ntdust3,ntdust4,ntdust5,ndust, &
      gq0,qgrs,duem, chem_opt_in,dust_opt_in, &
      dust_alpha_in,dust_gamma_in,pert_scale_dust, &
      emis_amp_dust, do_sppt_emis, sppt_wts, errmsg,errflg)

      implicit none

      integer, intent(in) :: im         ! horizontal_loop_extent
      integer, intent(in) :: kte        ! vertical_layer_dimension
      integer, intent(in) :: kme        ! number of vertical levels plus one
      integer, intent(in) :: ktau       ! current forecast iteration
      integer, intent(in) :: nsoil      ! soil vertical layer dimension
      integer, intent(in) :: ntrac      ! number of tracers
      integer, intent(in) :: ntdust1    ! index_for_dust_bin1      !!!! Hard coded for now....need to make this dynamic with the number of bins/tracers and give modal option
      integer, intent(in) :: ntdust2    ! index_for_dust_bin2      !!!! Hard coded for now....need to make this dynamic with the number of bins/tracers and give modal option
      integer, intent(in) :: ntdust3    ! index_for_dust_bin3      !!!! Hard coded for now....need to make this dynamic with the number of bins/tracers and give modal option
      integer, intent(in) :: ntdust4    ! index_for_dust_bin4      !!!! Hard coded for now....need to make this dynamic with the number of bins/tracers and give modal option
      integer, intent(in) :: ntdust5    ! index_for_dust_bin5      !!!! Hard coded for now....need to make this dynamic with the number of bins/tracers and give modal option
      integer, intent(in) :: ndust      ! number_of_dust_bins_for_diagnostics !!!! Hard coded for now....need to make this dynamic with the number of bins/tracers and give modal option

      real(kind_phys), intent(in) :: dt                   ! timestep_for_physics
      real(kind_phys), intent(in) :: emis_amp_dust        ! dust_emissions_perturbation_amplitude
      real(kind_phys), intent(in) :: pert_scale_dust      ! dust_emissions_scaling_factor

      ! sppt params 
      logical,         intent(in) :: do_sppt_emis                   ! flag_for_stochastic_emissions_perturbations
      real(kind=kind_phys), optional, intent(in) :: sppt_wts(:,:)   ! 

      ! surface description params
      integer, dimension(im), intent(in) :: land          ! landmask: sea/land/ice=0/1/2
      integer, dimension(im), intent(in) :: soiltyp       ! soil type at each grid cell
      integer, dimension(im), intent(in) :: vegtype       ! veg type at each grid cell

      integer, parameter :: ids=1, jds=1, jde=1, kds=1
      integer, parameter :: ims=1, jms=1, jme=1, kms=1
      integer, parameter :: its=1, jts=1, jte=1, kts=1

      real(kind=kind_phys), dimension(im), intent(in) :: garea  ! grid cell area

      real(kind=kind_phys), dimension(im, nsoil), intent(in) :: smc       ! volumetric fraction of soil moisture for lsm
      real(kind=kind_phys), dimension(im, ndust), intent(in) :: dust_in   ! fengsha dust input !!! need to think about this or the future

      ! physics params
      real(kind=kind_phys), dimension(im), intent(in) :: ustar ! 
      real(kind=kind_phys), dimension(ims:im, nsoil, jms:jme) :: smois
      real(kind_phys), dimension(im), intent(in) :: snow_cplchm
      real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_chem )  :: chem

      !>- dust & chemistry variables
      real(kind_phys), dimension(ims:im, jms:jme) :: ssm              ! SSM - in new scheme this is 1 always  
      real(kind_phys), dimension(ims:im, jms:jme) :: rdrag            ! drag parition (-)
      real(kind_phys), dimension(ims:im, jms:jme) :: uthr             ! dry threshold friction velocity    
      real(kind_phys), dimension(ims:im, jms:jme) :: snowh            ! snow height 
      real(kind_phys), dimension(ims:im, jms:jme) :: xland            ! landmask: sea/land/ice=0/1/2
      real(kind=kind_phys), dimension(ims:im, jms:jme) :: sandf       ! sand fraction 
      real(kind=kind_phys), dimension(ims:im, jms:jme) :: clayf       ! clay fraction 
      real(kind_phys), dimension(ims:im, jms:jme, 1:ndust) :: dust_emis
      integer,         dimension(ims:im, jms:jme) :: isltyp, ivgtyp
      real(kind_phys), dimension(ims:im, kms:kme, jms:jme, 1:num_moist)  :: moist

      real(kind_phys), parameter :: ugkg = 1.e-09_kind_phys !lzhang
      real(kind_phys), dimension(1:num_chem) :: ppm2ugkg
      real(kind_phys) :: dtstep

      ! Local Variables 
      integer :: i, j, k
      integer :: ide, ime, ite, kde
      real(kind_phys), dimension(ims:im,jms:jme) :: random_factor


      !==========================================================================================
      real(kind_phys), dimension(im,kme), intent(in) :: ph3d        ! geopotential at model layer interfaces (m2 s-2)
      real(kind_phys), dimension(im,kme), intent(in) :: pr3d        ! air pressure at model layer interfaces (Pa)
      real(kind_phys), dimension(im,kte), intent(in) :: phl3d       ! geopotential at model layer centers  (m2 s-2)
      real(kind_phys), dimension(im,kme), intent(in) :: prl3d       ! mean layer pressure (Pa)
      real(kind_phys), dimension(im,kme), intent(in) :: tk3d        ! updated temperature (K) 
      real(kind_phys), dimension(im,kme), intent(in) :: us3d        ! zonal wind updated by physics (m s-1)
      real(kind_phys), dimension(im,kme), intent(in) :: vs3d        ! meridional wind updated by physics (m s-1)
      real(kind_phys), dimension(im,kme), intent(in) :: spechum     ! water vapor specific humidity updated by physics (kg kg-1)
      real(kind_phys), dimension(im,kte,ntrac), intent(inout) :: gq0     ! tracer concentration updated by physics
      real(kind_phys), dimension(im,kte,ntrac), intent(inout) :: qgrs    ! model layer mean tracer concentration
      real(kind_phys), dimension(im,ndust    ), intent(inout) :: duem    ! instantaneous dust emission flux
      integer,        intent(in) :: chem_opt_in      ! CATchem option
      integer,        intent(in) :: dust_opt_in
      real(kind_phys),intent(in) :: dust_alpha_in
      integer,        intent(in) :: dust_gamma_in
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      real(kind_phys), dimension(1:im, 1:kme,jms:jme) :: t_phy, &
         p_phy, z_at_w, dz8w, p8w, t8w, rho_phy

      real(kind_phys), dimension(ims:im, jms:jme) :: ust, dxy
      !================================================================================


      errmsg = ''
      errflg = 0

      chem_opt          = chem_opt_in
      dust_opt          = dust_opt_in
      ! dust_calcdrag     = dust_calcdrag_in
      chem = 0.

      ! -- initialize dust emissions
      dust_emis = 0._kind_phys

      ! -- set domain
      ide=im
      ime=im
      ite=im
      kde=kte

      if(do_sppt_emis) then
         random_factor(:,jms) = pert_scale_dust*max( &
            min(1+(sppt_wts(:,kme/2)-1)*emis_amp_dust, 2.0_kind_phys), &
            0.0_kind_phys)
      else
         random_factor = 1.0
      endif

      ! -- volume to mass fraction conversion table (ppm -> ug/kg)
      ppm2ugkg         = 1._kind_phys
      ppm2ugkg(p_sulf) = 1.e+03_kind_phys * mw_so4_aer / mwdry

      ! -- compute accumulated large-scale and convective rainfall since last call
      if (ktau > 1) then
         dtstep = call_chemistry * dt
      else
         dtstep = dt
      end if

!>- get ready for chemistry run
      call catchem_chem_prep_dust(                                             &
         ktau,dtstep,                                                     &
         ustar,land,garea,                     &
         pr3d,ph3d,phl3d,tk3d,prl3d,spechum,                    &
         nsoil,smc,vegtype,soiltyp,                    &
         snow_cplchm,dust_in,                                      &
         ust,xland,dxy,                  &
         t_phy,p_phy,rho_phy,dz8w,p8w,t8w,z_at_w,         &
         ntdust1,ntdust2,ntdust3,ntdust4,ntdust5,                         & !!! HARD CODED Need a better way to do this 
         ntrac,gq0,num_chem, num_moist,ppm2ugkg,moist,chem,               &
         smois,ivgtyp,isltyp,                &
         snowh,clayf,rdrag,sandf,ssm,uthr,                           &
         ids,ide, jds,jde, kds,kde,                                       &
         ims,ime, jms,jme, kms,kme,                                       &
         its,ite, jts,jte, kts,kte)


      !-- compute dust
      !store_arrays = .false.
      select case (dust_opt)
       case (DUST_OPT_FENGSHA)
        !  dust_alpha    = dust_alpha_in  ! dust_alpha  !!! should be set in the chem side not in the physics side
        !  dust_gamma    = dust_gamma_in  ! dust_alpha  !!! should be set in the chem side not in the physics side
         do j=jts,jte
            do i=its,ite

               call DustEmissionFENGSHA(smois(i,1,j), clayf(i,j), sandf(i,j), ssm(i,j), &
                  rdrag(i,j), rho_phy(i,j,1), &
                  ust(i,j), uthr(i,j), &
                  dust_emis(i,j,:))
               dust_emis(i,j,:) = dust_emis(i,j,:) * dxy(i,j) * dt * random_factor(i,j) ! for sppt if activated

               ! modify concentrations 
               ! first convert concentrations to kg/kg 
               chem(i,kts,j,p_dust_1)=chem(i,kts,j,p_dust_1) * 1.0e-9
               chem(i,kts,j,p_dust_2)=chem(i,kts,j,p_dust_2) * 1.0e-9
               chem(i,kts,j,p_dust_3)=chem(i,kts,j,p_dust_3) * 1.0e-9
               chem(i,kts,j,p_dust_4)=chem(i,kts,j,p_dust_4) * 1.0e-9
               chem(i,kts,j,p_dust_5)=chem(i,kts,j,p_dust_5) * 1.0e-9

               ! now add dust emissions back to cocentration fields
               chem(i,kts,j,p_dust_1)=chem(i,kts,j,p_dust_1) + dust_emis(i,j,1) / rho_phy(i,j,1)
               chem(i,kts,j,p_dust_2)=chem(i,kts,j,p_dust_2) + dust_emis(i,j,2) / rho_phy(i,j,1)
               chem(i,kts,j,p_dust_3)=chem(i,kts,j,p_dust_3) + dust_emis(i,j,1) / rho_phy(i,j,1)
               chem(i,kts,j,p_dust_4)=chem(i,kts,j,p_dust_4) + dust_emis(i,j,1) / rho_phy(i,j,1)
               chem(i,kts,j,p_dust_5)=chem(i,kts,j,p_dust_5) + dust_emis(i,j,1) / rho_phy(i,j,1)

               ! convert concentrations back to ug/kg 
               chem(i,kts,j,p_dust_1)=chem(i,kts,j,p_dust_1) * 1.e9
               chem(i,kts,j,p_dust_2)=chem(i,kts,j,p_dust_2) * 1.e9
               chem(i,kts,j,p_dust_3)=chem(i,kts,j,p_dust_3) * 1.e9
               chem(i,kts,j,p_dust_4)=chem(i,kts,j,p_dust_4) * 1.e9
               chem(i,kts,j,p_dust_5)=chem(i,kts,j,p_dust_5) * 1.e9
            end do
         end do   

       case default
         errmsg = 'Logic error in catchem_chem_dust_wrapper_run: invalid dust_opt'
         errflg = 1
         return
         !store_arrays = .true.
      end select

      ! -- put chem stuff back into tracer array
      do k=kts,kte
         do i=its,ite
            gq0(i,k,ntdust1  )=ppm2ugkg(p_dust_1) * max(epsilc,chem(i,k,1,p_dust_1))
            gq0(i,k,ntdust2  )=ppm2ugkg(p_dust_2) * max(epsilc,chem(i,k,1,p_dust_2))
            gq0(i,k,ntdust3  )=ppm2ugkg(p_dust_3) * max(epsilc,chem(i,k,1,p_dust_3))
            gq0(i,k,ntdust4  )=ppm2ugkg(p_dust_4) * max(epsilc,chem(i,k,1,p_dust_4))
            gq0(i,k,ntdust5  )=ppm2ugkg(p_dust_5) * max(epsilc,chem(i,k,1,p_dust_5))
         enddo
      enddo

      do k=kts,kte
         do i=its,ite
            qgrs(i,k,ntdust1 )=gq0(i,k,ntdust1  )
            qgrs(i,k,ntdust2 )=gq0(i,k,ntdust2  )
            qgrs(i,k,ntdust3 )=gq0(i,k,ntdust3  )
            qgrs(i,k,ntdust4 )=gq0(i,k,ntdust4  )
            qgrs(i,k,ntdust5 )=gq0(i,k,ntdust5  )
         enddo
      enddo

      duem(:,:) = ugkg*dust_emis(:,:,1)  ! FIXME
!
   end subroutine catchem_chem_dust_wrapper_run
!> @}

   subroutine catchem_chem_prep_dust(                                      &
      ktau,dtstep,                     &
      ustar,land,garea,                     &
      pr3d,ph3d,phl3d,tk3d,prl3d,spechum,                &
      nsoil,smc,vegtype,soiltyp,                 &
      snow_cplchm,dust_in,                               &
      ust,xland,dxy,                          &
      t_phy,p_phy,rho_phy,dz8w,p8w,                  &
      t8w,                                              &
      z_at_w,                                              &
      ntdust1,ntdust2,ntdust3,ntdust4,ntdust5,               &
      ntrac,gq0,                                                     &
      num_chem, num_moist,                                &
      ppm2ugkg,                                             &
      moist,chem,                                   &
      smois,ivgtyp,isltyp,              &
      snowh,clayf,rdrag,sandf,ssm,uthr,            &
      ids,ide, jds,jde, kds,kde,                                     &
      ims,ime, jms,jme, kms,kme,                                     &
      its,ite, jts,jte, kts,kte)
      ! NOTE: currently many of the args above are unused

      !Chem input configuration
      integer, intent(in) :: ktau
      real(kind=kind_phys), intent(in) :: dtstep

      !catchem Chem variables
      integer,intent(in) ::  num_chem, num_moist
      integer,intent(in) ::  ids,ide, jds,jde, kds,kde,                      &
         ims,ime, jms,jme, kms,kme,                      &
         its,ite, jts,jte, kts,kte

      !FV3 input variables
      integer, intent(in) :: nsoil
      integer, dimension(ims:ime), intent(in) :: land
      integer, dimension(ims:ime), intent(in) :: vegtype
      integer, dimension(ims:ime), intent(in) :: soiltyp
      integer, intent(in) :: ntrac
      integer, intent(in) :: ntdust1
      integer, intent(in) :: ntdust2
      integer, intent(in) :: ntdust3
      integer, intent(in) :: ntdust4
      integer, intent(in) :: ntdust5

      real(kind=kind_phys), dimension(ims:ime), intent(in) :: ustar 
      real(kind=kind_phys), dimension(ims:ime), intent(in) :: garea 
      ! real(kind=kind_phys), dimension(ims:ime), intent(in) :: ts2d 
      ! real(kind=kind_phys), dimension(ims:ime), intent(in) :: sigmaf
      ! real(kind=kind_phys), dimension(ims:ime), intent(in) :: dswsfc
      real(kind=kind_phys), dimension(ims:ime), intent(in) :: snow_cplchm 
      ! real(kind=kind_phys), dimension(ims:ime), intent(in) :: hf2d
      ! real(kind=kind_phys), dimension(ims:ime), intent(in) :: pb2d
      real(kind=kind_phys), dimension(ims:ime, nsoil),   intent(in) :: smc
      real(kind=kind_phys), dimension(ims:ime, ndust),   intent(in) :: dust_in
      ! real(kind=kind_phys), dimension(ims:ime,    10),   intent(in) :: emi_in
      real(kind=kind_phys), dimension(ims:ime, kms:kme), intent(in) :: pr3d
      real(kind=kind_phys), dimension(ims:ime, kms:kme), intent(in) :: ph3d
      real(kind=kind_phys), dimension(ims:ime, kts:kte), intent(in) :: phl3d
      real(kind=kind_phys), dimension(ims:ime, kms:kme), intent(in) :: tk3d
      real(kind=kind_phys), dimension(ims:ime, kms:kme), intent(in) :: prl3d
      ! real(kind=kind_phys), dimension(ims:ime, kms:kme), intent(in) :: us3d
      ! real(kind=kind_phys), dimension(ims:ime, kms:kme), intent(in) :: vs3d
      real(kind=kind_phys), dimension(ims:ime, kms:kme), intent(in) ::spechum
      real(kind=kind_phys), dimension(ims:ime, kts:kte,ntrac), intent(in) :: gq0

      real(kind_phys), dimension(num_chem), intent(in) :: ppm2ugkg

      integer,dimension(ims:ime, jms:jme), intent(out) :: isltyp, ivgtyp
      ! real(kind_phys), dimension(ims:ime, jms:jme, 3), intent(inout) :: erod
      real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: t_phy
      real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: p_phy
      real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: rho_phy
      real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: dz8w
      real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: p8w
      real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: t8w
      real(kind_phys), dimension(ims:ime, jms:jme), intent(out) :: ust
      real(kind_phys), dimension(ims:ime, jms:jme), intent(out) :: xland
      real(kind_phys), dimension(ims:ime, jms:jme), intent(out) :: dxy
      real(kind_phys), dimension(ims:ime, jms:jme), intent(out) :: snowh
      real(kind_phys), dimension(ims:ime, jms:jme), intent(out) :: clayf
      real(kind_phys), dimension(ims:ime, jms:jme), intent(out) :: rdrag
      real(kind_phys), dimension(ims:ime, jms:jme), intent(out) :: sandf
      real(kind_phys), dimension(ims:ime, jms:jme), intent(out) :: ssm
      real(kind_phys), dimension(ims:ime, jms:jme), intent(out) :: uthr
      real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, num_moist), intent(out) :: moist
      real(kind_phys), dimension(ims:ime, kms:kme, jms:jme, num_chem),  intent(out) :: chem

      real(kind_phys), dimension(ims:ime, kms:kme, jms:jme), intent(out) :: z_at_w
      real(kind_phys), dimension(ims:ime, nsoil, jms:jme), intent(out) :: smois

      ! -- local variables
!   real(kind=kind_phys), dimension(ims:ime, kms:kme, jms:jme) :: p_phy
      integer i,ip,j,jp,k,kp,kk,kkp


      ! -- initialize output arrays
      isltyp         = 0._kind_phys
      ivgtyp         = 0._kind_phys
      ! rri            = 0._kind_phys
      t_phy          = 0._kind_phys
      ! u_phy          = 0._kind_phys
      ! v_phy          = 0._kind_phys
      p_phy          = 0._kind_phys
      rho_phy        = 0._kind_phys
      dz8w           = 0._kind_phys
      p8w            = 0._kind_phys
      t8w            = 0._kind_phys
      ! u10            = 0._kind_phys
      ! v10            = 0._kind_phys
      ust            = 0._kind_phys
      ! tsk            = 0._kind_phys
      xland          = 0._kind_phys
      ! xlat           = 0._kind_phys
      ! xlong          = 0._kind_phys
      dxy            = 0._kind_phys
      ! vegfrac        = 0._kind_phys
      ! rmol           = 0._kind_phys
      ! gsw            = 0._kind_phys
      ! znt            = 0._kind_phys
      ! hfx            = 0._kind_phys
      ! pbl            = 0._kind_phys
      snowh          = 0._kind_phys
      clayf          = 0._kind_phys
      rdrag          = 0._kind_phys
      sandf          = 0._kind_phys
      ssm            = 0._kind_phys
      uthr           = 0._kind_phys
      moist          = 0._kind_phys
      chem           = 0._kind_phys
      z_at_w         = 0._kind_phys


      do i=its,ite
         ust  (i,1)=ustar(i)
         dxy  (i,1)=garea(i)
         xland(i,1)=real(land(i))
         snowh(i,1)=snow_cplchm(i)*0.001
         clayf(i,1)=dust_in(i,1) !!! again need to think how this is coming in 
         rdrag(i,1)=dust_in(i,2) !!! again need to think how this is coming in
         sandf(i,1)=dust_in(i,3) !!! again need to think how this is coming in
         ssm  (i,1)=dust_in(i,4) !!! again need to think how this is coming in
         uthr (i,1)=dust_in(i,5) !!! again need to think how this is coming in
         ivgtyp (i,1)=vegtype(i) 
         isltyp (i,1)=soiltyp(i)
      enddo

      ! rmol=0.

      do k=1,nsoil
         do j=jts,jte
            do i=its,ite
               smois(i,k,j)=smc(i,k)
            enddo
         enddo
      enddo

      do j=jts,jte
         jp = j - jts + 1
         do i=its,ite
            ip = i - its + 1
            z_at_w(i,kts,j)=max(0._kind_phys, ph3d(ip,1)/g)
         enddo
      enddo

      do j=jts,jte
         jp = j - jts + 1
         do k=kts,kte
            kp = k - kts + 1
            do i=its,ite
               ip = i - its + 1
               dz8w(i,k,j)=abs(ph3d(ip,kp+1)-ph3d(ip,kp))/g
               z_at_w(i,k+1,j)=z_at_w(i,k,j)+dz8w(i,k,j)
            enddo
         enddo
      enddo

      do j=jts,jte
         jp = j - jts + 1
         do k=kts,kte+1
            kp = k - kts + 1
            do i=its,ite
               ip = i - its + 1
               p8w(i,k,j)=pr3d(ip,kp)
            enddo
         enddo
      enddo

      do j=jts,jte
         jp = j - jts + 1
         do k=kts,kte+1
            kk=min(k,kte)
            kkp = kk - kts + 1
            do i=its,ite
               ip = i - its + 1
               dz8w(i,k,j)=z_at_w(i,kk+1,j)-z_at_w(i,kk,j)
               t_phy(i,k,j)=tk3d(ip,kkp)
               p_phy(i,k,j)=prl3d(ip,kkp)
               ! u_phy(i,k,j)=us3d(ip,kkp)
               ! v_phy(i,k,j)=vs3d(ip,kkp)
               rho_phy(i,k,j)=p_phy(i,k,j)/(287.04*t_phy(i,k,j)*(1.+.608*spechum(ip,kkp)))
               ! rri(i,k,j)=1./rho_phy(i,k,j)
               moist(i,k,j,:)=0.
               moist(i,k,j,1)=gq0(ip,kkp,p_atm_shum)
               if (t_phy(i,k,j) > 265.) then
                  moist(i,k,j,2)=gq0(ip,kkp,p_atm_cldq)
                  moist(i,k,j,3)=0.
                  if (moist(i,k,j,2) < 1.e-8) moist(i,k,j,2)=0.
               else
                  moist(i,k,j,2)=0.
                  moist(i,k,j,3)=gq0(ip,kkp,p_atm_cldq)
                  if(moist(i,k,j,3) < 1.e-8)moist(i,k,j,3)=0.
               endif
               !--
            enddo
         enddo
      enddo

      do j=jts,jte
         do k=2,kte
            do i=its,ite
               t8w(i,k,j)=.5*(t_phy(i,k,j)+t_phy(i,k-1,j))
            enddo
         enddo
      enddo

      ! -- only used in phtolysis....
      do j=jts,jte
         do i=its,ite
            t8w(i,1,j)=t_phy(i,1,j)
            t8w(i,kte+1,j)=t_phy(i,kte,j)
         enddo
      enddo



      do k=kms,kte
         do i=ims,ime
            chem(i,k,jts,p_dust_1)=max(epsilc,gq0(i,k,ntdust1)/ppm2ugkg(p_dust_1))
            chem(i,k,jts,p_dust_2)=max(epsilc,gq0(i,k,ntdust2)/ppm2ugkg(p_dust_2))
            chem(i,k,jts,p_dust_3)=max(epsilc,gq0(i,k,ntdust3)/ppm2ugkg(p_dust_3))
            chem(i,k,jts,p_dust_4)=max(epsilc,gq0(i,k,ntdust4)/ppm2ugkg(p_dust_4))
            chem(i,k,jts,p_dust_5)=max(epsilc,gq0(i,k,ntdust5)/ppm2ugkg(p_dust_5))
         enddo
      enddo


      ! -- real-time application, keeping eruption constant

   end subroutine catchem_chem_prep_dust


!> @}
end module catchem_chem_dust_wrapper

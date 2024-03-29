module tphysbc_ml
    !----------------------------------------------- !
    ! Purpose:                                       !
    !                                                !
    ! Dropin replacement of the parametrized frow of !
    ! tphysbc,                                       !
    !                                                !
    !----------------------------------------------- !


implicit none
private
save

public :: tphysbc_ml_main

contains

subroutine tphysbc_ml_main (ztodt,   pblht,   tpert,   qpert,    &
                    fsns,    fsnt,    flns,    flnt,    state,   &
                    tend,    pbuf,    fsds,    landm,            &
                    cam_out, cam_in )


!----------------------------------------------------------------------- 
! 
! Purpose: 
! Evaluate and apply physical processes that are calculated BEFORE 
! coupling to land, sea, and ice models.  
!
! Processes currently included are: 
! dry adjustment, moist convection, stratiform, wet deposition, radiation
!
! Pass surface fields for separate surface flux calculations
! Dump appropriate fields to history file.
! 
! Method: 
!
! Each parameterization should be implemented with this sequence of calls:
!  1)  Call physics interface
!  2)  Check energy
!  3)  Call physics_update
! See Interface to Column Physics and Chemistry Packages 
!   http://www.ccsm.ucar.edu/models/atm-cam/docs/phys-interface/index.html
! 
! Author: CCM1, CMS Contact: J. Truesdale
!         modified by A. Gettelman and C. Craig Nov 2010 to separate micro/macro physics
! 
!-----------------------------------------------------------------------

   use shr_kind_mod,    only: r8 => shr_kind_r8
   use ppgrid
   use phys_buffer,     only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx, pbuf_get_fld_idx, pbuf_times
   use stratiform,      only: stratiform_tend
   use phys_control,    only: phys_getopts
   use microp_driver,   only: microp_driver_tend
   use macrop_driver,   only: macrop_driver_tend
   use physics_types,   only: physics_state, physics_tend, physics_ptend, physics_update, physics_ptend_init
   use cam_diagnostics, only: diag_conv_tend_ini, diag_phys_writeout, diag_conv, diag_export, diag_state_b4_phys_write
   use cam_history,     only: outfld
   use physconst,       only: cpair
   use constituents,    only: pcnst, qmin, cnst_get_ind
   use convect_deep,    only: convect_deep_tend, convect_deep_tend_2, deep_scheme_does_scav_trans
   use time_manager,    only: is_first_step, get_nstep
   use convect_shallow, only: convect_shallow_tend
   use check_energy,    only: check_energy_chng, check_energy_fix
   use check_energy,    only: check_tracers_data, check_tracers_init, check_tracers_chng
   use dycore,          only: dycore_is
   use aerosol_intr,    only: aerosol_wet_intr
   use camsrfexch_types,only: cam_out_t, cam_in_t
   use radiation,       only: radiation_tend
   use cloud_diagnostics, only: put_cloud_diagnostics
   use perf_mod
#ifdef MODAL_AERO
   use modal_aero_data, only: qneg3_worst_thresh_amode
#endif
   use mo_gas_phase_chemdr,only: map2chm
   use clybry_fam,         only: clybry_fam_adj
   use sslt_rebin,      only: sslt_rebin_adv
   use tropopause,      only: tropopause_output
   
   use spmd_utils,        only : masterproc
   use cam_logfile,       only : iulog
   use machine_learning_model, only: ml_model_tend
   use machine_learning_model_config
   use diag_ml,          only: ml_verify_bottom
   use export_state,    only: record_state_before, record_state_after

   implicit none

!
! Arguments
!
   real(r8), intent(in) :: ztodt                          ! 2 delta t (model time increment)
   real(r8), intent(in) :: pblht(pcols)                ! Planetary boundary layer height
   real(r8), intent(in) :: tpert(pcols)                ! Thermal temperature excess
   real(r8), intent(inout) :: qpert(pcols,pcnst)         ! Thermal humidity & constituent excess1
   real(r8), intent(inout) :: fsns(pcols)                   ! Surface solar absorbed flux
   real(r8), intent(inout) :: fsnt(pcols)                   ! Net column abs solar flux at model top
   real(r8), intent(inout) :: flns(pcols)                   ! Srf longwave cooling (up-down) flux
   real(r8), intent(inout) :: flnt(pcols)                   ! Net outgoing lw flux at model top
   real(r8), intent(inout) :: fsds(pcols)                   ! Surface solar down flux
   real(r8), intent(in) :: landm(pcols)                   ! land fraction ramp

   type(physics_state), intent(inout) :: state
   type(physics_tend ), intent(inout) :: tend
   type(pbuf_fld),      intent(inout) :: pbuf(pbuf_size_max)
   type(cam_out_t),     intent(inout) :: cam_out
   type(cam_in_t),      intent(in)    :: cam_in

!
!---------------------------Local workspace-----------------------------
!

   type(physics_state)   :: state_eq         ! equilibrium state variables
   type(physics_ptend)   :: ptend            ! indivdual parameterization tendencies

   integer :: nstep                          ! current timestep number
   integer      lat(pcols)                   ! current latitudes(indices)
   integer      lon(pcols)                   ! current longtitudes(indices)

   real(r8) :: net_flx(pcols)

   real(r8) :: zdu(pcols,pver)               ! detraining mass flux from deep convection
   real(r8) :: cmfmc(pcols,pverp)            ! Convective mass flux--m sub c

   real(r8) cmfcme(pcols,pver)                ! cmf condensation - evaporation
   real(r8) cmfmc2(pcols,pverp)               ! Moist convection cloud mass flux
   real(r8) coszrs(pcols)                     ! Cosine solar zenith angle
   real(r8) dlf(pcols,pver)                   ! Detraining cld H20 from shallow + deep convections
   real(r8) dlf2(pcols,pver)                  ! Detraining cld H20 from shallow convections
   real(r8) pflx(pcols,pverp)                 ! Conv rain flux thru out btm of lev
   real(r8) rtdt                              ! 1./ztodt

   integer lchnk                              ! chunk identifier
   integer ncol                               ! number of atmospheric columns

   integer  i,k,m                             ! Longitude, level, constituent indices
   integer :: ixcldice, ixcldliq              ! constituent indices for cloud liquid and ice water.
                                           

   real(r8) dellow(pcols)                     ! delta p for bottom three levels of model
   real(r8) tavg(pcols)                       ! mass weighted average temperature for 

! physics buffer fields to compute tendencies for stratiform package
   integer itim, ifld
   real(r8), pointer, dimension(:,:) :: cld        ! cloud fraction


! physics buffer fields for total energy and mass adjustment
   real(r8), pointer, dimension(:  ) :: teout
   real(r8), pointer, dimension(:,:) :: tini
   real(r8), pointer, dimension(:,:) :: qini
   real(r8), pointer, dimension(:,:) :: cldliqini
   real(r8), pointer, dimension(:,:) :: cldiceini
   real(r8), pointer, dimension(:,:) :: dtcore

   real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble

! convective precipitation variables
   real(r8) :: prec_zmc(pcols)                ! total precipitation from ZM convection
   real(r8) :: snow_zmc(pcols)                ! snow from ZM convection
   real(r8) :: prec_cmf(pcols)                ! total precipitation from Hack convection
   real(r8) :: snow_cmf(pcols)                ! snow from Hack convection

! stratiform precipitation variables
   real(r8) :: prec_str(pcols)    ! sfc flux of precip from stratiform (m/s)
   real(r8) :: snow_str(pcols)     ! sfc flux of snow from stratiform   (m/s)
   real(r8) :: prec_pcw(pcols)     ! total precip from prognostic cloud scheme
   real(r8) :: snow_pcw(pcols)     ! snow from prognostic cloud scheme
   real(r8) :: prec_sed(pcols)     ! total precip from cloud sedimentation
   real(r8) :: snow_sed(pcols)     ! snow from cloud ice sedimentation
   real(r8), pointer, dimension(:,:) :: cldo 

! energy checking variables
   real(r8) :: zero(pcols)                    ! array of zeros
   real(r8) :: rliq(pcols)                    ! vertical integral of liquid not yet in q(ixcldliq)
   real(r8) :: rliq2(pcols)                   ! vertical integral of liquid from shallow scheme
   real(r8) :: flx_cnd(pcols)
   real(r8) :: flx_heat(pcols)
   type(check_tracers_data):: tracerint             ! energy integrals and cummulative boundary fluxes
   real(r8) :: zero_tracers(pcols,pcnst)

!++ debug code to be removed after PBL scheme validation
   integer :: kmx
!-- debug code to be removed after PBL scheme validation

  real(r8)  :: cmeliq(pcols,pver)                      ! Rate of cond-evap of liq within the cloud
!  pass macro to micro

   real(r8) :: prect(pcols)     !
   real(r8) :: precc(pcols)     !
   real(r8) :: precl(pcols)     !
   real(r8) :: precsc(pcols)     !
   real(r8) :: precsl(pcols)     !
   real(r8) :: flds (pcols)                 ! downwelling lw flux at the surface
   real(r8) :: solin(pcols)                 ! incoming toa solar down flux

   real(r8), dimension(pcols) :: srfrad
   real(r8), dimension(pcols) :: soll
   real(r8), dimension(pcols) :: sols
   real(r8), dimension(pcols) :: solld
   real(r8), dimension(pcols) :: solsd


   type(physics_state) :: initial_state
   type(physics_tend ) :: initial_tend
   real(r8) relhum(pcols,pver) ! temporary workspace

   character(len=16) :: microp_scheme 
   call phys_getopts( microp_scheme_out = microp_scheme )

!-----------------------------------------------------------------------
   call t_startf('bc_init')

   zero = 0._r8
   zero_tracers(:,:) = 0._r8

   lchnk = state%lchnk
   ncol  = state%ncol

   rtdt = 1._r8/ztodt

   nstep = get_nstep()


! Associate pointers with physics buffer fields
   itim = pbuf_old_tim_idx()
   ifld = pbuf_get_fld_idx('CLD')
   cld  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

   ifld   = pbuf_get_fld_idx('TEOUT')
   teout  => pbuf(ifld)%fld_ptr(1,1:pcols,1,lchnk,itim)
   ifld   =  pbuf_get_fld_idx('QINI')
   qini   => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
   ifld   =  pbuf_get_fld_idx('CLDLIQINI')
   cldliqini => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
   ifld   =  pbuf_get_fld_idx('CLDICEINI')
   cldiceini => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
   ifld   =  pbuf_get_fld_idx('TINI')
   tini   => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1)
   ifld   =  pbuf_get_fld_idx('DTCORE')
   dtcore => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, itim)

   ifld    = pbuf_get_fld_idx('FRACIS')
   fracis  => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk, 1:pcnst)
   if(is_first_step())then
      ifld = pbuf_get_fld_idx('CLDO')
      do m=1,pbuf_times
         cldo => pbuf(ifld)%fld_ptr(1,1:pcols,1:pver,lchnk,m)
         cldo(:ncol,:) = 0
      enddo
   endif
!
! Set physics tendencies to 0
   tend %dTdt(:ncol,:pver)  = 0._r8
   tend %dudt(:ncol,:pver)  = 0._r8
   tend %dvdt(:ncol,:pver)  = 0._r8

   call physics_ptend_init (ptend) ! Initialize parameterization tendency structure
!
! Make sure that input tracers are all positive (probably unnecessary)
!
    
#ifdef MODAL_AERO
   call qneg3_modalx1( &
              'TPHYSBCb',lchnk  ,ncol    ,pcols   ,pver    , &
              1, pcnst, qmin  ,state%q, qneg3_worst_thresh_amode )
#else
   call qneg3('TPHYSBCb',lchnk  ,ncol    ,pcols   ,pver    , &
              1, pcnst, qmin  ,state%q )
#endif

   call clybry_fam_adj( ncol, lchnk, map2chm, state%q )

   fracis (:ncol,:,1:pcnst) = 1._r8
!
! Dump out "before physics" state
!
   call diag_state_b4_phys_write (state)

! compute mass integrals of input tracers state
   call check_tracers_init(state, tracerint)

   call t_stopf('bc_init')

!===================================================
! Global mean total energy fixer
!===================================================
   call t_startf('energy_fixer')

   !*** BAB's FV heating kludge *** save the initial temperature
   tini(:ncol,:pver) = state%t(:ncol,:pver)
   if (dycore_is('LR')) then
      call check_energy_fix(state, ptend, nstep, flx_heat)
      call physics_update(state, tend, ptend, ztodt)
      call check_energy_chng(state, tend, "chkengyfix", nstep, ztodt, zero, zero, zero, flx_heat)
   end if
! Save state for convective tendency calculations.
   call diag_conv_tend_ini(state,pbuf)

   call cnst_get_ind('CLDLIQ', ixcldliq)
   call cnst_get_ind('CLDICE', ixcldice)
   qini     (:ncol,:pver) = state%q(:ncol,:pver,       1)
   cldliqini(:ncol,:pver) = state%q(:ncol,:pver,ixcldliq)
   cldiceini(:ncol,:pver) = state%q(:ncol,:pver,ixcldice)

   call outfld('TEOUT', teout       , pcols, lchnk   )
   call outfld('TEINP', state%te_ini, pcols, lchnk   )
   call outfld('TEFIX', state%te_cur, pcols, lchnk   )

! set and output the dse change due to dynpkg
   if( nstep > pbuf_times-1 ) then
      do k = 1,pver
         dtcore(:ncol,k) = (state%s(:ncol,k) - dtcore(:ncol,k))/(cpair*ztodt)
      end do
      call outfld( 'DTCORE', dtcore, pcols, lchnk )
   end if

   call t_stopf('energy_fixer')
!
!===================================================
! Dry adjustment
! This code block is not a good example of interfacing a parameterization
!===================================================
   call t_startf('dry_adjustment')

! Copy state info for input to dadadj
! This is a kludge, so that dadadj does not have to be correctly reformulated in dry static energy

   ptend%s(:ncol,:pver)   = state%t(:ncol,:pver)
   ptend%q(:ncol,:pver,1) = state%q(:ncol,:pver,1)

   call dadadj (lchnk, ncol, state%pmid,  state%pint,  state%pdel,  &
                ptend%s, ptend%q(1,1,1))
   ptend%name  = 'dadadj'
   ptend%ls    = .TRUE.
   ptend%lq(1) = .TRUE.
   ptend%s(:ncol,:)   = (ptend%s(:ncol,:)   - state%t(:ncol,:)  )/ztodt * cpair
   ptend%q(:ncol,:,1) = (ptend%q(:ncol,:,1) - state%q(:ncol,:,1))/ztodt
   call physics_update (state, tend, ptend, ztodt)

   call t_stopf('dry_adjustment')

   
!
!===================================================
! ML model added code - for training purposes
! We repeat here the same call as in tphysbc_param.F90
! so that we can validate the inputs and outputs of
! a ML model.
!===================================================
   call record_state_before(state, cam_in, tend, initial_state, initial_tend, relhum)

!
!===================================================
! Moist convection
!===================================================
   !DF:CONV: call t_startf('moist_convection')
!
! Since the PBL doesn't pass constituent perturbations, they
! are zeroed here for input to the moist convection routine
!
   !DF:CONV: call t_startf ('convect_deep_tend')
   !DF:CONV: call convect_deep_tend(  prec_zmc,   &
   !DF:CONV:      pblht,    cmfmc,      cmfcme,             &
   !DF:CONV:      tpert,    dlf,        pflx,    zdu,       &
   !DF:CONV:      rliq,    &
   !DF:CONV:      ztodt,    snow_zmc,  &
   !DF:CONV:      state,   ptend, cam_in%landfrac, pbuf ) 
   !DF:CONV: call t_stopf('convect_deep_tend')

   !DF:CONV: call physics_update(state, tend, ptend, ztodt)

! Check energy integrals, including "reserved liquid"
   !DF:CONV: flx_cnd(:ncol) = prec_zmc(:ncol) + rliq(:ncol)
   !DF:CONV: call check_energy_chng(state, tend, "convect_deep", nstep, ztodt, zero, flx_cnd, snow_zmc, zero)

!
! Call Hack (1994) convection scheme to deal with shallow/mid-level convection
!
   !DF:CONV SH: call t_startf ('convect_shallow_tend')

   !DF:CONV SH: call convect_shallow_tend (ztodt   ,&
   !DF:CONV SH:      qpert     ,   &
   !DF:CONV SH:      pblht      ,     &
   !DF:CONV SH:      cmfmc      ,cmfmc2  ,  prec_cmf,   &
   !DF:CONV SH:      dlf        , dlf2,   rliq      , rliq2, & 
   !DF:CONV SH:      snow_cmf   , state, ptend,  pbuf       )
   !DF:CONV SH: call t_stopf ('convect_shallow_tend')
   
   !DF:CONV SH: call physics_update (state, tend, ptend, ztodt)

   !DF:CONV SH: flx_cnd(:ncol) = prec_cmf(:ncol) + rliq2(:ncol)
   !DF:CONV SH: call check_energy_chng(state, tend, "convect_shallow", nstep, ztodt, zero, flx_cnd, snow_cmf, zero)

   !DF:CONV SH: call check_tracers_chng(state, tracerint, "convect_shallow", nstep, ztodt, zero_tracers)

   !DF:CONV SH: call t_stopf('moist_convection')

   ! Rebin the 4-bin version of sea salt into bins for coarse and accumulation
   ! modes that correspond to the available optics data.  This is only necessary
   ! for CAM-RT.  But it's done here so that the microphysics code which is called
   ! from the stratiform interface has access to the same aerosols as the radiation
   ! code.
   !DF:STRATI: call sslt_rebin_adv(pbuf, state)

!if CAM4/RK microphysics

   !DF:STRATI: if( microp_scheme .eq. 'RK' ) then

!===================================================
! Calculate stratiform tendencey (sedimentation, detrain, cloud fraction and microphysics )
!===================================================
   !DF:STRATI:    call t_startf('stratiform_tend')

   !DF:STRATI:    call stratiform_tend(state, ptend, ztodt, &
   !DF:STRATI:         cam_in%icefrac, cam_in%landfrac, cam_in%ocnfrac, &
   !DF:STRATI:         landm, cam_in%snowhland, & ! sediment
   !DF:STRATI:         dlf, dlf2, & ! detrain
   !DF:STRATI:         rliq  , & ! check energy after detrain
   !DF:STRATI:         cmfmc,   cmfmc2, &
   !DF:STRATI:         cam_in%ts,      cam_in%sst,        zdu,  &
   !DF:STRATI:         prec_str, snow_str, prec_sed, snow_sed, prec_pcw, snow_pcw, & 
   !DF:STRATI:         pbuf, state_eq)

   !DF:STRATI:    call physics_update (state, tend, ptend, ztodt)
   !DF:STRATI:    call check_energy_chng(state, tend, "cldwat_tend", nstep, ztodt, zero, prec_str, snow_str, zero)
   
   !DF:STRATI:    call t_stopf('stratiform_tend')

   !DF:STRATI:  elseif( microp_scheme .eq. 'MG' ) then

!===================================================
! Calculate macrophysical tendencey (sedimentation, detrain, cloud fraction)
!===================================================

   !DF:STRATI:    call t_startf('macrop_tend')
    
   !DF:STRATI:    call macrop_driver_tend(state, ptend, ztodt, &
   !DF:STRATI:         cam_in%landfrac, cam_in%ocnfrac, &
   !DF:STRATI:         cam_in%snowhland, & ! sediment
   !DF:STRATI:         dlf, dlf2, & ! detrain
   !DF:STRATI:         cmfmc,   cmfmc2, &
   !DF:STRATI:         cam_in%ts,      cam_in%sst,        zdu,  &
   !DF:STRATI:          pbuf, state_eq, cmeliq)

!      call physics_update (state, tend, ptend, ztodt)
!      call check_energy_chng(state, tend, "macrop_tend", nstep, ztodt, zero, prec_str, snow_str, zero)

   !DF:STRATI:    call t_stopf('macrop_tend')

!===================================================
! Calculate cloud microphysics 
!===================================================

   !DF:STRATI:    call t_startf('microp_tend')

   !DF:STRATI:    call microp_driver_tend(state, ptend, ztodt, &
#ifdef MODAL_AERO
   !DF:STRATI:         cam_in%cflx, & ! constituent sources
#endif
   !DF:STRATI:         rliq  , & ! check energy after detrain
   !DF:STRATI:         prec_str, snow_str, prec_sed, snow_sed, prec_pcw, snow_pcw, & 
   !DF:STRATI:         pbuf, state_eq, cmeliq)

   !DF:STRATI:    call physics_update (state, tend, ptend, ztodt)
   !DF:STRATI:    call check_energy_chng(state, tend, "microp_tend", nstep, ztodt, zero, prec_str, snow_str, zero)

   !DF:STRATI:    call t_stopf('microp_tend')

   !DF:STRATI: endif

!end microphysics conditional.

   !DF:SCAVENG: if ( .not. deep_scheme_does_scav_trans() ) then

      !===================================================
      !  Aerosol wet chemistry determines scavenging fractions, and transformations
      !
      !
      !  Then do convective transport of all trace species except water vapor and
      !     cloud liquid and ice (we needed to do the scavenging first
      !     to determine the interstitial fraction) 
      !===================================================

   !DF:SCAVENG:    call t_startf('bc_aerosols')
   !DF:SCAVENG:    call aerosol_wet_intr (state, ptend, ztodt, pbuf, cam_out, dlf)
   !DF:SCAVENG:    call physics_update (state, tend, ptend, ztodt)

   !DF:SCAVENG:    call t_startf ('convect_deep_tend2')
   !DF:SCAVENG:    call convect_deep_tend_2( state,   ptend,  ztodt,  pbuf ) 
   !DF:SCAVENG:    call t_stopf ('convect_deep_tend2')

   !DF:SCAVENG:    call physics_update (state, tend, ptend, ztodt)

      ! check tracer integrals
   !DF:SCAVENG:    call check_tracers_chng(state, tracerint, "cmfmca", nstep, ztodt, ptend%cflx_srf)

   !DF:SCAVENG:    call t_stopf('bc_aerosols')

   !DF:SCAVENG: endif

   !===================================================
   ! Moist physical parameteriztions complete: 
   ! send dynamical variables, and derived variables to history file
   !===================================================

   !DF:CLD-DIAG: call t_startf('bc_history_write')
   !DF:CLD-DIAG: call diag_phys_writeout(state, cam_out%psl)
   !DF:CLD-DIAG: call diag_conv(state, ztodt,    &
   !DF:CLD-DIAG:      prec_zmc, snow_zmc, prec_cmf, snow_cmf, prec_sed, snow_sed, prec_pcw, snow_pcw)

   !DF:CLD-DIAG: call t_stopf('bc_history_write')

   !===================================================
   ! Write cloud diagnostics on history file
   !===================================================

   !DF:CLD-DIAG: call t_startf('bc_cld_diag_history_write')

   !DF:CLD-DIAG: call put_cloud_diagnostics(state, pbuf)

   !DF:CLD-DIAG: call t_stopf('bc_cld_diag_history_write')

   !===================================================
   ! ML code computations
   !===================================================
   call t_startf ('ml_model_tend')
   !if( masterproc ) write(iulog,*) 'tphysbc_ml about to call ml_model_tend'

   call ml_model_tend ( ztodt, state, ptend,  pbuf, &
      fsns,    fsnt,    flns,    flnt, fsds, flds, solin, &
      cam_in, cam_out, &
      prect, precc, precl, &
      precsc, precsl, & 
      srfrad, soll, solld, sols, solsd, &
      relhum )


   call physics_update (state, tend, ptend, ztodt)
   call check_energy_chng(state, tend, "mlmodel", nstep, ztodt, zero, zero, zero, net_flx)
   call t_stopf('ml_model_tend')

   call output_state(state, cam_in, &
              prect, precc, precl, precsc, precsl, &
              fsns, fsnt, flns, flnt, fsds, flds, solin, &
              srfrad, soll, solld, sols, solsd & 
              )
  
   call diag_phys_writeout(state, cam_out%psl)

   !===================================================
   ! Radiation computations
   !===================================================
   
   if (.not. is_radiation_from_ml()) then
      !call t_startf('radiation')

      !call radiation_tend(state,ptend,pbuf, &
      !     cam_out, cam_in, &
      !     cam_in%landfrac,landm,cam_in%icefrac, cam_in%snowhland, &
      !     fsns,    fsnt, flns,    flnt,  &
      !     fsds, net_flx)
      !! Set net flux used by spectral dycores
      !do i=1,ncol
      !   tend%flx_net(i) = net_flx(i)
      !end do
      !call physics_update(state, tend, ptend, ztodt)
      !call check_energy_chng(state, tend, "radheat", nstep, ztodt, zero, zero, zero, net_flx)

      !call t_stopf('radiation')
   end if !(.not. is_radiation_from_ml())

   ! Diagnose the location of the tropopause and its location to the history file(s).
   call t_startf('tropopause')
   call tropopause_output(state)
   call t_stopf('tropopause')

   ! Save atmospheric fields to force surface models
   !DF: SRFXFER: call t_startf('srfxfer')
   !DF: SRFXFER: call srfxfer (state,cam_out,prec_zmc,snow_zmc, &
   !DF: SRFXFER:      prec_cmf,snow_cmf,prec_sed,snow_sed, &
   !DF: SRFXFER:      prec_pcw,snow_pcw)
   !DF: SRFXFER: call t_stopf('srfxfer')

   ! Save atmospheric fields to force surface models
   call t_startf('srfxfer_ml')
   call srfxfer_ml (state,cam_out, &
           prect, precc, precl, &
           precsc, precsl&
           )
   call t_stopf('srfxfer_ml')
   !call ml_verify_bottom(state, cam_out, & 
   !        fsns, fsnt, flns, flnt, fsds, solin, .false.) 


   ! Write export state to history file
   call t_startf('diag_export')
   call diag_export(cam_out)
   call t_stopf('diag_export')

!
!===================================================
! ML model added code - for training purposes
! We repeat here the same call as in tphysbc_param.F90
! so that we can validate the inputs and outputs of
! a ML model.
!===================================================
   call record_state_after(state, tend, ztodt, initial_state, initial_tend)


end subroutine tphysbc_ml_main


subroutine srfxfer_ml(state,cam_out,prect,precc,precl, &
                precsc,precsl)
        !state,cam_out,prec_zmc,snow_zmc, &
        !           prec_cmf,snow_cmf,prec_sed,snow_sed, &
        !           prec_pcw,snow_pcw)

!-----------------------------------------------------------------------
!
! Purpose:
! Transfer atmospheric fields into necessary surface data structures
!
! Author: L. Bath  CMS Contact: M. Vertenstein
!
!-----------------------------------------------------------------------
   use shr_kind_mod,     only: r8 => shr_kind_r8
   use physics_types,    only: physics_state
   use ppgrid,           only: pver, pcols
   use cam_history,      only: outfld
   use comsrf,           only: psm1, srfrpdel, prcsnw
   use camsrfexch_types, only: cam_out_t
   use chem_surfvals,    only: chem_surfvals_get
   use co2_cycle,        only: co2_transport, c_i
   use physconst,        only: mwdry, mwco2
   use constituents,     only: pcnst
   use cam_control_mod,  only: rair
   
   implicit none

   !------------------------------Arguments--------------------------------
   !
   ! Input arguments
   !
   type(physics_state),  intent(in)    :: state
   type (cam_out_t),     intent(inout) :: cam_out
   real(r8), intent(in) :: prect(pcols)                ! total precipitation   from ZM convection
   real(r8), intent(in) :: precc(pcols)                ! total precipitation   from ZM convection
   real(r8), intent(in) :: precl(pcols)                ! total precipitation   from ZM convection
   real(r8), intent(in) :: precsc(pcols)                ! total precipitation   from ZM convection
   real(r8), intent(in) :: precsl(pcols)                ! total precipitation   from ZM convection

   !
   !---------------------------Local variables-----------------------------
   !
   integer :: i              ! Longitude index
   integer :: m              ! constituent index
   integer :: lchnk          ! Chunk index
   integer :: ncol
   !-----------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol
   do i=1,ncol
      cam_out%tbot(i)  = state%t(i,pver)
      cam_out%thbot(i) = state%t(i,pver) * state%exner(i,pver)
      cam_out%zbot(i)  = state%zm(i,pver)
      cam_out%ubot(i)  = state%u(i,pver)
      cam_out%vbot(i)  = state%v(i,pver)
      cam_out%pbot(i)  = state%pmid(i,pver)
      cam_out%rho(i)   = cam_out%pbot(i)/(rair*cam_out%tbot(i))
      cam_out%netsw(i) = cam_out%srfrad(i) - cam_out%flwds(i)
      psm1(i,lchnk)    = state%ps(i)
      srfrpdel(i,lchnk)= state%rpdel(i,pver)
   end do
   do m = 1, pcnst
     do i = 1, ncol
        cam_out%qbot(i,m) = state%q(i,pver,m)
     end do
   end do

   cam_out%co2diag(:ncol) = chem_surfvals_get('CO2VMR') * 1.0e+6_r8
   if (co2_transport()) then
      do i=1,ncol
         cam_out%co2prog(i) = state%q(i,pver,c_i(4)) * 1.0e+6_r8 *mwdry/mwco2
      end do
   end if
   !
   ! Precipation and snow rates from shallow convection, deep convection and stratiform processes.
   ! Compute total convective and stratiform precipitation and snow rates
   !
   do i=1,ncol
      cam_out%precc (i) = precc(i)
      cam_out%precl (i) = precl(i)
      cam_out%precsc(i) = precsc(i)
      cam_out%precsl(i) = precsl(i)

      ! jrm These checks should not be necessary if they exist in the parameterizations
      if (cam_out%precc(i) .lt.0._r8) cam_out%precc(i)=0._r8
      if (cam_out%precl(i) .lt.0._r8) cam_out%precl(i)=0._r8
      if (cam_out%precsc(i).lt.0._r8) cam_out%precsc(i)=0._r8
      if (cam_out%precsl(i).lt.0._r8) cam_out%precsl(i)=0._r8

      if (cam_out%precsc(i).gt.cam_out%precc(i)) cam_out%precsc(i)=cam_out%precc(i)
      if (cam_out%precsl(i).gt.cam_out%precl(i)) cam_out%precsl(i)=cam_out%precl(i)
      ! end jrm
   end do

   ! total snowfall rate: needed by slab ocean model
   prcsnw(:ncol,lchnk) = cam_out%precsc(:ncol) + cam_out%precsl(:ncol)


   !cam_out%bcphiwet = 0._r8
   !cam_out%bcphidry = 0._r8
   !cam_out%bcphodry = 0._r8
   !cam_out%ocphiwet = 0._r8
   !cam_out%ocphidry = 0._r8
   !cam_out%ocphodry = 0._r8
   !cam_out%dstwet1 = 0._r8
   !cam_out%dstdry1 = 0._r8
   !cam_out%dstwet2 = 0._r8
   !cam_out%dstdry2 = 0._r8
   !cam_out%dstwet3 = 0._r8
   !cam_out%dstdry3 = 0._r8
   !cam_out%dstwet4 = 0._r8
   !cam_out%dstdry4 = 0._r8

   !call ml_verify_bottom(state, cam_out, .false.) 

end subroutine srfxfer_ml


subroutine  output_state(state, cam_in, &
                prect, precc, precl, precsc, precsl, &
                fsns, fsnt, flns, flnt, fsds, flds, solin, &
                srfrad, soll, solld, sols, solsd )
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   use shr_kind_mod,    only: r8 => shr_kind_r8
   use ppgrid
   use physconst,     only: cpair
   use tidal_diag,    only: get_tidal_coeffs
   use constituents,  only: pcnst, cnst_name, cnst_longname, cnst_cam_outfld, ptendnam, dmetendnam, apcnst, bpcnst, &
                         cnst_get_ind
   use cam_history,     only: outfld
   use camsrfexch_types,only: cam_out_t, cam_in_t
   use physics_types,   only: physics_state, physics_tend, physics_ptend, physics_update, physics_ptend_init
   use machine_learning_model_config



! Arguments:

   !real(r8),            intent(in) :: ztodt   ! timestep for computing physics tendencies
   type(physics_state), intent(in) :: state
   type(cam_in_t),      intent(in) :: cam_in

! convective precipitation variables
   real(r8), intent(in) :: prect(pcols)                ! total precipitation   from ZM convection
   real(r8), intent(in) :: precc(pcols)                ! snow from ZM   convection
   real(r8), intent(in) :: precl(pcols)                ! total precipitation   from Hack convection
   real(r8), intent(in) :: precsc(pcols)               ! snow from   Hack   convection
   real(r8), intent(in) :: precsl(pcols)               ! total precipitation   from ZM convection

   real(r8), intent(in) :: fsns(pcols)                 ! Surface solar absorbed flux
   real(r8), intent(in) :: fsnt(pcols)                 ! Net column abs solar flux at model top
   real(r8), intent(in) :: flns(pcols)                 ! Srf longwave cooling (up-down) flux
   real(r8), intent(in) :: flnt(pcols)                 ! Net outgoing lw flux at model top
   real(r8), intent(in) :: fsds(pcols)                 ! Surface solar down flux
   real(r8), intent(in) :: flds(pcols)                 ! Surface lw down flux
   real(r8), intent(in) :: solin(pcols)                !

   real(r8), intent(in) :: srfrad(pcols)
   real(r8), intent(in) :: soll(pcols)
   real(r8), intent(in) :: sols(pcols)
   real(r8), intent(in) :: solld(pcols)
   real(r8), intent(in) :: solsd(pcols)


   integer :: i, k, m, lchnk, ncol

   lchnk = state%lchnk
   ncol  = state%ncol


   call outfld('PRECC   ', precc, pcols, lchnk )
   call outfld('PRECL   ', precl, pcols, lchnk )
   call outfld('PRECSC  ', precsc, pcols, lchnk )
   call outfld('PRECSL  ', precsl, pcols, lchnk )
   call outfld('PRECT   ', prect, pcols, lchnk )

   !call outfld('T       ',state%t , pcols   ,lchnk   )
   !call outfld('U       ',state%u , pcols   ,lchnk   )
   !call outfld('V       ',state%v , pcols   ,lchnk   )

   !call outfld('PS      ',state%ps, pcols   ,lchnk   )

   do m=1,pcnst
      if ( cnst_cam_outfld(m) ) then
         call outfld(cnst_name(m),state%q(1,1,m),pcols ,lchnk )
      end if
   end do

   call outfld('SHFLX',    cam_in%shf,       pcols, lchnk)
   call outfld('LHFLX',    cam_in%lhf,       pcols, lchnk)

   !call outfld('PTTEND',ftem3, pcols, lchnk )

   if (is_radiation_from_ml()) then
      call outfld('SOLIN   ',solin ,pcols,lchnk)
      call outfld('FSDS    ',fsds  ,pcols,lchnk)
      call outfld('FSNT    ',fsnt  ,pcols,lchnk)
      call outfld('FSNS    ',fsns  ,pcols,lchnk)

      call outfld('FLNS    ',flns  ,pcols,lchnk)
      call outfld('FLNT    ',flnt  ,pcols,lchnk)
      call outfld('FLDS    ',flds  ,pcols,lchnk)
      call outfld('SRFRAD  ',srfrad,pcols,lchnk)
      call outfld('SOLL    ',soll ,pcols,lchnk)
      call outfld('SOLLD   ',solld,pcols,lchnk)
      call outfld('SOLS    ',sols ,pcols,lchnk)
      call outfld('SOLSD   ',solsd,pcols,lchnk)
   end if

end subroutine output_state


end module tphysbc_ml

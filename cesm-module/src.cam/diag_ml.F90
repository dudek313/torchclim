module diag_ml
   use time_manager,    only: is_first_step
   use shr_kind_mod,      only : r8=>shr_kind_r8
   use phys_grid,       only: get_lat_p, get_lon_p 


   implicit none
   private
   save

   public :: &
             ml_diagnostics, &   !
             ml_diagnostics_pre, &   !
             ml_diagnostics_tend_post, &
             ml_verify_bottom, &
             ml_diagnostics_post_tend, &
             ml_diagnostics_pre_prognostic_vars


     !real(r8), dimension(1) :: lats = (/ -6._r8 /)
     !real(r8), dimension(1) :: lons = (/ 180._r8 /)

     !logical                  :: do_init_location = .true.
     !integer, dimension(1)    :: cols
     !integer, dimension(1)    :: chunks
     !integer, dimension(1)    :: owners
     !integer :: col

     logical             :: found_location = .false.
     integer             :: the_col = -1
     integer             :: chunk = -1


  contains

subroutine ml_diagnostics_tend_post(state, solin, pttend, pteq)

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
   !use machine_learning_model_config
   use cam_logfile,       only : iulog
   use abortutils,        only : endrun
   use spmd_utils,         only: masterproc

   implicit none

!
! Arguments
!
   type(physics_state), intent(in) :: state
   real(r8), dimension(pcols,pver) :: pteq
   real(r8), dimension(pcols,pver) :: pttend
   real(r8), dimension(pcols)      :: solin
   integer :: lchnk, col



   return


    col = find_coord(state)
    if( col < 0 ) then
        return
    end if

    write(iulog,*) '*************************************************'
    write(iulog,*) 'ml_model_tend predict_model solin  : \n', solin(col)
    write(iulog,*) 'ml_model_tend predict_model pttend  : \n', pttend(col,:)
    write(iulog,*) 'ml_model_tend predict_model pteq  : \n', pteq(col,:)
    write(iulog,*) '*************************************************'


 

end subroutine ml_diagnostics_tend_post



subroutine ml_diagnostics_pre (ztodt,   pblht,   tpert,   qpert,            &
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
   !use machine_learning_model_config
   use cam_logfile,       only : iulog
   use abortutils,        only : endrun
   use spmd_utils,         only: masterproc

   implicit none

!
! Arguments
!
   real(r8), intent(in) :: ztodt                          ! 2 delta t (model time increment)
   real(r8), intent(in) :: pblht(pcols)                ! Planetary boundary layer height
   real(r8), intent(in) :: tpert(pcols)                ! Thermal temperature excess
   real(r8), intent(in) :: qpert(pcols,pcnst)         ! Thermal humidity & constituent excess1
   real(r8), intent(in) :: fsns(pcols)                   ! Surface solar absorbed flux
   real(r8), intent(in) :: fsnt(pcols)                   ! Net column abs solar flux at model top
   real(r8), intent(in) :: flns(pcols)                   ! Srf longwave cooling (up-down) flux
   real(r8), intent(in) :: flnt(pcols)                   ! Net outgoing lw flux at model top
   real(r8), intent(in) :: fsds(pcols)                   ! Surface solar down flux
   real(r8), intent(in) :: landm(pcols)                   ! land fraction ramp

   !real(r8), intent(inout) :: solin(pcols)                 ! incoming toa solar down flux


   type(physics_state), intent(in) :: state
   type(physics_tend ), intent(in) :: tend
   type(pbuf_fld),      intent(in) :: pbuf(pbuf_size_max)
   type(cam_out_t),     intent(in) :: cam_out
   type(cam_in_t),      intent(in)    :: cam_in
   integer :: col




   return




   !if( .not. found_location ) then
   !   call init_coord(state)
      !do_init_location = .false. 
      !call get_chunk_coord_owner_p(1, lons, lats, chunks, cols, owners)
      !col = cols(1)
   !end if

    !lchnk = state%lchnk
    !ncol  = state%ncol 

    !if( masterproc ) then
    !if( found_location .and. (state%lchnk .eq. chunk) ) then
    
    col = find_coord(state)
    if( col > -1 ) then
      write(iulog,*) 'ml_diagnostics: start 0'
      write(iulog,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   
      write(iulog,*) 'tend%dtdt : ', tend%dtdt(col,:) 
      write(iulog,*) 'tend%dudt : ', tend%dudt(col,:) 
      write(iulog,*) 'tend%dvdt : ', tend%dvdt(col,:) 
    
      write(iulog,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(iulog,*) 'ml_diagnostics: end 0'

      write(iulog,*) 'ml_diagnostics: start 1'
      write(iulog,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   
      write(iulog,*) 'state%ps : ', state%ps(col) 
      write(iulog,*) 'state%t : ', state%t(col,:) 
      write(iulog,*) 'state%s : ', state%s(col,:) 
      write(iulog,*) 'state%u : ', state%u(col,:) 
      write(iulog,*) 'state%v : ', state%v(col,:) 
      write(iulog,*) 'state%omega : ', state%omega(col,:) 
      write(iulog,*) 'state%q : ', state%q(col,:,:)
      write(iulog,*) '...................................................................................'
      write(iulog,*) 'state%zm : ', state%zm(col,:)
      write(iulog,*) '...................................................................................'
      write(iulog,*) 'state%pint : ', state%pint(col,:)
      write(iulog,*) '...................................................................................'
      write(iulog,*) 'state%pmid : ', state%pmid(col,:)
      write(iulog,*) '...................................................................................'
      write(iulog,*) 'state%pdel : ', state%pdel(col,:)
      write(iulog,*) '...................................................................................'
      write(iulog,*) 'state%phis : ', state%phis(col)
      write(iulog,*) 'state%lat : ', state%lat(col) * 180._r8/3.14159265359_r8
      write(iulog,*) 'state%lon : ', state%lon(col) * 180._r8/3.14159265359_r8
    
      write(iulog,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(iulog,*) 'ml_diagnostics: end 1'

      write(iulog,*) 'ml_diagnostics: start 2'
      write(iulog,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
   
      write(iulog,*) 'cam_out%netsw : ', cam_out%netsw(col) 
      write(iulog,*) 'cam_out%srfrad : ', cam_out%srfrad(col) 
      write(iulog,*) 'cam_out%flwds : ', cam_out%flwds(col) 
      write(iulog,*) 'cam_out%precsc : ', cam_out%precsc(col) 
      write(iulog,*) 'cam_out%precsl : ', cam_out%precsl(col) 
      write(iulog,*) 'cam_out%precc : ', cam_out%precc(col) 
      write(iulog,*) 'cam_out%precl : ', cam_out%precl(col) 
      write(iulog,*) 'cam_out%soll : ', cam_out%soll(col) 
      write(iulog,*) 'cam_out%sols : ', cam_out%sols(col) 
      write(iulog,*) 'cam_out%solld : ', cam_out%solld(col) 
      write(iulog,*) 'cam_out%solsd : ', cam_out%solsd(col) 
      write(iulog,*) '...................................................................................'
      write(iulog,*) 'cam_out%tbot : ', cam_out%tbot(col) 
      write(iulog,*) 'cam_out%zbot : ', cam_out%zbot(col) 
      write(iulog,*) 'cam_out%ubot : ', cam_out%ubot(col) 
      write(iulog,*) 'cam_out%vbot : ', cam_out%vbot(col) 
      write(iulog,*) 'cam_out%qbot : ', cam_out%qbot(col, :) 
      write(iulog,*) 'cam_out%pbot : ', cam_out%pbot(col) 
      write(iulog,*) 'cam_out%rho  : ', cam_out%rho(col) 
      write(iulog,*) 'cam_out%thbot: ', cam_out%thbot(col) 
      !write(iulog,*) 'cam_out% : ', cam_out% 
      !write(iulog,*) 'cam_out% : ', cam_out% 
    
      write(iulog,*) '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
      write(iulog,*) 'ml_diagnostics: end 2'



      !call endrun('ml_diagnostics: printing diag and exiting...')
    end if

  end subroutine ml_diagnostics_pre



subroutine ml_diagnostics (ztodt,   pblht,   tpert,   qpert,            &
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
   !use machine_learning_model_config
   use cam_logfile,       only : iulog
   use abortutils,        only : endrun
   use spmd_utils,         only: masterproc

   implicit none

!
! Arguments
!
   real(r8), intent(in) :: ztodt                          ! 2 delta t (model time increment)
   real(r8), intent(in) :: pblht(pcols)                ! Planetary boundary layer height
   real(r8), intent(in) :: tpert(pcols)                ! Thermal temperature excess
   real(r8), intent(in) :: qpert(pcols,pcnst)         ! Thermal humidity & constituent excess1
   real(r8), intent(in) :: fsns(pcols)                   ! Surface solar absorbed flux
   real(r8), intent(in) :: fsnt(pcols)                   ! Net column abs solar flux at model top
   real(r8), intent(in) :: flns(pcols)                   ! Srf longwave cooling (up-down) flux
   real(r8), intent(in) :: flnt(pcols)                   ! Net outgoing lw flux at model top
   real(r8), intent(in) :: fsds(pcols)                   ! Surface solar down flux
   real(r8), intent(in) :: landm(pcols)                   ! land fraction ramp

   !real(r8), intent(inout) :: solin(pcols)                 ! incoming toa solar down flux


   type(physics_state), intent(in) :: state
   type(physics_tend ), intent(in) :: tend
   type(pbuf_fld),      intent(in) :: pbuf(pbuf_size_max)
   type(cam_out_t),     intent(in) :: cam_out
   type(cam_in_t),      intent(in)    :: cam_in
   integer :: col






   return








   !if( .not. found_location ) then
   !   call init_coord(state)
      !do_init_location = .false. 
      !call get_chunk_coord_owner_p(1, lons, lats, chunks, cols, owners)
      !col = cols(1)
   !end if

    !lchnk = state%lchnk
    !ncol  = state%ncol 

    !if( masterproc ) then
    !if( found_location .and. (state%lchnk .eq. chunk) ) then
    
    col = find_coord(state)
    if( col > -1 ) then
      write(iulog,*) 'ml_diagnostics: start 0'
      write(iulog,*) '===================================================================================='
   
      write(iulog,*) 'tend%dtdt : ', tend%dtdt(col,:) 
      write(iulog,*) 'tend%dudt : ', tend%dudt(col,:) 
      write(iulog,*) 'tend%dvdt : ', tend%dvdt(col,:) 
    
      write(iulog,*) '===================================================================================='
      write(iulog,*) 'ml_diagnostics: end 0'

      write(iulog,*) 'ml_diagnostics: start 1'
      write(iulog,*) '===================================================================================='
   
      write(iulog,*) 'state%ps : ', state%ps(col) 
      write(iulog,*) 'state%t : ', state%t(col,:) 
      write(iulog,*) 'state%s : ', state%s(col,:) 
      write(iulog,*) 'state%u : ', state%u(col,:) 
      write(iulog,*) 'state%v : ', state%v(col,:) 
      write(iulog,*) 'state%omega : ', state%omega(col,:) 
      write(iulog,*) 'state%q : ', state%q(col,:,:)
      write(iulog,*) '...................................................................................'
      write(iulog,*) 'state%zm : ', state%zm(col,:)
      write(iulog,*) '...................................................................................'
      write(iulog,*) 'state%pint : ', state%pint(col,:)
      write(iulog,*) '...................................................................................'
      write(iulog,*) 'state%pmid : ', state%pmid(col,:)
      write(iulog,*) '...................................................................................'
      write(iulog,*) 'state%pdel : ', state%pdel(col,:)
      write(iulog,*) '...................................................................................'
      write(iulog,*) 'state%phis : ', state%phis(col)
      write(iulog,*) 'state%lat : ', state%lat(col) * 180._r8/3.14159265359_r8
      write(iulog,*) 'state%lon : ', state%lon(col) * 180._r8/3.14159265359_r8
    
      write(iulog,*) '===================================================================================='
      write(iulog,*) 'ml_diagnostics: end 1'

      write(iulog,*) 'ml_diagnostics: start 2'
      write(iulog,*) '===================================================================================='
   
      write(iulog,*) 'cam_out%netsw : ', cam_out%netsw(col) 
      write(iulog,*) 'cam_out%srfrad : ', cam_out%srfrad(col) 
      write(iulog,*) 'cam_out%flwds : ', cam_out%flwds(col) 
      write(iulog,*) 'cam_out%precsc : ', cam_out%precsc(col) 
      write(iulog,*) 'cam_out%precsl : ', cam_out%precsl(col) 
      write(iulog,*) 'cam_out%precc : ', cam_out%precc(col) 
      write(iulog,*) 'cam_out%precl : ', cam_out%precl(col) 
      write(iulog,*) 'cam_out%soll : ', cam_out%soll(col) 
      write(iulog,*) 'cam_out%sols : ', cam_out%sols(col) 
      write(iulog,*) 'cam_out%solld : ', cam_out%solld(col) 
      write(iulog,*) 'cam_out%solsd : ', cam_out%solsd(col) 
      write(iulog,*) '...................................................................................'
      write(iulog,*) 'cam_out%tbot : ', cam_out%tbot(col) 
      write(iulog,*) 'cam_out%zbot : ', cam_out%zbot(col) 
      write(iulog,*) 'cam_out%ubot : ', cam_out%ubot(col) 
      write(iulog,*) 'cam_out%vbot : ', cam_out%vbot(col) 
      write(iulog,*) 'cam_out%qbot : ', cam_out%qbot(col, :) 
      write(iulog,*) 'cam_out%pbot : ', cam_out%pbot(col) 
      write(iulog,*) 'cam_out%rho  : ', cam_out%rho(col) 
      write(iulog,*) 'cam_out%thbot: ', cam_out%thbot(col) 
      !write(iulog,*) 'cam_out% : ', cam_out% 
      !write(iulog,*) 'cam_out% : ', cam_out% 
    
      write(iulog,*) '===================================================================================='
      write(iulog,*) 'ml_diagnostics: end 2'



      !call endrun('ml_diagnostics: printing diag and exiting...')
    end if

  end subroutine ml_diagnostics





subroutine ml_diagnostics_post_tend (ztodt,   pblht,   tpert,   qpert,            &
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
   use machine_learning_model_config, only: first_ml_time_step
   use cam_logfile,       only : iulog
   use abortutils,        only : endrun
   use spmd_utils,         only: masterproc

   implicit none

!
! Arguments
!
   real(r8), intent(in) :: ztodt                          ! 2 delta t (model time increment)
   real(r8), intent(in) :: pblht(pcols)                ! Planetary boundary layer height
   real(r8), intent(in) :: tpert(pcols)                ! Thermal temperature excess
   real(r8), intent(in) :: qpert(pcols,pcnst)         ! Thermal humidity & constituent excess1
   real(r8), intent(in) :: fsns(pcols)                   ! Surface solar absorbed flux
   real(r8), intent(in) :: fsnt(pcols)                   ! Net column abs solar flux at model top
   real(r8), intent(in) :: flns(pcols)                   ! Srf longwave cooling (up-down) flux
   real(r8), intent(in) :: flnt(pcols)                   ! Net outgoing lw flux at model top
   real(r8), intent(in) :: fsds(pcols)                   ! Surface solar down flux
   real(r8), intent(in) :: landm(pcols)                   ! land fraction ramp

   !real(r8), intent(inout) :: solin(pcols)                 ! incoming toa solar down flux


   type(physics_state), intent(in) :: state
   type(physics_tend ), intent(in) :: tend
   type(pbuf_fld),      intent(in) :: pbuf(pbuf_size_max)
   type(cam_out_t),     intent(in) :: cam_out
   type(cam_in_t),      intent(in) :: cam_in
   integer :: col, i, nstep
   integer lchnk                              ! chunk identifier
   integer ncol                               ! number of atmospheric columns
   real(r8) :: lon
   real(r8) :: lat

   real(r8), parameter :: pi = 3.14159265359_r8
   real(r8), parameter :: rad_to_deg = 180._r8/pi

   !return

   nstep = get_nstep()
   !
   !if( is_first_step() ) then
   !if( nstep .ne. 1 ) then
   if( nstep .ne. first_ml_time_step() ) then
       return
   end if

   lchnk = state%lchnk
   ncol  = state%ncol

   !nstep = get_nstep()

   do col=1,ncol
       lat = state%lat(col) * rad_to_deg
       lon = state%lon(col) * rad_to_deg

       call sleep(2) 
       write(iulog, '(A47,I,",",F,",",F,26(ES,:,","))') 'ml_diagnostics_post_tend: nstep lat lon, dtdt, ', nstep, lat, lon, tend%dtdt(col,:)
       call flush(iulog) 
   end do
  end subroutine ml_diagnostics_post_tend




subroutine ml_diagnostics_pre_prognostic_vars (ztodt,   pblht,   tpert,   qpert,            &
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
   use physconst,       only: cpair, rga
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
   use machine_learning_model_config, only: first_ml_time_step
   use ml_solin,        only : ml_calc_solin
   use cam_logfile,       only : iulog
   use abortutils,        only : endrun
   use spmd_utils,         only: masterproc

   implicit none

!
! Arguments
!
   real(r8), intent(in) :: ztodt                          ! 2 delta t (model time increment)
   real(r8), intent(in) :: pblht(pcols)                ! Planetary boundary layer height
   real(r8), intent(in) :: tpert(pcols)                ! Thermal temperature excess
   real(r8), intent(in) :: qpert(pcols,pcnst)         ! Thermal humidity & constituent excess1
   real(r8), intent(in) :: fsns(pcols)                   ! Surface solar absorbed flux
   real(r8), intent(in) :: fsnt(pcols)                   ! Net column abs solar flux at model top
   real(r8), intent(in) :: flns(pcols)                   ! Srf longwave cooling (up-down) flux
   real(r8), intent(in) :: flnt(pcols)                   ! Net outgoing lw flux at model top
   real(r8), intent(in) :: fsds(pcols)                   ! Surface solar down flux
   real(r8), intent(in) :: landm(pcols)                   ! land fraction ramp
   !real(r8), intent(inout) :: solin(pcols)                 ! incoming toa solar down flux


   type(physics_state), intent(in) :: state
   type(physics_tend ), intent(in) :: tend
   type(pbuf_fld),      intent(in) :: pbuf(pbuf_size_max)
   type(cam_out_t),     intent(in) :: cam_out
   type(cam_in_t),      intent(in) :: cam_in
   integer :: col, i, nstep
   integer lchnk                              ! chunk identifier
   integer ncol                               ! number of atmospheric columns
   real(r8) :: lon
   real(r8) :: lat
   
   real(r8) :: solin(pcols)                   ! Surface solar down flux



   real(r8), parameter :: pi = 3.14159265359_r8
   real(r8), parameter :: rad_to_deg = 180._r8/pi

   !return

   nstep = get_nstep()

   !if( is_first_step() ) then
   !if( nstep .ne. 1 ) then
   if( nstep .ne. first_ml_time_step() ) then
       return
   end if

   call ml_calc_solin(state, solin)

   lchnk = state%lchnk
   ncol  = state%ncol

   nstep = get_nstep()

   do col=1,ncol
       lat = state%lat(col) * rad_to_deg
       lon = state%lon(col) * rad_to_deg

       call sleep(2) 
       write(iulog, '(A47,I,",",F,",",F,",",26(ES,:,","))') 'ml_diagnostics_pre_prognostic_vars: q, ', nstep, lat, lon, state%q(col,:,1)
       call flush(iulog) 

       call sleep(2) 
       write(iulog, '(A47,I,",",F,",",F,",",26(ES,:,","))') 'ml_diagnostics_pre_prognostic_vars: temp, ', nstep, lat, lon, state%t(col,:)
       call flush(iulog) 

       call sleep(2) 
       write(iulog, '(A47,I,",",F,",",F,",",26(ES,:,","))') 'ml_diagnostics_pre_prognostic_vars: uwind, ', nstep, lat, lon, state%u(col,:)
       call flush(iulog) 

       call sleep(2) 
       write(iulog, '(A47,I,",",F,",",F,",",26(ES,:,","))') 'ml_diagnostics_pre_prognostic_vars: vwind, ', nstep, lat, lon, state%v(col,:)
       call flush(iulog) 

       call sleep(2) 
       write(iulog, '(A47,I,",",F,",",F,",",26(ES,:,","))') 'ml_diagnostics_pre_prognostic_vars: omega, ', nstep, lat, lon, state%omega(col,:)
       call flush(iulog) 

       call sleep(2) 
       write(iulog, '(A47,I,",",F,",",F,",",26(ES,:,","))') 'ml_diagnostics_pre_prognostic_vars: zm, ', nstep, lat, lon, &
       state%zm(col,:) + state%phis(col)*rga
       call flush(iulog) 

       call sleep(2) 
       write(iulog, '(A47,I,",",F,",",F,",",9(ES,:,","))') 'ml_diagnostics_pre_prognostic_vars: other, ', nstep, lat, lon, &
               state%ps(col), solin(col), cam_in%shf(col), cam_in%lhf(col), &
               fsns(col), flns(col), fsnt(col), flnt(col), fsds(col)
       call flush(iulog) 
   end do
  end subroutine ml_diagnostics_pre_prognostic_vars






  subroutine ml_verify_bottom(state, cam_out, &
                  fsns, fsnt, flns, flnt, fsds, solin, apply_fix)

     use phys_grid,         only: get_lat_p, get_lon_p 
     use physics_types,     only: physics_state
     use cam_logfile,       only: iulog
     use camsrfexch_types,  only: cam_out_t, cam_in_t
     use ppgrid,            only : pver, pcols, pverp

     implicit none

     type(physics_state), intent(in) :: state
     type(cam_out_t), intent(inout) :: cam_out
     logical, intent(in) :: apply_fix

     real(r8), intent(in) :: fsns(pcols)                   ! Surface solar absorbed flux
     real(r8), intent(in) :: fsnt(pcols)                   ! Net column abs solar flux at model top
     real(r8), intent(in) :: flns(pcols)                   ! Srf longwave cooling (up-down) flux
     real(r8), intent(in) :: flnt(pcols)                   ! Net outgoing lw flux at model top
     real(r8), intent(in) :: fsds(pcols)                   ! Surface solar down flux
     real(r8), intent(in) :: solin(pcols)                   ! Surface solar down flux


     integer lchnk                              ! chunk identifier
     integer ncol                               ! number of atmospheric columns
     integer i
     integer find_coord

     real(r8) :: lon
     real(r8) :: lat

     real(r8), parameter :: pi = 3.14159265359_r8
     real(r8), parameter :: rad_to_deg = 180._r8/pi
     real(r8), parameter :: max_ocean_temp = 271.5

     !return



     lchnk = state%lchnk
     ncol  = state%ncol

     do i=1,ncol
         lat = state%lat(i) * rad_to_deg 
         lon = state%lon(i) * rad_to_deg 

         !if(lat < 75._r8 .and. lat > -75._r8) then
         !if(lat < 65._r8 .and. lat > -65._r8) then
         !    cycle
         !end if

         !if(cam_out%tbot(i) <= -1.8_r8) then
         !if(cam_out%tbot(i) <= max_ocean_temp) then
         !    cycle
         !end if

         !write(iulog,*) 'df ml_verify_bottom, warn : ', lchnk, lat, lon, cam_out%tbot(i) 
         !if( apply_fix ) then
         !    cam_out%tbot(i) = 271.0_r8
         !    cam_out%thbot(i) = cam_out%tbot(i) * state%exner(i,pver)
         !end if
         write(iulog,*) 'df ml_verify_bottom, fsns   : ', lchnk, lat, lon, fsns(i)
         write(iulog,*) 'df ml_verify_bottom, fsnt   : ', lchnk, lat, lon, fsnt(i)
         write(iulog,*) 'df ml_verify_bottom, flns   : ', lchnk, lat, lon, flns(i)
         write(iulog,*) 'df ml_verify_bottom, flnt   : ', lchnk, lat, lon, flnt(i)
         write(iulog,*) 'df ml_verify_bottom, fsds   : ', lchnk, lat, lon, fsds(i)
         write(iulog,*) 'df ml_verify_bottom, solin  : ', lchnk, lat, lon, solin(i)


         write(iulog,*) 'df ml_verify_bottom, tbot   : ', lchnk, lat, lon, cam_out%tbot(i)
         write(iulog,*) 'df ml_verify_bottom, zbot   : ', lchnk, lat, lon, cam_out%zbot(i)
         write(iulog,*) 'df ml_verify_bottom, ubot   : ', lchnk, lat, lon, cam_out%ubot(i)
         write(iulog,*) 'df ml_verify_bottom, vbot   : ', lchnk, lat, lon, cam_out%vbot(i)
         write(iulog,*) 'df ml_verify_bottom, qbot   : ', lchnk, lat, lon, cam_out%qbot(i,:)
         write(iulog,*) 'df ml_verify_bottom, pbot   : ', lchnk, lat, lon, cam_out%pbot(i)
         write(iulog,*) 'df ml_verify_bottom, rho    : ', lchnk, lat, lon, cam_out%rho(i)
         write(iulog,*) 'df ml_verify_bottom, netsw  : ', lchnk, lat, lon, cam_out%netsw(i)
         write(iulog,*) 'df ml_verify_bottom, flwds  : ', lchnk, lat, lon, cam_out%flwds(i)
         write(iulog,*) 'df ml_verify_bottom, precsc : ', lchnk, lat, lon, cam_out%precsc(i)
         write(iulog,*) 'df ml_verify_bottom, precsl : ', lchnk, lat, lon, cam_out%precsl(i)
         write(iulog,*) 'df ml_verify_bottom, precc  : ', lchnk, lat, lon, cam_out%precc(i)
         write(iulog,*) 'df ml_verify_bottom, soll   : ', lchnk, lat, lon, cam_out%soll(i)
         write(iulog,*) 'df ml_verify_bottom, sols   : ', lchnk, lat, lon, cam_out%sols(i)
         write(iulog,*) 'df ml_verify_bottom, solld  : ', lchnk, lat, lon, cam_out%solld(i)
         write(iulog,*) 'df ml_verify_bottom, solsd  : ', lchnk, lat, lon, cam_out%solsd(i)
         write(iulog,*) 'df ml_verify_bottom, srfrad : ', lchnk, lat, lon, cam_out%srfrad(i)
         write(iulog,*) 'df ml_verify_bottom, thbot  : ', lchnk, lat, lon, cam_out%thbot(i)
         write(iulog,*) 'df ml_verify_bottom, co2prog: ', lchnk, lat, lon, cam_out%co2prog(i)
         write(iulog,*) 'df ml_verify_bottom, co2diag: ', lchnk, lat, lon, cam_out%co2diag(i)
         write(iulog,*) 'df ml_verify_bottom, psl    : ', lchnk, lat, lon, cam_out%psl(i)
         write(iulog,*) 'df ml_verify_bottom, bcphiwet: ', lchnk, lat, lon, cam_out%bcphiwet(i)
         write(iulog,*) 'df ml_verify_bottom, bcphidry: ', lchnk, lat, lon, cam_out%bcphidry(i)
         write(iulog,*) 'df ml_verify_bottom, bcphodry: ', lchnk, lat, lon, cam_out%bcphodry(i)
         write(iulog,*) 'df ml_verify_bottom, ocphiwet: ', lchnk, lat, lon, cam_out%ocphiwet(i)
         write(iulog,*) 'df ml_verify_bottom, ocphidry: ', lchnk, lat, lon, cam_out%ocphidry(i)
         write(iulog,*) 'df ml_verify_bottom, ocphodry: ', lchnk, lat, lon, cam_out%ocphodry(i)
         write(iulog,*) 'df ml_verify_bottom, dstwet1: ', lchnk, lat, lon, cam_out%dstwet1(i)
         write(iulog,*) 'df ml_verify_bottom, dstdry1: ', lchnk, lat, lon, cam_out%dstdry1(i)
         write(iulog,*) 'df ml_verify_bottom, dstwet2: ', lchnk, lat, lon, cam_out%dstwet2(i)
         write(iulog,*) 'df ml_verify_bottom, dstdry2: ', lchnk, lat, lon, cam_out%dstdry2(i)
         write(iulog,*) 'df ml_verify_bottom, dstwet3: ', lchnk, lat, lon, cam_out%dstwet3(i)
         write(iulog,*) 'df ml_verify_bottom, dstdry3: ', lchnk, lat, lon, cam_out%dstdry3(i)
         write(iulog,*) 'df ml_verify_bottom, dstwet4: ', lchnk, lat, lon, cam_out%dstwet4(i)
         write(iulog,*) 'df ml_verify_bottom, dstdry4: ', lchnk, lat, lon, cam_out%dstdry4(i)
     end do





  end subroutine

  function find_coord(state)
     use phys_grid,       only: get_lat_p, get_lon_p 
     use physics_types,   only: physics_state
     use cam_logfile,       only : iulog

     implicit none

     type(physics_state), intent(in) :: state

     integer lchnk                              ! chunk identifier
     integer ncol                               ! number of atmospheric columns
     integer i
     integer find_coord

     real(r8) :: lon
     real(r8) :: lat

     !tropics
     !real(r8), parameter :: the_lat = 6._r8
     !real(r8), parameter :: the_lon = 0._r8
     
     !nh greenland
     !real(r8), parameter :: the_lat = 78._r8
     !real(r8), parameter :: the_lon = 320._r8

     !the location of ice flux error
     real(r8), parameter :: the_lat = -67.2631578947369_r8
     real(r8), parameter :: the_lon = 85._r8


     real(r8), parameter :: pi = 3.14159265359_r8
     real(r8), parameter :: rad_to_deg = 180._r8/pi
     

     !real(r8), parameter :: the_lat = -6._r8
     !real(r8), parameter :: the_lon = 180._r8



     lchnk = state%lchnk
     ncol  = state%ncol

     if (found_location .and. lchnk == chunk) then
        find_coord = the_col
     end if



     do i=1,ncol
         lat = state%lat(i) * rad_to_deg 
         lon = state%lon(i) * rad_to_deg 

          !write(iulog,*) 'df testing : ', lchnk, lat, lon, the_lat, the_lon 

         !find the first coord within 10 degees of the desired lat/lon
         if ( (abs(lat-the_lat) < 1.) .and. (abs(lon-the_lon) < 1.) ) then
            found_location = .true.
            chunk = lchnk
            the_col = i 
            write(iulog,*) 'df testing found: ', lchnk, lat, lon, the_lat, the_lon, the_col 
            !found_location = .true.
            find_coord = the_col
            return
         end if
     end do

     find_coord = -1

  end function find_coord



  subroutine init_coord(state)
     use phys_grid,       only: get_lat_p, get_lon_p 
     use physics_types,   only: physics_state
     use cam_logfile,       only : iulog

     implicit none

     type(physics_state), intent(in) :: state

     integer lchnk                              ! chunk identifier
     integer ncol                               ! number of atmospheric columns
     integer i

     real(r8) :: lon
     real(r8) :: lat

     real(r8), parameter :: the_lat = 6._r8
     real(r8), parameter :: the_lon = 1._r8
     real(r8), parameter :: pi = 3.14159265359_r8
     real(r8), parameter :: rad_to_deg = 180._r8/pi
     

     !real(r8), parameter :: the_lat = -6._r8
     !real(r8), parameter :: the_lon = 180._r8



     if (found_location) then
        return
     end if

     lchnk = state%lchnk
     ncol  = state%ncol

     do i=1,ncol
         lat = state%lat(i) * rad_to_deg 
         lon = state%lon(i) * rad_to_deg 

          !write(iulog,*) 'df testing : ', lchnk, lat, lon, the_lat, the_lon 

         !find the first coord within 10 degees of the desired lat/lon
         if ( (abs(lat-the_lat) < 1.) .and. (abs(lon-the_lon) < 1.) ) then
            found_location = .true.
            chunk = lchnk
            the_col = i 
            write(iulog,*) 'df testing found: ', lchnk, lat, lon, the_lat, the_lon, the_col 

            exit 
         end if
     end do


  end subroutine init_coord


  subroutine old_init_coord(state)
     use phys_grid,       only: get_lat_p, get_lon_p 
     use physics_types,   only: physics_state
     use cam_logfile,       only : iulog

     implicit none

     type(physics_state), intent(in) :: state

     integer lchnk                              ! chunk identifier
     integer ncol                               ! number of atmospheric columns
     integer i

     real(r8) :: lon
     real(r8) :: lat

     !tropics
     !real(r8), parameter :: the_lat = 6._r8
     !real(r8), parameter :: the_lon = 1._r8

     real(r8), parameter :: the_lat = 78._r8
     real(r8), parameter :: the_lon = 320._r8

     !real(r8), parameter :: the_lat = -6._r8
     !real(r8), parameter :: the_lon = 180._r8



     if (found_location) then
        return
     end if

     lchnk = state%lchnk
     ncol  = state%ncol

     do i=1,ncol
         lat = get_lat_p(lchnk, i)
         lon = get_lon_p(lchnk, i)

          !write(iulog,*) 'df testing : ', lchnk, lat, lon, the_lat, the_lon 

         !find the first coord within 10 degees of the desired lat/lon
         if ( (abs(lat-the_lat) < 3.) .and. (abs(lon-the_lon) < 3.) ) then
            found_location = .true.
            chunk = lchnk
            the_col = i 
            write(iulog,*) 'df testing found: ', lchnk, lat, lon, the_lat, the_lon, the_col 

            exit 
         end if
     end do


  end subroutine old_init_coord


end module diag_ml

module physpkg
!-----------------------------------------------------------------------
! Purpose:
!
! Provides the interface to CAM physics package
!
! Revision history:
! Aug  2005,  E. B. Kluzek,  Creation of module from physpkg subroutine
! 2005-10-17  B. Eaton       Add contents of inti.F90 to phys_init().  Add
!                            initialization of grid info in phys_state.
! Nov 2010    A. Gettelman   Put micro/macro physics into separate routines
!-----------------------------------------------------------------------
  use shr_kind_mod,     only: r8 => shr_kind_r8
  use spmd_utils,       only: masterproc
  use physics_types,    only: physics_state, physics_tend, physics_state_set_grid, &
                              physics_ptend, physics_tend_init, physics_update,    &
                              physics_ptend_init, physics_type_alloc
  use phys_grid,        only: get_ncols_p
  use phys_gmean,       only: gmean_mass
  use ppgrid,           only: begchunk, endchunk, pcols
  use constituents,     only: pcnst, cnst_name, cnst_get_ind
  use camsrfexch_types, only: cam_out_t, cam_in_t
  use phys_buffer,      only: pbuf_old_tim_idx, pbuf_fld, pbuf_size_max, pbuf_get_fld_idx 
  use cam_control_mod,  only: ideal_phys, adiabatic
  use phys_control,     only: phys_do_flux_avg, waccmx_is
  use scamMod,          only: single_column, scm_crm_mode
  use flux_avg,         only: flux_avg_init
  use cldwat,           only: inimc
#ifdef SPMD
  use mpishorthand
#endif
   use perf_mod
   use cam_logfile,     only: iulog

   implicit none
   private

!  Physics buffer index
   integer :: teout_idx   = 0  
   save

! Public methods
   public phys_init   ! Public initialization method
   public phys_run1   ! First phase of the public run method
   public phys_run2   ! Second phase of the public run method
   public phys_final  ! Public finalization method
!
! Private module data
!

!======================================================================= 
contains
!======================================================================= 

subroutine phys_inidat( cam_out )
  use abortutils, only : endrun
  use buffer, only : tpert, qpert,pblht
  use phys_buffer,      only: pbuf, pbuf_times, pbuf_get_fld_idx
  use startup_initialconds, only : initial_file_get_id, topo_file_get_id
  use pio,            only : file_desc_t
  use ppgrid,         only: pver, pverp
  use ncdio_atm,          only: infld
  use dycore,            only: dycore_is
  use polar_avg,     only: polar_average
  use camsrfexch_types, only: cam_out_t
  use short_lived_species, only: initialize_short_lived_species
  use infnan, only : inf
  use comsrf,       only: landm, sgh, sgh30
  use cam_control_mod,    only: aqua_planet

  
  type(cam_out_t),         intent(inout) :: cam_out(begchunk:endchunk)
  integer :: lchnk, m, n, i, k, ncol
  type(file_desc_t), pointer :: ncid_ini, ncid_topo
  character(len=8) :: fieldname
  real(r8), pointer :: cldptr(:,:,:,:), convptr_3d(:,:,:,:)
  real(r8), pointer :: tptr(:,:), tptr3d(:,:,:)

  character*11 :: subname='phys_inidat' ! subroutine name

  logical :: found=.false., found2=.false.
  integer :: ierr
  character(len=4) :: dim1name
  integer :: ixcldice, ixcldliq
  nullify(tptr,tptr3d,cldptr,convptr_3d)

  ncid_ini=>initial_file_get_id()

!   dynamics variables are handled in dyn_init - here we read variables needed for physics 
!   but not dynamics

  if(dycore_is('UNSTRUCTURED')) then  
     dim1name='ncol'
  else
     dim1name='lon'
  end if
  if(aqua_planet) then
     sgh = 0._r8
     sgh30 = 0._r8
     landm = 0._r8
  else
     ncid_topo=>topo_file_get_id()
     call infld('SGH', ncid_topo, dim1name, 'lat', 1, pcols, begchunk, endchunk, &
          sgh, found, grid_map='PHYS')
     if(.not. found) call endrun('ERROR: SGH not found on topo file')
     
     call infld('SGH30', ncid_topo, dim1name, 'lat', 1, pcols, begchunk, endchunk, &
          sgh30, found, grid_map='PHYS')
     if(.not. found) then
        if (masterproc) write(iulog,*) 'Warning: Error reading SGH30 from topo file.'
        if (masterproc) write(iulog,*) 'The field SGH30 will be filled using data from SGH.'
        sgh30 = sgh
     end if
     
     call infld('LANDM_COSLAT', ncid_topo, dim1name, 'lat', 1, pcols, begchunk, endchunk, &
          landm, found, grid_map='PHYS')

     if(.not.found) call endrun(' ERROR: LANDM_COSLAT not found on topo dataset.')
  end if

  call infld('PBLH', ncid_ini, dim1name, 'lat', 1, pcols, begchunk, endchunk, &
       pblht, found, grid_map='PHYS')
  if(.not. found) then
     pblht(:,:) = 0._r8
     if (masterproc) write(iulog,*) 'PBLH initialized to 0.'
  end if

  call infld('TPERT', ncid_ini, dim1name, 'lat', 1, pcols, begchunk, endchunk, &
       tpert, found, grid_map='PHYS')
  if(.not. found) then
     tpert(:,:) = 0._r8
     if (masterproc) write(iulog,*) 'TPERT initialized to 0.'
  end if

  fieldname='QPERT'  
  qpert(:,:,:) = 0._r8
  allocate(tptr(1:pcols,begchunk:endchunk))
  call infld(fieldname, ncid_ini, dim1name, 'lat', 1, pcols, begchunk, endchunk, &
       tptr, found, grid_map='PHYS')
  if(found) then
     do lchnk=begchunk,endchunk
        qpert(:,1,lchnk) = tptr(:,lchnk)
     end do
  else
     if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
  end if
  deallocate(tptr)
   
  fieldname='CUSH'
  m = pbuf_get_fld_idx('cush')
  tptr3d => pbuf(m)%fld_ptr(1,1:pcols,1,begchunk:endchunk,1:pbuf_times)
  call infld(fieldname, ncid_ini, dim1name, 'lat', 1, pcols, begchunk, endchunk, &
       tptr3d(:,:,1), found, grid_map='PHYS')
  if(.not.found) then
     if(masterproc) write(iulog,*) trim(fieldname), ' initialized to 1000.'
     tptr3d=1000._r8
  else
     do n=2,pbuf_times
        tptr3d(:,:,n)=tptr3d(:,:,1)
     end do
  end if


  do lchnk=begchunk,endchunk
     cam_out(lchnk)%tbot(:) = inf
  end do


!
! 3-D fields
!
  fieldname='CLOUD'
  m = pbuf_get_fld_idx('CLD')
  if(associated(pbuf(m)%fld_ptr)) then
     cldptr => pbuf(m)%fld_ptr(1,:,:,begchunk:endchunk,:)
  else
     call endrun('pbuf not allocated in phys_inidat') 
  end if

  call infld(fieldname,ncid_ini,dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
       cldptr(:,:,:,1), found, grid_map='PHYS')

  if(.not. found) then
     if(masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
     cldptr(:,:,:,:) = 0._r8     
  else
     do n = 2, pbuf_times
        cldptr(:,:,:,n) = cldptr(:,:,:,1)
     end do
  end if

  fieldname = 'QCWAT'
  m = pbuf_get_fld_idx(fieldname)
  if(associated(pbuf(m)%fld_ptr)) then
     cldptr => pbuf(m)%fld_ptr(1,:,:,begchunk:endchunk,:)
  else
     call endrun('pbuf not allocated in phys_inidat') 
  end if
  call infld(fieldname,ncid_ini,dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
       cldptr(:,:,:,1), found, grid_map='PHYS')
  if(.not.found) then
     call infld('Q',ncid_ini,dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
          cldptr(:,:,:,1), found, grid_map='PHYS')
     if (found) then
        if (masterproc) write(iulog,*) trim(fieldname), ' initialized with Q'
	if(dycore_is('LR')) call polar_average(pver, cldptr(:,:,:,1)) 	
     else
        call endrun('  '//trim(subname)//' Error:  Q must be on Initial File')
     end if
  end if
  do n = 2, pbuf_times
     cldptr(:,:,:,n) = cldptr(:,:,:,1)
  end do 

  fieldname = 'ICCWAT'
  m = pbuf_get_fld_idx(fieldname)

  if(associated(pbuf(m)%fld_ptr)) then
     cldptr => pbuf(m)%fld_ptr(1,:,:,begchunk:endchunk,:)
  else
     call endrun('pbuf not allocated in phys_inidat') 
  end if

  cldptr=0._r8
  call infld(fieldname,ncid_ini,dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
       cldptr(:,:,:,1), found, grid_map='PHYS')
  if(.not. found) then
     allocate(tptr3d(pcols,pver,begchunk:endchunk))     
     call cnst_get_ind('CLDICE', ixcldice)
     call infld('CLDICE',ncid_ini,dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
          tptr3d, found, grid_map='PHYS')
     if(found) then
        cldptr(:,:,:,1)=tptr3d(:,:,:)
     end if
     if (masterproc) then
        if (found) then
           write(iulog,*) trim(fieldname), ' initialized with CLDICE'
        else
           write(iulog,*) trim(fieldname), ' initialized to 0.0'
        end if
     end if
     deallocate(tptr3d)
  end if
  do n = 2, pbuf_times
     cldptr(:,:,:,n) = cldptr(:,:,:,1)
  end do

  fieldname = 'LCWAT'
  m = pbuf_get_fld_idx(fieldname)

  if(associated(pbuf(m)%fld_ptr)) then
     cldptr => pbuf(m)%fld_ptr(1,:,:,begchunk:endchunk,:)
  else
     call endrun('pbuf not allocated in phys_inidat') 
  end if

  call infld(fieldname,ncid_ini,dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
       cldptr(:,:,:,1), found, grid_map='PHYS')
  if(.not. found) then
     allocate(tptr3d(pcols,pver,begchunk:endchunk))     
     call cnst_get_ind('CLDICE', ixcldice)
     call cnst_get_ind('CLDLIQ', ixcldliq)
     call infld('CLDICE',ncid_ini,dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
          tptr3d, found, grid_map='PHYS')
     if(found) then
        cldptr(:,:,:,1)=tptr3d(:,:,:)
     end if
     call infld('CLDLIQ',ncid_ini,dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
          tptr3d, found2, grid_map='PHYS')
     if(found2) then
        cldptr(:,:,:,1)=cldptr(:,:,:,1)+tptr3d(:,:,:)
     end if

     if (found .or. found2) then
        if (masterproc) write(iulog,*) trim(fieldname), ' initialized with CLDICE + CLDLIQ'
        if(dycore_is('LR')) call polar_average(pver, cldptr(:,:,:,1)) 	
     else
        if (masterproc)  write(iulog,*) trim(fieldname), ' initialized to 0.0'
        cldptr(:,:,:,1) = 0._r8
     end if
     deallocate(tptr3d)
  end if
  do n = 2, pbuf_times
     cldptr(:,:,:,n) = cldptr(:,:,:,1)
  end do

  fieldname = 'TCWAT'
  m = pbuf_get_fld_idx(fieldname)
  if(associated(pbuf(m)%fld_ptr)) then
     cldptr => pbuf(m)%fld_ptr(1,:,:,begchunk:endchunk,:)
  else
     call endrun('pbuf not allocated in phys_inidat') 
  end if
  call infld(fieldname,ncid_ini,dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
       cldptr(:,:,:,1), found, grid_map='PHYS')
  if(.not. found) then
     call infld('T',ncid_ini,dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
          cldptr(:,:,:,1), found, grid_map='PHYS')
     if(dycore_is('LR')) call polar_average(pver, cldptr(:,:,:,1)) 	
     if (masterproc) write(iulog,*) trim(fieldname), ' initialized with T'
  end if
  do n = 2, pbuf_times
     cldptr(:,:,:,n) = cldptr(:,:,:,1)
  end do

  fieldname = 'TKE'
  m = pbuf_get_fld_idx('tke')
  convptr_3d => pbuf(m)%fld_ptr(1,1:pcols,1:pverp,begchunk:endchunk,1:pbuf_times)
  call infld(fieldname, ncid_ini, dim1name, 'ilev', 'lat', 1, pcols, 1, pverp, begchunk, endchunk, &
       convptr_3d(:,:,:,1), found, grid_map='phys')
  if(.not. found) then
     convptr_3d(:,:,:,1) = 0.01_r8
     if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.01'
  end if
  do n = 2, pbuf_times
     convptr_3d(:,:,:,n) = convptr_3d(:,:,:,1)
  end do

  fieldname = 'KVM'
  m = pbuf_get_fld_idx('kvm')
  convptr_3d => pbuf(m)%fld_ptr(1,1:pcols,1:pverp,begchunk:endchunk,1:pbuf_times)
  call infld(fieldname, ncid_ini, dim1name, 'ilev', 'lat', 1, pcols, 1, pverp, begchunk, endchunk, &
       convptr_3d(:,:,:,1), found, grid_map='phys')
  if(.not. found) then
     convptr_3d(:,:,:,1) = 0._r8
     if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
  end if
  do n = 2, pbuf_times
     convptr_3d(:,:,:,n) = convptr_3d(:,:,:,1)
  end do

  fieldname = 'KVH'
  m = pbuf_get_fld_idx('kvh')
  convptr_3d => pbuf(m)%fld_ptr(1,1:pcols,1:pverp,begchunk:endchunk,1:pbuf_times)
  call infld(fieldname, ncid_ini, dim1name, 'ilev', 'lat', 1, pcols, 1, pverp, begchunk, endchunk, &
       convptr_3d(:,:,:,1), found, grid_map='phys')
  if(.not. found) then
     convptr_3d(:,:,:,1) = 0._r8
     if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
  end if
  do n = 2, pbuf_times
     convptr_3d(:,:,:,n) = convptr_3d(:,:,:,1)
  end do

  fieldname = 'CONCLD'
  m = pbuf_get_fld_idx('CONCLD')
  convptr_3d => pbuf(m)%fld_ptr(1,1:pcols,1:pver,begchunk:endchunk,1:pbuf_times)
  call infld(fieldname, ncid_ini, dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
       convptr_3d(:,:,:,1), found, grid_map='phys')
  if(.not. found) then
     convptr_3d(:,:,:,1) = 0.
     if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
  end if
  do n = 2, pbuf_times
     convptr_3d(:,:,:,n) = convptr_3d(:,:,:,1)
  end do

  call initialize_short_lived_species(ncid_ini)
end subroutine phys_inidat


subroutine phys_init( phys_state, phys_tend, pbuf, cam_out )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Initialization of physics package.
! 
!-----------------------------------------------------------------------

   use physconst,          only: rair, cpair, cpwv, gravit, stebol, epsilo, tmelt, &
                                 latvap, latice, rh2o, rhoh2o, pstd,  &
                                 karman, rhodair, physconst_init 
   use phys_control,        only: phys_ctl_init
   use ref_pres,           only: pref_edge, pref_mid

   use aerosol_intr,       only: prognostic_aerosol_initialize
   use cloud_rad_props,    only: cloud_rad_props_init
   use cam_control_mod,    only: nsrest  ! restart flag
   use check_energy,       only: check_energy_init
   use chemistry,          only: chem_init
   use prescribed_ozone,   only: prescribed_ozone_init
   use prescribed_ghg,     only: prescribed_ghg_init
   use prescribed_aero,    only: prescribed_aero_init
   use aerodep_flx,        only: aerodep_flx_init
   use aircraft_emit,      only: aircraft_emit_init
   use prescribed_volcaero,only: prescribed_volcaero_init
   use cloud_fraction,     only: cldfrc_init
   use co2_cycle,          only: co2_init, co2_transport
   use convect_deep,       only: convect_deep_init
   use convect_shallow,    only: convect_shallow_init
   use machine_learning_model,    only: ml_model_init
   use export_state,              only: export_state_init
   use cam_diagnostics,    only: diag_init
   use gw_drag,            only: gw_inti
   use cam3_aero_data,     only: cam3_aero_data_on, cam3_aero_data_init
   use cam3_ozone_data,    only: cam3_ozone_data_on, cam3_ozone_data_init
   use radheat,            only: radheat_init
   use radiation,          only: radiation_init
   use cloud_diagnostics,  only: cloud_diagnostics_init
   use stratiform,         only: stratiform_init
   use phys_control,       only: phys_getopts, waccmx_is
   use microp_driver,      only: microp_driver_init
   use macrop_driver,      only: macrop_driver_init
   use conv_water,         only: conv_water_init
   use tracers,            only: tracers_init
   use aoa_tracers,        only: aoa_tracers_init
   use rayleigh_friction,  only: rayleigh_friction_init
   use vertical_diffusion, only: vertical_diffusion_init
   use dycore,             only: dycore_is
   use phys_debug_util,    only: phys_debug_init
   use phys_buffer,        only: pbuf_size_max, pbuf_fld
   use rad_constituents,   only: rad_cnst_init
   use aer_rad_props,      only: aer_rad_props_init
#if ( defined WACCM_PHYS )
   use qbo,                only: qbo_init
   use iondrag,            only: iondrag_init
#endif
#if ( defined OFFLINE_DYN )
   use metdata,            only: metdata_phys_init
#endif
   use ionosphere,	   only: ionos_init  ! Initialization of ionosphere module (WACCM-X)
   use majorsp_diffusion,  only: mspd_init   ! Initialization of major species diffusion module (WACCM-X)
   use sslt_rebin,         only: sslt_rebin_init
   use tropopause,         only: tropopause_init
   use solar_data,         only: solar_data_init
   use rad_solar_var,      only: rad_solar_var_init

   ! Input/output arguments
   type(physics_state), pointer :: phys_state(:)
   type(physics_tend ), pointer :: phys_tend(:)
   type(pbuf_fld), intent(in), dimension(pbuf_size_max) :: pbuf  ! physics buffer
   type(cam_out_t),intent(inout)                        :: cam_out(begchunk:endchunk)

   ! local variables
   integer :: lchnk

! Get microphysics option

   character(len=16) :: microp_scheme 
   call phys_getopts( microp_scheme_out = microp_scheme )

   !-----------------------------------------------------------------------

   call physics_type_alloc(phys_state, phys_tend, begchunk, endchunk)

   do lchnk = begchunk,endchunk
      call physics_state_set_grid(lchnk, phys_state(lchnk))
   end do
   
   !-------------------------------------------------------------------------------------------
   ! Initialize any variables in physconst which are not temporally and/or spatially constant
   !------------------------------------------------------------------------------------------- 
   call physconst_init()

   !--------------------------------------
   ! Initialize phys_control variables 
   !--------------------------------------
   call phys_ctl_init

   ! Initialize debugging a physics column
   call phys_debug_init()

   ! diag_init makes addfld calls for dynamics fields that are output from
   ! the physics decomposition
   call diag_init()

   call check_energy_init()

   call tracers_init()

   ! age of air tracers
   call aoa_tracers_init()

   ! For adiabatic or ideal physics don't need to initialize any of the
   ! parameterizations below:
   if (adiabatic .or. ideal_phys) return

   if (nsrest .eq. 0) then
      call phys_inidat(cam_out) 
   end if

   ! CAM3 prescribed aerosols
   if (cam3_aero_data_on) call cam3_aero_data_init(phys_state(begchunk:endchunk))

   ! Initialize rad constituents and their properties
   call rad_cnst_init(pbuf, phys_state)
   call aer_rad_props_init()
   call cloud_rad_props_init()

   ! Initialize prognostic aerosols
   call prognostic_aerosol_initialize(phys_state(begchunk:endchunk))

   ! solar irradiance data modules
   call solar_data_init()

   ! Prognostic chemistry.
   call chem_init(phys_state(begchunk:endchunk))
 
   ! Prescribed tracers
   call prescribed_ozone_init()
   call prescribed_ghg_init()
   call prescribed_aero_init()
   call aerodep_flx_init()
   call aircraft_emit_init()
   call prescribed_volcaero_init()

   ! co2 cycle            
   if (co2_transport()) then
      call co2_init()
   end if

   ! CAM3 prescribed ozone
   if (cam3_ozone_data_on) call cam3_ozone_data_init(phys_state(begchunk:endchunk))

   call gw_inti(cpair, cpwv, gravit, rair, pref_edge)

   call rayleigh_friction_init()

   call vertical_diffusion_init

   if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 
     call mspd_init ()
     ! Initialization of ionosphere module if mode set to ionosphere
     if( waccmx_is('ionosphere') ) then
       call ionos_init()
     endif
   endif

   call tsinti(tmelt, latvap, rair, stebol, latice)

   call radiation_init
   call rad_solar_var_init()

   call cloud_diagnostics_init

   call radheat_init(pref_mid)

   call esinti(epsilo, latvap, latice, rh2o, cpair, tmelt)

   call convect_shallow_init(pref_edge)
   call ml_model_init(pref_edge)
   call export_state_init()

   call cldfrc_init

   call convect_deep_init(pref_edge)

   if( microp_scheme .eq. 'RK' ) then
      call stratiform_init
   elseif( microp_scheme .eq. 'MG' ) then 
      call macrop_driver_init
      call microp_driver_init
      call conv_water_init
   end if
   call inimc(tmelt, rhodair/1000.0_r8, gravit, rh2o)

#if ( defined WACCM_PHYS )
   call iondrag_init(pref_mid)
   call qbo_init
#endif

#if ( defined OFFLINE_DYN )
   call metdata_phys_init()
#endif
   call sslt_rebin_init
   call tropopause_init()

end subroutine phys_init

!
!-----------------------------------------------------------------------
!

subroutine phys_run1(phys_state, ztodt, phys_tend, pbuf, cam_in, cam_out)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! First part of atmospheric physics package before updating of surface models
! 
!-----------------------------------------------------------------------
   use ppgrid,         only: pver
   use time_manager,   only: get_nstep
   use cam_diagnostics,only: diag_allocate, diag_physvar_ic
   use check_energy,   only: check_energy_gmean
   use phys_buffer,    only: pbuf_allocate
   use buffer,         only: pblht, tpert, qpert
#if (defined BFB_CAM_SCAM_IOP )
   use cam_history,    only: outfld
#endif
   use comsrf,         only: fsns, fsnt, flns, flnt, landm, fsds
   use abortutils, only :endrun
!
! Input arguments
!
   real(r8), intent(in) :: ztodt            ! physics time step unless nstep=0
!
! Input/Output arguments
!
   type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
   type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend
   type(pbuf_fld),      intent(inout), dimension(pbuf_size_max)     :: pbuf
   type(cam_in_t),                     dimension(begchunk:endchunk) :: cam_in
   type(cam_out_t),                    dimension(begchunk:endchunk) :: cam_out
!-----------------------------------------------------------------------
!
!---------------------------Local workspace-----------------------------
!
   integer :: c                                 ! indices
   integer :: ncol                              ! number of columns
   integer :: nstep                             ! current timestep number
#if (! defined SPMD)
   integer  :: mpicom = 0
#endif

   call t_startf ('physpkg_st1')
   nstep = get_nstep()

   ! The following initialization depends on the import state (cam_in)
   ! being initialized.  This isn't true when cam_init is called, so need
   ! to postpone this initialization to here.
   if (nstep == 0 .and. phys_do_flux_avg()) call flux_avg_init(cam_in, pbuf)

   ! Compute total energy of input state and previous output state
   call t_startf ('chk_en_gmean')
   call check_energy_gmean(phys_state, pbuf, ztodt, nstep)
   call t_stopf ('chk_en_gmean')

   call t_stopf ('physpkg_st1')

   if ( adiabatic .or. ideal_phys )then
      call t_startf ('bc_physics')
      call phys_run1_adiabatic_or_ideal(ztodt, phys_state, phys_tend, pbuf)
      call t_stopf ('bc_physics')
   else
      call t_startf ('physpkg_st1')

      call pbuf_allocate('physpkg')
      call diag_allocate()

      !-----------------------------------------------------------------------
      ! Advance time information
      !-----------------------------------------------------------------------

      call advnce( phys_state, cam_out, pbuf)

      call t_stopf ('physpkg_st1')

#ifdef TRACER_CHECK
      call gmean_mass ('before tphysbc DRY', phys_state)
#endif


      !-----------------------------------------------------------------------
      ! Tendency physics before flux coupler invocation
      !-----------------------------------------------------------------------
      !

#if (defined BFB_CAM_SCAM_IOP )
      do c=begchunk, endchunk
         call outfld('Tg',cam_in(c)%ts,pcols   ,c     )
      end do
#endif

      call t_barrierf('sync_bc_physics', mpicom)
      call t_startf ('bc_physics')
      call t_adj_detailf(+1)

!$OMP PARALLEL DO PRIVATE (C)
      do c=begchunk, endchunk
         !
         ! Output physics terms to IC file
         !
         call t_startf ('diag_physvar_ic')
         call diag_physvar_ic ( c, pbuf, cam_out(c), cam_in(c) )
         call t_stopf ('diag_physvar_ic')

         call tphysbc (ztodt, pblht(1,c), tpert(1,c), qpert(1,1,c),                      &
                       fsns(1,c), fsnt(1,c), flns(1,c), flnt(1,c), phys_state(c),        &
                       phys_tend(c), pbuf,  fsds(1,c), landm(1,c),                       &
                       cam_out(c), cam_in(c) )

      end do

      call t_adj_detailf(-1)
      call t_stopf ('bc_physics')

! Don't call the rest in CRM mode
      if(single_column.and.scm_crm_mode) return

#ifdef TRACER_CHECK
      call gmean_mass ('between DRY', phys_state)
#endif
   end if

end subroutine phys_run1

!
!-----------------------------------------------------------------------
!

subroutine phys_run1_adiabatic_or_ideal(ztodt, phys_state, phys_tend, pbuf)
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Physics for adiabatic or idealized physics case.
! 
!-----------------------------------------------------------------------
   use time_manager,     only: get_nstep
   use cam_diagnostics,  only: diag_phys_writeout
   use check_energy,     only: check_energy_fix, check_energy_chng
   use dycore,           only: dycore_is
!
! Input arguments
!
   real(r8), intent(in) :: ztodt            ! physics time step unless nstep=0
!
! Input/Output arguments
!
   type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
   type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend
   type(pbuf_fld),      intent(inout), dimension(pbuf_size_max)     :: pbuf
!-----------------------------------------------------------------------
!---------------------------Local workspace-----------------------------
!
   integer             :: c               ! indices
   integer             :: nstep           ! current timestep number
   type(physics_ptend) :: ptend           ! indivdual parameterization tendencies
   real(r8)            :: flx_heat(pcols) ! effective sensible heat flux
   real(r8)            :: zero(pcols)     ! array of zeros

   ! physics buffer field for total energy
   integer itim, ncol
   real(r8), pointer, dimension(:) :: teout
   logical, SAVE :: first_exec_of_phys_run1_adiabatic_or_ideal  = .TRUE.
!-----------------------------------------------------------------------

   nstep = get_nstep()
   zero  = 0._r8

   ! Associate pointers with physics buffer fields
   itim = pbuf_old_tim_idx()
   if (first_exec_of_phys_run1_adiabatic_or_ideal) then
      teout_idx = pbuf_get_fld_idx('TEOUT')
      first_exec_of_phys_run1_adiabatic_or_ideal  = .FALSE.
   endif
!$OMP PARALLEL DO PRIVATE (C, NCOL, TEOUT, PTEND, FLX_HEAT)
   do c=begchunk, endchunk

      ncol = get_ncols_p(c)

      ! Associate pointers with physics buffer fields
      teout  => pbuf(teout_idx)%fld_ptr(1,1:pcols,1,c,itim)

      ! Initialize the physics tendencies to zero.
      call physics_tend_init(phys_tend(c))

      ! Dump dynamics variables to history buffers
      call diag_phys_writeout(phys_state(c))

      if (dycore_is('LR')) then
         call physics_ptend_init(ptend)
         call check_energy_fix(phys_state(c), ptend, nstep, flx_heat)
         call physics_update(phys_state(c), phys_tend(c), ptend, ztodt)
         call check_energy_chng(phys_state(c), phys_tend(c), "chkengyfix", nstep, ztodt, &
                                zero, zero, zero, flx_heat)
      end if

      if ( ideal_phys )then
         call t_startf('tphysidl')
         call tphysidl(ztodt, phys_state(c), phys_tend(c))
         call t_stopf('tphysidl')
      end if

      ! Save total enery after physics for energy conservation checks
      teout(:ncol) = phys_state(c)%te_cur(:ncol)

   end do

end subroutine phys_run1_adiabatic_or_ideal

!
!-----------------------------------------------------------------------
!

subroutine phys_run2(phys_state, ztodt, phys_tend, pbuf, cam_out, &
                     cam_in )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Second part of atmospheric physics package after updating of surface models
! 
!-----------------------------------------------------------------------
   use buffer,         only: pblht, tpert, qpert
   use mo_lightning,   only: lightning_no_prod
   use phys_buffer,    only: pbuf_update_tim_idx
   use phys_buffer,    only: pbuf_deallocate
   use cam_diagnostics,only: diag_deallocate, diag_surf
   use comsrf,         only: trefmxav, trefmnav, sgh, sgh30, fsds 
   use physconst,      only: stebol, latvap
!
! Input arguments
!
   real(r8), intent(in) :: ztodt                       ! physics time step unless nstep=0
!
! Input/Output arguments
!
   type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
   type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend
   type(pbuf_fld),      intent(inout), dimension(pbuf_size_max)     :: pbuf
   type(cam_out_t),     intent(inout), dimension(begchunk:endchunk) :: cam_out
   type(cam_in_t),      intent(inout), dimension(begchunk:endchunk) :: cam_in
!
!-----------------------------------------------------------------------
!---------------------------Local workspace-----------------------------
!
   integer :: c                                 ! chunk index
   integer :: ncol                              ! number of columns
#if (! defined SPMD)
   integer  :: mpicom = 0
#endif
   !
   ! If exit condition just return
   !

   if(single_column.and.scm_crm_mode) return

   if ( adiabatic .or. ideal_phys ) return
   !-----------------------------------------------------------------------
   ! Tendency physics after coupler 
   ! Not necessary at terminal timestep.
   !-----------------------------------------------------------------------
   !
   ! Set lightning production of NO
   call t_startf ('lightning_no_prod')
   call lightning_no_prod( phys_state, pbuf, cam_in )
   call t_stopf ('lightning_no_prod')

   call t_barrierf('sync_ac_physics', mpicom)
   call t_startf ('ac_physics')
   call t_adj_detailf(+1)

!$OMP PARALLEL DO PRIVATE (C, NCOL)

   do c=begchunk,endchunk
      ncol = get_ncols_p(c)
      !
      ! surface diagnostics for history files
      !
      call t_startf('diag_surf')
      call diag_surf(cam_in(c), cam_out(c), phys_state(c)%ps,trefmxav(1,c), trefmnav(1,c))
      call t_stopf('diag_surf')

      call tphysac(ztodt, pblht(1,c), qpert(1,1,c), tpert(1,c), cam_in(c),        &
                    sgh(1,c), sgh30(1,c), cam_out(c),                              &
                    phys_state(c), phys_tend(c), pbuf, fsds(1,c))
   end do                    ! Chunk loop

   call t_adj_detailf(-1)
   call t_stopf('ac_physics')

#ifdef TRACER_CHECK
   call gmean_mass ('after tphysac FV:WET)', phys_state)
#endif

   call t_startf ('physpkg_st2')
   call pbuf_deallocate('physpkg')
   call pbuf_update_tim_idx()
   call diag_deallocate()
   call t_stopf ('physpkg_st2')

end subroutine phys_run2

!
!----------------------------------------------------------------------- 
!

subroutine phys_final( phys_state, phys_tend )
  use chemistry, only : chem_final
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Finalization of physics package
! 
!-----------------------------------------------------------------------
   ! Input/output arguments
   type(physics_state), pointer :: phys_state(:)
   type(physics_tend ), pointer :: phys_tend(:)

   deallocate(phys_state)
   deallocate(phys_tend)
   call chem_final

end subroutine phys_final

end module physpkg

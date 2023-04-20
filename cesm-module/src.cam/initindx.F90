subroutine initindx
!----------------------------------------------------------------------- 
! 
! Purpose: Register constituents and physics buffer fields.
! 
! Author:    CSM Contact: M. Vertenstein, Aug. 1997
!            B.A. Boville, Oct 2001
!            A. Gettelman, Nov 2010 - put micro/macro physics into separate routines
! 
!-----------------------------------------------------------------------
  use shr_kind_mod,       only: r8 => shr_kind_r8
  use spmd_utils,         only: masterproc
  use constituents,       only: pcnst, cnst_add, cnst_chk_dim, cnst_name
  use phys_buffer,        only: pbuf_init, pbuf_print
  use cam_control_mod,    only: moist_physics
  use phys_control,       only: phys_do_flux_avg, phys_getopts, waccmx_is
  use chemistry,          only: chem_register
  use stratiform,         only: stratiform_register
  use microp_driver,      only: microp_driver_register
  use macrop_driver,      only: macrop_driver_register
  use conv_water,         only: conv_water_register
  use physconst,          only: mwdry, cpair, mwh2o, cpwv
  use tracers,            only: tracers_register
  use check_energy,       only: check_energy_register
  use aerosol_intr,       only: aerosol_register_cnst
  use cam3_aero_data,     only: cam3_aero_data_on, cam3_aero_data_register
  use cam3_ozone_data,    only: cam3_ozone_data_on, cam3_ozone_data_register
  use ghg_data,           only: ghg_data_register
  use vertical_diffusion, only: vd_register
  use convect_deep,       only: convect_deep_register
  use convect_shallow,    only: convect_shallow_register
  use machine_learning_model,    only: ml_model_register
  use radiation,          only: radiation_register
  use co2_cycle,          only: co2_register
  use flux_avg,           only: flux_avg_register
#if ( defined WACCM_PHYS )
  use exbdrift,           only: exbdrift_register
  use gw_drag,            only: gw_drag_register
#endif
  use iondrag,            only: iondrag_register
  use ionosphere,         only: ionos_register
  use string_utils,       only: to_lower
  use prescribed_ozone,   only: prescribed_ozone_register
  use prescribed_volcaero,only: prescribed_volcaero_register
  use prescribed_aero,    only: prescribed_aero_register
  use prescribed_ghg,     only: prescribed_ghg_register
  use sslt_rebin,         only: sslt_rebin_register
  use aoa_tracers,        only: aoa_tracers_register
  use aircraft_emit,      only: aircraft_emit_register
  use cam_diagnostics,    only: diag_register

  implicit none
!---------------------------Local variables-----------------------------
!
  integer  :: m        ! loop index
  integer  :: mm       ! constituent index 
!-----------------------------------------------------------------------

   character(len=16) :: microp_scheme 
   call phys_getopts( microp_scheme_out = microp_scheme )

  ! Initialize physics buffer
  call pbuf_init()

  ! Register water vapor.
  ! ***** N.B. ***** This must be the first call to cnst_add so that
  !                  water vapor is constituent 1.
  if (moist_physics) then
     call cnst_add('Q', mwh2o, cpwv, 1.E-12_r8, mm, &
                   longname='Specific humidity', readiv=.true.)
  else
     call cnst_add('Q', mwh2o, cpwv, 0.0_r8, mm, &
                   longname='Specific humidity', readiv=.false.)
  end if

  ! check energy package
  call check_energy_register

  ! If using an ideal/adiabatic physics option, the CAM physics parameterizations 
  ! aren't called.
   if (moist_physics) then

      ! register fluxes for saving across time
      if (phys_do_flux_avg()) call flux_avg_register()

      ! cloud water
      if( microp_scheme .eq. 'RK' ) then
         call stratiform_register()
      elseif( microp_scheme .eq. 'MG' ) then
         call macrop_driver_register()
         call microp_driver_register()
      end if
      call conv_water_register()

      ! chemical constituents
      call chem_register()

      ! co2 constituents
      call co2_register()

      ! register data model ozone with pbuf
      if (cam3_ozone_data_on) then
         call cam3_ozone_data_register()
      end if
      call prescribed_volcaero_register()
      call prescribed_ozone_register()
      call prescribed_aero_register()
      call prescribed_ghg_register()
      call sslt_rebin_register

      ! CAM3 prescribed aerosols
      if (cam3_aero_data_on) then
         call cam3_aero_data_register()
      end if

      ! register various data model gasses with pbuf
      call ghg_data_register()

#if ( defined WACCM_PHYS )
      ! Initialize e and b fields
      call exbdrift_register()
      ! waccm gravity wave drag
      call gw_drag_register()
#endif

      if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then 
        ! Register iondrag variables with pbuf
        call iondrag_register()
        ! Register ionosphere variables with pbuf if mode set to ionosphere
        if( waccmx_is('ionosphere') ) then
          call ionos_register()
        endif
      endif

      ! aerosols
      call aerosol_register_cnst()

      call aircraft_emit_register()

      ! deep convection
      call convect_deep_register

      !  shallow convection
      call convect_shallow_register
 
      call ml_model_register

      ! radiation
      call radiation_register

      ! vertical diffusion
      call vd_register()

   end if

  ! Register diagnostics PBUF
   call diag_register()

  ! Register age of air tracers
  call aoa_tracers_register()

  ! Register test tracers
  ! ***** N.B. ***** This is the last call to register constituents because
  !                  the test tracers fill the remaining available slots up
  !                  to constituent number PCNST -- regardless of what PCNST is set to.
  call tracers_register()

  ! All tracers registered, check that the dimensions are correct
  call cnst_chk_dim()

  ! ***NOTE*** No registering constituents after the call to cnst_chk_dim.

  ! Registration of physics buffer fields is also complete.
  if (masterproc) call pbuf_print()

end subroutine initindx


subroutine tphysac (ztodt,   pblh,    qpert,   tpert,  cam_in,  &
                    sgh,     sgh30,                                     &
                    cam_out,  state,   tend,    pbuf,           &
                    fsds    )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Tendency physics after coupling to land, sea, and ice models.
! Computes the following:
!   o Radon surface flux and decay (optional)
!   o Vertical diffusion and planetary boundary layer
!   o Multiple gravity wave drag
! 
! Method: 
! <Describe the algorithm(s) used in the routine.> 
! <Also include any applicable external references.> 
! 
! Author: CCM1, CMS Contact: J. Truesdale
! 
!-----------------------------------------------------------------------
   use shr_kind_mod,       only: r8 => shr_kind_r8
   use ppgrid,             only: pcols, pver, pverp
   use chemistry,          only: chem_is_active, chem_timestep_tend
   use cam_diagnostics,    only: diag_phys_tend_writeout
   use gw_drag,            only: gw_intr
   use vertical_diffusion, only: vertical_diffusion_tend
   use rayleigh_friction,  only: rayleigh_friction_tend
   use constituents,       only: cnst_get_ind
   use physics_types,      only: physics_state, physics_tend, physics_ptend, physics_update,    &
                                 physics_ptend_init, physics_dme_adjust, set_dry_to_wet
   use cam_logfile,   only : iulog
   use majorsp_diffusion,  only: mspd_intr  ! WACCM-X major diffusion
   use ionosphere,         only: ionos_intr ! WACCM-X ionosphere
   use phys_buffer,        only: pbuf_size_max, pbuf_fld, pbuf_old_tim_idx, pbuf_get_fld_idx
   use tracers,            only: tracers_timestep_tend
   use aoa_tracers,        only: aoa_tracers_timestep_tend
   use physconst,          only: zvir, gravit, rhoh2o, latvap,latice, cpair, rair
   use aerosol_intr,       only: aerosol_emis_intr, aerosol_drydep_intr
   use camsrfexch_types,   only: cam_out_t, cam_in_t     
   use check_energy,       only: check_energy_chng
   use check_energy,       only: check_tracers_data, check_tracers_init, check_tracers_chng
   use time_manager,       only: get_nstep
   use abortutils,         only: endrun
   use dycore,             only: dycore_is
   use cam_control_mod,    only: aqua_planet 
   use mo_gas_phase_chemdr,only: map2chm
   use clybry_fam,         only: clybry_fam_set
#if ( defined WACCM_PHYS )
   use charge_neutrality,  only: charge_fix
   use iondrag,            only: iondrag_calc, do_waccm_ions
   use qbo,                only: qbo_relax
#endif
   use perf_mod
   use phys_control,       only: phys_do_flux_avg, waccmx_is
   use flux_avg,           only: flux_avg_run

   use machine_learning_model_config


   implicit none

!
! Arguments
!
   real(r8), intent(in) :: ztodt                  ! Two times model timestep (2 delta-t)
   real(r8), intent(inout) :: pblh(pcols)         ! Planetary boundary layer height
   real(r8), intent(in) :: fsds(pcols)            ! down solar flux
   real(r8), intent(out) :: qpert(pcols)          ! Moisture/constit. perturbation (PBL)
   real(r8), intent(out) :: tpert(pcols)          ! Temperature perturbation (PBL)
   real(r8), intent(in) :: sgh(pcols)             ! Std. deviation of orography for gwd
   real(r8), intent(in) :: sgh30(pcols)           ! Std. deviation of 30s orography for tms

   type(cam_in_t),      intent(inout) :: cam_in
   type(cam_out_t),     intent(inout) :: cam_out
   type(physics_state), intent(inout) :: state
   type(physics_tend ), intent(inout) :: tend
   type(pbuf_fld),      intent(inout) :: pbuf(pbuf_size_max)

   type(check_tracers_data):: tracerint             ! tracer mass integrals and cummulative boundary fluxes

!
!---------------------------Local workspace-----------------------------
!
   type(physics_ptend)     :: ptend               ! indivdual parameterization tendencies

   integer  :: nstep                              ! current timestep number
   real(r8) :: zero(pcols)                        ! array of zeros

   integer :: lchnk                                ! chunk identifier
   integer :: ncol                                 ! number of atmospheric columns
   integer i,k,m                 ! Longitude, level indices
   integer :: yr, mon, day, tod       ! components of a date
   integer :: ixcldice, ixcldliq      ! constituent indices for cloud liquid and ice water.

   logical :: labort                            ! abort flag

   real(r8) tvm(pcols,pver)           ! virtual temperature
   real(r8) prect(pcols)              ! total precipitation
   real(r8) surfric(pcols)            ! surface friction velocity
   real(r8) obklen(pcols)             ! Obukhov length
   real(r8) :: fh2o(pcols)            ! h2o flux to balance source from methane chemistry
   real(r8) :: tmp_q     (pcols,pver) ! tmp space
   real(r8) :: tmp_cldliq(pcols,pver) ! tmp space
   real(r8) :: tmp_cldice(pcols,pver) ! tmp space
   real(r8) :: tmp_t     (pcols,pver) ! tmp space

! physics buffer fields for total energy and mass adjustment
   integer itim, ifld
   real(r8), pointer, dimension(:  ) :: teout
   real(r8), pointer, dimension(:,:) :: tini
   real(r8), pointer, dimension(:,:) :: cld
   real(r8), pointer, dimension(:,:) :: kvt
   real(r8), pointer, dimension(:,:) :: qini
   real(r8), pointer, dimension(:,:) :: cldliqini
   real(r8), pointer, dimension(:,:) :: cldiceini
   real(r8), pointer, dimension(:,:) :: dtcore
   real(r8), pointer, dimension(:,:) :: ast     ! relative humidity cloud fraction 
!
!-----------------------------------------------------------------------
!
   lchnk = state%lchnk
   ncol  = state%ncol

   nstep = get_nstep()


   if( .not. is_bypass_tphysac() ) then
      call tphysac_parametrized( &
                    ztodt,   pblh,    qpert,   tpert,  cam_in,  &
                    sgh,     sgh30,                                     &
                    cam_out,  state,   tend,    pbuf,           &
                    fsds &
              )
   end if

end subroutine tphysac

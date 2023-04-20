   module machine_learning_model

   !----------------------------------------------- !
   ! Purpose:                                       !
   !                                                !
   ! CAM interface to the shallow convection scheme !
   !                                                !
   ! Author: D.B. Coleman                           !
   !         Sungsu Park. Jan. 2010.                !
   !                                                !
   !----------------------------------------------- !

   use, intrinsic :: ieee_arithmetic,   only : ieee_is_finite, ieee_is_nan
   use shr_kind_mod,      only : r8=>shr_kind_r8
   use physconst,         only : cpair, cpairv, zvir, rga
   use ppgrid,            only : pver, pcols, pverp
   use zm_conv,           only : zm_conv_evap
   use cam_history,       only : outfld, addfld, add_default, phys_decomp
   use phys_control,      only : phys_getopts
   use torch_plugin,      only : predict_model, predict_model_v2, predict_model_v3, predict_model_v6
   use camsrfexch_types,only: cam_out_t, cam_in_t
   use cam_logfile,       only : iulog
   use phys_grid,       only: get_rlat_all_p, get_rlon_all_p
   use machine_learning_model_config

   implicit none
   private                 
   save

   public :: &
             ml_model_register,      & ! Register fields in physics buffer
             ml_model_init,          & ! Initialize shallow module
             ml_model_tend             ! Return tendencies

   ! The following namelist variable controls which shallow convection package is used.
   !        'Hack'   = Hack shallow convection (default)
   !        'UW'     = UW shallow convection by Sungsu Park and Christopher S. Bretherton
   !        'off'    = No shallow convection

   character(len=16) :: shallow_scheme      ! Default set in phys_control.F90, use namelist to change
   character(len=16) :: microp_scheme       ! Microphysics scheme
   logical           :: history_budget      ! Output tendencies and state variables for CAM4 T, qv, ql, qi
   integer           :: history_budget_histfile_num ! output history file number for budget fields

   ! Physics buffer indices 
   integer    ::     icwmrsh_idx    = 0  
   integer    ::      rprdsh_idx    = 0 
   integer    ::     rprdtot_idx    = 0 
   integer    ::      cldtop_idx    = 0 
   integer    ::      cldbot_idx    = 0 
   integer    ::        cush_idx    = 0 
   integer    :: nevapr_shcu_idx    = 0
   integer    ::       shfrc_idx    = 0 
   integer    ::         cld_idx    = 0 
   integer    ::      concld_idx    = 0
   integer    ::      rprddp_idx    = 0
   integer    ::         tke_idx    = 0

   integer    ::      solin_idx     = 0 


   integer :: & ! field index in physics buffer
      sh_flxprc_idx, &
      sh_flxsnw_idx, &
      sh_cldliq_idx, &
      sh_cldice_idx

   contains


  logical function is_finite(val)
        real(r8), intent(in) :: val

        is_finite = .true.

        !if(.not. ieee_is_finite( val )) then 
        if(1._r8/val == 0._r8) then 
            is_finite = .false.
            return
        end if

        !if(.not. ieee_is_nan( val )) then 
        if(.not.val<1._r8 .and. .not.val>1._r8 .and. .not.val==1._r8) then
            is_finite = .false.
            return
        end if
  end function



  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine ml_model_register

  !-------------------------------------------------- !
  ! Purpose : Register fields with the physics buffer !
  !-------------------------------------------------- !

  use phys_buffer,    only : pbuf_times, pbuf_add
  implicit none

  !call pbuf_add( 'SOLIN'   , 'physpkg' ,  1,  1,     1,        solin_idx )


  !call phys_getopts( shallow_scheme_out = shallow_scheme, microp_scheme_out = microp_scheme )

  !call pbuf_add( 'ICWMRSH'  , 'physpkg' ,  1,  pver,  1,       icwmrsh_idx )
  !call pbuf_add( 'RPRDSH'   , 'physpkg' ,  1,  pver,  1,        rprdsh_idx )
  !call pbuf_add( 'RPRDTOT'  , 'physpkg' ,  1,  pver,  1,       rprdtot_idx )
  !call pbuf_add( 'CLDTOP'   , 'physpkg' ,  1,  1,     1,        cldtop_idx )
  !call pbuf_add( 'CLDBOT'   , 'physpkg' ,  1,  1,     1,        cldbot_idx )
  !call pbuf_add( 'cush'     , 'global'  ,  1,  1,     pbuf_times, cush_idx ) 	
  !call pbuf_add( 'NEVAPR_SHCU', 'physpkg', 1,  pver,  1,   nevapr_shcu_idx )

  !if( shallow_scheme .eq. 'UW' ) then
  !    call pbuf_add( 'shfrc'  ,  'physpkg' ,  1,  pver,  1,  shfrc_idx )
  !endif

  !call pbuf_add('SH_FLXPRC', 'physpkg', 1, pverp, 1, sh_flxprc_idx)  ! shallow interface gbm flux_convective_cloud_rain+snow (kg/m2/s)
  !call pbuf_add('SH_FLXSNW', 'physpkg', 1, pverp, 1, sh_flxsnw_idx)  ! shallow interface gbm flux_convective_cloud_snow (kg/m2/s)
  !call pbuf_add('SH_CLDLIQ', 'physpkg', 1, pver,  1, sh_cldliq_idx)  ! shallow gbm cloud liquid water (kg/kg)
  !call pbuf_add('SH_CLDICE', 'physpkg', 1, pver,  1, sh_cldice_idx)  ! shallow gbm cloud ice water (kg/kg)

  end subroutine ml_model_register



  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine ml_model_init(pref_edge)

  !------------------------------------------------------------------------------- !
  ! Purpose : Declare output fields, and initialize variables needed by convection !
  !------------------------------------------------------------------------------- !

  use cam_history,       only : addfld, add_default, phys_decomp
  use ppgrid,            only : pcols, pver
  use hk_conv,           only : mfinti
  use uw_conv,           only : init_uw_conv
  use uwshcu,            only : init_uwshcu
  use physconst,         only : rair, gravit, latvap, rhoh2o, rh2o, zvir, tmelt, &
                                cappa, epsilo, latice, mwdry, mwh2o
  use pmgrid,            only : plev, plevp
  use spmd_utils,        only : masterproc
  use abortutils,        only : endrun
  use phys_control,      only : cam_physpkg_is
  use phys_buffer,       only : pbuf_get_fld_idx 

  implicit none

  real(r8), intent(in)       :: pref_edge(plevp)        ! Reference pressures at interfaces
    
  integer limcnv                                   ! Top interface level limit for convection
  integer k
  character(len=16)          :: eddy_scheme
    
  ! ------------------------------------------------- !
  ! Variables for detailed abalysis of UW-ShCu scheme !
  ! ------------------------------------------------- !

  !call addfld( 'qt_pre_Cu    ', 'kg/kg'   ,  pver ,  'I' , 'qt_preCU'                                         ,  phys_decomp )
  !call addfld( 'sl_pre_Cu    ', 'J/kg'    ,  pver ,  'I' , 'sl_preCU'                                         ,  phys_decomp )
  !call addfld( 'slv_pre_Cu   ', 'J/kg'    ,  pver ,  'I' , 'slv_preCU'                                        ,  phys_decomp )
  !call addfld( 'u_pre_Cu     ', 'm/s'     ,  pver ,  'I' , 'u_preCU'                                          ,  phys_decomp )
  !call addfld( 'v_pre_Cu     ', 'm/s'     ,  pver ,  'I' , 'v_preCU'                                          ,  phys_decomp )
  !call addfld( 'qv_pre_Cu    ', 'kg/kg'   ,  pver ,  'I' , 'qv_preCU'                                         ,  phys_decomp )
  !call addfld( 'ql_pre_Cu    ', 'kg/kg'   ,  pver ,  'I' , 'ql_preCU'                                         ,  phys_decomp )
  !call addfld( 'qi_pre_Cu    ', 'kg/kg'   ,  pver ,  'I' , 'qi_preCU'                                         ,  phys_decomp )
  !call addfld( 't_pre_Cu     ', 'K'       ,  pver ,  'I' , 't_preCU'                                          ,  phys_decomp )
  !call addfld( 'rh_pre_Cu    ', '%'       ,  pver ,  'I' , 'rh_preCU'                                         ,  phys_decomp )

  !call addfld( 'qt_aft_Cu    ', 'kg/kg'   ,  pver ,  'I' , 'qt_afterCU'                                       ,  phys_decomp )
  !call addfld( 'sl_aft_Cu    ', 'J/kg'    ,  pver ,  'I' , 'sl_afterCU'                                       ,  phys_decomp )
  !call addfld( 'slv_aft_Cu   ', 'J/kg'    ,  pver ,  'I' , 'slv_afterCU'                                      ,  phys_decomp )
  !call addfld( 'u_aft_Cu     ', 'm/s'     ,  pver ,  'I' , 'u_afterCU'                                        ,  phys_decomp )
  !call addfld( 'v_aft_Cu     ', 'm/s'     ,  pver ,  'I' , 'v_afterCU'                                        ,  phys_decomp )
  !call addfld( 'qv_aft_Cu    ', 'kg/kg'   ,  pver ,  'I' , 'qv_afterCU'                                       ,  phys_decomp )
  !call addfld( 'ql_aft_Cu    ', 'kg/kg'   ,  pver ,  'I' , 'ql_afterCU'                                       ,  phys_decomp )
  !call addfld( 'qi_aft_Cu    ', 'kg/kg'   ,  pver ,  'I' , 'qi_afterCU'                                       ,  phys_decomp )
  !call addfld( 't_aft_Cu     ', 'K'       ,  pver ,  'I' , 't_afterCU'                                        ,  phys_decomp )
  !call addfld( 'rh_aft_Cu    ', '%'       ,  pver ,  'I' , 'rh_afterCU'                                       ,  phys_decomp )

  !call addfld( 'tten_Cu      ', 'K/s'     ,  pver ,  'I' , 'Temprtaure tendency by cumulus convection'        ,  phys_decomp )
  !call addfld( 'rhten_Cu     ', '%/s'     ,  pver ,  'I' , 'RH tendency by cumumus convection'                ,  phys_decomp )

  ! ------------------------------------------- !
  ! Common Output for Shallow Convection Scheme !
  ! ------------------------------------------- !

  !call addfld( 'CMFDT   '     , 'K/s     ',  pver ,  'A' , 'T tendency - shallow convection'                           ,  phys_decomp )
  !call addfld( 'CMFDQ   '     , 'kg/kg/s ',  pver ,  'A' , 'QV tendency - shallow convection'                          ,  phys_decomp )
  !call addfld( 'CMFDLIQ '     , 'kg/kg/s ',  pver ,  'A' , 'Cloud liq tendency - shallow convection'                   ,  phys_decomp )
  !call addfld( 'CMFDICE '     , 'kg/kg/s ',  pver ,  'A' , 'Cloud ice tendency - shallow convection'                   ,  phys_decomp )
  !call addfld( 'CMFDQR  '     , 'kg/kg/s ',  pver ,  'A' , 'Q tendency - shallow convection rainout'                   ,  phys_decomp )
  !call addfld( 'EVAPTCM '     , 'K/s     ',  pver ,  'A' , 'T tendency - Evaporation/snow prod from Hack convection'   ,  phys_decomp )
  !call addfld( 'FZSNTCM '     , 'K/s     ',  pver ,  'A' , 'T tendency - Rain to snow conversion from Hack convection' ,  phys_decomp )
  !call addfld( 'EVSNTCM '     , 'K/s     ',  pver ,  'A' , 'T tendency - Snow to rain prod from Hack convection'       ,  phys_decomp )
  !call addfld( 'EVAPQCM '     , 'kg/kg/s ',  pver ,  'A' , 'Q tendency - Evaporation from Hack convection'             ,  phys_decomp )
  !call addfld( 'QC      '     , 'kg/kg/s ',  pver ,  'A' , 'Q tendency - shallow convection LW export'                 ,  phys_decomp )
  !call addfld( 'PRECSH  '     , 'm/s     ',  1,      'A' , 'Shallow Convection precipitation rate'                     ,  phys_decomp )
  !call addfld( 'CMFMC   '     , 'kg/m2/s ',  pverp,  'A' , 'Moist shallow convection mass flux'                        ,  phys_decomp )
  !call addfld( 'CMFSL   '     , 'W/m2    ',  pverp,  'A' , 'Moist shallow convection liquid water static energy flux'  ,  phys_decomp )
  !call addfld( 'CMFLQ   '     , 'W/m2    ',  pverp,  'A' , 'Moist shallow convection total water flux'                 ,  phys_decomp )
  !call addfld( 'CIN     '     , 'J/kg    ',  1    ,  'A' , 'Convective inhibition'                                     ,  phys_decomp )
  !call addfld( 'CBMF    '     , 'kg/m2/s ',  1    ,  'A' , 'Cloud base mass flux'                                      ,  phys_decomp )
  !call addfld( 'CLDTOP  '     , '1       ',  1    ,  'I' , 'Vertical index of cloud top'                               ,  phys_decomp )
  !call addfld( 'CLDBOT  '     , '1       ',  1    ,  'I' , 'Vertical index of cloud base'                              ,  phys_decomp )
  !call addfld( 'PCLDTOP '     , '1       ',  1    ,  'A' , 'Pressure of cloud top'                                     ,  phys_decomp )
  !call addfld( 'PCLDBOT '     , '1       ',  1    ,  'A' , 'Pressure of cloud base'                                    ,  phys_decomp )

  !call addfld( 'FREQSH '      , 'fraction',  1    ,  'A' , 'Fractional occurance of shallow convection'                ,  phys_decomp )
                                                                                                                    
  !call addfld( 'HKFLXPRC'     , 'kg/m2/s ',  pverp,  'A' , 'Flux of precipitation from HK convection'                  ,  phys_decomp )
  !call addfld( 'HKFLXSNW'     , 'kg/m2/s ',  pverp,  'A' , 'Flux of snow from HK convection'                           ,  phys_decomp )
  !call addfld( 'HKNTPRPD'     , 'kg/kg/s ',  pver ,  'A' , 'Net precipitation production from HK convection'           ,  phys_decomp )
  !call addfld( 'HKNTSNPD'     , 'kg/kg/s ',  pver ,  'A' , 'Net snow production from HK convection'                    ,  phys_decomp )
  !call addfld( 'HKEIHEAT'     , 'W/kg'    ,  pver ,  'A' , 'Heating by ice and evaporation in HK convection'           ,  phys_decomp )

  !if( shallow_scheme .eq. 'UW' ) then
  !   call addfld( 'UWFLXPRC'     , 'kg/m2/s ',  pverp,  'A' , 'Flux of precipitation from UW shallow convection'          ,  phys_decomp )
  !   call addfld( 'UWFLXSNW'     , 'kg/m2/s ',  pverp,  'A' , 'Flux of snow from UW shallow convection'                   ,  phys_decomp )
  !end if

  !call add_default( 'CMFDT   ', 1, ' ' )
  !call add_default( 'CMFDQ   ', 1, ' ' )
  !call add_default( 'CMFDQR  ', 1, ' ' )
  !call add_default( 'QC      ', 1, ' ' )
  !call add_default( 'PRECSH  ', 1, ' ' )
  !call add_default( 'CMFMC   ', 1, ' ' )
  !call add_default( 'FREQSH  ', 1, ' ' ) 

  !call phys_getopts( eddy_scheme_out = eddy_scheme, history_budget_out = history_budget, &
  !                   history_budget_histfile_num_out = history_budget_histfile_num)

  !if( history_budget ) then
  !    call add_default( 'CMFDLIQ  ', history_budget_histfile_num, ' ' )
  !    call add_default( 'CMFDICE  ', history_budget_histfile_num, ' ' )
  !    if( cam_physpkg_is('cam3') .or. cam_physpkg_is('cam4') ) then
  !       call add_default( 'EVAPQCM  ', history_budget_histfile_num, ' ' )
  !       call add_default( 'EVAPTCM  ', history_budget_histfile_num, ' ' )
  !    end if
  !    if( history_budget_histfile_num > 1 ) then
  !       call add_default( 'CMFDT   ', history_budget_histfile_num, ' ' )
  !       call add_default( 'CMFDQ   ', history_budget_histfile_num, ' ' )
  !    end if
  !end if

  !select case (shallow_scheme)

  !case('off')  ! None

  !   if( masterproc ) write(iulog,*) 'convect_shallow_init: shallow convection OFF'
  !   continue

  !case('Hack') ! Hack scheme

  !   if( masterproc ) write(iulog,*) 'convect_shallow_init: Hack shallow convection'
  ! ! Limit shallow convection to regions below 40 mb
  ! ! Note this calculation is repeated in the deep convection interface
  !   if( pref_edge(1) >= 4.e3_r8 ) then
  !       limcnv = 1
  !   else
  !       do k = 1, plev
  !          if( pref_edge(k) < 4.e3_r8 .and. pref_edge(k+1) >= 4.e3_r8 ) then
  !              limcnv = k
  !              goto 10
  !          end if
  !       end do
  !       limcnv = plevp
  !   end if
  !10 continue
  !   if( masterproc ) then
  !       write(iulog,*) 'MFINTI: Convection will be capped at intfc ', limcnv, ' which is ', pref_edge(limcnv), ' pascals'
  !   end if
  !   call mfinti( rair, cpair, gravit, latvap, rhoh2o, limcnv) ! Get args from inti.F90

  !case('UW') ! Park and Bretherton shallow convection scheme

  !   if( masterproc ) write(iulog,*) 'convect_shallow_init: UW shallow convection scheme (McCaa)'
  !   if( eddy_scheme .ne. 'diag_TKE' ) then
  !       write(iulog,*) 'ERROR: shallow convection scheme ', shallow_scheme, ' is incompatible with eddy scheme ', eddy_scheme
  !       call endrun( 'convect_shallow_init: shallow_scheme and eddy_scheme are incompatible' )
  !   endif
  !   call init_uwshcu( r8, latvap, cpair, latice, zvir, rair, gravit, mwh2o/mwdry )

  !   tke_idx =  pbuf_get_fld_idx('tke')

  !end select

  cld_idx     = pbuf_get_fld_idx('CLD')
  concld_idx  = pbuf_get_fld_idx('CONCLD')
  rprddp_idx  = pbuf_get_fld_idx('RPRDDP')
  rprdsh_idx  = pbuf_get_fld_idx('RPRDSH')
  rprdtot_idx = pbuf_get_fld_idx('RPRDTOT')

  end subroutine ml_model_init


  !=============================================================================== !
  !                                                                                !
  !=============================================================================== !

  subroutine ml_model_tend( ztodt, state   , ptend_all, pbuf, &
                            fsns,    fsnt,    flns,    flnt, fsds, flds, solin, &
                            cam_in, cam_out, &
                            prect, precc, precl, &
                            precsc, precsl, &
                            srfrad, soll, solld, sols, solsd, &
                            relhum )

   use cam_history,     only : outfld
   use physics_types,   only : physics_state, physics_ptend, physics_tend
   use physics_types,   only : physics_ptend_init, physics_tend_init, physics_update
   use physics_types,   only : physics_state_copy
   use physics_types,   only : physics_ptend_sum
   use phys_buffer,     only : pbuf_size_max, pbuf_fld, pbuf_old_tim_idx, pbuf_get_fld_idx
   use constituents,    only : pcnst, qmin, cnst_get_ind, cnst_get_type_byind
   use hk_conv,         only : cmfmca
   use uw_conv,         only : compute_uw_conv
   use uwshcu,          only : compute_uwshcu_inv
   use time_manager,    only : get_nstep, is_first_step
   use wv_saturation,   only : fqsatd, aqsat
   use physconst,       only : latice, latvap, rhoh2o
   use spmd_utils,      only : masterproc
   
   use diag_ml, only : ml_diagnostics_tend_post
   use ml_solin,        only : ml_calc_solin
   
   use check_energy,    only: check_energy_get_integrals
   use physconst,       only: cpair
   use dycore,          only: dycore_is
   use abortutils,        only : endrun


   implicit none

   ! ---------------------- !
   ! Input-Output Arguments !
   ! ---------------------- !

   type(physics_state), intent(inout)    :: state                           ! Physics state variables
   real(r8),            intent(in)    :: ztodt                           ! 2 delta-t  [ s ]
   !real(r8),            intent(in)    :: pblht(pcols)                    ! PBL height [ m ]

   type(pbuf_fld),      intent(inout), dimension(pbuf_size_max) :: pbuf  ! Physics buffer
   !real(r8),            intent(inout) :: qpert(pcols,pcnst)              ! PBL perturbation specific humidity
   !real(r8),            intent(inout) :: cmfmc(pcols,pverp)              ! Moist deep + shallow convection cloud mass flux [ kg/s/m2 ]
   !real(r8),            intent(inout) :: qc(pcols,pver)                  ! dq/dt due to export of cloud water into environment by shallow and deep convection [ kg/kg/s ]
   !real(r8),            intent(inout) :: rliq(pcols)                     ! Vertical integral of qc [ m/s ]
   
   real(r8),            intent(inout) :: fsns(pcols)                     ! 
   real(r8),            intent(inout) :: fsnt(pcols)                     ! 
   real(r8),            intent(inout) :: flns(pcols)                     ! 
   real(r8),            intent(inout) :: flnt(pcols)                     ! 
   real(r8),            intent(inout) :: fsds(pcols)                     ! 
   real(r8),            intent(inout) :: solin(pcols)                     ! 

   type(cam_out_t),     intent(inout) :: cam_out
   type(cam_in_t),      intent(in)    :: cam_in

   type(physics_ptend), intent(out)   :: ptend_all                       ! Indivdual parameterization tendencies
   real(r8),            intent(out)   :: prect(pcols)                    ! 
   real(r8),            intent(out)   :: precc(pcols)                    ! 
   real(r8),            intent(out)   :: precl(pcols)                    ! 
   real(r8),            intent(out)   :: precsc(pcols)                   ! 
   real(r8),            intent(out)   :: precsl(pcols)                   ! 

   real(r8), intent(out), dimension(pcols) :: flds                       ! 
   real(r8), intent(out), dimension(pcols) :: srfrad
   real(r8), intent(out), dimension(pcols) :: soll
   real(r8), intent(out), dimension(pcols) :: sols
   real(r8), intent(out), dimension(pcols) :: solld
   real(r8), intent(out), dimension(pcols) :: solsd

   real(r8), intent(in), dimension(pcols, pver) :: relhum


   ! --------------- !
   ! Local Variables ! 
   ! --------------- !
   integer  :: i, k, m
   integer  :: n, x
   integer  :: ilon                                                      ! Global longitude index of a column
   integer  :: ilat                                                      ! Global latitude  index of a column
   integer  :: lchnk                                                     ! Chunk identifier
   integer  :: ncol                                                      ! Number of atmospheric columns
   integer  :: nstep                                                     ! Current time step index
   integer  :: ixcldice, ixcldliq                                        ! Constituent indices for cloud liquid and ice water.
   integer  :: ixnumice, ixnumliq                                        ! Constituent indices for cloud liquid and ice number concentration

   !real(r8) :: ftem(pcols,pver)                                          ! Temporary workspace for outfld variables
   !real(r8) :: cnt2(pcols)                                               ! Top level of shallow convective activity
   !real(r8) :: cnb2(pcols)                                               ! Bottom level of convective activity
   !real(r8) :: tpert(pcols)                                              ! PBL perturbation theta
   !real(r8) :: ntprprd(pcols,pver)                                       ! Net precip production in layer
   !real(r8) :: ntsnprd(pcols,pver)                                       ! Net snow   production in layer
   !real(r8) :: tend_s_snwprd(pcols,pver)                                 ! Heating rate of snow production
   !real(r8) :: tend_s_snwevmlt(pcols,pver)                               ! Heating rate of evap/melting of snow
   !real(r8) :: slflx(pcols,pverp)                                        ! Shallow convective liquid water static energy flux
   !real(r8) :: qtflx(pcols,pverp)                                        ! Shallow convective total water flux
   !real(r8) :: cmfdqs(pcols, pver)                                       ! Shallow convective snow production
   real(r8) :: zero(pcols)                                               ! Array of zeros
   !real(r8) :: cbmf(pcols)                                               ! Shallow cloud base mass flux [ kg/s/m2 ]
   !real(r8) :: freqsh(pcols)                                             ! Frequency of shallow convection occurence
   !real(r8) :: pcnt(pcols)                                               ! Top    pressure level of shallow + deep convective activity
   !real(r8) :: pcnb(pcols)                                               ! Bottom pressure level of shallow + deep convective activity
   !real(r8) :: cmfsl(pcols,pverp )                                       ! Convective flux of liquid water static energy
   !real(r8) :: cmflq(pcols,pverp )                                       ! Convective flux of total water in energy unit
   
   !real(r8) :: ftem_preCu(pcols,pver)                                    ! Saturation vapor pressure after shallow Cu convection
   !real(r8) :: tem2(pcols,pver)                                          ! Saturation specific humidity and RH
   !real(r8) :: t_preCu(pcols,pver)                                       ! Temperature after shallow Cu convection
   !real(r8) :: tten(pcols,pver)                                          ! Temperature tendency after shallow Cu convection
   !real(r8) :: rhten(pcols,pver)                                         ! RH tendency after shallow Cu convection
   !real(r8) :: iccmr_UW(pcols,pver)                                      ! In-cloud Cumulus LWC+IWC [ kg/m2 ]
   !real(r8) :: icwmr_UW(pcols,pver)                                      ! In-cloud Cumulus LWC     [ kg/m2 ]
   !real(r8) :: icimr_UW(pcols,pver)                                      ! In-cloud Cumulus IWC     [ kg/m2 ]
   !real(r8) :: ptend_tracer(pcols,pver,pcnst)                            ! Tendencies of tracers
   !real(r8) :: sum1, sum2, sum3, pdelx 

   !real(r8), dimension(pcols,pver) :: sl, qt, slv
   !real(r8), dimension(pcols,pver) :: sl_preCu, qt_preCu, slv_preCu

   !type(physics_state) :: state1                                         ! Locally modify for evaporation to use, not returned
   type(physics_ptend) :: ptend_loc                                      ! Local tendency from processes, added up to return as ptend_all
   type(physics_tend ) :: tend                                           ! Physics tendencies ( empty, needed for physics_update call )

   integer itim, ifld

   real(r8), dimension(pcols,pver) :: spec_hum
   real(r8), dimension(pcols,pver) :: pteq
   real(r8), dimension(pcols,pver) :: ptecldliq
   real(r8), dimension(pcols,pver) :: ptecldice
   real(r8), dimension(pcols,pver) :: pttend
   real(r8), dimension(pcols,pver) :: dcq
   real(r8), dimension(pcols,pver) :: dtcond

   real(r8), dimension(pcols,pver) :: qrs
   real(r8), dimension(pcols,pver) :: qrl
   real(r8), dimension(pcols,pver) :: cld 
   real(r8), dimension(pcols,pver) :: concld
  


   real(r8) :: rtdt
   real(r8) :: heat_glob         ! global energy integral (FV only)
   real(r8) :: correction_term = 0._r8
   real(r8) :: corrected_pttend = 0._r8
   real(r8) :: qq = 0
   logical  :: debug_log = .true.
   logical  :: predict_cloud_tend = .false.

   !for debug purposes
   real(r8) :: fsns_pre
   real(r8) :: flns_pre
   real(r8) :: fsnt_pre
   real(r8) :: flnt_pre
   real(r8) :: fsds_pre


   ! ----------------------- !
   ! Main Computation Begins ! 
   ! ----------------------- !
   !if( masterproc ) write(iulog,*) 'ml_model_tend start...'

   rtdt = 1._r8/ztodt
   zero  = 0._r8
   nstep = get_nstep()
   lchnk = state%lchnk
   ncol  = state%ncol
  
   !call physics_state_copy( state, state1 )   ! Copy state to local state1.
   call physics_ptend_init( ptend_loc )       ! Initialize local ptend type
   call physics_ptend_init( ptend_all )       ! Initialize output ptend type
   call physics_tend_init( tend )             ! Tendency type here is a null place holder

   
   ! Associate pointers with physics buffer fields

   !itim   =  pbuf_old_tim_idx()
   !cld    => pbuf(cld_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

   !itim   =  pbuf_old_tim_idx()
   !concld => pbuf(concld_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,itim)

   !icwmr  => pbuf(icwmrsh_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,1)

   !rprddp => pbuf(rprddp_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,1)

   !rprdsh => pbuf(rprdsh_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,1)

   !evapcsh => pbuf(nevapr_shcu_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,1)

   !cnt    => pbuf(cldtop_idx)%fld_ptr(1,1:pcols,1,lchnk,1)

   !cnb    => pbuf(cldbot_idx)%fld_ptr(1,1:pcols,1,lchnk,1)

   !if( convect_shallow_use_shfrc() ) then ! Park-Bretherton UW Shallow Convection Schemes
   !    shfrc => pbuf(shfrc_idx)%fld_ptr(1,1:pcols,1:pver,lchnk,1)
   !endif

   !solin_buf    => pbuf(solin_idx)%fld_ptr(1,1:pcols,1,lchnk,1)


   ! Initialization

   !tpert(:ncol)         = 0._r8
   !qpert(:ncol,2:pcnst) = 0._r8

   call cnst_get_ind( 'CLDLIQ', ixcldliq )
   call cnst_get_ind( 'CLDICE', ixcldice )

   call ml_calc_solin(state, solin)


   !call sanitize_state(state, fsns, fsnt, flns, flnt, fsds, solin )

   do i=1,ncol
       if(.false.) then
           call predict_model_v2( & 
               state%q(i,:,1), &
               state%t(i,:), &
               state%u(i,:), &
               state%v(i,:), &
               state%omega(i,:), &
               state%zm(i,:pver) + state%phis(i)*rga, & !repeat whats done for Z3 in cam_diagnostics.F90
               state%ps(i), &
               cam_in%ts(i), & ! surface temp cam_in
                       !
               solin(i), & ! solin from direct calculation 
               cam_in%shf(i), & ! sensible heat flux from cam_in
               cam_in%lhf(i), & ! latend heat flux from cam_in

               !added for v2
               cam_in%landfrac(i), &
               cam_in%ocnfrac(i), &
               cam_in%icefrac(i), &

               !inout
               fsns(i),   & ! from tphysbc
               flns(i),   & ! from tphysbc
               fsnt(i),   & ! from tphysbc
               flnt(i),   & ! from tphysbc 
               fsds(i),   & ! from tphysbc 

               !out
               flds(i),   & ! 
               srfrad(i), & ! from tphysbc 
               soll(i),   & ! 
               sols(i),   & ! 
               solld(i),  & ! 
               solsd(i),  & ! 
                      ! 
               pteq(i,:),   & ! total physics moist tendencies
               pttend(i,:), & ! total physics heating tendencies
               !dcq(i,:),    & ! model v1 -- total moist tendencies from moist processes
               !dtcond(i,:), & ! model v1 -- total heating tendencies from moist processes
               !qrs(i,:),    & ! model v1 -- total heating tendencies due to SW radiation
               !qrl(i,:),    & ! model v1 -- total heating tendencies due to LW radiation
               !cld(i,:),    & ! model v1 -- total cloud fraction
               !concld(i,:), & ! model v1 -- total convective cloud fraction
           
               prect(i), & ! prect - not used aside from output
               precc(i), &
               precl(i), &
               precsc(i), &
               precsl(i), &
               .false. &
              )
       else if(.not. predict_cloud_tend) then
           call predict_model_v3( & 
               !relhum(i,:), &
               state%q(i,:,1), &
               state%t(i,:), &
               state%u(i,:), &
               state%v(i,:), &
               state%omega(i,:), &
               state%zm(i,:pver) + state%phis(i)*rga, & !repeat whats done for Z3 in cam_diagnostics.F90
               state%ps(i), &
               cam_in%ts(i), & ! surface temp cam_in
                       !
               solin(i), & ! solin from direct calculation 
               cam_in%shf(i), & ! sensible heat flux from cam_in
               cam_in%lhf(i), & ! latend heat flux from cam_in

               !added for v2
               !WARN: fixed land/ocean fraction
               !0.0_r8, & !cam_in%landfrac(i), &
               !1.0_r8, & !cam_in%ocnfrac(i), &
               cam_in%landfrac(i), &
               cam_in%ocnfrac(i), &
               cam_in%icefrac(i), &

               !out
               fsns(i),   & ! from tphysbc
               flns(i),   & ! from tphysbc
               fsnt(i),   & ! from tphysbc
               flnt(i),   & ! from tphysbc 
               fsds(i),   & ! from tphysbc 
               flds(i),   & ! 
               srfrad(i), & ! from tphysbc 
               soll(i),   & ! 
               sols(i),   & ! 
               solld(i),  & ! 
               solsd(i),  & ! 
                            ! 
               pteq(i,:),   & ! total physics moist tendencies
               pttend(i,:), & ! total physics heating tendencies
           
               prect(i), & ! prect - not used aside from output
               precc(i), &
               precl(i), &
               precsc(i), &
               precsl(i), &
               .false. &
              )
       else
           call predict_model_v6( & 
               !relhum(i,:), &
               state%q(i,:,1), &
               state%t(i,:), &
               state%u(i,:), &
               state%v(i,:), &
               state%omega(i,:), &
               state%zm(i,:pver) + state%phis(i)*rga, & !repeat whats done for Z3 in cam_diagnostics.F90
               state%ps(i), &
               cam_in%ts(i), & ! surface temp cam_in
                       !
               solin(i), & ! solin from direct calculation 
               cam_in%shf(i), & ! sensible heat flux from cam_in
               cam_in%lhf(i), & ! latend heat flux from cam_in

               !added for v2
               !WARN: fixed land/ocean fraction
               !0.0_r8, & !cam_in%landfrac(i), &
               !1.0_r8, & !cam_in%ocnfrac(i), &
               cam_in%landfrac(i), &
               cam_in%ocnfrac(i), &
               cam_in%icefrac(i), &

               !out
               fsns(i),   & ! from tphysbc
               flns(i),   & ! from tphysbc
               fsnt(i),   & ! from tphysbc
               flnt(i),   & ! from tphysbc 
               fsds(i),   & ! from tphysbc 
               flds(i),   & ! 
               srfrad(i), & ! from tphysbc 
               soll(i),   & ! 
               sols(i),   & ! 
               solld(i),  & ! 
               solsd(i),  & ! 
                            ! 
               pteq(i,:),   & ! total physics moist tendencies
               pttend(i,:), & ! total physics heating tendencies
               ptecldliq(i,:),   & ! total physics cloud liquid tendencies
               ptecldice(i,:),   & ! total physics cloud ice tendencies
           
               prect(i), & ! prect - not used aside from output
               precc(i), &
               precl(i), &
               precsc(i), &
               precsl(i), &
               .false. &
              )
 
       end if

   end do
 

   !call sanitize_tend(pttend, pteq, fsns, fsnt, flns, flnt, fsds, solin )


   ptend_loc%name  = 'ltorch'
   ptend_loc%ls    = .true.
   ptend_loc%lq(1) = .true.

   

   if (is_radiation_from_ml()) then
       do i=1,ncol
           !variables that are added to cam_out not in srfxfer - so we do the same here. ugly lack of isolation...
           !note that we are not using srfrad, it could be used for stability and validation...

           call check_swrad(solin(i), fsns(i), "fsns", cam_in%landfrac(i)) 
           call check_swrad(solin(i), fsnt(i), "fsnt", cam_in%landfrac(i)) 
           call check_swrad(solin(i), fsds(i), "fsds", cam_in%landfrac(i)) 

           !flds(i)  = min(flds(i), 800._r8)  ! seemingly a good guess of upper bound
           !flds(i)  = max(flds(i), 0._r8)  ! check sign of fsns
           
           !call check_swrad(solin(i), sols(i),  "sols", cam_in%landfrac(i)) 
           !call check_swrad(solin(i), soll(i),  "soll", cam_in%landfrac(i)) 
           !call check_swrad(solin(i), solsd(i), "solsd", cam_in%landfrac(i)) 
           !call check_swrad(solin(i), solld(i), "solld", cam_in%landfrac(i)) 
           
           call constrain_direct_diffuse(fsds(i), soll(i), solld(i), sols(i), solsd(i))
          
           !flds(i)  = max(flds(i), 0._r8)
           !flns(i)  = max(flns(i), 0._r8)

           !for consistency we set srfrad to be a sum rather than the predicted value
           srfrad(i) = fsns(i) + flds(i)
           cam_out%srfrad(i) = srfrad(i)
           cam_out%netsw(i)  = fsns(i)  
           cam_out%flwds(i)  = flds(i) 
           cam_out%sols(i)   = sols(i)
           cam_out%solsd(i)  = solsd(i)
           cam_out%soll(i)   = soll(i)
           cam_out%solld(i)  = solld(i)

           !if(.not. is_finite( flds(i) )) then
           !    call endrun( 'flds is not finite!')
           !end if 

           !if(.not. is_finite( flns(i) )) then
           !    flns(i) = 0._r8
           !    !call endrun( 'flns is not finite!')
           !end if 

           !if(.not. is_finite( fsns(i) )) then
           !    call endrun( 'fsns is not finite!')
           !end if 

           !if(.not. is_finite( soll(i) )) then
           !    call endrun( 'soll is not finite!')
           !end if 

           !if(.not. is_finite( solld(i) )) then
           !    call endrun( 'solld is not finite!')
           !end if 

           do k=1,pver


               ptend_loc%s(i,k)   = pttend(i,k) * cpairv(i,k,state%lchnk) !convert it back to heat, see physics_types.F90: physics_update 
               ptend_loc%q(i,k,1) = pteq(i,k) !DF: / rtdt --> note that in physics_update: state%q(i,k,m) = state%q(i,k,m) + ptend%q(i,k,m) * dt
 
               if(predict_cloud_tend) then
                   ptend_loc%q(i,k,ixcldliq) = ptecldliq(i,k) 
                   ptend_loc%q(i,k,ixcldice) = ptecldice(i,k) 
               end if
               !if(.not. is_finite( ptend_loc%q(i,k,1) )) then
               !    ptend_loc%q(i,k,1) = 0._r8
               !end if

               !if(.not. is_finite( ptend_loc%s(i,k) )) then
               !    ptend_loc%s(i,k) = 0._r8
               !end if

               !if(.not. is_finite( ptend_loc%s(i,k) )) then
               !    call flush_state(state, i, k, ptend_loc) 
               !    call endrun( 'ptend_loc%s(i,k) is not finite!')
               !end if
               !if(.not. is_finite( ptend_loc%q(i,k,1) )) then
               !    call flush_state(state, i, k, ptend_loc) 
               !    call endrun( 'ptend_loc%q(i,k,1) is not finite!')
               !end if
           end do

   
       end do
   else
       do i=1,ncol
           do k =1,pver
               !see: models/atm/cam/src/physics/cam/cam_diagnostics.F90: 
               ! in the following we do the inverse direction
               !1170          dtcond(i,k,lchnk) = (state%s(i,k) - dtcond(i,k,lchnk))*rtdt / cpair
               !1185             dqcond(i,k,m,lchnk) = (state%q(i,k,m) - dqcond(i,k,m,lchnk))*rtdt
               ptend_loc%s(i,k)   = dtcond(i,k) * cpair / rtdt
               ptend_loc%q(i,k,1) = dcq(i,k) / rtdt

               if(predict_cloud_tend) then
                   ptend_loc%q(i,k,ixcldliq) = ptecldliq(i,k) 
                   ptend_loc%q(i,k,ixcldice) = ptecldice(i,k) 
               end if
     
           end do
       end do
   end if


   !call ml_diagnostics_tend_post(state, solin, pttend, pteq)


   !Zero tendencies of constituents that are not handled by the ML model
   !
   !if(.not. predict_cloud_tend) then
   !    do i=2,pcnst
   !        ptend_loc%q(:ncol,:pver,i) = 0._r8
   !        state%q(:ncol,:pver,i) = qmin(i)
   !    end do
   !end if

   !ptend_loc%s(:ncol,:pver) = 0._r8
   !ptend_loc%s(:,1:5) = 0._r8
   !ptend_loc%s(:,24:26) = ptend_loc%s(:,24:26)
   
   !ptend_loc%q(:,1:7,1) = 0._r8

 
   !ptend_loc%lu = .TRUE.
   !ptend_loc%lv = .TRUE.
   ptend_all%name  = 'ml_model'
   ptend_all%ls    = .TRUE.
   ptend_all%lq(1) = .TRUE.


   call physics_ptend_sum( ptend_loc, ptend_all, state )

   return

  end subroutine ml_model_tend



  subroutine flush_state(state, col, k, ptend) 
      use cam_logfile,       only : iulog
      use physics_types,   only : physics_state, physics_ptend, physics_tend


  implicit none

      type(physics_state), intent(in) :: state  
      type(physics_ptend), intent(in) :: ptend
      integer, intent(in) :: col,k


      write(iulog, *) "ptend-s:"
      write(iulog, *) ptend%s(col, :)

      write(iulog, *) "ptend-q:"
      write(iulog, *) ptend%q(col, :, 1)

  end subroutine flush_state




  subroutine check_swrad(solin, rad_var, rad_var_id, landfrac) 
      use cam_logfile,       only : iulog


   implicit none

   character (*), intent(in) :: rad_var_id
   real(r8),      intent(inout) :: rad_var
   real(r8),      intent(in) :: solin
   real(r8),      intent(in) :: landfrac

   !if(solin < 150._r8) then 
   !     rad_var = 0._r8
   !     return
   !end if

   if(rad_var > solin) then        
        write(iulog, *) "ML HIGH EXCEEDANCE: ", rad_var_id, " > solin: ", rad_var, " > ", solin
        write(iulog, *) "ML HIGH EXCEEDANCE: solin: ", solin
        write(iulog, *) "ML HIGH EXCEEDANCE: landfrac: ", landfrac
        rad_var  = min(rad_var, solin)  ! check sign of fsns
   end if

   if(rad_var < 0._r8) then        
        write(iulog, *) "ML LOW EXCEEDANCE: ", rad_var_id, " < 0._r8: ", rad_var, " > 0._r8"
        rad_var  = max(rad_var, 0._r8)  ! check sign of fsns
   end if


  end subroutine check_swrad

  subroutine constrain_direct_diffuse(fsds, soll, solld, sols, solsd)
      use cam_logfile,       only : iulog


   implicit none

   real(r8),      intent(inout) :: soll
   real(r8),      intent(inout) :: solld
   real(r8),      intent(inout) :: sols
   real(r8),      intent(inout) :: solsd
   real(r8),      intent(in) :: fsds

   real(r8) :: sumall
   real(r8) :: scaler

   !if(solin .eq. 0._r8) then
   if(fsds .eq. 0._r8) then
       soll  = 0._r8
       solld = 0._r8
       sols  = 0._r8
       solsd = 0._r8
       return
   end if


   soll  = max(soll , 0._r8)
   solld = max(solld, 0._r8)

   sols  = max(sols , 0._r8)
   solsd = max(solsd, 0._r8)

   sumall = soll + solld + sols + solsd
   
   if(sumall <= 0._r8) then
       sumall = fsds
       sols  = fsds
       solsd = 0._r8
       soll  = 0._r8
       solld = 0._r8
       return
   end if

   !if(sumall .eq. fsds) then
   !    return
   !end if


   scaler = fsds/sumall 

   soll  = soll  * scaler
   solld = solld * scaler

   sols  = sols  * scaler
   solsd = solsd * scaler

  end subroutine constrain_direct_diffuse




           


  end module machine_learning_model

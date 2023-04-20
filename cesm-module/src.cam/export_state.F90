module export_state
   use time_manager,    only: is_first_step
   use shr_kind_mod,      only : r8=>shr_kind_r8
   use phys_grid,       only: get_lat_p, get_lon_p 
   use camsrfexch_types, only: cam_in_t, cam_out_t
   use physics_types, only: physics_state, physics_tend, physics_state_copy, physics_tend_init
   use ppgrid,        only: pcols, pver, pverp, begchunk, endchunk
   use cam_history,   only: outfld, write_inithist, hist_fld_active
   use constituents,  only: pcnst, cnst_name, cnst_longname, cnst_cam_outfld, ptendnam, dmetendnam, apcnst, bpcnst, &
                            cnst_get_ind

   use machine_learning_model_config, only: is_output_state_vars_from_ctrl_run 

   !----------------------------------------------- !
   ! Purpose:                                       !
   !                                                !
   ! This module assists with the initial learning  !
   ! of a ML/AI model. It allows the user to plug   !
   ! calls to output selected state variables to    !
   ! history files before and after the location    !
   ! where a ML/AI model will be placed. The data   !
   ! can then be used to boost a ML/AI model.       !
   !                                                !
   !----------------------------------------------- !

   implicit none
   private
   save

   integer :: ixcldice, ixcldliq              ! constituent indices for cloud liquid and ice water.
   integer :: ixnumice, ixnumliq



   public :: &
            export_state_init, &
            record_state_before, &
            record_state_after, &
            record_state_after_additional, &
            physics_tend_copy
  contains

  subroutine export_state_init()
      use cam_history,        only: addfld, add_default, phys_decomp

  !-------------------------------------------------- !
  ! Purpose : Register fields that will be added to   !
  ! CAM's history output for the sake of training ML  !
  ! models.                                           !
  !-------------------------------------------------- !


  implicit none
        integer  :: history_budget_histfile_num    ! output history file number for budget fields

        if( .not. is_output_state_vars_from_ctrl_run() ) then
                return
        end if

        call cnst_get_ind( 'CLDLIQ', ixcldliq )
        call cnst_get_ind( 'CLDICE', ixcldice )
        !call cnst_get_ind( 'NUMLIQ', ixnumliq )
        !call cnst_get_ind( 'NUMICE', ixnumice )


        !
        ! BEFORE:
        !
        call addfld ('B_T     ','K       ',pver, 'A','Temperature',phys_decomp)
        call addfld ('B_U     ','m/s     ',pver, 'A','Zonal wind',phys_decomp)
        call addfld ('B_V     ','m/s     ',pver, 'A','Meridional wind',phys_decomp)
        call addfld ('B_Q     ','kg/kg   ',pver, 'A',cnst_longname(1),phys_decomp)
        call addfld ('B_OMEGA ','Pa/s    ',pver, 'A','Vertical velocity (pressure)',phys_decomp)        
        call addfld ('B_Z3    ','m       ',pver, 'A','Geopotential Height (above sea level)',phys_decomp)
        call addfld ('B_RELHUM','percent ',pver, 'A','Relative humidity',phys_decomp)

        call addfld ('B_CLDLIQ     ','kg/kg   ',pver, 'A','Grid box averaged cloud liquid amount before ML', phys_decomp)
        call addfld ('B_CLDICE     ','kg/kg   ',pver, 'A','Grid box averaged cloud ice amount before ML', phys_decomp)

        call addfld ('B_NUMLIQ     ','1/kg   ',pver, 'A','Grid box averaged cloud liquid number before ML', phys_decomp)
        call addfld ('B_NUMICE     ','1/kg   ',pver, 'A','Grid box averaged cloud ice number before ML', phys_decomp)


        call addfld ('B_PS    ','Pa      ',1,    'A','Surface pressure',phys_decomp)
        call addfld ('B_TS    ','K       ',1,    'A','Surface temperature (radiative)',phys_decomp)

        !solin
        
        call addfld ('B_SHFLX ','W/m2    ',1,    'A','Surface sensible heat flux',phys_decomp)
        call addfld ('B_LHFLX ','W/m2    ',1,    'A','Surface latent heat flux',phys_decomp)


        ! defaults
        call add_default ('B_T     '  , 1, ' ')
        call add_default ('B_U     '  , 1, ' ')
        call add_default ('B_V     '  , 1, ' ')
        call add_default ('B_Q     '  , 1, ' ')
        call add_default ('B_CLDLIQ'  , 1, ' ')
        call add_default ('B_CLDICE'  , 1, ' ')
        call add_default ('B_NUMLIQ'  , 1, ' ')
        call add_default ('B_NUMICE'  , 1, ' ')
        call add_default ('B_OMEGA '  , 1, ' ')
        call add_default ('B_Z3    '  , 1, ' ')
        call add_default ('B_RELHUM'  , 1, ' ')
        call add_default ('B_PS    '  , 1, ' ')
        call add_default ('B_TS    '  , 1, ' ')
        call add_default ('B_SHFLX '  , 1, ' ')
        call add_default ('B_LHFLX '  , 1, ' ')
         


        !
        ! AFTER
        !
        call addfld ('A_PTTEND'   ,'K/s     ',pver, 'A','T total physics tendency',phys_decomp)
        
        !A_PTEQ
        call addfld ('A_'//ptendnam(       1),  'kg/kg/s ',pver, 'A',trim(cnst_name(       1))//' total physics tendency ',phys_decomp)

        !A_PTECLDLIQ, A_PTECLDICE
        call addfld ('A_PTECLDLIQ',  'kg/kg/s ',pver, 'A','CLDLIQ total physics tendency ',phys_decomp)
        call addfld ('A_PTECLDICE',  'kg/kg/s ',pver, 'A','CLDICE total physics tendency ',phys_decomp)
        !TODO: add somewhere NUMLIQ and NUMICE tendencies for CAM5 and beyond...

        call addfld ('A_TFIX    ','K/s     ',1,    'A'     ,'T fixer (T equivalent of Energy correction)',phys_decomp)

        call add_default ('A_PTTEND'          , 1, ' ')
        call add_default ('A_'//ptendnam(       1), 1, ' ')
        call add_default ('A_TFIX    '    , 1, ' ')



        !
        ! Additional diagnostics
        ! Used for cases where we want to output CAM's original param and ML 
        !
        call addfld ('A_PRECL   ','m/s     ',1,    'A','Large-scale (stable) precipitation rate (liq + ice)',phys_decomp)
        call addfld ('A_PRECC   ','m/s     ',1,    'A','Convective precipitation rate (liq + ice)',phys_decomp)
        call addfld ('A_PRECT   ','m/s     ',1,    'A','Total (convective and large-scale) precipitation rate (liq + ice)',phys_decomp)
        call addfld ('A_PRECSL  ','m/s     ',1,    'A','Large-scale (stable) snow rate (water equivalent)' ,phys_decomp)
        call addfld ('A_PRECSC  ','m/s     ',1,    'A','Convective snow rate (water equivalent)',phys_decomp)


        ! defaults
        call add_default ('A_PRECL   ', 1, ' ')
        call add_default ('A_PRECC   ', 1, ' ')
        call add_default ('A_PRECT   ', 1, ' ')
        call add_default ('A_PRECSL  ', 1, ' ')
        call add_default ('A_PRECSC  ', 1, ' ')




  end subroutine export_state_init


  subroutine record_state_before(state, cam_in, tend, state_copy, tend_copy, relhum)
        use physconst,          only: gravit, rga, rair, cpair, latvap, rearth, pi, cappa
        use wv_saturation,      only: aqsat
        implicit none

     
      
        !
        ! Arguments
        !
        type(physics_state), intent(in) :: state
        type(cam_in_t),  intent(in) :: cam_in
        type(physics_tend ), intent(in) :: tend
        
        type(physics_state), intent(out) :: state_copy
        type(physics_tend ), intent(out) :: tend_copy
        real(r8),            intent(out) :: relhum(pcols,pver) ! temporary workspace

        real(r8) z3(pcols,pver)   ! geo-potential height
        real(r8) ftem(pcols,pver) ! temporary workspace
        real(r8) tem2(pcols,pver) ! temporary workspace
        integer :: lchnk, ncol, k

        lchnk = state%lchnk
        ncol  = state%ncol


        call aqsat (state%t    ,state%pmid  ,tem2    ,ftem    ,pcols   , &
                ncol ,pver  ,1       ,pver    )
        relhum(:ncol,:) = state%q(:ncol,:,1)/ftem(:ncol,:)*100._r8
         

        if( .not. is_output_state_vars_from_ctrl_run() ) then
                return
        end if

         

        call outfld('B_T     ',state%t , pcols   ,lchnk   )
        call outfld('B_U     ',state%u , pcols   ,lchnk   )
        call outfld('B_V     ',state%v , pcols   ,lchnk   )
        call outfld('B_Q     ',state%q(:ncol,:pver,1), pcols ,lchnk )
        call outfld('B_CLDLIQ',state%q(:ncol,:pver,ixcldliq), pcols ,lchnk )
        call outfld('B_CLDICE',state%q(:ncol,:pver,ixcldice), pcols ,lchnk )
        !call outfld('B_NUMLIQ',state%q(:ncol,:pver,ixnumliq), pcols ,lchnk )
        !call outfld('B_NUMICE',state%q(:ncol,:pver,ixnumice), pcols ,lchnk )
        call outfld('B_OMEGA ',state%omega, pcols, lchnk  )

       call outfld ('B_RELHUM',relhum     , pcols, lchnk  )

        !
        ! Add height of surface to midpoint height above surface
        !
        do k = 1, pver
                z3(:ncol,k) = state%zm(:ncol,k) + state%phis(:ncol)*rga
        end do
        call outfld('B_Z3    ',z3,pcols,lchnk)
 
        call outfld('B_PS    ',state%ps, pcols   ,lchnk   )
        call outfld('B_TS    ',cam_in%ts,        pcols, lchnk)
        call outfld('B_SHFLX ',cam_in%shf,       pcols, lchnk)
        call outfld('B_LHFLX ',cam_in%lhf,       pcols, lchnk)



        !
        ! Record the state at this point in time
        !
        
        call physics_state_copy( state, state_copy )
        call physics_tend_init( tend_copy )
        call physics_tend_copy( tend, tend_copy )


  end subroutine record_state_before


  subroutine record_state_after(state, tend, ztodt, initial_state, initial_tend)
        use physconst,          only: gravit, rga, rair, cpair, latvap, rearth, pi, cappa
        use constituents,  only: pcnst, cnst_name, cnst_longname, cnst_cam_outfld, ptendnam, cnst_get_ind
        use dycore,          only: dycore_is
        use check_energy,    only: check_energy_get_integrals

        implicit none

     
      
        !
        ! Arguments
        !
        type(physics_state), intent(in) :: state
        type(physics_tend ), intent(in) :: tend
        real(r8)           , intent(in) :: ztodt            ! physics timestep
        type(physics_state), intent(in) :: initial_state
        type(physics_tend ), intent(in) :: initial_tend

        real(r8) :: ftem2(pcols     ) ! Temporary workspace for outfld variables
        real(r8) :: ftem3(pcols,pver) ! Temporary workspace for outfld variables
        real(r8) :: rtdt
        real(r8) :: heat_glob         ! global energy integral (FV only)
        integer  :: lchnk, ncol, k

        if( .not. is_output_state_vars_from_ctrl_run() ) then
                return
        end if

         
        lchnk = state%lchnk
        ncol  = state%ncol
        rtdt  = 1._r8/ztodt

        !
        ! T-tendency due to FV Energy fixer (remove from total physics tendency diagnostic)
        !
        !if (dycore_is('LR')) then
        !    call check_energy_get_integrals( heat_glob_out=heat_glob )
        !    ftem2(:ncol)  = heat_glob/cpair
        !    call outfld('A_TFIX', ftem2, pcols, lchnk   )
        !    ftem3(:ncol,:pver)  = tend%dtdt(:ncol,:pver) - heat_glob/cpair
        !else
        !    ftem3(:ncol,:pver)  = tend%dtdt(:ncol,:pver)
        !end if

        !
        ! Total physics tendency for Temperature
        !
        ftem3(:ncol,:pver)  = tend%dtdt(:ncol,:pver) - initial_tend%dtdt(:ncol,:pver) 
        call outfld('A_PTTEND',ftem3, pcols, lchnk )


        !
        ! Total physics tendency for moisture and other tracers
        !
        ftem3(:ncol,:pver) = (state%q(:ncol,:pver, 1) - initial_state%q(:ncol,:pver, 1) )*rtdt
        call outfld ('A_'//ptendnam(       1), ftem3, pcols, lchnk)

        ftem3(:ncol,:pver) = (state%q(:ncol,:pver, ixcldliq) - initial_state%q(:ncol,:pver, ixcldliq) )*rtdt
        call outfld ('A_PTECLDLIQ', ftem3, pcols, lchnk)

        ftem3(:ncol,:pver) = (state%q(:ncol,:pver, ixcldice) - initial_state%q(:ncol,:pver, ixcldice) )*rtdt
        call outfld ('A_PTECLDICE', ftem3, pcols, lchnk)

        !TODO: for cam5 and above

        !ftem3(:ncol,:pver) = (state%q(:ncol,:pver, ixnumliq) - initial_state%q(:ncol,:pver, ixnumliq) )*rtdt
        !call outfld ('A_PTENUMLIQ', ftem3, pcols, lchnk)

        !ftem3(:ncol,:pver) = (state%q(:ncol,:pver, ixnumice) - initial_state%q(:ncol,:pver, ixnumice) )*rtdt
        !call outfld ('A_PTENUMICE', ftem3, pcols, lchnk)

  end subroutine record_state_after


  subroutine record_state_after_additional(state, prect, precc, precl, precsc, precsl)
        use physconst,          only: gravit, rga, rair, cpair, latvap, rearth, pi, cappa
        use constituents,  only: pcnst, cnst_name, cnst_longname, cnst_cam_outfld, ptendnam, cnst_get_ind
        use dycore,          only: dycore_is
        use check_energy,    only: check_energy_get_integrals

        implicit none

     
      
        !
        ! Arguments
        !
        type(physics_state), intent(in) :: state
        real(r8),            intent(in) :: prect(pcols) 
        real(r8),            intent(in) :: precc(pcols) 
        real(r8),            intent(in) :: precl(pcols) 
        real(r8),            intent(in) :: precsc(pcols) 
        real(r8),            intent(in) :: precsl(pcols) 

        !real(r8) :: rtdt
        integer  :: lchnk, ncol, k

        if( .not. is_output_state_vars_from_ctrl_run() ) then
                return
        end if

         
        lchnk = state%lchnk
        ncol  = state%ncol
        !rtdt  = 1._r8/ztodt


        call outfld('A_PRECC   ', precc , pcols, lchnk )
        call outfld('A_PRECL   ', precl , pcols, lchnk )
        call outfld('A_PRECSC  ', precsc, pcols, lchnk )
        call outfld('A_PRECSL  ', precsl, pcols, lchnk )
        call outfld('A_PRECT   ', prect , pcols, lchnk )



  end subroutine record_state_after_additional



  subroutine physics_tend_copy(tend, tend_copy)

      implicit none

      !
      ! Arguments
      !
      type(physics_tend), intent(in) :: tend
      type(physics_tend), intent(out) :: tend_copy

      !
      ! Local variables
      !

      tend_copy%dtdt(:pcols,:pver) = tend%dtdt(:pcols,:pver)
      tend_copy%dudt(:pcols,:pver) = tend%dudt(:pcols,:pver)
      tend_copy%dvdt(:pcols,:pver) = tend%dvdt(:pcols,:pver)
      tend_copy%flx_net(:pcols)    = tend%flx_net(:pcols)
      tend_copy%te_tnd(:pcols)     = tend%te_tnd(:pcols)
      tend_copy%tw_tnd(:pcols)     = tend%tw_tnd(:pcols)

  end subroutine physics_tend_copy



end module export_state

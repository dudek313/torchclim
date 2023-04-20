module machine_learning_model_config
   use shr_kind_mod,      only : r8=>shr_kind_r8
   use time_manager,    only: is_first_step, get_nstep
   use physics_types,     only: physics_state

   !----------------------------------------------- !
   ! Purpose:                                       !
   !                                                !
   ! CAM interface to the shallow convection scheme !
   !                                                !
   ! Author: D.B. Coleman                           !
   !         Sungsu Park. Jan. 2010.                !
   !                                                !
   !----------------------------------------------- !

   implicit none
   private                 
   save

   public :: &
             is_ml_model_enabled,    &
             is_ml_model_standalone, &
             is_radiation_from_ml,   & 
             is_bypass_tphysac,      & 
             is_bypass_surface_drag, &   
             is_output_state_vars_from_ctrl_run, & 
             first_ml_time_step


  !logical, parameter :: ml_model_skip_first_steps = .false.
  logical, parameter :: ml_model_skip_first_steps = .true.
  
  logical, parameter :: ml_model_enabled = .false.
  !logical, parameter :: ml_model_enabled = .true.
 
  !logical, parameter :: ml_model_standalone = .false.
  logical, parameter :: ml_model_standalone = .true.


  logical, parameter :: radiation_from_ml = .true.

  !logical, parameter :: output_state_vars_from_ctrl_run = .false.
  logical, parameter :: output_state_vars_from_ctrl_run = .true.


  logical, parameter :: bypass_tphysac = .false.
  !logical, parameter :: bypass_tphysac = .true.
  logical, parameter :: bypass_surface_drag = .false.

  logical, parameter :: ml_only_in_lat_band = .false.

  real(r8), parameter :: pi = 3.14159265359_r8
  real(r8), parameter :: rad_to_deg = 180._r8/pi
  real(r8), parameter :: max_lat =  30._r8
  real(r8), parameter :: min_lat = -30._r8

  integer, parameter :: steps_to_skip = ((1*24*2) + 1) !days*hours*30min time step

  contains

  logical function is_ml_model_enabled(state)
        type(physics_state), intent(in) :: state

        if(.not. ml_model_enabled) then
            is_ml_model_enabled = .false.
            return
        end if

        !if (ml_model_skip_first_step .and. is_first_step()) then
        if (ml_model_skip_first_steps .and. is_skip_time_step()) then
            is_ml_model_enabled = .false.
            return
        end if

        if( .not. ml_only_in_lat_band) then 
            is_ml_model_enabled = .true.
            return 
        end if

        is_ml_model_enabled = is_in_lat_band(state)

  end function

  logical function is_ml_model_standalone()
        is_ml_model_standalone = ml_model_standalone !.and. is_ml_model_enabled(state)
  end function

  logical function is_radiation_from_ml()
        is_radiation_from_ml = radiation_from_ml !.and. is_ml_model_enabled(state)
  end function

  logical function is_bypass_tphysac()
        is_bypass_tphysac = bypass_tphysac !.and. is_ml_model_enabled(state)
  end function

  logical function is_bypass_surface_drag()
        is_bypass_surface_drag = bypass_surface_drag !.and. is_ml_model_enabled(state)
  end function

  logical function is_output_state_vars_from_ctrl_run()
         is_output_state_vars_from_ctrl_run = output_state_vars_from_ctrl_run
  end function

  integer function first_ml_time_step()
          first_ml_time_step = steps_to_skip
  end function

  !=============================================================================== !
  !                                                                                !
  ! PRIVATE                                                                                !
  !                                                                                !
  !=============================================================================== !


  logical function is_in_lat_band(state)
        type(physics_state), intent(in) :: state
        real(r8) lat

        !assume mid lat of the tile 
        lat = (minval(state%lat) + maxval(state%lat))/2._r8
        lat = lat * rad_to_deg

        if( lat < min_lat ) then
            is_in_lat_band = .false.
            return
        end if 

        if( lat > max_lat ) then
            is_in_lat_band = .false.
            return
        end if 

        is_in_lat_band = .true.
  end function

  logical function is_skip_time_step()
          !Designed to count steps after which the ML model runs
          !we do so in order to smooth noise of initial steps
          integer :: nstep
          nstep = get_nstep()

          is_skip_time_step = nstep < first_ml_time_step() 
  end function


  end module machine_learning_model_config

module ml_solin
   !
   ! This module distills the code to calculate solar insulation
   ! for more see radiation.F90 and radsw.F90
   !

   use shr_kind_mod,    only: r8 => shr_kind_r8
   use ppgrid

   implicit none
   private
   save

   public :: &
      ml_calc_solin 


   contains


   subroutine ml_calc_solin(state, solin)

      use time_manager,    only: get_curr_calday
      use phys_grid,       only: get_rlat_all_p, get_rlon_all_p
      use physics_types,   only: physics_state, physics_ptend
      !use camsrfexch_types,only: cam_out_t, cam_in_t

      type(physics_state), intent(in), target :: state
      !type(cam_in_t),      intent(in)         :: cam_in

      real(r8), intent(out) :: solin(pcols)     ! Incident solar flux
     
      integer lchnk, ncol
      real(r8) :: calday                        ! current calendar day
      real(r8) :: clat(pcols)                   ! current latitudes(radians)
      real(r8) :: clon(pcols)                   ! current longitudes(radians)
      real(r8) coszrs(pcols)                     ! Cosine solar zenith angle

      real(r8) eccf                 ! Earth/sun distance factor

      real(r8), parameter :: cgs2mks = 1.e-3_r8

      ! Gathered indicies of day and night columns
      !  chunk_column_index = IdxDay(daylight_column_index)
      integer :: Nday                      ! Number of daylight columns
      integer :: Nnite                     ! Number of night columns
      integer, dimension(pcols) :: IdxDay  ! Indicies of daylight coumns
      integer, dimension(pcols) :: IdxNite ! Indicies of night coumns

      integer  i,k,m 

      lchnk = state%lchnk
      ncol  = state%ncol

      calday = get_curr_calday()

      ! Cosine solar zenith angle for current time step
      call get_rlat_all_p(lchnk, ncol, clat)
      call get_rlon_all_p(lchnk, ncol, clon)
      call zenith (calday, clat, clon, coszrs, ncol)
 
      !call output_rad_data( pbuf, state, cam_in, landm, coszrs )

      ! Gather night/day column indices.
      Nday = 0
      Nnite = 0
      do i = 1, ncol
         if ( coszrs(i) > 0.0_r8 ) then
            Nday = Nday + 1
            IdxDay(Nday) = i
         else
            Nnite = Nnite + 1
            IdxNite(Nnite) = i
         end if
      end do

      call get_eccf(eccf)
 
      call radcswmx_ml(lchnk   ,ncol, &
         eccf    ,coszrs ,solin, &
         Nday    ,Nnite   ,IdxDay  ,IdxNite )
    
      !
      ! Convert units of shortwave fields needed by rest of model from CGS to MKS
      !
      do i=1,ncol
        solin(i) = solin(i)*cgs2mks
      end do

      !call outfld('SOLIN   ',solin ,pcols,lchnk)

   end subroutine !ml_calc_solin



   subroutine radcswmx_ml(lchnk   ,ncol, &
      eccf    ,E_coszrs ,solin, &
      Nday    ,Nnite   ,IdxDay  ,IdxNite )

      !use rad_solar_var,    only: get_variability
      use cmparray_mod,     only: CmpDayNite, ExpDayNite
      use solar_data,       only: sol_tsi, do_spctrl_scaling, ref_tsi
     

      !
      !
      !
      !
      !


 !
 ! Input arguments
 !
      integer, intent(in) :: lchnk             ! chunk identifier
      integer, intent(in) :: ncol              ! number of atmospheric columns
 
      integer,intent(in) :: Nday                      ! Number of daylight columns
      integer,intent(in) :: Nnite                     ! Number of night columns
      integer,intent(in), dimension(pcols) :: IdxDay  ! Indicies of daylight coumns
      integer,intent(in), dimension(pcols) :: IdxNite ! Indicies of night coumns

      real(r8), intent(in) :: eccf             ! Eccentricity factor (1./earth-sun dist^2)
      real(r8), intent(in) :: E_coszrs(pcols)    ! Cosine solar zenith angle

      real(r8), intent(out) :: solin(pcols)     ! Incident solar flux

      real(r8) :: coszrs(pcols)    ! Cosine solar zenith angle
      real(r8) :: tot_irrad

      integer i

      solin(1:ncol)    = 0.0_r8

 !
 ! If night everywhere, return:
 !
      if ( Nday == 0 ) then
         return
      endif

 !
 ! Rearrange input arrays
 !
      call CmpDayNite(E_coszrs, coszrs,    Nday, IdxDay, Nnite, IdxNite, 1, pcols)

      if ( do_spctrl_scaling ) then
         !call get_variability(sfac)
         tot_irrad = ref_tsi
      else
         tot_irrad = sol_tsi
      endif


 !
 ! Define solar incident radiation and interface pressures:
 !

 do i=1,Nday
      solin(i)  = tot_irrad*1.e3_r8*eccf*coszrs(i)
 end do
 !
 ! Rearrange output arrays
 !

      call ExpDayNite(solin,       Nday, IdxDay, Nnite, IdxNite, 1, pcols)

 
   end subroutine radcswmx_ml


   subroutine get_eccf(eccf)
      use shr_orb_mod
      use time_manager, only: get_curr_calday
      use cam_control_mod, only: lambm0, obliqr, mvelpp, eccen

      real(r8), intent(out) :: eccf                ! Earth-sun distance factor

      real(r8) :: calday       ! current calendar day
      real(r8) :: delta        ! Solar declination angle


      calday = get_curr_calday()
      call shr_orb_decl (calday  ,eccen     ,mvelpp  ,lambm0  ,obliqr  , &
         delta   ,eccf)
   end subroutine get_eccf

end module !ml_solin

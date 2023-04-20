subroutine ml_srfxfer(state,cam_out,prect, precc, precl, &
                   precsc, precsl)


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
   real(r8), intent(in) :: precc(pcols)                ! snow from ZM   convection
   real(r8), intent(in) :: precl(pcols)                ! total precipitation   from Hack convection
   real(r8), intent(in) :: precsc(pcols)                ! snow from   Hack   convection
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

end subroutine ml_srfxfer








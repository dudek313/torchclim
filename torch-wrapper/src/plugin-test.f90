program test_plugin
  use torch_plugin, only: predict_model, predict_model_v2
  use cam4_profile  

  implicit none


  real(r8) :: prect
  real(r8) :: precc
  real(r8) :: precl
  real(r8) :: precsc
  real(r8) :: precsl


  real(r8) :: flds
  real(r8) :: srfrad
  real(r8) :: psl
  real(r8) :: soll
  real(r8) :: sols
  real(r8) :: solld
  real(r8) :: solsd

  real(r8), dimension(26) :: pteq      ! total physics moist tendencies
  real(r8), dimension(26) :: pttend    ! total physics heating tendencies
  real(r8), dimension(26) :: dcq       ! total moist tendencies from moist processes
  real(r8), dimension(26) :: dtcond    ! total heating tendencies from moist processes
  real(r8), dimension(26) :: qrs       ! total heating tendencies due to SW radiation
  real(r8), dimension(26) :: qrl       ! total heating tendencies due to LW radiation
  real(r8), dimension(26) :: cld       ! total cloud fraction
  real(r8), dimension(26) :: concld    ! total convective cloud fraction


  integer :: i


  do i = 1,1

    if( .false. ) then

        call predict_model( &
               state_qv,&
               state_t, &
               state_u, &
               state_v, &
               state_omega, &
               state_z3, &
               state_ps, &
                       !
               solin, & ! solin from direct calculation
               cam_in_shf, & ! sensible heat flux from cam_in
               cam_in_lhf, & ! latend heat flux from cam_in

               !inout
               fsns,   & ! from tphysbc
               flns,   & ! from tphysbc
               fsnt,   & ! from tphysbc
               flnt,   & ! from tphysbc
               fsds,   & ! from tphysbc

               !out
               flds,   & !
               srfrad, & ! from tphysbc
               soll,   & !
               sols,   & !
               solld,  & !
               solsd,  & !
               psl  ,  & !
                      !
               pteq,   & ! total physics moist tendencies
               pttend, & ! total physics heating tendencies
               dcq,    & ! total moist tendencies from moist processes
               dtcond, & ! total heating tendencies from moist processes
               qrs,    & ! total heating tendencies due to SW radiation
               qrl,    & ! total heating tendencies due to LW radiation
               cld,    & ! total cloud fraction
               concld, & ! total convective cloud fraction

               prect, & ! prect - not used aside from output
               precc, &
               precl, &
               precsc, &
               precsl, &
              .true.)
    else
        call predict_model_v2( &
               state_qv,&
               state_t, &
               state_u, &
               state_v, &
               state_omega, &
               state_z3, &
               state_ps, &
               297._r8, & !TODO: should be TS!
                       !
               solin, & ! solin from direct calculation
               cam_in_shf, & ! sensible heat flux from cam_in
               cam_in_lhf, & ! latend heat flux from cam_in

               landfrac, & 
               ocnfrac, &
               icefrac, &

               !inout
               fsns,   & ! from tphysbc
               flns,   & ! from tphysbc
               fsnt,   & ! from tphysbc
               flnt,   & ! from tphysbc
               fsds,   & ! from tphysbc

               !out
               flds,   & !
               srfrad, & ! from tphysbc
               soll,   & !
               sols,   & !
               solld,  & !
               solsd,  & !
                      !
               pteq,   & ! total physics moist tendencies
               pttend, & ! total physics heating tendencies

               prect, & ! prect - not used aside from output
               precc, &
               precl, &
               precsc, &
               precsl, &
              .false.)

    end if

    print *,'test...'
  end do



  if(.true.) then
     write(*,*) 'test app printing results:'
     write(*,*) '=========================='

     !write(*,*) 'ml_model_tend, pttend  : ', pttend
     !write(*,*) 'ml_model_tend, pttend  : ', pteq
     !write(*,*) 'ml_model_tend, fsns  : ', fsns
     !write(*,*) 'ml_model_tend, flns  : ', flns
     !write(*,*) 'ml_model_tend, fsnt  : ', fsnt
     !write(*,*) 'ml_model_tend, flnt  : ', flnt
     !write(*,*) 'ml_model_tend, fsds  : ', fsds
     !write(*,*) 'ml_model_tend, flds  : ', flds


     write(*,*) 'pttend = [ '
     write(*,'(*(G0.7:", "))') pttend
     write(*,*) '    ]'
     write(*,*) ''
     write(*,*) 'pteq = [ '
     write(*,'(*(G0.7:", "))') pteq
     write(*,*) '    ]'
     write(*,*) ''
     write(*,*) 'precc = ', precc,''
     write(*,*) 'precl = ', precl,''
     write(*,*) ''
     write(*,*) 'fsns = ', fsns,''
     write(*,*) 'flns = ', flns,''
     write(*,*) 'fsnt = ', fsnt,''
     write(*,*) 'flnt = ', flnt,''
     write(*,*) 'fsds = ', fsds,''
     write(*,*) 'flds = ', flds,''
     write(*,*) 'srfrad = ', srfrad,''
     write(*,*) ''
     write(*,*) 'soll  = ', soll,''
     write(*,*) 'solld = ', solld,''
     write(*,*) 'sols  = ', sols,''
     write(*,*) 'solsd = ', solsd,''
  end if


  print *, 'the end'

end program

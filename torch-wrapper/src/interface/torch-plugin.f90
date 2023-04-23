module torch_plugin
    !DF1 use iso_c_binding
    use, intrinsic :: iso_c_binding

    private
    public create_model
    public predict_model
    public predict_model_v2
    public predict_model_v3
    public predict_model_v4
    public predict_model_v5
    public predict_model_v6

    ! include the c interface definitions for the private impl. 
    include "torch-wrap-cdef.f90"

    ! mimic the behavior of CESM
    integer,parameter :: R8 = selected_real_kind(12) ! 8 byte real -> SHR_KIND_R8
    integer,parameter :: R4 = selected_real_kind( 6) ! 4 byte real -> SHR_KIND_R4

contains

    function create_model(path)
        implicit none
        character(len=*), intent(in) :: path
        character(len=1, kind=C_CHAR) :: c_path(len_trim(path) + 1)
        type(c_ptr) :: create_model
        integer :: path_len, i

        ! copy fortran string to c string
        path_len = len_trim(path)
        do i = 1,path_len
            c_path(i) = path(i:i)
        end do

        c_path(path_len + 1) = C_NULL_CHAR

        ! create libtorch script model allocated on heap
        create_model = model_create_c(c_path)
    end function

    subroutine predict_model( &
            !inputs on vertical
            qv, t, u, v, omega, z3, &
            
            !inputs on 1d
            ps, solin, shflx, lhflx, & 
            
            !inout on 1d
            fsns, flns, fsnt, flnt, fsds, &
            
            ! 1d outputs related to radiation 
            flds, srfrad, soll, sols, solld, solsd, &
            psl, &

            !outputs on vertical
            pteq, pttend, &
            dqc, dtcond, &
            qrs, qrl, &
            cloud, concld, &

            ! outputs on 1d
            prect, precc, precl, &
            precsc, precsl, &

            debug_log &

            )
        implicit none
        REAL(r8), dimension(:), intent(in) :: qv
        REAL(r8), dimension(:), intent(in) :: t
        REAL(r8), dimension(:), intent(in) :: u
        REAL(r8), dimension(:), intent(in) :: v
        REAL(r8), dimension(:), intent(in) :: omega
        REAL(r8), dimension(:), intent(in) :: z3
        REAL(r8), intent(in) :: ps
        REAL(r8), intent(in) :: solin 
        REAL(r8), intent(in) :: shflx
        REAL(r8), intent(in) :: lhflx
        REAL(r8), intent(inout) :: fsns
        REAL(r8), intent(inout) :: flns
        REAL(r8), intent(inout) :: fsnt
        REAL(r8), intent(inout) :: flnt
        REAL(r8), intent(inout) :: fsds
        logical , intent(in) :: debug_log 

        REAL(r8), intent(out) :: flds
        REAL(r8), intent(out) :: srfrad
        REAL(r8), intent(out) :: soll
        REAL(r8), intent(out) :: sols
        REAL(r8), intent(out) :: solld
        REAL(r8), intent(out) :: solsd
        REAL(r8), intent(out) :: psl


        REAL(r8), dimension(:), intent(out) :: pteq
        REAL(r8), dimension(:), intent(out) :: pttend
        REAL(r8), dimension(:), intent(out) :: dqc
        REAL(r8), dimension(:), intent(out) :: dtcond
        REAL(r8), dimension(:), intent(out) :: qrs
        REAL(r8), dimension(:), intent(out) :: qrl
        REAL(r8), dimension(:), intent(out) :: cloud
        REAL(r8), dimension(:), intent(out) :: concld
        REAL(r8), intent(out) :: prect
        REAL(r8), intent(out) :: precc
        REAL(r8), intent(out) :: precl
        REAL(r8), intent(out) :: precsc
        REAL(r8), intent(out) :: precsl

        !note the implicit conversion from R8 to c_float - since the NN is based
        !on single prcision output from cesm cam, but is defied as capi type
        REAL(c_float), dimension((size(t) * 6) + 9) :: istacked
        REAL(c_float), dimension((size(t) * 8) + 5 + 7 + 5) :: ostacked
        integer :: sz = 0
        integer :: retval = 0
        integer :: output_size = 0

        if (debug_log) then
            write(*,*) "torch-plugin debug log, 3D inputs:"
            write(*,*) "----------------------------------"
            write(*,*) "temp : ", t
            write(*,*) "Qv   : ", qv
            write(*,*) "U    : ", u
            write(*,*) "V    : ", v
            write(*,*) "omega: ", omega
            write(*,*) "z3   : ", z3
            write(*,*) ""
            write(*,*) "torch-plugin debug log, 2D inputs:"
            write(*,*) "----------------------------------"
            write(*,*) "ps   : ", ps
            write(*,*) "solin: ", solin
            write(*,*) "shflx: ", shflx
            write(*,*) "lhflx: ", lhflx
            write(*,*) "fsns : ", fsns
            write(*,*) "flns : ", flns
            write(*,*) "fsnt : ", fsnt
            write(*,*) "flnt : ", flnt
            write(*,*) "fsds : ", fsds
            write(*,*) ""
         end if
        !write(*,*) "here 1"

        sz = size(t)
        output_size = size(ostacked)

        istacked((sz*0) + 1 : sz*1) = qv
        istacked((sz*1) + 1 : sz*2) = t
        istacked((sz*2) + 1 : sz*3) = u
        istacked((sz*3) + 1 : sz*4) = v
        istacked((sz*4) + 1 : sz*5) = omega
        istacked((sz*5) + 1 : sz*6) = z3

        istacked((sz*6) + 1) = ps
        istacked((sz*6) + 2) = solin
        istacked((sz*6) + 3) = shflx
        istacked((sz*6) + 4) = lhflx
        istacked((sz*6) + 5) = fsns
        istacked((sz*6) + 6) = flns
        istacked((sz*6) + 7) = fsnt
        istacked((sz*6) + 8) = flnt
        istacked((sz*6) + 9) = fsds


        call model_predict_c(                         &
                        input  = istacked,            &
                        input_size  = size(istacked), &
                        output = ostacked,            &
                        output_size = output_size,    &
                        retval = retval,              &
                        loopback = 0                  & !non-zero means test mode
                )


        pteq    = ostacked((sz*0) + 1 : sz*1)
        pttend  = ostacked((sz*1) + 1 : sz*2)
        dqc     = ostacked((sz*2) + 1 : sz*3)
        dtcond  = ostacked((sz*3) + 1 : sz*4)
        qrs     = ostacked((sz*4) + 1 : sz*5)
        qrl     = ostacked((sz*5) + 1 : sz*6)
        cloud   = ostacked((sz*6) + 1 : sz*7)
        concld  = ostacked((sz*7) + 1 : sz*8)

        fsns    = ostacked((sz*8) + 1)
        flns    = ostacked((sz*8) + 2)
        fsnt    = ostacked((sz*8) + 3)
        flnt    = ostacked((sz*8) + 4)
        fsds    = ostacked((sz*8) + 5)

        flds    = ostacked((sz*8) + 6)
        srfrad  = ostacked((sz*8) + 7)
        soll    = ostacked((sz*8) + 8)
        sols    = ostacked((sz*8) + 9)
        solld   = ostacked((sz*8) + 10)
        solsd   = ostacked((sz*8) + 11)
        psl     = ostacked((sz*8) + 12)

        prect   = ostacked((sz*8) + 13)
        precc   = ostacked((sz*8) + 14)
        precl   = ostacked((sz*8) + 15)
        precsc  = ostacked((sz*8) + 16)
        precsl  = ostacked((sz*8) + 17)


        if (debug_log) then
            !if (.true.) then
            write(*,*) "torch-plugin debug log, 3D outputs:"
            write(*,*) "-----------------------------------"
            write(*,*) "" 
            write(*,*) "ostacked : ", ostacked 
            write(*,*) "pttend: ", pttend
            write(*,*) "pteq  : ", pteq
            write(*,*) "" 
            write(*,*) "torch-plugin debug log, 2D outputs:"
            write(*,*) "-----------------------------------"
            write(*,*) "psl : ", psl
            write(*,*) "prect : ", prect
            write(*,*) "precc : ", precc
            write(*,*) "fsds  : ", fsds
            write(*,*) "flds  : ", flds
            write(*,*) "srfrad: ", srfrad
            write(*,*) "" 
            !write(*,*) "ostacked : ", ostacked 
         end if

    end subroutine

    subroutine predict_model_v2( &
            !inputs on vertical
            qv, t, u, v, omega, z3, &
            
            !inputs on 1d
            ps, ts, solin, shflx, lhflx, & 

            !added 1d
            landfrac, ocnfrac, icefrac, &

            !inout on 1d
            fsns, flns, fsnt, flnt, fsds, &
            
            ! 1d outputs related to radiation 
            flds, srfrad, soll, sols, solld, solsd, &

            !outputs on vertical
            pteq, pttend, &            !+

            ! outputs on 1d
            prect, precc, precl, &
            precsc, precsl, &

            debug_log &

            )
        implicit none
        REAL(r8), dimension(:), intent(in) :: qv
        REAL(r8), dimension(:), intent(in) :: t
        REAL(r8), dimension(:), intent(in) :: u
        REAL(r8), dimension(:), intent(in) :: v
        REAL(r8), dimension(:), intent(in) :: omega
        REAL(r8), dimension(:), intent(in) :: z3
        REAL(r8), intent(in) :: ps
        REAL(r8), intent(in) :: ts
        REAL(r8), intent(in) :: solin 
        REAL(r8), intent(in) :: shflx
        REAL(r8), intent(in) :: lhflx

        REAL(r8), intent(in) :: landfrac 
        REAL(r8), intent(in) :: ocnfrac
        REAL(r8), intent(in) :: icefrac

        REAL(r8), intent(inout) :: fsns
        REAL(r8), intent(inout) :: flns
        REAL(r8), intent(inout) :: fsnt
        REAL(r8), intent(inout) :: flnt
        REAL(r8), intent(inout) :: fsds
        logical , intent(in) :: debug_log 

        REAL(r8), intent(out) :: flds
        REAL(r8), intent(out) :: srfrad
        REAL(r8), intent(out) :: soll
        REAL(r8), intent(out) :: sols
        REAL(r8), intent(out) :: solld
        REAL(r8), intent(out) :: solsd


        REAL(r8), dimension(:), intent(out) :: pteq
        REAL(r8), dimension(:), intent(out) :: pttend
        REAL(r8), intent(out) :: prect
        REAL(r8), intent(out) :: precc
        REAL(r8), intent(out) :: precl
        REAL(r8), intent(out) :: precsc
        REAL(r8), intent(out) :: precsl

        !note the implicit conversion from R8 to c_float - since the NN is based
        !on single prcision output from cesm cam, but is defied as capi type
        !REAL(c_float), dimension((size(t) * 6) + 13) :: istacked
        REAL(c_float), dimension((size(t) * 6) + 12) :: istacked
        REAL(c_float), dimension((size(t) * 2) + 5 + 6 + 5) :: ostacked
        integer :: sz = 0
        integer :: retval = 0
        integer :: output_size = 0

        if (debug_log) then
            write(*,*) "torch-plugin debug log, 3D inputs:"
            write(*,*) "----------------------------------"
            write(*,*) "temp : ", t
            write(*,*) "Qv   : ", qv
            write(*,*) "U    : ", u
            write(*,*) "V    : ", v
            write(*,*) "omega: ", omega
            write(*,*) "z3   : ", z3
            write(*,*) ""
            write(*,*) "torch-plugin debug log, 2D inputs:"
            write(*,*) "----------------------------------"
            write(*,*) "ps   : ", ps
            write(*,*) "ts   : ", ts
            write(*,*) "solin: ", solin
            write(*,*) "shflx: ", shflx
            write(*,*) "lhflx: ", lhflx
            write(*,*) "fsns : ", fsns
            write(*,*) "flns : ", flns
            write(*,*) "fsnt : ", fsnt
            write(*,*) "flnt : ", flnt
            write(*,*) "fsds : ", fsds
            write(*,*) ""
         end if
        !write(*,*) "here 1"

        sz = size(t)
        output_size = size(ostacked)

        istacked((sz*0) + 1 : sz*1) = qv
        istacked((sz*1) + 1 : sz*2) = t
        istacked((sz*2) + 1 : sz*3) = u
        istacked((sz*3) + 1 : sz*4) = v
        istacked((sz*4) + 1 : sz*5) = omega
        istacked((sz*5) + 1 : sz*6) = z3

        istacked((sz*6) + 1) = ps
        istacked((sz*6) + 2) = solin
        istacked((sz*6) + 3) = shflx
        istacked((sz*6) + 4) = lhflx
        istacked((sz*6) + 5) = landfrac
        istacked((sz*6) + 6) = ocnfrac
        istacked((sz*6) + 7) = icefrac
        istacked((sz*6) + 8) = fsns
        istacked((sz*6) + 9) = flns
        istacked((sz*6) + 10) = fsnt
        istacked((sz*6) + 11) = flnt
        istacked((sz*6) + 12) = fsds
        !istacked((sz*6) + 13) = ts


        call model_predict_c(                         &
                        input  = istacked,            &
                        input_size  = size(istacked), &
                        output = ostacked,            &
                        output_size = output_size,    &
                        retval = retval,              &
                        loopback = 0                  & !non-zero means test mode
                )


        pttend  = ostacked((sz*0) + 1 : sz*1)
        pteq    = ostacked((sz*1) + 1 : sz*2)

        fsns    = ostacked((sz*2) + 1)
        flns    = ostacked((sz*2) + 2)
        fsnt    = ostacked((sz*2) + 3)
        flnt    = ostacked((sz*2) + 4)
        fsds    = ostacked((sz*2) + 5)

        flds    = ostacked((sz*2) + 6)
        srfrad  = ostacked((sz*2) + 7)
        soll    = ostacked((sz*2) + 8)
        sols    = ostacked((sz*2) + 9)
        solld   = ostacked((sz*2) + 10)
        solsd   = ostacked((sz*2) + 11)

        prect   = ostacked((sz*2) + 12)
        precc   = ostacked((sz*2) + 13)
        precl   = ostacked((sz*2) + 14)
        precsc  = ostacked((sz*2) + 15)
        precsl  = ostacked((sz*2) + 16)


        if (debug_log) then
            !if (.true.) then
            write(*,*) "torch-plugin debug log, 3D outputs:"
            write(*,*) "-----------------------------------"
            write(*,*) "" 
            write(*,*) "ostacked : ", ostacked 
            write(*,*) "pttend: ", pttend
            write(*,*) "pteq  : ", pteq
            write(*,*) "" 
            write(*,*) "torch-plugin debug log, 2D outputs:"
            write(*,*) "-----------------------------------"
            write(*,*) "prect : ", prect
            write(*,*) "precc : ", precc
            write(*,*) "fsds  : ", fsds
            write(*,*) "flds  : ", flds
            write(*,*) "srfrad: ", srfrad
            write(*,*) "" 
            !write(*,*) "ostacked : ", ostacked 
         end if
    end subroutine

    subroutine predict_model_v3( &
            !inputs on vertical
            qv, t, u, v, omega, z3, &
            
            !inputs on 1d
            ps, ts, solin, shflx, lhflx, & 

            !added 1d
            landfrac, ocnfrac, icefrac, &

            !out on 1d
            fsns, flns, fsnt, flnt, fsds, &
            
            ! 1d outputs related to radiation 
            flds, srfrad, soll, sols, solld, solsd, &

            !outputs on vertical
            pteq, pttend, &            !+

            ! outputs on 1d
            prect, precc, precl, &
            precsc, precsl, &

            debug_log &

            )
        implicit none
        REAL(r8), dimension(:), intent(in) :: qv
        REAL(r8), dimension(:), intent(in) :: t
        REAL(r8), dimension(:), intent(in) :: u
        REAL(r8), dimension(:), intent(in) :: v
        REAL(r8), dimension(:), intent(in) :: omega
        REAL(r8), dimension(:), intent(in) :: z3
        REAL(r8), intent(in) :: ps
        REAL(r8), intent(in) :: ts
        REAL(r8), intent(in) :: solin 
        REAL(r8), intent(in) :: shflx
        REAL(r8), intent(in) :: lhflx

        REAL(r8), intent(in) :: landfrac 
        REAL(r8), intent(in) :: ocnfrac
        REAL(r8), intent(in) :: icefrac

        REAL(r8), intent(out) :: fsns
        REAL(r8), intent(out) :: flns
        REAL(r8), intent(out) :: fsnt
        REAL(r8), intent(out) :: flnt
        REAL(r8), intent(out) :: fsds
        logical , intent(in) :: debug_log 

        REAL(r8), intent(out) :: flds
        REAL(r8), intent(out) :: srfrad
        REAL(r8), intent(out) :: soll
        REAL(r8), intent(out) :: sols
        REAL(r8), intent(out) :: solld
        REAL(r8), intent(out) :: solsd


        REAL(r8), dimension(:), intent(out) :: pteq
        REAL(r8), dimension(:), intent(out) :: pttend
        REAL(r8), intent(out) :: prect
        REAL(r8), intent(out) :: precc
        REAL(r8), intent(out) :: precl
        REAL(r8), intent(out) :: precsc
        REAL(r8), intent(out) :: precsl

        !note the implicit conversion from R8 to c_float - since the NN is based
        !on single prcision output from cesm cam, but is defied as capi type
        REAL(c_float), dimension((size(t) * 6) + 8) :: istacked
        REAL(c_float), dimension((size(t) * 2) + 5 + 6 + 5) :: ostacked
        integer :: sz = 0
        integer :: retval = 0
        integer :: output_size = 0

        if (debug_log) then
            write(*,*) "torch-plugin debug log, 3D inputs:"
            write(*,*) "----------------------------------"
            write(*,*) "temp : ", t
            write(*,*) "Qv   : ", qv
            write(*,*) "U    : ", u
            write(*,*) "V    : ", v
            write(*,*) "omega: ", omega
            write(*,*) "z3   : ", z3
            write(*,*) ""
            write(*,*) "torch-plugin debug log, 2D inputs:"
            write(*,*) "----------------------------------"
            write(*,*) "ps   : ", ps
            write(*,*) "ts   : ", ts
            write(*,*) "solin: ", solin
            write(*,*) "shflx: ", shflx
            write(*,*) "lhflx: ", lhflx
            write(*,*) ""
         end if
        !write(*,*) "here 1"

        sz = size(t)
        output_size = size(ostacked)

        istacked((sz*0) + 1 : sz*1) = qv
        istacked((sz*1) + 1 : sz*2) = t
        istacked((sz*2) + 1 : sz*3) = u
        istacked((sz*3) + 1 : sz*4) = v
        istacked((sz*4) + 1 : sz*5) = omega
        istacked((sz*5) + 1 : sz*6) = z3

        istacked((sz*6) + 1) = ps
        istacked((sz*6) + 2) = solin
        istacked((sz*6) + 3) = shflx
        istacked((sz*6) + 4) = lhflx
        istacked((sz*6) + 5) = landfrac
        istacked((sz*6) + 6) = ocnfrac
        istacked((sz*6) + 7) = icefrac
        istacked((sz*6) + 8) = ts


        call model_predict_c(                         &
                        input  = istacked,            &
                        input_size  = size(istacked), &
                        output = ostacked,            &
                        output_size = output_size,    &
                        retval = retval,              &
                        loopback = 0                  & !non-zero means test mode
                )


        pttend  = ostacked((sz*0) + 1 : sz*1)
        pteq    = ostacked((sz*1) + 1 : sz*2)

        fsns    = ostacked((sz*2) + 1)
        flns    = ostacked((sz*2) + 2)
        fsnt    = ostacked((sz*2) + 3)
        flnt    = ostacked((sz*2) + 4)
        fsds    = ostacked((sz*2) + 5)

        flds    = ostacked((sz*2) + 6)
        srfrad  = ostacked((sz*2) + 7)
        soll    = ostacked((sz*2) + 8)
        sols    = ostacked((sz*2) + 9)
        solld   = ostacked((sz*2) + 10)
        solsd   = ostacked((sz*2) + 11)

        prect   = ostacked((sz*2) + 12)
        precc   = ostacked((sz*2) + 13)
        precl   = ostacked((sz*2) + 14)
        precsc  = ostacked((sz*2) + 15)
        precsl  = ostacked((sz*2) + 16)


        if (debug_log) then
            !if (.true.) then
            write(*,*) "torch-plugin debug log, 3D outputs:"
            write(*,*) "-----------------------------------"
            write(*,*) "" 
            write(*,*) "ostacked : ", ostacked 
            write(*,*) "pttend: ", pttend
            write(*,*) "pteq  : ", pteq
            write(*,*) "" 
            write(*,*) "torch-plugin debug log, 2D outputs:"
            write(*,*) "-----------------------------------"
            write(*,*) "prect : ", prect
            write(*,*) "precc : ", precc
            write(*,*) "fsds  : ", fsds
            write(*,*) "flds  : ", flds
            write(*,*) "srfrad: ", srfrad
            write(*,*) "" 
            !write(*,*) "ostacked : ", ostacked 
         end if
    end subroutine

    subroutine predict_model_v4( &
            !inputs on vertical
            qv, t, u, v, omega, &
            
            !inputs on 1d
            ps, ts, solin, shflx, lhflx, & 

            !added 1d
            landfrac, ocnfrac, icefrac, &

            !out on 1d
            fsns, flns, fsnt, flnt, fsds, &
            
            ! 1d outputs related to radiation 
            flds, srfrad, soll, sols, solld, solsd, &

            !outputs on vertical
            pteq, pttend, &            !+

            ! outputs on 1d
            prect, precc, precl, &
            precsc, precsl, &

            debug_log &

            )
        implicit none
        REAL(r8), dimension(:), intent(in) :: qv
        REAL(r8), dimension(:), intent(in) :: t
        REAL(r8), dimension(:), intent(in) :: u
        REAL(r8), dimension(:), intent(in) :: v
        REAL(r8), dimension(:), intent(in) :: omega
        REAL(r8), intent(in) :: ps
        REAL(r8), intent(in) :: ts
        REAL(r8), intent(in) :: solin 
        REAL(r8), intent(in) :: shflx
        REAL(r8), intent(in) :: lhflx

        REAL(r8), intent(in) :: landfrac 
        REAL(r8), intent(in) :: ocnfrac
        REAL(r8), intent(in) :: icefrac

        REAL(r8), intent(out) :: fsns
        REAL(r8), intent(out) :: flns
        REAL(r8), intent(out) :: fsnt
        REAL(r8), intent(out) :: flnt
        REAL(r8), intent(out) :: fsds
        logical , intent(in) :: debug_log 

        REAL(r8), intent(out) :: flds
        REAL(r8), intent(out) :: srfrad
        REAL(r8), intent(out) :: soll
        REAL(r8), intent(out) :: sols
        REAL(r8), intent(out) :: solld
        REAL(r8), intent(out) :: solsd


        REAL(r8), dimension(:), intent(out) :: pteq
        REAL(r8), dimension(:), intent(out) :: pttend
        REAL(r8), intent(out) :: prect
        REAL(r8), intent(out) :: precc
        REAL(r8), intent(out) :: precl
        REAL(r8), intent(out) :: precsc
        REAL(r8), intent(out) :: precsl

        !note the implicit conversion from R8 to c_float - since the NN is based
        !on single prcision output from cesm cam, but is defied as capi type
        REAL(c_float), dimension((size(t) * 5) + 8) :: istacked
        REAL(c_float), dimension((size(t) * 2) + 5 + 6 + 5) :: ostacked
        integer :: sz = 0
        integer :: retval = 0
        integer :: output_size = 0

        if (debug_log) then
            write(*,*) "torch-plugin debug log, 3D inputs:"
            write(*,*) "----------------------------------"
            write(*,*) "temp : ", t
            write(*,*) "Qv   : ", qv
            write(*,*) "U    : ", u
            write(*,*) "V    : ", v
            write(*,*) "omega: ", omega
            write(*,*) ""
            write(*,*) "torch-plugin debug log, 2D inputs:"
            write(*,*) "----------------------------------"
            write(*,*) "ps   : ", ps
            write(*,*) "ts   : ", ts
            write(*,*) "solin: ", solin
            write(*,*) "shflx: ", shflx
            write(*,*) "lhflx: ", lhflx
            write(*,*) ""
         end if
        !write(*,*) "here 1"

        sz = size(t)
        output_size = size(ostacked)

        istacked((sz*0) + 1 : sz*1) = qv
        istacked((sz*1) + 1 : sz*2) = t
        istacked((sz*2) + 1 : sz*3) = u
        istacked((sz*3) + 1 : sz*4) = v
        istacked((sz*4) + 1 : sz*5) = omega

        istacked((sz*6) + 1) = ps
        istacked((sz*6) + 2) = solin
        istacked((sz*6) + 3) = shflx
        istacked((sz*6) + 4) = lhflx
        istacked((sz*6) + 5) = landfrac
        istacked((sz*6) + 6) = ocnfrac
        istacked((sz*6) + 7) = icefrac
        istacked((sz*6) + 8) = ts


        call model_predict_c(                         &
                        input  = istacked,            &
                        input_size  = size(istacked), &
                        output = ostacked,            &
                        output_size = output_size,    &
                        retval = retval,              &
                        loopback = 0                  & !non-zero means test mode
                )


        pttend  = ostacked((sz*0) + 1 : sz*1)
        pteq    = ostacked((sz*1) + 1 : sz*2)

        fsns    = ostacked((sz*2) + 1)
        flns    = ostacked((sz*2) + 2)
        fsnt    = ostacked((sz*2) + 3)
        flnt    = ostacked((sz*2) + 4)
        fsds    = ostacked((sz*2) + 5)

        flds    = ostacked((sz*2) + 6)
        srfrad  = ostacked((sz*2) + 7)
        soll    = ostacked((sz*2) + 8)
        sols    = ostacked((sz*2) + 9)
        solld   = ostacked((sz*2) + 10)
        solsd   = ostacked((sz*2) + 11)

        prect   = ostacked((sz*2) + 12)
        precc   = ostacked((sz*2) + 13)
        precl   = ostacked((sz*2) + 14)
        precsc  = ostacked((sz*2) + 15)
        precsl  = ostacked((sz*2) + 16)


        if (debug_log) then
            !if (.true.) then
            write(*,*) "torch-plugin debug log, 3D outputs:"
            write(*,*) "-----------------------------------"
            write(*,*) "" 
            write(*,*) "ostacked : ", ostacked 
            write(*,*) "pttend: ", pttend
            write(*,*) "pteq  : ", pteq
            write(*,*) "" 
            write(*,*) "torch-plugin debug log, 2D outputs:"
            write(*,*) "-----------------------------------"
            write(*,*) "prect : ", prect
            write(*,*) "precc : ", precc
            write(*,*) "fsds  : ", fsds
            write(*,*) "flds  : ", flds
            write(*,*) "srfrad: ", srfrad
            write(*,*) "" 
            !write(*,*) "ostacked : ", ostacked 
         end if

    end subroutine


    subroutine predict_model_v5( &
            !inputs on vertical
            relhum, qv, t, u, v, omega, z3, &
            
            !inputs on 1d
            ps, ts, solin, shflx, lhflx, & 

            !added 1d
            landfrac, ocnfrac, icefrac, &

            !out on 1d
            fsns, flns, fsnt, flnt, fsds, &
            
            ! 1d outputs related to radiation 
            flds, srfrad, soll, sols, solld, solsd, &

            !outputs on vertical
            pteq, pttend, &            !+

            ! outputs on 1d
            prect, precc, precl, &
            precsc, precsl, &

            debug_log &

            )
        implicit none
        REAL(r8), dimension(:), intent(in) :: relhum
        REAL(r8), dimension(:), intent(in) :: qv
        REAL(r8), dimension(:), intent(in) :: t
        REAL(r8), dimension(:), intent(in) :: u
        REAL(r8), dimension(:), intent(in) :: v
        REAL(r8), dimension(:), intent(in) :: omega
        REAL(r8), dimension(:), intent(in) :: z3
        REAL(r8), intent(in) :: ps
        REAL(r8), intent(in) :: ts
        REAL(r8), intent(in) :: solin 
        REAL(r8), intent(in) :: shflx
        REAL(r8), intent(in) :: lhflx

        REAL(r8), intent(in) :: landfrac 
        REAL(r8), intent(in) :: ocnfrac
        REAL(r8), intent(in) :: icefrac

        REAL(r8), intent(out) :: fsns
        REAL(r8), intent(out) :: flns
        REAL(r8), intent(out) :: fsnt
        REAL(r8), intent(out) :: flnt
        REAL(r8), intent(out) :: fsds
        logical , intent(in) :: debug_log 

        REAL(r8), intent(out) :: flds
        REAL(r8), intent(out) :: srfrad
        REAL(r8), intent(out) :: soll
        REAL(r8), intent(out) :: sols
        REAL(r8), intent(out) :: solld
        REAL(r8), intent(out) :: solsd


        REAL(r8), dimension(:), intent(out) :: pteq
        REAL(r8), dimension(:), intent(out) :: pttend
        REAL(r8), intent(out) :: prect
        REAL(r8), intent(out) :: precc
        REAL(r8), intent(out) :: precl
        REAL(r8), intent(out) :: precsc
        REAL(r8), intent(out) :: precsl

        !note the implicit conversion from R8 to c_float - since the NN is based
        !on single prcision output from cesm cam, but is defied as capi type
        REAL(c_float), dimension((size(t) * 7) + 8) :: istacked
        REAL(c_float), dimension((size(t) * 2) + 5 + 6 + 5) :: ostacked
        integer :: sz = 0
        integer :: retval = 0
        integer :: output_size = 0

        if (debug_log) then
            write(*,*) "torch-plugin debug log, 3D inputs:"
            write(*,*) "----------------------------------"
            write(*,*) "temp : ", t
            write(*,*) "Qv   : ", qv
            write(*,*) "U    : ", u
            write(*,*) "V    : ", v
            write(*,*) "omega: ", omega
            write(*,*) "z3   : ", z3
            write(*,*) ""
            write(*,*) "torch-plugin debug log, 2D inputs:"
            write(*,*) "----------------------------------"
            write(*,*) "ps   : ", ps
            write(*,*) "ts   : ", ts
            write(*,*) "solin: ", solin
            write(*,*) "shflx: ", shflx
            write(*,*) "lhflx: ", lhflx
            write(*,*) ""
         end if
        !write(*,*) "here 1"

        sz = size(t)
        output_size = size(ostacked)

        istacked((sz*0) + 1 : sz*1) = relhum
        istacked((sz*0) + 1 : sz*2) = qv
        istacked((sz*1) + 1 : sz*3) = t
        istacked((sz*2) + 1 : sz*4) = u
        istacked((sz*3) + 1 : sz*5) = v
        istacked((sz*4) + 1 : sz*6) = omega
        istacked((sz*5) + 1 : sz*7) = z3

        istacked((sz*7) + 1) = ps
        istacked((sz*7) + 2) = solin
        istacked((sz*7) + 3) = shflx
        istacked((sz*7) + 4) = lhflx
        istacked((sz*7) + 5) = landfrac
        istacked((sz*7) + 6) = ocnfrac
        istacked((sz*7) + 7) = icefrac
        istacked((sz*7) + 8) = ts


        call model_predict_c(                         &
                        input  = istacked,            &
                        input_size  = size(istacked), &
                        output = ostacked,            &
                        output_size = output_size,    &
                        retval = retval,              &
                        loopback = 0                  & !non-zero means test mode
                )


        pttend  = ostacked((sz*0) + 1 : sz*1)
        pteq    = ostacked((sz*1) + 1 : sz*2)

        fsns    = ostacked((sz*2) + 1)
        flns    = ostacked((sz*2) + 2)
        fsnt    = ostacked((sz*2) + 3)
        flnt    = ostacked((sz*2) + 4)
        fsds    = ostacked((sz*2) + 5)

        flds    = ostacked((sz*2) + 6)
        srfrad  = ostacked((sz*2) + 7)
        soll    = ostacked((sz*2) + 8)
        sols    = ostacked((sz*2) + 9)
        solld   = ostacked((sz*2) + 10)
        solsd   = ostacked((sz*2) + 11)

        prect   = ostacked((sz*2) + 12)
        precc   = ostacked((sz*2) + 13)
        precl   = ostacked((sz*2) + 14)
        precsc  = ostacked((sz*2) + 15)
        precsl  = ostacked((sz*2) + 16)


        if (debug_log) then
            !if (.true.) then
            write(*,*) "torch-plugin debug log, 3D outputs:"
            write(*,*) "-----------------------------------"
            write(*,*) "" 
            write(*,*) "ostacked : ", ostacked 
            write(*,*) "pttend: ", pttend
            write(*,*) "pteq  : ", pteq
            write(*,*) "" 
            write(*,*) "torch-plugin debug log, 2D outputs:"
            write(*,*) "-----------------------------------"
            write(*,*) "prect : ", prect
            write(*,*) "precc : ", precc
            write(*,*) "fsds  : ", fsds
            write(*,*) "flds  : ", flds
            write(*,*) "srfrad: ", srfrad
            write(*,*) "" 
            !write(*,*) "ostacked : ", ostacked 
         end if
    end subroutine


    subroutine predict_model_v6( &
            !inputs on vertical
            qv, t, u, v, omega, z3, &
            cldliq, cldice, &

            !inputs on 1d
            ps, ts, solin, shflx, lhflx, & 

            !added 1d
            landfrac, ocnfrac, icefrac, &

            !out on 1d
            fsns, flns, fsnt, flnt, fsds, &
            
            ! 1d outputs related to radiation 
            flds, srfrad, soll, sols, solld, solsd, &

            !outputs on vertical
            pteq, pttend, &            !+
            ptecldliq, ptecldice, &

            ! outputs on 1d
            prect, precc, precl, &
            precsc, precsl, &

            debug_log &

            )
        implicit none
        REAL(r8), dimension(:), intent(in) :: qv
        REAL(r8), dimension(:), intent(in) :: t
        REAL(r8), dimension(:), intent(in) :: u
        REAL(r8), dimension(:), intent(in) :: v
        REAL(r8), dimension(:), intent(in) :: omega
        REAL(r8), dimension(:), intent(in) :: z3
        REAL(r8), dimension(:), intent(in) :: cldliq
        REAL(r8), dimension(:), intent(in) :: cldice
        REAL(r8), intent(in) :: ps
        REAL(r8), intent(in) :: ts
        REAL(r8), intent(in) :: solin 
        REAL(r8), intent(in) :: shflx
        REAL(r8), intent(in) :: lhflx

        REAL(r8), intent(in) :: landfrac 
        REAL(r8), intent(in) :: ocnfrac
        REAL(r8), intent(in) :: icefrac

        REAL(r8), intent(out) :: fsns
        REAL(r8), intent(out) :: flns
        REAL(r8), intent(out) :: fsnt
        REAL(r8), intent(out) :: flnt
        REAL(r8), intent(out) :: fsds
        logical , intent(in) :: debug_log 

        REAL(r8), intent(out) :: flds
        REAL(r8), intent(out) :: srfrad
        REAL(r8), intent(out) :: soll
        REAL(r8), intent(out) :: sols
        REAL(r8), intent(out) :: solld
        REAL(r8), intent(out) :: solsd


        REAL(r8), dimension(:), intent(out) :: pteq
        REAL(r8), dimension(:), intent(out) :: ptecldliq
        REAL(r8), dimension(:), intent(out) :: ptecldice
        REAL(r8), dimension(:), intent(out) :: pttend
        REAL(r8), intent(out) :: prect
        REAL(r8), intent(out) :: precc
        REAL(r8), intent(out) :: precl
        REAL(r8), intent(out) :: precsc
        REAL(r8), intent(out) :: precsl

        !note the implicit conversion from R8 to c_float - since the NN is based
        !on single prcision output from cesm cam, but is defied as capi type
        REAL(c_float), dimension((size(t) * 8) + 8) :: istacked
        REAL(c_float), dimension((size(t) * 4) + 5 + 6 + 5) :: ostacked
        integer :: sz = 0
        integer :: retval = 0
        integer :: output_size = 0

        if (debug_log) then
            write(*,*) "torch-plugin debug log, 3D inputs:"
            write(*,*) "----------------------------------"
            write(*,*) "temp  : ", t
            write(*,*) "Qv    : ", qv
            write(*,*) "U     : ", u
            write(*,*) "V     : ", v
            write(*,*) "omega : ", omega
            write(*,*) "z3    : ", z3
            write(*,*) "cldliq: ", cldliq
            write(*,*) "cldice: ", cldice
            write(*,*) ""
            write(*,*) "torch-plugin debug log, 2D inputs:"
            write(*,*) "----------------------------------"
            write(*,*) "ps   : ", ps
            write(*,*) "ts   : ", ts
            write(*,*) "solin: ", solin
            write(*,*) "shflx: ", shflx
            write(*,*) "lhflx: ", lhflx
            write(*,*) ""
         end if
        !write(*,*) "here 1"

        sz = size(t)
        output_size = size(ostacked)

        istacked((sz*0) + 1 : sz*1) = qv
        istacked((sz*1) + 1 : sz*2) = t
        istacked((sz*2) + 1 : sz*3) = u
        istacked((sz*3) + 1 : sz*4) = v
        istacked((sz*4) + 1 : sz*5) = omega
        istacked((sz*5) + 1 : sz*6) = z3
        istacked((sz*6) + 1 : sz*7) = cldliq
        istacked((sz*7) + 1 : sz*8) = cldliq

        istacked((sz*8) + 1) = ps
        istacked((sz*8) + 2) = solin
        istacked((sz*8) + 3) = shflx
        istacked((sz*8) + 4) = lhflx
        istacked((sz*8) + 5) = landfrac
        istacked((sz*8) + 6) = ocnfrac
        istacked((sz*8) + 7) = icefrac
        istacked((sz*8) + 8) = ts


        call model_predict_c(                         &
                        input  = istacked,            &
                        input_size  = size(istacked), &
                        output = ostacked,            &
                        output_size = output_size,    &
                        retval = retval,              &
                        loopback = 0                  & !non-zero means test mode
                )


        pttend    = ostacked((sz*0) + 1 : sz*1)
        pteq      = ostacked((sz*1) + 1 : sz*2)
        ptecldliq = ostacked((sz*2) + 1 : sz*3)
        ptecldice = ostacked((sz*3) + 1 : sz*4)

        fsns    = ostacked((sz*4) + 1)
        flns    = ostacked((sz*4) + 2)
        fsnt    = ostacked((sz*4) + 3)
        flnt    = ostacked((sz*4) + 4)
        fsds    = ostacked((sz*4) + 5)

        flds    = ostacked((sz*4) + 6)
        srfrad  = ostacked((sz*4) + 7)
        soll    = ostacked((sz*4) + 8)
        sols    = ostacked((sz*4) + 9)
        solld   = ostacked((sz*4) + 10)
        solsd   = ostacked((sz*4) + 11)

        prect   = ostacked((sz*4) + 12)
        precc   = ostacked((sz*4) + 13)
        precl   = ostacked((sz*4) + 14)
        precsc  = ostacked((sz*4) + 15)
        precsl  = ostacked((sz*4) + 16)


        if (debug_log) then
            !if (.true.) then
            write(*,*) "torch-plugin debug log, 3D outputs:"
            write(*,*) "-----------------------------------"
            write(*,*) "" 
            write(*,*) "ostacked : ", ostacked 
            write(*,*) "pttend: ", pttend
            write(*,*) "pteq  : ", pteq
            write(*,*) "ptecldliq  : ", ptecldliq
            write(*,*) "ptecldice  : ", ptecldice
            write(*,*) "" 
            write(*,*) "torch-plugin debug log, 2D outputs:"
            write(*,*) "-----------------------------------"
            write(*,*) "prect : ", prect
            write(*,*) "precc : ", precc
            write(*,*) "fsds  : ", fsds
            write(*,*) "flds  : ", flds
            write(*,*) "srfrad: ", srfrad
            write(*,*) "" 
            !write(*,*) "ostacked : ", ostacked 
         end if
    end subroutine



end module

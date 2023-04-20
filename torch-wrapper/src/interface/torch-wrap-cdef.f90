
! C functions declaration
interface
        function model_create_c(path) bind(C, name="torch_wrap_create")
                use iso_c_binding
                implicit none
                type(c_ptr) :: model_create_c
                character(len=1, kind=C_CHAR), intent(in) :: path(*)
        end function

        subroutine model_predict_c(input, input_size, output, output_size, retval, loopback) bind(C, name="torch_wrap_predict")
                use iso_c_binding
                implicit none
                real(c_float), dimension(*), intent(in) ::  input
                integer(c_int), value, intent(in) :: input_size
                real(c_float), dimension(*), intent(out) ::  output
                integer(c_int), intent(inout) :: output_size
                integer(c_int), intent(out) :: retval
                integer(c_int), value, intent(in) :: loopback 
        end subroutine
                 

end interface

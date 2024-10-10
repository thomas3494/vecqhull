module fvqhull
    use, intrinsic :: iso_c_binding
    use iso_fortran_env, only: real64, int64

    interface
        integer(c_size_t) function VecQuickhull(n, Px, Py) &
                                                bind(C, name='VecQuickhull')
            use, intrinsic :: iso_c_binding
            integer (C_size_t), value :: n
            type (c_ptr), value :: Px, Py
        end function VecQuickhull

        integer(c_size_t) function VecQuickhullP(n, Px, Py) &
                                                bind(C, name='VecQuickhullP')
            use, intrinsic :: iso_c_binding
            integer (C_size_t), value :: n
            type (c_ptr), value :: Px
            type (c_ptr), value :: Py
        end function VecQuickhullP
    end interface

    contains
        function vqhull(n, P) result(m)
            integer(int64), intent(in) :: n
            real(real64), intent(inout), target :: P(:,:)
            m = VecQuickhull(n, C_loc(P(:, 1)), C_loc(P(:, 2)))
        end function vqhull

        function vqhullp(n, P) result(m)
            integer(int64), intent(in) :: n
            real(real64), intent(inout), target :: P(:,:)
            m = VecQuickhullP(n, C_loc(P(:, 1)), C_loc(P(:, 2)))
        end function vqhullp
end module

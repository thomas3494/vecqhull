program main
    use iso_fortran_env, only: real64, int64
    use fvqhull

    real(real64), allocatable :: P(:, :)
    integer(int64) :: n, m
    real(real64) :: two_pi_n

    n = 1e8
    two_pi_n = 8.d0 * datan(1.d0) / n

    allocate(P(n, 2))

    write(*, *) "Initialize circle with ", n, " points."

    do i = 1, n
        P(i, 1) = dcos(two_pi_n * i)
    end do
    do i = 1, n
        P(i, 2) = dsin(two_pi_n * i)
    end do

    write(*, *) "Start computing hull"

    m = vqhullp(n, P)

    write(*, *) "Done. Hull has ", m, " elements"

    deallocate(P)
end program

!model_natureをチェック
program maketruedata_binary
    implicit none

    real(8), parameter :: PI = 4*ATAN (1.d0)

    ! Constants for the two-way coupled system
    real(8), parameter :: h = 1.0d0
    real(8), parameter :: b = 10.0d0
    real(8), parameter :: c = 4.0d0

    ! Simulation parameters
    ! integer, parameter :: total_inits = 1000  ! アトラクタ上のデータの数
    integer, parameter :: data_span = 5
    integer, parameter :: n_timesteps = data_span*4*365*100  !6時間ごとのデータなので1日4つ。スピンアップは99年
    integer, parameter :: n_savedata = data_span*4*365
    integer, parameter :: kmax = 40
    integer, parameter :: jmax = 0
    real(8), parameter :: F_ex = 8.d0
    real(8), parameter :: dt = 1.0d-2
    real(8), parameter :: ADAY = data_span * 4 * dt


    ! Variables for storing x1, x2 data and selected random initial condition
    real(8) :: x1_true(kmax)
    ! real(8) :: x2_true(jmax*kmax)

    ! real(8) :: random_num
    ! integer :: random_index
    ! integer :: i, j, iostat, jmaxkmax
    integer :: l,m,n
    real(8) :: Ampli, Period, Freq

    ! Initialize state variables
    real(8) :: t

    ! File units
    integer, parameter :: iunit_x1 = 12
    ! integer, parameter :: iunit_x2 = 22
    integer, parameter :: iunit_t = 32
    integer, parameter :: iunit_x1t = 14
    ! integer, parameter :: iunit_x2t = 24

    ! ! Open binary files for input and output
    ! open(iunit_x1, file="1scale_x_initial_km40_F8.bin", form='unformatted', access='stream', action="read", iostat=iostat)
    ! if (iostat /= 0) then
    !     print *, "Error: Could not open x1 binary file."
    !     stop
    ! end if

    open(iunit_x1t, file="1scalePlusSinVarious_xtrue_timeseries_km40_F8_50span.bin", &
    & form='unformatted', access='stream', status="replace")

    open(iunit_t, file="1scalePlusSinVarious_t_truedata_fort_50span.bin", form='unformatted', access='stream', status="replace")

    do l = -1, 4
        do m = 1,15
            !! Initialize random seed and select a random index!!
            call random_seed()
            ! call random_number(random_num)
            ! random_index = int(random_num * total_inits) + 1  ! 1 〜 total_inits

            ! ! Rewind the file to the beginning before skipping!!
            ! rewind(iunit_x1)

            ! ! Move the file pointers to the random index !!
            ! call skip_records(iunit_x1, random_index - 1, kmax)

            ! ! Read the selected initial condition !!
            ! read(iunit_x1) x1_true  ! アトラクタ上の点 x1

            ! ! Output the selected random initial condition for verification !!
            ! print *, "Random initial condition (x1):", x1_true
            ! print *, "Random index:", random_index

            Ampli= 2.**l
            Period= 2.* m * dble(ADAY)
            Freq = 1./Period


            ! Time integration !
            t = 0.0d0
            call initialize(x1_true)

            do n = 1, n_timesteps

                ! Runge-Kutta 4th order integration step
                call RK4(t, x1_true, dt,Ampli,Freq)

                ! Update time
                t = t + dt

                if (n > n_timesteps - n_savedata .and. modulo(n, data_span) == 0) then
                    ! Write data to files
                    write(iunit_x1t) x1_true
                    write(iunit_t) t
                end if
            end do
        enddo
    enddo

    ! Close files
    close(iunit_x1)
    close(iunit_x1t)
    close(iunit_t)

contains


    ! Function to initialize x1 and x2
    subroutine initialize(x1)
        real(8), intent(out) :: x1(kmax)

        ! Initialize x1 and x2 with random values
        x1 = stdnorm_err(kmax)
    end subroutine initialize


    ! Subroutine to skip a given number of records in a binary file
    subroutine skip_records(unit, num_records, record_size)
        integer, intent(in) :: unit, num_records, record_size
        real(8) :: temp(record_size)
        integer :: i, iostat

        do i = 1, num_records
            read(unit, iostat=iostat) temp
            if (iostat /= 0) then
                print *, "Error: Failed to skip records in binary file."
                stop
            end if
        end do
    end subroutine skip_records

    ! Runge-Kutta 4th order integration step
    subroutine RK4(t, x1, dt, A, B)
        real(8), intent(out) :: x1(kmax)
        real(8), intent(in) :: t, dt
        real(8), intent(in) :: A, B
        real(8) :: a1(kmax), a2(kmax), a3(kmax), a4(kmax)
        real(8) :: temp_x1(kmax)

        call model_nature_1(t, x1, a1, A, B)

        temp_x1 = x1 + 0.5d0 * dt * a1
        call model_nature_1(t, temp_x1, a2, A, B)

        temp_x1 = x1 + 0.5d0 * dt * a2
        call model_nature_1(t, temp_x1, a3, A, B)

        temp_x1 = x1 + dt * a3
        call model_nature_1(t, temp_x1, a4, A, B)

        x1 = x1 + dt / 6.0d0 * (a1 + 2.0d0 * (a2 + a3) + a4)
    end subroutine RK4

    ! Model for x1
    ! A: Amplitude, B: Frequency
    subroutine model_nature_1(t, X1, model, A, B)
        real(8), intent(in) :: t, X1(kmax)
        real(8), intent(out) :: model(kmax)
        real(8), intent(in) :: A, B
        integer :: k

        do k = 1, kmax
            model(k) = X1(modulo(k - 2, kmax) + 1) * (X1(modulo(k, kmax) + 1) - X1(modulo(k - 3, kmax) + 1))  &
                       - X1(k) + F_ex + A * sin(2.0d0* PI* B * t)
        end do
    end subroutine model_nature_1


    function stdnorm_err(n) result(z1)
        implicit none
        integer, intent(in) :: n
        real(8) :: z1(n)
        real(8) :: u1, u2
        integer :: i

        ! Box-Muller method for generating normal random numbers
        do i = 1, n
            call random_number(u1)
            call random_number(u2)
            z1(i) = sqrt(-2.0d0 * log(u1)) * cos(2.0d0 * acos(-1.0d0) * u2)
        end do
    end function stdnorm_err
end program maketruedata_binary

      
!   This program uses internal subprograms for the Cahn-Hilliard eq. 
!            
!         
!   Author :
!              Shahid Maqbool
! 
!   Modified :
!                13 Feb. 2023
!
!   To compile and run :
!                          Check ReadMe
!              
!------------------------------------------------------------------------------


program fd_ch_test
  implicit none

  ! ===========================================================================
  !                                parameters
  ! ===========================================================================

  ! simulation cell 

  integer ( kind = 4 ), parameter :: Nx = 64
  integer ( kind = 4 ), parameter :: Ny = 64
  integer ( kind = 4 ), parameter :: dx = 1
  integer ( kind = 4 ), parameter :: dy = 1

  ! time integration 

  integer ( kind = 4 ), parameter :: no_of_steps = 10000
  integer ( kind = 4 ), parameter :: frequency   = 1000
  integer ( kind = 4 )            :: step 
  real ( kind = 8 )   , parameter :: dt = 0.01
  real ( kind = 8 )               :: start, finish

  ! material specific 

  real ( kind = 8 )   , parameter :: initial_con = 0.4
  real ( kind = 8 )   , parameter :: mobility = 1.0
  real ( kind = 8 )   , parameter :: grad_coef = 0.5

  ! microstructure 

  real ( kind = 8 )   , parameter :: noise = 0.02
  real ( kind = 8 )   , parameter :: A  = 1.0
  integer ( kind = 4 )            :: i, j, jp, jm, ip, im
  real ( kind = 8 )   , dimension ( Nx, Ny ) :: con, lap_con, dfdcon
  real ( kind = 8 )   , dimension ( Nx, Ny ) :: dummy_con, lap_dummy


  call cpu_time ( start )



  ! ===========================================================================
  !                            initial microstucture
  ! ===========================================================================



  call Introduce_fluctuation ( con, initial_con, noise )



  ! ===========================================================================
  !                         evolution of microstructure 
  ! ===========================================================================


  time_loop: do step = 1, no_of_steps
  

     spatial_loop:  do concurrent ( j = 1 : Nx , i = 1 : Ny )
	 

        call Set_boundary_conditions (i, j, jp, jm, ip, im )

        call Compute_derivative_free_energy ( A, con, dfdcon, i, j )

        call Evaluate_laplacian ( lap_con, con, dx, dy, dummy_con, &
             & dfdcon, grad_coef, lap_dummy, i, j, ip, jp, im, jm )

        call Perform_time_integration ( con, dt, mobility, lap_dummy, i, j )
		

     end do spatial_loop


     ! adjust concentration in range

     where ( con >= 0.99999 )  con = 0.99999
     where ( con <  0.00001 )  con = 0.00001

     ! print steps

     if ( mod ( step, frequency ) .eq. 0 ) print *, 'Done steps  = ', step


  end do time_loop


  call cpu_time ( finish )



  ! ===========================================================================
  !                                  Output 
  ! ===========================================================================



  call Write_input_parameters_on_file
  call Output_concentration_on_file


contains



  ! ===========================================================================
  !                           internal subprograms 
  ! ===========================================================================



  subroutine Introduce_fluctuation ( con_, initial_con_, noise_ )
    implicit none

    real ( kind = 8 ), intent ( in )                        :: noise_
    real ( kind = 8 ), intent ( in )                        :: initial_con_
    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( out ) :: con_
    real ( kind = 8 ), dimension ( Nx, Ny )                 :: r_

    call random_number ( r_ )

    con_ = initial_con_ + noise_*( 0.5 - r_ )

  end subroutine Introduce_fluctuation



  ! ---------------------------------------------------------------------------



  pure subroutine Set_boundary_conditions ( i_, j_, jp_, jm_, ip_, im_ )
    implicit none

    integer ( kind = 4 ), intent ( in )  :: i_, j_
    integer ( kind = 4 ), intent ( out ) :: jp_, jm_, ip_, im_

    jp_ = j_ + 1
    jm_ = j_ - 1

    ip_ = i_ + 1
    im_ = i_ - 1

    if ( im_ == 0 ) im_ = Nx
    if ( ip_ == ( Nx + 1 ) ) ip_ = 1
    if ( jm_ == 0 ) jm_ = Ny
    if ( jp_ == ( Ny + 1 ) ) jp_ = 1


  end subroutine Set_boundary_conditions



  ! ---------------------------------------------------------------------------



  pure subroutine Compute_derivative_free_energy ( A_, con_, dfdcon_, i_, j_ )
    implicit none


    real ( kind = 8 ), intent ( in )                        :: A_
    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( in )  :: con_
    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( out ) :: dfdcon_
    integer ( kind = 4 ), intent ( in )                     :: i_, j_


    dfdcon_(i_,j_) = A_*( 2.0*con_(i_,j_)*( 1.0 - con_(i_,j_) )**2 - &
         & 2.0*con_(i_,j_)**2*( 1.0 - con_(i_,j_) ) )


  end subroutine Compute_derivative_free_energy



  ! ---------------------------------------------------------------------------



  pure subroutine Evaluate_laplacian ( lap_con_, con_, dx_, dy_, dummy_con_, &
       & dfdcon_, grad_coef_, lap_dummy_, i_ , j_, ip_, jp_, im_, jm_ )


    real ( kind = 8 ), intent ( in )                       :: grad_coef_
    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( out ):: lap_con_
    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( out ):: dummy_con_
    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( in ) :: con_, dfdcon_
    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( out ):: lap_dummy_
    integer ( kind = 4 ), intent ( in )                    :: dx_, dy_, i_, j_
    integer ( kind = 4 ), intent ( in )                    :: ip_, jp_
    integer ( kind = 4 ), intent ( in )                    :: im_, jm_


    lap_con_(i_,j_) = ( con(ip_,j_) + con(im_,j_) + con(i_,jm_) + con(i_,jp_) &
         & - 4.0*con(i_,j_) ) / ( dx_*dy_ )

    dummy_con_(i_,j_) = dfdcon_(i_,j_) - grad_coef_*lap_con_(i_,j_)

    lap_dummy_(i_,j_) = ( dummy_con_(ip_,j_) + dummy_con_(im_,j_) & 
         & + dummy_con_(i_,jm_) + dummy_con_(i_,jp_) &
         & - 4.0*dummy_con_(i_,j_) ) / ( dx_*dy_ )


  end subroutine Evaluate_laplacian



  ! ---------------------------------------------------------------------------



  pure subroutine Perform_time_integration ( con_, dt_, mobility_, &
       & lap_dummy_, i_, j_ )
    implicit none


    real ( kind = 8 ), intent ( in )                           :: dt_
    real ( kind = 8 ), intent ( in )                           :: mobility_
    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( in )     :: lap_dummy_
    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( in out ) :: con_
    integer ( kind = 4 ), intent ( in )                        :: i_, j_


    con_(i_,j_)  = con_(i_,j_) + dt_*mobility_*lap_dummy_(i_,j_)


  end subroutine Perform_time_integration



  ! ---------------------------------------------------------------------------



  subroutine Write_input_parameters_on_file
    open ( unit = 1, file = 'parameters.txt', status = 'replace' )

    write ( 1, * ) ""
    write ( 1, * ) '------------------------------------------------'
    write ( 1, * ) '                 Parameters                     '
    write ( 1, * ) '------------------------------------------------'
    write ( 1, * ) ""
    write ( 1, 100 ) Nx
    write ( 1, 200 ) Ny
    write ( 1, 300 ) no_of_steps
    write ( 1, 400 ) dt
    write ( 1, 500 ) initial_con
    write ( 1, 600 ) mobility
    write ( 1, 700 ) grad_coef
    write ( 1, 800 ) noise
    write ( 1, 900 ) A
    write ( 1, 1000 ) finish - start
    write ( 1, * ) ""
    write ( 1, * ) '------------------------------------------------'

    close ( 1 )

100 format ( 6x, 'Nx', 23x, ' = ', i2 )
200 format ( 6x, 'Ny', 23x, ' = ', i2 )
300 format ( 6x, 'number of steps', 10x, ' = ', i5 )
400 format ( 6x, 'time difference', 10x, ' = ',f4.2 )
500 format ( 6x, 'initial concentration', 4x, ' = ', f3.1 )
600 format ( 6x, 'mobility', 17x, ' = ', f3.1 )
700 format ( 6x, 'gradient coefficient', 5x, ' = ', f3.1 )
800 format ( 6x, 'noise', 20x' = ', f4.2 )
900 format ( 6x, 'barrier height', 11x' = ', f3.1 )
1000 format( 6x, 'computed time in seconds',1x' = ', f4.2 )

    
  end subroutine Write_input_parameters_on_file



  ! ---------------------------------------------------------------------------



  subroutine Output_concentration_on_file
    implicit none

    character ( len = 80 ) :: filename 

    write ( filename, '( "con_", i0.3,".dat" )' ) no_of_steps
    open ( 2, file = filename, status = 'replace' )

    do i = 1, Nx    
       write ( 2, * ) ( con(i,j), j = 1, Ny )
    end do

    close ( 2 )

  end subroutine Output_concentration_on_file

  

end program fd_ch_test

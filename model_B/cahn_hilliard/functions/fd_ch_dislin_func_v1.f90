      
!   This program uses internal subprograms for the Cahn-Hilliard eq. 
!            
!         
!   Author :
!              Shahid Maqbool
! 
!   Modified :
!                13 Feb. 2023, 10 August 2023
!
!   To compile and run :
!                          Check ReadMe
!              
!------------------------------------------------------------------------------


program fd_ch_test
  use Dislin
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
  real ( kind = 8 )   , dimension ( Nx, Ny ) :: dummy_con, laplace


  call cpu_time ( start )


  ! ===========================================================================
  !                            initial microstructure
  ! ===========================================================================



  con =  Introduce_fluctuation( initial_con, noise )



  ! ===========================================================================
  !                         evolution of microstructure 
  ! ===========================================================================


  time_loop: do step = 1, no_of_steps
  

     column: do j = 1, Ny
        row: do i = 1, Ny
		

           dfdcon = Deriv_free_energy ( A, con, i, j )

           laplace = Laplacian ( lap_con, con, dx, dy, dummy_con, &
                & dfdcon, grad_coef, i, j, ip, jp, im, jm )

           con(i,j)  = con(i,j) + dt*mobility*laplace(i,j)


        end do row
     end do column



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
  call Dislin_color_plot ( )


contains



  ! ===========================================================================
  !                          internal subprograms 
  ! ===========================================================================



  function Introduce_fluctuation (  initial_con_, noise_ )
    implicit none

    real ( kind = 8 ), intent ( in )        :: noise_
    real ( kind = 8 ), intent ( in )        :: initial_con_
    real ( kind = 8 ), dimension ( Nx, Ny ) :: Introduce_fluctuation 
    real ( kind = 8 ), dimension ( Nx, Ny ) :: r

    call random_number ( r )

    Introduce_fluctuation = initial_con_ + noise_*( 0.5 - r )

    
  end function Introduce_fluctuation



  ! ---------------------------------------------------------------------------



  pure function deriv_free_energy ( A_, con_, i_, j_ )
    implicit none


    real ( kind = 8 ), dimension ( Nx, Ny ) :: deriv_free_energy
    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( in )  :: con_
    real ( kind = 8 ), intent ( in )                        :: A_    
    integer ( kind = 4 ), intent ( in )                     :: i_, j_


    Deriv_free_energy(i_,j_) = A_*( 2.0*con_(i_,j_)*( 1.0 - con_(i_,j_) )**2 - &
         & 2.0*con_(i_,j_)**2*( 1.0 - con_(i_,j_) ) )


  end function deriv_free_energy



  ! ---------------------------------------------------------------------------



  function Laplacian ( lap_con_, con_, dx_, dy_, dummy_con_, &
       & dfdcon_, grad_coef_, i_ , j_, ip_, jp_, im_, jm_ )
    implicit none

    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( in out ):: lap_con_
    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( in out ):: dummy_con_
    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( in out ):: dfdcon_
    real ( kind = 8 ), intent ( in )                       :: grad_coef_
    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( in ) :: con_ 
    real ( kind = 8 ), dimension ( Nx, Ny )                :: laplacian
    integer ( kind = 4 ), intent ( in )                    :: dx_, dy_, i_, j_
    integer ( kind = 4 ), intent ( inout )                 :: ip_, jp_
    integer ( kind = 4 ), intent ( inout )                 :: im_, jm_


    jp_ = j_ + 1
    jm_ = j_ - 1

    ip_ = i_ + 1
    im_ = i_ - 1

    if ( im_ == 0 ) im_ = Nx
    if ( ip_ == ( Nx + 1 ) ) ip_ = 1
    if ( jm_ == 0 ) jm_ = Ny
    if ( jp_ == ( Ny + 1 ) ) jp_ = 1


    lap_con_(i_,j_) = ( con_(ip_,j_) + con_(im_,j_) + con_(i_,jm_) + con_(i_,jp_) &
         & - 4.0*con_(i_,j_) ) / ( dx_*dy_ )

    dummy_con_(i_,j_) = dfdcon_(i_,j_) - grad_coef_*lap_con_(i_,j_)

    Laplacian(i_,j_) = ( dummy_con_(ip_,j_) + dummy_con_(im_,j_) & 
         & + dummy_con_(i_,jm_) + dummy_con_(i_,jp_) &
         & - 4.0*dummy_con_(i_,j_) ) / ( dx_*dy_ )


  end function laplacian



  ! ---------------------------------------------------------------------------



  subroutine Write_input_parameters_on_file
    open ( unit = 1, file = 'parameters.txt', status = 'replace' )

    write ( 1, * ) ""
    write ( 1, * ) '------------------------------------------------'
    write ( 1, * ) '  Input Parameters for Cahn-Hilliard Equation   '
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


  
  ! ---------------------------------------------------------------------------

  

  subroutine Dislin_color_plot ( )

    call scrmod ( 'REVERS' )
    call metafl ( 'png' )

    call disini ( )

    call hwfont ( )
    call titlin ( 'Color Plot', 4 )

    call name ( 'Nx', 'X' )
    call name ( 'Ny', 'Y' )
    call name ( 'Concentration', 'Z' )

    call intax ( )
    call autres ( Nx, Ny )
    call axspos ( 350, 1700 )
    call ax3len ( 1400, 1400, 1400 )

    call labdig ( 2, 'Z' )        
    call graf3 ( 0.d0, 64.d0, 0.d0, 16.d0, 0.d0, 64.d0, 0.d0, 16.d0, &
         & 0.05d0, 1.0d0, 0.05d0, 0.05d0 )
    call crvmat ( con, Nx, Ny, 1, 1 )

    call height ( 50 )
    call title ( )
    call mpaepl ( 3 )

    call disfin ( )

    stop

  end subroutine Dislin_color_plot



  
end program fd_ch_test

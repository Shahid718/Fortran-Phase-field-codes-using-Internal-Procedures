      
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

  ! time integration 

  integer ( kind = 4 ), parameter :: no_of_steps = 10000
  integer ( kind = 4 ), parameter :: frequency   = 1000
  integer ( kind = 4 )            :: step 
  real ( kind = 8 )   , parameter :: dt = 0.01
  real ( kind = 8 )               :: start, finish

  ! material specific 

  real ( kind = 8 )   , parameter :: initial_con = 0.4
  real ( kind = 8 )   , parameter :: mobility = 1.0

  ! microstructure 

  real ( kind = 8 )   , dimension ( Nx, Ny ) :: con, lap_con, dfdcon
  real ( kind = 8 )   , dimension ( Nx, Ny ) :: dummy_con, laplace
  integer ( kind = 4 ) :: i, j 


  call cpu_time ( start )


  ! ===========================================================================
  !                            initial microstucture
  ! ===========================================================================



  con =  Introduce_fluctuation( initial_con )



  ! ===========================================================================
  !                         evolution of microstructure 
  ! ===========================================================================


  time_loop: do step = 1, no_of_steps


     con  = con + dt*mobility*Laplacian( con, dfdcon, lap_con, dummy_con )


     ! adjust concentration in range

     where ( con >= 0.99999 )  con = 0.99999
     where ( con <  0.00001 )  con = 0.00001


     ! print steps on the console

     if ( mod ( step, frequency ) .eq. 0 ) print *, 'Done steps  = ', step


  end do time_loop


  call cpu_time ( finish )



  ! ===========================================================================
  !                                  Output 
  ! ===========================================================================



  call Output_files


contains



  ! ===========================================================================
  !                          internal subprograms 
  ! ===========================================================================


  ! ---------------------------------------------------------------------------



  function Introduce_fluctuation ( initial_con_ )
    implicit none

    real ( kind = 8 )  :: noise = 0.02
    real ( kind = 8 ), intent ( in )        :: initial_con_
    real ( kind = 8 ), dimension ( Nx, Ny ) :: Introduce_fluctuation 
    real ( kind = 8 ), dimension ( Nx, Ny ) :: r

    call random_number ( r )

    Introduce_fluctuation = initial_con_ + noise*( 0.5 - r )

  end function Introduce_fluctuation



  ! ---------------------------------------------------------------------------
  


  function Laplacian ( con_, dfdcon_, lap_con_, dummy_con_ )
    implicit none

    real ( kind = 8 ) , dimension ( Nx, Ny ), intent (in out) :: con_
    real ( kind = 8 ) , dimension ( Nx, Ny ), intent ( out )  :: dfdcon_
    real ( kind = 8 ) , dimension ( Nx, Ny ), intent ( out )  :: lap_con_
    real ( kind = 8 ) , dimension ( Nx, Ny ), intent ( out )  :: dummy_con_
    real ( kind = 8 ) , dimension ( Nx, Ny )  :: laplacian
    real ( kind = 8 )    :: grad_coef = 0.5
    real ( kind = 8 )    :: dt = 0.01
    real ( kind = 8 )    :: A  = 1.0
    integer ( kind = 4 ) :: i, j, jp, jm, ip, im, dx = 1, dy = 1


    do concurrent ( i = 1:Nx, j =1:Ny )


       ! free energy derivative

       dfdcon_(i,j) = A*( 2.0*con_(i,j)*( 1.0 - con_(i,j) )**2 &
            - 2.0*con_(i,j)**2*( 1.0 - con_(i,j) ) )

       ! laplace evaluation

       jp = j + 1
       jm = j - 1

       ip = i + 1
       im = i - 1

       if ( im == 0 ) im = Nx
       if ( ip == ( Nx + 1) ) ip = 1
       if ( jm == 0 ) jm = Ny
       if ( jp == ( Ny + 1) ) jp = 1

       lap_con_(i,j)   = ( con_(ip,j) + con_(im,j) + con_(i,jm) + &
            & con_(i,jp) - 4.0*con_(i,j) ) /( dx*dy )

       dummy_con_(i,j) = dfdcon_(i,j) - grad_coef*lap_con_(i,j)

       Laplacian(i,j) = ( dummy_con_(ip,j) + dummy_con_(im,j) + &
            & dummy_con_(i,jm) + dummy_con_(i,jp) - &
            & 4.0*dummy_con_(i,j) ) / ( dx*dy )


    end do

end function laplacian


  ! ---------------------------------------------------------------------------
  


  subroutine Output_files
    implicit none

    ! parameters file

    real ( kind = 8 )      :: noise = 0.02
    real ( kind = 8 )      :: initial_con = 0.4
    real ( kind = 8 )      :: mobility = 1.0
    real ( kind = 8 )      :: grad_coef = 0.5
    real ( kind = 8 )      :: dt = 0.01
    real ( kind = 8 )      :: A  = 1.0
    character ( len = 80 ) :: filename
    integer ( kind = 4 )   :: i, j

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


    ! con write file
    
    write ( filename, '( "con_", i0.3,".dat" )' ) no_of_steps
    open ( 2, file = filename, status = 'replace' )

    do i = 1, Nx    
       write ( 2, * ) ( con(i,j), j = 1, Ny )
    end do

    close ( 2 )

  end subroutine Output_files



end program fd_ch_test

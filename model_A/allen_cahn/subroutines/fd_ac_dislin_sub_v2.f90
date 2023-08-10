
!   This program uses internal subprograms for Allen-Cahn Equation.
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


program fd_ac_test
  use Dislin
  implicit none


  ! ===========================================================================
  !                                parameters
  ! ===========================================================================


  ! simulation cell 

  integer ( kind = 4 ), parameter :: Nx = 64
  integer ( kind = 4 ), parameter :: Ny = 64

  ! time integration

  integer ( kind = 4 ), parameter :: no_of_steps = 1500
  integer ( kind = 4)             :: frequency = 100
  integer (kind = 4 )             :: step 
  real ( kind = 8 )               :: start, finish

  ! microstructure

  real ( kind = 8 ), dimension ( Nx, Ny ) :: phi, dfdphi
  real ( kind = 8 ), dimension ( Nx, Ny ) :: lap_phi, dummy_phi


  call cpu_time ( start )



  ! ===========================================================================
  !                            initial microstructure
  ! ===========================================================================



  call Introduce_fluctuation ( phi)



  ! ===========================================================================
  !                         evolution of microstructure 
  ! ===========================================================================



  time_loop: do step = 1, no_of_steps


     call Perform_evolution ( phi , dfdphi, lap_phi, dummy_phi )


     ! adjust order parameter in range

     where ( phi >= 0.99999 ) phi = 0.99999
     where ( phi < 0.00001  ) phi = 0.00001


     ! print steps on console

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



  subroutine Introduce_fluctuation ( phi_ )
    implicit none

    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( out ) :: phi_
    real ( kind = 8 ), dimension ( Nx, Ny ) :: r
    real ( kind = 8 )                       :: noise = 0.02
    real ( kind = 8 )                       :: initial_phi = 0.5

    call random_number ( r )

    phi_ = initial_phi + noise*( 0.5 - r )

  end subroutine Introduce_fluctuation



  ! ---------------------------------------------------------------------------



  subroutine Perform_evolution ( phi_, dfdphi_, lap_phi_, dummy_phi_ )
    implicit none

    real ( kind = 8 ) , dimension ( Nx, Ny ), intent (in out) :: phi_
    real ( kind = 8 ) , dimension ( Nx, Ny ), intent ( out )  :: dfdphi_
    real ( kind = 8 ) , dimension ( Nx, Ny ), intent ( out )  :: lap_phi_
    real ( kind = 8 ) , dimension ( Nx, Ny ), intent ( out )  :: dummy_phi_
    real ( kind = 8 )    :: mobility = 1.0
    real ( kind = 8 )    :: grad_coef = 1.0
    real ( kind = 8 )    :: dt = 0.01
    real ( kind = 8 )    :: A  = 1.0
    integer ( kind = 4 ) :: i, j, jp, jm, ip, im, dx = 2, dy = 2


    do concurrent ( j = 1:Ny, i = 1:Nx )


       ! free energy derivative

       dfdphi_(i,j) = A*( 2.0*phi_(i,j)*( 1.0 - phi_(i,j) )**2 &
            *( 1.0 - 2*phi_(i,j) ) )

       ! laplace evaluation

       jp = j + 1
       jm = j - 1

       ip = i + 1
       im = i - 1

       if ( im == 0 ) im = Nx
       if ( ip == ( Nx + 1 ) ) ip = 1
       if ( jm == 0 ) jm = Ny
       if ( jp == ( Ny + 1 ) ) jp = 1

       lap_phi_(i,j) = ( phi_(ip,j) + phi_(im,j) + phi_(i,jm) + &
            phi_(i,jp) - 4.0*phi_(i,j)) /( dx*dy )              

       ! time integration

       dummy_phi_(i,j) = dfdphi_(i,j) - grad_coef*lap_phi_(i,j)

       phi_(i,j) = phi_(i,j) - dt*mobility*dummy_phi(i,j)


    end do


  end subroutine Perform_evolution
  

  
  ! ---------------------------------------------------------------------------
  


  subroutine Output_files
    implicit none

    ! parameters file

    real ( kind = 8 )      :: noise = 0.02
    real ( kind = 8 )      :: initial_phi = 0.5
    real ( kind = 8 )      :: mobility = 1.0
    real ( kind = 8 )      :: grad_coef = 1.0
    real ( kind = 8 )      :: dt = 0.01
    real ( kind = 8 )      :: A  = 1.0
    character ( len = 80 ) :: filename
    integer ( kind = 4 )   :: i, j

    open ( unit = 1, file = 'parameters.txt', status = 'replace' )


    write ( 1, * ) ""
    write ( 1, * ) '------------------------------------------------'
    write ( 1, * ) '   Input Parameters for Allen-Cahn Equation     '
    write ( 1, * ) '------------------------------------------------'
    write ( 1, * ) ""
    write ( 1, 100 ) Nx
    write ( 1, 200 ) Ny
    write ( 1, 300 ) no_of_steps
    write ( 1, 400 ) dt
    write ( 1, 500 ) initial_phi
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
300 format ( 6x, 'number of steps', 10x, ' = ', i4 )
400 format ( 6x, 'time difference', 10x, ' = ',f4.2 )
500 format ( 6x, 'initial phi', 14x, ' = ', f3.1 )
600 format ( 6x, 'mobility', 17x, ' = ', f3.1 )
700 format ( 6x, 'gradient coefficient', 5x, ' = ', f3.1 )
800 format ( 6x, 'noise', 20x' = ', f4.2 )
900 format ( 6x, 'barrier height', 11x' = ', f3.1 )
1000 format( 6x, 'computed time in seconds',1x' = ', f4.2 )


    ! phi file

    write ( filename, '( "phi_", i0.3,".dat" )' ) no_of_steps
    open ( 2, file = filename, status = 'replace' )

    do i = 1, Nx    
       write ( 2, * ) ( phi(i,j), j = 1, Ny )
    end do

    close ( 2 )


    ! dislin plot

    call scrmod ( 'REVERS' )
    call metafl ( 'png' )

    call disini ( )

    call hwfont ( )
    call titlin ( 'Color Plot', 4 )

    call name ( 'Nx', 'X' )
    call name ( 'Ny', 'Y' )
    call name ( 'phi', 'Z' )

    call intax ( )
    call autres ( Nx, Ny )
    call axspos ( 350, 1700 )
    call ax3len ( 1400, 1400, 1400 )

    call labdig ( 2, 'Z' )        
    call graf3 ( 0.d0, 64.d0, 0.d0, 16.d0, 0.d0, 64.d0, 0.d0, 16.d0, &
         & 0.05d0, 1.0d0, 0.05d0, 0.05d0 )
    call crvmat ( phi, Nx, Ny, 1, 1 )

    call height ( 50 )
    call title ( )
    call mpaepl ( 3 )

    call disfin ( )

    stop

  end subroutine Output_files



end program fd_ac_test

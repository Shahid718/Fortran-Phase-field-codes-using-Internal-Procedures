      
!    This program uses internal subprograms for the Allen-Cahn eq.
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


program fd_ac_test
  implicit none


  ! ===========================================================================
  !                                parameters
  ! ===========================================================================


  ! simulation cell 

  integer ( kind = 4 ), parameter :: Nx = 64
  integer ( kind = 4 ), parameter :: Ny = 64
  integer ( kind = 4 ), parameter :: dx = 2
  integer ( kind = 4 ), parameter :: dy = 2

  ! time integration

  integer ( kind = 4 ), parameter :: no_of_steps = 1500
  integer ( kind = 4)             :: frequency = 100
  integer (kind = 4 )             :: step 
  real ( kind = 8 )   , parameter :: dt = 0.01
  real ( kind = 8 )               :: start, finish

  ! material specific 

  real ( kind = 8 )   , parameter :: initial_phi = 0.5
  real ( kind = 8 )   , parameter :: mobility = 1.0
  real ( kind = 8 )   , parameter :: grad_coef = 1.0

  ! microstructure

  real ( kind = 8 )   , parameter :: noise = 0.02
  real ( kind = 8 )   , parameter :: A  = 1.0
  integer ( kind = 4 )           :: i, j, jp, jm, ip, im
  real ( kind = 8 )   , dimension ( Nx, Ny ) :: r, phi, dfdphi
  real ( kind = 8 )   , dimension ( Nx, Ny ) :: lap_phi, laplace


  call cpu_time ( start )

  

  ! ===========================================================================
  !                            initial microstucture
  ! ===========================================================================



  phi =  Introduce_fluctuation( initial_phi, noise )



  ! ===========================================================================
  !                         evolution of microstructure 
  ! ===========================================================================


  time_loop: do step = 1, no_of_steps


     column: do j = 1, Ny
        row:do i = 1, Ny


           dfdphi = Deriv_free_energy ( A, phi, i, j )

           laplace = Laplacian ( lap_phi, phi, dx, dy, &
                &  i, j, ip, jp, im, jm )

           phi(i,j) = phi(i,j) - dt*mobility*( dfdphi(i,j) - &
                & grad_coef*laplace(i,j) )


        end do row
     end do column


     ! adjust order parameter in range

     if ( phi(i,j) >= 0.99999 ) phi(i,j) = 0.99999
     if ( phi(i,j) < 0.00001 )  phi(i,j) = 0.00001


     ! print steps on the console

     if ( mod ( step, frequency ) .eq. 0 ) print *, 'Done steps  = ', step


  end do time_loop


  call cpu_time ( finish )



  ! ===========================================================================
  !                                  Output 
  ! ===========================================================================



  call Write_input_parameters_on_file
  call Output_phi_on_file


contains



  ! ===========================================================================
  !                          internal subprograms 
  ! ===========================================================================



  function Introduce_fluctuation ( initial_phi_, noise_ )
    implicit none

    
    real ( kind = 8 ), intent ( in )        :: noise_
    real ( kind = 8 ), intent ( in )        :: initial_phi_
    real ( kind = 8 ), dimension ( Nx, Ny ) :: Introduce_fluctuation
    real ( kind = 8 ), dimension ( Nx, Ny ) :: r

    
    call random_number ( r )

    Introduce_fluctuation = initial_phi + noise_*( 0.5 - r )


  end function Introduce_fluctuation



  ! ---------------------------------------------------------------------------



  pure function Deriv_free_energy ( A_, phi_, i_, j_ )
    implicit none


    real ( kind = 8 ), dimension ( Nx, Ny ) :: Deriv_free_energy
    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( in ) :: phi_
    real ( kind = 8 ), intent ( in )                       :: A_
    integer ( kind = 4 ), intent ( in )                    :: i_, j_


    Deriv_free_energy(i,j) = A_*( 2.0*phi_(i_,j_)*( 1.0 - phi_(i_,j_) )**2 &
         *( 1.0 - 2*phi_(i_,j_) ) )


  end function Deriv_free_energy



  ! ---------------------------------------------------------------------------



  function Laplacian ( lap_phi_, phi_, dx_, dy_, &
       & i_ , j_, ip_, jp_, im_, jm_ )


    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( in out ):: lap_phi_
    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( in ) :: phi_
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


    Laplacian(i_,j_) = ( phi_(ip_,j_) + phi_(im_,j_) + phi_(i_,jm_) + &
         phi_(i_,jp_) - 4.0*phi_(i_,j_) ) / ( dx_*dy_ )              


  end function Laplacian



  ! ---------------------------------------------------------------------------



  subroutine Write_input_parameters_on_file
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
400 format ( 6x, 'time difference', 10x, ' = ', f4.2 )
500 format ( 6x, 'initial phi', 14x, ' = ', f3.1 )
600 format ( 6x, 'mobility', 17x, ' = ', f3.1 )
700 format ( 6x, 'gradient coefficient', 5x, ' = ', f3.1 )
800 format ( 6x, 'noise', 20x' = ', f4.2 )
900 format ( 6x, 'barrier height', 11x' = ', f3.1 )
1000 format( 6x, 'computed time in seconds',1x' = ', f4.2 )


  end subroutine Write_input_parameters_on_file



  ! ---------------------------------------------------------------------------



  subroutine Output_phi_on_file
    implicit none

    character ( len = 80 ) :: filename 

    write ( filename, '( "phi_", i0.3,".dat" )' ) no_of_steps
    open ( 2, file = filename, status = 'replace' )

    do i = 1, Nx    
       write ( 2, * ) ( phi(i,j), j = 1, Ny )
    end do

    close ( 2 )


  end subroutine Output_phi_on_file



end program fd_ac_test

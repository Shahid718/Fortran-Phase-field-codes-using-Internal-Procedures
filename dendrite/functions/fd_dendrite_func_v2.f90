      
!   This program uses internal subprograms for the dendrite model. 
!            
!         
!   Author :
!              Shahid Maqbool
! 
!   Modified :
!               13 Feb. 2023, 10 August 2023
!
!   To compile and run :
!                            check ReadMe file
!
!-------------------------------------------------------------------------------



program fd_Kobayashi_model_test
  implicit none


  ! ============================================================================
  !                                parameters
  ! ============================================================================


  ! simulation cell

  integer ( kind = 4 ), parameter :: Nx = 300
  integer ( kind = 4 ), parameter :: Ny = 300

  ! time integeration

  integer ( kind = 4 ) :: no_of_steps = 2000
  integer ( kind = 4 ) :: frequency = 100
  integer ( kind = 4 ) :: step 
  real ( kind = 8 )    :: start, finish

  ! initial nuclei and evolution 

  real ( kind = 8 )                       :: seed  = 5.0
  real ( kind = 8 ) , dimension( Nx, Ny ) :: phi, tempr
  real ( kind = 8 ) , dimension( Nx, Ny ) :: lap_phi, lap_tempr
  real ( kind = 8 ) , dimension( Nx, Ny ) :: phidx, phidy
  real ( kind = 8 ) , dimension( Nx, Ny ) :: epsil, epsilon_deriv
  real ( kind = 8 )                       :: phi_old, term1, term2
  real ( kind = 8 )                       :: theta, m


  call cpu_time ( start )


  ! ============================================================================
  !                          initial microstucture
  ! ============================================================================



  phi = 0.0
  tempr = 0.0

  phi = Initial_microstructure ( seed )



  ! ============================================================================
  !                       evolution of microstructure 
  ! ============================================================================



  time_loop: do step = 1, no_of_steps


     call Perform_evolution ( phi, tempr, lap_phi, lap_tempr, phidx, &
          & phidy, theta, epsil, epsilon_deriv, phi_old, term1, term2, m )


     ! print steps on the console

     if ( mod ( step, frequency ) .eq. 0 ) print *, 'Done steps  = ', step


  end do time_loop


  call cpu_time ( finish )



  ! ===========================================================================
  !                                  Output 
  ! ===========================================================================



  call Output_files

contains



  ! ============================================================================
  !                            internal subprograms 
  ! ============================================================================



  function Initial_microstructure ( seed_ )
    implicit none

    real ( kind = 8 ), dimension ( Nx, Ny ) :: Initial_microstructure
    real ( kind = 8 ), intent ( in )        :: seed_
    integer ( kind = 4 )                    :: i, j


    do i = 1, Nx
       do j = 1, Ny
          if ( (i - Nx/2.0)*(i - Nx/2.0) + (j - Ny/2.0)*(j - Ny/2.0)&
               & < seed_ ) then
             Initial_microstructure(i,j) = 1.0
          end if
       end do
    end do


  end function Initial_microstructure



  ! ---------------------------------------------------------------------------



  subroutine Perform_evolution ( phi_, tempr_, lap_phi_, lap_tempr_, phidx_, &
       & phidy_, theta_, epsil_, epsilon_deriv_, phi_old_, term1_, term2_, &
       & m_ ) 
    implicit none

    integer ( kind = 4 ) :: i, j, jp, jm, ip, im
    real ( kind = 8 )    :: dx = 0.03, dy = 0.03
    real ( kind = 8 )    :: dtime  = 1.0e-4
    real ( kind = 8 )    :: tau   = 0.0003
    real ( kind = 8 )    :: epsilonb = 0.01
    real ( kind = 8 )    :: kappa = 1.8
    real ( kind = 8 )    :: delta = 0.02
    real ( kind = 8 )    :: aniso = 6.0
    real ( kind = 8 )    :: alpha = 0.9
    real ( kind = 8 )    :: gama  = 10.0
    real ( kind = 8 )    :: teq   = 1.0
    real ( kind = 8 )    :: theta0 = 0.2 
    real ( kind = 8 )    :: pix = 4.0*atan(1.0)

    real ( kind = 8 ), dimension( Nx, Ny ), intent ( in out ) :: phi_, tempr_
    real ( kind = 8 ), dimension( Nx, Ny ), intent ( in out ) :: lap_phi_
    real ( kind = 8 ), dimension( Nx, Ny ), intent ( in out ) :: lap_tempr_
    real ( kind = 8 ), dimension( Nx, Ny ), intent ( in out ) :: phidx_, phidy_
    real ( kind = 8 ), dimension( Nx, Ny ), intent ( in out ) :: epsil_
    real ( kind = 8 ), dimension( Nx, Ny ), intent ( in out ) :: epsilon_deriv_
    real ( kind = 8 )                                         :: term1_, term2_
    real ( kind = 8 ), intent ( in out )                      :: phi_old_
    real ( kind = 8 ), intent ( in out )                      :: theta_, m_


    do concurrent ( j=1:Ny, i=1:Nx )

       jp = j + 1
       jm = j - 1

       ip = i + 1
       im = i - 1

       if ( im == 0 ) im = Nx
       if ( ip == ( Nx + 1) ) ip = 1
       if ( jm == 0 ) jm = Ny
       if ( jp == ( Ny + 1) ) jp = 1

       ! laplacian

       lap_phi_(i,j) = ( phi_(ip,j) + phi_(im,j) + phi_(i,jm) + phi_(i,jp)&
            & - 4.0*phi_(i,j)) / ( dx*dy )
       lap_tempr_(i,j) = ( tempr_(ip,j) + tempr_(im,j) + tempr_(i,jm) + &
            & tempr_(i,jp) - 4.0*tempr_(i,j)) / ( dx*dy )

       ! gradients

       phidx_(i,j) = ( phi_(ip,j) - phi_(im,j) ) / dx
       phidy_(i,j) = ( phi_(i,jp) - phi_(i,jm) ) / dy

       ! angle

       theta_  = atan2( phidy_(i,j), phidx_(i,j) )

       ! epsilon and its derivative

       epsil_(i,j) = epsilonb*( 1.0 + delta*cos(aniso* &
            & ( theta - theta0 ) ) )
       epsilon_deriv_(i,j) = -epsilonb*aniso*delta*sin &
            & ( aniso*( theta - theta0 ) )


    end do

    do concurrent ( j=1:Ny, i=1:Nx )

       
       jp = j + 1
       jm = j - 1

       ip = i + 1
       im = i - 1

       if ( im == 0 ) im = Nx
       if ( ip == ( Nx + 1) ) ip = 1
       if ( jm == 0 ) jm = Ny
       if ( jp == ( Ny + 1) ) jp = 1

       phi_old_ = phi_(i,j)

       ! term1 and term2

       term1_ = ( epsil_(i,jp)*epsilon_deriv_(i,jp)*phidx_(i,jp)&
            & - epsil_(i,jm)*epsilon_deriv_(i,jm)*phidx_(i,jm) ) / dy
       term2_ = -( epsil_(ip,j)*epsilon_deriv_(ip,j)*phidy_(ip,j)&
            & - epsil_(im,j)*epsilon_deriv_(im,j)*phidy_(im,j) ) / dx

       ! factor m

       m_ = alpha/pix*atan( gama*( teq - tempr_(i,j) ) )

       ! time integration

       phi_(i,j) = phi_(i,j) + ( dtime/tau )*( term1_ + term2_ + &
            & epsil_(i,j)**2*lap_phi_(i,j) ) + &
            & phi_old_*( 1.0 - phi_old_ )*( phi_old_ -0.5 + m_ )
       tempr_(i,j) = tempr_(i,j) + dtime*lap_tempr_(i,j) &
            & + kappa*( phi_(i,j) - phi_old_ )

    end do


  end subroutine Perform_evolution



  ! ---------------------------------------------------------------------------



  subroutine output_files
    implicit none

    integer ( kind = 4 ), parameter :: Nx = 300
    integer ( kind = 4 ), parameter :: Ny = 300
    integer ( kind = 4 )   :: no_of_steps = 2000
    real ( kind = 8 )      :: dtime = 1.0e-4
    real ( kind = 8 )      :: seed = 5.0
    real ( kind = 8 )      :: delta = 0.02
    real ( kind = 8 )      :: kappa = 1.8
    real ( kind = 8 )      :: aniso = 6.0
    real ( kind = 8 )      :: gama  = 10.0
    integer ( kind = 4 )   :: i, j
    character ( len = 80 ) :: filename1, filename2 

    open ( unit = 1, file = 'parameters.txt', status = 'replace' )

    write ( 1, * ) ""
    write ( 1, * ) '------------------------------------------------'
    write ( 1, * ) '  Input Parameters of Dendritic Solidification  '
    write ( 1, * ) '------------------------------------------------'
    write ( 1, * ) ""
    write ( 1, 100 ) Nx
    write ( 1, 200 ) Ny
    write ( 1, 300 ) no_of_steps
    write ( 1, 400 ) dtime
    write ( 1, 500 ) seed
    write ( 1, 600 ) delta
    write ( 1, 700 ) kappa
    write ( 1, 800 ) aniso
    write ( 1, 900 ) gama
    write ( 1, 1000 ) finish - start
    write ( 1, * ) ""
    write ( 1, * ) '------------------------------------------------'

    close ( 1 )

100 format ( 6x, 'Nx', 22x, ' = ', i3 )
200 format ( 6x, 'Ny', 22x, ' = ', i3 )
300 format ( 6x, 'number of steps', 9x, ' = ', i4 )
400 format ( 6x, 'time difference', 9x, ' = ', f6.4 )
500 format ( 6x, 'initial seed', 12x, ' = ', f3.1 )
600 format ( 6x, 'delta', 19x, ' = ', f4.2 )
700 format ( 6x, 'gradient coefficient', 4x, ' = ', f3.1 )
800 format ( 6x, 'anisotropy', 14x' = ', f3.1 )
900 format ( 6x, 'gamma', 19x' = ', f4.1 )
1000 format( 6x, 'computed time in seconds = ', f5.2 )


    ! output field files

    write ( filename1, '( "phi_", i0.3,".dat" )' ) no_of_steps
    write ( filename2, '( "temperature_", i0.3,".dat" )' ) no_of_steps

    open ( 2, file = filename1, status = 'replace' )
    open ( 3, file = filename2, status = 'replace' )

    do i = 1, Nx    
       write ( 2, * ) ( phi(i,j), j = 1, Ny )
       write ( 3, * ) ( tempr(i,j), j = 1, Ny )
    end do

    close ( 2 )
    close ( 3 )


  end subroutine output_files



end program fd_Kobayashi_model_test

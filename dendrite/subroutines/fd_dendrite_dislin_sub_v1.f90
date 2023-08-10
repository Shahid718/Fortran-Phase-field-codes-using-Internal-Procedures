
!   This program uses internal subprograms for the dendrite model. 
!            
!         
!   Author :
!              Shahid Maqbool
! 
!   Modified :
!                13 Feb. 2023, 10 August 2023
!
!   To compile and run :
!                            check ReadMe file
!
!-------------------------------------------------------------------------------


program fd_Kobayashi_model_test
  use Dislin
  implicit none


  ! ============================================================================
  !                                parameters
  ! ============================================================================


  ! simulation cell

  integer ( kind = 4 ), parameter :: Nx = 300
  integer ( kind = 4 ), parameter :: Ny = 300
  real ( kind = 8 )               :: dx = 0.03
  real ( kind = 8 )               :: dy = 0.03

  ! time integeration

  integer (kind = 4 ) :: no_of_steps = 2000
  integer (kind = 4 ) :: frequency = 100
  integer (kind = 4 ) :: step 
  real ( kind = 8 )   :: dtime  = 1.0e-4
  real ( kind = 8 )   :: start, finish

  ! material specific 

  real ( kind = 8 )   :: tau   = 0.0003
  real ( kind = 8 )   :: epsilonb = 0.01
  real ( kind = 8 )   :: kappa = 1.8
  real ( kind = 8 )   :: delta = 0.02
  real ( kind = 8 )   :: aniso = 6.0
  real ( kind = 8 )   :: alpha = 0.9
  real ( kind = 8 )   :: gama  = 10.0
  real ( kind = 8 )   :: teq   = 1.0
  real ( kind = 8 )   :: theta0= 0.2 

  real ( kind = 8 )   :: pix   = 4.0*atan(1.0)

  ! initial nuclei and evolution 

  real ( kind = 8 )                       :: seed  = 5.0
  real ( kind = 8 ) , dimension( Nx, Ny ) :: phi, tempr
  real ( kind = 8 ) , dimension( Nx, Ny ) :: lap_phi, lap_tempr
  real ( kind = 8 ) , dimension( Nx, Ny ) :: phidx, phidy
  real ( kind = 8 ) , dimension( Nx, Ny ) :: epsil, epsilon_deriv
  real ( kind = 8 )                       :: phi_old, term1, term2
  real ( kind = 8 )                       :: theta, m
  integer ( kind = 4 )                    :: i, j, ip, im, jp, jm


  call cpu_time ( start )



  ! ============================================================================
  !                          initial microstructure
  ! ============================================================================



  call Initial_microstructure ( phi, tempr, seed, i, j )



  ! ============================================================================
  !                   Setting Initial Dislin routines for animation 
  ! ============================================================================



  call metafl ( 'png' )
  call scrmod ( 'revers' )
  call disini ( )



  ! ============================================================================
  !                       evolution of microstructure 
  ! ============================================================================



  time_loop: do step = 1, no_of_steps


     
     first_spatial_loop:  do concurrent ( j = 1 : Nx , i = 1 : Ny )
        

        call Set_boundary_conditions (i, j, jp, jm, ip, im )

        call Evaluate_laplacian ( phi, tempr, lap_phi, lap_tempr,&
             & dx, dy, i , j, ip, jp, im, jm )

        call Calculate_gradients ( phidx, phidy, phi, dx, dy, i, j, im, ip,&
             & jm, jp )

        call Find_angle ( theta, phidx, phidy, i, j )

        call Compute_epsilon_and_its_derivative ( epsil, epsilonb, delta, &
             & aniso, theta, theta0, epsilon_deriv, i, j )

        
     end do first_spatial_loop

     

     second_spatial_loop:  do concurrent ( j = 1 : Nx , i = 1 : Ny )
        

        call Set_boundary_conditions (i, j, jp, jm, ip, im )

        call Update_phi ( phi_old, phi, i, j )

        call Measure_term1_and_term2 ( term1, term2, epsil, epsilon_deriv,&
             &  phidx, phidy, dy, dx, i, j, ip, im, jp, jm )

        call Determine_factor_m ( m, alpha, pix, gama, teq, tempr, i ,j )

        call Perform_time_integration ( phi, dtime, tau, term1, term2,&
             & epsil, lap_phi, phi_old, m, tempr, lap_tempr, kappa, i, j)

        
     end do second_spatial_loop


     ! print steps on the console

     if ( mod ( step, frequency ) .eq. 0 ) print *, 'Done steps  = ', step

     ! save multiplot

     call Dislin_color_multi_plot ( )


  end do time_loop


  call cpu_time ( finish )


  ! ===========================================================================
  !                                  Output 
  ! ===========================================================================



  call Write_input_parameters_on_file
  call Output_phi_and_temperature_on_files


contains



  ! ============================================================================
  !                            internal subprograms 
  ! ============================================================================



  subroutine Initial_microstructure ( phi_, tempr_, seed_, i_, j_ )
    implicit none


    real ( kind = 8 ), intent ( in )                        :: seed_
    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( out ) :: phi_, tempr_
    integer ( kind = 4 ), intent (in out )                  :: i_, j_


    phi_ = 0.0
    tempr_ = 0.0

    do i_ = 1, Nx
       do j_ = 1, Ny
          if ( (i_ - Nx/2.0)*(i_ - Nx/2.0) + (j_ - Ny/2.0)*(j_ - Ny/2.0)&
               & < seed_ ) then
             phi_(i_,j_) = 1.0
          end if
       end do
    end do

  end subroutine Initial_microstructure



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



  pure subroutine Evaluate_laplacian ( phi_, tempr_, lap_phi_, lap_tempr_, &
       & dx_, dy_, i_ , j_, ip_, jp_, im_, jm_ )
    implicit none


    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( in ) :: phi_,tempr_
    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( out ):: lap_phi_
    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( out ):: lap_tempr_
    real ( kind = 8 ), intent ( in )                       :: dx_, dy_
    integer ( kind = 4 ), intent ( in )                    :: i_, j_, ip_, jp_
    integer ( kind = 4 ), intent ( in )                    :: im_, jm_


    lap_phi_(i_,j_) = ( phi_(ip_,j_) + phi_(im_,j_) + phi_(i_,jm_) &
         + phi_(i_,jp_) - 4.0*phi_(i_,j_)) / ( dx_*dy_ )

    lap_tempr_(i_,j_) = ( tempr_(ip_,j_) + tempr_(im_,j_) + tempr_(i_,jm_) + &
         & tempr_(i_,jp_) - 4.0*tempr_(i_,j_)) / ( dx_*dy_ )


  end subroutine Evaluate_laplacian



  ! ---------------------------------------------------------------------------



  pure subroutine Calculate_gradients ( phidx_, phidy_, phi_, dx_, &
       & dy_, i_, j_, im_, ip_, jm_, jp_ )
    implicit none


    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( in ) :: phi_
    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( out ):: phidx_, phidy_
    real ( kind = 8 ), intent ( in )                       :: dx_, dy_
    integer ( kind = 4 ), intent ( in )                    :: i_, j_, ip_, jp_
    integer ( kind = 4 ), intent ( in )                    :: im_, jm_


    phidx_(i_,j_) = ( phi_(ip_,j_) - phi_(im_,j_) ) / dx_
    phidy_(i_,j_) = ( phi_(i_,jp_) - phi_(i_,jm_) ) / dy_


  end subroutine Calculate_gradients



  ! ---------------------------------------------------------------------------



  pure subroutine Find_angle ( theta_, phidx_, phidy_, i_, j_ )
    implicit none


    real ( kind = 8 )                      , intent ( out ) :: theta_
    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( out ) :: phidx_, phidy_
    integer ( kind = 4 ), intent ( in )                     :: i_, j_


    theta_ = atan2( phidy_(i_,j_),phidx_(i_,j_) )


  end subroutine Find_angle



  ! ---------------------------------------------------------------------------



  pure subroutine Compute_epsilon_and_its_derivative ( epsil_, epsilonb_, & 
       & delta_, aniso_, theta_, theta0_, epsilon_deriv_, i_, j_ )
    implicit none


    real ( kind = 8 ), intent ( in )                        :: delta_, aniso_
    real ( kind = 8 ), intent ( in )                        :: theta_, theta0_
    real ( kind = 8 ), intent ( in )                        :: epsilonb_
    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( out ) :: epsil_
    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( out ) :: epsilon_deriv_
    integer ( kind = 4 ), intent ( in )                     :: i_, j_


    epsil_(i_,j_) = epsilonb_*( 1.0 + delta_*cos(aniso_*&
         & ( theta_ - theta0_ ) ) )
    epsilon_deriv_(i_,j_) = -epsilonb_*aniso_*delta_*sin&
         & ( aniso_*( theta_ - theta0_ ) )


  end subroutine Compute_epsilon_and_its_derivative



  ! ---------------------------------------------------------------------------



  pure subroutine Update_phi ( phi_old_, phi_, i_, j_ )
    implicit none


    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( in ) :: phi_
    real ( kind = 8 ), intent ( out )                      :: phi_old_
    integer ( kind = 4 ), intent ( in )                    :: i_, j_


    phi_old_ = phi_(i_,j_)


  end subroutine Update_phi



  ! ---------------------------------------------------------------------------



  pure subroutine Measure_term1_and_term2 ( term1_, term2_, epsil_, &
       & epsilon_deriv_, phidx_, phidy_, dy_, dx_, i_, j_, ip_, im_, jp_, jm_ )
    implicit none


    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( in )  :: epsil_
    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( in )  :: epsilon_deriv_
    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( in )  :: phidx_, phidy_
    real ( kind = 8 ), intent ( out )                       :: term1_, term2_
    integer ( kind = 4 ), intent ( in )                     :: i_, j_, ip_, im_
    integer ( kind = 4 ), intent ( in )                     :: jp_, jm_
    real ( kind = 8 ), intent ( in )                        :: dx_, dy_


    term1_ = ( epsil_(i_,jp_)*epsilon_deriv_(i_,jp_)*phidx_(i_,jp_)&
         & - epsil_(i_,jm_)*epsilon_deriv_(i_,jm_)*phidx_(i_,jm_) ) / dy_

    term2_ = -( epsil_(ip_,j_)*epsilon_deriv_(ip_,j_)*phidy_(ip_,j_)&
         & - epsil_(im_,j_)*epsilon_deriv_(im_,j_)*phidy_(im_,j_) ) / dx_


  end subroutine Measure_term1_and_term2



  ! ---------------------------------------------------------------------------



  pure subroutine determine_factor_m ( m_, alpha_, pix_, gama_, teq_, &
       & tempr_, i_, j_ )
    implicit none


    real ( kind = 8 ), intent ( in )                       :: alpha_, pix_ 
    real ( kind = 8 ), intent ( in )                       :: gama_, teq_
    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( in ) :: tempr_
    real ( kind = 8 ), intent ( out )                      :: m_
    integer ( kind = 4 ), intent ( in )                    :: i_, j_


    m_ = alpha_/pix_*atan( gama_*( teq_ - tempr_(i_,j_) ) )


  end subroutine determine_factor_m



  ! ---------------------------------------------------------------------------



  pure subroutine Perform_time_integration ( phi_, dtime_, tau_, term1_, term2_&
       & ,epsil_, lap_phi_, phi_old_, m_, tempr_, lap_tempr_, kappa_, i_, j_)
    implicit none


    real ( kind = 8 ), intent ( in )                         :: dtime_, tau_
    real ( kind = 8 ), intent ( in )                         :: kappa_
    real ( kind = 8 ), intent ( in )                         :: term1_, term2_
    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( in )   :: epsil_, lap_phi_
    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( in )   :: lap_tempr_
    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( inout ):: phi_, tempr_
    real ( kind = 8 ), intent ( in )                         :: m_, phi_old_
    integer ( kind = 4 ), intent ( in )                      :: i_, j_


    phi_(i_,j_) = phi_(i_,j_) + ( dtime_/tau_ )*( term1_ + term2_ +&
         & epsil_(i_,j_)**2*lap_phi_(i_,j_) ) + &
         & phi_old_*( 1.0 - phi_old_ )*( phi_old_ - 0.5 + m_ )

    tempr_(i_,j_) = tempr_(i_,j_) + dtime_*lap_tempr_(i_,j_) &
         & + kappa_*( phi(i_,j_) - phi_old_ )


  end subroutine Perform_time_integration



  ! ---------------------------------------------------------------------------



  subroutine Write_input_parameters_on_file
    open ( unit = 1, file = 'input_parameters.txt', status = 'replace' )

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
400 format ( 6x, 'time difference', 9x, ' = ', e7.1 )
500 format ( 6x, 'initial seed', 12x, ' = ', f3.1 )
600 format ( 6x, 'delta', 19x, ' = ', e7.1 )
700 format ( 6x, 'gradient coefficient', 4x, ' = ', f3.1 )
800 format ( 6x, 'anisotropy', 14x' = ', f3.1 )
900 format ( 6x, 'gamma', 19x' = ', f4.1 )
1000 format( 6x, 'computed time in seconds = ', f5.2 )


  end subroutine Write_input_parameters_on_file



  ! ---------------------------------------------------------------------------



  subroutine Output_phi_and_temperature_on_files
    implicit none


    character ( len = 80 ) :: filename1, filename2 

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


  end subroutine Output_phi_and_temperature_on_files



  ! ---------------------------------------------------------------------------



  subroutine Dislin_color_multi_plot ( )

    call complx ( )
    call intax ( )
    call labdig ( 2, 'Z' )

    if ( step == 1 ) then

       call color ( 'Blue')
       call messag ( 'phi', 650, 100 )
       call color ( 'Fore' )
       call axspos ( 400, 800 )
       call ax3len ( 600, 600, 600 )
       call graf3 ( 0.d0, 300.d0, 0.d0, 100.d0, 0.d0, 300.d0,&
            & 0.d0, 100.d0, 0.0d0, 1.2d0, 0.0d0, 0.2d0 )
       call crvmat ( phi, Nx, Ny, 1, 1 )
       call endgrf
       call messag ( 'step = 1', 600, 950 )

       call color ( 'Blue')
       call messag ( 'temperature', 1800, 100 )
       call color ( 'Fore' )
       call axspos ( 1700, 800 )
       call ax3len ( 600, 600, 600 )
       call graf3 ( 0.d0, 300.d0, 0.d0, 100.d0, 0.d0, 300.d0,&
            & 0.d0, 100.d0, 0.0d0, 1.2d0, 0.0d0, 0.2d0 )
       call crvmat ( tempr, Nx, Ny, 1, 1 )
       call endgrf
       call messag ( 'step = 1', 1900, 950 )

    else if ( step == 2000 ) then

       call axspos ( 400, 1700 )
       call ax3len ( 600, 600, 600 )
       call graf3 ( 0.d0, 300.d0, 0.d0, 100.d0, 0.d0, 300.d0,&
            & 0.d0, 100.d0, 0.0d0, 1.2d0, 0.0d0, 0.2d0 )
       call crvmat ( phi, Nx, Ny, 1, 1 )
       call endgrf
       call messag ( 'step = 2000', 550, 1850 )

       call axspos ( 1700, 1700 )
       call ax3len ( 600, 600, 600 )
       call graf3 ( 0.d0, 300.d0, 0.d0, 100.d0, 0.d0, 300.d0,&
            & 0.d0, 100.d0, 0.0d0, 1.2d0, 0.0d0, 0.2d0 )
       call crvmat ( tempr, Nx, Ny, 1, 1 )
       call endgrf
       call messag ( 'step = 2000', 1850, 1850 )

       call disfin ( )

    end if


  end subroutine Dislin_color_multi_plot




end program fd_Kobayashi_model_test

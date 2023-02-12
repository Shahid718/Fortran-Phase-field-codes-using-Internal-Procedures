
!   This program uses internal subprograms for the dendrite model. 
!            
!         
!   Author :
!              Shahid Maqbool
! 
!   Modified :
!               13 Feb. 2023
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
  real ( kind = 8 )   :: mu    = 1.0
  real ( kind = 8 )   :: kappa = 1.8
  real ( kind = 8 )   :: delta = 0.02
  real ( kind = 8 )   :: aniso = 6.0
  real ( kind = 8 )   :: alpha = 0.9
  real ( kind = 8 )   :: gama  = 10.0
  real ( kind = 8 )   :: teq   = 1.0
  real ( kind = 8 )   :: theta0= 0.2 

  real ( kind = 8 )   :: pix   = 4.0*atan(1.0)

  ! initial nuclei and evolution 

  real ( kind = 8 )                      :: seed  = 5.0
  real ( kind = 8 ), dimension( Nx, Ny ) :: phi, tempr
  real ( kind = 8 ), dimension( Nx, Ny ) :: lap_phi, lap_tempr
  real ( kind = 8 ), dimension( Nx, Ny ) :: phidx, phidy
  real ( kind = 8 ), dimension( Nx, Ny ) :: epsil, epsilon_deriv
  real ( kind = 8 )                      :: phi_old, term1, term2
  real ( kind = 8 )                      :: theta, m
  integer ( kind = 4 )                   :: i, j, istep, ip, im, jp, jm


  call cpu_time ( start )



  ! ============================================================================
  !                          initial microstucture
  ! ============================================================================



  tempr = 0.0
  phi = 0.0

  phi = Initial_microstructure ( seed )



  ! ============================================================================
  !                   Setting Initial Dislin routines for multiplot 
  ! ============================================================================



  call metafl ( 'png' )
  call scrmod ( 'revers' )
  call disini ( )



  ! ============================================================================
  !                       evolution of microstructure 
  ! ============================================================================

  

  time_loop: do step = 1, no_of_steps


     do concurrent ( i = 1:Nx, j = 1:Ny )

        ! boundary conditions

        jp = j + 1
        jm = j - 1

        ip = i + 1
        im = i - 1

        if ( im == 0 ) im = Nx
        if ( ip == ( Nx + 1 ) ) ip = 1
        if ( jm == 0 ) jm = Ny
        if ( jp == ( Ny + 1 ) ) jp = 1

        ! invoking laplacian function

        lap_phi = Laplacian ( phi, dx, dy, i , j, ip, jp, im, jm )
        lap_tempr = Laplacian (tempr, dx, dy, i , j, ip, jp, im, jm )

        ! invoking gradient function

        phidx = Gradient_x ( phi, dx, i, j, im, ip, jm, jp )
        phidy = Gradient_y ( phi, dy, i, j, im, ip, jm, jp )

        ! invoking angle function

        theta = Theta_angle ( phidy, phidx, i, j )

        ! epsilon and its derivative

        epsil(i,j) = epsilonb*( 1.0 + delta*cos(aniso*&
             & ( theta - theta0 ) ) )
        epsilon_deriv(i,j) = -epsilonb*aniso*delta*sin&
             & ( aniso*( theta - theta0 ) )

     end do

     do concurrent ( i = 1:Nx, j = 1:Ny )

        ! boundary conditions

        jp = j + 1
        jm = j - 1

        ip = i + 1
        im = i - 1

        if ( im == 0 ) im = Nx
        if ( ip == ( Nx + 1 ) ) ip = 1
        if ( jm == 0 ) jm = Ny
        if ( jp == ( Ny + 1 ) ) jp = 1

        ! update phi

        phi_old = phi(i,j)

        ! term1 and term2

        term1 = ( epsil(i,jp)*epsilon_deriv(i,jp)*phidx(i,jp) &
             & - epsil(i,jm)*epsilon_deriv(i,jm)*phidx(i,jm) ) / dy

        term2 = - ( epsil(ip,j)*epsilon_deriv(ip,j)*phidy(ip,j) &
             & - epsil(im,j)*epsilon_deriv(im,j)*phidy(im,j) ) / dx

        ! invoking function 

        m = Factor_m ( alpha, pix, gama, teq, tempr, i, j )

        ! time integration

        phi(i,j) = phi(i,j) + ( dtime/tau )*( term1 + term2 + &
             & epsil(i,j)**2*lap_phi(i,j) ) + &
             & phi_old*( 1.0 - phi_old )*( phi_old - 0.5 + m )

        tempr(i,j) = tempr(i,j) + dtime*lap_tempr(i,j) &
             & + kappa*( phi(i,j) - phi_old )

     end do


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



  pure function Laplacian ( Order_parameter, dx_, dy_, i_, j_, ip_, jp_, im_, jm_ )
    implicit none

    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( in ) :: Order_parameter
    real ( kind = 8 ), dimension ( Nx, Ny ) :: laplacian
    real ( kind = 8 ), intent ( in )   :: dx_, dy_
    integer ( kind = 4 ),intent ( in ) :: i_, j_
    integer ( kind = 4 ),intent ( in ) :: jp_, jm_, ip_, im_


    Laplacian(i_,j_) = ( Order_parameter(ip_,j_) + Order_parameter(im_,j_) + &
         & Order_parameter(i_,jm_) +  Order_parameter(i_,jp_) - &
         & 4.0* Order_parameter(i_,j_) ) / ( dx_*dy_ )


  end function Laplacian



  ! ---------------------------------------------------------------------------


  pure function Gradient_x ( phi_, dx_, i_, j_, im_, ip_, jm_, jp_ )
    implicit none


    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( in ) :: phi_
    real ( kind = 8 ), dimension ( Nx, Ny ):: gradient_x
    real ( kind = 8 ), intent ( in )       :: dx_
    integer ( kind = 4 ), intent ( in )    :: i_, j_, ip_, im_, jp_, jm_


    Gradient_x(i_,j_) = ( phi_(ip_,j_) - phi_(im_,j_) ) / dx_


  end function Gradient_x


  pure function Gradient_y ( phi_, dy_, i_, j_, im_, ip_, jm_, jp_ )
    implicit none


    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( in ) :: phi_
    real ( kind = 8 ), dimension ( Nx, Ny ):: gradient_y
    real ( kind = 8 ), intent ( in )       :: dy_
    integer ( kind = 4 ), intent ( in )    :: i_, j_, ip_, im_, jp_, jm_


    Gradient_y(i_,j_) = ( phi_(i_,jp_) - phi_(i_,jm_) ) / dy_


  end function Gradient_y



  ! ---------------------------------------------------------------------------



  pure function Theta_angle (  phidy_, phidx_, i_, j_ )
    implicit none


    real ( kind = 8 ) :: theta_angle
    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( in ) :: phidx_, phidy_
    integer ( kind = 4 ), intent ( in )                     :: i_, j_


    Theta_angle = atan2( phidy_(i_,j_), phidx_(i_,j_) )


  end function Theta_angle



  ! ---------------------------------------------------------------------------



  pure function Factor_m ( alpha_, pix_, gama_, teq_, tempr_, i_, j_ )
    implicit none


    real ( kind = 8 )                   :: factor_m
    real ( kind = 8 ), intent ( in )    :: alpha_, pix_ 
    real ( kind = 8 ), intent ( in )    :: gama_, teq_
    integer ( kind = 4 ), intent ( in ) :: i_, j_
    real ( kind = 8 ), dimension ( Nx, Ny ), intent ( in ) :: tempr_


    Factor_m  = alpha_ / pix_*atan( gama_*( teq_ - tempr_(i_,j_) ) )


  end function Factor_m



  ! ---------------------------------------------------------------------------



  subroutine Write_input_parameters_on_file
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

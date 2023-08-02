!=====================================================================================================
!    Siwei Zhu
!    China University of Gesciences, Wuhan, China
!    selwey1996@gmail.com
!
!    The code is a parallel Fortran program targeted at the 1D forward modeling specifically designed
!    for overhead transmission lines.
!
!    The code is a free program. You can redistribute it and/or modify it under the terms of the GNU
!    General Public License.
!
!    Feel free to contact the author for any kind of problem related to this code.
!=====================================================================================================

module parameter_set

!============================================================================================
!This file sets parameters used during forward modeling.
!============================================================================================
    
    implicit none

!---------------------------------------------
! Global constants
!---------------------------------------------

    real(8), parameter :: pi = 3.141592653589793d0
    real(8), parameter :: mu0 = 4d-7*pi                       ! The permeability of air
    real(8), parameter :: eps0 = 1d0/(mu0*299792458.0d0**2)   ! The permittivity of air
    real(8), parameter :: sig0 = 1d-8

!---------------------------------------------
! Power-line parameters
!---------------------------------------------
    ! Geometric parameters of overhead power lines
    integer, public :: num_circuit                         ! Number of circuits, each circuit contains three phase concuctors
    integer, public :: num_span                            ! Number of spans
    integer, public, allocatable :: bundle_num(:)          ! Bundle number per phase conductor 
    real(8), public, allocatable :: x_coor_suspension(:,:) ! X coordinates of suspension points
    real(8), public, allocatable :: y_coor_suspension(:,:) ! Y coordinates of suspension points
    real(8), public, allocatable :: z_coor_suspension(:,:) ! Z coordinates of suspension points
    real(8), public, allocatable :: q(:,:)                 ! Catenary coefficient
    ! Current loaded on overhead power lines
    real(8), public :: freq_I                  ! Frequency
    real(8), public, allocatable :: amp_I(:,:) ! Amplitude
    real(8), public, allocatable :: ph_I(:,:)  ! Initial phase

!---------------------------------------------
! Receiver parameters
!---------------------------------------------
    integer, public :: num_o    ! Number of observation points
    real(8), public :: xos, xoe ! X cordinates of the starting and ending point for observation, and point spacing
    real(8), public :: yos, yoe ! Y cordinates of the starting and ending point for observation, and point spacing
    real(8), public :: zos, zoe ! Z cordinates of the starting and ending point for observation, and point spacing

!---------------------------------------------
! 1D model parameters
!---------------------------------------------
    integer, public :: num_layer                                    ! Number of model layers
    real(8), public, allocatable :: sig_layer(:), eps_layer(:)      ! Conductivity and permittivity of layerd model
    real(8), public, allocatable :: mu_layer(:), thickness_layer(:) ! Permeability and thickness of layerd model

!---------------------------------------------
! Parameter for Gauss quadrature integration
! and Hankel transform
!---------------------------------------------
    integer, public :: numIntegPts      ! Number of points to use for Gauss quadrature integration
    character(32), public :: HTmethod1D ! Filter options listed in FilterModules.f90
!---------------------------------------------
! Contained subroutines:
!---------------------------------------------
    contains
!==============================================================================!
!======================================================== power_line_parameter !
!==============================================================================!

    subroutine power_line_parameter

    ! Geometric parameters of overhead power lines

        num_circuit = 1
        num_span = 5
        allocate(bundle_num(num_circuit))
        bundle_num = [1]

        ! Value in row i and column j: x/y/z coordinate of the ith suspension point for the j bundled conductor
        allocate(x_coor_suspension(num_span+1,3*sum(bundle_num)))
        allocate(y_coor_suspension(num_span+1,3*sum(bundle_num)))
        allocate(z_coor_suspension(num_span+1,3*sum(bundle_num)))
        x_coor_suspension = reshape([-10.0,-10.0,-10.0,-10.0,-10.0,-10.0,&
                                     0.0,0.0,0.0,0.0,0.0,0.0,&
                                     10.0,10.0,10.0,10.0,10.0,10.0]&
                                    ,(/num_span+1,3*sum(bundle_num)/))
        y_coor_suspension = reshape([-1250.0,-750.0,-250.0,250.0,750.0,1250.0,&
                                     -1250.0,-750.0,-250.0,250.0,750.0,1250.0,&
                                     -1250.0,-750.0,-250.0,250.0,750.0,1250.0]&
                                    ,(/num_span+1,3*sum(bundle_num)/))
        z_coor_suspension = reshape([-30.0,-30.0,-30.0,-30.0,-30.0,-30.0,&
                                     -30.0,-30.0,-30.0,-30.0,-30.0,-30.0,&
                                     -30.0,-30.0,-30.0,-30.0,-30.0,-30.0]&
                                    ,(/num_span+1,3*sum(bundle_num)/))

        
        ! Value in row i and column j: Catenary coefficient for the ith catenary for the jth phase conductor
        allocate(q(num_span,3*num_circuit))
        
        q = reshape([1.0/3.198295e-4,1.0/3.198295e-4,1.0/3.198295e-4,1.0/3.198295e-4,1.0/3.198295e-4,&
                     1.0/3.198295e-4,1.0/3.198295e-4,1.0/3.198295e-4,1.0/3.198295e-4,1.0/3.198295e-4,&
                     1.0/3.198295e-4,1.0/3.198295e-4,1.0/3.198295e-4,1.0/3.198295e-4,1.0/3.198295e-4]&
                     ,(/num_span,3*num_circuit/))

        
    ! Current loaded on overhead power lines
        freq_I = 50.0
        allocate(amp_I(3,num_circuit))
        allocate(ph_I(3,num_circuit))
        ! Value in row i and column j: Amplitude of loading current on each bundled conductor
        ! of the ith phase conductor in the jth circuit
        amp_I = reshape([1.0, 1.0, 1.0]&
                         ,(/3,num_circuit/))
        ! Value in row i and column j: Initial phase of loading current on each bundled conductor
        ! of the ith phase conductor in the jth circuit                 
        ph_I = reshape([-120.0, 0.0, 120.0]&
                        ,(/3,num_circuit/))

    end subroutine power_line_parameter

!==============================================================================!
!========================================================== receiver_parameter !
!==============================================================================!

    subroutine receiver_parameter

        num_o = 101

        xos = -1000.0
        xoe = 1000.0

        yos = -100.0
        yoe = -100.0

        zos = 0.0
        zoe = 0.0

    end subroutine receiver_parameter

!==============================================================================!
!============================================================= model_parameter !
!==============================================================================!

    subroutine model_parameter

        ! Number of model layers (excluding air layer)
        num_layer = 3
        
        allocate(sig_layer(num_layer),mu_layer(num_layer))
        allocate(eps_layer(num_layer),thickness_layer(num_layer-1))

        ! Conductivity of each layer (excluding air layer)
        sig_layer = [0.01,0.1,0.001]
        ! Permeability of each layer (excluding air layer), usually set as the value of air layer
        mu_layer = [mu0,mu0,mu0]
        ! Permittivity of each layer (excluding air layer), usually set as the value of air layer
        eps_layer = [eps0,eps0,eps0]

        ! Thickness of each layer (excluding air layer and the last layer)
        ! Please ensure that the array thickness_layer is not empty.
        ! If the number of model layers is less than 2, you could set a number arbitrarily.
        thickness_layer = [200.0,100.0]

    end subroutine model_parameter

!==============================================================================!
!============================================================= GS_HK_parameter !
!==============================================================================!

    subroutine GS_HK_parameter

        numIntegPts = 100        ! The number of points must be greater than 0

        HTmethod1D = 'kk_ht_201' ! Use 201 point HT digital filters.

    end subroutine GS_HK_parameter

end module parameter_set
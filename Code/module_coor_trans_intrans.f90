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

module coor_trans_intrans
!============================================================================================
!This file is used for coordinate translation and inverse coordinate translation.
!============================================================================================

    contains

!==============================================================================!
!====================================================== Coordinate translation !
!==============================================================================!     
    subroutine coor_trans(P_start,P_end,q_coef,xo,yo,zo,alpha,h1,h2,l,a,xo_tr,yo_tr,zo_tr)
        implicit none
        real(8), intent(in) :: P_start(:), P_end(:), q_coef
        real(8), intent(in) :: xo(:), yo(:), zo(:)
        real(8), intent(inout) :: alpha, h1, h2, l, a
        real(8), intent(inout) :: xo_tr(:), yo_tr(:), zo_tr(:)

        real(8) :: start_1(3), end_1(3)
        real(8), allocatable :: xo_1(:), yo_1(:), zo_1(:)

        ! Coordinate translation
        ! The starting point will be on the z axis after transformation
        start_1 = [P_start(1:2)-P_start(1:2),P_start(3)] ! Starting point of the catenary
        end_1 = [P_end(1:2)-P_start(1:2),P_end(3)]       ! Ending point of the catenary
        allocate(xo_1(size(xo)),yo_1(size(yo)),zo_1(size(zo)))
        xo_1 = xo-P_start(1)                             ! x coordinates of observation points
        yo_1 = yo-P_start(2)                             ! y coordinates of observation points
        zo_1 = zo                                        ! z coordinates of observation points

        ! Coordinate rotation
        ! The catenary will be in the xOz plane after transformation
        alpha = -atan2(end_1(2),end_1(1))           ! Rotation angle
        h1 = start_1(3)                             ! Z coordinate of starting point after coordinate transformation
        h2 = end_1(3)                               ! Z coordinate of ending point after coordinate transformation
        l = end_1(1)*cos(alpha)-end_1(2)*sin(alpha) ! X coordinate of ending point after coordinate transformation   
		a = l/2+q_coef*asinh((h2-h1)/2/q_coef/sinh(l/2/q_coef)) ! X coordinate of the lowest point on the catenary
        xo_tr = xo_1*cos(alpha)-yo_1*sin(alpha)     ! X coordinates of observation points after coordinate transformation
		yo_tr = yo_1*cos(alpha)+xo_1*sin(alpha)     ! Y coordinates of observation points after coordinate transformation
		zo_tr = zo_1                                ! Z coordinates of observation points after coordinate transformation
        
    end subroutine coor_trans

!==============================================================================!
!================================= Inverse Coordinate translation of responses !
!==============================================================================!         
    subroutine in_trans(Ex_in,Ey_in,Ez_in,Hx_in,Hy_in,Hz_in,alpha,Ex_out,Ey_out,Ez_out,Hx_out,Hy_out,Hz_out)

        complex(8), intent(in) :: Ex_in(:), Ey_in(:), Ez_in(:), Hx_in(:), Hy_in(:), Hz_in(:)
        real(8), intent(in) :: alpha
        complex(8), intent(inout) :: Ex_out(:), Ey_out(:), Ez_out(:), Hx_out(:), Hy_out(:), Hz_out(:)

        Ex_out = Ex_in*cos(-alpha)-Ey_in*sin(-alpha)     ! Ex after inverse coordinate transformation
		Ey_out = Ey_in*cos(-alpha)+Ex_in*sin(-alpha)     ! Ey after inverse coordinate transformation
        Hx_out = Hx_in*cos(-alpha)-Hy_in*sin(-alpha)     ! Hx after inverse coordinate transformation
		Hy_out = Hy_in*cos(-alpha)+Hx_in*sin(-alpha)     ! Hy after inverse coordinate transformation
        Ez_out = Ez_in
        Hz_out = Hz_in

    end subroutine in_trans

!==============================================================================!
!====================================================== Compute function asinh !
!==============================================================================!         
    function asinh(x)
        real(8), intent(in) :: x
        real(8) :: asinh
        asinh = log(x+sqrt(x**2+1))
    end function


end module coor_trans_intrans
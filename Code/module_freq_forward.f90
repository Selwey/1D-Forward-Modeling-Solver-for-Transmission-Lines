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

module freq_forward
!============================================================================================
!This file culculates forward responses.
!============================================================================================

    use parameter_set
    use GS_Hk
    !$ use omp_lib

    implicit none

    contains

    subroutine freq_sounding(h1, h2, l, q_coef, a, I_current, w, xo_tr, yo_tr, zo_tr,&
                             Ex_t, Ey_t, Ez_t, Hx_t, Hy_t, Hz_t)

        implicit none
        
        real(8), intent(in) :: h1, h2, l, q_coef, a, w
        complex(8), intent(in) :: I_current
        real(8), intent(in) :: xo_tr(:), yo_tr(:), zo_tr(:)
        complex(8), intent(inout) :: Ex_t(:), Ey_t(:), Ez_t(:), Hx_t(:), Hy_t(:), Hz_t(:)

        real(8), allocatable :: gx(:), gw(:), gz(:)
        integer :: i, j

        real(8) :: X_G, Y, P_G
        
        complex(8) :: impe0, admi0 ! Impedance and admittance of air
        real(8), allocatable :: lam_G(:)
        complex(8) :: k0
        complex(8), allocatable :: u0_G(:)
        complex(8), allocatable :: rm_G(:), re_G(:)
        complex(8), allocatable :: eZ_nega_G(:), eZ_posi_G(:)

        complex(8) :: Ex_h1, Ex_h2, Ex_h3, Ex_h4
        complex(8) :: Ey_h1, Ey_h2, Ey_h3
        complex(8) :: Ez_h1, Ez_h2
        complex(8) :: Hx_h1, Hx_h2, Hx_h3
        complex(8) :: Hy_h1, Hy_h2, Hy_h3, Hy_h4
        complex(8) :: Hz_h1
        complex(8) :: Ex_G, Ey_G, Ez_G, Hx_G, Hy_G, Hz_G

        ! Impedance and admittance of air
        impe0 = dcmplx(0.d0,mu0*w)
        admi0 = dcmplx(sig0,eps0*w)
        k0 = sqrt(-impe0*admi0)
        
        ! Move Gauss integral weights and integral points from [-1,1] to [0,l]
        allocate(gx(size(xint)),gw(size(weights)),gz(size(xint)))
        gx = xint*l/2.d0+l/2.d0 
        gw = weights*l/2.d0
        gz = q_coef*(cosh(a/q_coef)-cosh((gx-a)/q_coef))+h1

        ! Allocate lamda and u0 used for Hankel transformation
        allocate(lam_G(ndhtfc),u0_G(ndhtfc))

        ! Allocate re and rm used for Hankel transformation
        allocate(rm_G(size(base)),re_G(size(base)))

        ! Allocate eZ_nega_G and eZ_posi_G used for Hankel transformation
        allocate(eZ_nega_G(size(base)),eZ_posi_G(size(base)))
        
        !$OMP PARALLEL DO private(Y, Ex_G, Ey_G, Ez_G, Hx_G, Hy_G, Hz_G)
        do i = 1, num_o ! Loop over measuring points

            Y = yo_tr(i)-0.d0

            Ex_G = (0.d0,0.d0)
            Ey_G = (0.d0,0.d0)
            Ez_G = (0.d0,0.d0)
            Hx_G = (0.d0,0.d0)
            Hy_G = (0.d0,0.d0)
            Hz_G = (0.d0,0.d0)
            
            
            !$OMP  PARALLEL DO private(X_G, P_G, lam_G, u0_G, eZ_posi_G, eZ_nega_G, rm_G, re_G, Ex_h1, &
            !$OMP                      Ex_h2, Ex_h3, Ex_h4, Ey_h1, Ey_h2, Ey_h3, Ez_h1, Ez_h2, Hx_h1, Hx_h2, &
            !$OMP                      Hx_h3, Hy_h1, Hy_h2, Hy_h3, Hy_h4, Hz_h1) &
            !$OMP              reduction(+:Ex_G,Ey_G,Ez_G,Hx_G,Hy_G,Hz_G)
            do j = 1, size(gx) ! Loop over Gauss points
            
                X_G = xo_tr(i)-gx(j)
                P_G = sqrt(X_G**2+Y**2)

                lam_G = 1.d0/P_G*base
                u0_G = sqrt(lam_G**2-k0**2)

                eZ_posi_G = exp(u0_G*(zo_tr(i)+gz(j)))
                eZ_nega_G = exp(-u0_G*abs(zo_tr(i)-gz(j)))

                ! Reflection coefficient
                rm_G = RTM(num_layer,sig_layer,mu_layer,eps_layer,thickness_layer,u0_G,admi0,lam_G,w)
                re_G = RTE(num_layer,sig_layer,mu_layer,eps_layer,thickness_layer,u0_G,impe0,lam_G,w)
            
                ! Compute the Hankel integral items

                Ex_h1 = 1.d0/P_G*sum(((impe0/u0_G-u0_G/admi0)*eZ_nega_G+&
                                      (u0_G/admi0*rm_G+impe0/u0_G*re_G)*eZ_posi_G)*htj1)
                Ex_h2 = 1.d0/P_G*sum(lam_G*((impe0/u0_G-u0_G/admi0)*eZ_nega_G+&
                                            (u0_G/admi0*rm_G+impe0/u0_G*re_G)*eZ_posi_G)*htj0)
                Ex_h3 = 1.d0/P_G*sum(lam_G/u0_G*(eZ_nega_G+re_G*eZ_posi_G)*htj0)
                Ex_h4 = 1.d0/P_G*sum(lam_G**2*(-eZ_nega_G+rm_G*eZ_posi_G)*htj1)

                Ey_h1 = Ex_h1
                Ey_h2 = Ex_h2
                Ey_h3 = Ex_h4

                Ez_h1 = 1.d0/P_G*sum(lam_G**2*(-eZ_nega_G-rm_G*eZ_posi_G)*htj1)
                Ez_h2 = 1.d0/P_G*sum(lam_G**3/u0_G*(eZ_nega_G+rm_G*eZ_posi_G)*htj0)

                Hx_h1 = 1.d0/P_G*sum((rm_G+re_G)*eZ_posi_G*htj1)
                Hx_h2 = 1.d0/P_G*sum(lam_G*(rm_G+re_G)*eZ_posi_G*htj0)
                Hx_h3 = 1.d0/P_G*sum(lam_G**2/u0_G*(eZ_nega_G+rm_G*eZ_posi_G)*htj1)

                Hy_h1 = 1.d0/P_G*sum(lam_G*(-eZ_nega_G+re_G*eZ_posi_G)*htj0)
                Hy_h2 = Hx_h1
                Hy_h3 = Hx_h2
                Hy_h4 = Hx_h3

                Hz_h1 = 1.d0/P_G*sum(lam_G**2/u0_G*(eZ_nega_G+re_G*eZ_posi_G)*htj1)

                if (zo_tr(i)-gz(j)<0.d0) then

                    Ex_h4 = 1.d0/P_G*sum(lam_G**2*(eZ_nega_G+rm_G*eZ_posi_G)*htj1)
                    Ey_h3 = Ex_h4
                    Ez_h1 = 1.d0/P_G*sum(lam_G**2*(eZ_nega_G-rm_G*eZ_posi_G)*htj1)
                    Hy_h1 = 1.d0/P_G*sum(lam_G*(eZ_nega_G+re_G*eZ_posi_G)*htj0)

                end if

                Ex_G = Ex_G+I_current/4.d0/pi*gw(j)*&
                            ((Y**2-X_G**2)/P_G**3*Ex_h1+X_G**2/P_G**2*Ex_h2-&
                             impe0*Ex_h3+sinh((gx(j)-a)/q_coef)*X_G/P_G/admi0*Ex_h4)

                Ey_G = Ey_G+I_current/4.d0/pi*gw(j)*&
                            (-2.d0*Y*X_G/P_G**3*Ey_h1+X_G*Y/P_G**2*Ey_h2+&
                             sinh((gx(j)-a)/q_coef)*Y/P_G/admi0*Ey_h3)

                Ez_G = Ez_G+I_current/4.d0/pi*gw(j)*&
                            (-X_G/P_G/admi0*Ez_h1-sinh((gx(j)-a)/q_coef)/admi0*Ez_h2)

                Hx_G = Hx_G+I_current/4.d0/pi*gw(j)*&
                            (-2.d0*Y*X_G/P_G**3*Hx_h1+X_G*Y/P_G**2*Hx_h2+&
                             sinh((gx(j)-a)/q_coef)*Y/P_G*Hx_h3)

                Hy_G = Hy_G+I_current/4.d0/pi*gw(j)*&
                            (-(Y**2-X_G**2)/P_G**3*Hy_h2-X_G**2/P_G**2*Hy_h3+Hy_h1-&
                             sinh((gx(j)-a)/q_coef)*X_G/P_G*Hy_h4)

                Hz_G = Hz_G+I_current/4.d0/pi*gw(j)*(Y/P_G*Hz_h1)

            end do
            !$OMP END PARALLEL DO
            

            Ex_t(i) = Ex_G
            Ey_t(i) = Ey_G
            Ez_t(i) = Ez_G
            Hx_t(i) = Hx_G
            Hy_t(i) = Hy_G
            Hz_t(i) = Hz_G

        end do
        !$OMP END PARALLEL DO
    end subroutine freq_sounding

!==============================================================================!
!============================================= Computes reflection coefficient !
!==============================================================================!     
    ! For TM mode
    function RTM(layernum,sig,mu,eps,h,upper_u,upper_admi,lamda,w)
        
        implicit none

        integer, intent(in) :: layernum
        real(8), intent(in) :: sig(:), mu(:), eps(:), h(:)
        complex(8), intent(in) :: upper_u(:), upper_admi
        real(8), intent(in) :: lamda(:), w
        complex(8) :: RTM(size(lamda))

        integer :: j
        complex(8) :: Z_n(size(lamda)), Z_n_in(size(lamda))
        complex(8) :: Z0(size(lamda))
        complex(8) :: impe_n, admi_n, k_n, u_n(size(lamda))

        Z0 = upper_u/upper_admi ! Intrinsic impedance of the air layer

        impe_n = dcmplx(0.d0,mu(layernum)*w)
        admi_n = dcmplx(sig(layernum),eps(layernum)*w)
        k_n = sqrt(-impe_n*admi_n)
        u_n = sqrt(lamda**2-k_n**2)
        Z_n = u_n/admi_n ! Intrinsic impedance of the last layer
        Z_n_in = Z_n     ! Interface impedance of the last layer

        do j = layernum-1, 1, -1
            impe_n = dcmplx(0.d0,mu(j)*w)
            admi_n = dcmplx(sig(j),eps(j)*w)
            k_n = sqrt(-impe_n*admi_n)
            u_n = sqrt(lamda**2-k_n**2)
            ! Intrinsic impedance of the jth layer
            Z_n = u_n/admi_n 
            ! Interface impedance of the jth layer's upper surface
            Z_n_in = Z_n*(Z_n_in+Z_n*tanh(u_n*h(j)))/&
                         (Z_n+Z_n_in*tanh(u_n*h(j)))
        end do

        RTM = (Z0-Z_n_in)/(Z0+Z_n_in) ! Reflection coefficient of the air-ground interface

        return

    end function RTM

    ! For TE mode
    function RTE(layernum,sig,mu,eps,h,upper_u,upper_impe,lamda,w)
        
        implicit none

        integer, intent(in) :: layernum
        real(8), intent(in) :: sig(:), mu(:), eps(:), h(:)
        complex(8), intent(in) :: upper_u(:), upper_impe
        real(8), intent(in) :: lamda(:), w
        complex(8) :: RTE(size(lamda))

        integer :: j
        complex(8) :: Y_n(size(lamda)), Y_n_in(size(lamda))
        complex(8) :: Y0(size(lamda))
        complex(8) :: impe_n, admi_n, k_n, u_n(size(lamda))

        Y0 = upper_u/upper_impe ! Intrinsic admittance of the air layer

        impe_n = dcmplx(0.d0,mu(layernum)*w)
        admi_n = dcmplx(sig(layernum),eps(layernum)*w)
        k_n = sqrt(-impe_n*admi_n)
        u_n = sqrt(lamda**2-k_n**2)
        Y_n = u_n/impe_n ! Intrinsic admittance of the last layer
        Y_n_in = Y_n     ! Interface admittance of the last layer

        do j = layernum-1, 1, -1
            impe_n = dcmplx(0.d0,mu(j)*w)
            admi_n = dcmplx(sig(j),eps(j)*w)
            k_n = sqrt(-impe_n*admi_n)
            u_n = sqrt(lamda**2-k_n**2)
            ! Intrinsic admittance of the jth layer
            Y_n = u_n/impe_n 
            ! Interface admittance of the jth layer's upper surface
            Y_n_in = Y_n*(Y_n_in+Y_n*tanh(u_n*h(j)))/&
                         (Y_n+Y_n_in*tanh(u_n*h(j)))
        end do

        RTE = (Y0-Y_n_in)/(Y0+Y_n_in) ! Reflection coefficient of the air-ground interface

        return

    end function RTE

end module freq_forward
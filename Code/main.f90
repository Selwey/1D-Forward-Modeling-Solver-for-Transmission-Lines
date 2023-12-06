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

program main

!============================================================================================
!The main program, which is used for data preparation, module invocation, and result output.
!============================================================================================

    use parameter_set
    use GS_Hk
    use coor_trans_intrans
    use freq_forward
    !$ use omp_lib
    
    implicit none

    real(8) :: w
    real(8), allocatable :: xo(:), yo(:), zo(:), xo_tr(:), yo_tr(:), zo_tr(:)
    real(8) :: P_start(3), P_end(3)
    real(8) :: q_coef, alpha, h1, h2, l, a
    complex(8) :: I_current
    complex(8), allocatable :: Ex_t(:), Ey_t(:), Ez_t(:), Hx_t(:), Hy_t(:), Hz_t(:)
    complex(8), allocatable :: Ex_intr(:), Ey_intr(:), Ez_intr(:),&
                               Hx_intr(:), Hy_intr(:), Hz_intr(:)
    complex(8), allocatable :: Ex(:), Ey(:), Ez(:), Hx(:), Hy(:), Hz(:)
    
    integer :: i, j
    integer :: num_catenary, mth_line, nth_point, oth_circuit
    integer :: pth_line_each_circuit, qth_ph_each_circuit
    integer, allocatable :: num_line_circuit(:)

	character(len=300) :: currentpath
    character(len=10) :: ch

    real :: t1, t2
    
    !$ real :: t3, t4

    !$ integer :: num_cup
    
    !$ num_cup = omp_get_num_procs()
    !$ call omp_set_num_threads(num_cup)

    call cpu_time(t1)
    !$ t3 = omp_get_wtime()

    call getcwd(currentpath)  ! The working path
    call power_line_parameter ! Read power-line parameters
    call receiver_parameter   ! Read observation parameters
    call model_parameter      ! Read model parameters
    call Hankel_parameter     ! Read parameters used for Hankel transformaion
    call Gauss_parameter      ! Read parameters used for Gauss integral

    ! Angular frequency
    w = 2*pi*freq_I
    ! Coordinates of observation points
    allocate(xo(num_o),yo(num_o),zo(num_o))
    allocate(xo_tr(num_o),yo_tr(num_o),zo_tr(num_o))
    xo(1) = xos
    yo(1) = yos
    zo(1) = zos
    if (num_o>1) then
        do i = 2, num_o
            xo(i) = xos+(i-1)*(xoe-xos)/(num_o-1)
            yo(i) = yos+(i-1)*(yoe-yos)/(num_o-1)
            zo(i) = zos+(i-1)*(zoe-zos)/(num_o-1)
        end do
    end if

    ! Total number of catenaries
    num_catenary = 3*sum(bundle_num)*num_span

    ! num_line_circuit(i): Total number of bundled conductors in the first i circuits
    allocate(num_line_circuit(num_circuit))
    do i = 1, num_circuit
        num_line_circuit(i) = 3*sum(bundle_num(1:i))
    end do

    ! Responses of each catenary at all measuring points (after coordinate transformation)
    allocate(Ex_t(num_o),Ey_t(num_o),Ez_t(num_o),Hx_t(num_o),Hy_t(num_o),Hz_t(num_o))
    ! Responses of each catenary at all measuring points
    allocate(Ex_intr(num_o),Ey_intr(num_o),Ez_intr(num_o),&
             Hx_intr(num_o),Hy_intr(num_o),Hz_intr(num_o))
    ! Responses at all measuring points
    allocate(Ex(num_o),Ey(num_o),Ez(num_o),Hx(num_o),Hy(num_o),Hz(num_o))
    Ex = (0.d0,0.d0)
    Ey = (0.d0,0.d0)
    Ez = (0.d0,0.d0)
    Hx = (0.d0,0.d0)
    Hy = (0.d0,0.d0)
    Hz = (0.d0,0.d0)

    oth_circuit = 1
    do i = 1, num_catenary

        if (mod(i,num_span) > 0) then
            ! The mth bundled conductor for all circuits to which the current catenary belongs
            mth_line = int(i/num_span)+1 
            ! The nth catenary for the mth bundled conductor which the current catenary is
            nth_point = mod(i,num_span)  
        else
            mth_line = int(i/num_span)
            nth_point = mod(i,num_span)+num_span
        end if

        ! The oth circuit to which the current catenary belongs
        if (mth_line > num_line_circuit(oth_circuit)) then
            oth_circuit = oth_circuit+1
        end if

        ! The pth bundled conductor for the oth circuit to which the current catenary belongs 
        pth_line_each_circuit=3*bundle_num(oth_circuit)+&
                              mth_line-num_line_circuit(oth_circuit)
        
        ! The qth phase for the oth circuit to which the current catenary belongs 
        qth_ph_each_circuit=int(pth_line_each_circuit/bundle_num(oth_circuit))+1
        if (mod(pth_line_each_circuit,bundle_num(oth_circuit))==0) then
            qth_ph_each_circuit=int(pth_line_each_circuit/bundle_num(oth_circuit))
        end if

        ! Coordinate of the starting point of the ith catenary
        P_start = [x_coor_suspension(nth_point,mth_line),&
                   y_coor_suspension(nth_point,mth_line),&
                   z_coor_suspension(nth_point,mth_line)]
        ! Coordinate of the ending point of the ith catenary
        P_end = [x_coor_suspension(nth_point+1,mth_line),&
                 y_coor_suspension(nth_point+1,mth_line),&
                 z_coor_suspension(nth_point+1,mth_line)]
        ! Catenary coefficient of the ith catenary
        q_coef =  q(nth_point,3*(oth_circuit-1)+qth_ph_each_circuit)

        ! Current loading on the ith catenary
        I_current = amp_I(qth_ph_each_circuit,oth_circuit)*&
                    dcmplx(cos(pi*ph_I(qth_ph_each_circuit,oth_circuit)/180.0),&
                           sin(pi*ph_I(qth_ph_each_circuit,oth_circuit)/180.0))
        ! write(*,*) I_current
        ! Coordinate transformation
        call coor_trans(P_start,P_end,q_coef,xo,yo,zo,alpha,h1,h2,l,a,xo_tr,yo_tr,zo_tr)
        if (q_coef*(cosh(a/q_coef)-1)+h1>0.0) then
            write(*,*) 'Warning! Please ensure the lowest point of the',i,'th catenary is above the ground!'
            stop
        end if
        
        ! Forward modeling in frequency domain
        call freq_sounding(h1, h2, l, q_coef, a, I_current, w, xo_tr, yo_tr, zo_tr,&
                           Ex_t, Ey_t, Ez_t, Hx_t, Hy_t, Hz_t)

        ! Inverse coordinate transformation of responses
        call in_trans(Ex_t, Ey_t, Ez_t, Hx_t, Hy_t, Hz_t, alpha,&
                      Ex_intr, Ey_intr, Ez_intr, Hx_intr, Hy_intr, Hz_intr)

        ! Final responses at all measuring points
        Ex = Ex + Ex_intr
        Ey = Ey + Ey_intr
        Ez = Ez + Ez_intr
        Hx = Hx + Hx_intr
        Hy = Hy + Hy_intr
        Hz = Hz + Hz_intr

        ! Output temporary files
        ! write(ch,'(I10)') i
        ! open(unit=16, file=trim(currentpath)//'/single_catenary_'//trim(adjustl(ch))//'.txt')
        ! do j = 1,num_o
        !     write(16,'(3F12.3,12E16.7)') xo(j),yo(j),zo(j),&
        !                                  real(Ex_intr(j)),aimag(Ex_intr(j)),&
        !                                  real(Ey_intr(j)),aimag(Ey_intr(j)),&
        !                                  real(Ez_intr(j)),aimag(Ez_intr(j)),&
        !                                  real(Hx_intr(j)),aimag(Hx_intr(j)),&
        !                                  real(Hy_intr(j)),aimag(Hy_intr(j)),&
        !                                  real(Hz_intr(j)),aimag(Hz_intr(j))
        ! enddo
        ! close(16)
        
    end do

    open(unit=16, file=trim(currentpath)//'/response.txt')
        do j = 1,num_o
            write(16,'(3F12.3,12E16.7)') xo(j),yo(j),zo(j),&
                                         real(Ex(j)),aimag(Ex(j)),&
                                         real(Ey(j)),aimag(Ey(j)),&
                                         real(Ez(j)),aimag(Ez(j)),&
                                         real(Hx(j)),aimag(Hx(j)),&
                                         real(Hy(j)),aimag(Hy(j)),&
                                         real(Hz(j)),aimag(Hz(j))
        enddo
    close(16)
    
    call cpu_time(t2)
    !$ t4 = omp_get_wtime()
    write(*,*) 'The cpu time is:', t2-t1, 's'
    !$ write(*,*) 'The run time is:', t4-t3, 's'

end program main
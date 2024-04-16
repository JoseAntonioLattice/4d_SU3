module starts

  use iso_fortran_env, only : dp => real64, i4 => int32
  implicit none

  private
  public :: cold_start, hot_start, cross, det

contains


  subroutine cold_start(U)
    use data_types_observables, only : link_variable
    type(link_variable), intent(out), dimension(:,:,:,:) :: U
    integer(i4) :: Lx,Lt,x,y,z,w,mu

    Lx = size(U(:,1,1,1))
    Lt = size(U(1,1,1,:)) 
    do x = 1, Lx
       do y = 1, Lx
          do z = 1, Lx
             do w = 1, Lt
                do mu = 1, 4
                   U(x,y,z,w)%link(mu)%matrix = reshape([1.0_dp,0.0_dp,0.0_dp, 0.0_dp,1.0_dp,0.0_dp, 0.0_dp,0.0_dp,1.0_dp], [3,3])
                end do
             end do
          end do
       end do
    end do

  end subroutine cold_start

  subroutine hot_start(U)
    use data_types_observables, only : link_variable
    type(link_variable), intent(out), dimension(:,:,:,:) :: U
    integer(i4) :: x,y,z,w,mu, Lx,Lt
    real(dp), dimension(6) :: r1, r2
    complex(dp), dimension(3) :: uu, vv

    Lx = size(U(:,1,1,1))
    Lt = size(U(1,1,1,:))
    do x = 1, Lx
       do y = 1, Lx
          do z = 1, Lx
             do w = 1, Lt
                do mu = 1, 4
                   call random_number(r1); call random_number(r2)
                   r1 = r1 - 0.5_dp ; r2 = r2 - 0.5_dp;
                   r1 = r1/norm2(r1); r2 = r2/norm2(r2)
                   uu = [cmplx(r1(1),r1(2),dp), cmplx(r1(3),r1(4),dp), cmplx(r1(5), r1(6), dp)]
                   vv = [cmplx(r2(1),r2(2),dp), cmplx(r2(3),r2(4),dp), cmplx(r2(5), r2(6), dp)]
                   vv = vv - uu * dot_product(uu,vv)
                   vv = vv / sqrt((vv(1)%re)**2 + (vv(1)%im)**2 + (vv(2)%re)**2 + (vv(2)%im)**2 +  (vv(3)%re)**2 + (vv(3)%im)**2)
                   U(x,y,z,w) % link(mu) % matrix = transpose( reshape([uu,vv,cross(conjg(uu), conjg(vv))], [3,3]) )
                end do
             end do
          end do
       end do
    end do

  end subroutine hot_start

  function cross(u,v)
    complex(dp), dimension(3), intent(in) :: u, v
    complex(dp), dimension(3) :: cross

    cross(1) = u(2) * v(3) - u(3) * v(2) 
    cross(2) =-u(1) * v(3) + u(3) * v(1)
    cross(3) = u(1) * v(2) - u(2) * v(1)
  end function cross


  function det(A)

    complex(dp), dimension(3,3), intent(in) :: A
    integer(i4) :: i, info
    integer :: ipiv(3)
    real(dp) :: sgn
    complex(dp) :: det

    ipiv = 0

    call zgetrf(3,3,A,3,ipiv,info)

    det = 1.0_dp

    do i = 1, 3
       det = det*A(i,i)
    end do

    sgn = 1.0_dp
    do i = 1, 3
       if(ipiv(i) /= i) sgn = -sgn
    end do

    det = sgn*det
    
  end function det
  
end module starts

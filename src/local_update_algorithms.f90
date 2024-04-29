module local_update_algorithms

  use data_types_observables
  use matrix_operations
  !use dynamics, only : staples
  use iso_fortran_env, only : dp => real64, i4 => int32
  use periodic_boundary_conditions_mod, only : ip_func, im_func
  use get_index_mod

  implicit none

contains

  subroutine metropolis(Delta_S,U,Up)
    real(dp), intent(in) :: Delta_S
    complex(dp), dimension(:,:), intent(inout) :: U
    complex(dp), dimension(:,:), intent(in) :: Up
    real(dp) :: r, prob

    prob = min(1.0_dp,exp(-Delta_S))

    call random_number(r)
    if( prob >= r )then
       U = Up
    end if

  end subroutine metropolis


 subroutine glauber(Delta_S,U,Up)
    real(dp), intent(in) :: Delta_S
    complex(dp), dimension(:,:), intent(inout) :: U
    complex(dp), dimension(:,:), intent(in) :: Up
    real(dp) :: r

    call random_number(r)
    if ( 1/(exp(Delta_S) + 1.0_dp) > r)then
       U = Up
    end if

  end subroutine glauber

  subroutine heatbath(U,x,mu,beta)
    type(link_variable), intent(inout) :: U(:,:,:,:)
    integer(i4), intent(in) :: x(4)
    integer(i4), intent(in) :: mu
    real(dp), intent(in) :: beta
    type(complex_3x3_matrix) :: R, S, T, W, A, Up
    complex(dp), dimension(2,2) :: R2, S2, T2, W2

    
    A = dagger(staples(U,x,mu))
    W = U(x(1),x(2),x(3),x(4))%link(mu) * A

    W2 = turn_2x2(W,1)
    call heatbath_SU2(W2,R2,2*beta/3)
    R = turn_to_SU3(R2,1)
    call heatbath_SU2(turn_2x2(R*W,2),S2,2*beta/3)
    S = turn_to_SU3(S2,2)
    call heatbath_SU2(turn_2x2(S*R*W,3),T2,2*beta/3)
    T = turn_to_SU3(T2,3)
    
    U(x(1),x(2),x(3),x(4))%link(mu) = T * S * R * U(x(1),x(2),x(3),x(4))%link(mu)
    
  end subroutine heatbath


  function turn_2x2(A,n) result(X)
    type(complex_3x3_matrix) :: A
    integer(i4), intent(in) :: n
    complex(dp), dimension(2,2) :: X

    select case(n)
    case(1)
       X(1,1) = A%matrix(1,1)
       X(1,2) = A%matrix(1,2)
       X(2,1) = A%matrix(2,1)
       X(2,2) = A%matrix(2,2)
    case(2)
       X(1,1) = A%matrix(1,1)
       X(1,2) = A%matrix(1,3)
       X(2,1) = A%matrix(3,1)
       X(2,2) = A%matrix(3,3)
    case(3)
       X(1,1) = A%matrix(2,2)
       X(1,2) = A%matrix(2,3)
       X(2,1) = A%matrix(3,2)
       X(2,2) = A%matrix(3,3)
    end select

  
  end function turn_2x2
  
  function turn_to_SU3(A,n) result(X)
    complex(dp), dimension(2,2) :: A
    integer(i4) :: n
    type(complex_3x3_matrix) :: X

    X%matrix = 0.0_dp

    select case(n)
    case(1)
       X%matrix(1,1) = A(1,1)
       X%matrix(1,2) = A(1,2)
       X%matrix(2,1) = A(2,1)
       X%matrix(2,2) = A(2,2)
       X%matrix(3,3) = 1.0_dp
    case(2)
       X%matrix(1,1) = A(1,1)
       X%matrix(1,3) = A(1,2)
       X%matrix(3,1) = A(2,1)
       X%matrix(3,3) = A(2,2)
       X%matrix(2,2) = 1.0_dp
    case(3)
       X%matrix(2,2) = A(1,1)
       X%matrix(2,3) = A(1,2)
       X%matrix(3,2) = A(2,1)
       X%matrix(3,3) = A(2,2)
       X%matrix(1,1) = 1.0_dp
    end select
    
  end function turn_to_SU3
  
  subroutine heatbath_SU2(W,R,beta)
    complex(dp), dimension(2,2), intent(in)  :: W
    complex(dp), dimension(2,2), intent(out) :: R
    real(dp), intent(in) :: beta

    complex(dp), dimension(2,2) :: V

    real(dp) :: k, a, b, s, u, x0,x(3)

    logical :: accept
    
    R = turn_to_SU2_matrix(W)
    k = sqrt(det2(R))

    V = R/k

    a = exp(-2*k*beta)
    b = 1.0_dp

    accept = .false.
    do while( accept .eqv. .false. )
       u = random_uniform(a,b)
       x0 = 1.0_dp + log(u)/(beta*k)
       call random_number(s)
       if( s > 1.0_dp - sqrt(1.0_dp - x0**2)) accept = .true.
    end do
    
    x = sqrt(1.0_dp - x0**2) * random_vector()

    R = SU2_matrix(cmplx(x0,x(1),dp),cmplx(x(2),x(3),dp))

    R = matmul(R,conjg(transpose(V)))
    
  end subroutine heatbath_SU2

  function turn_to_SU2_matrix(M) result(X)
    complex(dp), dimension(2,2), intent(in) :: M
    complex(dp), dimension(2,2) :: X

    complex(dp) :: a, b

    a = 0.5 * ( M(1,1) + conjg(M(2,2)) )
    b = 0.5 * ( M(1,2) - conjg(M(2,1)) )

    X = SU2_matrix(a,b)
    
  end function turn_to_SU2_matrix

  
  function det2(W)
    complex(dp), dimension(2,2) :: W
    real(dp) :: det2
    
    det2 = (W(1,1)%re)**2 + (W(1,1)%im)**2 + (W(1,2)%re)**2 + (W(1,2)%im)**2
        
  end function det2
  
  
  function random_vector() result(y)
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), dimension(3) :: y
    real(dp) :: r, n

    call random_number(r)

    
    y(1) = random_uniform(-1.0_dp, 1.0_dp)
    n = sqrt(1.0_dp - (y(1))**2)
    y(2) = n * cos(2*pi*r)
    y(3) = n * sin(2*pi*r)

  end function random_vector

  function SU2_matrix(a,b) result(matrix)
    complex(dp), intent(in)  :: a, b
    complex(dp), dimension(2,2) :: matrix
      matrix(1,1) = a
      matrix(1,2) = b
      matrix(2,1) = -conjg(b)
      matrix(2,2) =  conjg(a)
  end function SU2_matrix

  subroutine generate_lambdasq(det_A,beta,lambdasq,s)
    real(dp), intent(in) :: det_A,beta
    real(dp), intent(out):: lambdasq, s
    real(dp) :: r(3)
    real(dp), parameter :: pi = acos(-1.0_dp)

    call random_number(r)
    r = 1.0_dp - r
    lambdasq = -(1/(2*det_A*beta)) * ( log(r(1)) + (cos(2*pi*r(2)))**2 * log(r(3)) )
    call random_number(s)
  end subroutine generate_lambdasq


  subroutine create_complex_numbers(a,b)
    complex(dp), intent(out) :: a, b
    real(dp), dimension(0:3) :: r, x
    real(dp), parameter :: eps = 0.1_dp
    real(dp) :: norm_r

    call random_number(r)
    r = r - 0.5_dp
    norm_r = sqrt(r(1)**2 + r(2)**2 + r(3)**2)

    x(1:3) = eps*r(1:3)/norm_r
    x(0) = sgn(r(0)) * sqrt(1.0_dp - eps**2)

    a = cmplx(x(0),x(1),dp)
    b = cmplx(x(2),x(3),dp)

  end subroutine create_complex_numbers

  subroutine create_biased_update(X)
    type(complex_3x3_matrix), intent(out) :: X
    complex(dp), dimension(3) :: uu, vv
    real(dp), dimension(6) :: r1,r2
    type(complex_3x3_matrix) :: R, S, T
    complex(dp) :: a, b

    R%matrix = reshape([1.0_dp,0.0_dp,0.0_dp, 0.0_dp,1.0_dp,0.0_dp, 0.0_dp,0.0_dp,1.0_dp], [3,3])
    S = R
    T = R

    call create_complex_numbers(a,b)
    R%matrix(1,1) = a; R%matrix(1,2) = b; R%matrix(2,1)= -conjg(b); R%matrix(2,2) = conjg(a)

    
    call create_complex_numbers(a,b)
    S%matrix(1,1) = a; S%matrix(1,3) = b; S%matrix(3,1)= -conjg(b); S%matrix(3,3) = conjg(a)

    
    call create_complex_numbers(a,b)
    T%matrix(2,2) = a; T%matrix(2,3) = b; T%matrix(3,2)= -conjg(b); T%matrix(3,3) = conjg(a)

    X = R * S * T
    
  end subroutine create_biased_update

  pure function sgn(x)
    real(dp), intent(in) :: x
    integer(i4) :: sgn

    if( x > 0.0_dp )then
      sgn = 1
    elseif( x < 0.0_dp)then
      sgn = -1
    else
      sgn = 0
    end if

  end function sgn


  function cross(u,v)
    complex(dp), dimension(3), intent(in) :: u, v
    complex(dp), dimension(3) :: cross
    
    cross(1) = u(2) * v(3) - u(3) * v(2) 
    cross(2) =-u(1) * v(3) + u(3) * v(1)
    cross(3) = u(1) * v(2) - u(2) * v(1)
  end function cross
    
  function DS(U,mu,Up,x,beta,N)
    type(link_variable), dimension(:,:,:,:), intent(in) :: U
    integer(i4), intent(in) :: x(4), mu
    real(dp), intent(in) :: beta
    integer(i4), intent(in) :: N
    type(complex_3x3_matrix), intent(out) :: Up
    real(dp) :: DS

    call create_biased_update(Up)
    Up = Up * U(x(1),x(2),x(3),x(4))%link(mu)
    !call create_update(Up); Up = Up * U(x,y)%link(mu)

    DS = - (beta/N) * real( tr( (Up - U(x(1),x(2),x(3),x(4))%link(mu)) * dagger(staples(U,x,mu)) ) )

  end function DS

      function staples(U,x,mu) result(A)

    type(link_variable), dimension(:,:,:,:), intent(in) :: U
    integer(i4), intent(in) :: x(4), mu
    integer(i4) :: nu
    integer(i4), parameter :: d = 4
    type(complex_3x3_matrix) :: A

    integer(i4), dimension(4) :: ipx_mu, ipx_nu, imx_nu, ipx_mu_imx_nu

    A%matrix = 0.0_dp
    ipx_mu = ip_func(x,mu)
    do nu = 1, d
       if(nu .ne. mu)then
          
          ipx_nu = ip_func(x,nu)
          imx_nu = im_func(x,nu)
          ipx_mu_imx_nu = ip_func(imx_nu,mu)

          A = A +    U(   x(1)  ,   x(2)  ,   x(3)  ,   x(4)  )%link(nu)  &
                   * U(ipx_nu(1),ipx_nu(2),ipx_nu(3),ipx_nu(4))%link(mu)  &
            * dagger(U(ipx_mu(1),ipx_mu(2),ipx_mu(3),ipx_mu(4))%link(nu)) &
            + dagger(U(imx_nu(1),imx_nu(2),imx_nu(3),imx_nu(4))%link(nu)) &
                   * U(imx_nu(1),imx_nu(2),imx_nu(3),imx_nu(4))%link(mu)  &
                   * U(ipx_mu_imx_nu(1),ipx_mu_imx_nu(2),ipx_mu_imx_nu(3),ipx_mu_imx_nu(4))%link(nu)
       end if
    end do
  end function staples


  function random_uniform(a,b) result(y)
    real(dp), intent(in) :: a, b
    real(dp) :: y, r

    call random_number(r)

    y = a + r * ( b - a )

  end function random_uniform

end module local_update_algorithms

module dynamics

  use iso_fortran_env, only: dp => real64, i4 => int32
  use data_types_observables, only : link_variable, complex_3x3_matrix
  use matrix_operations
  use local_update_algorithms
  use periodic_boundary_conditions_mod
  use get_index_mod
  use create_files

  implicit none
  integer(i4), dimension(4,4,4,4) :: levi_civita
  real(dp), parameter :: pi = acos(-1.0_dp)
  !private !:: dp, i4, link_variable
  !public :: sweeps, create_update, sgn, drand, take_measurements, dagger, tr, gauge_transformation, DS

contains

  subroutine set_levi_civita()

    levi_civita = 0
    
    levi_civita(1,2,3,4) = 1
    levi_civita(1,3,4,2) = 1
    levi_civita(1,4,2,3) = 1
    levi_civita(2,1,4,3) = 1
    levi_civita(2,3,1,4) = 1
    levi_civita(2,4,3,1) = 1
    levi_civita(3,1,2,4) = 1
    levi_civita(3,2,4,1) = 1
    levi_civita(3,4,1,2) = 1
    levi_civita(4,1,3,2) = 1
    levi_civita(4,2,1,3) = 1
    levi_civita(4,3,2,1) = 1

    levi_civita(1,2,4,3) = -1
    levi_civita(1,3,2,4) = -1
    levi_civita(1,4,3,2) = -1
    levi_civita(2,1,3,4) = -1
    levi_civita(2,3,4,1) = -1
    levi_civita(2,4,1,3) = -1
    levi_civita(3,1,4,2) = -1
    levi_civita(3,2,1,4) = -1
    levi_civita(3,4,2,1) = -1
    levi_civita(4,1,2,3) = -1
    levi_civita(4,2,3,1) = -1
    levi_civita(4,3,1,2) = -1

    !print*, levi_civita
  end subroutine set_levi_civita
  

  subroutine equilibrium_dynamics(U,Lx,Lt,beta,N,d,algorithm,N_thermalization,N_measurements,N_skip,equilibrium)
    use starts
    use statistics
    type(link_variable), intent(inout), dimension(:,:,:,:) :: U
    integer(i4), intent(in)  :: Lx,Lt, N, d
    real(dp), intent(in), dimension(:) :: beta
    character(*), intent(in) :: algorithm
    integer(i4), intent(in) :: N_thermalization, N_measurements, N_skip
    logical, intent(in) :: equilibrium
    integer(i4) :: i_beta
    
    !allocate(E_p%array(N_measurements))

    call set_levi_civita
    call cold_start(U)

    do i_beta = 1, size(beta)
       call thermalization(U,Lx,Lt,beta(i_beta),N,d,algorithm,N_thermalization)
       call create_measurements_file(Lx,Lt,beta(i_beta),algorithm,equilibrium)
       !call measurements_sweeps(U,Lx,Lt,beta(i_beta),N,d,algorithm,N_measurements,N_skip)
       call wilson_flow_euler(U)
    end do
  end subroutine equilibrium_dynamics
  
  subroutine thermalization(U,Lx,Lt,beta,N,d,algorithm,N_thermalization)
    type(link_variable), intent(inout), dimension(:,:,:,:) :: U
    integer(i4), intent(in)  :: Lx,Lt, N, d
    real(dp), intent(in) :: beta
    character(*), intent(in) :: algorithm
    integer(i4), intent(in) :: N_thermalization
    integer(i4) :: i

    do i = 1, N_thermalization
       call sweeps(U,Lx,Lt,beta,N,d,algorithm)
       !if(mod(i,10) == 0) call normalization(U,Lx,Lt)
    end do
   end subroutine thermalization


   subroutine measurements_sweeps(U,Lx,Lt,beta,N,d,algorithm,N_measurements,N_skip)
     type(link_variable), intent(inout), dimension(:,:,:,:) :: U
     integer(i4), intent(in)  :: Lx,Lt, N, d
     real(dp), intent(in) :: beta
     character(*), intent(in) :: algorithm
     integer(i4), intent(in) :: N_measurements, N_skip
     real(dp) :: E_p
     complex(dp) :: avr_polyakov_loop
     !real(dp) :: correlation_polyakov_loop(Lx/2-1)
     complex(dp) :: correlation_polyakov_loop(Lx/2-1)
     integer(i4) :: i

     do i = 1, N_measurements*N_skip
        call sweeps(U,Lx,Lt,beta,N,d,algorithm)
        if( mod(i,N_skip) == 0)then
           call take_measurements(U,Lx,Lt,E_p,avr_polyakov_loop,correlation_polyakov_loop)
           write(100,*) E_p,avr_polyakov_loop,correlation_polyakov_loop 
        end if
        if(mod(i,10) == 0) call normalization(U,Lx,Lt)
     end do 

   end subroutine measurements_sweeps

  subroutine sweeps(U,Lx,Lt,beta,N,d,algorithm)
    type(link_variable), intent(inout), dimension(:,:,:,:) :: U
    integer(i4), intent(in)  :: Lx,Lt, N, d
    real(dp), intent(in) :: beta
    character(*), intent(in) :: algorithm
    integer(i4) :: x, y,z,w, mu

    if ( N == 3 .and. algorithm == 'heatbath' )then
       do x = 1, Lx
          do y = 1, Lx
             do z = 1, Lx
                do w = 1, Lt
                   do mu = 1, d
                      call heatbath(U,[x,y,z,w],mu,beta)
                   end do
                end do
             end do
          end do
       end do
    else if( N == 3 .and. algorithm == 'overrelaxation' )then
       do x = 1, Lx
          do y = 1, Lx
             do z = 1, Lx
                do w = 1, Lt
                   do mu = 1, d
                      call overrelaxation(U,[x,y,z,w],mu)
                   end do
                end do
             end do
          end do
       end do
    end if
  end subroutine sweeps


  subroutine take_measurements(U,Lx,Lt,Ep,avr_polyakov_loop,correlation_polyakov_loop)
    type(link_variable), dimension(:,:,:,:), intent(in) :: U
    integer(i4), intent(in) ::  Lx,Lt
    real(dp), intent(out) :: Ep
    complex(dp) :: polyakov_loop_array(Lx,Lx,Lx), polyakov_loop_proj(Lx)
    complex(dp), intent(out) :: avr_polyakov_loop
    complex(dp), intent(out) :: correlation_polyakov_loop(Lx/2-1)
    integer(i4) :: x,y,z,w,t,mu,nu, xp, yp, zp, zpp
    integer(i4), parameter :: d = 4, number_of_planes = d*(d-1)/2
    complex(dp) :: avg_poly(Lx)
    

        Ep = 0.0_dp
    do x = 1, Lx
       do y = 1, Lx
          do z = 1, Lx
             polyakov_loop_array(x,y,z) = polyakov_loop(U,[x,y,z],Lt)
             do w = 1, Lt
                do mu = 1, d - 1
                   do nu = mu + 1, d
                      Ep = Ep + real(tr(plaquette(U,[x,y,z,w],mu,nu)),dp)
                   end do
                end do
             end do
          end do
       end do
    end do

    
    do x = 1, Lx
       polyakov_loop_proj(x) = sum(polyakov_loop_array(x,:,:))/Lx**2
    end do


    
    avr_polyakov_loop = sum(polyakov_loop_array)/Lx**3 

    correlation_polyakov_loop = (0.0_dp,0.0_dp)

    do t = 1, Lx/2 - 1
       correlation_polyakov_loop(t) = polyakov_loop_proj(1)*conjg(polyakov_loop_proj(t))
       !do x = 1, Lx
          !xp = mod(x+t,Lx); if(xp == 0) xp = Lx
          !do y = 1, Lx
             !yp = mod(y+t,Lx); if(yp == 0) yp = Lx
             !do z = 1, Lx
                !zp = mod(z+t,Lx); if(zp == 0) zp = Lx
                !do xp = 1, Lx
                   !do yp = 1,Lx
                      !do zp = 1,Lx
                         !zpp = mod(zp+t,Lx); if(zpp == 0) zpp = Lx
                         !correlation_polyakov_loop(t) = correlation_polyakov_loop(t) + &
                         !     polyakov_loop_array(x,y,zp) * &
                         !     conjg( polyakov_loop_array(xp,yp,zpp)) !+ &
                         !polyakov_loop_array(x,yp,z) + &
                         !polyakov_loop_array(x,y,zp) &
                         !)
                         !                correlation_polyakov_loop(t) = correlation_polyakov_loop(t) + &
                         !                     wilson_loop(U,[x,y,z],1,t,Lt,Lx) + &
                         !                     wilson_loop(U,[x,y,z],2,t,Lt,Lx) + &
                         !
                         !wilson_loop(U,[x,y,z],3,t,Lt,Lx)
                      !end do
                   !end do
             !end do
          !end do
       !end do
    end do
    
    !correlation_polyakov_loop = correlation_polyakov_loop/(Lx**3)

    !do x = 1, Lx
    !avg_poly(x) = sum(polyakov_loop_array(:,:,x))/Lx**2
    !end do

    !do t=1,Lx/2-1
    !   do x=1,Lx
    !   xp = mod(x+t,Lx); if(xp == 0) xp = Lx
    !   correlation_polyakov_loop(t) = correlation_polyakov_loop(t) + avg_poly(x)*conjg(avg_poly(xp))
    !   end do
    !end do

    !correlation_polyakov_loop = correlation_polyakov_loop/Lx
    
    Ep =  Ep/(3*number_of_planes*Lx**3*Lt)
    
  end subroutine take_measurements


  
  function polyakov_loop(U,x,L)
    type(link_variable), dimension(:,:,:,:), intent(in) :: U
    integer(i4), dimension(3), intent(in) :: x 
    integer(i4), intent(in) :: L
    type(complex_3x3_matrix) :: product
    complex(dp) :: polyakov_loop
    integer(i4) :: t

    product = U(x(1),x(2),x(3),1)%link(4)
    do t = 2, L
       product = product*U(x(1),x(2),x(3),t)%link(4)
    end do

    polyakov_loop = tr(product)
    
  end function polyakov_loop


  function wilson_loop(U,x,mu,nx,nt,Lx)
    type(link_variable), dimension(:,:,:,:), intent(in) :: U
    integer(i4), dimension(3), intent(in) :: x 
    integer(i4), intent(in) :: mu,nx,nt,Lx
    type(complex_3x3_matrix) :: product1,product2,product3,product4
    complex(dp) :: wilson_loop
    integer(i4) :: i
    integer(i4), dimension(3) :: y,xp

    y = x
    xp = x
    product1%matrix = one
    product3%matrix = one
    y(mu) = mod(y(mu)+nx,Lx); if(y(mu) == 0) y(mu) = Lx
    do i = 1,nt
       product1 = product1 * U(x(1),x(2),x(3),i)%link(4)
       product3 = product3 * U(y(1),y(2),y(3),i)%link(4)
    end do

    product2%matrix = one
    product4%matrix = one
    do i = 1,nx
       xp(mu) = mod(xp(mu)+i,Lx); if(xp(mu) == 0) xp(mu) = Lx
       product2 = product2 * U(xp(1),xp(2),xp(3),nt)%link(mu)
       product4 = product4 * U(xp(1),xp(2),xp(3),1 )%link(mu) 
    end do
    
    wilson_loop = tr(product1*product2*dagger(product3)*dagger(product4))
    
  end function wilson_loop
  
  function inv_polyakov_loop(U,x,L)
    type(link_variable), dimension(:,:,:,:), intent(in) :: U
    integer(i4), dimension(3), intent(in) :: x 
    integer(i4), intent(in) :: L
    type(complex_3x3_matrix) :: product
    complex(dp) :: inv_polyakov_loop
    integer(i4) :: t

    product = dagger(U(x(1),x(2),x(3),L)%link(4))
    do t = L-1, 1, -1
       product = product*dagger(U(x(1),x(2),x(3),t)%link(4))
    end do

    inv_polyakov_loop = tr(product)
    
  end function inv_polyakov_loop
  
  
  function action(U,Lx,Lt,beta_N)

    type(link_variable), dimension(:,:,:,:), intent(in) :: U
    integer(i4), intent(in) ::  Lx,Lt
    
    real(dp) :: action, beta_N
    integer(i4) :: x,y,z,w,mu,nu
    integer(i4), parameter :: d = 4, number_of_planes = d*(d-1)/2
    
    
    action = 0.0_dp

    do x = 1, Lx
       do y = 1, Lx
          do z = 1, Lx
             do w = 1, Lt
                do mu = 1, d - 1
                   
                   do nu = mu + 1, d
                      action = action + real(tr(plaquette(U,[x,y,z,w],mu,nu)),dp)
                   end do
                   
                end do
             end do
          end do
       end do
    end do
    action =  - beta_N * action/number_of_planes

  end function action

    function energy(U)

    type(link_variable), dimension(:,:,:,:), intent(in) :: U
        
    real(dp) :: energy
    integer(i4) :: x,y,z,w,mu,nu
    integer(i4), parameter :: d = 4, number_of_planes = d*(d-1)/2
    
    
    energy = 0.0_dp

    do x = 1, size(U(:,1,1,1))
       do y = 1, size(U(1,:,1,1))
          do z = 1, size(U(1,1,:,1))
             do w = 1, size(U(1,1,1,:))
                do mu = 1, d - 1
                   
                   do nu = mu + 1, d
                      energy = energy + real(tr(plaquette(U,[x,y,z,w],mu,nu)),dp)
                   end do
                   
                end do
             end do
          end do
       end do
    end do
    !action =  - beta_N * action/number_of_planes

  end function energy

  subroutine normalization(U,Lx,Lt)
    type(link_variable), dimension(:,:,:,:), intent(inout) :: U
    integer(i4), intent(in) :: Lx,Lt
    integer(i4) :: x,y,z,w, mu
    complex(dp), dimension(3) :: u_vec, v_vec
    real(dp) :: norm
    
    do x = 1, Lx
       do y = 1, Lx
          do z = 1, Lx
             do w = 1, Lt
                do mu = 1, 4
                   u_vec = U(x,y,z,w)%link(mu)%matrix(1,:)
                   norm = sqrt( (u_vec(1)%re)**2 + (u_vec(2)%re)**2 + (u_vec(3)%re)**2 +  &
                                (u_vec(1)%im)**2 + (u_vec(2)%im)**2 + (u_vec(3)%im)**2)
                   u_vec = u_vec/norm
                   U(x,y,z,w)%link(mu)%matrix(1,:) = u_vec
                   v_vec = U(x,y,z,w)%link(mu)%matrix(2,:)
                   norm = sqrt( (v_vec(1)%re)**2 + (v_vec(2)%re)**2 + (v_vec(3)%re)**2 +  &
                                (v_vec(1)%im)**2 + (v_vec(2)%im)**2 + (v_vec(3)%im)**2)
                   v_vec = v_vec/norm
                   U(x,y,z,w)%link(mu)%matrix(2,:) = v_vec
                   U(x,y,z,w)%link(mu)%matrix(3,:) = cross_3d(conjg(u_vec),conjg(v_vec))
                end do
             end do
          end do
       end do
    end do
    
    
  end subroutine normalization

  function plaquette(U,x,mu,nu)
    type(link_variable), dimension(:,:,:,:), intent(in) :: U
    integer(i4), intent(in) :: x(4), mu, nu
    type(complex_3x3_matrix) :: plaquette
    integer(i4), dimension(4) :: ipx_mu, ipx_nu
    
    !             x, mu           x + mu, nu                   x + nu, mu                    x, nu

    ipx_mu = ip_func(x,mu)
    ipx_nu = ip_func(x,nu)
    plaquette = U(x(1),x(2),x(3),x(4))%link(mu) * U(ipx_mu(1),ipx_mu(2),ipx_mu(3),ipx_mu(4))%link(nu) * &
         dagger(U(ipx_nu(1),ipx_nu(2),ipx_nu(3),ipx_nu(4))%link(mu)) * dagger(U(x(1),x(2),x(3),x(4))%link(nu))
  end function plaquette



  function zeta(U,x,mu)
    type(link_variable), dimension(:,:,:,:), intent(in) :: U
    integer(i4), intent(in) :: x(4), mu
    type(complex_3x3_matrix) :: zeta
    zeta = TA(U(x(1),x(2),x(3),x(4))%link(mu)*dagger(staples(U,x,mu)))
    !zeta%matrix = -zeta%matrix
  end function zeta

  function TA(W)
    type(complex_3x3_matrix), intent(in) :: W
    type(complex_3x3_matrix) :: TA
    TA = (W - dagger(W))
    Ta%matrix = TA%matrix/2
    TA%matrix = TA%matrix - one*tr(W - dagger(W))/6
    
  end function TA

  subroutine wilson_flow_euler(U)
    type(link_variable), dimension(:,:,:,:), intent(inout) :: U
    type(link_variable), dimension(size(U(:,1,1,1)),size(U(1,:,1,1)),size(U(1,1,:,1)),size(U(1,1,1,:))) :: V
    type(complex_3x3_matrix) :: B
    integer(i4) :: x(4), mu
    real(dp) :: epsilon = 0.01_dp
    integer :: i, x1, x2, x3, x4
    integer, parameter :: n = 100
    complex(dp) :: Q
    
    !print*, 0, topological_density(U,x)
    print*, "inside wilson flow"
    do i = 1, n
       Q = 0.0_dp
       do x1 = 1, size(U(:,1,1,1))
          do x2 = 1, size(U(1,:,1,1))
             do x3 = 1, size(U(1,1,:,1))
                do x4 = 1, size(U(1,1,1,:))
                   x = [x1,x2,x3,x4]
                   do mu = 1, 4
                      B = Zeta(U,x,mu)
                      B%matrix = B%matrix*epsilon
                      V(x(1),x(2),x(3),x(4))%link(mu)%matrix = matmul(my_exp(B%matrix) , U(x(1),x(2),x(3),x(4))%link(mu)%matrix)
                   end do
                   Q = Q + topological_density(U,x)
                end do
             end do
          end do
       end do
       U = V
       print*, i,energy(U),-real(Q/(32*pi**2))
    end do
  end subroutine wilson_flow_euler

  subroutine wilson_flow_rk4(U)
    type(link_variable), dimension(:,:,:,:), intent(inout) :: U
    type(link_variable), dimension(size(U(:,1,1,1)),size(U(1,:,1,1)),size(U(1,1,:,1)),size(U(1,1,1,:))) :: V
    type(complex_3x3_matrix) :: B
    integer(i4) :: x(4), mu
    real(dp) :: epsilon = 0.1_dp
    integer :: i, x1, x2, x3, x4
    integer, parameter :: n = 100
    complex(dp) :: Q
    
    !print*, 0, topological_density(U,x)
    print*, "inside wilson flow"
    do i = 1, n
       Q = 0.0_dp
       do x1 = 1, size(U(:,1,1,1))
          do x2 = 1, size(U(1,:,1,1))
             do x3 = 1, size(U(1,1,:,1))
                do x4 = 1, size(U(1,1,1,:))
                   x = [x1,x2,x3,x4]
                   do mu = 1, 4
                      B = Zeta(U,x,mu)
                      B%matrix = B%matrix*epsilon
                      V(x(1),x(2),x(3),x(4))%link(mu)%matrix = matmul(my_exp(B%matrix) , U(x(1),x(2),x(3),x(4))%link(mu)%matrix)
                   end do
                   Q = Q + topological_density(U,x)
                end do
             end do
          end do
       end do
       U = V
       !print*, i, action(U)!-Q/(32*pi**2)
    end do
  end subroutine wilson_flow_rk4

   function my_exp(X) result(expX)
    ! Computes exp(X) for X in the su(3) Lie algebra using Cayley-Hamilton theorem.
    !
    ! By Cayley-Hamilton, any X in su(3) (traceless, antihermitian) satisfies:
    !   X^3 = -t*X + d*I
    ! where:
    !   t = -1/2 * tr(X^2)
    !   d =  det(X)
    !
    ! Therefore exp(X) = q0*I + q1*X + q2*X^2, with coefficients q0,q1,q2
    ! computed via Horner's method on the Taylor series:
    !   exp(X) = sum_{n=0}^{inf} X^n / n!
    !
    ! Recurrence (running from n=K down to n=0):
    !   q2_new = q1_old
    !   q1_new = q0_old - t * q2_old
    !   q0_new = 1/(i+1)! + d * q2_old
    complex(dp), dimension(3,3), intent(in) :: X
    complex(dp), dimension(3,3) :: expX, B
    integer, parameter :: K = 20
    complex(dp), dimension(3,3) :: Id
    complex(dp) :: q0, q1, q2, q0old, q1old, q2old
    complex(dp) :: d, t, trB
    integer :: i
    
    Id = 0.0_dp
    Id(1,1) = 1.0_dp
    Id(2,2) = 1.0_dp
    Id(3,3) = 1.0_dp
 
    ! X^2 and its trace (needed for t)
    B  = matmul(X, X)
    trB = B(1,1) + B(2,2) + B(3,3)
 
    ! Cayley-Hamilton coefficients
    ! For traceless X: characteristic poly is lambda^3 + t*lambda - d = 0
    t = -0.5_dp * trB          ! coefficient of lambda (no ii factor needed)
    d = determinant(3, X)      ! constant term: det(X), no extra ii factor
 
    ! Horner recurrence initialised at step K
    q0old = 1.0_dp / gamma(1.0_dp*(K+1))
    q1old = (0.0_dp, 0.0_dp)
    q2old = (0.0_dp, 0.0_dp)
 
    do i = K-1, 0, -1
       q0 = 1.0_dp/gamma(1.0_dp*(i+1)) + d * q2old   ! X^3 = -t*X + d*I  => +d (not -ii*d)
       q1 = q0old - t * q2old
       q2 = q1old
       q0old = q0
       q1old = q1
       q2old = q2
    end do
 
    expX = q0*Id + q1*X + q2*B
 
  end function my_exp
 


  pure recursive function determinant(n, a) result(det)
    implicit none
    integer, intent(in) :: n
    complex(dp), dimension(3,3), intent(in) :: a
    complex(dp) :: det
    integer :: i, sgn
    complex(dp), dimension(n-1, n-1) :: b
    
    if (n == 1) then
       det = a(1,1)
    else
       det = 0.0
       sgn = 1
       do i = 1, n
          ! Extract submatrix
          b(:, :(i-1)) = a(2:, :i-1)
          b(:, i:) = a(2:, i+1:)
          det = det + sgn * a(1, i) * determinant(n-1, b)
          sgn = -sgn
       end do
    end if
  end function determinant
  
  ! -- Returns the inverse of a general squared matrix A
  function inv(A) result(Ainv)
    implicit none
    type(complex_3x3_matrix)::  A
    type(complex_3x3_matrix) :: Ainv
    complex(dp)            :: work(3)            ! work array for LAPACK
    integer         :: n,info,ipiv(3)     ! pivot indices
    
    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv%matrix = A%matrix
    n = 3
    ! SGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call zGETRF(n,n,Ainv%matrix,n,ipiv,info)
    if (info.ne.0) stop 'Matrix is numerically singular!'
    ! SGETRI computes the inverse of a matrix using the LU factorization
    ! computed by SGETRF.
    call zGETRI(n,Ainv%matrix,n,ipiv,work,n,info)
    if (info.ne.0) stop 'Matrix inversion failed!'
  end function inv
  
  function F(U,x,mu,nu)
    type(link_variable), dimension(:,:,:,:), intent(in) :: U
    integer(i4), intent(in) :: x(4), mu, nu
    type(complex_3x3_matrix) :: F,Q
    integer(i4), dimension(4) :: ipx_mu, ipx_nu,imx_mu, imx_nu, &
         x_im_mu_ip_nu, x_im_mu_im_nu,  x_ip_mu_im_nu


    ipx_mu = ip_func(x,mu)
    ipx_nu = ip_func(x,nu)
    imx_mu = im_func(x,mu)
    imx_nu = im_func(x,nu)

    x_im_mu_ip_nu = im_func(ipx_nu,mu)
    x_im_mu_im_nu = im_func(imx_nu,mu)
    x_ip_mu_im_nu = ip_func(imx_nu,mu)

    
    Q = U(x(1),x(2),x(3),x(4))%link(mu) * &
         U(ipx_mu(1),ipx_mu(2),ipx_mu(3),ipx_mu(4))%link(nu) * &
         dagger(U(ipx_nu(1),ipx_nu(2),ipx_nu(3),ipx_nu(4))%link(mu)) *&
         dagger(U(x(1),x(2),x(3),x(4))%link(nu)) + &
         U(x(1),x(2),x(3),x(4))%link(nu) * &
         dagger(U(x_im_mu_ip_nu(1),x_im_mu_ip_nu(2),x_im_mu_ip_nu(3),x_im_mu_ip_nu(4))%link(mu)) * &
         dagger(U(imx_mu(1),imx_mu(2),imx_mu(3),imx_mu(4))%link(nu)) &
         * U(imx_mu(1),imx_mu(2),imx_mu(3),imx_mu(4))%link(mu) + &
         dagger(U(imx_mu(1),imx_mu(2),imx_mu(3),imx_mu(4))%link(mu)) &
         * dagger(U(x_im_mu_im_nu(1),x_im_mu_im_nu(2),x_im_mu_im_nu(3),x_im_mu_im_nu(4))%link(nu)) * &
         U(x_im_mu_im_nu(1),x_im_mu_im_nu(2),x_im_mu_im_nu(3),x_im_mu_im_nu(4))%link(mu) &
         * U(imx_nu(1),imx_nu(2),imx_nu(3),imx_nu(4))%link(nu) + &
         dagger(U(imx_nu(1),imx_nu(2),imx_nu(3),imx_nu(4))%link(nu))&
         * U(imx_nu(1),imx_nu(2),imx_nu(3),imx_nu(4))%link(mu) * &
         U(x_ip_mu_im_nu(1),x_ip_mu_im_nu(2),x_ip_mu_im_nu(3),x_ip_mu_im_nu(4))%link(nu) &
         * dagger(U(x(1),x(2),x(3),x(4))%link(mu))


    F = (Q - dagger(Q))/8.0_dp
    
  end function F

  function topological_density(U,x)
    type(link_variable), dimension(:,:,:,:), intent(in) :: U
    integer(i4), intent(in) :: x(4)
    integer(i4) :: mu, nu, rho, sigma
    complex(dp) :: topological_density

    topological_density = 0.0_dp
    do mu = 1, 4
       do nu = 1, 4
          do rho = 1, 4
             do sigma = 1, 4
                topological_density = topological_density - levi_civita(mu,nu,rho,sigma)*tr(F(U,x,mu,nu)*F(U,x,rho,sigma))
             end do
          end do
       end do
    end do
    !topological_density = -topological_density/(32*pi**2)
  end function topological_density

  
end module dynamics

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
       call measurements_sweeps(U,Lx,Lt,beta(i_beta),N,d,algorithm,N_measurements,N_skip)
       
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
       if(mod(i,10) == 0) call normalization(U,Lx,Lt)
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
    complex(dp) :: polyakov_loop_array(Lx,Lx,Lx)
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

    avr_polyakov_loop = sum(polyakov_loop_array)/Lx**3 

    correlation_polyakov_loop = (0.0_dp,0.0_dp)

    do t = 1, Lx/2 - 1
       do x = 1, Lx
          !xp = mod(x+t,Lx); if(xp == 0) xp = Lx
          do y = 1, Lx
             !yp = mod(y+t,Lx); if(yp == 0) yp = Lx
             !do z = 1, Lx
                !zp = mod(z+t,Lx); if(zp == 0) zp = Lx
                do xp = 1, Lx
                   do yp = 1,Lx
                      do zp = 1,Lx
                         zpp = mod(zp+t,Lx); if(zpp == 0) zpp = Lx
                         correlation_polyakov_loop(t) = correlation_polyakov_loop(t) + &
                              polyakov_loop_array(x,y,zp) * &
                              conjg( polyakov_loop_array(xp,yp,zpp)) !+ &
                         !polyakov_loop_array(x,yp,z) + &
                         !polyakov_loop_array(x,y,zp) &
                         !)
                         !                correlation_polyakov_loop(t) = correlation_polyakov_loop(t) + &
                         !                     wilson_loop(U,[x,y,z],1,t,Lt,Lx) + &
                         !                     wilson_loop(U,[x,y,z],2,t,Lt,Lx) + &
                         !
                         !wilson_loop(U,[x,y,z],3,t,Lt,Lx)
                      end do
                   end do
             end do
          end do
       end do
    end do
    
    correlation_polyakov_loop = correlation_polyakov_loop/(Lx**3)

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
    zeta%matrix = -zeta%matrix
  end function zeta

  function TA(W)
    type(complex_3x3_matrix), intent(in) :: W
    type(complex_3x3_matrix) :: TA
    TA = (W - dagger(W))
    Ta%matrix = TA%matrix/2
    TA%matrix = TA%matrix - one*tr(W - dagger(W))/6
    
  end function TA

  subroutine wilson_flow_euler(U,V,x,mu)
    type(link_variable), dimension(:,:,:,:), intent(inout) :: U
    type(complex_3x3_matrix), intent(out) :: V
    type(complex_3x3_matrix) :: B
    integer(i4), intent(in) :: x(4), mu
    real(dp) :: epsilon = 0.1_dp
    integer :: i
    integer, parameter :: n = 30
    !type(complex_3x3_matrix), dimension(n), intent(out) :: V

    V = U(x(1),x(2),x(3),x(4))%link(mu)
    do i = 1, n
       B = Zeta(U,x,mu)
       B%matrix = B%matrix*epsilon
       U(x(1),x(2),x(3),x(4))%link(mu) = my_exp(B) * U(x(1),x(2),x(3),x(4))%link(mu)
       !V(i) = U(x(1),x(2),x(3),x(4))%link(mu)
    end do
  end subroutine wilson_flow_euler

  function my_exp(W) result(res)
    type(complex_3x3_matrix), intent(in) :: W
    type(complex_3x3_matrix) :: res, B, C
    
    integer, parameter :: n = 3
    integer, parameter :: lda = 3
    integer, parameter :: ldvl  = n
    integer, parameter :: ldvr = n
    integer, parameter :: lwork = 2*n
    
    complex(dp), dimension(n,n) :: A
    complex(dp), dimension(lwork) :: work
    complex(dp), dimension(n) :: eigenv
    complex(dp), dimension(ldvl,n) :: vl
    complex(dp), dimension(ldvr,n) :: vr
    real(dp), dimension(2*n) :: rwork
    integer :: info

    A = W%matrix
    call zgeev('N', 'V', n, A, lda,eigenv, vl, ldvl, vr, ldvr, WORK, lwork, rwork,INFO)

    C%matrix = vr
    B%matrix = one
    B%matrix(1,1) = exp(eigenv(1))
    B%matrix(2,2) = exp(eigenv(2))
    B%matrix(3,3) = exp(eigenv(3))
    
    res = C*B*inv(C)

  end function my_exp

  
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

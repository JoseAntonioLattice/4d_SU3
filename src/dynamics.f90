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
    
    call set_levi_civita

    ! ---- elegir start ----
    ! Para probar el instanton seed, reemplaza cold_start por:
    !   call instanton_start(U, rho=4.0_dp, charge=1_i4)
    ! y pon N_thermalization=0 para correr el flujo directamente.
    call cold_start(U)
    ! call instanton_start(U, rho=4.0_dp, charge=1_i4)

    do i_beta = 1, size(beta)
       call thermalization(U,Lx,Lt,beta(i_beta),N,d,algorithm,N_thermalization)
       call create_measurements_file(Lx,Lt,beta(i_beta),algorithm,equilibrium)
       !call measurements_sweeps(U,Lx,Lt,beta(i_beta),N,d,algorithm,N_measurements,N_skip)
       call wilson_flow_rk4(U)
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

  ! NOTE: the original `energy(U)` function was removed here. It had no
  ! variable declarations at all (Ep, x,y,z,w, Lx, Lt, mu, nu, d, ...),
  ! so it could not compile under `implicit none`, and it was not called
  ! from anywhere else in the code. Its functionality is already covered
  ! by `take_measurements` (Ep) and `action`.
  
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
    ! Z = -T_A(U * Sigma^dagger).
    ! With S = -(beta/N) Re tr(U Sigma^dagger) and the su(3) inner product
    ! (X,Y) = -2 tr(XY) (Luscher App. A.6), the gradient-DESCENT direction
    ! is Z = -T_A(Omega) (the minus sign IS needed - confirmed empirically:
    ! without it, the action increases monotonically toward a maximum
    ! instead of decreasing toward a minimum).
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
    integer(i4) :: x(4), mu, Lx, Lt
    real(dp) :: epsilon = 0.1_dp
    integer :: i, x1, x2, x3, x4
    integer, parameter :: n = 100
    real(dp) :: S

    Lx = size(U(:,1,1,1))
    Lt = size(U(1,1,1,:))

    print*, "inside wilson flow Euler"
    print '(A6,2X,A14,2X,A18,2X,A18)', "# step", "t", "action", "Q"
    call write_flow_observables(U, 0, 0.0_dp, 'euler')
    S = action(U,Lx,Lt,1.0_dp)
    print '(I6,2X,F14.6,2X,F18.10,2X,F18.10)', 0, 0.0_dp, S, topological_charge_clover(U)

    do i = 1, n
       do x1 = 1, Lx
          do x2 = 1, Lx
             do x3 = 1, Lx
                do x4 = 1, Lt
                   x = [x1,x2,x3,x4]
                   do mu = 1, 4
                      B = Zeta(U,x,mu)
                      B%matrix = B%matrix*epsilon
                      V(x1,x2,x3,x4)%link(mu)%matrix = matmul(my_exp(B%matrix), U(x1,x2,x3,x4)%link(mu)%matrix)
                   end do
                end do
             end do
          end do
       end do
       U = V
       S = action(U,Lx,Lt,1.0_dp)
       call write_flow_observables(U, i, i*epsilon, 'euler')
       print '(I6,2X,F14.6,2X,F18.10,2X,F18.10)', i, i*epsilon, S, topological_charge_clover(U)
    end do
  end subroutine wilson_flow_euler

  subroutine wilson_flow_rk4(U)
    ! Luscher's RK4 integrator for the Wilson/gradient flow
    ! (M. Luscher, "Properties and uses of the Wilson flow in lattice QCD",
    ! JHEP 1008:071 (2010), eq. (4.4)):
    !
    !   W0 = V(t)
    !   W1 = exp( (1/4) Z0 ) W0,                       Z0 = eps * Z(W0)
    !   W2 = exp( (8/9) Z1 - (17/36) Z0 ) W1,          Z1 = eps * Z(W1)
    !   W3 = exp( (3/4) Z2 - (8/9) Z1 + (17/36) Z0 ) W2,  Z2 = eps * Z(W2)
    !   V(t+eps) = W3
    !
    ! This is 3rd-order accurate in eps, so a much larger step size can be
    ! used than for wilson_flow_euler (1st order) at the same accuracy.
    type(link_variable), dimension(:,:,:,:), intent(inout) :: U
    type(link_variable), dimension(size(U(:,1,1,1)),size(U(1,:,1,1)),size(U(1,1,:,1)),size(U(1,1,1,:))) :: W, Unew
    type(complex_3x3_matrix), dimension(size(U(:,1,1,1)),size(U(1,:,1,1)),size(U(1,1,:,1)),size(U(1,1,1,:)),4) :: Z0, Z1
    type(complex_3x3_matrix) :: Ztmp, A
    integer(i4) :: x(4), mu, Lx, Lt
    real(dp) :: epsilon = 0.1_dp
    integer :: i, x1, x2, x3, x4
    integer, parameter :: n = 100
    real(dp) :: S

    Lx = size(U(:,1,1,1))
    Lt = size(U(1,1,1,:))

    print*, "inside wilson flow RK4"
    print '(A6,2X,A14,2X,A18,2X,A18)', "# step", "t", "action", "Q"
    call write_flow_observables(U, 0, 0.0_dp, 'rk4')
    S = action(U,Lx,Lt,1.0_dp)
    print '(I6,2X,F14.6,2X,F18.10,2X,F18.10)', 0, 0.0_dp, S, topological_charge_clover(U)

    do i = 1, n

       ! ---- W0 = U ----
       W = U

       ! ---- Z0 = eps * Z(W0) ----
       do x1 = 1, Lx
          do x2 = 1, Lx
             do x3 = 1, Lx
                do x4 = 1, Lt
                   x = [x1,x2,x3,x4]
                   do mu = 1, 4
                      Ztmp = Zeta(W,x,mu)
                      Ztmp%matrix = Ztmp%matrix * epsilon
                      Z0(x1,x2,x3,x4,mu) = Ztmp
                   end do
                end do
             end do
          end do
       end do

       ! ---- W1 = exp(Z0/4) W0  (purely local: each link only needs its own Z0) ----
       do x1 = 1, Lx
          do x2 = 1, Lx
             do x3 = 1, Lx
                do x4 = 1, Lt
                   do mu = 1, 4
                      A%matrix = Z0(x1,x2,x3,x4,mu)%matrix * 0.25_dp
                      W(x1,x2,x3,x4)%link(mu)%matrix = &
                           matmul(my_exp(A%matrix), W(x1,x2,x3,x4)%link(mu)%matrix)
                   end do
                end do
             end do
          end do
       end do

       ! ---- Z1 = eps * Z(W1), evaluated on the (now complete) field W1 ----
       do x1 = 1, Lx
          do x2 = 1, Lx
             do x3 = 1, Lx
                do x4 = 1, Lt
                   x = [x1,x2,x3,x4]
                   do mu = 1, 4
                      Ztmp = Zeta(W,x,mu)
                      Ztmp%matrix = Ztmp%matrix * epsilon
                      Z1(x1,x2,x3,x4,mu) = Ztmp
                   end do
                end do
             end do
          end do
       end do

       ! ---- W2 = exp(8/9 Z1 - 17/36 Z0) W1  (purely local) ----
       do x1 = 1, Lx
          do x2 = 1, Lx
             do x3 = 1, Lx
                do x4 = 1, Lt
                   do mu = 1, 4
                      A%matrix = (8.0_dp/9.0_dp) * Z1(x1,x2,x3,x4,mu)%matrix &
                               - (17.0_dp/36.0_dp) * Z0(x1,x2,x3,x4,mu)%matrix
                      W(x1,x2,x3,x4)%link(mu)%matrix = &
                           matmul(my_exp(A%matrix), W(x1,x2,x3,x4)%link(mu)%matrix)
                   end do
                end do
             end do
          end do
       end do

       ! ---- W3 = exp(3/4 Z2 - 8/9 Z1 + 17/36 Z0) W2,  Z2 = eps * Z(W2) ----
       ! Z2 is computed from W (= W2, not yet modified in this pass) and the
       ! result is written into Unew, so the staples used for Z2 always see
       ! the consistent, unmodified W2 field.
       do x1 = 1, Lx
          do x2 = 1, Lx
             do x3 = 1, Lx
                do x4 = 1, Lt
                   x = [x1,x2,x3,x4]
                   do mu = 1, 4
                      Ztmp = Zeta(W,x,mu)
                      A%matrix = 0.75_dp * Ztmp%matrix * epsilon &
                               - (8.0_dp/9.0_dp)  * Z1(x1,x2,x3,x4,mu)%matrix &
                               + (17.0_dp/36.0_dp) * Z0(x1,x2,x3,x4,mu)%matrix
                      Unew(x1,x2,x3,x4)%link(mu)%matrix = &
                           matmul(my_exp(A%matrix), W(x1,x2,x3,x4)%link(mu)%matrix)
                   end do
                end do
             end do
          end do
       end do

       U = Unew

       S = action(U,Lx,Lt,1.0_dp)
       call write_flow_observables(U, i, i*epsilon, 'rk4')
       print '(I6,2X,F14.6,2X,F18.10,2X,F18.10)', i, i*epsilon, S, topological_charge_clover(U)
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
    type(complex_3x3_matrix) :: Fmunu(4,4)

    ! F(U,x,mu,nu) is expensive (index lookups + ~12 matmuls/daggers).
    ! The naive quadruple loop below calls F(...) twice per iteration,
    ! i.e. up to 4^4 * 2 = 512 evaluations per lattice site, even though
    ! there are only 16 distinct (mu,nu) pairs. Precompute them once.
    do mu = 1, 4
       do nu = 1, 4
          Fmunu(mu,nu) = F(U,x,mu,nu)
       end do
    end do

    topological_density = 0.0_dp
    do mu = 1, 4
       do nu = 1, 4
          do rho = 1, 4
             do sigma = 1, 4
                ! Skip the ~93% of terms where levi_civita = 0
                ! (only 24 of the 256 (mu,nu,rho,sigma) combinations
                ! are nonzero), avoiding an unnecessary 3x3 matmul + tr.
                if (levi_civita(mu,nu,rho,sigma) /= 0) then
                   topological_density = topological_density &
                        - levi_civita(mu,nu,rho,sigma)*tr(Fmunu(mu,nu)*Fmunu(rho,sigma))
                end if
             end do
          end do
       end do
    end do
    !topological_density = -topological_density/(32*pi**2)
  end function topological_density

  ! ============================================================
  ! Standalone field-strength tensor F_{mu nu}(x)
  ! ============================================================
  !
  ! Returns the anti-Hermitian traceless clover-leaf estimate of the
  ! continuum field-strength tensor at site x in the (mu,nu) plane:
  !
  !   F_{mu nu}(x) = (1/8) [ Q(x,mu,nu) - Q^dagger(x,mu,nu) ]
  !
  ! where Q is the sum of the 4 oriented plaquettes (clover) around x
  ! in the (mu,nu) plane.
  !
  ! This function is IDENTICAL in result to the private function F(U,x,mu,nu)
  ! used internally by topological_density, but exposed here as a public
  ! standalone routine so it can be called, measured, and printed
  ! independently of any flow or charge calculation.
  !
  ! Example use:
  !   Fmunu = F_field(U, [x1,x2,x3,x4], 1, 2)   ! F_{12} at (x1,x2,x3,x4)
  !   print*, tr(Fmunu * Fmunu)                    ! action density contribution
  function F_field(U, x, mu, nu) result(Fmunu)
    type(link_variable), dimension(:,:,:,:), intent(in) :: U
    integer(i4), intent(in) :: x(4), mu, nu
    type(complex_3x3_matrix) :: Fmunu, Q, tmp
    integer(i4), dimension(4) :: ipx_mu, ipx_nu, imx_mu, imx_nu, &
         x_im_mu_ip_nu, x_im_mu_im_nu, x_ip_mu_im_nu

    ipx_mu = ip_func(x,mu)
    ipx_nu = ip_func(x,nu)
    imx_mu = im_func(x,mu)
    imx_nu = im_func(x,nu)
    x_im_mu_ip_nu = im_func(ipx_nu,mu)
    x_im_mu_im_nu = im_func(imx_nu,mu)
    x_ip_mu_im_nu = ip_func(imx_nu,mu)

    ! Clover: sum of 4 oriented plaquettes around x in the (mu,nu) plane.
    ! Plaquette 1: forward mu, forward nu  (P_{++})
    Q = U(x(1),x(2),x(3),x(4))%link(mu) * &
        U(ipx_mu(1),ipx_mu(2),ipx_mu(3),ipx_mu(4))%link(nu) * &
        dagger(U(ipx_nu(1),ipx_nu(2),ipx_nu(3),ipx_nu(4))%link(mu)) * &
        dagger(U(x(1),x(2),x(3),x(4))%link(nu))
    ! Plaquette 2: backward mu, forward nu  (P_{-+})
    tmp = U(x(1),x(2),x(3),x(4))%link(nu) * &
          dagger(U(x_im_mu_ip_nu(1),x_im_mu_ip_nu(2),x_im_mu_ip_nu(3),x_im_mu_ip_nu(4))%link(mu)) * &
          dagger(U(imx_mu(1),imx_mu(2),imx_mu(3),imx_mu(4))%link(nu)) * &
          U(imx_mu(1),imx_mu(2),imx_mu(3),imx_mu(4))%link(mu)
    Q%matrix = Q%matrix + tmp%matrix
    ! Plaquette 3: backward mu, backward nu (P_{--})
    tmp = dagger(U(imx_mu(1),imx_mu(2),imx_mu(3),imx_mu(4))%link(mu)) * &
          dagger(U(x_im_mu_im_nu(1),x_im_mu_im_nu(2),x_im_mu_im_nu(3),x_im_mu_im_nu(4))%link(nu)) * &
          U(x_im_mu_im_nu(1),x_im_mu_im_nu(2),x_im_mu_im_nu(3),x_im_mu_im_nu(4))%link(mu) * &
          U(imx_nu(1),imx_nu(2),imx_nu(3),imx_nu(4))%link(nu)
    Q%matrix = Q%matrix + tmp%matrix
    ! Plaquette 4: forward mu, backward nu  (P_{+-})
    tmp = dagger(U(imx_nu(1),imx_nu(2),imx_nu(3),imx_nu(4))%link(nu)) * &
          U(imx_nu(1),imx_nu(2),imx_nu(3),imx_nu(4))%link(mu) * &
          U(x_ip_mu_im_nu(1),x_ip_mu_im_nu(2),x_ip_mu_im_nu(3),x_ip_mu_im_nu(4))%link(nu) * &
          dagger(U(x(1),x(2),x(3),x(4))%link(mu))
    Q%matrix = Q%matrix + tmp%matrix

    Fmunu = (Q - dagger(Q)) / 8.0_dp

  end function F_field

  ! ============================================================
  ! Standalone topological density at a single site x
  ! ============================================================
  !
  ! Returns the (unnormalized) topological density at site x:
  !
  !   q(x) = - sum_{mu,nu,rho,sigma} eps_{mu nu rho sigma}
  !               tr[ F_{mu nu}(x) F_{rho sigma}(x) ]
  !
  ! To get the physical charge density, normalize by 1/(32 pi^2):
  !
  !   q_phys(x) = q(x) / (-32 pi^2)
  !
  ! and then  Q = sum_x q_phys(x).
  !
  ! This is equivalent to the private topological_density(U,x) used
  ! internally by the flow routines, but uses F_field() so it can be
  ! called as a fully standalone measurement.
  !
  ! Example use:
  !   q = topological_density_at(U, [x1,x2,x3,x4])
  !   print*, "q(x) =", real(-q/(32*pi**2), dp)
  function topological_density_at(U, x) result(qx)
    type(link_variable), dimension(:,:,:,:), intent(in) :: U
    integer(i4), intent(in) :: x(4)
    complex(dp) :: qx
    type(complex_3x3_matrix) :: Fmunu(4,4)
    integer(i4) :: mu, nu, rho, sigma

    ! Precompute all 16 distinct (mu,nu) pairs once (see optimization note
    ! in topological_density: avoids 512 evaluations of F_field per site).
    do mu = 1, 4
       do nu = 1, 4
          Fmunu(mu,nu) = F_field(U, x, mu, nu)
       end do
    end do

    qx = 0.0_dp
    do mu = 1, 4
       do nu = 1, 4
          do rho = 1, 4
             do sigma = 1, 4
                if (levi_civita(mu,nu,rho,sigma) /= 0) then
                   qx = qx - levi_civita(mu,nu,rho,sigma) * tr(Fmunu(mu,nu) * Fmunu(rho,sigma))
                end if
             end do
          end do
       end do
    end do

  end function topological_density_at

  ! ============================================================
  ! Standalone topological charge density over the full lattice
  ! ============================================================
  !
  ! Returns a real 4D array q_phys(x1,x2,x3,x4) with the properly
  ! normalized topological charge density at every site:
  !
  !   q_phys(x) = - tr[ F_{mu nu}(x) F_{rho sigma}(x) ] * eps / (32 pi^2)
  !
  ! so that   Q = sum_{x} q_phys(x)   gives the total topological charge.
  !
  ! The full array is useful for:
  !   - Writing the density field to disk for visualization
  !   - Identifying the location and size of individual instantons
  !   - Checking that the charge is localized (instantons/anti-instantons
  !     visible as peaks/troughs) after sufficient gradient flow
  !
  ! Example use:
  !   real(dp), allocatable :: qdens(:,:,:,:)
  !   allocate(qdens(Lx,Lx,Lx,Lt))
  !   qdens = topological_charge_density(U)
  !   print*, "Q =", sum(qdens)
  !   print*, "max q(x) =", maxval(qdens)
  !   print*, "min q(x) =", minval(qdens)
  function topological_charge_density(U) result(qdens)
    type(link_variable), dimension(:,:,:,:), intent(in) :: U
    real(dp), dimension(size(U(:,1,1,1)),size(U(1,:,1,1)),size(U(1,1,:,1)),size(U(1,1,1,:))) :: qdens
    integer(i4) :: x1, x2, x3, x4, Lx, Lt
    integer(i4) :: x(4)

    Lx = size(U(:,1,1,1))
    Lt = size(U(1,1,1,:))

    do x1 = 1, Lx
       do x2 = 1, Lx
          do x3 = 1, Lx
             do x4 = 1, Lt
                x = [x1,x2,x3,x4]
                qdens(x1,x2,x3,x4) = real(-topological_density_at(U,x) / (32*pi**2), dp)
             end do
          end do
       end do
    end do

  end function topological_charge_density

  function topological_charge(U)
    ! Computes the (lattice-discretized) topological charge of the
    ! configuration U, summing topological_density(U,x) over the whole
    ! lattice and applying the standard normalization:
    !
    !   Q = - (1/32 pi^2) * sum_x topological_density(U,x)
    !
    ! This is completely independent of the flow routines: it can be
    ! called on ANY configuration -- a cold start (U = I everywhere,
    ! where it must return exactly 0), a thermalized configuration,
    ! or a configuration at any point during/after a flow -- simply by
    ! passing U. Useful as a standalone measurement / sanity check.
    type(link_variable), dimension(:,:,:,:), intent(in) :: U
    real(dp) :: topological_charge
    complex(dp) :: Qsum
    integer(i4) :: x1, x2, x3, x4, Lx, Lt
    integer(i4) :: x(4)

    Lx = size(U(:,1,1,1))
    Lt = size(U(1,1,1,:))

    Qsum = 0.0_dp
    do x1 = 1, Lx
       do x2 = 1, Lx
          do x3 = 1, Lx
             do x4 = 1, Lt
                x = [x1,x2,x3,x4]
                Qsum = Qsum + topological_density(U,x)
             end do
          end do
       end do
    end do

    topological_charge = real(-Qsum/(32*pi**2), dp)

  end function topological_charge

  ! ============================================================
  ! Raw clover Q_{mu nu}(x) tensor  (ec. 9 del documento)
  ! ============================================================
  !
  ! Devuelve la suma cruda de las 4 plaquetas orientadas alrededor de x
  ! en el plano (mu,nu), SIN el factor 1/8i. Esta es la cantidad Q_{mu nu}
  ! que aparece directamente en la formula optimizada ec. (19):
  !
  !   Q_{mu nu}(x) = P_{++} + P_{-+} + P_{--} + P_{+-}
  !
  ! donde cada P es una plaqueta orientada. En terminos de F^clover:
  !
  !   F^clover_{mu nu}(x) = (1/8i) [ Q_{mu nu}(x) - Q_{nu mu}(x) ]
  !                       = (Q_{mu nu}(x) - Q^dagger_{mu nu}(x)) / 8
  !
  ! Nota: Q_{nu mu}(x) = Q^dagger_{mu nu}(x) para matrices unitarias,
  ! por lo que no hace falta calcularlo por separado.
  function clover_Q(U, x, mu, nu) result(Qmunu)
    type(link_variable), dimension(:,:,:,:), intent(in) :: U
    integer(i4), intent(in) :: x(4), mu, nu
    type(complex_3x3_matrix) :: Qmunu, tmp
    integer(i4), dimension(4) :: ipx_mu, ipx_nu, imx_mu, imx_nu, &
         x_im_mu_ip_nu, x_im_mu_im_nu, x_ip_mu_im_nu

    ipx_mu = ip_func(x, mu)
    ipx_nu = ip_func(x, nu)
    imx_mu = im_func(x, mu)
    imx_nu = im_func(x, nu)
    x_im_mu_ip_nu = im_func(ipx_nu, mu)
    x_im_mu_im_nu = im_func(imx_nu, mu)
    x_ip_mu_im_nu = ip_func(imx_nu, mu)

    ! Plaqueta 1: forward mu, forward nu  (P_{++})
    Qmunu = U(x(1),x(2),x(3),x(4))%link(mu) * &
            U(ipx_mu(1),ipx_mu(2),ipx_mu(3),ipx_mu(4))%link(nu) * &
            dagger(U(ipx_nu(1),ipx_nu(2),ipx_nu(3),ipx_nu(4))%link(mu)) * &
            dagger(U(x(1),x(2),x(3),x(4))%link(nu))
    ! Plaqueta 2: backward mu, forward nu (P_{-+})
    tmp = U(x(1),x(2),x(3),x(4))%link(nu) * &
          dagger(U(x_im_mu_ip_nu(1),x_im_mu_ip_nu(2),x_im_mu_ip_nu(3),x_im_mu_ip_nu(4))%link(mu)) * &
          dagger(U(imx_mu(1),imx_mu(2),imx_mu(3),imx_mu(4))%link(nu)) * &
          U(imx_mu(1),imx_mu(2),imx_mu(3),imx_mu(4))%link(mu)
    Qmunu%matrix = Qmunu%matrix + tmp%matrix
    ! Plaqueta 3: backward mu, backward nu (P_{--})
    tmp = dagger(U(imx_mu(1),imx_mu(2),imx_mu(3),imx_mu(4))%link(mu)) * &
          dagger(U(x_im_mu_im_nu(1),x_im_mu_im_nu(2),x_im_mu_im_nu(3),x_im_mu_im_nu(4))%link(nu)) * &
          U(x_im_mu_im_nu(1),x_im_mu_im_nu(2),x_im_mu_im_nu(3),x_im_mu_im_nu(4))%link(mu) * &
          U(imx_nu(1),imx_nu(2),imx_nu(3),imx_nu(4))%link(nu)
    Qmunu%matrix = Qmunu%matrix + tmp%matrix
    ! Plaqueta 4: forward mu, backward nu  (P_{+-})
    tmp = dagger(U(imx_nu(1),imx_nu(2),imx_nu(3),imx_nu(4))%link(nu)) * &
          U(imx_nu(1),imx_nu(2),imx_nu(3),imx_nu(4))%link(mu) * &
          U(x_ip_mu_im_nu(1),x_ip_mu_im_nu(2),x_ip_mu_im_nu(3),x_ip_mu_im_nu(4))%link(nu) * &
          dagger(U(x(1),x(2),x(3),x(4))%link(mu))
    Qmunu%matrix = Qmunu%matrix + tmp%matrix

  end function clover_Q

  ! ============================================================
  ! Densidad topologica optimizada (ec. 19 del documento)
  ! ============================================================
  !
  ! Implementa la formula reducida que explota las 3 simetrias:
  !
  !   eps_{mu nu rho sigma} tr(F_{mn} F_{rs}) = eps_{nu mu rho sigma} tr(F_{nm} F_{rs})
  !   eps_{mu nu rho sigma} tr(F_{mn} F_{rs}) = eps_{mu nu sigma rho} tr(F_{mn} F_{sr})
  !   eps_{mu nu rho sigma} tr(F_{mn} F_{rs}) = eps_{rho sigma mu nu} tr(F_{rs} F_{mn})
  !
  ! reduciendo los 24 terminos no-nulos del Levi-Civita a solo 3:
  !
  !   q(x) = - (1/128 pi^2) * sum_{[mnrs] in {[1234],[1324],[1423]}}
  !               eps_{mnrs} tr( Q_{mn}(x) [Q_{rs}(x) - Q^dag_{rs}(x)] )
  !
  ! Relacion con tr(F F): tr(F^cl_{mn} F^cl_{rs}) = (1/32) tr(Q_{mn}(Q_{rs}-Q^dag_{rs}))
  ! (ec. 18), y la suma reducida equivale a sumar los 24 terminos completos,
  ! con el factor total 1/(4 * 32 * pi^2) = 1/(128 pi^2).
  !
  ! Las 3 permutaciones independientes y sus signos (eps):
  !   [1,2,3,4] -> eps = +1
  !   [1,3,2,4] -> eps = -1   (intercambio de 2 y 3)
  !   [1,4,2,3] -> eps = +1   (transposicion ciclica)
  function topological_density_clover(U, x) result(qx)
    type(link_variable), dimension(:,:,:,:), intent(in) :: U
    integer(i4), intent(in) :: x(4)
    real(dp) :: qx
    type(complex_3x3_matrix) :: Q12, Q34, Q13, Q24, Q14, Q23
    complex(dp) :: term

    ! Solo calculamos los 6 clovers que aparecen en las 3 permutaciones
    Q12 = clover_Q(U, x, 1, 2)
    Q34 = clover_Q(U, x, 3, 4)
    Q13 = clover_Q(U, x, 1, 3)
    Q24 = clover_Q(U, x, 2, 4)
    Q14 = clover_Q(U, x, 1, 4)
    Q23 = clover_Q(U, x, 2, 3)

    ! [1,2,3,4]: eps_{1234} = +1
    !   tr( Q12 * (Q34 - Q34^dag) )
    term = tr( Q12 * (Q34 - dagger(Q34)) )

    ! [1,3,2,4]: eps_{1324} = -1  (intercambio de indices 2<->3)
    !   - tr( Q13 * (Q24 - Q24^dag) )
    term = term - tr( Q13 * (Q24 - dagger(Q24)) )

    ! [1,4,2,3]: eps_{1423} = +1  (permutacion ciclica 2->4->3->2)
    !   tr( Q14 * (Q23 - Q23^dag) )
    term = term + tr( Q14 * (Q23 - dagger(Q23)) )

    qx = real(-term / (128.0_dp * pi**2), dp)

  end function topological_density_clover

  ! ============================================================
  ! Carga topologica optimizada (suma sobre todo el lattice)
  ! ============================================================
  !
  ! Versión optimizada de topological_charge usando la formula reducida
  ! de 3 terminos (ec. 19). Es algebraicamente identica al resultado
  ! de topological_charge pero ~8x mas rapida, ya que:
  !   - Solo calcula 6 clovers por sitio (en vez de 16)
  !   - Solo evalua 3 trazas de matrices 3x3 (en vez de 24)
  !
  ! Uso:
  !   Q = topological_charge_clover(U)
  !   print*, "Q =", Q
  function topological_charge_clover(U) result(Q)
    type(link_variable), dimension(:,:,:,:), intent(in) :: U
    real(dp) :: Q
    integer(i4) :: x1, x2, x3, x4, Lx, Lt
    integer(i4) :: x(4)

    Lx = size(U(:,1,1,1))
    Lt = size(U(1,1,1,:))

    Q = 0.0_dp
    do x1 = 1, Lx
       do x2 = 1, Lx
          do x3 = 1, Lx
             do x4 = 1, Lt
                x = [x1,x2,x3,x4]
                Q = Q + topological_density_clover(U, x)
             end do
          end do
       end do
    end do

  end function topological_charge_clover

  ! ============================================================
  ! Escribe accion y carga topologica a archivo en data/
  ! ============================================================
  !
  ! Guarda en data/gradient_flow_<label>.dat una tabla con columnas:
  !   step   t=step*epsilon   action   Q
  !
  ! Uso desde fuera del modulo:
  !   call write_flow_observables(U, 0, 0.0_dp, 'rk4')
  !   ! (despues de cada paso del flujo)
  !   call write_flow_observables(U, i, i*epsilon, 'rk4')
  !
  ! El archivo se crea/sobreescribe en el primer llamado (step=0)
  ! y se agrega una linea en los llamados posteriores.
  subroutine write_flow_observables(U, step, t, label)
    type(link_variable), dimension(:,:,:,:), intent(in) :: U
    integer, intent(in) :: step
    real(dp), intent(in) :: t
    character(*), intent(in) :: label
    integer, parameter :: iunit = 42
    integer :: Lx, Lt
    real(dp) :: S, Q
    character(len=256) :: filename

    Lx = size(U(:,1,1,1))
    Lt = size(U(1,1,1,:))
    S  = action(U, Lx, Lt, 1.0_dp)
    Q  = topological_charge_clover(U)

    write(filename, '(A,A,A)') 'data/gradient_flow_', trim(label), '.dat'

    if (step == 0) then
       ! Primer paso: crear/sobreescribir el archivo y escribir cabecera
       open(unit=iunit, file=trim(filename), status='replace', action='write')
       write(iunit, '(A6,3X,A14,3X,A18,3X,A18)') &
            '# step', 't', 'action', 'Q'
    else
       open(unit=iunit, file=trim(filename), status='old', &
            action='write', position='append')
    end if

    write(iunit, '(I6,3X,F14.6,3X,F18.10,3X,F18.10)') step, t, S, Q
    close(iunit)

  end subroutine write_flow_observables

  ! ============================================================
  ! Instanton seed: configura un (anti-)instanton BPST exacto
  ! embebido en el bloque SU(2) superior izquierdo de SU(3)
  ! ============================================================
  !
  ! Construye los links U_mu(x) a partir del campo de gauge continuo
  ! del instanton BPST discretizado en el punto medio del link:
  !
  !   U_mu(x) = exp( i Gmu(x + mu_hat/2) )
  !
  ! donde el potencial de gauge del instanton (en gauge regular) es:
  !
  !   Gmu^a(x) = 2 * eta^a_{mu nu} * (x-x0)_nu / [(x-x0)^2 + rho^2]
  !
  ! con eta^a_{mu nu} el tensor de 't Hooft (BPST 1975):
  !   eta^a_{mu nu} = eps_{a mu nu}         para mu,nu en {1,2,3}
  !   eta^a_{mu  4} = +delta_{a,mu}         para mu en {1,2,3}
  !   eta^a_{4  nu} = -delta_{a,nu}         para nu en {1,2,3}
  !   eta^a_{mu nu} = -eta^a_{nu mu}        (antisimetrico)
  !
  ! Para anti-instanton (Q=-1): eta -> eta_bar con signo opuesto
  ! en los terminos que involucran la direccion 4.
  !
  ! El link U_mu(x) se obtiene exponenciando exactamente el algebra su(2)
  ! usando la formula cerrada exp(i n_a sigma_a/2) y luego embebiendo
  ! en SU(3) via el bloque 2x2 superior izquierdo.
  !
  ! Argumentos:
  !   U      : campo gauge a inicializar
  !   rho    : tamano del instanton en unidades de a (recomendado: 4.0 - 6.0)
  !   charge : +1 (instanton, Q=+1) o -1 (anti-instanton, Q=-1)
  subroutine instanton_start(U, rho, charge)
    type(link_variable), intent(inout), dimension(:,:,:,:) :: U
    real(dp),    intent(in) :: rho
    integer(i4), intent(in) :: charge

    integer(i4) :: Lx, Lt, ix1, ix2, ix3, ix4, imu, inu, ic
    integer(i4) :: xvec(4)
    real(dp)    :: x0(4)
    real(dp)    :: xmid(4), dxvec(4), r2
    ! Potencial de gauge Gmu^a (ic=color, imu=direccion): renombrado de A para evitar
    ! conflicto con el indice 'a' (Fortran es case-insensitive)
    real(dp)    :: Gfield(4,3)
    ! Tensor de 't Hooft: eta(mu,nu,color)
    real(dp)    :: eta(4,4,3)
    real(dp)    :: nvec, c0, c1
    complex(dp) :: E2(2,2)

    Lx = size(U(:,1,1,1))
    Lt = size(U(1,1,1,:))

    ! Centro geometrico de la red
    x0 = [0.5_dp*(Lx+1), 0.5_dp*(Lx+1), 0.5_dp*(Lx+1), 0.5_dp*(Lt+1)]

    ! ---- Tabla completa del tensor de 't Hooft eta^a_{mu nu} ----
    ! Indice de almacenamiento: eta(mu, nu, color_a)
    ! Definicion BPST: antisimetrico, eta^a_{mu nu} = eps_{a,mu,nu} para mu,nu<4
    !                  eta^a_{mu,4} = +delta_{a,mu},  eta^a_{4,nu} = -delta_{a,nu}
    eta = 0.0_dp

    ! -- Parte espacial: eps_{a,mu,nu} para mu,nu in {1,2,3} --
    eta(2,3,1) = +1.0_dp;  eta(3,2,1) = -1.0_dp
    eta(3,1,2) = +1.0_dp;  eta(1,3,2) = -1.0_dp
    eta(1,2,3) = +1.0_dp;  eta(2,1,3) = -1.0_dp

    ! Nota de convencion: con el estimador clover de este codigo,
    ! eta_bar (antidual, signo temporal negativo) da Q=+1,
    ! y eta (autodual, signo temporal positivo) da Q=-1.
    ! Esto es consistente con el signo del estimador clover usado
    ! en topological_density_clover.
    !
    ! charge=+1 (instanton): usar eta_bar -> signo temporal NEGATIVO
    ! charge=-1 (anti-instanton): usar eta  -> signo temporal POSITIVO
    if (charge == +1_i4) then
      eta(1,4,1) = -1.0_dp;  eta(4,1,1) = +1.0_dp
      eta(2,4,2) = -1.0_dp;  eta(4,2,2) = +1.0_dp
      eta(3,4,3) = -1.0_dp;  eta(4,3,3) = +1.0_dp
    else
      eta(1,4,1) = +1.0_dp;  eta(4,1,1) = -1.0_dp
      eta(2,4,2) = +1.0_dp;  eta(4,2,2) = -1.0_dp
      eta(3,4,3) = +1.0_dp;  eta(4,3,3) = -1.0_dp
    end if

    ! ---- Construir los links ----
    do ix4 = 1, Lt
      do ix3 = 1, Lx
        do ix2 = 1, Lx
          do ix1 = 1, Lx
            xvec = [ix1, ix2, ix3, ix4]

            do imu = 1, 4

              ! Punto medio del link (en unidades de a)
              xmid      = real(xvec, dp)
              xmid(imu) = xmid(imu) + 0.5_dp

              ! Desplazamiento al centro con imagen periodica mas cercana
              dxvec(1) = xmid(1) - x0(1);  dxvec(1) = dxvec(1) - Lx*nint(dxvec(1)/Lx)
              dxvec(2) = xmid(2) - x0(2);  dxvec(2) = dxvec(2) - Lx*nint(dxvec(2)/Lx)
              dxvec(3) = xmid(3) - x0(3);  dxvec(3) = dxvec(3) - Lx*nint(dxvec(3)/Lx)
              dxvec(4) = xmid(4) - x0(4);  dxvec(4) = dxvec(4) - Lt*nint(dxvec(4)/Lt)

              r2 = sum(dxvec**2) + rho**2

              ! Potencial de gauge: Gfield(imu, ic) = 2 * sum_nu eta^ic_{imu,nu} * dxvec(nu) / r2
              do ic = 1, 3
                Gfield(imu, ic) = 0.0_dp
                do inu = 1, 4
                  Gfield(imu, ic) = Gfield(imu, ic) + eta(imu, inu, ic) * dxvec(inu)
                end do
                Gfield(imu, ic) = Gfield(imu, ic) * 2.0_dp / r2
              end do

              ! Exponenciar el algebra su(2) exactamente:
              ! X = i * Gfield^a * sigma_a/2  =>  exp(X) = c0*I + i*c1*(G^a sigma_a)
              ! donde c0 = cos(|G|/2),  c1 = sin(|G|/2)/|G|,  |G| = sqrt(G^a G^a)
              nvec = sqrt(Gfield(imu,1)**2 + Gfield(imu,2)**2 + Gfield(imu,3)**2)

              if (nvec < 1.0e-14_dp) then
                E2(1,1) = (1.0_dp, 0.0_dp);  E2(1,2) = (0.0_dp, 0.0_dp)
                E2(2,1) = (0.0_dp, 0.0_dp);  E2(2,2) = (1.0_dp, 0.0_dp)
              else
                c0 = cos(0.5_dp * nvec)
                c1 = sin(0.5_dp * nvec) / nvec
                ! exp(X) con X = i*G^a*sigma_a/2:
                ! sigma_1 = [[0,1],[1,0]], sigma_2 = [[0,-i],[i,0]], sigma_3 = [[1,0],[0,-1]]
                ! => (1,1): c0 + i*c1*G^3
                !    (1,2): i*c1*G^1 + c1*G^2  (= c1*(G^2 + i*G^1))
                !    (2,1): i*c1*G^1 - c1*G^2  (= c1*(-G^2 + i*G^1))
                !    (2,2): c0 - i*c1*G^3
                E2(1,1) = cmplx( c0,              c1*Gfield(imu,3), dp)
                E2(1,2) = cmplx( c1*Gfield(imu,2), c1*Gfield(imu,1), dp)
                E2(2,1) = cmplx(-c1*Gfield(imu,2), c1*Gfield(imu,1), dp)
                E2(2,2) = cmplx( c0,             -c1*Gfield(imu,3), dp)
              end if

              ! Embeber SU(2) en SU(3): bloque 2x2 superior + (3,3) = 1
              U(ix1,ix2,ix3,ix4)%link(imu)%matrix         = (0.0_dp, 0.0_dp)
              U(ix1,ix2,ix3,ix4)%link(imu)%matrix(1:2,1:2) = E2
              U(ix1,ix2,ix3,ix4)%link(imu)%matrix(3,3)     = (1.0_dp, 0.0_dp)

            end do  ! imu
          end do
        end do
      end do
    end do

    write(*,'(A,F6.2,A,I2,A,I0,A,I0)') &
      '  instanton_start: rho/a=', rho, '  Q=', charge, &
      '  red=', Lx, '^3 x ', Lt

  end subroutine instanton_start

  
end module dynamics

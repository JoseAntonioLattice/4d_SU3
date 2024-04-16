module dynamics

  use iso_fortran_env, only: dp => real64, i4 => int32
  use data_types_observables, only : link_variable, complex_3x3_matrix
  use matrix_operations
  use local_update_algorithms
  use periodic_boundary_conditions_mod
  use get_index_mod
  use create_files
  
  implicit none

  !private !:: dp, i4, link_variable
  !public :: sweeps, create_update, sgn, drand, take_measurements, dagger, tr, gauge_transformation, DS

contains

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
       !if(mod(i,10) == 0) call normalization(U,L)
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
     real(dp) :: correlation_polyakov_loop(Lx/2-1)
     integer(i4) :: i

     do i = 1, N_measurements*N_skip
        call sweeps(U,Lx,Lt,beta,N,d,algorithm)
        if( mod(i,N_skip) == 0)then
           call take_measurements(U,Lx,Lt,E_p,avr_polyakov_loop,correlation_polyakov_loop)
           write(100,*) E_p,avr_polyakov_loop%re, avr_polyakov_loop%im,correlation_polyakov_loop 
        end if
        !if(mod(i,10) == 0) call normalization(U,L)
     end do 

   end subroutine measurements_sweeps

  subroutine sweeps(U,Lx,Lt,beta,N,d,algorithm)
    type(link_variable), intent(inout), dimension(:,:,:,:) :: U
    integer(i4), intent(in)  :: Lx,Lt, N, d
    real(dp), intent(in) :: beta
    character(*), intent(in) :: algorithm
    type(complex_3x3_matrix) :: Up
    integer(i4) :: x, y,z,w, mu
    real(dp) :: Delta_S

    do x = 1, Lx
       do y = 1, Lx
          do z = 1, Lx
             do w = 1, Lt
                do mu = 1, d
                   !Delta_S = DS(U,mu,Up,[x,y,z,w],beta,N)
                   !call metropolis(Delta_S,U(x,y,z,w)%link(mu)%matrix,Up%matrix)
                   !print*, x,y,z,w,mu
                   call heatbath(U,[x,y,z,w],mu,beta)
                end do
             end do
          end do
       end do
    end do
  end subroutine sweeps


  subroutine take_measurements(U,Lx,Lt,Ep,avr_polyakov_loop,correlation_polyakov_loop)
    type(link_variable), dimension(:,:,:,:), intent(in) :: U
    integer(i4), intent(in) ::  Lx,Lt
    real(dp), intent(out) :: Ep
    complex(dp) :: polyakov_loop_array(Lx,Lx,Lx)
    complex(dp), intent(out) :: avr_polyakov_loop
    real(dp):: correlation_polyakov_loop(Lx/2-1)
    integer(i4) :: x,y,z,w,t,mu,nu, xp, yp, zp
    integer(i4), parameter :: d = 4, number_of_planes = d*(d-1)/2
    
    
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
    correlation_polyakov_loop = 0.0_dp

    do t = 1, Lx/2 - 1
       do x = 1, Lx
          xp = mod(x+t,Lx); if(xp == 0) xp = Lx
          do y = 1, Lx
             yp = mod(y+t,Lx); if(yp == 0) yp = Lx
             do z = 1, Lx
                zp = mod(z+t,Lx); if(zp == 0) zp = Lx
                                
                correlation_polyakov_loop(t) = correlation_polyakov_loop(t) + &
                     real(polyakov_loop_array(x,y,z) * conjg( polyakov_loop_array(xp,y,z) + &
                     polyakov_loop_array(x,yp,z) + &
                     polyakov_loop_array(x,y,zp)))
             end do
          end do
       end do
    end do
    
    Ep =  Ep/(3*number_of_planes*Lx**3*Lt)
    !print*, ep!, inv_polyakov_loop_array(1,1,1), conjg(polyakov_loop_array(1,1,1))
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

end module dynamics

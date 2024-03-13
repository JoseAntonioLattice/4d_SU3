module dynamics

  use iso_fortran_env, only: dp => real64, i4 => int32
  use data_types_observables, only : link_variable, complex_3x3_matrix
  use matrix_operations
  use local_update_algorithms
  use periodic_boundary_conditions_mod
  use get_index_mod

  implicit none

  !private !:: dp, i4, link_variable
  !public :: sweeps, create_update, sgn, drand, take_measurements, dagger, tr, gauge_transformation, DS

contains

  subroutine equilibrium_dynamics(U,L,beta,N,d,algorithm,N_thermalization,N_measurements,N_skip)
    use starts
    type(link_variable), intent(inout), dimension(:,:,:,:) :: U
    integer(i4), intent(in)  :: L, N, d
    real(dp), intent(in), dimension(:) :: beta
    character(*), intent(in) :: algorithm
    integer(i4), intent(in) :: N_thermalization, N_measurements, N_skip
    
    integer(i4) :: i_beta
    
    allocate(E_p%array(N_measurements))
    
    call hot_start(U)

    do i_beta = 1, size(beta)
       call thermalization(U,L,beta(i_beta),N,d,algorithm,N_thermalization)
       call measurements_sweeps(U,L,beta(i_beta),N,d,algorithm,N_measurements,N_skip,E_p%array)
       !write(100,*) E_p%array(i_beta)
       !print*, E_p%array
    end do
  end subroutine equilibrium_dynamics
  
  subroutine thermalization(U,L,beta,N,d,algorithm,N_thermalization)
    type(link_variable), intent(inout), dimension(:,:,:,:) :: U
    integer(i4), intent(in)  :: L, N, d
    real(dp), intent(in) :: beta
    character(*), intent(in) :: algorithm
    integer(i4), intent(in) :: N_thermalization

    integer(i4) :: i

    do i = 1, N_thermalization
       call sweeps(U,L,beta,N,d,algorithm)
    end do
   end subroutine thermalization


   subroutine measurements_sweeps(U,L,beta,N,d,algorithm,N_measurements,N_skip,E_p)
     type(link_variable), intent(inout), dimension(:,:,:,:) :: U
     integer(i4), intent(in)  :: L, N, d
     real(dp), intent(in) :: beta
     character(*), intent(in) :: algorithm
     integer(i4), intent(in) :: N_measurements, N_skip
     real(dp), intent(out) :: E_p(:)
     integer(i4) :: i

     do i = 1, N_measurements*N_skip
        call sweeps(U,L,beta,N,d,algorithm)
        if( mod(i,N_skip) == 0)then
           E_p(i/N_skip) = action(U,L,-1.0_dp/N)/(L**d)
        end if
     end do

   end subroutine measurements_sweeps

  subroutine sweeps(U,L,beta,N,d,algorithm)
    type(link_variable), intent(inout), dimension(:,:,:,:) :: U
    integer(i4), intent(in)  :: L, N, d
    real(dp), intent(in) :: beta
    character(*), intent(in) :: algorithm
    type(complex_3x3_matrix) :: Up
    integer(i4) :: x, y,z,w, mu
    real(dp) :: Delta_S

    do x = 1, L
       do y = 1, L
          do z = 1, L
             do w = 1, L
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

  function action(U,L,beta_N)

    type(link_variable), dimension(:,:,:,:), intent(in) :: U
    integer(i4), intent(in) ::  L
    real(dp) :: action, beta_N
    integer(i4) :: x,y,z,w,mu,nu
    integer(i4), parameter :: d = 4, number_of_planes = d*(d-1)/2
    
    
    action = 0.0_dp
    do x = 1, L
       do y = 1, L
          do z = 1, L
             do w = 1, L
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

 function plaquette(U,x,mu,nu)
   type(link_variable), dimension(:,:,:,:), intent(in) :: U
    integer(i4), intent(in) :: x(4), mu, nu
    type(complex_3x3_matrix) :: plaquette
    integer(i4) :: L
    integer(i4), dimension(4) :: ipx_mu, ipx_nu
    L = size(U(:,1,1,1))
    !             x, mu           x + mu, nu                   x + nu, mu                    x, nu

    ipx_mu = ip_func(x,mu)
    ipx_nu = ip_func(x,nu)
    plaquette = U(x(1),x(2),x(3),x(4))%link(mu) * U(ipx_mu(1),ipx_mu(2),ipx_mu(3),ipx_mu(4))%link(nu) * &
         dagger(U(ipx_nu(1),ipx_nu(2),ipx_nu(3),ipx_nu(4))%link(mu)) * dagger(U(x(1),x(2),x(3),x(4))%link(nu))
  end function plaquette

end module dynamics

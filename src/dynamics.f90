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

  subroutine equilibrium_dynamics(U,L,beta,N,d,algorithm,N_thermalization,N_measurements,N_skip,equilibrium)
    use starts
    use statistics
    type(link_variable), intent(inout), dimension(:,:,:,:) :: U
    integer(i4), intent(in)  :: L, N, d
    real(dp), intent(in), dimension(:) :: beta
    character(*), intent(in) :: algorithm
    integer(i4), intent(in) :: N_thermalization, N_measurements, N_skip
    logical, intent(in) :: equilibrium
    integer(i4) :: i_beta
    
    !allocate(E_p%array(N_measurements))
    
    call hot_start(U)

    do i_beta = 1, size(beta)
       call thermalization(U,L,beta(i_beta),N,d,algorithm,N_thermalization)
       call create_measurements_file(L,beta(i_beta),algorithm,equilibrium)
       call measurements_sweeps(U,L,beta(i_beta),N,d,algorithm,N_measurements,N_skip)
       
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


   subroutine measurements_sweeps(U,L,beta,N,d,algorithm,N_measurements,N_skip)
     type(link_variable), intent(inout), dimension(:,:,:,:) :: U
     integer(i4), intent(in)  :: L, N, d
     real(dp), intent(in) :: beta
     character(*), intent(in) :: algorithm
     integer(i4), intent(in) :: N_measurements, N_skip
     real(dp) :: E_p
     complex(dp) :: correlation_polyakov_loop(0:L-1)
     integer(i4) :: i

     do i = 1, N_measurements*N_skip
        call sweeps(U,L,beta,N,d,algorithm)
        if( mod(i,N_skip) == 0)then
           call take_measurements(U,L,E_p,correlation_polyakov_loop)
           write(100,*) E_p,correlation_polyakov_loop
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


  subroutine take_measurements(U,L,Ep,correlation_polyakov_loop)
    type(link_variable), dimension(:,:,:,:), intent(in) :: U
    integer(i4), intent(in) ::  L
    real(dp), intent(out) :: Ep
    complex(dp) :: correlation_polyakov_loop(0:L-1)
    complex(dp) :: poly_loop1
    integer(i4) :: x,y,z,w,mu,nu
    integer(i4), parameter :: d = 4, number_of_planes = d*(d-1)/2
    
    
    Ep = 0.0_dp
    poly_loop1 = polyakov_loop(U,[1,1,1],L)
    do x = 1, L
       correlation_polyakov_loop(x-1) = poly_loop1*conjg(polyakov_loop(U,[1,1,x],L))
       do y = 1, L
          do z = 1, L
             do w = 1, L
                do mu = 1, d - 1
                   do nu = mu + 1, d
                      Ep = Ep + real(tr(plaquette(U,[x,y,z,w],mu,nu)),dp)
                   end do
                end do
             end do
          end do
       end do
    end do
    !print*, correlation_polyakov_loop(0)
    Ep =  Ep/(4*number_of_planes*L**4)

  end subroutine take_measurements

  function polyakov_loop(U,x,L)
    type(link_variable), dimension(:,:,:,:), intent(in) :: U
    integer(i4), dimension(3), intent(in) :: x 
    integer(i4), intent(in) :: L
    type(complex_3x3_matrix) :: product
    complex(dp) :: polyakov_loop
    integer(i4) :: t

    product%matrix = 0.0_dp
    product%matrix(1,1) = 1.0_dp
    product%matrix(2,2) = 1.0_dp
    product%matrix(3,3) = 1.0_dp

    do t = 1, L
       product = product*U(x(1),x(2),x(3),t)%link(4)
    end do

    polyakov_loop = tr(product)
    
  end function polyakov_loop

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

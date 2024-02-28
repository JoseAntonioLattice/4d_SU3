module dynamics

  use iso_fortran_env, only: dp => real64, i4 => int32
  use data_types_observables, only : link_variable, complex_3x3_matrix
  use matrix_operations
  !use local_update_algorithms
  use periodic_boundary_conditions_mod
  !use get_index_mod

  implicit none

  !private !:: dp, i4, link_variable
  !public :: sweeps, create_update, sgn, drand, take_measurements, dagger, tr, gauge_transformation, DS

contains

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
           !E_p(i/N_skip) = action(U,-1.0_dp/N,d)/(L**d)
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
                  ! Delta_S = (beta/N) * DS(U,mu,Up,[x,y,z,w])
                  ! call metropolis(Delta_S,U(x,y,z,w)%link(mu)%matrix,Up%matrix)
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
                      !print*, "Inside action. Inside loop", x, mu, nu
                      action = action + real(tr(plaquette(U,[x,y,z,w],mu,nu)),dp)
                   end do
                end do
             end do
          end do
       end do
    end do
    !print*, "Inside action. outside loop"
    action =  - beta_N * action/number_of_planes
  end function action

  !subroutine gauge_transformation(U)


    !use periodic_boundary_contidions_mod, only : ip
    !use parameters, only : L

   ! type(link_variable), dimension(L), intent(inout) :: U
  !  type(complex_2x2_matrix), dimension(L,L) :: V

 !   integer(i4) :: x, y

!    do x = 1, L
    !   do y = 1, L
    !      call create_update(V(x,y))
    !   end do
    !end do

   ! do x = 1, L
   !    do y = 1, L
    !      U(x,y)%link(1) = V(x,y)*U(x,y)%link(1)*dagger(V(ip(x),y))
   !       U(x,y)%link(2) = V(x,y)*U(x,y)%link(2)*dagger(V(x,ip(y)))
   !    end do
   ! end do

  !end subroutine gauge_transformation


 function plaquette(U,x,mu,nu)
   type(link_variable), dimension(:,:,:,:), intent(in) :: U
    integer(i4), intent(in) :: x(4), mu, nu
    type(complex_3x3_matrix) :: plaquette
    integer(i4) :: L
    integer(i4), dimension(4) :: ipx_mu, ipx_nu
    L = size(U(:,1,1,1))
    !             x, mu           x + mu, nu                   x + nu, mu                    x, nu
    !print*, "Inside plaquette"
    ipx_mu = ip_func(x,mu)
    ipx_nu = ip_func(x,nu)
    plaquette = U(x(1),x(2),x(3),x(4))%link(mu) * U(ipx_mu(1),ipx_mu(2),ipx_mu(3),ipx_mu(4))%link(nu) * &
         dagger(U(ipx_nu(1),ipx_nu(2),ipx_nu(3),ipx_nu(4))%link(mu)) * dagger(U(x(1),x(2),x(3),x(4))%link(nu))
  end function plaquette

!  function DS2(U,x,mu,beta_N,d,Up)
!    type(link_variable), dimension(:,:), intent(inout) :: U
!    integer(i4), intent(in) :: x, mu, d
!    real(dp), intent(in) :: beta_N
    !type(complex_2x2_matrix), intent(out) :: Up
   ! type(complex_2x2_matrix) :: Uold
  !  real(dp) :: DS2, Sold, Snew

 !   call create_unbiased_update(Up)

!    Uold = U(x)%link(mu)

    !Sold = action(U,beta_N,d)
    !U(x)%link(mu) = Up
    !Snew = action(U,beta_N,d)
    !U(x)%link(mu) = Uold

   ! DS2 = Snew - Sold

  !end function DS2

end module dynamics

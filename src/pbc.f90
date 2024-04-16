module periodic_boundary_conditions_mod

  use iso_fortran_env, only : i4 => int32

  implicit none

  integer, allocatable, dimension(:) :: ip, im, ip_t, im_t

contains

  subroutine set_periodic_bounds(L,Lt)

    integer(i4), intent(in) :: L, Lt

    allocate(ip(L),im(L),ip_t(L),im_t(L))
    call initialize(L,Lt)

  end subroutine set_periodic_bounds

  subroutine initialize(L,Lt)

    integer(i4), intent(in) :: L,Lt
    integer(i4) :: i

    do i = 1, L
       ip(i) = i + 1
       im(i) = i - 1
    end do
    ip(L) = 1
    im(1) = L

    
    do i = 1, Lt
       ip_t(i) = i + 1
       im_t(i) = i - 1
    end do
    ip_t(Lt) = 1
    im_t(1) = Lt
    
  end subroutine initialize

  function ip_func(vector,mu)
    integer, dimension(:), intent(in) :: vector
    integer, intent(in) :: mu

    integer, dimension(size(vector)) :: ip_func


    ip_func = vector

    ip_func(mu) = ip(vector(mu))
    if (mu == 4) ip_func(mu) = ip_t(vector(mu))
    
  end function ip_func


  function im_func(vector,mu)

    integer, dimension(:), intent(in) :: vector
    integer, intent(in) :: mu

    integer, dimension(size(vector)) :: im_func


    im_func = vector

    im_func(mu) = im(vector(mu))
    if (mu == 4) im_func(mu) = im_t(vector(mu))
    
  end function im_func

end module periodic_boundary_conditions_mod

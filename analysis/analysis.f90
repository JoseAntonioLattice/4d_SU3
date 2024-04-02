program analysis

  use iso_fortran_env, only : dp => real64, i4 => int32
  use statistics
  use number2string_mod
  implicit none

  integer(i4) :: L
  integer(i4) :: N_thermalization
  integer(i4) :: N_measurements
  integer(i4) :: N_skip
  character(25) :: algorithm
  logical :: equilibrium
  
  namelist /input_parameters/ L, N_thermalization, N_measurements, N_skip, algorithm, equilibrium

  integer(i4) :: inunit
  character(256) :: data_file, data 
  
  type observable
     real(dp), allocatable :: array(:)
     real(dp) :: avr
     real(dp) :: err
  end type observable

  type(observable) :: Ep

  integer(i4) :: i, i_beta, bins
  real(dp), allocatable, dimension(:) :: beta, auto_correlation

  beta = [(i*0.1_dp,i=1,80)]
  
  call read_input()


  allocate(Ep%array(N_measurements),auto_correlation(N_measurements))

  data = "data/L="//trim(int2str(L))//"_equilibrium_"//trim(algorithm)//".dat"
  open(unit = 666, file = trim(data), status = "unknown")
  do i_beta = 1,size(beta)
     data_file = "data/L="//trim(int2str(L))//"/"//"equilibrium"//"/"//trim(algorithm)&
          //"/beta="//real2str(beta(i_beta))//"/observables.dat"
     !print*, trim(data_file)
     open( newunit = inunit, file = trim(data_file) )
     
     do i = 1, N_measurements
        read(inunit,*) Ep%array(i)
     end do
     close(inunit)
     
     open(unit = 69, file = "data/L="//trim(int2str(L))//"/"//"equilibrium"//"/"//trim(algorithm)&
          //"/beta="//real2str(beta(i_beta))//"/auto_correlation.dat")
     auto_correlation = acf(Ep%array)
     do i = 1, size(Ep%array)
        write(69,*) auto_correlation(i)
     end do
     
     call max_jackknife_error_2(Ep%array,Ep%avr,Ep%err,bins)
     write(666,*) beta(i_beta),Ep%avr, Ep%err,bins
  end do
contains

  function acf(x)
    real(dp), intent(in), dimension(:) :: x
    real(dp), dimension(0:size(x)-1) :: acf
    integer(i4) :: k, n, t
    real(dp) :: avr, var, suma
    
    n = size(x)

    avr = sum(x)/n
    var = (sum(x**2)-n*avr**2)/(n-1)
    
    do t = 0, n-1
       suma = 0.0_dp
       do k = 1, n-t
          suma = suma + (x(k)-avr)*(x(k+t)-avr)
       end do
       acf(t) = 1/((n-t)*var)*suma
    end do
    
  end function acf
  
  subroutine read_input()
    integer(i4) :: in_unit
    character(20) :: parameters_file

    read(*,*) parameters_file
    open(newunit = in_unit, file = parameters_file)
    read(in_unit, nml = input_parameters)
    close(in_unit)
    write(*,nml = input_parameters)
    
  end subroutine read_input

  
end program analysis

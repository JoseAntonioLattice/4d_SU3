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
  real(dp), allocatable, dimension(:,:) :: corr_poly
  real(dp), allocatable, dimension(:) :: avr_corr_poly, err_corr_poly
  real(dp), allocatable, dimension(:,:) :: correlation_polyakov_loop
  integer(i4) :: i,j, i_beta, bins1, bins2
  real(dp), allocatable, dimension(:) :: beta, auto_correlation

  beta = [5.7_dp]![(i*0.1_dp,i=1,80)]
  
  call read_input()

 

  allocate(Ep%array(N_measurements),auto_correlation(N_measurements), corr_poly(N_measurements,L/2-1))
  allocate(correlation_polyakov_loop(N_measurements,L/2-1))
  allocate(avr_corr_poly(L/2-1), err_corr_poly(L/2-1))
  data = "data/L="//trim(int2str(L))//"_equilibrium_"//trim(algorithm)//".dat"
  open(unit = 666, file = trim(data), status = "unknown")
  do i_beta = 1,size(beta)
     data_file = "data/L="//trim(int2str(L))//"/"//"equilibrium"//"/"//trim(algorithm)&
          //"/beta="//real2str(beta(i_beta))//"/observables.dat"
     !print*, trim(data_file)
     open( newunit = inunit, file = trim(data_file) )
     
     do i = 1, N_measurements
        read(inunit,*) Ep%array(i), correlation_polyakov_loop(i,:)
        corr_poly(i,:) = correlation_polyakov_loop(i,:)
     end do
     close(inunit)
     
     open(unit = 69, file = "data/L="//trim(int2str(L))//"/"//"equilibrium"//"/"//trim(algorithm)&
          //"/beta="//real2str(beta(i_beta))//"/auto_correlation.dat")
     auto_correlation = acf(Ep%array,N_measurements)
     do i = 1, size(Ep%array)
        write(69,*) auto_correlation(i)
     end do
     
     call max_jackknife_error_2(Ep%array,Ep%avr,Ep%err,bins1)
     do j = 1,L/2-1
        call max_jackknife_error_2(corr_poly(:,j),avr_corr_poly(j),err_corr_poly(j),bins2)
        avr_corr_poly(j) = avr_corr_poly(j)/(3*L**3)
        err_corr_poly(j) = err_corr_poly(j)/(3*L**3)
     end do
     write(666,*) beta(i_beta),Ep%avr, Ep%err,bins1,&!abs(avr_corr_poly)!, err_corr_poly, &
          -log(abs(avr_corr_poly))/L!, abs(err_corr_poly/(avr_corr_poly*L)) 
  end do
contains

  function acf(x,n)
    integer(i4), intent(in) :: n
    real(dp), intent(in), dimension(0:n-1) :: x
    real(dp), dimension(0:size(x)-1) :: acf
    integer(i4) :: k, t
    real(dp) :: avr, xim,d, suma
    


    avr = sum(x)/n
    do t = 0, n-1
       suma = 0.0_dp
       d = 0.0_dp
       do k = 0, n-1
          xim = x(k) - avr
          suma = suma + xim*(x(mod(k+t,n)) - avr)
          d = d + xim*xim
       end do
       acf(t) = suma/d
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

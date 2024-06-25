program analysis

  use iso_fortran_env, only : dp => real64, i4 => int32
  use statistics
  use number2string_mod
  implicit none

  integer(i4) :: Lx, Lt
  integer(i4) :: N_thermalization
  integer(i4) :: N_measurements
  integer(i4) :: N_skip
  character(25) :: algorithm
  logical :: equilibrium
  
  namelist /input_parameters/ Lx,Lt, N_thermalization, N_measurements, N_skip, algorithm, equilibrium
  namelist /analysis_parameters/ Lx,Lt, algorithm, equilibrium
  
  integer(i4) :: inunit
  character(256) :: data_file, data 
  
  type observable
     real(dp), allocatable :: array(:)
     real(dp) :: avr
     real(dp) :: err
  end type observable

  type(observable) :: Ep,abs_poly
  real(dp), allocatable, dimension(:,:) :: corr_poly
  complex(dp), allocatable, dimension(:) :: avr_corr_poly, err_corr_poly
  complex(dp), allocatable, dimension(:,:) :: correlation_polyakov_loop
  integer(i4) :: i,j,k,l, i_beta, bins1, bins2
  real(dp), allocatable, dimension(:) :: beta, auto_correlation
  complex(dp), allocatable :: poly(:)
  complex(dp) :: poly_avr, poly_err
  integer(i4), allocatable :: n_ms(:)
  logical :: ex

  beta = [real(dp) ::]
  beta = [beta, 5.7_dp]![(i*0.1_dp,i=1,80)]
  call read_input()  
  n_ms = [integer ::]
  i = 30
  do
     i = i + 1
     data_file = "data/Lx="//trim(int2str(Lx))//"/Lt="//trim(int2str(Lt))//"/"//"equilibrium"//"/"//trim(algorithm)&
          //"/beta="//real2str(beta(1))//"/observables_"//trim(int2str(i))//".dat"
     
     inquire(file = trim(data_file), exist = ex)
     if(ex .eqv. .false.)exit
     print*, data_file
     open(unit = 420, file = trim(data_file))
     read(420, nml = input_parameters)
     write(*, nml = input_parameters)
     n_ms = [n_ms, N_measurements]
  end do
  print*, n_ms,'holas'
  !stop
 
  allocate(Ep%array(sum(n_ms)),auto_correlation(sum(n_ms)), corr_poly(sum(n_ms),Lx/2-1))
  allocate(abs_poly%array(sum(n_ms)), poly(sum(n_ms)) )
  allocate(correlation_polyakov_loop(sum(n_ms),Lx/2-1))
  allocate(avr_corr_poly(Lx/2-1), err_corr_poly(Lx/2-1))
  data = "data/Lx="//trim(int2str(Lx))//"_Lt="//trim(int2str(Lt))//"_equilibrium_"//trim(algorithm)//".dat"
  open(unit = 666, file = trim(data), status = "unknown")
  open(unit = 777, file = 'data_Lx='//trim(int2str(Lx))//"_Lt="//trim(int2str(Lt))//'_correlation_polyakovloop.dat' &
       , status = "unknown")
  !print*, '1 ok'
  do i_beta = 1,size(beta)
     k = 30
     l = 1
     do
     k = k + 1   
     data_file = "data/Lx="//trim(int2str(Lx))//"/Lt="//trim(int2str(Lt))//"/"//"equilibrium"//"/"//trim(algorithm)&
          //"/beta="//real2str(beta(i_beta))//"/observables_"//trim(int2str(k))//".dat"
     
     inquire(file = trim(data_file), exist = ex)
     if(ex .eqv. .false.)exit
     open( newunit = inunit, file = trim(data_file) )
     print*, trim(data_file)
     read(inunit, nml = input_parameters)
     
     do i = (l-1)*n_ms(l)+1, sum(n_ms(1:l))
     
        read(inunit,*) Ep%array(i),poly(i), correlation_polyakov_loop(i,:)
       
        abs_poly%array(i) = abs(poly(i))
     end do
     close(inunit)
     l = l + 1
     end do
     
     open(unit = 69, file = "data/Lx="//trim(int2str(Lx))//"/Lt="//trim(int2str(Lt))//"/"//"equilibrium"//"/"//trim(algorithm)&
          //"/beta="//real2str(beta(i_beta))//"/auto_correlation.dat")
     auto_correlation = acf(Ep%array,sum(N_ms))
     
     do i = 1, size(Ep%array)
        write(69,*) auto_correlation(i)
     end do
     
     print*,'size E array', size(Ep%array)
     call max_jackknife_error_2(poly%re,poly_avr%re,poly_err%re,bins1)
     call max_jackknife_error_2(poly%im,poly_avr%im,poly_err%im,bins1)
     call max_jackknife_error_2(Ep%array,Ep%avr,Ep%err,bins1)
     call max_jackknife_error_2(abs_poly%array,abs_poly%avr,abs_poly%err,bins1)
     !abs_poly%avr = abs_poly%avr/Lx**3
     !abs_poly%err = abs_poly%err/Lx**3
     !poly_avr = poly_avr/Lx**3
     !poly_err = poly_err/Lx**3
     do j = 1,Lx/2-1
        call max_jackknife_error_2(correlation_polyakov_loop(:,j)%re,avr_corr_poly(j)%re,err_corr_poly(j)%re,bins2)
        call max_jackknife_error_2(correlation_polyakov_loop(:,j)%im,avr_corr_poly(j)%im,err_corr_poly(j)%im,bins2)
        !avr_corr_poly(j) = avr_corr_poly(j)/(3*Lx**3)
        !err_corr_poly(j) = err_corr_poly(j)/(3*Lx**3)
        !avr_corr_poly(j) = avr_corr_poly(j) - abs_poly%avr**2
        !err_corr_poly(j) = sqrt((err_corr_poly(j))**2 + 4*(abs_poly%avr*abs_poly%err)**2)
     end do
     write(666,*) beta(i_beta),Ep%avr, Ep%err,bins1
     do j = 1, Lx/2 -1
        write(777,*) -abs(log(avr_corr_poly(j)) )/Lt, abs(err_corr_poly(j)/(avr_corr_poly(j)*Lt))
     end do
  end do
  close(666)
  close(777)
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
    character(100) :: analysis_file

    read(*,*) analysis_file
    open(newunit = in_unit, file = trim(analysis_file))
    read(in_unit, nml = analysis_parameters)
    close(in_unit)
    write(*,nml = analysis_parameters)
    
  end subroutine read_input


  function complex_avr(x)
    complex(dp), intent(in) :: x(:)
    complex(dp) :: complex_avr

    complex_avr = sum(x)/size(x)

    
  end function complex_avr
end program analysis

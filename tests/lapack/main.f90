program main

  use iso_fortran_env, only : dp => real64
  implicit none

  integer, parameter :: n = 2
  integer, parameter :: lda = 2
  integer, parameter :: ldvl  = n
  integer, parameter :: ldvr = n
  integer, parameter :: lwork = 2*n
  
  complex(dp), dimension(n,n) :: A
  complex(dp), dimension(lwork) :: work
  complex(dp), dimension(n) :: w
  complex(dp), dimension(ldvl,n) :: vl
  complex(dp), dimension(ldvr,n) :: vr
  real(dp), dimension(2*n) :: rwork
  integer :: info
  
  
  !A = reshape([4,0,5,3,-5,4,-3,0,0,-3,4,5,3,-5,0,4],shape(A))

  A = reshape([(0.0_dp,1.0_dp),(1.0_dp,0.0_dp),&
               (2.5_dp,0.0_dp),(0.0_dp,0.0_dp)],shape(A))
  !print*, A
  call zgeev('N', 'V', n, A, lda, W, vl, ldvl, vr, ldvr, WORK, lwork, rwork,INFO)

  print*, info
  print*, w(1)
  print*, w(2)
  print*, vr(:,1)
  print*, vr(:,2) 
  
  
end program main

program main

  use iso_fortran_env, only : dp => real64
  implicit none

  integer, parameter :: n = 2
  integer, parameter :: lda = 2
  integer, parameter :: ldvl  = n
  integer, parameter :: ldvr = n
  integer, parameter :: lwork = 4*n
  
  real(dp), dimension(lda,n) :: A, B
  real(dp), dimension(lwork) :: work
  real(dp), dimension(n) :: wr,wi
  real(dp), dimension(ldvl,n) :: vl
  real(dp), dimension(ldvr,n) :: vr
  real(dp), dimension(4*n) :: rwork
  integer, dimension(n) :: ipiv
  integer :: info, info2
  integer, parameter :: lwork2 = n
  real(dp), dimension(lwork2) :: work2
  
  complex(dp), dimension(n,n) :: C
  complex(dp), dimension(n) :: eigenval
  complex(dp), dimension(n,n) :: reigenvec,leigenvec
  complex(dp), dimension(lwork) :: cwork
  real(dp), dimension(2*n) :: crwork

  
  !A = reshape([4,0,5,3,-5,4,-3,0,0,-3,4,5,3,-5,0,4],shape(A))

  A = reshape([3,1,5,-1],shape(A))
  !print*, A
  B = reshape([5.0_dp,1.0_dp,-1.0_dp,1.0_dp],shape(B))
  call dgeev('N', 'V', n, A, lda, Wr,wi, vl, ldvl, vr, ldvr, WORK, lwork,INFO)

  C = reshape([2.0_dp,2.0_dp,-5.0_dp,-4.0_dp],shape(B))
  call zgeev('N', 'V', n, C, lda, eigenval, leigenvec, ldvl, reigenvec, ldvr, cWORK, lwork,crwork,INFO)

  print*, eigenval
  leigenvec = 0.0_dp
  leigenvec(1,1) = exp(eigenval(1))
  leigenvec(2,2) = exp(eigenval(2))
  !reshape([1/6,-1,5,-1],shape(A))
  !call dgetri(n,B,lda,ipiv,work2,lwork2,info2)
  !vl(1,1) = exp(wr(1))
  !vl(2,2) = exp(wr(2))
  print*, matmul(reigenvec,matmul(leigenvec,cinv(reigenvec)))
  !B = reshape([w(1),(0.0_dp,0.0_dp),(0.0_dp,0.0_dp),w(2)],shape(B))
  !print*, info
  !print*, matmul(vr,matmul(vl,inv(vr)))
  !print*, info2
  !print*,B
  !print*, matmul(vl,vr)
  !print*, wr
  !print*, A
  !print*, w(1)
  !print*, w(2)
  !print*, vl(:,1)
  !print*, vl(:,2) 
  !print*, matmul(B,vr(:,1)) 
  !print*, w(1)*vr(:,1)

contains
  
  ! -- Returns the inverse of a general squared matrix A
  function inv(A) result(Ainv)
    implicit none
    real(dp),intent(in) :: A(:,:)
    real(dp)            :: Ainv(size(A,1),size(A,2))
    real(dp)            :: work(size(A,1))            ! work array for LAPACK
    integer         :: n,info,ipiv(size(A,1))     ! pivot indices
    
    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)
    ! SGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call dGETRF(n,n,Ainv,n,ipiv,info)
    if (info.ne.0) stop 'Matrix is numerically singular!'
    ! SGETRI computes the inverse of a matrix using the LU factorization
    ! computed by SGETRF.
    call dGETRI(n,Ainv,n,ipiv,work,n,info)
    if (info.ne.0) stop 'Matrix inversion failed!'
  end function inv

    ! -- Returns the inverse of a general squared matrix A
  function cinv(A) result(Ainv)
    implicit none
    complex(dp),intent(in) :: A(:,:)
    complex(dp)            :: Ainv(size(A,1),size(A,2))
    complex(dp)            :: work(size(A,1))            ! work array for LAPACK
    integer         :: n,info,ipiv(size(A,1))     ! pivot indices
    
    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)
    ! SGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call zGETRF(n,n,Ainv,n,ipiv,info)
    if (info.ne.0) stop 'Matrix is numerically singular!'
    ! SGETRI computes the inverse of a matrix using the LU factorization
    ! computed by SGETRF.
    call zGETRI(n,Ainv,n,ipiv,work,n,info)
    if (info.ne.0) stop 'Matrix inversion failed!'
  end function cinv
  
end program main

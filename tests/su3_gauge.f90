!=============================================================================
! SU(3) pure gauge theory in 4 dimensions
!
! Contents
! --------
!  1. Parameters & global types
!  2. SU(3) algebra: matmul, dagger, tr, det, TA, exp
!  3. Lattice geometry: periodic boundary conditions
!  4. Observables: plaquette, action, staples
!  5. Heatbath update  (Kennedy-Pendleton algorithm for SU(2) subgroups)
!  6. Gradient flow    (Luscher RK4 integrator)
!  7. Topological charge  (clover estimator, reduced 3-term formula)
!  8. Main program
!=============================================================================
module su3_types
  use iso_fortran_env, only: dp => real64, i4 => int32
  implicit none

  ! ---- lattice parameters (change here) ----
  integer(i4), parameter :: Lx = 16        ! spatial extent
  integer(i4), parameter :: Lt = 16        ! temporal extent
  integer(i4), parameter :: Nd = 4        ! number of dimensions
  real(dp),    parameter :: beta = 6.0_dp ! Wilson beta
  integer(i4), parameter :: N_therm = 200 ! thermalization sweeps
  integer(i4), parameter :: N_meas  = 20  ! configurations to measure
  integer(i4), parameter :: N_skip  = 10  ! sweeps between measurements
  integer(i4), parameter :: N_flow  = 150 ! gradient flow steps per config
  real(dp),    parameter :: eps_flow = 0.05_dp ! flow step size
  integer(i4), parameter :: N_print = 10  ! print Q every N_print flow steps
  ! ------------------------------------------

  real(dp), parameter :: pi = acos(-1.0_dp)
  integer(i4), parameter :: Nc = 3         ! number of colors

  ! 3x3 complex matrix
  type :: mat3
    complex(dp) :: m(3,3)
  end type mat3

  ! The gauge field: U(x1,x2,x3,x4,mu)
  type(mat3), allocatable :: U(:,:,:,:,:)

  ! Levi-Civita symbol
  integer(i4) :: eps4(4,4,4,4)

end module su3_types

!=============================================================================
! SU(3) matrix operations
!=============================================================================
module su3_algebra
  use su3_types
  implicit none

contains

  ! ---------- matrix multiply ----------
  function mm(A,B) result(C)
    type(mat3), intent(in) :: A, B
    type(mat3) :: C
    integer :: i,j,k
    C%m = 0.0_dp
    do i=1,3; do j=1,3; do k=1,3
      C%m(i,j) = C%m(i,j) + A%m(i,k)*B%m(k,j)
    end do; end do; end do
  end function mm

  ! ---------- hermitian conjugate ----------
  function dag(A) result(B)
    type(mat3), intent(in) :: A
    type(mat3) :: B
    B%m = conjg(transpose(A%m))
  end function dag

  ! ---------- trace ----------
  function tr3(A) result(t)
    type(mat3), intent(in) :: A
    complex(dp) :: t
    t = A%m(1,1) + A%m(2,2) + A%m(3,3)
  end function tr3

  ! ---------- determinant (3x3) ----------
  function det3(A) result(d)
    type(mat3), intent(in) :: A
    complex(dp) :: d
    associate(m => A%m)
      d = m(1,1)*(m(2,2)*m(3,3)-m(2,3)*m(3,2)) &
        - m(1,2)*(m(2,1)*m(3,3)-m(2,3)*m(3,1)) &
        + m(1,3)*(m(2,1)*m(3,2)-m(2,2)*m(3,1))
    end associate
  end function det3

  ! ---------- identity matrix ----------
  function eye3() result(I)
    type(mat3) :: I
    I%m = 0.0_dp
    I%m(1,1) = 1.0_dp; I%m(2,2) = 1.0_dp; I%m(3,3) = 1.0_dp
  end function eye3

  ! ---------- traceless anti-Hermitian projection T_A ----------
  ! T_A(W) = (W - W†)/2 - tr(W - W†)/6 * I
  function TA(W) result(X)
    type(mat3), intent(in) :: W
    type(mat3) :: X, AH, Id
    complex(dp) :: tr_AH
    AH%m = (W%m - conjg(transpose(W%m))) * 0.5_dp
    tr_AH = tr3(AH)
    Id = eye3()
    X%m = AH%m - (tr_AH/3.0_dp) * Id%m
  end function TA

  ! ---------- matrix exponential via Cayley-Hamilton (Luscher 2009) ----------
  ! For X in su(3), exp(X) = q0*I + q1*X + q2*X^2
  ! Horner recursion with K=20 terms (sufficient for ||X||<1).
  function mat_exp(X) result(E)
    type(mat3), intent(in) :: X
    type(mat3) :: E, X2, Id
    integer, parameter :: K = 20
    complex(dp) :: q0, q1, q2, q0o, q1o, q2o
    complex(dp) :: t, d, trX2
    integer :: n

    X2%m = matmul(X%m, X%m)
    trX2 = tr3(X2)
    t = -0.5_dp * trX2           ! Cayley-Hamilton invariant t
    d = (0.0_dp,1.0_dp)*det3(X)  ! Cayley-Hamilton invariant d = i*det(X)

    q0o = 1.0_dp / gamma(real(K+1,dp))
    q1o = (0.0_dp, 0.0_dp)
    q2o = (0.0_dp, 0.0_dp)

    do n = K-1, 0, -1
      q0 = 1.0_dp/gamma(real(n+1,dp)) - (0.0_dp,1.0_dp)*d*q2o
      q1 = q0o - t*q2o
      q2 = q1o
      q0o = q0; q1o = q1; q2o = q2
    end do

    Id = eye3()
    E%m = q0*Id%m + q1*X%m + q2*X2%m
  end function mat_exp

  ! ---------- reunitarize by Gram-Schmidt + det correction ----------
  subroutine reunitarize(A)
    type(mat3), intent(inout) :: A
    complex(dp) :: r1(3), r2(3), r3(3)
    complex(dp) :: proj, nrm, d
    integer :: i

    ! Row 1
    r1 = A%m(1,:)
    nrm = sqrt(real(dot_product(r1,r1),dp))
    r1 = r1 / nrm

    ! Row 2: subtract projection onto row 1
    r2 = A%m(2,:)
    proj = dot_product(r1, r2)          ! <r1|r2>
    r2 = r2 - proj*r1
    nrm = sqrt(real(dot_product(r2,r2),dp))
    r2 = r2 / nrm

    ! Row 3: cross product of conj(r1) x conj(r2)  -> ensures det = +1
    r3(1) = conjg(r1(2))*conjg(r2(3)) - conjg(r1(3))*conjg(r2(2))
    r3(2) = conjg(r1(3))*conjg(r2(1)) - conjg(r1(1))*conjg(r2(3))
    r3(3) = conjg(r1(1))*conjg(r2(2)) - conjg(r1(2))*conjg(r2(1))

    A%m(1,:) = r1
    A%m(2,:) = r2
    A%m(3,:) = r3
  end subroutine reunitarize

end module su3_algebra

!=============================================================================
! Periodic boundary conditions
!=============================================================================
module lattice_bc
  use su3_types
  implicit none

contains

  ! Periodic index shift +1 in direction mu (mu=1..4 -> dims x,y,z,t)
  function ip(x, mu) result(xp)
    integer(i4), intent(in) :: x(4), mu
    integer(i4) :: xp(4), L(4)
    L = [Lx, Lx, Lx, Lt]
    xp = x
    xp(mu) = mod(x(mu), L(mu)) + 1
  end function ip

  ! Periodic index shift -1 in direction mu
  function im(x, mu) result(xm)
    integer(i4), intent(in) :: x(4), mu
    integer(i4) :: xm(4), L(4)
    L = [Lx, Lx, Lx, Lt]
    xm = x
    xm(mu) = mod(x(mu) - 2 + L(mu), L(mu)) + 1
  end function im

end module lattice_bc

!=============================================================================
! Observables: plaquette, action, staples, clover F_{mu nu}
!=============================================================================
module observables
  use su3_types
  use su3_algebra
  use lattice_bc
  implicit none

contains

  ! ---- single plaquette U_mu(x) U_nu(x+mu) U_mu†(x+nu) U_nu†(x) ----
  function plaquette(V, x, mu, nu) result(P)
    type(mat3), intent(in) :: V(:,:,:,:,:)
    integer(i4), intent(in) :: x(4), mu, nu
    type(mat3) :: P
    integer(i4) :: xpmu(4), xpnu(4)
    xpmu = ip(x,mu); xpnu = ip(x,nu)
    P = mm(mm(mm(V(x(1),x(2),x(3),x(4),mu), &
                  V(xpmu(1),xpmu(2),xpmu(3),xpmu(4),nu)), &
               dag(V(xpnu(1),xpnu(2),xpnu(3),xpnu(4),mu))), &
            dag(V(x(1),x(2),x(3),x(4),nu)))
  end function plaquette

  ! ---- Wilson action S = (beta/3) * sum_{x,mu<nu} Re tr(1 - P_{mu nu}) ----
  ! Returns plaquette average P = <Re tr P_{mu nu}> / 3
  function wilson_action(V) result(Sp)
    type(mat3), intent(in) :: V(:,:,:,:,:)
    real(dp) :: Sp
    integer(i4) :: x(4), mu, nu, x1,x2,x3,x4
    integer(i4), parameter :: Npl = Nd*(Nd-1)/2
    real(dp) :: S

    S = 0.0_dp
    do x4=1,Lt; do x3=1,Lx; do x2=1,Lx; do x1=1,Lx
      x = [x1,x2,x3,x4]
      do mu=1,Nd-1
        do nu=mu+1,Nd
          S = S + real(tr3(plaquette(V,x,mu,nu)), dp)
        end do
      end do
    end do; end do; end do; end do
    Sp = S / (3.0_dp * Npl * Lx**3 * Lt)
  end function wilson_action

  ! ---- sum of staples for link U_mu(x) ----
  ! Sigma_mu(x) = sum_{nu /= mu} [ U_nu(x) U_mu(x+nu) U_nu†(x+mu)
  !                               + U_nu†(x-nu) U_mu(x-nu) U_nu(x+mu-nu) ]
  function staples(V, x, mu) result(Sig)
    type(mat3), intent(in) :: V(:,:,:,:,:)
    integer(i4), intent(in) :: x(4), mu
    type(mat3) :: Sig
    type(mat3) :: tmp
    integer(i4) :: nu, xpmu(4), xpnu(4), xmnu(4), xpmumnu(4)

    Sig%m = 0.0_dp
    xpmu = ip(x, mu)
    do nu = 1, Nd
      if (nu == mu) cycle

      xpnu    = ip(x, nu)
      xmnu    = im(x, nu)
      xpmumnu = im(xpmu, nu)

      ! Forward staple: U_nu(x) * U_mu(x+nu) * U_nu†(x+mu)
      tmp = mm(mm(V(x(1),x(2),x(3),x(4),nu), &
                  V(xpnu(1),xpnu(2),xpnu(3),xpnu(4),mu)), &
               dag(V(xpmu(1),xpmu(2),xpmu(3),xpmu(4),nu)))
      Sig%m = Sig%m + tmp%m

      ! Backward staple: U_nu†(x-nu) * U_mu(x-nu) * U_nu(x+mu-nu)
      tmp = mm(mm(dag(V(xmnu(1),xmnu(2),xmnu(3),xmnu(4),nu)), &
                   V(xmnu(1),xmnu(2),xmnu(3),xmnu(4),mu)), &
                V(xpmumnu(1),xpmumnu(2),xpmumnu(3),xpmumnu(4),nu))
      Sig%m = Sig%m + tmp%m
    end do
  end function staples

  ! ---- clover Q_{mu nu}(x): sum of 4 oriented plaquettes (raw, no 1/8i) ----
  function clover_Q(V, x, mu, nu) result(Q)
    type(mat3), intent(in) :: V(:,:,:,:,:)
    integer(i4), intent(in) :: x(4), mu, nu
    type(mat3) :: Q, tmp
    integer(i4) :: xpmu(4), xpnu(4), xmmu(4), xmnu(4)
    integer(i4) :: xmmupnu(4), xmmuMnu(4), xpmuMnu(4)

    xpmu    = ip(x,mu);    xpnu    = ip(x,nu)
    xmmu    = im(x,mu);    xmnu    = im(x,nu)
    xmmupnu = im(xpnu,mu); xmmuMnu = im(xmnu,mu); xpmuMnu = ip(xmnu,mu)

    ! P_{++}
    Q = mm(mm(mm(V(x(1),x(2),x(3),x(4),mu), &
                  V(xpmu(1),xpmu(2),xpmu(3),xpmu(4),nu)), &
               dag(V(xpnu(1),xpnu(2),xpnu(3),xpnu(4),mu))), &
            dag(V(x(1),x(2),x(3),x(4),nu)))
    ! P_{-+}
    tmp = mm(mm(mm(V(x(1),x(2),x(3),x(4),nu), &
                    dag(V(xmmupnu(1),xmmupnu(2),xmmupnu(3),xmmupnu(4),mu))), &
                 dag(V(xmmu(1),xmmu(2),xmmu(3),xmmu(4),nu))), &
              V(xmmu(1),xmmu(2),xmmu(3),xmmu(4),mu))
    Q%m = Q%m + tmp%m
    ! P_{--}
    tmp = mm(mm(mm(dag(V(xmmu(1),xmmu(2),xmmu(3),xmmu(4),mu)), &
                    dag(V(xmmuMnu(1),xmmuMnu(2),xmmuMnu(3),xmmuMnu(4),nu))), &
                 V(xmmuMnu(1),xmmuMnu(2),xmmuMnu(3),xmmuMnu(4),mu)), &
              V(xmnu(1),xmnu(2),xmnu(3),xmnu(4),nu))
    Q%m = Q%m + tmp%m
    ! P_{+-}
    tmp = mm(mm(mm(dag(V(xmnu(1),xmnu(2),xmnu(3),xmnu(4),nu)), &
                    V(xmnu(1),xmnu(2),xmnu(3),xmnu(4),mu)), &
                 V(xpmuMnu(1),xpmuMnu(2),xpmuMnu(3),xpmuMnu(4),nu)), &
              dag(V(x(1),x(2),x(3),x(4),mu)))
    Q%m = Q%m + tmp%m
  end function clover_Q

  ! ---- topological charge density at site x (3-term formula, Garcia Hernandez) ----
  ! q(x) = - (1/128 pi^2) * sum_{[mnrs] in {[1234],[1324],[1423]}}
  !             eps_{mnrs} tr( Q_{mn} [Q_{rs} - Q_{rs}†] )
  function topo_density(V, x) result(qx)
    type(mat3), intent(in) :: V(:,:,:,:,:)
    integer(i4), intent(in) :: x(4)
    real(dp) :: qx
    type(mat3) :: Q12, Q34, Q13, Q24, Q14, Q23
    complex(dp) :: s

    Q12 = clover_Q(V,x,1,2);  Q34 = clover_Q(V,x,3,4)
    Q13 = clover_Q(V,x,1,3);  Q24 = clover_Q(V,x,2,4)
    Q14 = clover_Q(V,x,1,4);  Q23 = clover_Q(V,x,2,3)

    ! eps_{1234}=+1: tr(Q12*(Q34-Q34†))
    s  =  tr3(mm(Q12, mat3_sub_dag(Q34)))
    ! eps_{1324}=-1: -tr(Q13*(Q24-Q24†))
    s  = s - tr3(mm(Q13, mat3_sub_dag(Q24)))
    ! eps_{1423}=+1: tr(Q14*(Q23-Q23†))
    s  = s + tr3(mm(Q14, mat3_sub_dag(Q23)))

    qx = real(-s / (128.0_dp * pi**2), dp)
  end function topo_density

  ! helper: A - A†
  function mat3_sub_dag(A) result(B)
    type(mat3), intent(in) :: A
    type(mat3) :: B
    B%m = A%m - conjg(transpose(A%m))
  end function mat3_sub_dag

  ! ---- total topological charge Q = sum_x q(x) ----
  function topo_charge(V) result(Q)
    type(mat3), intent(in) :: V(:,:,:,:,:)
    real(dp) :: Q
    integer(i4) :: x1,x2,x3,x4
    Q = 0.0_dp
    do x4=1,Lt; do x3=1,Lx; do x2=1,Lx; do x1=1,Lx
      Q = Q + topo_density(V,[x1,x2,x3,x4])
    end do; end do; end do; end do
  end function topo_charge

end module observables

!=============================================================================
! Heatbath update (Kennedy-Pendleton algorithm for SU(2) subgroups)
!=============================================================================
module heatbath_mod
  use su3_types
  use su3_algebra
  use lattice_bc
  use observables
  implicit none

contains

  ! ---- random number in [0,1) ----
  function rand01() result(r)
    real(dp) :: r
    call random_number(r)
  end function rand01

  ! ---- draw SU(2) matrix from heat-bath distribution ----
  ! Generates X ~ exp(alpha * beta_N * x0) following Kennedy-Pendleton.
  subroutine heatbath_su2(alpha, beta_N, x0, x1, x2, x3)
    real(dp), intent(in)  :: alpha, beta_N
    real(dp), intent(out) :: x0, x1, x2, x3
    real(dp) :: r1, r2, r3, lambda2, tmp, phi, ct, st
    real(dp) :: a   ! = alpha * beta_N / Nc
    logical :: accepted

    a = alpha * beta_N / real(Nc, dp)

    ! Kennedy-Pendleton: generate x0 with distribution ~ sqrt(1-x0^2) exp(a*x0)
    accepted = .false.
    do while (.not. accepted)
      r1 = rand01(); r2 = rand01(); r3 = rand01()
      ! Avoid log(0)
      if (r1 < 1e-15_dp) r1 = 1e-15_dp
      if (r2 < 1e-15_dp) r2 = 1e-15_dp
      if (r3 < 1e-15_dp) r3 = 1e-15_dp
      lambda2 = -(1.0_dp/(2.0_dp*a)) * (log(r1) + cos(2.0_dp*pi*r2)**2 * log(r3))
      tmp = rand01()
      if (tmp**2 <= 1.0_dp - lambda2) then
        accepted = .true.
      end if
    end do
    x0 = 1.0_dp - 2.0_dp*lambda2

    ! Random unit vector (x1,x2,x3) on sphere of radius sqrt(1-x0^2)
    tmp = sqrt(max(1.0_dp - x0**2, 0.0_dp))
    ct  = 2.0_dp*rand01() - 1.0_dp          ! cos(theta)
    phi = 2.0_dp*pi*rand01()
    st  = sqrt(max(1.0_dp - ct**2, 0.0_dp))
    x1 = tmp * st * cos(phi)
    x2 = tmp * st * sin(phi)
    x3 = tmp * ct
  end subroutine heatbath_su2

  ! ---- embed SU(2) update into SU(3) for a given subgroup ij ----
  ! sub=1 -> (1,2), sub=2 -> (1,3), sub=3 -> (2,3)
  subroutine su2_update(W, Sig, sub, beta_N)
    type(mat3), intent(inout) :: W
    type(mat3), intent(in)    :: Sig
    integer,    intent(in)    :: sub
    real(dp),   intent(in)    :: beta_N

    integer  :: i, j
    type(mat3) :: V, Xmat, Xinv
    complex(dp) :: a11, a12, a21, a22
    real(dp)    :: alpha, x0, x1, x2, x3
    complex(dp) :: u0, u1, u2, u3

    ! Row/column indices for this SU(2) subgroup
    select case(sub)
      case(1); i=1; j=2
      case(2); i=1; j=3
      case(3); i=2; j=3
    end select

    ! V = W * Sig  -- project onto the (i,j) SU(2) block
    V = mm(W, Sig)
    a11 = V%m(i,i); a12 = V%m(i,j)
    a21 = V%m(j,i); a22 = V%m(j,j)

    ! alpha = sqrt(det of the 2x2 block)
    alpha = sqrt(max(real(a11*a22 - a12*a21, dp), 0.0_dp))
    if (alpha < 1e-14_dp) return  ! degenerate, skip

    ! Draw new SU(2) element X = x0*I + i*(x1*s1 + x2*s2 + x3*s3)
    call heatbath_su2(alpha, beta_N, x0, x1, x2, x3)

    ! Build the SU(3) matrix that replaces the (i,j) block
    ! X_new operates on rows/columns i,j
    Xmat = eye3()
    u0 = cmplx(x0,  0.0_dp, dp)
    u1 = cmplx(0.0_dp,  x3, dp)
    u2 = cmplx(x2,   x1, dp)   ! x2 + i*x1
    u3 = cmplx(-x2,  x1, dp)   ! -x2 + i*x1

    ! The new link update: replace W by X * V^{-1}/alpha * W
    ! Equivalently, the SU(2) matrix acting on the (i,j) subspace is:
    !   new_block = (1/alpha) * [x0 -x3    ] * [a11† a21†]
    !                           [x3  x0    ]   [a12† a22†]
    ! We construct the full SU(3) update matrix directly.
    Xmat%m = 0.0_dp
    Xmat%m(i,i) =  cmplx(x0,  x3, dp)
    Xmat%m(i,j) =  cmplx(x2,  x1, dp)
    Xmat%m(j,i) = -cmplx(x2, -x1, dp)
    Xmat%m(j,j) =  cmplx(x0, -x3, dp)
    ! Set the remaining diagonal element (the one not in {i,j}) to 1
    block
      integer :: k
      k = 6 - i - j   ! 6-1-2=3, 6-1-3=2, 6-2-3=1
      Xmat%m(k,k) = 1.0_dp
    end block

    ! The SU(2) submatrix must multiply by (a†/alpha) to form X_normalized
    ! W_new = X_normalized * W
    ! Build Xinv = (normalized block of Sig) embedded in SU(3)
    Xinv = eye3()
    Xinv%m(i,i) = conjg(a11)/alpha
    Xinv%m(i,j) = conjg(a21)/alpha
    Xinv%m(j,i) = conjg(a12)/alpha
    Xinv%m(j,j) = conjg(a22)/alpha

    W = mm(mm(Xmat, Xinv), W)
    call reunitarize(W)
  end subroutine su2_update

  ! ---- full heatbath update for link U_mu(x) ----
  subroutine hb_link(x, mu)
    integer(i4), intent(in) :: x(4), mu
    type(mat3) :: Sig
    integer :: sub

    Sig = staples(U, x, mu)
    ! Cycle over the 3 SU(2) subgroups
    do sub = 1, 3
      call su2_update(U(x(1),x(2),x(3),x(4),mu), Sig, sub, beta)
    end do
  end subroutine hb_link

  ! ---- one heatbath sweep over the whole lattice ----
  subroutine hb_sweep()
    integer(i4) :: x1,x2,x3,x4,mu
    do x4=1,Lt; do x3=1,Lx; do x2=1,Lx; do x1=1,Lx
      do mu=1,Nd
        call hb_link([x1,x2,x3,x4], mu)
      end do
    end do; end do; end do; end do
  end subroutine hb_sweep

end module heatbath_mod

!=============================================================================
! Gradient flow  (Luscher RK4, 3rd order in epsilon)
!=============================================================================
module gradient_flow
  use su3_types
  use su3_algebra
  use lattice_bc
  use observables
  implicit none

contains

  ! ---- zeta: the Lie-algebra force Z_mu(x) = -T_A(U_mu * Sigma†) ----
  function zeta(V, x, mu) result(Z)
    type(mat3), intent(in) :: V(:,:,:,:,:)
    integer(i4), intent(in) :: x(4), mu
    type(mat3) :: Z, Sig
    Sig = staples(V, x, mu)
    Z = TA(mm(V(x(1),x(2),x(3),x(4),mu), dag(Sig)))
    Z%m = -Z%m
  end function zeta

  ! ---- one RK4 step of duration epsilon ----
  ! W0 = V(t)
  ! W1 = exp((1/4) Z0) W0               Z0 = eps*Z(W0)
  ! W2 = exp((8/9) Z1 - (17/36) Z0) W1  Z1 = eps*Z(W1)
  ! W3 = exp((3/4) Z2 - (8/9)Z1 + (17/36)Z0) W2  Z2=eps*Z(W2)
  ! V(t+eps) = W3
  subroutine flow_step(V, eps)
    type(mat3), intent(inout) :: V(:,:,:,:,:)
    real(dp),   intent(in)    :: eps
    type(mat3), allocatable :: W(:,:,:,:,:)
    type(mat3), allocatable :: Z0(:,:,:,:,:), Z1(:,:,:,:,:)
    type(mat3) :: Ztmp, A
    integer(i4) :: x1,x2,x3,x4,mu
    integer(i4) :: x(4)

    allocate(W(Lx,Lx,Lx,Lt,Nd))
    allocate(Z0(Lx,Lx,Lx,Lt,Nd))
    allocate(Z1(Lx,Lx,Lx,Lt,Nd))

    ! --- stage 1: compute Z0, advance W0 -> W1 ---
    W = V
    do x4=1,Lt; do x3=1,Lx; do x2=1,Lx; do x1=1,Lx
      x = [x1,x2,x3,x4]
      do mu=1,Nd
        Ztmp = zeta(W, x, mu)
        Z0(x1,x2,x3,x4,mu)%m = eps * Ztmp%m
      end do
    end do; end do; end do; end do
    do x4=1,Lt; do x3=1,Lx; do x2=1,Lx; do x1=1,Lx
      do mu=1,Nd
        A%m = 0.25_dp * Z0(x1,x2,x3,x4,mu)%m
        W(x1,x2,x3,x4,mu) = mm(mat_exp(A), W(x1,x2,x3,x4,mu))
      end do
    end do; end do; end do; end do

    ! --- stage 2: compute Z1, advance W1 -> W2 ---
    do x4=1,Lt; do x3=1,Lx; do x2=1,Lx; do x1=1,Lx
      x = [x1,x2,x3,x4]
      do mu=1,Nd
        Ztmp = zeta(W, x, mu)
        Z1(x1,x2,x3,x4,mu)%m = eps * Ztmp%m
      end do
    end do; end do; end do; end do
    do x4=1,Lt; do x3=1,Lx; do x2=1,Lx; do x1=1,Lx
      do mu=1,Nd
        A%m = (8.0_dp/9.0_dp)  * Z1(x1,x2,x3,x4,mu)%m &
            - (17.0_dp/36.0_dp)* Z0(x1,x2,x3,x4,mu)%m
        W(x1,x2,x3,x4,mu) = mm(mat_exp(A), W(x1,x2,x3,x4,mu))
      end do
    end do; end do; end do; end do

    ! --- stage 3: compute Z2, advance W2 -> W3 = V(t+eps) ---
    do x4=1,Lt; do x3=1,Lx; do x2=1,Lx; do x1=1,Lx
      x = [x1,x2,x3,x4]
      do mu=1,Nd
        Ztmp = zeta(W, x, mu)
        A%m = 0.75_dp * eps * Ztmp%m &
            - (8.0_dp/9.0_dp)  * Z1(x1,x2,x3,x4,mu)%m &
            + (17.0_dp/36.0_dp)* Z0(x1,x2,x3,x4,mu)%m
        V(x1,x2,x3,x4,mu) = mm(mat_exp(A), W(x1,x2,x3,x4,mu))
      end do
    end do; end do; end do; end do

    deallocate(W, Z0, Z1)
  end subroutine flow_step

  ! ---- run N_flow steps of gradient flow, print/save observables ----
  subroutine run_flow(V, iconf, funit)
    type(mat3), intent(inout) :: V(:,:,:,:,:)
    integer,    intent(in)    :: iconf   ! configuration index (for output)
    integer,    intent(in)    :: funit   ! open file unit to write to
    type(mat3), allocatable :: Vflow(:,:,:,:,:)
    real(dp) :: t, S, Q
    integer :: k

    allocate(Vflow(Lx,Lx,Lx,Lt,Nd))
    Vflow = V   ! work on a copy, don't modify the MC configuration

    do k = 0, N_flow
      t = k * eps_flow
      if (mod(k, N_print) == 0) then
        S = wilson_action(Vflow)
        Q = topo_charge(Vflow)
        write(funit, '(I4,2X,I5,2X,F10.4,2X,F14.8,2X,F14.6)') &
             iconf, k, t, S, Q
        write(*,'(A,I4,A,I4,A,F6.2,A,F10.6,A,F10.4)') &
             '  conf=', iconf, '  step=', k, '  t=', t, &
             '  <P>=', S, '  Q=', Q
      end if
      if (k < N_flow) call flow_step(Vflow, eps_flow)
    end do

    deallocate(Vflow)
  end subroutine run_flow

end module gradient_flow

!=============================================================================
! Main program
!=============================================================================
program su3_pure_gauge
  use su3_types
  use su3_algebra
  use lattice_bc
  use observables
  use heatbath_mod
  use gradient_flow
  implicit none

  integer(i4) :: x1,x2,x3,x4,mu, isw, iconf
  real(dp)    :: S
  integer     :: funit_flow = 10, funit_plaq = 11
  character(len=64) :: fname

  ! ---- initialize Levi-Civita ----
  call init_eps4()

  ! ---- allocate gauge field ----
  allocate(U(Lx,Lx,Lx,Lt,Nd))

  ! ---- hot start: all links = random SU(3) ----
  call random_seed()
  do x4=1,Lt; do x3=1,Lx; do x2=1,Lx; do x1=1,Lx
    do mu=1,Nd
      call random_su3(U(x1,x2,x3,x4,mu))
    end do
  end do; end do; end do; end do

  ! ---- open output files ----
  call system('mkdir -p data')
  open(funit_plaq, file='data/plaquette.dat', status='replace', action='write')
  write(funit_plaq,'(A)') '# sweep    <P>'
  open(funit_flow,  file='data/gradient_flow.dat', status='replace', action='write')
  write(funit_flow,'(A)') &
    '# conf   flow_step    t           <P>           Q'

  ! ---- thermalization ----
  write(*,*) '=== Thermalization ==='
  write(*,'(A,I0,A,F6.3,A,I0,A,I0)') &
    '  Lx=', Lx, '  beta=', beta, &
    '  N_therm=', N_therm, '  Nc=', Nc
  do isw = 1, N_therm
    call hb_sweep()
    if (mod(isw,20)==0) then
      S = wilson_action(U)
      write(*,'(A,I4,A,F10.6)') '  sweep ', isw, '   <P>=', S
      write(funit_plaq,'(I6,2X,F14.10)') isw, S
    end if
  end do

  ! ---- measurement loop ----
  write(*,*) ''
  write(*,*) '=== Measurement + Gradient Flow ==='
  do iconf = 1, N_meas
    ! skip sweeps between measurements
    do isw = 1, N_skip
      call hb_sweep()
    end do

    S = wilson_action(U)
    write(*,'(A,I3,A,F10.6)') 'Config ', iconf, '   <P>=', S
    write(funit_plaq,'(I6,2X,F14.10)') N_therm + (iconf-1)*N_skip + N_skip, S

    ! apply gradient flow and measure Q
    call run_flow(U, iconf, funit_flow)
    write(*,*)
  end do

  close(funit_plaq)
  close(funit_flow)
  write(*,*) 'Done. Output in data/plaquette.dat and data/gradient_flow.dat'
  deallocate(U)

contains

  ! ---- generate a random SU(3) matrix via random complex matrix + GS ----
  subroutine random_su3(A)
    type(mat3), intent(out) :: A
    real(dp) :: re(3,3), im(3,3)
    call random_number(re); call random_number(im)
    A%m = cmplx(re - 0.5_dp, im - 0.5_dp, dp)
    call reunitarize(A)
  end subroutine random_su3

  subroutine init_eps4()
    integer :: a,b,c,d, perm(4), i, sgn
    integer :: perms(24,4), signs(24)
    integer :: idx

    eps4 = 0
    ! Generate all 24 permutations of (1,2,3,4) with their signs
    idx = 0
    call gen_perms([1,2,3,4], 4, perms, signs, idx)
    do i = 1, 24
      eps4(perms(i,1),perms(i,2),perms(i,3),perms(i,4)) = signs(i)
    end do
  end subroutine init_eps4

  recursive subroutine gen_perms(arr, n, perms, signs, idx)
    integer, intent(in)    :: arr(4), n
    integer, intent(inout) :: perms(24,4), signs(24), idx
    integer :: arr2(4), i, tmp
    if (n == 1) then
      idx = idx + 1
      perms(idx,:) = arr
      signs(idx) = perm_sign(arr)
      return
    end if
    do i = 1, n
      arr2 = arr
      tmp = arr2(i); arr2(i) = arr2(n); arr2(n) = tmp
      call gen_perms(arr2, n-1, perms, signs, idx)
    end do
  end subroutine gen_perms

  function perm_sign(p) result(s)
    integer, intent(in) :: p(4)
    integer :: s, i, j, inv
    inv = 0
    do i = 1, 4
      do j = i+1, 4
        if (p(i) > p(j)) inv = inv + 1
      end do
    end do
    s = 1 - 2*mod(inv,2)
  end function perm_sign

end program su3_pure_gauge

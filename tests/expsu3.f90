module su3_exponential
    use, intrinsic :: iso_fortran_env, only: wp => real64
    implicit none
    
    complex(wp), parameter :: ii = (0.0_wp, 1.0_wp)

contains

    subroutine exp_su3(X, expX)
        ! X es una matriz de 3x3 en su(3): antihermítica y de traza cero
        complex(wp), intent(in)  :: X(3,3)
        complex(wp), intent(out) :: expX(3,3)
        
        complex(wp) :: X2(3,3), Identity(3,3)
        real(wp)    :: detX_im, c0, c1
        complex(wp) :: detX
        real(wp)    :: q, r, theta
        complex(wp) :: la(3) ! Autovalores
        complex(wp) :: f0, f1, f2
        integer     :: i
        
        ! 1. Definir la matriz Identidad
        Identity = reshape([ (1.0_wp, 0.0_wp), (0.0_wp, 0.0_wp), (0.0_wp, 0.0_wp), &
                             (0.0_wp, 0.0_wp), (1.0_wp, 0.0_wp), (0.0_wp, 0.0_wp), &
                             (0.0_wp, 0.0_wp), (0.0_wp, 0.0_wp), (1.0_wp, 0.0_wp) ], [3,3])
        
        ! X^2 es necesaria para la combinación lineal
        X2 = matmul(X, X)
        
        ! 2. Invariantes algebraicos para una matriz antihermítica de traza cero
        ! c1 = -1/2 * Tr(X^2)  --> Como X es antihermítica, c1 siempre es real y >= 0
        c1 = -0.5_wp * real(X2(1,1) + X2(2,2) + X2(3,3))
        
        ! c0 = det(X). Para su(3), el determinante es puramente imaginario.
        detX = X(1,1)*(X(2,2)*X(3,3) - X(2,3)*X(3,2)) - &
               X(1,2)*(X(2,1)*X(3,3) - X(2,3)*X(3,1)) + &
               X(1,3)*(X(2,1)*X(3,2) - X(2,2)*X(3,1))
        detX_im = aimag(detX)
        
        ! Caso trivial: Si la matriz es cero
        if (c1 < 1.0e-12_wp) then
            expX = Identity
            return
        end if
        
        ! 3. Encontrar los autovalores analíticamente (son puramente imaginarios: lambda = i * tr_lambda)
        ! Usamos el método de Cardano para la ecuación cúbica: t^3 - c1*t - detX_im = 0
        q = c1 / 3.0_wp
        r = detX_im / 2.0_wp
        
        ! Evitar errores de redondeo en el acos
        if (abs(r / sqrt(q3(q))) > 1.0_wp) then
            theta = 0.0_wp
        else
            theta = acos(r / sqrt(q**3))
        end if
        
        ! Los tres autovalores de X son i * la_k
        la(1) = ii * 2.0_wp * sqrt(q) * cos(theta / 3.0_wp)
        la(2) = ii * 2.0_wp * sqrt(q) * cos((theta + 2.0_wp*3.141592653589793_wp) / 3.0_wp)
        la(3) = ii * 2.0_wp * sqrt(q) * cos((theta - 2.0_wp*3.141592653589793_wp) / 3.0_wp)
        
        ! 4. Resolver las componentes de la expansión: exp(X) = f0*I + f1*X + f2*X^2
        ! Usando interpolación de Lagrange para matrices
        f0 = (la(2)*la(3)*exp(la(1))) / ((la(1)-la(2))*(la(1)-la(3))) + &
             (la(1)*la(3)*exp(la(2))) / ((la(2)-la(1))*(la(2)-la(3))) + &
             (la(1)*la(2)*exp(la(3))) / ((la(3)-la(1))*(la(3)-la(2)))
             
        f1 = -(la(2)+la(3))*exp(la(1)) / ((la(1)-la(2))*(la(1)-la(3))) - &
              (la(1)+la(3))*exp(la(2)) / ((la(2)-la(1))*(la(2)-la(3))) - &
              (la(1)+la(2))*exp(la(3)) / ((la(3)-la(1))*(la(3)-la(2)))
              
        f2 = exp(la(1)) / ((la(1)-la(2))*(la(1)-la(3))) + &
             exp(la(2)) / ((la(2)-la(1))*(la(2)-la(3))) + &
             exp(la(3)) / ((la(3)-la(1))*(la(3)-la(2)))
        
        ! 5. Reconstruir la matriz final
        expX = f0 * Identity + f1 * X + f2 * X2
        
    contains
        ! Función auxiliar para evitar problemas de potencias flotantes
        pure real(wp) function q3(val)
            real(wp), intent(in) :: val
            q3 = val * val * val
        end function q3
    end subroutine exp_su3

end module su3_exponential

! --- Programa de prueba ---
program test_su3
    use su3_exponential
    implicit none
    
    complex(wp) :: X(3,3), expX(3,3)
    integer :: i
    
    ! Ejemplo: Una matriz antihermítica generada por las matrices de Gell-Mann (ej. i * lambda_1 + i * lambda_8)
    ! Debe cumplir que X^dagger = -X y Tr(X) = 0
    X = reshape([ (0.0_wp, 0.0_wp), (0.0_wp, 1.0_wp), (0.0_wp, 0.0_wp), &
                  (0.0_wp, 1.0_wp), (0.0_wp, 0.0_wp), (0.0_wp, 0.0_wp), &
                  (0.0_wp, 0.0_wp), (0.0_wp, 0.0_wp), (0.0_wp, 0.0_wp) ], [3,3])
                  
    ! Forzar traza cero y corregir un elemento para que sea su(3) real
    X(3,3) = (0.0_wp, 0.0_wp) 
    
    print *, "Matriz original X en su(3):"
    do i = 1, 3
        print '(3("(",F5.2,",",F5.2,") "))', X(i,:)
    end do
    
    call exp_su3(X, expX)
    
    print *, ""
    print *, "Exponencial exp(X) en SU(3):"
    do i = 1, 3
       print '(3("(",F10.7,",",F10.7,") "))', expX(i,:)
    end do

    expX = taylor_exp(X)
    print *, ""
    print *, "Exponencial exp(X) en SU(3):"
    do i = 1, 3
       print '(3("(",F10.7,",",F10.7,") "))', expX(i,:)
    end do

  contains

    function taylor_exp(X)
      complex(wp), dimension(3,3), intent(in) :: X
      complex(wp), dimension(3,3) :: taylor_exp, prodX
      integer :: i
      taylor_exp = (0.0_wp,0.0_wp)

      prodX = (0.0_wp,0.0_wp)
      prodX(1,1) = (1.0_wp,0.0_wp)
      prodX(2,2) = (1.0_wp,0.0_wp)
      prodX(3,3) = (1.0_wp,0.0_wp)
      Taylor_exp = prodX
      do i = 1, 50
         prodX = matmul(x,prodX)/i 
         taylor_exp = taylor_exp + prodX
      end do
      
    end function taylor_exp
     
end program test_su3

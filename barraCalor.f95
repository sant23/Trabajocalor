module calor

    use lineal 
    
    implicit none

    real(8), parameter:: pos1 = 493, pos2 = 293, km = 2.9, ka = 16.3

    contains

    subroutine caso1 (n, sol, h)

        implicit none

        integer, intent (in) :: n 
        real(8), intent (out) :: sol (n), h 
        real(8) :: D (n, n), A (n, n), B (n) 
        integer :: i

        D = 0d0; !Establecemos todos los valores a 0
        h = 1d0/(n-1) 
        do i = 1, n-1 !Derivadas centrales
            D (i, i+1) = 0.5d0
            D (i+1, i) = -0.5d0
        end do
        D (1,1) = -1d0; D (n, n) = 1d0; D (1, 2) = 1d0; D (n, n-1) = -1d0 ! Derivadas progresiva y regresiva
       
        D = ka*matmul (D, D)/h**2; A = 0d0 !Obtención de la matriz (por factor común se extrae h fuera)
        
        B = 0d0; A (2:n-1, 1:n) = D (2:n-1, :)!Añadimos todas las filas a la matriz Gaussiana menos la primera y la ultima
        A (1, 1) = 1d0; A (n, n) = 1d0; B (1) = pos1; B (n) = pos2 !Añadimos las condiciones de contorno en la primera y última 

        call gauss (A, B, sol) !Resolvemos el sistema

    end subroutine caso1


    subroutine caso2 (n, solDef, h)

        implicit none

        integer, intent (in) :: n !Numero de puntos por cada material
        real(8), intent (out) :: solDef (2*n-1), h !Solución con el número de puntos totales
        real(8) :: D (n, n), A (2*n, 2*n), B (2*n), sol (2*n) !Intérvalo, D y M, matriz A que contiene a las dos M, B del sistema, solucion inicial del sistema (con punto doble)
        integer :: i

        D = 0d0 !Establecemos todos los valores a 0
        h = 1d0/(2*n-2) 
        do i = 1, n-1 !Diagonales menores de la matriz D (Derivadas centrales)
            D (i, i+1) = 0.5d0 !Superior
            D (i+1, i) = -0.5d0 !Inferior
        end do
        
        D (1,1) = -1d0; D (n, n) = 1d0; D (1, 2) = 1d0; D (n, n-1) = -1d0 !Derivadas progresiva y regresiva
        
        D = matmul (D, D)/h**2   !Sacamos la matriz general M
        !Establecemos todo a 0, excepto los dos cuadrantes de la diagonal principal, cada uno con su k correspondiente, excepto la primera, la última fila y las dos del medio
        A = 0d0; B = 0d0
        A (2:n-1, 1:n) = ka*D (2:n-1, :); A (n+2:2*n-1, n+1:2*n) = km*D (2:n-1, :)
        A (n+1, n-1) = ka; A (n+1, n) = -ka; A (n+1, n+1) = -km; A (n+1, n+2) = km !El pto doble tiene que tener la misma derivada
        A (1, 1) = 1d0; A (2*n, 2*n) = 1d0; A (n, n) = 1d0; A (n, n+1) = -1d0; !Condiciones de contorno en A y Temperatura del punto medio es igual, T(n) = T(n+1) 
        B (1) = pos1; B (2*n) = pos2 !Condiciones de contorno
        
        call gauss (A, B, sol) !Resolvemos sistema
        do i = 1, n
            solDef (i) = sol (i) !Una vez obtenida la solución, eliminamos el punto medio, ya que este lo hemos establecido nosotros
            solDef (2*n-i) = sol (2*n-i+1)
        end do

    end subroutine caso2

end module calor
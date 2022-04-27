program trabajo

    use calor

    implicit none 
    
    real(8), allocatable :: Sol (:), Sol2(:)
    real(8) :: h
    integer :: n, i

    write (*,*) 'Introduce un numero impar de intervalos'
    read (*,*) n
    allocate (Sol (n), Sol2 (n))

    write(*,*)'la salida se generara en dos ficheros de texto'

    call caso1 (n, Sol, h)
    call caso2 (ceiling (n/2d0), Sol2, h)
    
    

    open (23, file = 'Temperaturas1Mat.txt') 
    write (23,*) 'x,   y'
    do i = 1, n
        write (23,*) (i-1)*h, ',',  Sol (i)
    end do
    close (23) 

    open (23, file = 'Temperaturas2Mat.txt') 
    write (23,*) 'x,   y'
    do i = 1, n
        write (23,*) (i-1)*h, ',',  Sol2 (i)
    end do
    close (23) 

    
end program trabajo
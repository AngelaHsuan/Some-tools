!--------------------------------------------------
subroutine Gaussian_eliminate (Mat, Nrow, Ncol,s,x)
    implicit NONE
    
    integer :: Nrow, Ncol, i, j, k
    real*8 :: Eli
    real*8, dimension (Nrow) :: s,x
    real*8, dimension (Nrow, Ncol) :: Mat
    real*8, dimension (Nrow, Ncol+1) :: Mat_new, Mat_old
    
    !Mat_old * x = s
    
    !---------------> first row
    !---------------> second row
    ! |       |
    ! |       |
    ! |       |
    !£¾      £¾
    !first  second
    !column column
    
    ! Nrow indicates the quantity of the matrix row
    ! Ncol indicates the quantity of the matrix column
    
    ! Do matrix collision
    
    do j = 1,Nrow,1
        do k = 1,Ncol,1
            Mat_old(j,k) = Mat(j,k)
        end do
        Mat_old(j,Ncol+1) = s(j)
    end do
    
    ! Do the elimination
    Mat_new = Mat_old
    do k = 1,Nrow-1,1
        do j = k+1,Nrow,1
            do i = k,Ncol+1,1
                Mat_new(j,i) = Mat_new(j,i) - Mat_old(k,i)*Mat_old(j,k)/Mat_old(k,k)
            end do
            Mat_old = Mat_new
        end do
    end do
    
    !do j = 1,Nrow
    !    do k = 1,Ncol+1
    !        write(*,*)Mat_old(j,k)
    !    end do
    !end do
    !pause
    
    ! Solve
    do j = 1,Nrow,1
        do k = 1,Ncol,1
            Mat(j,k) = Mat_old(j,k)
        end do
        s(j) = Mat_old(j,Ncol+1)
    end do
    
    x(Nrow) = s(Nrow)/Mat(Nrow,Ncol)
    do j = Nrow-1,1,-1
        Eli = 0.d0
        do i = j,Ncol-1,1
            Eli = Eli + Mat(j,i+1)*x(i+1)
        end do
        x(j) = (s(j)-Eli)/Mat(j,j)
    end do
    
end subroutine Gaussian_eliminate
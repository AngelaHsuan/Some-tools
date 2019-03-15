subroutine Matrix_multi(Mat_A,Mat_B,n,m,l,Mat_result)
    implicit NONE
    
    integer :: n,m,l,i,j,k  ! Matrix A has dimension n*m ; Matrix B has dimension m*l
    real*8, dimension (n,m) :: Mat_A
    real*8, dimension (m,l) :: Mat_B
    real*8, dimension (n,l) :: Mat_result
    
    Mat_result(:,:) = 0
    do i = 1,n,1
        do j = 1,l,1
            do k = 1,m,1
                Mat_result(i,j) = Mat_result(i,j) + Mat_A(i,k)*Mat_B(k,j)
            end do
        end do
    end do
    
end subroutine Matrix_multi
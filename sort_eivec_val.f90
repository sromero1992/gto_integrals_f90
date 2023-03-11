subroutine sort_eigvec_val(eig,eigv,nbas)
  real(8),intent(inout)  :: eig(nbas), eigv(nbas,nbas)
  integer,intent(in)     :: nbas  
  real(8),allocatable    :: eig_tmp(:),eigv_tmp(:,:)
  integer,allocatable    :: indx_array(:) 
  integer :: nbas,i  
  
  allocate(eig_tmp(nbas))
  allocate(eigv_tmp(nbas,nbas))
  allocate(indx_array(nbas))
  eig_tmp = eig
  eigv_tmp= eigv

  !write(*,*)'index=i, eig(i) before call' 
  do i=1,nbas
     indx_array(i)=i
     !write(*,*) indx_array(i), eig(indx_array(i))
  end do
 
  call quicksort2(eig_tmp,indx_array,1,nbas)
  
  !write(*,*)'i, indx_array, eig(i)  after call' 
  do i=1,nbas
     !write(*,*) i,indx_array(i), eig(indx_array(i))
     eig(i) = eig_tmp(i)
     eigv(:,i) = eigv_tmp(:, indx_array(i) )
  end do
 

end subroutine

recursive subroutine quicksort2(array,indx_array, first, last)
  implicit none
  real(8)  :: array(*), x, t
  integer  :: first, last, indx_array(*), i, j, t2

  x = array( (first+last) / 2 )
  i = first
  j = last
  do
     do while (array(i) < x)
        i=i+1
     end do
     do while (x < array(j))
        j=j-1
     end do
     if (i >= j) exit
     t = array(i);  array(i) = array(j);  array(j) = t
     t2 = indx_array(i); indx_array(i)=indx_array(j); indx_array(j)=t2
     i=i+1
     j=j-1
  end do
  if (first < i-1) call quicksort2(array,indx_array, first, i-1)
  if (j+1 < last)  call quicksort2(array,indx_array, j+1, last)
end subroutine 


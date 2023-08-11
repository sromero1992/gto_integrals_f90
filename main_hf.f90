program basis_test
  use class_basis_function
  use integral
  use module_com
  use module_g
  use mat_build
  use diag
  implicit none
  real(8), allocatable    :: S_TMP(:,:)

  jobz='V'
  uplo='U'


  del          = 1.0d-7
  e_tol        = 1.0d-6
  chk_nrlmol   = .FALSE.
  debug        = .FALSE.
  ext_bas      = .FALSE.
  f_bol        = .FALSE. ! Activate F functions
  SPDF(:)      = (/1, 3, 6, 10/) 
  SPDF_SUM(:)  = (/0, 1, 4, 10, 20/)    


 
  call shell_gen(shells)
  call read_info
   
  if (f_bol) then 
     l_max = 4
     nbas = sum(g(:)%sf + 3*g(:)%pf + 6*g(:)%df + 10*g(:)%ff )
  else 
     l_max = 3
     nbas = sum(g(:)%sf + 3*g(:)%pf + 6*g(:)%df )
  end if

  !call basis_moments
  call okven_matbuild

  allocate(S_TMP(nbas,nbas))
  !call diagon_fortran(1,.TRUE.,S_TMP) 
  call diagon(1,.FALSE.,S_TMP) 

  !call hf_scf(ener_final)
  !write(*,*) 'FINAL ENERGY : ', ener_final

end program 


!#########################################################################
subroutine shell_gen(shells)
   ! shells define the angular momentum
   ! shells(li, max_nfunc,xyz)
   ! li angular momentum li = 1, 2, 3, 4  (S, P, D, F)
   ! max_nfunc = max( 1, 3, 6, 10) depending on basis e.g. up to D, max_nfunc = 6
   ! xyz =>  1, 2, 3
   integer,intent(inout)   :: shells(4,10,3)
   !integer,intent(inout)   :: shells(4,10,3)
   integer                 :: ii
   character(10)           :: chr_tmp
   logical                 :: read_file

   read_file = .FALSE.

   shells(:,:,:) = 0

  ! P-shell Px => (1,0,0) 
   shells(2,1,:) = (/1,0,0/)
  ! P-shell Py => (0,1,0)
   shells(2,2,:) = (/0,1,0/)
  ! P-shell Py => (0,0,1)
   shells(2,3,:) = (/0,0,1/)


  ! D-shell Dxx => (2,0,0)
   shells(3,1,:) = (/2,0,0/)
  ! D-shell Dyy => (0,2,0)
   shells(3,2,:) = (/0,2,0/)
  ! D-shell Dzz => (0,0,2)
   shells(3,3,:) = (/0,0,2/)
  ! D-shell Dxy => (1,1,0)
   shells(3,4,:) = (/1,1,0/)
  ! D-shell Dxz => (1,0,1)
   shells(3,5,:) = (/1,0,1/)
  ! D-shell Dyz => (0,1,1)
   shells(3,6,:) = (/0,1,1/)
  

  ! F-shell Dxxx => (3,0,0)
   shells(4,1,:) = (/3,0,0/)
  ! D-shell Dyyy => (0,3,0)
   shells(4,2,:) = (/0,3,0/)
  ! D-shell Dzzz => (0,0,3)
   shells(4,3,:) = (/0,0,3/)
  ! D-shell Dxxy => (2,1,0)
   shells(4,4,:) = (/2,1,0/)
  ! D-shell Dxxz => (2,0,1)
   shells(4,5,:) = (/2,0,1/)
  ! D-shell Dxyy => (1,2,0)
   shells(4,6,:) = (/1,2,0/)
  ! D-shell Dyyz => (0,2,1)
   shells(4,7,:) = (/0,2,1/)
  ! D-shell Dxxz => (1,0,2)
   shells(4,8,:) = (/1,0,2/)
  ! D-shell Dxyy => (0,1,2)
   shells(4,9,:) = (/0,1,2/)
  ! D-shell Dxyy => (1,1,1)
   shells(4,10,:) = (/1,1,1/)
   

   if (read_file) then
      open(21,file='shells.dat')
      read(21,*) chr_tmp
      do ii = 1, 3
         read(21,*) shells(2,ii,1), shells(2,ii,2), shells(2,ii,3) 
      end do
      read(21,*) chr_tmp
      do ii = 1, 6
         read(21,*) shells(3,ii,1), shells(3,ii,2), shells(3,ii,3) 
      end do
      read(21,*) chr_tmp
      do ii = 1, 10
         read(21,*) shells(4,ii,1), shells(4,ii,2), shells(4,ii,3) 
      end do
   end if
end subroutine




subroutine AV_VD_test(V,D,A,n)
   !V = eigV
   !D = eig 
   !A = F_MAT
   implicit none
   integer             :: i, j, k, n
   real(8), intent(in) :: V(n,n), D(n), A(n,n)
   real(8)             :: max_dif, S1(n,n), S2(n,n), E(n,n)
  
   write(*,*)
   write(*,*) 'MATRIX V = eigV :'
   do i = 1, n
      write(*,*) V(i,:)
   end do 
   write(*,*)

   write(*,*)
   write(*,*) 'VECTOR D = eig :'
   do i = 1, n
      write(*,*) D(i)
   end do 
   write(*,*)
   
   write(*,*)
   write(*,*) 'MATRIX A = S/F MAT :'
   do i = 1, n
      write(*,*) A(i,:)
   end do 
   write(*,*)
 
   S1 = 0.0d0 
   do i = 1, n
      do j = 1, n
         do k = 1, n
         S1(i,j) = S1(i,j) +  A(i,k)*V(k,j)
         end do
      end do
   end do   
   write(*,*)
   write(*,*) 'S1 = AV :'
   do i = 1, n
      write(*,*) S1(i,:)
   end do  
   write(*,*)

   write(*,*)
   S2 = 0.0d0 
   do i = 1, n
      S2(:,i) = S2(:,i) + V(:,i)*D(i)
   end do   
   write(*,*) 'S2 = VD :'
   do i = 1, n
      write(*,*) S2(i,:)
   end do  
   write(*,*)



   E=S1-S2
   max_dif = 0 
   write(*,*)
   write(*,*) 'HISTORY OF MAX_DIF'
   do i = 1, n
      do j = 1, n
         if ( abs( E(i,j) ) .GT. abs(max_dif) ) then
            max_dif = E(i,j)
            write(*,*) max_dif
         end if
      end do
   end do
   write(*,*) 'MAX MAT DIFFERENCE:'
   write(*,*) max_dif
   write(*,*)
end subroutine

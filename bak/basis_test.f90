program basis_test
  use class_basis_function
  use integral
  use com_vars
  use mat_build
  use diag
  implicit none
  ! type(gaussian)         :: g_tmp, ga, gb, gc, gd     ! Declare gaussian basis 
  !type(integ)            :: Ven ! Declaration of integral variables that uses g objects
  real(8), allocatable   :: S_TMP(:,:), X_TMP(:,:)
  real(8)                :: ener(30), nuc_attraction, ener_rep 


  integer                :: mx_scf, n_scf, ia
  logical                :: conv

  jobz='V'
  uplo='U'


  del          = 1.0d-7
  ext_bas      = .FALSE.
  f_bol        = .FALSE. ! Activate F functions
  SPDF(:)      = (/1, 3, 6, 10/) 
  SPDF_SUM(:)  = (/1, 4, 10, 20/)    


 
  call shell_gen(shells)
  call read_info
   
  if (f_bol) then 
     l_max = 4
     nbas = sum(g(:)%sf + 3*g(:)%pf + 6*g(:)%df + 10*g(:)%ff )
  else 
     l_max = 3
     nbas = sum(g(:)%sf + 3*g(:)%pf + 6*g(:)%df )
  end if

  call okven_matbuild


  !allocate(TMP_MAT(nbas,nbas))
  
  allocate(S_TMP(nbas,nbas))
  call diagon(1,.FALSE.,S_TMP) 
  !S_TMP = S_MAT
  !call diagon(2,.FALSE.,S_TMP)
  !allocate( COUL_MAT( nbas, nbas) )
  !call coul_build  

  !################# preparing SCF cycle ########################3
  write(*,*)
  write(*,*)' PREPATATION FOR SCF CYCLE'

  !S transformation
  allocate(X_TMP(nbas,nbas))
  call s_transform(1,.FALSE.,X_TMP)
  !call s_transform(2,.FALSE.,X_TMP)

  !write(*,*)
  !write(*,*) 'SUBROUTINE MAT-TEST FOR AV = VD'
  !call AV_VD_test(eigV,eig,S_MAT,nbas)




  allocate(H0_MAT(nbas,nbas))
  allocate( F_MAT(nbas,nbas))
  allocate( Fp_MAT(nbas,nbas))
  allocate( F_TMP(nbas,nbas))
  allocate( F_eig(nbas))
  allocate( F_eigV(nbas,nbas))
  allocate( COUL_MAT( nbas, nbas) )
  allocate( COUL_VEC( nbas*(nbas+1)/2) )

  H0_MAT = T_MAT + Ven_MAT 
  F_MAT  = H0_MAT

  call diagon(3,.TRUE.,F_MAT)

  write(*,*)
  write(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%% SCF CYCLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
  conv   = .FALSE.
  n_scf  = 1
  mx_scf = 30
  do while ( (n_scf < mx_scf) .AND. .NOT. conv )
    write(*,*) 'ITERATION :',n_scf
    n_scf = n_scf + 1
    call coul_build  
    write(*,*)
    !write(*,*) 'Vee MATRIX :'
    !do i = 1, nbas
    !   write(*,*) COUL_MAT(i,:)
    !end do 
    F_MAT = H0_MAT + COUL_MAT
    !F' rotation of Fock operatior to molecular basis 
     
    call diagon(3,.TRUE.,F_MAT)

    !calculate energy
    call get_energy(ener(n_scf))
    write(*,*) 'ENERGY IS : ', ener(n_scf)
    !ener(n_scf) = n_scf
    if ( (n_scf > 2) .AND. ( abs( ener(n_scf) - ener(n_scf-1) ) .LT. 1.0d-6 )  ) then
        write(*,*) 'ENERGY CONVERGET AT ITERATION : ', n_scf
        write(*,*) 'ENERGY : ',ener(n_scf)
        conv = .TRUE.
    end if
  end do

  call nuc_rep(ener_rep) 

  ener(n_scf) = ener(n_scf) + ener_rep
  write(*,*) 'FINAL ENERGY : ', ener(n_scf)  

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

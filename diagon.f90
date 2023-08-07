subroutine  diagon(mode,verbose,D_TMP) 
  !mode 1 calculates density for S_MAT
  !mode 2 calculates density for matrix D_TMP 
  !mode 3 calculates rotation for molecular orbitals F' 
  use module_com
  use diag 
  implicit none
  real(8),intent(inout) :: D_TMP(nbas,nbas)
  integer               :: mode
  logical               :: verbose
  jobz='V'
  uplo='U'
  
  ! This will be replaced by sewell subroutines
  ! DEGSYM
  allocate(WORK(1))
  if ( .not. allocated(eig     ) ) allocate(    eig(nbas))
  if ( .not. allocated(eig_tmp ) ) allocate(eig_tmp(nbas))
  if ( .not. allocated(eigV    ) ) allocate(    eigV(nbas,nbas))
  if ( .not. allocated(eigV_tmp) ) allocate(eigV_tmp(nbas,nbas))
  if ( .not. allocated(D_mat   ) ) allocate(   D_MAT(nbas,nbas))


  if (mode .EQ. 1) then
     write(*,*) 'DIAGONALIZE S_MAT'

     eigv = S_MAT
     LWORk=-1
    
     write(*,*) 'SIZE OF S_MAT:',size(S_MAT)
     CALL DSYEV( JOBZ, UPLO, NBAS, eigV, NBAS, eig, WORK, LWORK, INFO )
     LWORK= INT(WORK(1))
     write(*,*)'WORK QUERY (LWORK)',LWORK
     deallocate(WORK)
     allocate(WORK(LWORK),STAT=ierr)
   
     CALL DSYEV( JOBZ, UPLO, NBAS, eigV, NBAS, eig, WORK, LWORK, INFO )
   
     write(*,*)'INFO : ',info
     if (info .EQ. 0 ) then
       write(*,*)'DSYEV finished successfully'
     else if (info .GT. 0) then
       write(*,*)'DSYEV faild to converge'
     else
       write(*,*)'DSYEV info=-i, argument had an illegal value'
     end if
   
   
     call sort_eigvec_val(eig, eigV, nbas)
     
     if (verbose) then  
        write(*,*)'EVALUES OF S_MAT : '
        do i = 1, nbas
           write(*,*) i, eig(i)
        end do
        write(*,*)
        write(*,*)'EVECTORS OF S_MAT (COLUMN WISE) :'
        do i = 1, nbas
           write(*,*)eigV(i,:)
           write(*,*)
        end do 
     end if 
     deallocate(work)
     write(*,*) 'END OF DIAG FOR S_MAT'

  else if ( mode .EQ. 2) then

     write(*,*) 'DIAGONALIZE INPUT MATRIX D_TMP '
     eigV_tmp = D_TMP
     LWORk=-1
    
     write(*,*) 'SIZE OF D_TMP MATRIX :',size(D_TMP)
     CALL DSYEV( JOBZ, UPLO, NBAS, eigV_tmp, NBAS, eig_tmp, WORK, LWORK, INFO )
     LWORK= INT(WORK(1))
     write(*,*)'WORK QUERY (LWORK)',LWORK
     deallocate(WORK)
     allocate(WORK(LWORK),STAT=ierr)
   
     CALL DSYEV( JOBZ, UPLO, NBAS, eigV_tmp, NBAS, eig_tmp, WORK, LWORK, INFO )
   
     write(*,*)'INFO : ',info
     if (info .EQ. 0 ) then
       write(*,*)'DSYEV finished successfully'
     else if (info .GT. 0) then
       write(*,*)'DSYEV faild to converge'
     else
       write(*,*)'DSYEV info=-i, argument had an illegal value'
     end if
   
   
     call sort_eigvec_val(eig_tmp, eigV_tmp, nbas)
   
     if (verbose) then 
        write(*,*)'EVALUES OF D_TMP MATRIX : '
        do i = 1, nbas
           write(*,*) i, eig_tmp(i)
        end do
        write(*,*)
        write(*,*)'EVECTORS OF D_TMP MATRIX : '
        do i = 1, nbas
           write(*,*)eigV_tmp(i,:)
        end do 
     end if
     deallocate(work)
     write(*,*) 'END OF DIAG  D_TMP' 
  
  else if (mode .EQ. 3) then


   !Transformation to get molecular orbitals
   !F'=X^(T)*F*X
   if (verbose) then 
      write(*,*)
      write(*,*) '################### F_MAT #######################'
      do i = 1, nbas
         write(*,*)F_MAT(i,:)
      end do
   end if

   Fp_MAT  = 0.0d0
   do i = 1, nbas
      do j = 1, nbas
         do l = 1, nbas
            do k =1, nbas
               Fp_MAT(i,j) = Fp_MAT(i,j) + X_MAT(l,i)*F_MAT(l,k)*X_MAT(k,j)
            end do
         end do
      end do
   end do

   if (verbose) then 
      write(*,*)
      write(*,*) '################### F_MAT transformed (Fp_MAT) #######################'
      do i = 1, nbas
         write(*,*)Fp_MAT(i,:)
      end do
   end if

   F_eigV = 0.0d0
   F_eigV = Fp_MAT
   LWORK=-1
   !allocate(WORK(1))
   write(*,*)
   write(*,*) 'DSYEV TO FIND EVALUES AND EVECTORS:'
   write(*,*) 'SIZE OF Fp_MAT, NBAS:',size(Fp_MAT),NBAS
   CALL DSYEV( JOBZ, UPLO, NBAS, F_eigV, NBAS, F_eig, WORK, LWORK, INFO )
   LWORK= INT(WORK(1))
   write(*,*)'WORK QUERY (LWORK)', LWORK
   deallocate(WORK)
   allocate(WORK(LWORK),STAT=ierr)
 
   CALL DSYEV( JOBZ, UPLO, NBAS, F_eigV, NBAS, F_eig, WORK, LWORK, INFO )
   deallocate(WORK)
 
   write(*,*)'INFO : ',info
   if (info .EQ. 0 ) then
     write(*,*)'DSYEV finished successfully'
   else if (info .GT. 0) then
     write(*,*)'DSYEV faild to converge'
   else
     write(*,*)'DSYEV info=-i, argument had an illegal value'
   end if
 
   if (verbose) then
      write(*,*)
      write(*,*)'################### EVALUES OF Fp_MAT ####################: '
      do i = 1, nbas
         write(*,*) i, F_eig(i)
      end do
      write(*,*)
      !Transformed eigen vectors of Fp_MAT to F_MAT space
      !write(*,*)
      !write(*,*) '################## EIGENVECTORS OF Fp_MAT (Cp)  ####################'
      !do i = 1, nbas
      !   write(*,*) F_eigV(i,:)
      !end do 
   end if 

   F_TMP=0.0d0
   do i = 1, nbas
      do j =1, nbas
         do k = 1, nbas
            F_TMP(i,j) = F_TMP(i,j) + X_MAT(i,k)*F_eigV(k,j)
         end do
      end do
   end do
   !Molecular orbitals
   F_eigV = F_TMP

   if (verbose) then
      write(*,*)
      write(*,*) '################## EVECTORS OF Fp_MAT (C)  ####################' 
      do i = 1, nbas
         write(*,*)F_eigV(i,:)
      end do
      write(*,*)
   end if

   !Build density here for SCF cycle
   call dens(1, verbose, F_eigV, D_MAT)
   !write(*,*)
   !write(*,*) 'SUBROUTINE MAT-TEST FOR AV = VD'
   !call AV_VD_test(F_eigV,F_eig,Fp_MAT,nbas)
   !write(*,*)
   
  end if

end subroutine 


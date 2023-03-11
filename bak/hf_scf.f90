subroutine hf_scf(ener_fin)
   !use class_basis_function
   !use integral
   use com_vars
   use mat_build
   use diag
   real(8), intent(inout)  :: ener_fin
   real(8), allocatable    :: X_TMP(:,:)
   real(8)                 :: ener(30), ener_rep
   integer                 :: mx_scf, n_scf
   logical                 :: conv

    
   !################# preparing SCF cycle ########################3
   write(*,*)
   write(*,*)' PREPATATION FOR SCF CYCLE'
 
   !S transformation
   allocate(X_TMP(nbas,nbas))
   call s_transform(1,.FALSE.,X_TMP)
   !call s_transform(2,.FALSE.,X_TMP)
 
   !write(*,*)
   !write(*,*) 'SUBROUTINE MAT-TEST (EIGENVALUE TEST) FOR AV = VD'
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
 
   call diagon(3,.FALSE.,F_MAT)
 


   write(*,*)
   write(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%% SCF CYCLE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
   conv   = .FALSE.
   n_scf  = 1
   mx_scf = 30
   do while ( (n_scf < mx_scf) .AND. .NOT. conv )
     write(*,*) 'ITERATION :',n_scf
     n_scf = n_scf + 1
     write(*,*) 'BUILDING COULOM INIT'
     call coul_build

     F_MAT = H0_MAT + COUL_MAT
     !F' rotation of Fock operatior to molecular basis

     call diagon(3,.FALSE.,F_MAT)

     !calculate energy
     call get_energy(ener(n_scf))
     write(*,*) 'ENERGY IS : ', ener(n_scf)
     if (debug) exit
     if ( (n_scf > 2) .AND. ( abs( ener(n_scf) - ener(n_scf-1) ) .LT. del )  ) then
         write(*,*) 'ENERGY CONVERGET AT ITERATION : ', n_scf
         write(*,*) 'ENERGY : ',ener(n_scf)
         conv = .TRUE.
     end if
   end do

   call nuc_rep(ener_rep)
   ener(n_scf) = ener(n_scf) + ener_rep
   write(*,*) 'FINAL ENERGY : ', ener(n_scf)
   ener_fin = ener(n_scf)

end subroutine

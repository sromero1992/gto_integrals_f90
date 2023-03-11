program read_isymgen
implicit none
integer              :: tot_natom, i, j, k, l, m,  ai, bj, kk, nuc_chrg, elec_chrg, sup_sf, &
                        sup_pf, sup_df, nb, ib, jb, ext_base, nbas, nbas2, nrec, ibm, jbm
integer, allocatable :: atm_type(:), natm_type(:), ng(:), sf(:), pf(:), df(:)
real(8), allocatable :: alp(:,:),cg(:,:,:),r_atm(:,:), ovl_mat(:,:),ovlbaby(:)
real(8)              :: elec(2), ovl_prod(3), ovl_prodSP(3), nrec2, gp3d, PijAB(3), pij, &
                        dum3(3), del
character (7)        :: all_elec, atm_type_str, end_str, atm_str, dum_str
logical              :: le, verbose, not_ext_bas

verbose   = .TRUE.
not_ext_bas   = .TRUE.
del = 1.0d-7

open(11,file='ISYMGEN')
read(11,*) tot_natom !Equivalent atoms
allocate(atm_type(tot_natom))
allocate(natm_type(tot_natom))
allocate(ng(tot_natom))
allocate(sf(tot_natom))
allocate(pf(tot_natom))
allocate(df(tot_natom))
!allocate(nb(tot_natom))


do i=1,tot_natom
   read(11,*) elec_chrg, nuc_chrg !of current atom
   read(11,*) all_elec
   read(11,*) natm_type(i)
   read(11,*) atm_type_str
   !this reads ALL-ATM & EXTRA BASIS label
   if (natm_type(i) .GT. 1) then
     do j =1, natm_type(i)
        read(11,*) 
     end do
   end if
   read(11,*) ng(i) !Number of gaussians
   read(11,*) sf(i), pf(i), df(i) !Gaussian functions S,P,D
   read(11,*) sup_sf, sup_pf, sup_df !Extra basis/supplemental of S,P,D
   !allocate(nb(tot_natom)) you have different number of basis for each different atom
   nb =sf(i)+pf(i) !*3!+6*df(i)!+sup_sf+sup_pf+sup_df !Total number of functions
   !nb =sf(i)+pf(i)+df(i)+sup_sf+sup_pf+sup_df !Total number of functions
   allocate( alp( ng(i), natm_type(i) ) ) 
   allocate( cg( ng(i), nb, natm_type(i) ))
   !Reading all alpha exponents
   do j=1,ng(i)/3
      read(11,*) alp(3*j-2,natm_type(i)), alp(3*j-1,natm_type(i)), alp(3*j,natm_type(i))
   end do
   !Reading all gaussian contraction coefficients Cij
   do ib=1,nb
      do j=1,ng(i)/3
         read(11,*) cg(3*j-2,ib,natm_type(i)), cg(3*j-1,ib,natm_type(i)), cg(3*j,ib,natm_type(i))
      end do
   end do
   !primitive gaussian functions psi_i = sum_j (Cij*phi_j)
end do
close(11)
!read(11,*)
!read(11,*) end_str
!if (end_str == "ELECTRO") close(11)

!Print out all the read
if (verbose) then
write(*,*) 'Verbose ISYMGEN file'
i=1
!do i=1,tot_natom
  write(*,*) tot_natom
  write(*,*) elec_chrg, nuc_chrg
  write(*,*) all_elec
  write(*,*) natm_type(i)
  write(*,*) atm_type_str
  write(*,*) ng(i)
  write(*,*) sf(i), pf(i), df(i) 
  write(*,*) sup_sf, sup_pf, sup_df 
  !do j=1,ng(i)/3
  !   write(*,*) alp(3*j-2,natm_type(i)), alp(3*j-1,natm_type(i)), alp(3*j,natm_type(i))
  !end do
  !do ib=1,nb !i-th function or the whole basis function (sum_i Ci*exp(-alp_i*|r-A|^2))
  !  do j=1,ng(i)/3
  !     write(*,*) cg(3*j-2,ib,natm_type(i)), cg(3*j-1,ib,natm_type(i)), cg(3*j,ib,natm_type(i))
  !  end do
  !end do
  write(*,*) 'Gaussian exponents'
  write(*,*) alp(:,natm_type(i))
  write(*,*)
  write(*,*) 'Contraction coefficients'
  do ib=1,nb
       write(*,*) cg(:,ib,natm_type(i))
       write(*,*)
  end do
!end do
end if

!nb =sf(i)+pf(i)*3!+6*df(i)!+sup_sf+sup_pf+sup_df !Total number of functions
!alp are the exponent arguments and cg are the gaussian coefficients of primitive Gaussians

!Atom positions and SYMBOL information
allocate(r_atm(3,natm_type(i)))
open(12,file='SYMBOL')
write(*,*) 'Atom positions and number of electrons'
do i=1,10
   read(12,*)
end do
j=1
do i=1, natm_type(j)
   read(12,*)  atm_str, dum_str, r_atm(1,i), r_atm(2,i), r_atm(3,i), dum_str
end do
read(12,*) dum_str, dum_str, elec(1), elec(2)
read(12,*) dum_str, dum_str, ext_base
close(12)

!Print read file info
if (verbose) then
  write(*,*)
  write(*,*) r_atm(:,1)
  write(*,*) r_atm(:,2)
  write(*,*) elec(:)
  write(*,*) ext_base
  write(*,*)
end if

!Calc OVERLAP
write(*,*) 'OVL matrix'
l=1 !Atom type 1
!do l=1,natm_type(1) 
   write(*,*) 'values of s and p bare gaussians :',sf(l),pf(l)
   nb = sf(l)+3*pf(l)
   write(*,*) 'basis functions :',nb
   allocate( ovl_mat( nb, nb) ) !sf -> nb nbas, so nb = sum_i { [sf(i)+ pf(i) + df(i)]+ suplemental}
   ovl_mat(:,:) =0.0d0
   do ai = 1, natm_type(l) !atom position A: nbasis in A
    do bj = ai+1, natm_type(l) !atom position B: nbasis in B
         do ib = 1, nb ! basis i
          do jb = 1, nb ! basis j
             !Loop to get the integral of bare gaussians
             do i = 1, ng(l) !basis j of bare gaussians
              do j = 1, ng(l) !basis k of bare gaussians
                 gp3d = 0.0d0
                 call D1_gprod( alp(i, natm_type(l) ), alp(j, natm_type(l) ), r_atm(:,ai), r_atm(:,bj), ovl_prod(:))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            S-S OVL              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 if ( ( ib .LE. sf(l) ) .AND. ( jb .LE. sf(l) ) ) then
                     gp3d = ovl_prod(1)*ovl_prod(2)*ovl_prod(3)*cg(i, ib, natm_type(l) )* cg(j, jb, natm_type(l))
                 end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         S-P/P-S  OVL              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   if ( ( ib .LE. sf(l) ) .AND. ( jb .GT. sf(l)) ) then !loop over the p functions S-Pi
                   jbm = mod(jb-sf(l),pf(l))
                   if (jbm .EQ. 0) jbm=3 
                   gp3d = ovl_prod(1)*ovl_prod(2)*ovl_prod(3)*cg(i, ib, natm_type(l) )* cg(j, jbm, natm_type(l))
                   pij = alp(i, natm_type(l) ) + alp(j, natm_type(l) )
                   PijAB(:) = ( alp(i, natm_type(l) )*r_atm(:,ai) + alp(j, natm_type(l) )*r_atm(:,bj) )/pij

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            S-P OVL              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                   if (jb .LE. sf(l) + pf(l) ) then
                     gp3d = gp3d*(PijAB(1)-r_atm(1,bj)) !for S-Px functions
                   else if ( ( jb .LE. sf(l) + 2*pf(l) ) ) then !.AND. ( jb .GT. sf(l)+pf(l)) ) then
                     gp3d = gp3d*(PijAB(2)-r_atm(2,bj)) !for S-Py functions
                   else if ( ( jb .LE. sf(l) + 3*pf(l) ) ) then ! .AND. ( jb .GT. sf(l)+2*pf(l)) ) then
                     gp3d = gp3d*(PijAB(3)-r_atm(3,bj)) !for S-Pz functions
                   end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            P-S OVL              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 else if ( (jb .LE. sf(l) ) .AND. (ib .GT. sf(l)) ) then !loop over the p functions in Pi-S 
                   ibm = mod(ib-sf(l),pf(l))
                   if (ibm .EQ. 0) ibm=3 
                   gp3d = ovl_prod(1)*ovl_prod(2)*ovl_prod(3)*cg(i, ibm, natm_type(l) )* cg(j, jb, natm_type(l))
                   pij = alp(i, natm_type(l) ) + alp(j, natm_type(l) )
                   PijAB(:) = ( alp(i, natm_type(l) )*r_atm(:,ai) + alp(j, natm_type(l) )*r_atm(:,bj) )/pij
                   !!!!!!! Px-S  !!!!!
                   if ( ib .LE. sf(l) + pf(l) ) then
                     gp3d = gp3d*(PijAB(1)-r_atm(1,ai)) !for Px-S functions 
                   !!!!!!! Py-S  !!!!!
                   else if ( ( ib .LE. sf(l) + 2*pf(l) ) ) then
                     gp3d = gp3d*(PijAB(2)-r_atm(2,ai)) !for Py-S functions
                   !!!!!!! Pz-S  !!!!!
                   else if ( ( ib .LE. sf(l) + 3*pf(l) ) ) then
                     gp3d = gp3d*(PijAB(3)-r_atm(3,ai)) !for Pz-S functions
                   end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            P-P OVL              !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                 else if ( (jb .GT. sf(l) ) .and. (ib .GT. sf(l) ) ) then
                   ibm = mod(ib-sf(l),pf(l))
                   if (ibm .EQ. 0) ibm=3 
                   jbm = mod(jb-sf(l),pf(l))
                   if (jbm .EQ. 0) jbm=3 
                   gp3d = ovl_prod(1)*ovl_prod(2)*ovl_prod(3)*cg(i, ibm, natm_type(l) )* cg(j, jbm, natm_type(l))
                   pij = alp(i, natm_type(l) ) + alp(j, natm_type(l) )
                   PijAB(:) = ( alp(i, natm_type(l) )*r_atm(:,ai) + alp(j, natm_type(l) )*r_atm(:,bj) )/pij
                   !!!!! Px-Pi !!!!!!!
                   if ( ib .LE. sf(l) + pf(l) ) then
                       if (jb .LE. sf(l) + pf(l) ) then
                         gp3d = gp3d*(0.5d0/pij + PijAB(1)**2.0d0 +r_atm(1,ai)*r_atm(1,bj) ) !for Px-Px functions  
                       else if (jb .LE. sf(l) + 2*pf(l) ) then
                         gp3d = gp3d*( PijAB(1)-r_atm(1,ai) )*( PijAB(2)-r_atm(2,bj) ) !for Px-Py functions
                       else if (jb .LE. sf(l) + 3*pf(l) ) then
                         gp3d = gp3d*( PijAB(1)-r_atm(1,ai) )*( PijAB(3)-r_atm(3,bj) ) !for Px-Pz functions
                       end if

                   !!!!! Py-Pi !!!!!!!
                   else if ( ib .LE. sf(l) + 2*pf(l) ) then
                       if (jb .LE. sf(l) + pf(l) ) then
                         gp3d = gp3d*( PijAB(2)-r_atm(2,ai) )*( PijAB(1)-r_atm(1,bj) ) !for Py-Px functions  
                       else if (jb .LE. sf(l) + 2*pf(l) ) then
                         gp3d = gp3d*(0.5d0/pij + PijAB(2)**2.0d0 +r_atm(2,ai)*r_atm(2,bj) ) !for Py-Py functions  
                       else if (jb .LE. sf(l) + 3*pf(l) ) then
                         gp3d = gp3d*( PijAB(2)-r_atm(2,ai) )*( PijAB(3)-r_atm(3,bj) ) !for Py-Pz functions
                       end if

                   !!!!! Pz-Pi !!!!!!!
                   else if ( ib .LE. sf(l) + 3*pf(l) ) then
                       if (jb .LE. sf(l) + pf(l) ) then
                         gp3d = gp3d*( PijAB(3)-r_atm(3,ai) )*( PijAB(1)-r_atm(1,bj) ) !for Pz-Px functions  
                       else if (jb .LE. sf(l) + 2*pf(l) ) then
                         gp3d = gp3d*( PijAB(3)-r_atm(3,ai) )*( PijAB(2)-r_atm(2,bj) ) !for Pz-Py functions 
                       else if (jb .LE. sf(l) + 3*pf(l) ) then
                         gp3d = gp3d*(0.5d0/pij + PijAB(3)**2.0d0 +r_atm(3,ai)*r_atm(3,bj) ) !for Pz-Pz functions  
                       end if

                   end if
                 end if
                 ovl_mat(ib,jb) = ovl_mat(ib,jb) + gp3d
              end do
             end do
          end do
          !write(*,*) ovl_mat(ib,:)
          !write(*,*)
          !write(*,*)
         end do
    end do
   end do
!end do


!Writting overlap results
write(*,*)
inquire(file='ovl.txt',exist=le)
if (le) then
  open(13,file='ovl.txt')
  close(13,status='delete')
end if
i=1
open(13,file='ovl.txt',status='unknown')
do ib = 1,nb
    write(13,*) ovl_mat(ib,:)
end do
close(13)



!Check against ovlbaby
open(14,file='OVLBABY',form='UNFORMATTED')
read(14) nrec
write(*,*) 'nrec file :',nrec
allocate(ovlbaby(nrec))
read(14) (ovlbaby(j),j=1,nrec)
close(14)
open(15,file='ovl_ovlbaby.txt')

l=0
m=0
do j=1, nrec
   if (ovlbaby(j) .NE. 0.0d0) then 
      l = l+1
      do i=1, nb
         do k=1, nb
            if (ovl_mat(i,k) .NE. 0.0d0) then  
               if (  abs( ovlbaby(j) - ovl_mat(i,k) ) .LE.  del  ) then
                  write(15,*)ovlbaby(j),j,ovl_mat(i,k),i,k, '     X'
                  write(*,*) ovlbaby(j),j
                  m = m+1
               else
                  write(15,*)ovlbaby(j),j,ovl_mat(i,k),i,k
               end if
            end if
         end do
      end do
   end if
end do
close(15)
write(*,*) 'OVLBABY total elements :',nrec
write(*,*) 'OVLBABY non-zero elements :', l 
write(*,*) 'ovl similar non-zero elements elements :',m
write(*,*)
write(*,*)

if ( .false. ) then
l = 0
write(*,*)'OVL matrix non zero elements:'
write(*,*)' ovl_mat, i, j'
do i=1, nb
   do k=1, nb
      if (ovl_mat(i,k) .NE. 0.0d0) then 
          write(*,*) ovl_mat(i,k),i,k
          l = l+1
      end if
   end do
end do

write(*,*) 'Total non-zero elements of ovl matrix: ', l
end if

end program



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine D1_gprod(alp1,beta1,A,B,govl)
implicit none
real(8)             :: p, pc(3), qc(3), q, q2(3)
real(8),intent(in)  :: alp1, beta1, A(3), B(3)
real(8),intent(out) :: govl(3)
real(8), parameter  :: pi=4*atan(1.0d0)


p = alp1+beta1 ! gaussian coefficient
pc(:) = (alp1*A(:) + beta1*B(:) )/p
q = alp1*beta1/p
qc(:) = A(:)-B(:)
q2(1) = (A(1)-B(1))**2.0d0
q2(2) = (A(2)-B(2))**2.0d0
q2(3) = (A(3)-B(3))**2.0d0
govl(:) = sqrt(pi/p)*exp(-q*q2(:))

return 
end subroutine



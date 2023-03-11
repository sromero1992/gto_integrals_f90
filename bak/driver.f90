program basis_test
  use class_basis_function
  use integral
  use com_vars
  implicit none
  type(gaussian)             :: g_tmp, ga, gb, gc, gd     ! Declare gaussian basis 
  type(integ)                :: S, T, Ven, Vee, Veex! Declaration of integral variables that uses g objects
  integer,allocatable        :: idx(:,:,:), jdx(:,:,:), N_SPDF(:,:), N_SPDF_SUM(:,:)
  integer                    :: SPDF(4), SPDF_SUM(4), idx_tmp, l_max, &
                                li, lj, lli, llj, iSPDF, jSPDF, iloc, jloc, & 
                                isa, isb, isc, isd, li2, lj2, lli2, llj2, iSPDF2, jSPDF2 , iloc2, jloc2 
  real(8)                    :: test_val
  real(8), allocatable       :: G_MAT(:,:,:,:) 


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



  write(*,*) 'NBAS : ',nbas
  all_atms = sum(n_atm(:))
  allocate( idx(all_atms, l_max, SPDF(l_max) ) )
  allocate( jdx(all_atms, l_max, SPDF(l_max) ) )
  !Check this definition for more classes such as H2O
  allocate( N_SPDF(4,all_atms))
  allocate( N_SPDF_SUM(4,all_atms))
 
  idx = 0
  jdx = 0
  idx_tmp = 0 
  do isa = 1, all_atms
     N_SPDF(:,isa) = (/ g(isa)%sf, g(isa)%pf, g(isa)%df, g(isa)%ff /) 
     N_SPDF_SUM(:,isa) = (/ 0, N_SPDF(1,isa),  sum(N_SPDF(1:2,isa)), sum(N_SPDF(1:3,isa)), sum(N_SPDF(1:4,isa))/)
     do li = 1, l_max
        do lli = 1, SPDF(li) 
            idx(isa, li, lli) = idx_tmp 
            jdx(isa, li, lli) = idx_tmp
            idx_tmp = idx_tmp + N_SPDF(li,isa) 
        end do
     end do
  end do 
 
 
  allocate(   S_MAT( nbas, nbas) )
  allocate(   T_MAT( nbas, nbas) )
  allocate( Ven_MAT( nbas, nbas) )
  S_MAT = 0.0d0
  T_MAT = 0.0d0
  Ven_MAT = 0.0d0
  do isa = 1, all_atms
     do isb = isa, all_atms
        do li = 1, l_max
           do lj =1, l_max
              do iSPDF = 1, N_SPDF(li,isa) 
                 do jSPDF = 1, N_SPDF(lj,isb) 
                    do lli = 1, SPDF(li)  
                       do llj = 1, SPDF(lj)  
                          iloc = idx(isa,li,lli)  + iSPDF
                          jloc = jdx(isb,lj,llj)  + jSPDF
                          g(isa)%shell = shells(li, lli, :) 
                          g_tmp = g(isa)
                          g(isb)%shell = shells(lj, llj, :)
                          call S%S_int( g_tmp, g(isb), iSPDF + N_SPDF_SUM(li,isa), jSPDF + N_SPDF_SUM(lj,isb) )      
                          call T%T_int( g_tmp, g(isb), iSPDF + N_SPDF_SUM(li,isa), jSPDF + N_SPDF_SUM(lj,isb) )      
                          call Ven%Ven_int( g_tmp, g(isb), iSPDF + N_SPDF_SUM(li,isa), jSPDF + N_SPDF_SUM(lj,isb), g(isa)%origin)      
                          S_MAT( iloc, jloc)    = S_MAT( iloc, jloc )   + S%int_val
                          T_MAT( iloc, jloc)    = T_MAT( iloc, jloc )   + T%int_val
                          Ven_MAT( iloc, jloc)  = Ven_MAT( iloc, jloc ) + Ven%int_val
                       end do
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do


  !Overlap matrix
  open(13,file='ovl2.txt',status='unknown')
  do ib = 1,nbas
      write(13,*) S_MAT(ib,:)
  end do
  close(13)
  !Kinetic energy matrix
  open(13,file='kin2.txt',status='unknown')
  do ib = 1,nbas
      write(13,*) T_MAT(ib,:)
  end do
  close(13)
  !Ven matrix
  open(13,file='coul_en.txt',status='unknown')
  do ib = 1,nbas
      write(13,*) Ven_MAT(ib,:)
  end do
  close(13)

 
  
  write(*,*)
  write(*,*) 'NBAS( NBAS + 1 )/2 :', nbas * ( nbas + 1 )/2
  allocate( S_VEC( nbas * ( nbas + 1 )/2 ) )
  allocate( T_VEC( nbas * ( nbas + 1 )/2 ) )
  allocate( Ven_VEC( nbas * ( nbas + 1 )/2 ) )
  l=0
  m=0
  do i=1, nbas
    do k=i, nbas
       l = l + 1
       S_VEC(l)   = S_MAT(i,k)
       T_VEC(l)   = T_MAT(i,k)
       Ven_VEC(l) = Ven_MAT(i,k)
       if (S_VEC(l) .NE. 0.0d0) then
           m = m+1
           !write(*,*) S_VEC(l)
       end if
    end do
  end do
  write(*,*) 'Number of elements in S_VEC : ', l  
  write(*,*) 'Number of non-zero elements in S_VEC : ', m 
   
  
  
  !Reading OVLBABY 1D array
  open(14,file='OVLBABY',form='UNFORMATTED')
  read(14) nrec
  write(*,*) 'nrec file :',nrec
  allocate(ovlbaby(nrec))
  read(14) (ovlbaby(j),j=1,nrec)
  close(14)
  
  
  !Check against ovlbaby which elements do we have
  open(15,file='ovl_ovlbaby2.txt') 
  l=0
  m=0
   do j = 1, nrec   
      if (ovlbaby(j) .NE. 0.0d0 ) then
         l = l+1
         do i = 1, nbas*(nbas + 1)/2 
            if ( (ovlbaby(j) .NE. 0.0d0) .AND. ( S_VEC(i) .NE. 0.0d0) ) then
               if (  abs( ovlbaby(j) - S_VEC(i) ) .LE.  del  ) then
  
                  write(15,*) ovlbaby(j), j, S_VEC(i), i ,'     X'
                  !write(*,*) ovlbaby(j), j, S_VEC(i), i
                  ovlbaby(j) = -1000.0d0
                  S_VEC(i)   = 1000.0d0
                  m = m+1
                  exit 
               else if ( (ovlbaby(j) .NE. -1000.0d0) .AND. (S_VEC(i) .NE. 1000.0d0)) then 
                  write(15,*) ovlbaby(j), j, S_VEC(i), i 
               end if
            end if
         end do
      end if
  end do               
  close(15)
  write(*,*) 'OVLBABY total elements :',nrec
  write(*,*) 'OVLBABY non-zero elements :', l
  write(*,*) 'S_VEC-OVLBABY similar non-zero elements elements :',m


!###################################### Density construction #######################################

  call dens


!########################################## COULOMB Vee INTEGRAL ###################################################
  allocate( COUL_MAT( nbas, nbas) )
  allocate( G_MAT( nbas, nbas, nbas, nbas) )
  COUL_MAT = 0.0d0
  G_MAT    = 0.0d0

  do isa = 1, all_atms
     ga = g(isa)
     do isb = 1, all_atms
        gb = g(isb)
        do li = 1, l_max
           do lj =1, l_max
              do iSPDF = 1, N_SPDF(li,isa) 
                 do jSPDF = 1, N_SPDF(lj,isb) 
                    do lli = 1, SPDF(li)  
                       do llj = 1, SPDF(lj)  
                          iloc = idx(isa,li,lli)  + iSPDF
                          jloc = jdx(isb,lj,llj)  + jSPDF
                          ga%shell = shells(li, lli, :) 
                          gb%shell = shells(lj, llj, :)
                          do isc = 1, all_atms
                             gc = g(isc)
                             do isd = 1, all_atms
                                gd = g(isd)
                                do li2 = 1, l_max
                                   do lj2 =1, l_max
                                      do iSPDF2 = 1, N_SPDF(li2,isc) 
                                         do jSPDF2 = 1, N_SPDF(lj2,isd) 
                                            !write(*,*)'isa,isb,isc,isd :',isa,isb,isc,isd
                                            !write(*,*)'iSPDF,jSPDF,iSPDF2,jSPDF2:',iSPDF,jSPDF,iSPDF2,jSPDF2
                                            do lli2 = 1, SPDF(li2)  
                                               do llj2 = 1, SPDF(lj2)  
                                                  iloc2 = idx(isc,li2,lli2)  + iSPDF2
                                                  jloc2 = jdx(isd,lj2,llj2)  + jSPDF2
                                                  gc%shell = shells(li2, lli2, :) 
                                                  gd%shell = shells(lj2, llj2, :)

                                                  call Vee%Vee_int( ga, gb, gc, gd, &
                                                     iSPDF + N_SPDF_SUM(li,isa) , jSPDF + N_SPDF_SUM(lj,isb), &
                                                     iSPDF2+ N_SPDF_SUM(li2,isc), jSPDF2+ N_SPDF_SUM(lj2,isd) )

                                                  !call Veex%Vee_int( ga, gc, gb, gd, &
                                                  !iSPDF + N_SPDF_SUM(li,isa), iSPDF2 + N_SPDF_SUM(lj2,isc), &
                                                  !jSPDF + N_SPDF_SUM(li,isb), jSPDF2 + N_SPDF_SUM(lj2,isd) )
                                                  
                                                  !COUL_MAT(iloc,jloc) = COUL_MAT(iloc,jloc) + D_MAT(iloc,jloc)* &
                                                  !                     ( Vee%int_val - Veex%int_val )
                                                  G_MAT(iloc,jloc,iloc2,jloc2) = Vee%int_val !Create G_mat to save doble calculations! 
                                                  !coul=vee-veex
                                               end do
                                            end do
                                         end do
                                      end do
                                   end do
                                end do
                             end do
                          end do
                       end do
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do

  do isa = 1, all_atms
     do isb = 1, all_atms
        do li = 1, l_max
           do lj =1, l_max
              do iSPDF = 1, N_SPDF(li,isa) 
                 do jSPDF = 1, N_SPDF(lj,isb) 
                    do lli = 1, SPDF(li)  
                       do llj = 1, SPDF(lj)  
                          iloc = idx(isa,li,lli)  + iSPDF
                          jloc = jdx(isb,lj,llj)  + jSPDF
                          do isc = 1, all_atms
                             do isd = 1, all_atms
                                do li2 = 1, l_max
                                   do lj2 =1, l_max
                                      do iSPDF2 = 1, N_SPDF(li2,isc) 
                                         do jSPDF2 = 1, N_SPDF(lj2,isd) 
                                            !write(*,*)'isa,isb,isc,isd :',isa,isb,isc,isd
                                            !write(*,*)'iSPDF,jSPDF,iSPDF2,jSPDF2:',iSPDF,jSPDF,iSPDF2,jSPDF2
                                            do lli2 = 1, SPDF(li2)  
                                               do llj2 = 1, SPDF(lj2)  
                                                  iloc2 = idx(isc,li2,lli2)  + iSPDF2
                                                  jloc2 = jdx(isd,lj2,llj2)  + jSPDF2

                                                  
                                                  COUL_MAT(iloc,jloc) = COUL_MAT(iloc,jloc) + D_MAT(iloc,jloc)* &
                                                        ( G_MAT(iloc,jloc,iloc2,jloc2) - G_MAT(iloc,iloc2,jloc,jloc2) )
                                                  !G_MAT(iloc,jloc,iloc2,jloc2) = Vee%int_val !Create G_mat to save doble calculations! 
                                                  !COUL_MAT(iloc,jloc) = COUL_MAT(iloc,jloc)+ &
                                                  !D(c,d)*( G_MAT(iloc, jloc, iloc2, jloc2) - G_MAT(iloc, iloc2, jloc, jloc2) )
                                               end do
                                            end do
                                         end do
                                      end do
                                   end do
                                end do
                             end do
                          end do
                       end do
                    end do
                 end do
              end do
           end do
        end do
     end do
  end do
  deallocate(G_MAT)

  !Ven matrix
  open(13,file='coul_ee.txt',status='unknown')
  do ib = 1,nbas
      write(13,*) COUL_MAT(ib,:)
  end do
  close(13)

  write(*,*)
  allocate( COUL_VEC( nbas * ( nbas + 1 )/2 ) )
  l=0
  m=0
  do i=1, nbas
    do k=i, nbas
       l = l + 1
       COUL_VEC(l) = COUL_MAT(i,k)
       if (COUL_VEC(l) .NE. 0.0d0) then
           m = m+1
       end if
    end do
  end do
  write(*,*) 'Number of elements in COUL_VEC : ', l  
  write(*,*) 'Number of non-zero elements in COUL_VEC : ', m 
   
  
  
 ! !Reading COULOMB 1D array
 ! open(14,file='COULOMB',form='UNFORMATTED')
 ! read(14) nrec
 ! write(*,*) 'nrec file :',nrec
 ! allocate(COULOMB(nrec))
 ! read(14) (COULOMB(j),j=1,nrec)
 ! close(14)
 ! 
 ! 
 ! !Check against ovlbaby which elements do we have
 ! open(15,file='coulomb_coul.txt') 
 ! l=0
 ! m=0
 !  do j = 1, nrec   
 !     if (COULOMB(j) .NE. 0.0d0 ) then
 !        l = l+1
 !        do i = 1, nbas*(nbas + 1)/2 
 !           if ( (COULOMB(j) .NE. 0.0d0) .AND. ( COUL_VEC(i) .NE. 0.0d0) ) then
 !              if (  abs( COULOMB(j) - COUL_VEC(i) ) .LE.  del  ) then
 ! 
 !                 write(15,*) COULOMB(j), j, COUL_VEC(i), i ,'     X'
 !                 !write(*,*) ovlbaby(j), j, S_VEC(i), i
 !                 ovlbaby(j) = -1000.0d0
 !                 S_VEC(i)   = 1000.0d0
 !                 m = m+1
 !                 exit 
 !              else if ( (COULOMB(j) .NE. -1000.0d0) .AND. (COUL_VEC(i) .NE. 1000.0d0)) then 
 !                 write(15,*) COULOMB(j), j, COUL_VEC(i), i 
 !              end if
 !           end if
 !        end do
 !     end if
 ! end do               
 ! close(15)
 ! write(*,*) 'COULOMB total elements :',nrec
 ! write(*,*) 'COULOMB non-zero elements :', l
 ! write(*,*) 'COUL_VEC-COULOMB similar non-zero elements elements :',m





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
  

  ! D-shell Dxxx => (3,0,0)
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


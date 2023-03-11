subroutine okven_matbuild
  use class_basis_function
  use integral
  use com_vars
  use mat_build
  implicit none
  type(gaussian)             :: ga, gb    ! Declare gaussian basis 
  type(integ)                :: S, T, Ven! Declaration of integral variables that uses g objects
  logical                    :: verbose
  integer                    :: ii, jj 
  verbose = .FALSE.

  if (f_bol) then 
     l_max = 4
     nbas = sum(g(:)%sf + 3*g(:)%pf + 6*g(:)%df + 10*g(:)%ff )
  else 
     l_max = 3
     nbas = sum(g(:)%sf + 3*g(:)%pf + 6*g(:)%df )
  end if



  write(*,*) 'NBAS : ',nbas
  all_atms = sum(site_type(:))
  allocate( idx(all_atms, l_max, SPDF(l_max) ) )
  !Check this definition for more classes such as H2O
  allocate( N_SPDF(4,all_atms))
  allocate( N_SPDF_SUM(4,all_atms))
 
  idx     = 0
  idx_tmp = 0
  do isite = 1, nsite 
     do isa = 1, site_type(isite) 
        N_SPDF(:,isa)     = (/ g(isa)%sf, g(isa)%pf, g(isa)%df, g(isa)%ff /) 
        N_SPDF_SUM(:,isa) = (/ 0, N_SPDF(1,isa),  sum(N_SPDF(1:2,isa)), sum(N_SPDF(1:3,isa)), sum(N_SPDF(1:4,isa))/)
        do li = 1, l_max
           do lli = 1, N_SPDF(li,isa) 
               idx(isa, li, lli) = idx_tmp 
               idx_tmp = idx_tmp + SPDF(li) 
           end do
        end do
     end do 
  end do 
 
 
  allocate(   S_MAT( nbas, nbas) )
  allocate(   T_MAT( nbas, nbas) )
  allocate( Ven_MAT( nbas, nbas) )
  S_MAT   = 0.0d0
  T_MAT   = 0.0d0
  Ven_MAT = 0.0d0

  do isite = 1, nsite
     do isa = 1, site_type(isite)
        
        ga = g(isa)
        do jsite = 1, nsite
           do isb = 1, site_type(jsite)

              gb = g(isb) 
              do li = 1, l_max

                 ii = N_SPDF_SUM(li,isa)
                 do lj =1, l_max

                    jj = N_SPDF_SUM(lj,isb)
                    do lli = 1, N_SPDF(li,isa)
                       do llj = 1, N_SPDF(lj,isb) 
                          do iSPDF = 1, SPDF(li) 

                             iloc = idx(isa,li,lli)  + iSPDF
                             ga%shell = shells(li, iSPDF, :) 
                             do jSPDF = 1, SPDF(lj)

                                jloc = idx(isb,lj,llj)  + jSPDF
                                gb%shell = shells(lj, jSPDF, :)

                                call S%S_int( ga, gb, ii + lli , jj + llj )      
                                call T%T_int( ga, gb, ii + lli , jj + llj )      
                                S_MAT( iloc, jloc)    = S_MAT( iloc, jloc )   + S%int_val
                                T_MAT( iloc, jloc)    = T_MAT( iloc, jloc )   + T%int_val
                                do i = 1, all_atms
                                   call Ven%Ven_int( ga, gb, ii + lli , jj + llj ,g(i)%origin)      
                                   Ven_MAT( iloc, jloc)  = Ven_MAT( iloc, jloc ) + 1.0d0*g(isa)%qn*Ven%int_val
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

  !Overlap matrix
  write(*,*) 
  write(*,*) 'S MATRIX CREATED :'
  open(13,file='ovl2.txt',status='unknown')
  do ib = 1,nbas
      write(13,*) S_MAT(ib,:)
      if (verbose) write(*,*) S_MAT(ib,:)
  end do
  close(13)
  !Kinetic energy matrix
  write(*,*) 
  write(*,*) 'T MATRIX CREATED:'
  open(13,file='kin2.txt',status='unknown')
  do ib = 1,nbas
      write(13,*) T_MAT(ib,:)
      if (verbose) write(*,*) T_MAT(ib,:)
  end do
  close(13)
  !Ven matrix
  write(*,*) 
  write(*,*) 'Ven MATRIX CREATED :'
  open(13,file='coul_en.txt',status='unknown')
  do ib = 1,nbas
      write(13,*) Ven_MAT(ib,:)
      if (verbose) write(*,*) Ven_MAT(ib,:)
  end do
  close(13)
  write(*,*) 
 
  
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
       end if
    end do
  end do
  write(*,*) 'Number of elements in S_VEC : ', l  
  write(*,*) 'Number of non-zero elements in S_VEC : ', m 
   
  
  
  if (chk_nrlmol) then 
     !Reading OVLBABY 1D array
     write(*,*) 'READING OVLBABY '
     open(14,file='OVLBABY',form='UNFORMATTED')
     read(14) nrec
     write(*,*) 'nrec file :',nrec
     allocate(ovlbaby(nrec))
     read(14) (ovlbaby(j),j=1,nrec)
     close(14)
      
     !Reading OVLBABY 1D array
     write(*,*) 'READING HAMBABY '
     open(14,file='HAMBABY',form='UNFORMATTED')
     read(14) nrec
     write(*,*) 'nrec file :',nrec
     allocate(hambaby(nrec))
     read(14) (hambaby(j),j=1,nrec)
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
     write(*,*) 'S_VEC---OVLBABY similar non-zero elements elements :',m
     write(*,*) 'MAXVAL OVL-SVEC :', MAXVAL( abs(S_VEC-ovlbaby) )  

     !Check against hambaby which elements do we have
     open(15,file='T_hambaby2.txt') 
     l=0
     m=0
      do j = 1, nrec   
         if (hambaby(j) .NE. 0.0d0 ) then
            l = l+1
            do i = 1, nbas*(nbas + 1)/2 
               if ( (hambaby(j) .NE. 0.0d0) .AND. ( T_VEC(i) .NE. 0.0d0) ) then
                  if (  abs( hambaby(j) - T_VEC(i) ) .LE.  del  ) then
     
                     !write(*,*) hambaby(j), j, T_VEC(i), i
                     write(15,*) hambaby(j), j, T_VEC(i), i ,'     X'
                     hambaby(j) = -1000.0d0
                     T_VEC(i)   = 1000.0d0
                     m = m+1
                     exit 
                  else if ( (hambaby(j) .NE. -1000.0d0) .AND. (T_VEC(i) .NE. 1000.0d0)) then 
                     write(15,*) hambaby(j), j, T_VEC(i), i 
                  end if
               end if
            end do
         end if
     end do               
     close(15)
     write(*,*) 'HAMBABY total elements :',nrec
     write(*,*) 'HAMBABY non-zero elements :', l
     write(*,*) 'T_VEC---HAMBABY similar non-zero elements elements :',m
     write(*,*) 'MAXVAL KIN-TVEC :', MAXVAL( abs(T_VEC-hambaby) )  
 
     deallocate(ovlbaby)
     deallocate(hambaby)
  end if
end subroutine


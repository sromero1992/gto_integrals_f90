subroutine okven_matbuild
  use class_basis_function
  use integral
  use module_com
  use module_g
  use mat_build
  implicit none
  type(gaussian)             :: ga, gb    ! Declare gaussian basis 
  type(integ)                :: S, T, Ven! Declaration of integral variables that uses g objects
  logical                    :: verbose
  real(8)                    :: atm_chg 
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
  !do isite = 1, nsite 
     !do isa = 1, site_type(isite) 
     do isa = 1, all_atms 
        N_SPDF(:,isa)     = (/ g(isa)%sf, g(isa)%pf, g(isa)%df, g(isa)%ff /) 
        N_SPDF_SUM(:,isa) = (/ 0, N_SPDF(1,isa),  sum(N_SPDF(1:2,isa)), sum(N_SPDF(1:3,isa)), sum(N_SPDF(1:4,isa))/)
        !This is identifying the spdf multiplicity later
        do li = 1, l_max
           do ifunc = 1, N_SPDF(li,isa) 
               !index depending on site, angular quantum number and number of functions of that angular 
               idx(isa, li, ifunc) = idx_tmp
               idx_tmp = idx_tmp + SPDF(li) 
           end do
        end do
     end do 
  !end do 
 
 
  allocate(   S_MAT( nbas, nbas) )
  allocate(   T_MAT( nbas, nbas) )
  allocate( Ven_MAT( nbas, nbas) )
  S_MAT   = 0.0d0
  T_MAT   = 0.0d0
  Ven_MAT = 0.0d0
  !do isite = 1, nsite !type sites i
     !do isa = 1, site_type(isite) !site a
     do isa = 1, all_atms !site a
        
        ga = g(isa)
        atm_chg = 1.0d0*ga%qn
        !do jsite = 1, nsite !types site j
           !do isb = 1, site_type(jsite) !site b
           !do isb = isa, site_type(jsite) !site b
           do isb = 1, all_atms !site b
           

              gb = g(isb) 
              do li = 1, l_max !angular momentum i

                 do lj = 1, l_max ! angular momentum j
                 !do lj =li, l_max ! angular momentum j

                    allocate( S%MAT( SPDF(li), SPDF(lj)) ) 
                    allocate( T%MAT( SPDF(li), SPDF(lj)) ) 
                    allocate( Ven%MAT( SPDF(li), SPDF(lj)) ) 

                    do ifunc = 1, N_SPDF(li,isa) ! number of bare SPDF functions

                       iloc = idx(isa,li,ifunc) 
                       do jfunc = 1, N_SPDF(lj,isb)
                       !do jfunc = ifunc, N_SPDF(lj,isb)

                          jloc = idx(isb,lj,jfunc)
                          !----------------------------------------------------------------------------------------------
                          !Those are integrating bare gaussians for provided function, and do it for all angular momentum 
                          !----------------------------------------------------------------------------------------------
                          call S%S_int( ga, gb, li, lj, ifunc , jfunc)      
                          S_MAT( iloc+1:iloc+SPDF(li), jloc+1:jloc+SPDF(lj) ) = S%MAT(1:SPDF(li),1:SPDF(lj))
                          !S_MAT( jloc+1:jloc+SPDF(lj), iloc+1:iloc+SPDF(li) ) = transpose( S%MAT(1:SPDF(li),1:SPDF(lj)) )

                          call T%T_int( ga, gb, li, lj, ifunc , jfunc)      
                          T_MAT( iloc+1:iloc+SPDF(li), jloc+1:jloc+SPDF(lj) ) = T%MAT(1:SPDF(li),1:SPDF(lj))
                          !T_MAT( jloc+1:jloc+SPDF(lj), iloc+1:iloc+SPDF(li) ) = transpose( T%MAT(1:SPDF(li),1:SPDF(lj)) )

                          !do i = 1, all_atms
                             call Ven%Ven_int( ga, gb, li, lj, ifunc , jfunc, ga%origin, atm_chg)      
                             Ven_MAT( iloc+1:iloc+SPDF(li), jloc+1:jloc+SPDF(lj) ) = Ven%MAT(1:SPDF(li),1:SPDF(lj)) !+ &
                                                                               !Ven_MAT( iloc+1:iloc+SPDF(li), jloc+1:jloc+SPDF(lj) )                                                 
                          !end do
                          !Ven_MAT( jloc+1:jloc+SPDF(lj), iloc+1:iloc+SPDF(li) ) = transpose(Ven_MAT( jloc+1:jloc+SPDF(lj),iloc+1:iloc+SPDF(li) )) 
                       end do
                    end do
                    deallocate( S%MAT, T%MAT, Ven%MAT)
                 end do
              end do
           end do
        !end do
     end do
  !end do

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
     do j = 1, nrec   
        write(15,*) j, S_vec(j),ovlbaby(j), S_Vec(j)-ovlbaby(j)
     end do    
     close(15)
     write(*,*) 'OVLBABY total elements :',nrec
     write(*,*) 'OVLBABY non-zero elements :', l
     write(*,*) 'S_VEC---OVLBABY similar non-zero elements elements :',m
     write(*,*) 'MAXVAL OVL-SVEC :', MAXVAL( abs(S_VEC-ovlbaby) )  

     !Check against hambaby which elements do we have
     open(15,file='T_hambaby2.txt') 
     do j = 1, nrec   
        write(15,*) j, T_vec(j),hambaby(j), T_Vec(j)-hambaby(j)
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


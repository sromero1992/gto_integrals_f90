program basis_test
  use class_basis_function
  use integral
  implicit none
  type(gaussian),allocatable :: g(:),g_tmp ! Declare gaussian basis 
  type(integ)                :: S, T, Ven, Vee ! Declare gaussian basis 
  real(8)                    :: x0(3,2), temp_exps(3), & !, boys1, R, RC(3), nuc_attraction, kinetic, overlap, &
                                el(2),  tmp1, tmp2, del
  integer                    :: n2, nb, nbas(3), tot_natm_type, atms_of_type_n,all_atms, i, j, k, kk, l, m, mm, ia, ib, bi, ai, & !jb, &
                                n_contracted_g, sf_tmp(2), pf_tmp(2), df_tmp(2), g_index, g_upper_bound, jm, im, nrec, &
                                nbasl(6),ima,jmb
  integer, allocatable       :: qe(:), qn(:), n_atm(:), sf(:), pf(:), df(:), sf_sup(:), pf_sup(:), df_sup(:), &
                                nb_tmp(:), nbas_tmp(:)
  real(8), allocatable       :: S_MAT(:,:), T_MAT(:,:), ovlbaby(:), S_VEC(:)
  character(10)              :: dum_str
  logical                    :: ext_bas, le


  del = 1.0d-7
  nb = 3
  tot_natm_type = 1
  atms_of_type_n = 2
  allocate(g(atms_of_type_n))
  !x0(xi,no.atom) 
  x0(1,1)=0.2d0
  x0(2,1)=0.0d0
  x0(3,1)=0.0d0
  !!!!!!!!!!!!!!!!
  !x0(:,2)=x0(:,1)
  x0(1,2)=0.0d0
  x0(2,2)=0.7d0
  x0(3,2)=0.4d0
  !!!!!!!!!!!!!



  g(1)%origin = x0(:,1)
  g(2)%origin = x0(:,2)
  g(1)%nb   = 3
  g(2)%nb   = 3
  g(1)%nbas = 1
  g(2)%nbas = 1

  call g(1)%allocation()
  call g(2)%allocation()

!  temp_exps(1)= 3.42525091d0
!  temp_exps(2)= 0.62391373d0
!  temp_exps(3)= 0.16885540d0
!  g(1)%exps    = temp_exps
!  g(2)%exps  = temp_exps
!
!  g(1)%coef(1) = 0.15432897d0
!  g(1)%coef(2) = 0.53532814d0
!  g(1)%coef(3) = 0.44463454d0
!  g(2)%coef  = g(1)%coef 
!
!  g(1)%shell(1)= 0
!  g(1)%shell(2)= 0
!  g(1)%shell(3)= 0
!  g(2)%shell = g(1)%shell
!
!  call g(1)%normalization()
!  call g(2)%normalization()
!  call g(1)%print_gaussian_info()  ! Call a class subroutine
!  call g(2)%print_gaussian_info()  ! Call a class subroutine
!
!  call S%S_int(g(1), g(2))
!  write(*,*) ' S(a,b) is : ', S%int_val 
!  call T%T_int(g(1), g(2))
!  write(*,*) ' T(a,b) is : ', T%int_val 
!  write(*,*) 'Boys1(1,0.5)', boys1(1,0.5d0)
!  write(*,*) 'R_{tuv}^n(p*RPC^2)', R(1,0,0,1,0.5d0,x0(:,1),1.0d0)
!  write(*,*) 'kinetic :', kinetic(g(1)%exps(1), g(1)%shell, g(1)%origin, g(2)%exps(1), g(2)%shell, g(2)%origin)  
!  write(*,*) 'overlap :', overlap(g(1)%exps(1), g(1)%shell, g(1)%origin, g(2)%exps(1), g(2)%shell, g(2)%origin)  
!
!  RC(1) =0.5d0
!  RC(2) =0.0d0
!  RC(3) =0.3d0
!  write(*,*) ' nuc_attraction(a, lmn1, RA, b, lmn2, RB, RC)', &
!             nuc_attraction(g(1)%exps(1), g(1)%shell, x0(:,1), g(2)%exps(1), g(2)%shell, x0(:,2), RC)
!  call Ven%Ven_int(g(1), g(2),RC)
!  write(*,*) ' Ven(a,b) is : ', Ven%int_val 
!  call Vee%Vee_int(g(1)%exps(1), g(1)%shell, g(1)%origin, g(1)%exps(1), g(1)%shell, g(1)%origin, &
!                   g(2)%exps(1), g(2)%shell, g(2)%origin, g(2)%exps(1), g(2)%shell, g(2)%origin)
!  write(*,*) ' Vee(a,b,a,b) is : ', Vee%int_val 
 
  deallocate(g)

  
  ext_bas = .FALSE.
  open(11,file='ISYMGEN')
  read(11,*) tot_natm_type  

  allocate(   n_atm(tot_natm_type))
  allocate(nbas_tmp(tot_natm_type))
  allocate(  nb_tmp(tot_natm_type))
  !Reading mainly total number of atoms in calculation to construct the objects
  do i = 1, tot_natm_type
     read(11,*) tmp1, tmp2 !qn is the atom e.g. qn=1 is hydrogen 
     read(11,*) dum_str
     read(11,*) n_atm(i) !number of atoms of type i
     do j = 1, n_atm(i) + 1
         read(11,*) dum_str
     end do
     read(11,*) nb_tmp(i)!Number of bare gaussians
     read(11,*) sf_tmp(1),pf_tmp(1),df_tmp(1) 
     read(11,*) sf_tmp(2),pf_tmp(2),df_tmp(2) 
     !The following turns off suplemental S P D functions
     if ( ext_bas ) then
        nbas_tmp(i) =  sum(sf_tmp(:)) + sum(pf_tmp(:)) + sum(df_tmp(:)) 
     else
        nbas_tmp(i) =  sf_tmp(1) + pf_tmp(1) + df_tmp(1) 
     end if
     write(*,*) 'number of S P D functions (includes supplemental): ',nbas_tmp(i)
     do j = 1, nbas_tmp(i) + 1 !+1 for the gaussian exponent
        do k = 1, nb_tmp(i) / 3
           read(11,*) dum_str 
        end do
     end do
     if (i .EQ. tot_natm_type) then
        read(11,*) dum_str
        if ( dum_str  .EQ. 'ELECTRONS') then
            write(*,*) 'Successful allocation of needed arrays for objects'
        end if 
     end if
  end do  
  write(*,*) 'functions S and supplemental: ',sf_tmp(1) ,',',sf_tmp(2)
  write(*,*) 'functions P and supplemental: ',pf_tmp(1) ,',',pf_tmp(2)
  write(*,*) 'functions D and supplemental: ',df_tmp(1) ,',',df_tmp(2)

  n_contracted_g = sum(n_atm(:))
  write(*,*) 'total number of contracted gaussians: ', n_contracted_g

  ! Allocation for contracted gaussians
  allocate(      g (n_contracted_g)) !Number of contracted gaussian basis 
  ! Allocation for ith atom type
  allocate(      qe (tot_natm_type)) !Number of electrons
  allocate(      qn (tot_natm_type)) !Nuclei number, for each gaussian
  !S P D functions and supplemental
  allocate(      sf (tot_natm_type)) 
  allocate(      pf (tot_natm_type)) 
  allocate(      df (tot_natm_type)) 
  allocate(  sf_sup (tot_natm_type)) 
  allocate(  pf_sup (tot_natm_type)) 
  allocate(  df_sup (tot_natm_type)) 
  write(*,*) 'allocation completed'

  rewind(unit = 11)
  read(11,*) tot_natm_type  
  do i = 1, tot_natm_type 
     read(11,*) qe(i), qn(i) !qn is the atom e.g. qn=1 is hydrogen 
     read(11,*) dum_str
     read(11,*) n_atm(i) !number of atoms of type i
     do j = 1, n_atm(i) + 1
         read(11,*) dum_str
     end do
     !Think ontranslation from g(1:n_atm(i)) and next type g(sum(1:n_atm(i-1))+1:n_atm(i))
     if ( i .EQ. 1 ) then
        g_index = 1  !first gaussian
     else 
        g_index = sum( n_atm( 1: (i-1) ) ) + 1
     end if
     g(g_index)%n_atm = qn(i) ! Atomic number
     read(11,*) g(g_index)%nb!Number of bare gaussians
     write(*,*) 'total number of bare gaussian (',1,') : ', g(g_index)%nb
     !Read S P D functions
     read(11,*) g(g_index)%sf    , g(g_index)%pf    , g(g_index)%df
     read(11,*) g(g_index)%sf_sup, g(g_index)%pf_sup, g(g_index)%df_sup 
     if (ext_bas) then
        !nbas in gaussian structure is the number of basis functions
        g(g_index)%nbas =  g(g_index)%sf + g(g_index)%pf + g(g_index)%df + &
                           g(g_index)%sf_sup + g(g_index)%pf_sup +  g(g_index)%df_sup 
     else
        g(g_index)%nbas = g(g_index)%sf + g(g_index)%pf + g(g_index)%df  
     end if 
     write(*,*) 'total number of basis functions(',g_index,') : ', g(g_index)%nbas
     call g(g_index)%allocation() !Allocates nb in exponents and coefficients
     write(*,*) 'size of gaussian basis coefficients C(ib,ibas) :', size(g(g_index)%coef)
     !Reading all alpha exponents
     do j = 1, g(g_index)%nb/3
        read(11,*) g(g_index)%exps(3*j-2), g(g_index)%exps(3*j-1), g(g_index)%exps(3*j)
     end do
     if (ext_bas) then
        write(*,*) 'NBAS : ',g(g_index)%nbas
        do ib = 1, g(g_index)%nbas !nbas is the number of basis functions in g(i) basis
           do j = 1, g(g_index)%nb/3 !nb is the number of bare gaussians
              read(11,*) g(g_index)%coef(3*j-2, ib),  g(g_index)%coef(3*j-1, ib),  g(g_index)%coef(3*j, ib) 
           end do
        end do
     else
        !Reading all gaussian contraction coefficients Cij
        write(*,*) 'Number of Gaussian basis functions: ',g(g_index)%nbas ! This is contains the number of basis functions
        do ib = 1, g(g_index)%nbas
           !This ifs are for skipping supplemental basis functions
           if      ( ( ib .EQ. g(g_index)%sf+1 ) .AND. (g(g_index)%sf_sup .NE. 0 ) ) then
               n2 = g(g_index)%sf_sup
           else if ( ( ib .EQ. g(g_index)%sf+g(g_index)%pf+1 ) .AND. (g(g_index)%pf_sup .NE. 0)) then
               n2 = g(g_index)%pf_sup
           else if ( ( ib .GT. g(g_index)%sf+g(g_index)%pf+g(g_index)%df+1 ) .AND. (g(g_index)%df_sup .NE. 0) ) then
               n2 = g(g_index)%df_sup
           else 
               n2 = 0
           end if
           do k = 1, n2
              do j = 1, g(g_index)%nb/3
                 read(11,*) temp_exps(1),  temp_exps(1),  temp_exps(1) 
              end do
           end do
           do j = 1, g(g_index)%nb/3
              read(11,*) g(g_index)%coef(3*j-2,ib),  g(g_index)%coef(3*j-1,ib),  g(g_index)%coef(3*j,ib) 
           end do
           !write(*,*) g(g_index)%coef(:,ib)
        end do
     end if 
!Pass all the read to all the contracted gaussians of same type!
     if ( (n_atm(i) .GT. 1) ) then
        if (i .GT. 1) then
           g_upper_bound = sum(n_atm(1:i))
        else
           g_upper_bound = n_atm(1)
        end if
        do k = g_index + 1, g_upper_bound 
           g(k) = g(g_index)
        end do
     end if
  end do
  close(11)

  !Read tom positions from  SYMBOL 
  open(12,file='SYMBOL')
  do i=1,10
     read(12,*)
  end do
  do i = 1, tot_natm_type
     do j = 1, n_atm(i)
        read(12,*)  dum_str, dum_str, g(j)%origin(1), g(j)%origin(2), g(j)%origin(3), dum_str
     end do
  end do
  read(12,*) dum_str, dum_str, el(1), el(2) ! number of electrons up and dn
  close(12)

  deallocate(  nb_tmp)
  deallocate(nbas_tmp)
  deallocate(      sf) 
  deallocate(      pf) 
  deallocate(      df) 
  deallocate(  sf_sup) 
  deallocate(  pf_sup) 
  deallocate(  df_sup) 
 



  !For my ISYMGEN and SYMBOL, I have to same hydrogens
  g(1)%shell(1) = 0
  g(1)%shell(2) = 0
  g(1)%shell(3) = 0
  g(2)%shell = g(1)%shell
  
  !call g(1)%print_gaussian_info() 
  !call g(2)%print_gaussian_info() 
  !call S%S_int( g(1), g(2))
  !write(*,*) ' S(a,b) is : ', S%int_val 

  !construct S_MAT !using nbas here made problems ...
  if (ext_bas) then
      nbas = sum ( g(:)%sf ) !+ g(:)%sf_sup+ g(:)%pf*3 + g(:)%pf_sup*3 + g(:)%df*6 +g(:)%df_sup*6 )
  else
      !nbas = sum ( g(:)%sf ) !+ g(:)%pf*3 + g(:)%df*6)
      nbas(1) = sum ( g(:)%sf )
      nbas(2) = sum ( g(:)%sf + g(:)%pf*3)
      nbas(3) = sum ( g(:)%sf + g(:)%pf*3 + g(:)%df*6)
  end if
  write(*,*) 'NBAS last:', nbas(3)
  allocate( S_MAT(nbas(3),nbas(3)))
  allocate( T_MAT(nbas(3),nbas(3)))
  all_atms = sum(n_atm(:))
  l=1 !Has to be over all the different types (l) of atoms
  !write(*,*) 'Overlap matrix'
  S_MAT = 0.0d0
  do ai = 1, all_atms   ! atom sites
     do bi = ai, all_atms
        nbasl(1) = g(ai)%sf
        nbasl(2) = g(ai)%sf + g(ai)%pf*3
        nbasl(3) = g(ai)%sf + g(ai)%pf*3 + g(ai)%df*6
        nbasl(4) = g(bi)%sf
        nbasl(5) = g(bi)%sf + g(bi)%pf*3
        nbasl(6) = g(bi)%sf + g(bi)%pf*3 + g(bi)%df*6
        do i = 1 + (ai-1)*nbasl(3), ai* nbasl(3) !nbas !Loop over number of basis functions
           do j = max(1 + (bi-1)*nbasl(6), i), bi*nbasl(6) !nbas
!Beging of S-X overlaps
              g(ai)%shell = 0
              g(bi)%shell = 0
              ! ima and jmb are variables to do a correct blocking
              ima = mod(i, nbasl(3) )
              if (ima .EQ. 0) ima = nbasl(3)
              jmb = mod(j, nbasl(6) )
              if (jmb .EQ. 0) jmb = nbasl(6)
              im = mod(i , nbasl(3)) !This is used for indexing of coefficients in gaussian basis
              if (im .EQ. 0) im = nbasl(3)

              if ( ima .LE. g(ai)%sf  ) then  
                 g_tmp = g(ai) !Temp gaussian basis set to keep a different angular momentum
!Calculate S-S overlap
                 if  ( jmb .LE. g(bi)%sf ) then !bi are the columns of S_MAT
                     jm = mod(j , g(bi)%sf)
                     if ( j .GT. nbasl(6) ) jm = mod(j - nbasl(6) , g(bi)%sf)
                     if (jm .EQ. 0) jm = g(bi)%sf
                     call S%S_int( g_tmp, g(bi), im, jm)      
                     call T%T_int( g_tmp, g(bi), im, jm)      
                     !S_MAT(i,j) = S_MAT(i,j) + S%int_val
                 end if
!Loop approach for  S-Pi overlap
                 do k = 1, 3
                    if ( ( jmb .LE. g(bi)%sf + g(bi)%pf*k ) .AND. (jmb .GT. g(bi)%sf + g(bi)%pf*(k-1) )  )  then  !Maybe a loop to compact this
                       g(bi)%shell = 0
                       g(bi)%shell(k) = 1 !Shell to activate, e.g. (1,0,0) => Px to calc overlap S-Px
                       jm = mod(j - g(bi)%sf , g(bi)%pf)
                       if ( j .GT. nbasl(6) ) jm = mod(j - g(bi)%sf - nbasl(6) , g(bi)%pf)
                       if (jm .EQ. 0) jm = g(bi)%pf 
                       call S%S_int( g_tmp, g(bi), im, g(bi)%sf + jm) !g_tmp is g(ai), when ai=bi this is used    
                       call T%T_int( g_tmp, g(bi), im, g(bi)%sf + jm)     
                       !S_MAT(i,j) = S_MAT(i,j) + S%int_val
                       exit
                    end if
                 end do
!Loop approach for S-Di overlap    ************** working now!!
                 do kk = 1 , 6
                    if ( ( jmb .LE. g(bi)%sf + g(bi)%pf*3+ g(bi)%df*kk ) .AND. &
                         ( jmb .GT. g(bi)%sf + g(bi)%pf*3 + g(bi)%df*(kk-1) )  )  then  !Maybe a loop to compact this
                       !Shell combinations: (2,0,0) => dxx ;  (0,0,2) => dzz  ;   (1,0,1) => dxz
                       !                    (0,2,0) => dyy ;  (1,1,0) => dxy  ;   (0,1,1) => dyz
                       g(bi)%shell = 0
                       if ( kk .LE. 3) then
                          g(bi)%shell(kk) = 2 ! kk=1  =>  dxx (2,0,0) ; kk=2  => dyy =(0,2,0) ; kk=3  =>  dzz =(0,0,3)
                       else if ( kk .GT. 3) then
                          if ( kk .EQ. 4 ) then  ! (1,1,0) dxy
                             g(bi)%shell(1) = 1 
                             g(bi)%shell(2) = 1
                          else if ( kk .EQ. 5 ) then  ! (1,0,1) dxz
                             g(bi)%shell(1) = 1 
                             g(bi)%shell(3) = 1
                          else                   ! (0,1,1) dyz
                             g(bi)%shell(2) = 1 
                             g(bi)%shell(3) = 1
                          end if
                       end if 
                       jm = mod(j - g(bi)%sf -3*g(bi)%pf, g(bi)%df)
                       if ( j .GT. nbasl(6) ) jm = mod(j - g(bi)%sf -3*g(bi)%pf - nbasl(6) , g(bi)%df)
                       if (jm .EQ. 0) jm = g(bi)%df
                       call S%S_int( g_tmp, g(bi), im, g(bi)%sf + g(bi)%pf + jm)      
                       call T%T_int( g_tmp, g(bi), im, g(bi)%sf + g(bi)%pf + jm)      
                       !S_MAT(i,j) = S_MAT(i,j) + S%int_val
                       exit !This exit is done for avoiding more j evaluations
                    end if
                 end do
              end if 
!    End S-X overlaps

!Begin of overlap P-X -> X=S,Pi,Di
              do m = 1, 3 
                 if ( (ima .LE. g(ai)%sf + g(ai)%pf*m) .AND. (ima .GT. g(ai)%sf + g(ai)%pf*(m-1)  )) then 
                    g(ai)%shell = 0
                    im = mod(i - g(ai)%sf , g(ai)%pf)
                    if ( i .GT. nbasl(3) ) im = mod(i - g(ai)%sf - nbasl(3) , g(ai)%pf)
                    if (im .EQ. 0) im = g(ai)%pf 
                    g(ai)%shell(m) = 1
                    g_tmp = g(ai)  !Temp gaussian basis set to keep angular momentum
                    g(bi)%shell = 0
!Overlap Pi-S
                    if  ( jmb .LE. g(bi)%sf ) then !bi are the columns of S_MAT
                       jm = mod(j , g(bi)%sf)
                       if ( j .GT. nbasl(6) ) jm = mod(j - nbasl(6) , g(bi)%sf)
                       if (jm .EQ. 0) jm = g(bi)%sf
                       call S%S_int( g_tmp, g(bi), g(ai)%sf +im, jm)      
                       call T%T_int( g_tmp, g(bi), g(ai)%sf +im, jm)      
                       !S_MAT(i,j) = S_MAT(i,j) + S%int_val
                       exit
                    end if

!Loop approach for Pi-Pj overlap
                    do k = 1, 3  
                       if ( ( jmb .LE. g(bi)%sf + g(bi)%pf*k ) .AND. (jmb .GT. g(bi)%sf + g(bi)%pf*(k-1) )  )  then  !Maybe a loop to compact this
                          g(bi)%shell = 0
                          g(bi)%shell(k) = 1 !Shell to activate, e.g. (1,0,0) => Px to calc overlap S-Px
                          jm = mod(j - g(bi)%sf , g(bi)%pf)
                          if ( j .GT. nbasl(6) ) jm = mod(j - g(bi)%sf - nbasl(6) , g(bi)%pf)
                          if (jm .EQ. 0) jm = g(bi)%pf 
                          call S%S_int( g_tmp, g(bi), g(ai)%sf + im, g(bi)%sf + jm)      
                          call T%T_int( g_tmp, g(bi), g(ai)%sf + im, g(bi)%sf + jm)      
                          !S_MAT(i,j) = S_MAT(i,j) + S%int_val
                          exit
                       end if
                    end do
!Loop approach for Pi-Dj overlap      ****This section is working!!!
                    do kk = 1 , 6
                       if ( ( jmb .LE. g(bi)%sf + g(bi)%pf*3+ g(bi)%df*kk ) .AND. &
                            (jmb .GT. g(bi)%sf + g(bi)%pf*3 + g(bi)%df*(kk-1) )  )  then  !Maybe a loop to compact this
                          !Shell combinations: (2,0,0) => dxx ;  (0,0,2) => dzz  ;   (1,0,1) => dxz
                          !                    (0,2,0) => dyy ;  (1,1,0) => dxy  ;   (0,1,1) => dyz
                          g(bi)%shell = 0
                          if ( kk .LE. 3) then
                             g(bi)%shell(kk) = 2 ! kk=1  =>  dxx (2,0,0) ; kk=2  => dyy =(0,2,0) ; kk=3  =>  dzz =(0,0,3)
                          else if ( kk .GT. 3) then
                             if ( kk .EQ. 4 ) then  ! (1,1,0) dxy
                                g(bi)%shell(1) = 1 
                                g(bi)%shell(2) = 1
                             else if ( kk .EQ. 5 ) then  ! (1,0,1) dxz
                                g(bi)%shell(1) = 1 
                                g(bi)%shell(3) = 1
                             else                   ! (0,1,1) dyz
                                g(bi)%shell(2) = 1 
                                g(bi)%shell(3) = 1
                             end if
                          end if 
                          jm = mod(j - g(bi)%sf -3*g(bi)%pf, g(bi)%df)
                          if ( j .GT. nbasl(6) ) jm = mod(j - g(bi)%sf -3*g(bi)%pf - nbasl(6) , g(bi)%df)
                          if (jm .EQ. 0) jm = g(bi)%df
                          call S%S_int( g_tmp, g(bi), g(ai)%sf + im, g(bi)%sf + g(bi)%pf + jm)      
                          call T%T_int( g_tmp, g(bi), g(ai)%sf + im, g(bi)%sf + g(bi)%pf + jm)      
                          !S_MAT(i,j) = S_MAT(i,j) + S%int_val
                        !  write(*,*) '----------------------- Pi-Dj -----------------------------'
                        !  write(*,*) '**********Shell ai ,  bi: ',g_tmp%shell,g(bi)%shell
                        !  write(*,*) '***** ITERATION VALS  m, k, im, jm: ', m, k, im, jm
                        !  write(*,*) '****Ingredients of integral', i, j, im+g(ai)%sf, jm+ g(bi)%sf+g(bi)%pf, S%int_val
                          exit !This exit is done for avoiding more j evaluations
                       end if
                    end do
                 end if 
              end do
!   End Pi-X overlaps

!Begin of overlap Di-X  *************************This section should be re-worked
              do k = 1, 6 
                  if ( ( ima .LE. g(ai)%sf + g(ai)%pf*3 + g(ai)%df*k ) .AND. &
                       ( ima .GT. g(ai)%sf + g(ai)%pf*3 + g(ai)%df*(k-1) )  )  then
                     !Shell combinations: (2,0,0) => dxx ;  (0,0,2) => dzz  ;   (1,0,1) => dxz
                     !                    (0,2,0) => dyy ;  (1,1,0) => dxy  ;   (0,1,1) => dyz
                     g(ai)%shell = 0
                     g(bi)%shell = 0
                     if ( k .LE. 3) then
                        g(ai)%shell(k) = 2 ! kk=1  =>  dxx (2,0,0) ; kk=2  => dyy =(0,2,0) ; kk=3  =>  dzz =(0,0,3)
                     else if ( k .GT. 3) then
                        if ( k .EQ. 4 ) then  ! (1,1,0) dxy
                           g(ai)%shell(1) = 1 
                           g(ai)%shell(2) = 1
                        else if ( k .EQ. 5 ) then  ! (1,0,1) dxz
                           g(ai)%shell(1) = 1 
                           g(ai)%shell(3) = 1
                        else                   ! (0,1,1) dyz
                           g(ai)%shell(2) = 1 
                           g(ai)%shell(3) = 1
                        end if
                     end if
                     g_tmp = g(ai)
                     g(bi)%shell = 0
                     im = mod(i- g(ai)%sf -3*g(ai)%pf, g(ai)%df)
                     if ( j .GT. nbasl(6) ) jm = mod(j - g(bi)%sf -3*g(bi)%pf - nbasl(6) , g(bi)%df)
                     if (im .EQ. 0) im = g(ai)%df
!!Overlap for Di-S
                     if  ( jmb .LE. g(bi)%sf ) then !bi are the columns of S_MAT
                        jm = mod(j , g(bi)%sf)
                        if ( j .GT. nbasl(6) ) jm = mod(j - nbasl(6) , g(bi)%sf)
                        if (jm .EQ. 0) jm = g(bi)%sf
                        call S%S_int( g_tmp, g(bi), g(ai)%sf +g(ai)%pf + im, jm)      
                        call T%T_int( g_tmp, g(bi), g(ai)%sf +g(ai)%pf + im, jm)      
                        !S_MAT(i,j) = S_MAT(i,j) + S%int_val
                        !write(*,*) '------------------ Di-S -------------------------------'
                        !write(*,*) '****Ingredients of integral', i, j, im+g(ai)%sf+g(ai)%pf, jm+ g(bi)%sf, S%int_val
                        exit
                     end if
!!Loop approach for Di-Pj
                     do kk = 1, 3  
                        if ( ( jmb .LE. g(bi)%sf + g(bi)%pf*kk ) .AND. (jmb .GT. g(bi)%sf + g(bi)%pf*(kk-1) )  )  then  !Maybe a loop to compact this
                           g(bi)%shell = 0
                           g(bi)%shell(kk) = 1 !Shell to activate, e.g. (1,0,0) => Px to calc overlap S-Px
                           jm = mod(j - g(bi)%sf , g(bi)%pf)
                           if ( j .GT. nbasl(6) ) jm = mod(j - g(bi)%sf - nbasl(6) , g(bi)%pf)
                           if (jm .EQ. 0) jm = g(bi)%pf 
                           call S%S_int( g_tmp, g(bi), g(ai)%sf + g(ai)%pf + im, g(bi)%sf + jm)      
                           call T%T_int( g_tmp, g(bi), g(ai)%sf + g(ai)%pf + im, g(bi)%sf + jm)      
                           !S_MAT(i,j) = S_MAT(i,j) + S%int_val
                           !write(*,*) '------------------ Di-Pj ----------------------------'
                           !write(*,*) '**********Shell ai ,  bi: ',g_tmp%shell,g(bi)%shell
                           !write(*,*) '***** ITERATION VALS  m, k, im, jm: ', m, k, im, jm
                           !write(*,*) '****Ingredients of integral', i, j, im+g(ai)%sf+g(ai)%pf, jm+ g(bi)%sf, S%int_val
                           exit
                        end if
                     end do
!Loop approach for Di-Dj overlap
                     do kk = 1 , 6
                        if ( ( jmb .LE. g(bi)%sf + g(bi)%pf*3+ g(bi)%df*kk ) .AND. &
                             ( jmb .GT. g(bi)%sf + g(bi)%pf*3 + g(bi)%df*(kk-1) )  )  then  !Maybe a loop to compact this
                           !Shell combinations: (2,0,0) => dxx ;  (0,0,2) => dzz  ;   (1,0,1) => dxz
                           !                    (0,2,0) => dyy ;  (1,1,0) => dxy  ;   (0,1,1) => dyz
                           g(bi)%shell = 0
                           if ( kk .LE. 3) then
                              g(bi)%shell(kk) = 2 ! kk=1  =>  dxx (2,0,0) ; kk=2  => dyy =(0,2,0) ; kk=3  =>  dzz =(0,0,3)
                           else if ( kk .GT. 3) then
                              if ( kk .EQ. 4 ) then  ! (1,1,0) dxy
                                 g(bi)%shell(1) = 1 
                                 g(bi)%shell(2) = 1
                              else if ( kk .EQ. 5 ) then  ! (1,0,1) dxz
                                 g(bi)%shell(1) = 1 
                                 g(bi)%shell(3) = 1
                              else                   ! (0,1,1) dyz
                                 g(bi)%shell(2) = 1 
                                 g(bi)%shell(3) = 1
                              end if
                           end if 
                           jm = mod(j - g(bi)%sf -3*g(bi)%pf, g(bi)%df)
                           if ( j .GT. nbasl(6) ) jm = mod(j - g(bi)%sf -3*g(bi)%pf - nbasl(6) , g(bi)%df)
                           if (jm .EQ. 0) jm = g(bi)%df
                           call S%S_int( g_tmp, g(bi), g(ai)%sf +g(ai)%pf + im, g(bi)%sf + g(bi)%pf + jm)      
                           call T%T_int( g_tmp, g(bi), g(ai)%sf +g(ai)%pf + im, g(bi)%sf + g(bi)%pf + jm)      
                           !S_MAT(i,j) = S_MAT(i,j) + S%int_val
                           !write(*,*) '-------------------- Di-Dj -------------------------------'
                           !write(*,*) '**********Shell ai ,  bi: ',g_tmp%shell,g(bi)%shell
                           !write(*,*) '****Ingredients of integral', i, j, im+g(ai)%sf+g(ai)%pf, jm+ g(bi)%sf, S%int_val
                           exit !This exit is done for avoiding more j evaluations
                        end if
                     end do
                  end if 
              end do
              S_MAT(i,j) = S_MAT(i,j) + S%int_val
              T_MAT(i,j) = T_MAT(i,j) + T%int_val
           end do
        end do
     end do
  end do
  
 write(*,*)
 inquire(file='ovl2.txt',exist=le)
 if (le) then
   open(13,file='ovl2.txt')
   close(13,status='delete')
 end if
 i=1
!OVERLAP matrix
 open(13,file='ovl2.txt',status='unknown')
 do ib = 1,nbas(3)
     write(13,*) S_MAT(ib,:)
 end do
 close(13)
!Kinetic energy matrix
open(13,file='kin2.txt',status='unknown')
 do ib = 1,nbas(3)
     write(13,*) T_MAT(ib,:)
 end do
 close(13)



 if ( .FALSE. ) then
    l = 0
    write(*,*)' ovl_mat, i, j'
    do i=1, nbas(3)
       write(*,*) S_MAT(i,:)
       write(*,*)
    end do
 end if
 

 write(*,*) 'nbas(nbas+1)/2 :',nbas(3)*(nbas(3)+1)/2
 allocate(S_VEC(nbas(3)*(nbas(3)+1)/2))
 l=0
 m=0
 write(*,*) 'Non- zero elements of S_VEC:' 
 do i=1, nbas(3)
   do k=i, nbas(3)
      l = l + 1
      S_VEC(l) = S_MAT(i,k)
      if (S_VEC(l) .NE. 0.0d0) then
          m = m+1
          !write(*,*) S_VEC(l)
      end if
   end do
 end do
write(*,*) 'Number of elements in S_VEC : ', l  
write(*,*) 'Number of non-zero elements in S_VEC : ', m 
 
 !Reading OVLBABY to vector
 open(14,file='OVLBABY',form='UNFORMATTED')
 read(14) nrec
 write(*,*) 'nrec file :',nrec
 allocate(ovlbaby(nrec))
 read(14) (ovlbaby(j),j=1,nrec)
 close(14)

! call quicksort(ovlbaby, 1, nrec)
! call quicksort(S_VEC, 1, nrec)

 open(15,file='ovl_ovlbaby.txt')
  do j = 1, nrec
     if ( ovlbaby(j) .NE. 0.0d0) then
        if (   abs( ovlbaby(j) - S_VEC(j) ) .LE.  del ) then
           write(15,*) ovlbaby(j), S_VEC(j), j,'   X'
        else 
           write(15,*) ovlbaby(j), S_VEC(j), j
        end if
     end if
 end do
 close(15)

 !Check against ovlbaby which elements do we have
 open(15,file='ovl_ovlbaby2.txt') 
 l=0
 m=0
  do j = 1, nrec   
     if (ovlbaby(j) .NE. 0.0d0 ) then
        l = l+1
        do i = 1, nbas(3)*(nbas(3) + 1)/2 
           if ( (ovlbaby(j) .NE. 0.0d0) .AND. ( S_VEC(i) .NE. 0.0d0) ) then
              if (  abs( ovlbaby(j) - S_VEC(i) ) .LE.  del  ) then

                 write(15,*) ovlbaby(j), j, S_VEC(i), i ,'     X'
                 write(*,*) ovlbaby(j), j, S_VEC(i), i
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
 write(*,*)
 write(*,*)





!do im = 1, 6
!   g(1)%shell = 0
!   if ( im .LE. 3) then
!      g(1)%shell(im) = 2 ! kk=1  =>  dxx (2,0,0) ; kk=2  => dyy =(0,2,0) ; kk=3  =>  dzz =(0,0,3)
!   else if ( im .GT. 3) then
!      if ( im .EQ. 4 ) then  ! (1,1,0) dxy
!         g(1)%shell(1) = 1 
!         g(1)%shell(2) = 1
!      else if ( im .EQ. 5 ) then  ! (1,0,1) dxz
!         g(1)%shell(1) = 1 
!         g(1)%shell(3) = 1
!      else                   ! (0,1,1) dyz
!         g(1)%shell(2) = 1 
!         g(1)%shell(3) = 1
!      end if
!   end if 
!   do jm = 1, 6 
!      g(2)%shell = 0
!      if ( jm .LE. 3) then
!         g(2)%shell(jm) = 2 ! kk=1  =>  dxx (2,0,0) ; kk=2  => dyy =(0,2,0) ; kk=3  =>  dzz =(0,0,3)
!      else if ( jm .GT. 3) then
!         if ( jm .EQ. 4 ) then  ! (1,1,0) dxy
!            g(2)%shell(1) = 1 
!            g(2)%shell(2) = 1
!         else if ( jm .EQ. 5 ) then  ! (1,0,1) dxz
!            g(2)%shell(1) = 1 
!            g(2)%shell(3) = 1
!         else                   ! (0,1,1) dyz
!            g(2)%shell(2) = 1 
!            g(2)%shell(3) = 1
!         end if
!      end if 
!      do i = 8, 8
!         do j = 8, 8
!            call S%S_int( g(1), g(2),i,j)
!            write(*,*) '**********S_int ingridients (i,j,shell) : ', i, j, g(1)%shell, g(2)%shell,  S%int_val
!         end do
!      end do
!    end do
!end do

!########################################## COULOMB INTEGRALS ###################################################
!Electron-nuclei












end program 



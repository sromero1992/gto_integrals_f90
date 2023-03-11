program basis_test
  use class_basis_function
  use integral
  use com_vars
  implicit none
  type(gaussian)             :: g_tmp     ! Declare gaussian basis 
  type(integ)                :: S, T, Ven, Vee ! Declaration of integral variables that uses g objects
  
  del = 1.0d-7
  ext_bas = .FALSE.
  f_bol = .TRUE.

! in this subroutine you can re-define the order of shell integrals 
! such as dxx dyy dzz dxy dxz dyz of D-shell
  call shell_gen(s_shells, p_shells, d_shells, f_shells)

  call read_info

  !construct S_MAT !using nbas here made problems ...
  if (ext_bas) then
      if (f_bol) then
         nbas = sum ( g(:)%sf + g(:)%sf_sup+ 3*(g(:)%pf + g(:)%pf_sup) + & 
                 6*(g(:)%df +g(:)%df_sup) + 10*(g(:)%ff + g(:)%ff_sup) )
      else
         nbas = sum ( g(:)%sf + g(:)%sf_sup+ 3*(g(:)%pf + g(:)%pf_sup) + 6*(g(:)%df +g(:)%df_sup ) )
      end if
  else
      if (f_bol) then
         nbas = sum ( g(:)%sf  + 3*g(:)%pf + 6*g(:)%df + 10*g(:)%ff )
      else
         nbas = sum ( g(:)%sf  + 3*g(:)%pf + 6*g(:)%df )
      end if
  end if
  write(*,*) 'NBAS last:', nbas
  allocate(   S_MAT( nbas, nbas) )
  allocate(   T_MAT( nbas, nbas) )
  allocate( Ven_MAT( nbas, nbas) )
  all_atms = sum(n_atm(:))
  !write(*,*) 'Overlap matrix'
  S_MAT = 0.0d0
  T_MAT = 0.0d0
  Ven_MAT = 0.0d0
  do ai = 1, all_atms   ! atom sites
     do bi = ai, all_atms
        if ( f_bol) then
           nbasl(1) = g(ai)%sf + 3*g(ai)%pf + 6*g(ai)%df +10*g(ai)%df  !nbas for loop
           nbasl(2) = g(bi)%sf + 3*g(bi)%pf + 6*g(bi)%df +10*g(bi)%df 
        else 
           nbasl(1) = g(ai)%sf + 3*g(ai)%pf + 6*g(ai)%df !nbas for loop
           nbasl(2) = g(bi)%sf + 3*g(bi)%pf + 6*g(bi)%df  
        end if
        do i = 1 + (ai-1)*nbasl(1), ai* nbasl(1) !nbas !Loop over number of basis functions
           do j = max(1 + (bi-1)*nbasl(2), i), bi*nbasl(2) !nbas
!Beging of S-X overlaps
              g(ai)%shell = 0
              g(bi)%shell = 0
              ! ima and jmb are variables to do a correct blocking
              ima = mod(i, nbasl(1) )
              if (ima .EQ. 0) ima = nbasl(1)
              jmb = mod(j, nbasl(2) )
              if (jmb .EQ. 0) jmb = nbasl(2)
              im = mod(i , nbasl(1)) !This is used for indexing of coefficients in gaussian basis
              if (im .EQ. 0) im = nbasl(1)

              if ( ima .LE. g(ai)%sf  ) then  
                 g_tmp = g(ai) !Temp gaussian basis set to keep a different angular momentum
!Calculate S-S overlap
                 if  ( jmb .LE. g(bi)%sf ) then !bi are the columns of S_MAT
                     jm = mod(j , g(bi)%sf)
                     if ( j .GT. nbasl(2) ) jm = mod(j - nbasl(2) , g(bi)%sf)
                     if (jm .EQ. 0) jm = g(bi)%sf
                     call S%S_int( g_tmp, g(bi), im, jm)      
                     call T%T_int( g_tmp, g(bi), im, jm)      
                     call Ven%Ven_int( g_tmp, g(bi), im, jm, g(ai)%origin)     
                     !S_MAT(i,j) = S_MAT(i,j) + S%int_val
                 end if
!Loop approach for  S-Pi overlap
                 do k = 1, 3
                    if ( ( jmb .LE. g(bi)%sf + g(bi)%pf*k ) .AND. (jmb .GT. g(bi)%sf + g(bi)%pf*(k-1) )  )  then  !Maybe a loop to compact this
                       g(bi)%shell = p_shells(k,:)
                       jm = mod(j - g(bi)%sf , g(bi)%pf)
                       if ( j .GT. nbasl(2) ) jm = mod(j - g(bi)%sf - nbasl(2) , g(bi)%pf)
                       if (jm .EQ. 0) jm = g(bi)%pf 
                       call S%S_int( g_tmp, g(bi), im, g(bi)%sf + jm) !g_tmp is g(ai), when ai=bi this is used    
                       call T%T_int( g_tmp, g(bi), im, g(bi)%sf + jm)     
                       call Ven%Ven_int( g_tmp, g(bi), im, g(bi)%sf + jm, g(ai)%origin)      
                       !S_MAT(i,j) = S_MAT(i,j) + S%int_val
                       exit
                    end if
                 end do
!Loop approach for S-Di overlap   
                 do kk = 1 , 6
                    if ( ( jmb .LE. g(bi)%sf + 3*g(bi)%pf + g(bi)%df*kk ) .AND. &
                         ( jmb .GT. g(bi)%sf + 3*g(bi)%pf + g(bi)%df*(kk-1) )  )  then  !Maybe a loop to compact this
                       g(bi)%shell = d_shells(kk,:)
                       jm = mod(j - g(bi)%sf -3*g(bi)%pf, g(bi)%df)
                       if ( j .GT. nbasl(2) ) jm = mod(j - g(bi)%sf -3*g(bi)%pf - nbasl(2) , g(bi)%df)
                       if (jm .EQ. 0) jm = g(bi)%df
                       call S%S_int( g_tmp, g(bi), im, g(bi)%sf + g(bi)%pf + jm)      
                       call T%T_int( g_tmp, g(bi), im, g(bi)%sf + g(bi)%pf + jm)      
                       call Ven%Ven_int( g_tmp, g(bi), im, g(bi)%sf + g(bi)%pf + jm, g(ai)%origin)      
                       !S_MAT(i,j) = S_MAT(i,j) + S%int_val
                       exit !This exit is done for avoiding more j evaluations
                    end if
                 end do
!Loop approach for S-Fi overlap  
                 if ( f_bol ) then 
                    do kk = 1 , 10
                       if ( ( jmb .LE. g(bi)%sf + 3*g(bi)%pf + 6*g(bi)%df + g(bi)%ff*kk ) .AND. &
                            ( jmb .GT. g(bi)%sf + 3*g(bi)%pf + 6*g(bi)%df + g(bi)%ff*(kk-1) )  )  then  !Maybe a loop to compact this
                          g(bi)%shell = f_shells(kk,:)
                          jm = mod(j - g(bi)%sf -3*g(bi)%pf -6*g(bi)%df, g(bi)%ff)
                          if ( j .GT. nbasl(2) ) jm = mod(j - g(bi)%sf -3*g(bi)%pf - 6*g(bi)%df - nbasl(2) , g(bi)%ff)
                          if (jm .EQ. 0) jm = g(bi)%ff
                          call S%S_int( g_tmp, g(bi), im, g(bi)%sf + g(bi)%pf + g(bi)%df + jm)      
                          call T%T_int( g_tmp, g(bi), im, g(bi)%sf + g(bi)%pf + g(bi)%df + jm)      
                          call Ven%Ven_int( g_tmp, g(bi), im, g(bi)%sf + g(bi)%pf + g(bi)%df + jm, g(ai)%origin)      
                          !S_MAT(i,j) = S_MAT(i,j) + S%int_val
                          exit !This exit is done for avoiding more j evaluations
                       end if
                    end do
                 end if
              end if 
!    End S-X overlaps

!Begin of overlap P-X -> X = S, Pi, Di
              do m = 1, 3 
                 if ( ( ima .LE. g(ai)%sf + g(ai)%pf*m ) .AND. &
                      ( ima .GT. g(ai)%sf + g(ai)%pf*(m-1)  ) ) then 
                    im = mod(i - g(ai)%sf , g(ai)%pf)
                    if ( i .GT. nbasl(1) ) im = mod(i - g(ai)%sf - nbasl(1) , g(ai)%pf)
                    if (im .EQ. 0) im = g(ai)%pf 
                    g(ai)%shell = p_shells(m,:)
                    g_tmp = g(ai)  !Temp gaussian basis set to keep angular momentum
!Overlap Pi-S
                    if  ( jmb .LE. g(bi)%sf ) then !bi are the columns of S_MAT
                       g(bi)%shell = s_shells
                       jm = mod(j , g(bi)%sf)
                       if ( j .GT. nbasl(2) ) jm = mod(j - nbasl(2) , g(bi)%sf)
                       if (jm .EQ. 0) jm = g(bi)%sf
                       call S%S_int( g_tmp, g(bi), g(ai)%sf +im, jm)      
                       call T%T_int( g_tmp, g(bi), g(ai)%sf +im, jm)      
                       call Ven%Ven_int( g_tmp, g(bi), g(ai)%sf +im, jm ,g(ai)%origin) 
                       !S_MAT(i,j) = S_MAT(i,j) + S%int_val
                       exit
                    end if

!Loop approach for Pi-Pj overlap
                    do k = 1, 3  
                       if ( ( jmb .LE. g(bi)%sf + g(bi)%pf*k ) .AND. & 
                            ( jmb .GT. g(bi)%sf + g(bi)%pf*(k-1) )  )  then  !Maybe a loop to compact this
                          g(bi)%shell = p_shells(k,:) 
                          jm = mod(j - g(bi)%sf , g(bi)%pf)
                          if ( j .GT. nbasl(2) ) jm = mod(j - g(bi)%sf - nbasl(2) , g(bi)%pf)
                          if (jm .EQ. 0) jm = g(bi)%pf 
                          call S%S_int( g_tmp, g(bi), g(ai)%sf + im, g(bi)%sf + jm)      
                          call T%T_int( g_tmp, g(bi), g(ai)%sf + im, g(bi)%sf + jm)      
                          call Ven%Ven_int( g_tmp, g(bi), g(ai)%sf + im, g(bi)%sf + jm, g(ai)%origin)      
                          !S_MAT(i,j) = S_MAT(i,j) + S%int_val
                          exit
                       end if
                    end do
!Loop approach for Pi-Dj overlap    
                    do kk = 1 , 6
                       if ( ( jmb .LE. g(bi)%sf + 3*g(bi)%pf + g(bi)%df*kk ) .AND. &
                            ( jmb .GT. g(bi)%sf + 3*g(bi)%pf + g(bi)%df*(kk-1) )  )  then  !Maybe a loop to compact this
                          g(bi)%shell = d_shells(kk,:)
                          jm = mod(j - g(bi)%sf -3*g(bi)%pf, g(bi)%df)
                          if ( j .GT. nbasl(2) ) jm = mod(j - g(bi)%sf -3*g(bi)%pf - nbasl(2) , g(bi)%df)
                          if (jm .EQ. 0) jm = g(bi)%df
                          call S%S_int( g_tmp, g(bi), g(ai)%sf + im, g(bi)%sf + g(bi)%pf + jm)      
                          call T%T_int( g_tmp, g(bi), g(ai)%sf + im, g(bi)%sf + g(bi)%pf + jm)      
                          call Ven%Ven_int( g_tmp, g(bi), g(ai)%sf + im, g(bi)%sf + g(bi)%pf + jm, g(ai)%origin)      
                          !S_MAT(i,j) = S_MAT(i,j) + S%int_val
                          exit !This exit is done for avoiding more j evaluations
                       end if
                    end do
!Loop approach for Pi-Fj overlap   
                    if ( f_bol ) then 
                       do kk = 1 , 10
                          if ( ( jmb .LE. g(bi)%sf + 3*g(bi)%pf + 6*g(bi)%df + g(bi)%ff*kk ) .AND. &
                               ( jmb .GT. g(bi)%sf + 3*g(bi)%pf + 6*g(bi)%df + g(bi)%ff*(kk-1) )  )  then  !Maybe a loop to compact this
                             g(bi)%shell = f_shells(kk,:)
                             jm = mod(j - g(bi)%sf - 3*g(bi)%pf - 6*g(bi)%df, g(bi)%ff)
                             if ( j .GT. nbasl(2) ) jm = mod(j - g(bi)%sf -3*g(bi)%pf - 6*g(bi)%df - nbasl(2) , g(bi)%ff)
                             if (jm .EQ. 0) jm = g(bi)%ff
                             call S%S_int( g_tmp, g(bi), g(ai)%sf + im, g(bi)%sf + g(bi)%pf + g(bi)%df + jm)      
                             call T%T_int( g_tmp, g(bi), g(ai)%sf + im, g(bi)%sf + g(bi)%pf + g(bi)%df + jm)      
                             call Ven%Ven_int( g_tmp, g(bi), g(ai)%sf + im, g(bi)%sf + g(bi)%pf + g(bi)%df + jm, g(ai)%origin)      
                             !S_MAT(i,j) = S_MAT(i,j) + S%int_val
                             exit !This exit is done for avoiding more j evaluations
                          end if
                       end do
                    end if
                 end if 
              end do
!   End Pi-X overlaps

!Begin of overlap Di-X  *************************This section should be re-worked
              do k = 1, 6 
                 if ( ( ima .LE. g(ai)%sf + 3*g(ai)%pf + g(ai)%df*k ) .AND. &
                      ( ima .GT. g(ai)%sf + 3*g(ai)%pf + g(ai)%df*(k-1) )  )  then
                     g(ai)%shell = d_shells(k,:)
                     g_tmp = g(ai)
                     im = mod(i- g(ai)%sf -3*g(ai)%pf, g(ai)%df)
                     if ( i .GT. nbasl(1) ) im = mod(i - g(ai)%sf -3*g(ai)%pf - nbasl(1) , g(ai)%df)
                     if (im .EQ. 0) im = g(ai)%df
!!Overlap for Di-S
                     if  ( jmb .LE. g(bi)%sf ) then !bi are the columns of S_MAT
                        g(bi)%shell = s_shells
                        jm = mod(j , g(bi)%sf)
                        if ( j .GT. nbasl(2) ) jm = mod(j - nbasl(2) , g(bi)%sf)
                        if (jm .EQ. 0) jm = g(bi)%sf
                        call S%S_int( g_tmp, g(bi), g(ai)%sf +g(ai)%pf + im, jm)      
                        call T%T_int( g_tmp, g(bi), g(ai)%sf +g(ai)%pf + im, jm)      
                        call Ven%Ven_int( g_tmp, g(bi), g(ai)%sf +g(ai)%pf + im, jm, g(ai)%origin)      
                        !S_MAT(i,j) = S_MAT(i,j) + S%int_val
                        exit
                     end if
!!Loop approach for Di-Pj
                     do kk = 1, 3  
                        if ( ( jmb .LE. g(bi)%sf + g(bi)%pf*kk ) .AND. &
                             ( jmb .GT. g(bi)%sf + g(bi)%pf*(kk-1) )  )  then  
                           g(bi)%shell = p_shells(kk,:)
                           jm = mod(j - g(bi)%sf , g(bi)%pf)
                           if ( j .GT. nbasl(2) ) jm = mod(j - g(bi)%sf - nbasl(2) , g(bi)%pf)
                           if (jm .EQ. 0) jm = g(bi)%pf 
                           call S%S_int( g_tmp, g(bi), g(ai)%sf + g(ai)%pf + im, g(bi)%sf + jm)      
                           call T%T_int( g_tmp, g(bi), g(ai)%sf + g(ai)%pf + im, g(bi)%sf + jm)      
                           call Ven%Ven_int( g_tmp, g(bi), g(ai)%sf + g(ai)%pf + im, g(bi)%sf + jm, g(ai)%origin)      
                           !S_MAT(i,j) = S_MAT(i,j) + S%int_val
                           exit
                        end if
                     end do
!Loop approach for Di-Dj overlap
                     do kk = 1 , 6
                        if ( ( jmb .LE. g(bi)%sf + 3*g(bi)%pf + g(bi)%df*kk ) .AND. &
                             ( jmb .GT. g(bi)%sf + 3*g(bi)%pf + g(bi)%df*(kk-1) )  )  then  
                           g(bi)%shell = d_shells(kk,:)
                           jm = mod(j - g(bi)%sf -3*g(bi)%pf, g(bi)%df)
                           if ( j .GT. nbasl(2) ) jm = mod(j - g(bi)%sf -3*g(bi)%pf - nbasl(2) , g(bi)%df)
                           if (jm .EQ. 0) jm = g(bi)%df
                           call S%S_int( g_tmp, g(bi), g(ai)%sf +g(ai)%pf + im, g(bi)%sf + g(bi)%pf + jm)      
                           call T%T_int( g_tmp, g(bi), g(ai)%sf +g(ai)%pf + im, g(bi)%sf + g(bi)%pf + jm)      
                           call Ven%Ven_int( g_tmp, g(bi), g(ai)%sf +g(ai)%pf + im, g(bi)%sf + g(bi)%pf + jm, g(ai)%origin)  
                           !S_MAT(i,j) = S_MAT(i,j) + S%int_val
                           !write(*,*) '-------------------- Di-Dj -------------------------------'
                           !write(*,*) '**********Shell ai ,  bi: ',g_tmp%shell,g(bi)%shell
                           !write(*,*) '****Ingredients of integral', i, j, im+g(ai)%sf+g(ai)%pf, jm+ g(bi)%sf, S%int_val
                           exit !This exit is done for avoiding more j evaluations
                        end if
                     end do
!Loop approach for Di-Fj overlap
                     if ( f_bol ) then 
                        do kk = 1 , 10
                           if ( ( jmb .LE. g(bi)%sf + 3*g(bi)%pf + 6*g(bi)%df + g(bi)%ff*kk ) .AND. &
                                ( jmb .GT. g(bi)%sf + 3*g(bi)%pf + 6*g(bi)%df + g(bi)%ff*(kk-1) )  )  then  !Maybe a loop to compact this
                              g(bi)%shell = f_shells(kk,:)
                              jm = mod(j - g(bi)%sf -3*g(bi)%pf -6*g(bi)%df, g(bi)%ff)
                              if ( j .GT. nbasl(2) ) jm = mod(j - g(bi)%sf -3*g(bi)%pf - 6*g(bi)%df - nbasl(2) , g(bi)%ff)
                              if (jm .EQ. 0) jm = g(bi)%ff
                              call S%S_int( g_tmp, g(bi), g(ai)%sf + g(ai)%pf + im, g(bi)%sf + g(bi)%pf + g(bi)%df + jm)      
                              call T%T_int( g_tmp, g(bi), g(ai)%sf + g(ai)%pf + im, g(bi)%sf + g(bi)%pf + g(bi)%df + jm)      
                              call Ven%Ven_int( g_tmp, g(bi), g(ai)%sf + g(ai)%pf + im, g(bi)%sf + g(bi)%pf + g(bi)%df + &
                                    jm, g(ai)%origin)      
                              !S_MAT(i,j) = S_MAT(i,j) + S%int_val
                              exit !This exit is done for avoiding more j evaluations
                           end if
                        end do
                     end if
                 end if 
              end do

!Begin of overlap Fi-X  *************************This section should be re-worked
              if ( f_bol) then
                 do k = 1, 10
                    if ( ( ima .LE. g(ai)%sf + 3*g(ai)%pf + 6*g(ai)%df + g(ai)%ff*k ) .AND. &
                         ( ima .GT. g(ai)%sf + 3*g(ai)%pf + 6*g(ai)%df + g(ai)%ff*(k-1) )  )  then
                        g(ai)%shell = f_shells(k,:)
                        g_tmp = g(ai)
                        im = mod(i- g(ai)%sf -3*g(ai)%pf -6*g(ai)%df, g(ai)%ff)
                        if ( i .GT. nbasl(1) ) im = mod(i - g(ai)%sf -3*g(ai)%pf - 6*g(ai)%df - nbasl(1) , g(ai)%ff)
                        if (im .EQ. 0) im = g(ai)%ff
 !!Overlap for Fi-S
                        if  ( jmb .LE. g(bi)%sf ) then !bi are the columns of S_MAT
                           g(bi)%shell = s_shells
                           jm = mod(j , g(bi)%sf)
                           if ( j .GT. nbasl(2) ) jm = mod(j - nbasl(2) , g(bi)%sf)
                           if (jm .EQ. 0) jm = g(bi)%sf
                           call S%S_int( g_tmp, g(bi), g(ai)%sf + g(ai)%pf + g(ai)%df + im, jm)      
                           call T%T_int( g_tmp, g(bi), g(ai)%sf + g(ai)%pf + g(ai)%df + im, jm)      
                           call Ven%Ven_int( g_tmp, g(bi), g(ai)%sf + g(ai)%pf + g(ai)%df + im, jm, g(ai)%origin)      
                           !S_MAT(i,j) = S_MAT(i,j) + S%int_val
                           exit
                        end if
!!Loop approach for Fi-Pj
                        do kk = 1, 3  
                           if ( ( jmb .LE. g(bi)%sf + g(bi)%pf*kk ) .AND. &
                                ( jmb .GT. g(bi)%sf + g(bi)%pf*(kk-1) )  )  then  
                              g(bi)%shell = p_shells(kk,:)
                              jm = mod(j - g(bi)%sf , g(bi)%pf)
                              if ( j .GT. nbasl(2) ) jm = mod(j - g(bi)%sf - nbasl(2) , g(bi)%pf)
                              if (jm .EQ. 0) jm = g(bi)%pf 
                              call S%S_int( g_tmp, g(bi), g(ai)%sf + g(ai)%pf + g(ai)%df + im, &
                                            g(bi)%sf + jm)      
                              call T%T_int( g_tmp, g(bi), g(ai)%sf + g(ai)%pf + g(ai)%df + im, &
                                            g(bi)%sf + jm)      
                              call Ven%Ven_int( g_tmp, g(bi), g(ai)%sf + g(ai)%pf + g(ai)%df + im, &
                                                g(bi)%sf + jm, g(ai)%origin)      
                              !S_MAT(i,j) = S_MAT(i,j) + S%int_val
                              exit
                           end if
                        end do
!Loop approach for Fi-Dj overlap
                        do kk = 1 , 6
                           if ( ( jmb .LE. g(bi)%sf + 3*g(bi)%pf + g(bi)%df*kk ) .AND. &
                                ( jmb .GT. g(bi)%sf + 3*g(bi)%pf + g(bi)%df*(kk-1) )  )  then  
                              g(bi)%shell = d_shells(kk,:)
                              jm = mod(j - g(bi)%sf -3*g(bi)%pf, g(bi)%df)
                              if ( j .GT. nbasl(2) ) jm = mod(j - g(bi)%sf -3*g(bi)%pf - nbasl(2) , g(bi)%df)
                              if (jm .EQ. 0) jm = g(bi)%df
                              call S%S_int( g_tmp, g(bi), g(ai)%sf + g(ai)%pf + g(ai)%df + im, &
                                            g(bi)%sf + g(bi)%pf + jm)      
                              call T%T_int( g_tmp, g(bi), g(ai)%sf + g(ai)%pf + g(ai)%df + im, &
                                            g(bi)%sf + g(bi)%pf + jm)      
                              call Ven%Ven_int( g_tmp, g(bi), g(ai)%sf + g(ai)%pf + g(ai)%df + im, &
                                                g(bi)%sf + g(bi)%pf + jm, g(ai)%origin)      
                              !S_MAT(i,j) = S_MAT(i,j) + S%int_val
                              !write(*,*) '-------------------- Di-Dj -------------------------------'
                              !write(*,*) '**********Shell ai ,  bi: ',g_tmp%shell,g(bi)%shell
                              !write(*,*) '****Ingredients of integral', i, j, im+g(ai)%sf+g(ai)%pf, jm+ g(bi)%sf, S%int_val
                              exit !This exit is done for avoiding more j evaluations
                           end if
                        end do
!Loop approach for Fi-Fj overlap
                        do kk = 1 , 10
                           if ( ( jmb .LE. g(bi)%sf + 3*g(bi)%pf + 6*g(bi)%df + g(bi)%ff*kk ) .AND. &
                                ( jmb .GT. g(bi)%sf + 3*g(bi)%pf + 6*g(bi)%df + g(bi)%ff*(kk-1) )  )  then  !Maybe a loop to compact this
                              g(bi)%shell = f_shells(kk,:)
                              jm = mod(j - g(bi)%sf -3*g(bi)%pf -6*g(bi)%df, g(bi)%ff)
                              if ( j .GT. nbasl(2) ) jm = mod(j - g(bi)%sf -3*g(bi)%pf - 6*g(bi)%df - nbasl(2) , g(bi)%ff)
                              if (jm .EQ. 0) jm = g(bi)%ff
                              call S%S_int( g_tmp, g(bi), g(ai)%sf + g(ai)%pf + g(ai)%df + im, &
                                            g(bi)%sf + g(bi)%pf + g(bi)%df + jm)      
                              call T%T_int( g_tmp, g(bi), g(ai)%sf + g(ai)%pf + g(ai)%df + im, &
                                            g(bi)%sf + g(bi)%pf + g(bi)%df + jm)      
                              call Ven%Ven_int( g_tmp, g(bi), g(ai)%sf + g(ai)%pf + g(ai)%df + im, &
                                                g(bi)%sf + g(bi)%pf + g(bi)%df + jm, g(ai)%origin)      
                             ! write(*,*) '-------------------- Fi-Fj -------------------------------'
                             ! write(*,*) '**********Shell ai ,  bi: ',g_tmp%shell,g(bi)%shell
                             ! write(*,*) '****Ingredients of integral ai, bi, i, j, im, jm,', ai, bi, i, j, &
                             !            im + g(ai)%sf + g(ai)%pf, jm + g(bi)%sf + g(bi)%pf , S%int_val
                              !S_MAT(i,j) = S_MAT(i,j) + S%int_val
                              exit !This exit is done for avoiding more j evaluations
                           end if
                        end do
                    end if
                 end do
              end if
              S_MAT(i,j) = S_MAT(i,j) + S%int_val
              T_MAT(i,j) = T_MAT(i,j) + T%int_val
              Ven_MAT(i,j) = Ven_MAT(i,j) + Ven%int_val*g(ai)%n_atm
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



 write(*,*) 'NBAS( NBAS + 1 )/2 :', nbas * ( nbas + 1 )/2
 allocate( S_VEC( nbas * ( nbas + 1 )/2 ) )
 l=0
 m=0
 !write(*,*) 'Non- zero elements of S_VEC:' 
 do i=1, nbas
   do k=i, nbas
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
 write(*,*)
 write(*,*)





!do im = 1, 10
!   !g(1)%shell = 0
!   g(1)%shell = f_shells(im,:)
!   do jm = 1, 10 
!      g(2)%shell = f_shells(jm,:)
!     do i = 9, 9
!         do j = 9, 9
!            call S%S_int( g(1), g(2),i,j)
!            write(*,*) '**********S_int ingridients (i,j,shell) : ', i, j, g(1)%shell, g(2)%shell,  S%int_val
!         end do
!      end do
!   end do
!end do

!########################################## COULOMB Vee INTEGRAL ###################################################
!Electron-nuclei
!write(*,*) 'nuclear attraction 1: ',nuc_attraction(g(1)%exps(1), g(1)%shell, g(1)%origin, g(1)%exps(2), &
!                                       g(1)%shell, g(1)%origin, g(1)%origin)







end program 




!#########################################################################
subroutine shell_gen(s_shells, p_shells, d_shells, f_shells)
   !implicit none
   integer,intent(inout)   :: s_shells(3), p_shells(3,3), d_shells(6,3), f_shells(10,3)
   integer                 :: ii
   character(10)           :: chr_tmp
   logical                 :: read_file

   read_file = .FALSE.

   s_shells(:) = 0

  ! P-shell Px => (1,0,0) 
   p_shells(1,1) = 1
   p_shells(1,2) = 0
   p_shells(1,3) = 0
  ! P-shell Py => (0,1,0)
   p_shells(2,1) = 0
   p_shells(2,2) = 1
   p_shells(2,3) = 0
  ! P-shell Py => (0,0,1)
   p_shells(3,1) = 0
   p_shells(3,2) = 0
   p_shells(3,3) = 1



  ! D-shell Dxx => (2,0,0)
   d_shells(1,1) = 2
   d_shells(1,2) = 0
   d_shells(1,3) = 0
  ! D-shell Dyy => (0,2,0)
   d_shells(2,1) = 0
   d_shells(2,2) = 2
   d_shells(2,3) = 0
  ! D-shell Dzz => (0,0,2)
   d_shells(3,1) = 0
   d_shells(3,2) = 0
   d_shells(3,3) = 2
  ! D-shell Dxy => (1,1,0)
   d_shells(4,1) = 1
   d_shells(4,2) = 1
   d_shells(4,3) = 0
  ! D-shell Dxz => (1,0,1)
   d_shells(5,1) = 1
   d_shells(5,2) = 0
   d_shells(5,3) = 1
  ! D-shell Dyz => (0,1,1)
   d_shells(6,1) = 0
   d_shells(6,2) = 1
   d_shells(6,3) = 1
  


  ! D-shell Dxxx => (3,0,0)
   f_shells(1,1) = 3
   f_shells(1,2) = 0
   f_shells(1,3) = 0
  ! D-shell Dyyy => (0,3,0)
   f_shells(2,1) = 0
   f_shells(2,2) = 3
   f_shells(2,3) = 0
  ! D-shell Dzzz => (0,0,3)
   f_shells(3,1) = 0
   f_shells(3,2) = 0
   f_shells(3,3) = 3
  ! D-shell Dxxy => (2,1,0)
   f_shells(4,1) = 2
   f_shells(4,2) = 1
   f_shells(4,3) = 0
  ! D-shell Dxxz => (2,0,1)
   f_shells(5,1) = 2
   f_shells(5,2) = 0
   f_shells(5,3) = 1
  ! D-shell Dxyy => (1,2,0)
   f_shells(6,1) = 1
   f_shells(6,2) = 2
   f_shells(6,3) = 0
  ! D-shell Dyyz => (0,2,1)
   f_shells(7,1) = 0
   f_shells(7,2) = 2
   f_shells(7,3) = 1
  ! D-shell Dxxz => (1,0,2)
   f_shells(8,1) = 1
   f_shells(8,2) = 0
   f_shells(8,3) = 2
  ! D-shell Dxyy => (0,1,2)
   f_shells(9,1) = 0
   f_shells(9,2) = 1
   f_shells(9,3) = 2
  ! D-shell Dxyy => (1,1,1)
   f_shells(10,1) = 1
   f_shells(10,2) = 1
   f_shells(10,3) = 1
   

   if (read_file) then
      open(21,file='shells.dat')
      read(21,*) chr_tmp
      do ii = 1, 3
         read(21,*) p_shells(ii,1), p_shells(ii,2), p_shells(ii,3) 
      end do
      read(21,*) chr_tmp
      do ii = 1, 6
         read(21,*) d_shells(ii,1), d_shells(ii,2), d_shells(ii,3) 
      end do
      read(21,*) chr_tmp
      do ii = 1, 10
         read(21,*) f_shells(ii,1), f_shells(ii,2), f_shells(ii,3) 
      end do
   end if
end subroutine


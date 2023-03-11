subroutine coul_build
   use class_basis_function
   use integral
   use mat_build
   use com_vars 
   type(gaussian)         :: ga, gb, gc, gd     ! Declare gaussian basis 
   type(integ)            :: Vee ! Declaration of integral variables that uses g objects
   real(8), allocatable   :: Gab(:,:)
   real(8)                :: coul_val, t0, t1, t2, t3, t4 !, tau2 


!########################################## COULOMB Vee INTEGRAL ###################################################
   call cpu_time(t0)
   tau = 10.0d-10 
   allocate( G_MAT( nbas, nbas, nbas, nbas) )
   allocate(   Gab(             nbas, nbas) )
   Gab(      :,:) = 0.0d0
   COUL_MAT( :,:) = 0.0d0
   G_MAT(:,:,:,:) = 0.0d0

   call cpu_time(t1)
   write(*,*) '  ALLOCATION-DECLARATION IN COULOMB : ', t1-t0, ' sec'

   !PRE-SCREENING BY CAUCHY-SCHWARZ INEQUALITY || g(a,b,c,d) || <= sqrt( g(a,b,a,b) ) * sqrt( g(c,d,c,d) )
   !CALCULATE Gab = sqrt( g(a,b,a,b))
   do isa = 1, all_atms
      ga = g(isa)
      do isb = isa, all_atms
         gb = g(isb)
         do li = 1, l_max
            do lj = 1, l_max
                do iSPDF = 1, N_SPDF(li,isa)
                  do jSPDF = 1, N_SPDF(lj,isb)
                     do lli = 1, SPDF(li)
                        do llj = 1, SPDF(lj)
                           iloc = idx(isa,li,lli)  + iSPDF
                           jloc = jdx(isb,lj,llj)  + jSPDF
                           ga%shell = shells(li, lli, :)
                           gb%shell = shells(lj, llj, :)

                           if (Gab(iloc,jloc) .EQ. 0.0d0) then

                              call Vee%Vee_int( ga, gb, ga, gb, &
                                   iSPDF + N_SPDF_SUM(li,isa) , jSPDF + N_SPDF_SUM(lj,isb), &
                                   iSPDF + N_SPDF_SUM(li,isa) , jSPDF + N_SPDF_SUM(lj,isb) )

                              Gab(iloc,jloc) = dsqrt(Vee%int_val)
                              Gab(jloc,iloc) = dsqrt(Vee%int_val)

                           end if
                           !write(*,*) 'iloc, jloc, Gab', iloc, jloc, Gab(iloc,jloc)

                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
   end do

   call cpu_time(t2)
   write(*,*) '  TOTAL TIME IN Gab LOOP (SCREENING): ', t2-t1, ' sec'

   !CALCULATION OF Gabcd 
   do isa = 1, all_atms
      ga = g(isa)
      do isb = isa, all_atms
         gb = g(isb)
         do li = 1, l_max
            do lj = 1, l_max
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
                              do isd = isc, all_atms
                                 gd = g(isd)
                                 do li2 = 1, l_max
                                    do lj2 = 1, l_max
                                        do iSPDF2 = 1, N_SPDF(li2,isc)
                                          do jSPDF2 = 1, N_SPDF(lj2,isd)
                                             do lli2 = 1, SPDF(li2)
                                                do llj2 = 1, SPDF(lj2)
                                                   iloc2 = idx(isc,li2,lli2)  + iSPDF2
                                                   jloc2 = jdx(isd,lj2,llj2)  + jSPDF2
                                                   gc%shell = shells(li2, lli2, :)
                                                   gd%shell = shells(lj2, llj2, :)
                                                   !write(*,*) 'iloc, jloc, iloc2, jloc2', &
                                                   !        iloc, jloc, iloc2, jloc2
    

                                                   if ( ( G_MAT(iloc,jloc,iloc2,jloc2) .EQ. 0.0d0 ) .AND. &
                                                        ( Gab(iloc,jloc)*Gab(iloc2,jloc2) .GT. tau ) ) then

                                                      call Vee%Vee_int( ga, gb, gc, gd, &
                                                         iSPDF + N_SPDF_SUM(li,isa) , jSPDF + N_SPDF_SUM(lj,isb), &
                                                         iSPDF2+ N_SPDF_SUM(li2,isc), jSPDF2+ N_SPDF_SUM(lj2,isd) )

                                                      G_MAT(iloc,jloc,iloc2,jloc2) = Vee%int_val
                                                      G_MAT(iloc,jloc,jloc2,iloc2) = Vee%int_val
                                                      G_MAT(jloc,iloc,iloc2,jloc2) = Vee%int_val
                                                      G_MAT(jloc,iloc,jloc2,iloc2) = Vee%int_val

                                                      G_MAT(iloc2,jloc2,iloc,jloc) = Vee%int_val
                                                      G_MAT(iloc2,jloc2,jloc,iloc) = Vee%int_val
                                                      G_MAT(jloc2,iloc2,iloc,jloc) = Vee%int_val
                                                      G_MAT(jloc2,iloc2,jloc,iloc) = Vee%int_val

                                                   end if

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
   call cpu_time(t3)

   write(*,*) '  TOTAL TIME IN Gabcd LOOP: ', t3-t2, ' sec'
   !CALCULATION OF COUL = COUL +  Dcd*(2.0*Gabcd -Gacbd)
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
                           coul_val = 0.0d0
                           do isc = 1, all_atms
                              do isd = 1, all_atms
                                 do li2 = 1, l_max
                                    do lj2 =1, l_max
                                       do iSPDF2 = 1, N_SPDF(li2,isc)
                                          do jSPDF2 = 1, N_SPDF(lj2,isd)
                                             do lli2 = 1, SPDF(li2)
                                                do llj2 = 1, SPDF(lj2)
                                                   iloc2 = idx(isc,li2,lli2)  + iSPDF2
                                                   jloc2 = jdx(isd,lj2,llj2)  + jSPDF2
 
                                                   !COUL_MAT(iloc,jloc) = COUL_MAT(iloc,jloc) + D_MAT(iloc2,jloc2)* &
                                                   !      ( 2.0d0*G_MAT(iloc,jloc,iloc2,jloc2) - G_MAT(iloc,iloc2,jloc,jloc2) )
                                                   coul_val = coul_val + D_MAT(iloc2,jloc2)* &
                                                         ( 2.0d0*G_MAT(iloc,jloc,iloc2,jloc2) - G_MAT(iloc,iloc2,jloc,jloc2) )

                                                end do
                                             end do
                                          end do
                                       end do
                                    end do
                                 end do
                              end do
                              COUL_MAT(iloc,jloc) = coul_val
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
   call cpu_time(t4)

   write(*,*) '  TOTAL TIME IN COUL-LOOP : ', t4-t3, ' sec'
 
   !Ven matrix
   open(13,file='coul_ee.txt',status='unknown')
   do ib = 1,nbas
       write(13,*) COUL_MAT(ib,:)
   end do
   close(13)
 
   write(*,*)
   !allocate( COUL_VEC( nbas * ( nbas + 1 )/2 ) )
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
   !write(*,*) 'Number of elements in COUL_VEC : ', l
   !write(*,*) 'Number of non-zero elements in COUL_VEC : ', m
 
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

  call cpu_time(t5)

  write(*,*) '  TOTAL TIME IN COULOMB : ', t5-t0, ' sec'

end subroutine

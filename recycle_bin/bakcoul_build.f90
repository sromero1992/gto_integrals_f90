subroutine coul_build
   use class_basis_function
   use integral
   use mat_build
   use com_vars 
   type(gaussian)         :: ga, gb, gc, gd     ! Declare gaussian basis 
   type(integ)            :: Vee ! Declaration of integral variables that uses g objects




!########################################## COULOMB Vee INTEGRAL ###################################################
   allocate( G_MAT( nbas, nbas, nbas, nbas) )
   COUL_MAT(:,:)  = 0.0d0
   G_MAT(:,:,:,:) = 0.0d0
   
   !open(20,file='gabcd.txt') 
   do isa = 1, all_atms
      ga = g(isa)
      do isb = 1, all_atms
         gb = g(isb)
         do li = 1, l_max
            do lj =1, l_max
               do iSPDF = 1, N_SPDF(li,isa)
                  do jSPDF = iSPDF, N_SPDF(lj,isb)
                  !do jSPDF = 1, N_SPDF(lj,isb)
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
                                          do jSPDF2 = iSPDF2, N_SPDF(lj2,isd)
                                          !do jSPDF2 = 1, N_SPDF(lj2,isd)
                                             do lli2 = 1, SPDF(li2)
                                                do llj2 = 1, SPDF(lj2)
                                                   iloc2 = idx(isc,li2,lli2)  + iSPDF2
                                                   jloc2 = jdx(isd,lj2,llj2)  + jSPDF2
                                                   write(*,*) 'iloc, jloc, iloc2, jloc2', iloc, jloc, iloc2, jloc2
                                                   gc%shell = shells(li2, lli2, :)
                                                   gd%shell = shells(lj2, llj2, :)
 
                                                   call Vee%Vee_int( ga, gb, gc, gd, &
                                                      iSPDF + N_SPDF_SUM(li,isa) , jSPDF + N_SPDF_SUM(lj,isb), &
                                                      iSPDF2+ N_SPDF_SUM(li2,isc), jSPDF2+ N_SPDF_SUM(lj2,isd) )
                                                   G_MAT(iloc,jloc,iloc2,jloc2) = Vee%int_val  
                                                   G_MAT(iloc,jloc,jloc2,iloc2) = Vee%int_val  
                                                   G_MAT(jloc,iloc,iloc2,jloc2) = Vee%int_val  
                                                   G_MAT(jloc,iloc,jloc2,iloc2) = Vee%int_val  
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
 
                                                   COUL_MAT(iloc,jloc) = COUL_MAT(iloc,jloc) + D_MAT(iloc2,jloc2)* &
                                                         ( 2.0d0*G_MAT(iloc,jloc,iloc2,jloc2) - G_MAT(iloc,iloc2,jloc,jloc2) )

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


end subroutine

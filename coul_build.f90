subroutine coul_build
   use class_basis_function
   use integral
   use mat_build
   use com_vars 
   type(gaussian)         :: ga, gb, gc, gd     ! Declare gaussian basis 
   type(integ)            :: Vee ! Declaration of integral variables that uses g objects
   real(8), allocatable   :: Gab(:,:)
   real(8)                :: coul_val, t0, t1, t2, t3, t4 !, tau2 
   integer                :: ii, jj, ii2, jj2

!########################################## COULOMB Vee INTEGRAL ###################################################
   call cpu_time(t0)
   tau = 1e-10 
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

            ii = N_SPDF_SUM(li,isa)
            do lj = 1, l_max
               
                jj = N_SPDF_SUM(lj,isb)

                do  lli = 1, N_SPDF(li,isa)
                   do llj = 1, N_SPDF(lj,isb)
                      do iSPDF = 1, SPDF(li)

                         iloc = idx(isa,li,lli)  + iSPDF
                         ga%shell = shells(li, iSPDF, :)
                         do jSPDF = 1, SPDF(lj)

                           jloc = idx(isb,lj,llj)  + jSPDF
                           gb%shell = shells(lj, jSPDF, :)

                           if (Gab(iloc,jloc) .EQ. 0.0d0) then

                              call Vee%Vee_int( ga, gb, ga, gb, &
                                            ii + lli, jj + llj, &
                                            ii + lli, jj + llj  )

                              Gab(iloc,jloc) = dsqrt(Vee%int_val)
                              Gab(jloc,iloc) = dsqrt(Vee%int_val)

                           end if

                         end do !End iSPDF

                     end do !End jSPDF
                  end do !End llj
               end do !End lli

            end do !End lj
         end do !End li
  
      end do !End isb
   end do !End isa

   call cpu_time(t2)
   write(*,*) '  TOTAL TIME IN Gab LOOP (SCREENING): ', t2-t1, ' sec'

   !CALCULATION OF Gabcd 
   do isa = 1, all_atms

      ga = g(isa)

      do isb = isa, all_atms

         gb = g(isb)

         do li = 1, l_max

            ii = N_SPDF_SUM(li,isa)

            do lj = 1, l_max

                jj = N_SPDF_SUM(lj,isb)

                do lli = 1, N_SPDF(li,isa)
                   do llj = 1, N_SPDF(lj,isb)
                      do iSPDF = 1, SPDF(li)

                         iloc = idx(isa,li,lli)  + iSPDF
                         ga%shell = shells(li, iSPDF, :)

                         do jSPDF = 1, SPDF(lj)

                           jloc = idx(isb,lj,llj)  + jSPDF
                           gb%shell = shells(lj, jSPDF, :)

                           do isc = 1, all_atms

                              gc = g(isc)

                              do isd = isc, all_atms

                                 gd = g(isd)

                                 do li2 = 1, l_max

                                    ii2 = N_SPDF_SUM(li2,isc)

                                    do lj2 = 1, l_max 
                                        
                                       jj2 = N_SPDF_SUM(lj2,isd)  

                                       do lli2 = 1,  N_SPDF(li2,isc)
                                          do llj2 = 1, N_SPDF(lj2,isd)
                                             do iSPDF2 = 1, SPDF(li2)
                                                
                                                iloc2    = idx(isc,li2,lli2)  + iSPDF2
                                                gc%shell = shells(li2, iSPDF2, :)

                                                do jSPDF2 = 1, SPDF(lj2)

                                                   jloc2    = idx(isd,lj2,llj2)  + jSPDF2
                                                   gd%shell = shells(lj2, jSPDF2, :)

                                                   if ( ( G_MAT(iloc,jloc,iloc2,jloc2) .EQ. 0.0d0 ) .AND. &
                                                        ( Gab(iloc,jloc)*Gab(iloc2,jloc2) .GT. tau ) ) then

                                                      call Vee%Vee_int(    ga   ,    gb   ,     gc    ,     gd   ,&
                                                                        ii + lli, jj + llj, ii2 + lli2, jj2 + llj2)   

                                                      

                                                      G_MAT(iloc,jloc,iloc2,jloc2) = Vee%int_val
                                                      G_MAT(iloc,jloc,jloc2,iloc2) = Vee%int_val
                                                      G_MAT(jloc,iloc,iloc2,jloc2) = Vee%int_val
                                                      G_MAT(jloc,iloc,jloc2,iloc2) = Vee%int_val

                                                      G_MAT(iloc2,jloc2,iloc,jloc) = Vee%int_val
                                                      G_MAT(iloc2,jloc2,jloc,iloc) = Vee%int_val
                                                      G_MAT(jloc2,iloc2,iloc,jloc) = Vee%int_val
                                                      G_MAT(jloc2,iloc2,jloc,iloc) = Vee%int_val

                                                   end if

                                                end do ! jSPDF2 number of ang. momentum s =>1, p = > 3 (px,py,pz)
                                             end do    ! iSPDF2 number of ang. momentum s =>1, p = > 3 (px,py,pz)
                                          end do       ! bare gaussians of llj2 (site bi, ang.mom lj2)
                                       end do          ! bare gaussians of lli2 (site ai, ang.mom li2)

                                    end do ! lj2 ang.mom for site di
                                 end do    ! li2 ang.mom for site ci

                              end do ! di site
                           end do    ! ci site

                        end do ! jSPDF number of ang. momentum s =>1, p = > 3 (px,py,pz)
                     end do    ! iSPDF number of ang. momentum s =>1, p = > 3 (px,py,pz) 
                  end do       ! bare gaussians of llj (site bi, ang.mom lj)
               end do          ! bare gaussians of lli (site ai, ang.mom li)

            end do ! lj ang.mom for site bi
         end do    ! li  ang.mom for site ai

      end do ! bi site 
   end do    ! ai site

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
                           jloc = idx(isb,lj,llj)  + jSPDF
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
                                                   jloc2 = idx(isd,lj2,llj2)  + jSPDF2
 
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

  call cpu_time(t5)

  write(*,*) '  TOTAL TIME IN COULOMB : ', t5-t0, ' sec'

end subroutine

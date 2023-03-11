subroutine basis_moments
  use class_basis_function
  use integral
  use com_vars
  use mat_build
  implicit none
  integer                    :: lmax2, SPD(3), SPDsum(4), Nmom, rmm_ord, ii, jj, ispd, jspd
  integer, allocatable       :: idx2(:,:,:), N_SPD(:,:)
  real(8), allocatable       :: basisMomMat(:,:,:)
  real(8)                    :: Ai(3), Aj(3), alp, bet, cci, ccj, Rm_tmp(10,10,10)


  SPD      = (/ 1 , 3, 6/)
  SPDsum   = (/ 0, 1 , 4, 10/)
  all_atms = sum(n_atm(:))
  lmax2    = 3

  allocate(  idx2( all_atms, lmax2, 6) ) 
  allocate( N_SPD( 3, all_atms)        ) 

  write(*,*) 'Begin of basis moment subroutine'
  
  idx2    = 0
  NBAS    = 0
  do isa = 1, all_atms
     N_SPD(:,isa) = (/ g(isa)%sf, g(isa)%pf, g(isa)%df /)      
     do li = 1, lmax2 !li=1 => S, li=2 => P, li=3 => D
        do lli = 1, N_SPD(li,isa) !Number of S,P,D functions 
           idx2(isa, li, lli) = idx_tmp
           idx_tmp            = idx_tmp + SPD(li) 
        end do
        NBAS = NBAS + N_SPD(li,isa) * SPD(li)
    end do
  end do 
  
  write(*,*) 'NBAS :', NBAS

  Nmom = 3
  allocate( basisMomMat( NBAS, NBAS, Nmom ) )

  basisMomMat = 0.0d0
  do isa = 1, all_atms

     Ai(:) = g(isa)%origin 
     do isb = 1, all_atms

        Aj(:) = g(isb)%origin
        do ia = 1, size(g(isa)%exps)

           alp = g(isa)%exps(ia) 
           do ib = 1, size(g(isb)%exps) 
              bet = g(isb)%exps(ib) 
              rmm_ord = 2
              call getRmoment( alp, bet, Ai, Aj, Rm_tmp, rmm_ord)
              !Looping over Orbitals

              do li = 1, lmax2

                 ii = SPDsum(li)
                 do lj =1, lmax2

                    jj = SPDsum(lj)
                    do lli = 1, N_SPD(li,isa) 

                      cci = g(isa)%coef(ia,lli) 
                      do llj = 1, N_SPD(lj,isb)  

                        ccj = g(isb)%coef(ib,llj) 
                        do iSPD = 1,SPD(li) 

                           iloc = idx2(isa,li,lli)  + iSPD
                           do jSPD = 1, SPD(lj) 

                              jloc = idx2(isb,lj,llj)  + jSPD
                              basisMommat( iloc, jloc, :) = basisMommat( iloc, jloc, :) + &
                                           cci*ccj*Rm_tmp( ii + ispd, jj + jspd,1:Nmom) 
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
 
  write(*,*) 'Writting moments to files'
  
  open(22,file='Rmoments.txt')
  do k = 1, 1
     !write(22,*) 'Rmoment matrix of moment:',k 
     do i = 1, NBAS
        write(22,*) (basisMommat( i, j, k), j=1, NBAS)
     end do
  end do
  close(22)


  write(*,*) 'End of basis moments'
end subroutine


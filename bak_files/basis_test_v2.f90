program basis_test
  use class_basis_function
  implicit none
  real(8)              :: x0(3,2), temp_exps(3), el(2)
  type(gaussian),allocatable :: g(:)! Declare gaussian basis 
  integer  n2, fact2, nb, nbas, n_atm_type, i, j, k, l, ib, jb
  integer, allocatable :: qe(:), qn(:), n_atm(:), sf(:), pf(:), df(:), sf_sup(:), pf_sup(:), df_sup(:) 
  character(10)        :: dum_str
  logical              :: ext_bas
  
  ext_bas = .FALSE.
  open(11,file='ISYMGEN')
  read(11,*) n_atm_type  
  allocate(g(n_atm_type)) !Number of gaussian basis => should be equals to number of atoms/number of contracted gaussians
  allocate(qe(n_atm_type)) !Number of electrons
  allocate(qn(n_atm_type)) ! Number of nuclei
  allocate(n_atm(n_atm_type)) !Number of atoms of n-type
  !S P D functions and supplemental
  allocate(sf(n_atm_type)) 
  allocate(pf(n_atm_type)) 
  allocate(df(n_atm_type)) 
  allocate(sf_sup(n_atm_type)) 
  allocate(pf_sup(n_atm_type)) 
  allocate(df_sup(n_atm_type)) 


  do i = 1, n_atm_type
     read(11,*) qe(i), qn(i) !qn is the atom e.g. qn=1 is hydrogen 
     read(11,*) dum_str
     read(11,*) n_atm(i)
     do j = 1, n_atm(i) + 1
         read(11,*) dum_str
     end do
     read(11,*) g(i)%nb!Number of bare gaussians
     read(11,*) sf(i),pf(i),df(i) 
     read(11,*) sf_sup(i),pf_sup(i),df_sup(i) 
     call g(i)%allocation() !Allocates nb in exponents and coefficients
     !Reading all alpha exponents
     do j = 1, g(i)%nb/3
        read(11,*) g(i)%exps(3*j-2), g(i)%exps(3*j-1), g(i)%exps(3*j)
     end do

!!!!!!!!!!!!!!!!!!!!!!!!!!! Check from here ...
     !turn on extra basis??
     if (ext_bas) then
        nbas = sf(i) + sf_sup(i) + pf(i) + pf_sup(i) + df(i) + df_sup(i) !not complete nbas since it should be sf + 3*pf + 6*df
        do ib = 1, nbas
           do j = 1, g(i)%nb/3
              !read(11,*) g(i)%coef(3*j-2, ib),  g(i)%coef(3*j-1, ib),  g(i)%coef(3*j, ib) 
           end do
        end do
     else
        nbas = sf(i) + pf(i) + df(i) !not complete nbas since it should be sf + 3*pf + 6*df
        !Reading all gaussian contraction coefficients Cij
        do ib = 1, nbas
           do j = 1, g(i)%nb/3
              !read(11,*) g(i)%coef(3*j-2,ib),  g(i)%coef(3*j-1,ib),  g(i)%coef(3*j,ib) 
           end do
           !This ifs are for skipping supplemental basis functions
           if      ( ( nbas .GT. sf(i) ) .AND. (sf_sup(i) .NE. 0 ) ) then
               n2 = sf_sup(i)
           else if ( ( nbas .GT. sf(i)+sf_sup(i)+pf(i) ) .AND. (pf_sup(i) .NE. 0)) then
               n2 = pf_sup(i)
           else if ( ( nbas .GT. sf(i)+sf_sup(i)+pf(i)+pf_sup(i)+df(i)) .AND. (df_sup(i) .NE. 0) ) then
               n2 = df_sup(i)
           else 
               n2 = 0
           end if
           do k = 1, n2
              do j = 1, g(i)%nb/3
                 read(11,*) temp_exps(1),  temp_exps(1),  temp_exps(1) 
              end do
           end do
        end do
     end if 
     !primitive gaussian functions psi_i = sum_j (Cij*phi_j)
  end do
  close(11)

  !Atom positions and SYMBOL information
   open(12,file='SYMBOL')
   do i=1,10
      read(12,*)
   end do
   do i = 1, n_atm_type
      do j = 1, n_atm(i)
         read(12,*)  dum_str, dum_str, g(i)%origin(1), g(i)%origin(2), g(i)%origin(3), dum_str
      end do
   end do
   read(12,*) dum_str, dum_str, el(1), el(2) ! number of electrons up and dn
   close(12)
   !For my ISYMGEN and SYMBOL, I have to same hydrogens



  g(1)%shell(1)= 0
  g(1)%shell(2)= 0
  g(1)%shell(3)= 0
  g(2)%shell = g(1)%shell
  !call g(1)%normalization()
  !call g(1)%print_gaussian_info()  ! Call a class subroutine
  !call g(2)%print_gaussian_info()  ! Call a class subroutine
end program 


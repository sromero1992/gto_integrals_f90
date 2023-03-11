subroutine read_info
   use class_basis_function
   use com_vars
   integer                    :: gs_tmp, li, lmax, nnuc
   real(8)                    :: dum_r(50)
   character(10)              :: dum_str
   integer, allocatable       :: sf(:), pf(:), df(:), sf_sup(:), &
                                 pf_sup(:), df_sup(:), spdf_f(:,:,:) 
   integer                    :: ff_tmp(2)

   lmax = 4
   open(11,file='ISYMGEN')
      read(11,*) nsite
      write(*,*) "Total number of atom type:", nsite 
      !allocation of atom types
      allocate(   site_type(nsite))
      allocate(nbas_tmp(nsite))
      allocate(  nb_tmp(nsite))
      allocate(spdf_f(lmax,2,nsite)) ! 1 = basis, 2 for supplemental basis
      spdf_f(:,:,:) = 0

      !Reading total number of atoms and atom types to construct the objects
      do i = 1, nsite

         read(11,*) tmp1, tmp2   !qn is the atom e.g. qn=1 is hydrogen 
         read(11,*) dum_str
         read(11,*) site_type(i)     !number of atoms of type i
         do j = 1, site_type(i) 
             read(11,*) dum_str
         end do
         read(11,*)
         read(11,*) nb_tmp(i)    !Number of bare gaussians
         write(*,*) 'NUMBER OF BARE GAUSSIANS: ', nb_tmp(i) 

         if (f_bol ) then
            read(11,*) ( spdf_f(li,1,i), li = 1, lmax ) 
            read(11,*) ( spdf_f(li,2,i), li = 1, lmax )
            ! The following turns off suplemental S P D F functions
            if ( ext_bas ) then
               nbas_tmp(i) =  sum( spdf_f(:,:,i) )
               write(*,*) 'number of S P D F functions (includes supplemental): ',nbas_tmp(i)
            else
               nbas_tmp(i) =  sum( spdf_f(:,1,i) )
               write(*,*) 'number of S P D F functions : ',nbas_tmp(i)
            end if

         else

            read(11,*)  ( spdf_f(li,1,i), li = 1, lmax-1)  
            read(11,*)  ( spdf_f(li,2,i), li = 1, lmax-1) 
            nbas_tmp(i) =  sum( spdf_f(:,:,i) )
            write(*,*) 'number of S P D functions (includes supplemental): ',nbas_tmp(i)

         end if
         read(11,*) (dum_r(k), k=1,nb_tmp(i))
         do j = 1, nbas_tmp(i)
            read(11,*) (dum_r(k), k=1,nb_tmp(i))
            read(11,*)
         end do 
       
         write(*,*) 'Atomic type : ',i
         write(*,*) 'functions S and supplemental: ',spdf_f(1,1,i) ,',',spdf_f(1,2,i)
         write(*,*) 'functions P and supplemental: ',spdf_f(2,1,i) ,',',spdf_f(2,2,i)
         write(*,*) 'functions D and supplemental: ',spdf_f(3,1,i) ,',',spdf_f(3,2,i)


         if (i .EQ. nsite) then

            read(11,*) dum_str
            if ( dum_str  .EQ. 'ELECTRONS') then
                write(*,*) 'Successful read of SPDF and total number of equivalent atoms '
            else
                write(*,*) 'You didnt reach the end of ISYMGEN'
            end if

         end if

      end do

      if ( f_bol ) write(*,*) 'functions F and supplemental: ',ff_tmp(1) ,',',ff_tmp(2)
 
      n_contracted_g = sum(site_type(:))
      write(*,*) 'Total number of gaussians (superposition): ', n_contracted_g
    
      ! Allocation for contracted gaussians
      allocate(      g (n_contracted_g)) !Number of contracted gaussian basis 
      ! Allocation for ith atom type
      allocate(      qe (nsite)) !Number of electrons
      allocate(      qn (nsite)) !Nuclei number, for each gaussian
      !S P D functions and supplemental
      allocate(      sf (nsite))
      allocate(      pf (nsite))
      allocate(      df (nsite))
      allocate(  sf_sup (nsite))
      allocate(  pf_sup (nsite))
      allocate(  df_sup (nsite))
      write(*,*) 'Allocation completed'
    
    !Reading everything into objects
      rewind(unit = 11)
      read(11,*) nsite

      do i = 1, nsite

         read(11,*) qe(i), qn(i)   !qn is the atom e.g. qn=1 is hydrogen 
         read(11,*) dum_str
         read(11,*) site_type(i)       !number of atoms of type i
         do j = 1, site_type(i) 
             read(11,*) dum_str
         end do
         read(11,*)
 
         if ( i .EQ. 1 ) then
            g_index = 1            !first contracted gaussian info
         else
            !g_index = sum( site_type( 1: (i-1) ) ) + 1
            g_index = g_index + 1
         end if
         g(g_index)%n_atm = qn(i)  !Atomic number
         read(11,*) g(g_index)%nb  !Number of bare gaussians
         write(*,*) 'total number of bare gaussian (',g_index,') : ', g(g_index)%nb
 
         if ( f_bol ) then
            !Read S P D F functions
            read(11,*) g(g_index)%sf    , g(g_index)%pf    , g(g_index)%df, g(g_index)%ff
            read(11,*) g(g_index)%sf_sup, g(g_index)%pf_sup, g(g_index)%df_sup, g(g_index)%ff_sup

            if (ext_bas) then
               !nbas in gaussian structure is the number of basis functions
               g(g_index)%nbas =  g(g_index)%sf + g(g_index)%pf + g(g_index)%df + g(g_index)%ff + &
                                  g(g_index)%sf_sup + g(g_index)%pf_sup + g(g_index)%df_sup + g(g_index)%ff_sup
               write(*,*) 'total number of basis functions(',g_index,') includes sup func : ', g(g_index)%nbas

            else
               g(g_index)%nbas = g(g_index)%sf + g(g_index)%pf + g(g_index)%df + g(g_index)%ff
               write(*,*) 'total number of basis functions(',g_index,') : ', g(g_index)%nbas
               gs_tmp =  g(g_index)%sf_sup + g(g_index)%pf_sup + g(g_index)%df_sup + g(g_index)%ff_sup
               write(*,*) 'total number of supplemental basis functions(',g_index,') : ', gs_tmp

           end if

         else

            !Read S P D functions
            read(11,*) g(g_index)%sf    , g(g_index)%pf    , g(g_index)%df
            read(11,*) g(g_index)%sf_sup, g(g_index)%pf_sup, g(g_index)%df_sup

            if (ext_bas) then
               !nbas in gaussian structure is the number of basis functions
               g(g_index)%nbas =  g(g_index)%sf + g(g_index)%pf + g(g_index)%df + &
                                  g(g_index)%sf_sup + g(g_index)%pf_sup +  g(g_index)%df_sup
               write(*,*) 'total number of basis functions(',g_index,') includes sup func : ', g(g_index)%nbas

            else
               g(g_index)%nbas = g(g_index)%sf + g(g_index)%pf + g(g_index)%df
               write(*,*) 'total number of basis functions(',g_index,') : ', g(g_index)%nbas
               gs_tmp =  g(g_index)%sf_sup + g(g_index)%pf_sup + g(g_index)%df_sup
               write(*,*) 'total number of supplemental basis functions(',g_index,') : ', gs_tmp

            end if

         end if

      !Allocates nb in exponents and coefficients
      call g(g_index)%allocation() 
      write(*,*) 'size of gaussian basis coefficients C(ib,ibas) :', size(g(g_index)%coef)

      !Reading all alpha exponents
      read(11,*) (g(g_index)%exps(j), j = 1,  g(g_index)%nb  )
      !write(*,*) 'alphas: ',g(g_index)%exps(:)

      if (ext_bas) then

         write(*,*) 'NBAS : ',g(g_index)%nbas
         do ib = 1, g(g_index)%nbas !nbas is the number of basis functions in g(i) basis
            read(11,*) ( g(g_index)%coef(j, ib), j = 1,  g(g_index)%nb )
         end do

      else

         !Reading all gaussian contraction coefficients Cij
         write(*,*) 'Number of Gaussian basis functions: ',g(g_index)%nbas ! This is contains the number of basis functions
         do ib = 1, g(g_index)%nbas  
            !This ifs are for skipping supplemental basis function

            if      ( ( ib .EQ. g(g_index)%sf ) .AND. (g(g_index)%sf_sup .NE. 0 ) ) then
                n2 = g(g_index)%sf_sup
            else if ( ( ib .EQ. g(g_index)%sf+g(g_index)%pf ) .AND. (g(g_index)%pf_sup .NE. 0)) then
                n2 = g(g_index)%pf_sup
            else if ( ( ib .EQ. g(g_index)%sf+g(g_index)%pf+g(g_index)%df ) .AND. (g(g_index)%df_sup .NE. 0) ) then
                n2 = g(g_index)%df_sup
            else
                n2 = 0
            end if

            read(11,*) (g(g_index)%coef(j,ib), j = 1, g(g_index)%nb )
            !write(*,*) (g(g_index)%coef(j,ib), j = 1, g(g_index)%nb )
            do k = 1, n2
               read(11,*) ( temp_exps(1), j = 1, g(g_index)%nb ) 
            end do

         end do

      end if
      g(g_index)%qn = qn(i) 
      g(g_index)%qe = qe(i) 
 
!!Pass all the read basis to any same kind of atom 
      if ( (site_type(i) .GT. 1) ) then

         if (i .GT. 1) then
            g_upper_bound = sum(site_type(1:i))
         else
            g_upper_bound = site_type(1)
         end if

         do k = g_index + 1, g_upper_bound
            g(k) = g(g_index)
         end do

      end if

   end do
   close(11)
 
!!Read atom positions from  SYMBOL 
   open(12,file='SYMBOL')

   do i = 1, 6
      read(12,*)
   end do

   read(12,*) nnuc
   do i = 1, nnuc + 1
      read(12,*)
   end do

   do i = 1, nnuc
       read(12,*)  dum_str, dum_str, g(i)%origin(1), g(i)%origin(2), g(i)%origin(3), dum_str
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

   write(*,*) 'End of read_info'
end subroutine read_info

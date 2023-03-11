module com_vars
   use class_basis_function
   implicit none
   ! module of common variables for ovl,kin,ven matrices.
   type(gaussian),allocatable :: g(:)

   real(8)                    :: temp_exps(3), el(2),  tmp1, tmp2, del, ener_final, e_tol, tau
   
   integer                    :: i, j, k, l, m, n2, ia, ib, nbas, tot_natm_type, all_atms, & 
                                 n_contracted_g, g_index, g_upper_bound, nrec, nelec, &
                                 shells(4,10,3), nsite, isite, jsite 
 
   integer, allocatable       :: qe(:), qn(:), n_atm(:), nb_tmp(:), nbas_tmp(:), site_type(:)
 
   real(8), allocatable       :: S_MAT(:,:), T_MAT(:,:), Ven_MAT(:,:), COUL_MAT(:,:),&
                                 ovlbaby(:), S_VEC(:), T_VEC(:), Ven_VEC(:), D_MAT(:,:), &
                                 eig(:), eigV(:,:), COULOMB(:), COUL_VEC(:), G_MAT(:,:,:,:),&
                                 X_MAT(:,:), eig_tmp(:), eigV_tmp(:,:), H0_MAT(:,:), F_MAT(:,:),&
                                 Fp_MAT(:,:),F_TMP(:,:), F_eigV(:,:), F_eig(:), MAT_TMP(:,:), &
                                 hambaby(:)
 
   logical                    :: ext_bas, le, f_bol, debug, chk_nrlmol

end module com_vars
!##############################################################################################
module mat_build 

  integer,allocatable        :: idx(:,:,:), jdx(:,:,:), N_SPDF(:,:), N_SPDF_SUM(:,:)

  integer                    :: SPDF(4), SPDF_SUM(5), idx_tmp, l_max, &
                                li, lj, lli, llj, iSPDF, jSPDF, iloc, jloc, &
                                isa, isb, isc, isd, li2, lj2, lli2, llj2, iSPDF2, &
                                jSPDF2 , iloc2, jloc2

end module

!#############################################################################################
module diag

  real(8),allocatable   :: work(:),iwork(:)

  character             :: jobz, uplo

  integer               :: info, lwork, liwork, ierr

end module

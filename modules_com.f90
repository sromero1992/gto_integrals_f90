module module_com
   implicit none

   real(8)                    :: temp_exps(3), el(2),  tmp1, tmp2, del, ener_final, e_tol, tau
   
   integer                    :: i, j, k, l, m, n2, ia, ib, ifunc, jfunc, nbas, nfunc, tot_natm_type, all_atms, & 
                                 n_contracted_g, g_index, g_upper_bound, nrec, nelec, ifunc2, jfunc2, &
                                 shells(4,10,3), nsite, isite, jsite 
 
   integer, allocatable       :: qe(:), qn(:), n_atm(:), nbg_tmp(:), nfunc_tmp(:), site_type(:)
 
   real(8), allocatable       :: S_MAT(:,:), T_MAT(:,:), Ven_MAT(:,:), COUL_MAT(:,:),&
                                 ovlbaby(:), S_VEC(:), T_VEC(:), Ven_VEC(:), D_MAT(:,:), &
                                 eig(:), eigV(:,:), COULOMB(:), COUL_VEC(:), G_MAT(:,:,:,:),&
                                 X_MAT(:,:), eig_tmp(:), eigV_tmp(:,:), H0_MAT(:,:), F_MAT(:,:),&
                                 Fp_MAT(:,:),F_TMP(:,:), F_eigV(:,:), F_eig(:), MAT_TMP(:,:), &
                                 hambaby(:)
 
   logical                    :: ext_bas, le, f_bol, debug, chk_nrlmol

end module module_com


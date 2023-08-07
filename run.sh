rm -v test_bas  
gfortran   modules_com.f90 modules_diag.f90 diagon.f90 dens.f90 sort_eivec_val.f90 quick_sort.f90 basis_functions.f90 read_info.f90  okven_matbuild.f90 coul_build.f90 s_transform.f90 get_energy.f90 nuc_rep.f90 hf_scf.f90 basis_moments.f90  getRmoment.f90 main_hf.f90  -o test_bas -llapack -lblas -fbacktrace -fbounds-check -Wall -std=f2008 -fall-intrinsics -g -ffpe-trap=zero,overflow,underflow -O3  #-funroll-loops 
#time ./test_bas > out.txt
#more out.txt 
time ./test_bas

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

program loop_spd
   implicit none
   integer    SPDF(4), SPDF_SUM(4)
   
   SPDF(:)      = [1 3  6 10]
   SPDF_SUM(:)  = [1 4 10 20]
   
   
   idx_tmp = 0
   do isite = 1, n_atm
      do l = 1, l_max
         do ll = 1, nsf(l) 
             idx(isite, l, ll) =  idx_tmp 
             idx_tmp = idx_tmp +   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
    do isite = 1, n_atm
      do l = 1, l_max
         do ll = 1, nsf(l) 
             do ispd = 1, SPD(l) 
                 loc = idx(isite,l,ll)  + ispd
   
   
end program  

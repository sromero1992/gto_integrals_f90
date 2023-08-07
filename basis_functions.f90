module class_basis_function
  implicit none
  private
  public     :: gaussian 
  !real(8)    :: pi = 4.0d0*datan(1.0d0) ! wide class constant


  type gaussian
      !real(8), allocatable :: coef,exps
      real(8), pointer :: coef(:,:) => null()
      real(8), pointer :: exps(:)   => null()
      real(8), pointer :: norm(:)   => null()
      real(8) :: origin(3)
      integer :: nbg,nbg_sup      !nbg   => number of bare gaussians
      integer :: nfunc,nfunc_sup  !nfunc => number of basis functions
      integer :: qe, qn           !q? number of electrons and protos
      integer :: n_atm            !atomic number
      integer :: shell(3)         !shell for angular momentum numbers 
      integer :: nbas             !nbas  => n_sf + 3 n_pf + 6 n_df   (number of basis SPDF functions) 
      integer :: sf , pf, df, ff  !basis functions for S,P,D and F 
      integer :: sf_sup, pf_sup, df_sup, ff_sup !basis functions for supplemental S,P,D and F
      !store all basis function values here for easy access  
      integer :: func(4)       
      integer :: func_sup(4) 
  contains 
      procedure, public :: print_gaussian_info 
      procedure, public :: allocation
      !procedure, public :: normalization
  end type gaussian 
contains


  subroutine print_gaussian_info(this)
     class(gaussian), intent(inout) :: this
             integer                      :: ifunc
       write(*,*) '#########################################################################' 
       write(*,*) 'Atomic number            : ',this%n_atm 
       write(*,*) 'Atomic coordinates       : ',this%origin(:)   
       write(*,*) 'Number of basis functions: ',this%nbg   
       write(*,*) 'Bare Gaussian exponents    : '
       write(*,*)  this%exps(1:this%nbg)  
       write(*,*) 'Bare Gaussian coefficients : '
       do ifunc = 1, this%nfunc
       write(*,*) this%coef(1:this%nbg,ifunc)   
       end do
       !write(*,*) 'Normalization: ',this%norm   
       write(*,*) 'Electronic shell: ',this%shell(:)  
       write(*,*) '#########################################################################' 
  end subroutine
  subroutine allocation(this)
     !this%nbg    => number of bare gaussians (exponents)
     !this%nfunc  => number of basis functions 
     class(gaussian), intent(inout) :: this
       allocate(this%coef(this%nbg,this%nfunc))
       allocate(this%exps(this%nbg))
       allocate(this%norm(this%nbg))
       this%shell(:) = 0
  end subroutine
  !subroutine normalization(this)  
  !   class(gaussian), intent(inout) :: this
  !     integer :: l, m, n, ll, fact2, nb, ia, ib, i
  !     real(8) :: p1, Nor, p ! p1 prefactor
  !     l = this%shell(1)
  !     m = this%shell(2)
  !     n = this%shell(3)
  !     ll= l + m + n
  !     !this%norm is an array of equal lenght to the number of primitive gaussians
  !     !normalization of primitive Gaussians:
  !     this%norm = (2**(2*ll+1.5d0) * this%exps**(ll+1.5d0) / (fact2(2*l-1) * fact2(2*m-1) * fact2(2*n-1) *pi**1.5d0 )  )**0.5d0
  !     p1 = pi**1.5d0 * fact2( 2*l-1 ) * fact2( 2*m-1 )*fact2( 2*n-1 )/2.0d0**ll

  !     do i = 1, this%nbas
  !     Nor = 0.0d0
  !     nb = this%nb
  !     do ia = 1, nb
  !        do ib = 1, nb
  !           p = this%exps(ia) + this%exps(ib)
  !           !Nor = Nor + this%norm(ia) * this%norm(ib) * this%coef(ia) * this%coef(ib) / p**1.5d0
  !           Nor = Nor + this%norm(ia) * this%norm(ib) * this%coef(ia,i) * this%coef(ib,i) / p**1.5d0
  !        end do
  !     end do
  !     write(*,*) ' Nor :',Nor
  !     write(*,*) ' p1  :',p1
  !     Nor = Nor * p1
  !     Nor = Nor**(-0.5d0)
  !     do ia = 1, nb
  !        this%coef(ia,i) = this%coef(ia,i) * Nor 
  !     end do
  !     end do
  !end subroutine
end module class_basis_function
!************************************************************************************

!************************************************************************************
module  integral 
  !This module/class is defines to introduce two gaussian classes into a function
  use class_basis_function 
  implicit none
  private
  public          :: integ 
  !type(gaussian)  :: g1, g2
  !real(8)         :: pi = 4.0d0*datan(1.0d0)

  type integ
      real(8)              :: int_val
      real(8), allocatable :: Mat(:,:)
  contains 
      procedure, public :: S_int
      procedure, public :: T_int
      procedure, public :: Ven_int
      procedure, public :: Vee_int
  end type integ 
contains

  !subroutine S_int(this, g1, g2, i, j)
  subroutine S_int(this, ga, gb, li, lj, ifunc, jfunc, max_iSPDF, max_jSPDF)
     ! Evaluates the overlap between two contracted Gaussians
     use mat_build, only             : SPDF, N_SPDF, N_SPDF_SUM
     use module_com,  only             : shells

     class(integ),    intent(inout) :: this !This is polymorphic
     class(gaussian), intent(inout) :: ga, gb
       !integer,       intent(in)    :: i, j
       integer,       intent(in)    :: li, lj, ifunc, jfunc, max_iSPDF, max_jSPDF
       integer                      :: iSPDF, jSPDF, inbg, jnbg
       real(8)                      :: overlap, s_val

       !s_val = 0.0d0       
       !do ia = 1, g1%nbg
       !   cof(1) = g1%coef(ia,i) 
       !   if ( cof(1) .NE. 0.0d0  ) then
       !      do ib = 1, g2%nbg
       !         cof(2) = g2%coef(ib,j)
       !         if ( cof(2) .NE. 0.0d0  ) then
       !            !s_val = s_val + g1%norm(ia) * g2%norm(ib) * g1%coef(ia,i) * g2%coef(ib,j) * & 
       !            s_val = s_val + cof(1) * cof(2) * &
       !                    overlap( g1%exps(ia), g1%shell, g1%origin, g2%exps(ib), g2%shell, g2%origin) 
       !         end if 
       !      end do
       !   end if
       !end do

       if ( allocated(this%Mat) ) deallocate(this%Mat)
       allocate(this%Mat(max_iSPDF,max_jSPDF) )
       this%Mat(:,:) = 0.0d0

       do iSPDF = 1, SPDF(li)

          ga%shell = shells(li,iSPDF,:)

          do jSPDF = 1, SPDF(lj)

             gb%shell = shells(lj,jSPDF,:)

             s_val = 0.0d0
             do inbg = 1, ga%nbg

                if ( ga%coef( inbg,ifunc)== 0) cycle

                do jnbg = 1, gb%nbg

                   if ( gb%coef( jnbg,jfunc) == 0) cycle
             
                   s_val =  s_val + ga%coef(inbg,ifunc) * gb%coef(jnbg,jfunc) * & 
                          overlap( ga%exps(inbg), ga%shell, ga%origin, gb%exps(jnbg), gb%shell, gb%origin) 

                end do
             end do   
             this%Mat(iSPDF,jSPDF) = s_val
          end do
       end do

  end subroutine

  subroutine T_int(this, g1, g2, i, j)
    ! Evaluates the kinetic energy between two contracted Gaussians
     class(integ),    intent(inout) :: this
     class(gaussian), intent(inout) :: g1, g2
       integer,       intent(in)    :: i, j
       integer                      :: ia, ib !, i, j
       real(8)                      :: kinetic, kin_val, cof(2)

       kin_val = 0.0d0 
       do ia = 1, g1%nbg
          cof(1) =  g1%coef(ia,i)
          if ( cof(1)  .NE. 0.0d0  ) then
             do ib = 1, g2%nbg
                cof(2) = g2%coef(ib,j) 
                if ( cof(2)  .NE. 0.0d0  ) then
                   !kin_val = kin_val + g1%norm(ia) * g2%norm(ib) * g1%coef(ia,i) * g2%coef(ib,j) * &
                   kin_val = kin_val + cof(1) * cof(2) * &
                             kinetic( g1%exps(ia), g1%shell, g1%origin, g2%exps(ib), g2%shell, g2%origin)
                end if
             end do
          end if
       end do
       this%int_val = kin_val
  end subroutine

  subroutine Ven_int(this, g1, g2, i, j, C) 
     !This subroutine evaluates the coulomb potential integral within electron-nucleus
     !a  : contracted Gaussian a (Basis function object) 
     !b  : contracted Gaussian b (Basis function object) 
     !C  : center of nucleus
     !Self-note, still needs to be multiplied by Z_n
     class(integ),    intent(inout) :: this
     class(gaussian), intent(inout) :: g1, g2
       real(8),       intent(in)    :: C(3)
       real(8)                      :: V_col, nuc_attraction, cof(2)
       integer,       intent(in)    :: i, j
       integer                      :: ia, ib 

       V_col = 0.0d0
       do ia = 1, g1%nbg
          cof(1) = g1%coef(ia,i)
          if (cof(1) .NE. 0.0d0) then 
             do ib = 1, g2%nbg
                cof(2) = g2%coef(ib,j)
                if (   cof(2)  .NE. 0.0d0  ) then
                   !V_col = V_col + g1%norm(ia) * g2%norm(ib) * g1%coef(ia,i) * g2%coef(ib,j) * &
                   V_col = V_col + cof(1) * cof(2) * & 
                        nuc_attraction( g1%exps(ia), g1%shell, g1%origin, g2%exps(ib),  g2%shell, g2%origin, C)
                end if
             end do
          end if
       end do
       this%int_val = V_col
  end subroutine

!459 function elec_repulsion( a, lmn1, RA, b, lmn2, RB, c, lmn3, RC, d, lmn4, RD) result(val)
  subroutine Vee_int(this, ga, gb, gc, gd, i, j, k, l)
     !This subroutine evaluates the four center electron repulsion 
     !a is the contracted gaussian a (Basis function object)
     !b is the contracted gaussian b (Basis function object)
     !c is the contracted gaussian c (Basis function object)
     !d is the contracted gaussian d (Basis function object)
     !i,j,k,l are the index of basis set
     class(integ),    intent(inout) :: this
     class(gaussian), intent(inout) :: ga, gb, gc, gd
       integer,       intent(in)    :: i, j, k, l
       integer                      :: ia, ib, ic, id
       real(8)                      :: val, elec_repulsion, c(4), tau

       tau = 1.0d-10      
       val = 0.0d0
       do ia = 1, ga%nbg
          c(1) =  ga%coef(ia,i) 
          if ( c(1) .NE. 0.0d0 ) then 
             do ib = 1, gb%nbg
                c(2) = gb%coef(ib,j)
                if ( c(2) .NE. 0.0d0 ) then 
                   do ic = 1, gc%nbg
                      c(3) = gc%coef(ic,k)
                      if ( c(3) .NE. 0.0d0 ) then 
                         do id = 1, gd%nbg
                            c(4) = gd%coef(id,l)
                            if ( (   c(4)  .NE. 0.0d0 ) .AND. &
                                 ( ( ga%exps(ia) + gb%exps(ib) ) * ( gc%exps(ic) + gd%exps(id) ) .GT. tau )     ) then  

                               val  = val + c(1) * c(2) * c(3) * c(4)  * &
                                      elec_repulsion( ga%exps(ia), ga%shell, ga%origin, &
                                                      gb%exps(ib), gb%shell, gb%origin, &
                                                      gc%exps(ic), gc%shell, gc%origin, &
                                                      gd%exps(id), gd%shell, gd%origin   )
                            end if
                         end do
                      end if
                   end do
                end if
             end do
          end if
       end do
       this%int_val = val
  end subroutine

end module integral
!***********************************************************************************

!******************************************************************
recursive function fact2(n) result(n2)
   !this function calculates n!!
   implicit none
   integer, intent(in), value :: n
   integer                    :: n2
   if ( ( n .LE. 1 ) .and. (n .GE. -1) ) then 
       n2 = 1 
   else if ( n .LT. -1 ) then
       n2 = 0 
   else
       n2 = n*fact2(n-2)
   end if
end function

!******************************************************************
recursive function E(i, j, t, Qx, a, b) result(res)
   !This function calculates the recursive value of Hermite Gaussian coefficients
   !a   => exponent of Gaussian A
   !b   => exponent of Gaussian B
   !i,j => orbital angular momentum of Gaussian A and gassian B, e.g. XA^i *XB^j
   !t   => number of nodes in Hermite, this depends on the kind of integral e.g. overlap t=0, dipole moment t=1
   !Qx  => distance within Gaussian centers (XA-XB)
   implicit none
   integer, intent(in), value :: i, j, t
   real(8), intent(in), value :: a, b, Qx
   real(8)                    :: p, q, res 
   p = a + b
   q = a*b/p
   if ( (t .LT. 0) .OR. (t .GT. (i+j) ) ) then
      res = 0.0d0
   else if ( (i .EQ. 0) .AND. (j .EQ. 0 ) .AND. (t .EQ. 0) ) then
      res = exp( -q * Qx**2 )
   else if ( j .EQ. 0 ) then
      res = E(i-1, j, t-1, Qx, a, b)/(2*p) - q*Qx/a* E(i-1, j, t, Qx, a, b) + (t+1)* E(i-1, j, t+1, Qx, a, b)
   else
      res = E(i, j-1, t-1, Qx, a, b)/(2*p) + q*Qx/b* E(i, j-1, t, Qx, a, b) + (t+1)* E(i, j-1, t+1, Qx, a, b)
   end if
end function

!*********************************************************************
function  overlap(a,lmn1, RA, b, lmn2, RB) result(ovl)
   !Evaluates the overlap integral between gaussians A and B
   !a    => Gaussian exponent of A
   !b    => Gaussian exponent of B
   !RA   => center of Gaussian A
   !RB   => center of Gaussian B
   !lmn1 => orbital angular momentum (e.g. (1,0,0) -> (l,m,n) for Gaussian A
   !lmn2 => orbital angular momentum for Gaussian B
   implicit none
   interface 
     function E(i, j, t, Qx, a, b) result(res)
       integer, intent(in), value :: i, j, t
       real(8), intent(in), value :: a, b, Qx
       real(8)                    :: p, q, res
     end function E
   end interface
   real(8), intent(in) :: a, RA(3), b, RB(3)
   integer, intent(in) :: lmn1(3), lmn2(3)
   integer             :: l(2), m(2), n(2) 
   real(8)             :: Sx, Sy, Sz, Q(3), pi, p, ovl

   pi = 4.0d0*datan(1.0d0)
   l = (/lmn1(1), lmn2(1) /)   
   m = (/lmn1(2), lmn2(2) /)   
   n = (/lmn1(3), lmn2(3) /)   

   p  = a + b
   Q  = RA - RB
   Sx = E(l(1), l(2), 0, Q(1), a, b)
   Sy = E(m(1), m(2), 0, Q(2), a, b)
   Sz = E(n(1), n(2), 0, Q(3), a, b)
   ovl= Sx*Sy*Sz*(pi/p)**1.5d0

end function
!***********************************************************************
function kinetic(a, lmn1, RA, b, lmn2, RB) result(kin)
   !Evaluates the kinetic energy integral between gaussians A and B
   !a    => Gaussian exponent of A
   !b    => Gaussian exponent of B
   !RA   => center of Gaussian A
   !RB   => center of Gaussian B
   !lmn1 => orbital angular momentum (e.g. (1,0,0) -> (l,m,n) for Gaussian A
   !lmn2 => orbital angular momentum for Gaussian B
   implicit none
   interface
     function E(i, j, t, Qx, a, b) result(res)
       integer, intent(in), value :: i, j, t
       real(8), intent(in), value :: a, b, Qx
       real(8)                    :: p, q, res
     end function E
   end interface
   real(8), intent(in) :: a, RA(3), b, RB(3)
   integer, intent(in) :: lmn1(3), lmn2(3)
   integer             :: l(2), m(2), n(2), lmn1_tmp(3), lmn2_tmp(3)
   real(8)             :: overlap , Tab0, Tab1, Tab2, kin ,tmp

   l  = (/ lmn1(1), lmn2(1) /)
   m  = (/ lmn1(2), lmn2(2) /)
   n  = (/ lmn1(3), lmn2(3) /)
   lmn1_tmp = lmn1
   lmn2_tmp = lmn2

   ! Tab = -0.5 ( Dij^2*Skl*Smn + Sij*Dkl^2*Smn + Sij*Skl*Dmn^2)
   ! Dij^2 =j(j-1) *S(i,j-2) - 2b*(2j+1)* S(i,j) + 4b^2*S(i,j+2) 

   ! Tab0 definition (j*(j-1)*S(i,j-2))
   lmn2_tmp(1) = l(2) - 2
   Tab0 =  l(2) * ( l(2) - 1.0d0 )*overlap(a, lmn1_tmp, RA, b,lmn2_tmp, RB)     
   lmn2_tmp(1) = lmn2(1)
   lmn2_tmp(2) = m(2) - 2
   Tab0 = Tab0 + m(2) * ( m(2) - 1.0d0 )*overlap(a, lmn1_tmp, RA, b, lmn2_tmp, RB)     
   lmn2_tmp(2) = lmn2(2)
   lmn2_tmp(3) = n(2) - 2
   Tab0 = Tab0 + n(2) * ( n(2) - 1.0d0 )*overlap(a, lmn1_tmp, RA, b, lmn2_tmp, RB)     
   lmn2_tmp(3) = lmn2(3)
   Tab0 = -0.5d0 * Tab0

   ! Tab1 definition of 2b*(2j+1)* S(i,j)

   !IF ANY INCONSISTENCY IN T MATRIX, CHANGE TO tmp
   tmp = overlap(a, lmn1_tmp, RA, b,lmn2_tmp, RB)
   !write(*,*)'overlap before Tab1:', tmp
   Tab1 = tmp ! I had to do this because Tab1 store junk...
   !Tab1 = overlap(a, lmn1_tmp, RB, b, lmn2_tmp, RB)
   !write(*,*) 'overlap after Tab1 :', Tab1
   Tab1 = b * ( 2.0d0 * ( l(2) + m(2) + n(2) ) + 3.0d0 ) * Tab1
   !Tab1 = b * ( 2.0d0 * ( l(2) + m(2) + n(2) ) + 3.0d0 ) * &
   !       overlap(a, lmn1_tmp, RA, b,lmn2_tmp, RB)

   ! Tab2 definition of 4b^2*S(i,j+2)
   lmn2_tmp(1) = l(2) + 2
   Tab2 =  overlap(a, lmn1_tmp, RA, b,lmn2_tmp, RB)*1.0d0     
   lmn2_tmp(1) = lmn2(1)
   lmn2_tmp(2) = m(2) + 2
   Tab2 = Tab2 + overlap(a, lmn1_tmp, RA, b,lmn2_tmp, RB)*1.0d0     
   lmn2_tmp(2) = lmn2(2)
   lmn2_tmp(3) = n(2) + 2
   Tab2 = Tab2 + overlap(a, lmn1_tmp, RA, b,lmn2_tmp, RB)*1.0d0     
   lmn2_tmp(3) = lmn2(3)
   Tab2 = -2.0d0 * b**2.0d0 * Tab2 

   kin = Tab0 + Tab1 + Tab2
end function 


!*********************************************************************************
recursive function boys1(n,T) result(boys_val)
   !This function calculates the boys function integral for Coulomb
   !n  : Boys' function index
   !T  : Boys' function variable
   implicit none
   integer, intent(in), value :: n
   real(8), intent(in), value :: T
   real(8)                    :: boys_val, pi
   
   pi =4.0d0*atan(1.0d0)
   if ( (T .GT. 1e-8) .AND. (n .EQ. 0) ) then
      boys_val = ( pi/(4.0d0*T) )**0.5d0*derf( T**0.5d0 )
   else if ( T .LE. 1e-8) then
      !First two terms of Taylor expansion to avoid division over zero
      boys_val = 1.0d0/(2.0d0*n+1.0d0) - T/(2.0d0*n+3.0d0) 
   else 
      boys_val = ( ( 2.0d0*n - 1.0d0 )*boys1(n-1,T) - dexp(-T) )/ (2.0d0*T)
   end if

end function



!******************************************************************************************
recursive function R(t, u, v, n, p, RPC1, RPC2, RPC3) result(val) 
   !This function computes the auxiliary  Hermite Coulomb integral function R^n_{tuv}(p,RP,RC)
   !t,u,v:  are order of the Coulomb Hermite derivative in x,y,z
   !n    :  order of Boys function 
   !RPC  :  vector distance between P and C        
   !F_n(Tb)
   implicit none
   interface
     function boys1(n,T) result(boys_val)
       integer, intent(in), value :: n
       real(8), intent(in), value :: T
       real(8)                    :: boys_val, pi
     end function boys1
   end interface 
   integer, intent(in), value :: t, u, v, n   
   real(8), intent(in), value :: p, RPC1, RPC2, RPC3
   real(8)                    :: val, Tb !, boys1 !, RPCnorm
   
   Tb = p*( RPC1**2.0d0 + RPC2**2.0d0 + RPC3**2.0d0) !p*RCP^2
   val = 0.0d0
   if ( (t .EQ. 0) .AND. (u .EQ. 0) .AND. (v .EQ. 0) ) then
        val = val + ( -2.0d0*p )**n * boys1(n,Tb) !R^n_{000}=(-2p)^n*F_n(p*Rpc^2)
   else if ( (t .EQ. 0) .AND. (u .EQ. 0) ) then
        if ( v .GT. 1 ) val = val + ( v - 1.0d0 )*R(t, u, v-2, n+1, p, RPC1, RPC2, RPC3)
        val = val + RPC3 * R( t, u, v-1, n+1, p, RPC1, RPC2, RPC3)
   else if ( t .EQ. 0 ) then
        if ( u .GT. 1 ) val = val + ( u - 1.0d0 )*R(t, u-2, v, n+1, p, RPC1, RPC2, RPC3) 
        val = val + RPC2 * R( t, u-1, v, n+1, p, RPC1, RPC2, RPC3)
   else 
        if ( t .GT. 1 ) val = val + ( t - 1.0d0 )*R(t-2, u, v, n+1, p, RPC1, RPC2, RPC3)
        val = val + RPC1 * R( t-1, u, v, n+1, p, RPC1, RPC2, RPC3)
   end if
end function 

!*********************************************************************************

function nuc_attraction(a, lmn1, RA, b, lmn2, RB, RC) result(val)
   !This function evaluates the nuclear attraction within two Gaussians
   !a    :  exponent of Gaussian A
   !b    :  exponent of Gaussian B
   !lmn1 :  integer vector(3) with the angular momentum of Ax^l Ay^m Az^n
   !lmn2 :  integer vector(3) with the angular momentum 
   !RA   :  origin of Gaussian A
   !RB   :  origin of Gaussian B
   !RC    :  origin of nuclear center C 
   implicit none
   interface 
     function E(i, j, t, Qx, a, b) result(res)
       integer, intent(in), value :: i, j, t
       real(8), intent(in), value :: a, b, Qx
       real(8)                    :: p, q, res
     end function E
     function R(t, u, v, n, p, RPC1, RPC2, RPC3) result(val)
       integer, intent(in), value :: t, u, v, n   
       real(8), intent(in), value :: p, RPC1, RPC2, RPC3
       real(8)                    :: val, Tb, boys1, RPCnorm
     end function R
   end interface
   integer, intent(in)   :: lmn1(3), lmn2(3)
   real(8), intent(in)   :: a, b, RA(3), RB(3), RC(3)
   integer               :: l(2), m(2), n(2), t, u, v
   real(8)               :: p, RP(3), RPC(3), RAB(3), val, pi,E1(3)

   pi = 4.0d0*datan(1.0d0)
   l = (/ lmn1(1), lmn2(1) /)
   m = (/ lmn1(2), lmn2(2) /)
   n = (/ lmn1(3), lmn2(3) /)
   p  = a + b 
   RAB = RA - RB
   RP  = (a*RA + b*RB)/p !Gaussian product center
   RPC = RP - RC

   val = 0.0d0
   do t = 0, (l(1) + l(2) + 1) 
      E1(1) = E(l(1), l(2), t, RAB(1), a, b)
      if ( E1(1) .NE. 0.0d0) then
         do u = 0, (m(1) + m(2) +1)
            E1(2) = E(m(1), m(2), u, RAB(2), a, b)
               if ( E1(2) .NE. 0.0d0) then
               do v = 0, (n(1) + n(2) +1)
                  E1(3) = E(n(1), n(2), v, RAB(3), a, b)
                  if ( E1(3) .NE. 0.0d0) then

                     val = val + E1(1) * E1(2) * E1(3) * &
                                 R(t, u, v, 0, p, RPC(1), RPC(2), RPC(3))
               
                  end if
               end do
            end if
         end do 
      end if
   end do
   val = -val*2.0d0*pi/p
end function 



function elec_repulsion( a, lmn1, RA, b, lmn2, RB, c, lmn3, RC, d, lmn4, RD) result(val)
   !This function evaluates the coulomb potential integral within electron-electron (repulsion)
   !a, b, c, d  :  are the respective gaussian exponents of gaussians A, B, C and D
   !lmni        :  are the respective respective angular momentum for gaussians A, B, C and D (e.g. i=1 <=> a%shell =(0,0,0))
   !RA, RB, RC, RD : are the respective gaussian center or each gaussian A, B, C and D
     implicit none
     interface 
       function E(i, j, t, Qx, a, b) result(res)
         integer, intent(in), value :: i, j, t
         real(8), intent(in), value :: a, b, Qx
         real(8)                    :: p, q, res
       end function E
       function R(t, u, v, n, p, RPC1, RPC2, RPC3) result(val)
         integer, intent(in), value :: t, u, v, n   
         real(8), intent(in), value :: p, RPC1, RPC2, RPC3
         real(8)                    :: val, Tb, boys1, RPCnorm
       end function R
     end interface
     real(8),       intent(in) :: a, b, c, d, RA(3), RB(3), RC(3), RD(3)
     integer,       intent(in) :: lmn1(3), lmn2(3), lmn3(3), lmn4(3)
     real(8)                   :: RP1(3), RP2(3), RP1P2(3), p1, p2, alp, &
                                  RAB(3), RCD(3), pi, val, E1(3), E2(3), val0
     integer                   :: t, u, v, tau, nu, phi, l(4), m(4), n(4)

     pi = 4.0d0*datan(1.0d0)
     l(1:4) = (/lmn1(1), lmn2(1), lmn3(1), lmn4(1)/) 
     m(1:4) = (/lmn1(2), lmn2(2), lmn3(2), lmn4(2)/)
     n(1:4) = (/lmn1(3), lmn2(3), lmn3(3), lmn4(3)/)
     p1     = a + b            ! composite exponent of gaussians A and B
     p2     = c + d            ! composite exponent of gaussians C and D
     RP1    = (a*RA + b*RB)/p1 ! composite center of gaussians A and B                  
     RP2    = (c*RC + d*RD)/p2 ! composite center of gaussians A and B                  
     RP1P2  = RP1 - RP2
     alp    = p1*p2/( p1 + p2 )        
     RAB    = RA - RB
     RCD    = RC - RD

     val = 0.0d0
     do tau = 0, ( l(3) + l(4) +1 )
        E2(1) = E(l(3), l(4), tau, RCD(1), c, d)
        if ( E2(1) .NE. 0.0d0 ) then
           do nu = 0, ( m(3) + m(4) +1 )
              E2(2) = E(m(3), m(4), nu , RCD(2), c, d)
              if ( E2(2) .NE. 0.0d0 ) then
                 do phi = 0, ( n(3) + n(4) +1 )
                    E2(3) = E(n(3), n(4), phi, RCD(3), c, d)
                    if ( E2(3) .NE. 0.0d0 ) then
                       val0 = (-1.0d0)**( tau + nu + phi ) * E2(1) * E2(2) * E2(3)
                       do t = 0, ( l(1) + l(2) + 1 )
                          E1(1) = E(l(1), l(2), t  , RAB(1), a, b)
                          if ( E1(1) .NE. 0.0d0 ) then
                             do u = 0, ( m(1) + m(2) + 1 )
                                E1(2) = E(m(1), m(2), u  , RAB(2), a, b)
                                if ( E1(2) .NE. 0.0d0 ) then
                                   do v = 0, ( n(1) + n(2) +1 )
                                      E1(3) = E(n(1), n(2), v  , RAB(3), a, b)

                                      if ( E1(3) .NE. 0.0d0 ) then
                                         val = val + E1(1) * E1(2) * E1(3) * val0 * &
                                                    R(t+tau, u+nu, v+phi, 0, alp, RP1P2(1), RP1P2(2), RP1P2(3))

                                      end if
                                   end do
                                end if
                             end do
                          end if
                       end do
                    end if
                 end do
              end if 
           end do
        end if
     end do
     val = val*2.0d0*pi**2.5d0 / (p1*p2*(p1 + p2)**0.5d0)
end function



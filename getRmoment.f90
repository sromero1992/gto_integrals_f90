subroutine getRmoment(alp,bet,A,B,W,rmm_ord) 

implicit none 

real(8),  intent(in) :: alp, bet, A(3), B(3)                 
integer,  intent(in) :: rmm_ord                              
! Generated for lmax=2 (orb_max=10), rmm_ord=2 (Nmom=10) 

!                      W(orb_max,orb_max,Nmom)               
real(8), intent(out) :: W(10,10,10)    

! parameters                                    
integer :: lmax=2                               
integer :: orb_max=10                           
integer :: Nmom=10                              
integer :: ML(3,10) !ML(3,orb_max)              
integer :: L3(3,10) !L3(3,Nmom)                 


!         F(0:lmax, 0:lmax¸ 0:rmm_ord)   
real(8) :: F(0:2, 0:2, 0:2)               

real(8) :: AT, xs, P, Fn(3)                                
real(8) :: Qn(7) !Qn(lmax+rmm_ord)                         
real(8) :: V(3,3,3,3) !V(lmax+1,lmax+1,Nmom ,3) 
! looping indices    
integer :: i,j,k,ix  
! The order of ML and L3 needs to be dicussed  !!!!!! ML
ML=reshape([     &        
    0,   0,   0, & 
    1,   0,   0, & 
    0,   1,   0, & 
    0,   0,   1, & 
    2,   0,   0, & 
    0,   2,   0, & 
    0,   0,   2, & 
    1,   1,   0, & 
    1,   0,   1, & 
    0,   1,   1 ], shape(ML))

L3=reshape([     &        
    0,   0,   0, & 
    1,   0,   0, & 
    0,   1,   0, & 
    0,   0,   1, & 
    2,   0,   0, & 
    0,   2,   0, & 
    0,   0,   2, & 
    1,   1,   0, & 
    1,   0,   1, & 
    0,   1,   1 ], shape(L3))


AT=1.0d0/(alp+bet)          
! Eq. Q_{n}=\int_{-\infty}^{\infty}x^{2n}e^{-\alpha x^{2}}\,dx=\frac{(2n-1)!!}{\alpha^{n}2^{n+1}}\sqrt{\frac{\pi}{\alpha}}
! P=alp+bet                         
! (ix+jx-1)!!/(2*P)                 
!Q1=1/(2*(alp+bet)= 0.5d0*AT        
!Q2=3*Q1*Q1                         
!Q3=5*Q1*Q2                         
!Q4=7*Q1*Q3                         
!Qn can be calculated recrusively   
Qn(1)=0.5d0*AT                      
do i=2,lmax+rmm_ord                 
  Qn(i)=(2.0*i-1)*Qn(1)*Qn(i-1)     
enddo                             

xs=0.0       
do ix=1,3    

  xs=xs+(alp*bet*((A(ix)-B(ix))**2))*AT  
  P=(alp*A(ix)+bet*B(ix))*AT  
  Fn(1)=P-A(ix)               
  Fn(2)=P-B(ix)               
  Fn(3)=P                     
  
! This part can be unrolled to reduce the number of ops   
  do i=0,lmax                                                 
    do j=0,lmax                                               
      do k=0,rmm_ord                                          
        F(i,j,k) = Fn(1)**i *Fn(2)**j *Fn(3)**k               
      enddo                                                   
    enddo                                                     
  enddo                                                       



  !(x-A)^0 (x-B)^0 x^0 
  V(1,1,1,ix) =         1 

  !(x-A)^0 (x-B)^1 x^0 
  V(1,2,1,ix) =         F(0,1,0)       

  !(x-A)^0 (x-B)^2 x^0 
  V(1,3,1,ix) =         F(0,2,0)        +                 Qn(1)

  !(x-A)^1 (x-B)^0 x^0 
  V(2,1,1,ix) =         F(1,0,0)       

  !(x-A)^1 (x-B)^1 x^0 
  V(2,2,1,ix) =         F(1,1,0)        +                 Qn(1)

  !(x-A)^1 (x-B)^2 x^0 
  V(2,3,1,ix) =         F(1,2,0)        +       F(1,0,0) *Qn(1) +  2.0 *F(0,1,0) *Qn(1)
               

  !(x-A)^2 (x-B)^0 x^0 
  V(3,1,1,ix) =         F(2,0,0)        +                 Qn(1)

  !(x-A)^2 (x-B)^1 x^0 
  V(3,2,1,ix) =         F(2,1,0)        +  2.0 *F(1,0,0) *Qn(1) +       F(0,1,0) *Qn(1)
               

  !(x-A)^2 (x-B)^2 x^0 
  V(3,3,1,ix) =         F(2,2,0)        +       F(2,0,0) *Qn(1) +  4.0 *F(1,1,0) *Qn(1)  &
                +       F(0,2,0) *Qn(1) +                 Qn(2)

  !(x-A)^0 (x-B)^0 x^1 
  V(1,1,2,ix) =         F(0,0,1)       

  !(x-A)^0 (x-B)^1 x^1 
  V(1,2,2,ix) =         F(0,1,1)        +                 Qn(1)

  !(x-A)^0 (x-B)^2 x^1 
  V(1,3,2,ix) =         F(0,2,1)        +  2.0 *F(0,1,0) *Qn(1) +       F(0,0,1) *Qn(1)
               

  !(x-A)^1 (x-B)^0 x^1 
  V(2,1,2,ix) =         F(1,0,1)        +                 Qn(1)

  !(x-A)^1 (x-B)^1 x^1 
  V(2,2,2,ix) =         F(1,1,1)        +       F(1,0,0) *Qn(1) +       F(0,1,0) *Qn(1)  &
                +       F(0,0,1) *Qn(1)

  !(x-A)^1 (x-B)^2 x^1 
  V(2,3,2,ix) =         F(1,2,1)        +  2.0 *F(1,1,0) *Qn(1) +       F(1,0,1) *Qn(1)  &
                +       F(0,2,0) *Qn(1) +  2.0 *F(0,1,1) *Qn(1) +                 Qn(2)
               

  !(x-A)^2 (x-B)^0 x^1 
  V(3,1,2,ix) =         F(2,0,1)        +  2.0 *F(1,0,0) *Qn(1) +       F(0,0,1) *Qn(1)
               

  !(x-A)^2 (x-B)^1 x^1 
  V(3,2,2,ix) =         F(2,1,1)        +       F(2,0,0) *Qn(1) +  2.0 *F(1,1,0) *Qn(1)  &
                +  2.0 *F(1,0,1) *Qn(1) +       F(0,1,1) *Qn(1) +                 Qn(2)
               

  !(x-A)^2 (x-B)^2 x^1 
  V(3,3,2,ix) =         F(2,2,1)        +  2.0 *F(2,1,0) *Qn(1) +       F(2,0,1) *Qn(1)  &
                +  2.0 *F(1,2,0) *Qn(1) +  4.0 *F(1,1,1) *Qn(1) +  2.0 *F(1,0,0) *Qn(2)  &
                +       F(0,2,1) *Qn(1) +  2.0 *F(0,1,0) *Qn(2) +       F(0,0,1) *Qn(2)
               

  !(x-A)^0 (x-B)^0 x^2 
  V(1,1,3,ix) =         F(0,0,2)        +                 Qn(1)

  !(x-A)^0 (x-B)^1 x^2 
  V(1,2,3,ix) =         F(0,1,2)        +       F(0,1,0) *Qn(1) +  2.0 *F(0,0,1) *Qn(1)
               

  !(x-A)^0 (x-B)^2 x^2 
  V(1,3,3,ix) =         F(0,2,2)        +       F(0,2,0) *Qn(1) +  4.0 *F(0,1,1) *Qn(1)  &
                +       F(0,0,2) *Qn(1) +                 Qn(2)

  !(x-A)^1 (x-B)^0 x^2 
  V(2,1,3,ix) =         F(1,0,2)        +       F(1,0,0) *Qn(1) +  2.0 *F(0,0,1) *Qn(1)
               

  !(x-A)^1 (x-B)^1 x^2 
  V(2,2,3,ix) =         F(1,1,2)        +       F(1,1,0) *Qn(1) +  2.0 *F(1,0,1) *Qn(1)  &
                +  2.0 *F(0,1,1) *Qn(1) +       F(0,0,2) *Qn(1) +                 Qn(2)
               

  !(x-A)^1 (x-B)^2 x^2 
  V(2,3,3,ix) =         F(1,2,2)        +       F(1,2,0) *Qn(1) +  4.0 *F(1,1,1) *Qn(1)  &
                +       F(1,0,2) *Qn(1) +       F(1,0,0) *Qn(2) +  2.0 *F(0,2,1) *Qn(1)  &
                +  2.0 *F(0,1,2) *Qn(1) +  2.0 *F(0,1,0) *Qn(2) +  2.0 *F(0,0,1) *Qn(2)
               

  !(x-A)^2 (x-B)^0 x^2 
  V(3,1,3,ix) =         F(2,0,2)        +       F(2,0,0) *Qn(1) +  4.0 *F(1,0,1) *Qn(1)  &
                +       F(0,0,2) *Qn(1) +                 Qn(2)

  !(x-A)^2 (x-B)^1 x^2 
  V(3,2,3,ix) =         F(2,1,2)        +       F(2,1,0) *Qn(1) +  2.0 *F(2,0,1) *Qn(1)  &
                +  4.0 *F(1,1,1) *Qn(1) +  2.0 *F(1,0,2) *Qn(1) +  2.0 *F(1,0,0) *Qn(2)  &
                +       F(0,1,2) *Qn(1) +       F(0,1,0) *Qn(2) +  2.0 *F(0,0,1) *Qn(2)
               

  !(x-A)^2 (x-B)^2 x^2 
  V(3,3,3,ix) =         F(2,2,2)        +       F(2,2,0) *Qn(1) +  4.0 *F(2,1,1) *Qn(1)  &
                +       F(2,0,2) *Qn(1) +       F(2,0,0) *Qn(2) +  4.0 *F(1,2,1) *Qn(1)  &
                +  4.0 *F(1,1,2) *Qn(1) +  4.0 *F(1,1,0) *Qn(2) +  4.0 *F(1,0,1) *Qn(2)  &
                +       F(0,2,2) *Qn(1) +       F(0,2,0) *Qn(2) +  4.0 *F(0,1,1) *Qn(2)  &
                +       F(0,0,2) *Qn(2) +                 Qn(3)
enddo !ix

do j=1, orb_max  ! X Y Z 1                               
 do i=1, orb_max                                         
  do k=1, Nmom                                           
    W(i,j,k)=V( ML(1,i)+1, ML(1,j)+1, L3(1,k)+1, 1 )   & 
            *V( ML(2,i)+1, ML(2,j)+1, L3(2,k)+1, 2 )   & 
            *V( ML(3,i)+1, ML(3,j)+1, L3(3,k)+1, 3 )     
  enddo                                                  
 enddo                                                   
enddo                                                  

W=W*exp(-xs)*(sqrt(3.14159265358979324D0*AT)**3)       

end subroutine


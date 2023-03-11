module class_basis_function
  implicit none
  private
  public  :: gaussian !, normalize
  real    :: pi = 3.1415926535897931d0 ! Class-wide private constant

  type gaussian
      !real(8), allocatable :: coef,exps
      real(8), pointer :: coef(:,:) => null()
      real(8), pointer :: exps(:) => null()
      real(8), pointer :: norm(:) => null()
      real(8) :: origin(3)
      integer :: shell(3), nb, nbas !nb = number of bare gaussians, nbas number of basis functions (S,P,D) functions without moment
  contains 
      procedure, public :: print_gaussian_info 
      procedure, public :: allocation
      procedure, public :: normalization
      !procedure, public :: read_gaussian
  end type gaussian 
contains
  subroutine print_gaussian_info(this)
     class(gaussian), intent(inout) :: this
       integer :: ii
       write(*,*) '#######################################' 
       write(*,*) 'Atomic coordinates: ',this%origin(:)   
       write(*,*) 'Number of bare gaussians: ',this%nb   
       write(*,*) 'Number of basis functions: ',this%nbas   
       write(*,*) 'Gaussian exponents: ',this%exps(:)  
       !write(*,*) 'Gaussian coefficients: ', this%coef(:)   
       !write(*,*) 'Normalization: ',this%norm   
       !write(*,*) 'Electronic shell: ',this%shell(:)  
       write(*,*) '#######################################' 
  end subroutine
  subroutine allocation(this)
     class(gaussian), intent(inout) :: this
       allocate( this%coef( this%nb, this%nbas ))
       allocate( this%exps( this%nb ))
       allocate( this%norm( this%nb ))
       !write(*,*)  'Size of exponents pointer array: ',size(this%exps)
  end subroutine
!  subroutine read_gaussian(this)
!     class(gaussian), intent(inout) :: this

!  end subroutine
  subroutine normalization(this)  
     class(gaussian), intent(inout) :: this
       integer :: l, m, n, ll, fact2, nb, ia, ib, ibas
       real(8) :: p1, Nor, p ! p1 prefactor
       l = this%shell(1)
       m = this%shell(2)
       n = this%shell(3)
       ll= l + m + n
       !this%norm is an array of equal lenght to the number of primitive gaussians
       !it is normalized first the primivitve bare gaussian functions 
       this%norm = (2**(2*ll+1.5d0) * this%exps**(ll+1.5d0) / (fact2(2*l-1) * fact2(2*m-1) * fact2(2*n-1) *pi**1.5d0 )  )**0.5d0
       p1 = pi**1.5d0 * fact2( 2*l-1 ) * fact2( 2*m-1 )*fact2( 2*n-1 )/2**ll

       !Check this ibas loopi*******
       do ibas = 1, this%nbas
          Nor = 0.0d0
          nb = this%nb
          do ia = 1, nb
             do ib = 1, nb
                   p = this%exps(ia) + this%exps(ib)
!                   Nor = Nor + this%norm(ia) * this%norm(ib) * this%coef(ia,ibas) * this%coef(ib) / p**1.5d0
!                   Nor = Nor + this%norm(ia) * this%norm(ib) * this%coef(ia) * this%coef(ib) / p**1.5d0
             end do
          end do
          Nor = Nor * p1
          Nor = Nor**(-0.5d0)
         do ia = 1, nb
!               this%coef(ia) = this%coef(ia) * Nor 
         end do
       end do
  end subroutine
end module class_basis_function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive function fact2(n) result(n2)
   integer, intent(in)  :: n
   integer              :: n2
   if ( ( n .LE. 1 ) .and. (n .GE. -1) ) then 
       n2 = 1 
   else if ( n .LT. -1 ) then
       n2 = 0 
   else
       n2 = n*fact2(n-2)
   end if
end function

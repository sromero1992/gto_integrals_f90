program geo_opt

implicit none

real(8) :: pi, dh, err,  pot, optp, rad(4),R_fix(3),d_fix(3),min_rad,rad_test
real(8), dimension(:,:), allocatable :: R, D, Rop, Dop, R_fod_fix 
real(8), dimension(:), allocatable :: th, ph, phold,thold
integer :: i, j, k, m, ios, coun, N, coun2, Nup, Ndn, AN, N_fod_fix
integer :: selec

pi = 4.d0 * atan(1.d0)


open(20,file = 'FOD_INPUT.DAT')
read(20,*) R_fix(1), R_fix(2), R_fix(3)
read(20,*) N
read(20,*) N_fod_fix
allocate( R_fod_fix(N_fod_fix,3) )
do i = 1,N_fod_fix
    read(20,*) R_fod_fix(i,1) ,R_fod_fix(i,2) ,R_fod_fix(i,3)
    R_fod_fix(i,:) = R_fod_fix(i,:) - R_fix(:)
end do
close(20)

!Allocation of the arrrays
allocate( R(N,3) )
allocate( Rop(N,3) )
allocate( D( N*(N-1)/2 ,3) )
allocate( Dop( N*(N-1)/2 ,3) )
allocate( th(N) )
allocate( ph(N) )
allocate( phold(N) )
allocate( thold(N) )


!Radius of the shells, the center is out of the shells
rad(1) = 0.9
rad(2) = 1.5
rad(3) = 2.2
rad(4) = 3.3
dh = 1 !initial step size
R(:,:) = 0.0d0
Rop(:,:) = 100.0d0 !optimum positions
optp = 10000.0d0 !optimum potential
Dop(:,:) = 0.0d0 !Optimum distance
thold(:) = 0.0d0
phold(:) = 0.0d0
coun = 0
coun2 =  0
!Begin of the monte carlo
call srand( time() )
write(*,*) ' Your radius should be smaller than : '
min_rad = 10.0
do i=1,N_fod_fix
   rad_test = norm2(R_fix(:)-R_fod_fix(i,:))
   if ( min_rad > rad_test ) then
       min_rad = rad_test  
   end if
end do
write(*,*) min_rad
write(*,*) ' radius = '
read(*,*) rad(1)
selec = 1

call get_R(R,rad(selec), th,thold,ph,phold,pi,N,dh)

do m = 1,5000
   pot = 0 
   k = 0 
   D(:,:) = 1e10
   do i=1,N-1
       do j = i+1, N
       !try to get the distance for spin ups, 1:Nup-1 (one core elec) and the other Nup:N
           k = k+1 
           D(k,:) = R(j,:) - R(i,:)!D(1,:) is the difference between 2 and 1
           pot = pot + 1/norm2( D(k,:) )
       end do
   end do
!   do i=1,N
!       do j=1,N_fod_fix
!          d_fix = R(i,:) - R_fod_fix(j,:)
!          pot = pot + 1/norm2(d_fix(:))
!       end do
!   end do

   !Testing if the potential improves
   if (pot .LT. optp) then
        optp = pot 
        !Dop = D
        Rop = R 
        thold = th
        phold = ph
        coun2 = 0 
       !write(*,*) 'The potential improved is : ',pot,' at iter ',m 
   else
        coun2 = coun2 + 1
        coun = coun + 1
        if (coun .GE. 100) then
            dh = dh*0.85
            coun = 0
        end if
   end if
! stoping criteria
  if (coun2 .GE. 300) then
      write(*,*) ' No improvement reached, exit at iter ',m
      exit
  end if
  call get_R(R,rad(selec), th,thold,ph,phold,pi,N,dh)
end do
write(*,*) 'The optimal potential is : ',optp
write(*,*) 'Iteration No. : ',m
write(*,*) 'The set of optimal positions (no core electrons)  are : '
!do m = 1, N
       !write(*,*) Rop(m,:)
!end do

do i=1,N
    Rop(i,:) = Rop(i,:) +R_fix(:)
end do
do i=1,N_fod_fix
    R_fod_fix(i,:) = R_fod_fix(i,:) +R_fix(:)
end do
call print_fod(Rop,N,R_fod_fix ,N_fod_fix, R_fix)
!call system('jmol plot.xyz')


end

!###############################################################

subroutine get_R( V_pos, radii, th_get, thold_get, ph_get, phold_get, pi, N_get, dh)
real(8) :: radii, pi, dh, thold_get(N_get), phold_get(N_get)
real(8) :: V_pos(N_get,3), th_get(N_get), ph_get(N_get)

do i = 1,N_get
    th_get(i) = thold_get(i) + pi*(rand()-0.5)*2*dh
    ph_get(i) = phold_get(i) + 2*pi*(rand()-0.5)*2*dh
    !These are X,Y and Z positions
    V_pos(i,1) = radii*sin(th_get(i))*cos(ph_get(i))
    V_pos(i,2) = radii*sin(th_get(i))*sin(ph_get(i))
    V_pos(i,3) = radii*cos(th_get(i))
end do

end subroutine

!############################################################
subroutine print_fod(Rop,N,R_fix,N_fix,R_fix_atm)!,Nup,Ndn)
implicit none
integer :: m,N,N_fix,ios
real(8) :: Rop(N,3),R_fix(N_fix,3), R_fix_atm(3)

do m = 1, N 
        write(*,*) Rop(m,:)
end do
do m = 1,N_fix
        write(*,*) R_fix(m,:)
end do
 

open(15,file = 'FRMIDT',access = 'append')
write(15,*) N + N_fix 
open(10,file = 'plot.xyz')
write(10,*,iostat=ios) N + N_fix + 1 
write(10,*,iostat=ios) 'plot points '

do m = 1,N_fix
          write(15,*) R_fix(m,:)
          write(10,*,iostat=ios)'Hs',R_fix(m,1),R_fix(m,2),R_fix(m,3),0.001
end do

do m = 1, N 
           write(15,*) Rop(m,:)
           write(10,*,iostat=ios)'Hs',Rop(m,1),Rop(m,2),Rop(m,3),0.001
end do
write(10,*,iostat=ios)'Pu',R_fix_atm(1),R_fix_atm(2),R_fix_atm(3),0.001
close(10)
close(15) 

end subroutine

!
!    geom_util is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    geom_util is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with geom_util.  If not, see <https://www.gnu.org/licenses/>.
!
! Copyright (C) 2019 - 2020 Holger Kruse, Institute of Biophysics of the CAS, Czechia

program geom_util
  implicit none

integer ( kind = 4 ), parameter :: dim_num = 3
real ( kind = 8 ) centroid(dim_num)

real(8), allocatable:: v(:,:), xyz(:,:)
integer, allocatable :: al(:)
integer nat,i,j,ial
character(120) cname,arg,ns,frmt,mode

! FOR PLANE DISTANCE
integer nat1,nat2,k
real(8), allocatable :: coord1(:,:), coord2(:,:)
!integer, allocatable :: ic1(:),ic2(:)

integer n ! polygon size

print*,''
print*,'x-------------------------x'
print*,'x       geom_util         x'
print*,'x-------------------------x'
print*,''
print*,'version: 0.4b'
print*,''
print*,' computes controid of a 3d polygon'
print*,' e.g. for a ring in a molecule'
print*,''
print*,' also distance between fitted planes'
print*,''
print*,'usage:'
print*,'Geometry centroid <file xyz> <int ring size> <index list of ring>'
print*,'Geometry plane <file xyz> <index list plane1> <index list plane2>'
print*,''
print*,''
print*,' <file xyz> = xmol or tmol, prints centroid in given units!'
print*,''
print*,'!! atom numbers must be given in order of the ring !!'
print*,'!! (clock -or anticlock-wise does not matter)      !!'
print*,''


! file with coordinates
call getarg(1,mode)
call getarg(2,cname)

call get_nat(cname,nat)
write(6,'(a,I4,a)')' --> found ',nat,' atoms!'
! now read coordinates
allocate(xyz(3,nat))
call rdfile(cname,nat,xyz)



if(trim(mode)=='centroid') then
! ring size
call getarg(3,arg)
read(arg,*) n
write(6,'(a,I2)')' --> ring size: ', n

! ring info
allocate(al(n))
call getarg(4,arg)
call atlist(arg,al,ial)

allocate(v(3,n))

! order atom list
!call int_bsort(n,al)
write(ns,*) n
write(frmt,'(a)') "(a,2x,"//trim(adjustl(ns))//"(I1,x))"
write(6,frmt),' --> ring atoms:',al

! assign polygon vertices
do j=1,n
 do i=1,nat
  if(i==al(j)) then
   v(1:3,j)=xyz(1:3,i)
  endif
 enddo
enddo


call r8mat_transpose_print ( dim_num, n, v, '  ring coordinates:' )

! compute!
call polygon_centroid_3d ( n, v, centroid )
print*, 'CENTROID'
print*, centroid

elseif(trim(mode)=='plane') then


allocate(al(100))


! ATOMS OF PLANE 1
al=0
call getarg(3,arg)
call atlist(arg,al,ial)

nat1=0
do i=1,nat
 if(any(al==i)) nat1=nat1+1
enddo
print*,'Atoms plane 1', nat1
allocate(coord1(3,nat1))

k=0
do i=1,nat
 if(any(al==i)) then
   k=k+1
   coord1(1:3,k)=xyz(1:3,i)
 endif
enddo

! ATOMS OF PLANE 1
al=0
call getarg(4,arg)
call atlist(arg,al,ial)

nat2=0
do i=1,nat
 if(any(al==i)) nat2=nat2+1
enddo
print*,'Atoms plane 2', nat2
allocate(coord2(3,nat2))

k=0
do i=1,nat
 if(any(al==i)) then
   k=k+1
   coord2(1:3,k)=xyz(1:3,i)
 endif
enddo


print*,coord1
print*,''
print*,coord2

call dist_plane(nat1,nat2,coord1,coord2)
endif






end

! count atoms in XMOL or TMOL file
subroutine get_nat(filen,nat)
implicit none
integer nat,io
character(*) filen
character(120) aa
logical fstr

nat=0
open(newunit=io,file=filen)
do
  read(io,'(a)',end=555) aa
  if( fstr(aa,'$coord') ) then
  do
   read(io,'(a)') aa
   nat=nat+1
   if(fstr(aa,'$')) goto 666
  enddo
  else
    nat=nat+1
  endif
enddo
555 nat=nat-2
666 close(io)
end

! very simple extraction of the coordinates from TMOL or XMOL files
! no elements
subroutine rdfile(filen,nat,xyz)
implicit none
character(*) filen
character(120) aa
real(8) xyz(3,nat)
integer io,nat,i,j
logical fstr


open(newunit=io,file=filen)
do
  read(io,'(a)') aa
  if( fstr(aa,'$coord') ) then
  do i=1,nat
   read(io,*) xyz(1:3,i),aa
  enddo
  xyz=xyz*0.52917720859d0
  exit
  else
   read(io,'(a)') aa ! title
   do i=1,nat
    read(io,*) aa,xyz(1:3,i)
   enddo
   exit
  endif
enddo
close(io)
end


subroutine bsort(n,e)
! bubble sort the vector e
implicit none
integer i,j,k,l,nn,n,ii
real(8) e(n),tt
character(80) cc
logical order

nn=n
order=.false.
do
if(order) exit
order=.true.
 do i=1,nn-1
    if (e(i).gt.e(i+1) ) then ! swap
      tt=e(i)
      e(i)=e(i+1)
      e(i+1)=tt
      order = .false.
     endif
 enddo
nn=nn-1
enddo

return
end subroutine


subroutine int_bsort(n,e)
! bubble sort the vector e
implicit none
integer i,j,k,l,nn,n,ii
integer e(n),tt
character(80) cc
logical order

nn=n
order=.false.
do
if(order) exit
order=.true.
 do i=1,nn-1
    if (e(i).gt.e(i+1) ) then ! swap
      tt=e(i)
      e(i)=e(i+1)
      e(i+1)=tt
      order = .false.
     endif
 enddo
nn=nn-1
enddo

return
end subroutine



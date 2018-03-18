program Cavity_Flow
	IMPLICIT NONE
    character(len=7) outputname
    character(len=7) inputname
	integer :: i,j,k,l,m,n
    integer :: Nx=129
    integer :: Ny=129
    real :: s,err
    real :: delx
    real :: dely
    real :: delt=0.001
    real :: Re
    real :: hx
    real :: hy
	real :: vor0(129,129)
    real :: vor(129,129)
    real :: fin0(129,129)
    real :: fin(129,129)
    real :: u(129,129)
    real :: v(129,129)

    real :: a0
    real :: a1
    real :: a2
    real :: a3
    real :: a4

    real :: check

50	format(2x,E20.10,E20.10,E20.10,E20.10)
    delx=1.0/(Nx-1.0)
    dely=1.0/(Ny-1.0)
    hx=1.0/(delx*delx)
    hy=1.0/(dely*dely)

!l control the reynald number
l=8

!initialization
if(l.eq.1) then
    do j=1,Ny
      do i=1,Nx
        vor0(i,j)=0.0
        fin0(i,j)=0.0
        vor(i,j)=0.0
        fin(i,j)=0.0
        if(j==Ny) then
          u(i,j)=1.0
        else
          u(i,j)=0.0
        endif
        v(i,j)=0.0
      enddo
    enddo
else
  	write(inputname,"(i3,a4)")99+l,'.txt'
	open(21,file=inputname) 
	do j=1,Ny
      do i=1,Nx
        read(21,50)fin0(i,j),vor0(i,j),u(i,j),v(i,j)
        write(*,*)i,j,fin0(i,j),vor0(i,j),u(i,j),v(i,j)
        fin(i,j)=fin0(i,j)
        vor(i,j)=vor0(i,j)
      enddo
    enddo
	close(21)
endif

write(*,*)delx,dely,hx,hy,delt

!Reynald number

	if(l==1) then
  		Re=0.1
	else if(l==2) then
  		Re=1.0
    else if(l==3) then
      	Re=10.0
    else if(l==4) then
      	Re=100.0
    else if(l==5) then
      	Re=200.0
    else if(l==6) then
      	Re=400.0
    else if(l==7) then
      	Re=1000.0
    else if(l==8) then
      	Re=3200
    endif 

  	k=0
do while(.true.)
	k=k+1
  
!vortity and stream function inside the cavity 
	do j=2,(Ny-1)
        s=0.0
    	do i=2,(Nx-1)
        	a0=1.0/delt+(2.0*(hx+hy))/Re
            a1=-u(i+1,j)/(2.0*delx)+hx/Re
            a2=u(i-1,j)/(2.0*delx)+hx/Re
            a3=-v(i,j+1)/(2.0*dely)+hy/Re
            a4=v(i,j-1)/(2.0*dely)+hy/Re

            vor(i,j)=(a1*vor0(i+1,j)+a2*vor(i-1,j)+a3*vor0(i,j+1)+a4*vor(i,j-1)+vor0(i,j)/delt)/a0
            fin(i,j)=(vor(i,j)+hx*(fin0(i+1,j)+fin(i-1,j))+hx*(fin0(i,j+1)+fin(i,j-1)))/(2.0*(hx+hy))
            
            s=s+abs(vor(i,j)-vor0(i,j))
        enddo
        if(j.eq.2) then
          err=s
        else
          if(s.gt.err) then
            err=s
          else
            continue
          endif
        endif
    enddo

!vorticity on the wall
    do i=2,(Nx-1)
    	vor(i,1)=-2.0*hy*(fin(i,2))
        vor(i,Ny)=-2.0*hy*(fin(i,(Ny-1))+1.0*dely)
    enddo

    do i=2,(Ny-1)
		vor(1,i)=-2.0*hx*fin(2,i)
        vor(Nx,i)=-2.0*hx*fin((Nx-1),i)
    enddo

!velocity inside the cavity
    do j=2,(Ny-1)
    	do i=2,(Nx-1)
			u(i,j)=(fin(i,j+1)-fin(i,j-1))/(2.0*dely)
            v(i,j)=-(fin(i+1,j)-fin(i-1,j))/(2.0*delx)           
        enddo
    enddo

!check convergence(L-inf norm)
    check=(1.0/(Nx-2))*err
    if(k.eq.1) then
      continue
    else
      if(check<0.000001) exit
    endif

!Last moment=Next moment        
	do j=1,Ny
    	do i=1,Nx
			vor0(i,j)=vor(i,j)
           	fin0(i,j)=fin(i,j)           
        enddo
    enddo  
enddo

!write to txt file
    write(outputname,"(i3,a4)")100+l,'.txt'
	open(20,file=outputname)  
    do j=1,Ny
      do i=1,Nx
    	write(20,50)fin(i,j),vor(i,j),u(i,j),v(i,j)
      enddo
    enddo
    write(20,*)k
    close(20)


write(*,*)delx,dely,hx,hy,delt

    
end program Cavity_Flow
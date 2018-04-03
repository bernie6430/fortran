program Cavity_Flow_refine
	IMPLICIT NONE
    character(len=7) outputname
    character(len=7) inputname
	integer :: i,j,k,l
    integer :: Nx=257
    integer :: Ny=257
    integer :: Nox
    integer :: Noy
    real :: s,err
    real :: delx
    real :: dely
    real :: delt=0.001
    real :: Re
    real :: hx
    real :: hy
	real :: vor0(257,257)
    real :: vor(257,257)
    real :: fin0(257,257)
    real :: fin(257,257)
    real :: u(257,257)
    real :: v(257,257)

    real :: a0
    real :: a1
    real :: a2
    real :: a3
    real :: a4

    real :: b0(129,129)
    real :: b1(129,129)
    real :: b2(129,129)
    real :: b3(129,129)


    real :: check

50	format(2x,E20.10,E20.10,E20.10,E20.10)
    delx=1.0/(Nx-1.0)
    dely=1.0/(Ny-1.0)
    hx=1.0/(delx*delx)
    hy=1.0/(dely*dely)

!l control the reynald number
l=9

!initialization


	Nox=129
    Noy=129
  	write(inputname,"(i3,a4)")100+1,'.txt'
	open(21,file=inputname) 
	do j=1,Noy
      do i=1,Nox
        read(21,50)b0(i,j),b1(i,j),b2(i,j),b3(i,j)
      enddo
    enddo
	close(21)
    
    b1(1,1)=-2.0*hx*b0(2,1)
    b1(Nox,1)=-2.0*hx*b0((Nox-1),1)
    b1(1,Noy)=-2.0*hx*b0(2,(Noy-1))
    b1(Nox,Noy)=-2.0*hx*b0((Nox-1),Noy)

    do j=1,Ny
      do i=1,Nx
        if((Mod(j,2).eq.0).and.(Mod(i,2).eq.0)) then
           fin0(i,j)=0.25*(b0(i/2,j/2)+b0(i/2+1,j/2)+b0(i/2,j/2+1)+b0(i/2+1,j/2+1))
   		   vor0(i,j)=0.25*(b1(i/2,j/2)+b1(i/2+1,j/2)+b1(i/2,j/2+1)+b1(i/2+1,j/2+1))
           u(i,j)=0.25*(b2(i/2,j/2)+b2(i/2+1,j/2)+b2(i/2,j/2+1)+b2(i/2+1,j/2+1))
           v(i,j)=0.25*(b3(i/2,j/2)+b3(i/2+1,j/2)+b3(i/2,j/2+1)+b3(i/2+1,j/2+1))
        else if((Mod(j,2).ne.0).and.(Mod(i,2).eq.0)) then
           fin0(i,j)=0.5*(b0(i/2,(j+1)/2)+b0(i/2+1,(j+1)/2))
           vor0(i,j)=0.5*(b1(i/2,(j+1)/2)+b1(i/2+1,(j+1)/2))
           u(i,j)=0.5*(b2(i/2,(j+1)/2)+b2(i/2+1,(j+1)/2))
           v(i,j)=0.5*(b3(i/2,(j+1)/2)+b3(i/2+1,(j+1)/2))
        else if((Mod(j,2).eq.0).and.(Mod(i,2).ne.0)) then
           fin0(i,j)=0.5*(b0((i+1)/2,j/2)+b0((i+1)/2,j/2+1))
           vor0(i,j)=0.5*(b1((i+1)/2,j/2)+b1((i+1)/2,j/2+1))
           u(i,j)=0.5*(b2((i+1)/2,j/2)+b2((i+1)/2,j/2+1))
           v(i,j)=0.5*(b3((i+1)/2,j/2)+b3((i+1)/2,j/2+1))
        else if((Mod(j,2).ne.0).and.(Mod(i,2).ne.0)) then
           fin0(i,j)=b0((i+1)/2,(j+1)/2)
           vor0(i,j)=b1((i+1)/2,(j+1)/2)
           u(i,j)=b2((i+1)/2,(j+1)/2)
           v(i,j)=b3((i+1)/2,(j+1)/2)
        endif
        fin(i,j)=fin0(i,j)
   		vor(i,j)=vor0(i,j)
      enddo
    enddo


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
    else if(l==9) then
      	Re=5000
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
    write(outputname,"(i3,a4)")200+l,'.txt'
	open(20,file=outputname)  
    do j=1,Ny
      do i=1,Nx
    	write(20,50)fin(i,j),vor(i,j),u(i,j),v(i,j)
      enddo
    enddo
    write(20,*)k
    close(20)

end program Cavity_Flow_refine
c--------------------------------------------------------------------------c
c     Forward gravity anomaly of 3D Fault Model
c     Using Spectral method by FAST FOURIER TRANSFORM 
c--------------------------------------------------------------------------c
c     input/output Grid file: Goldensoftware Surfer's Ascii format
c     gravity anomaly unit: g.u.
c
c     parameter define: please set in a .txt file
c     input grid: set in the parameters file
c     output: see the parameters file 
c
c         -----------------
c        /\  θ            \ 
c       /  \    FAULT 3D    \
c      /    \    Ver 1.3     \ 2b
c     \      \    Sep.2017    \
c      \      \      2a        \
c       \     |-----------------| h1
c        \    /                 /
c         \  /                 /
c          \/  α             /
c          |------------------|   h2
c--------------------------------------------------------------------------c
c     Code by Dr. Shi CHEN  E-mail:chenshi80@gmail.com
c     Date:2013-11-28
c     Last coding:2017-09-29
c--------------------------------------------------------------------------c
      PROGRAM Sp_Fault3D
	PARAMETER  (pi=3.14159267,G=6.672)
       
      real,ALLOCATABLE::x(:),gb(:,:),xi(:),gi(:),temp_q(:),dg(:,:)
      real,allocatable::bodies(:,:),gall(:,:)
      character*20,allocatable::btype(:,:)
      character*40 inname,fgrid,fbname,fgname,fbou1,fbou2
	Integer m,n,k,method,datatype,expendtype,calum,iftype,infft
	real z,h,x0,d,lou,id,ilou,ix0,iz,ih,Xmin,Xmax,Ymin,Ymax,Zmin,Zmax
      real dx,dy,a,b,length,width,h1,h2,alfa,pointx,pointy,topw,botw
      real tx1,tx2,tx3,rx,ry,rms,xita,fac,facg
      real t1,t2,tt1,tt2               ! this is for time check      
c      inname='cmd1.txt'
	call readparams(inname)
      write(*,*) '========================================'
      write(*,*) 'c         -----------------            c'
      write(*,*) 'c        /\  θ            \           c'
      write(*,*) 'c       /  \    FAULT 3D    \          c'
      write(*,*) 'c      /    \    Ver 1.3     \ 2b      c'
      write(*,*) 'c     \      \    Sep.2017    \        c'
      write(*,*) 'c      \      \      2a        \       c'
      write(*,*) 'c       \     |-----------------| h1   c'
      write(*,*) 'c        \    /                 /      c'
      write(*,*) 'c         \  /    WELCOME!!!   /       c'
      write(*,*) 'c          \/  α             /        c'
      write(*,*) 'c          |------------------|   h2   c'
      write(*,*) '========================================'

      call ReadCmd(inname,fgrid,fbname,fgname,ilen,iunit
     +,infft,igtype,isalldata,isoutbou,mlen,rou0,cx0,cy0)
       
      ALLOCATE (bodies(mlen,9),btype(mlen,2))
      call Readbodies(fbname,bodies,btype,mlen)
      call ReadGRD(fgrid,m,n,Xmin,Xmax,Ymin,Ymax,Zmin,Zmax)
      if (ilen.eq. 2) then
       fac=1000.0
      else
       fac=1.0
      endif
      if (iunit.eq. 1) then
       facg=0.1
      else
       facg=1.0
      endif
      ALLOCATE (gb(m,n),dg(m,n),gall(m,n))

	dx=(Xmax-Xmin)/(M-1)
	dy=(Ymax-Ymin)/(N-1)

	write(*,*) 'Dx,Dy',dx,dy
	write(*,*) 'M,RB,X0,Y0',mlen,rou0,cx0,cy0
      !mGal=10 g.u.
      !km=1000 m
	call cpu_time(t1)
c	tt1=secnds(0.0)
	gall=0.0
      do i=1,mlen
      rou=bodies(i,1)-rou0
      h1= bodies(i,2)
      h2= bodies(i,3)
      length = bodies(i,4)
      width = bodies(i,5)
      alfa = bodies(i,6)
      xita = bodies(i,7)
      pointx = bodies(i,8) !cx0
      pointy = bodies(i,9) !cy0
      if ((abs(alfa)<1.0E-30).or.(abs(alfa-180)<1.0E-30)) then
       STOP 'The value of alfa is not 0 or 180!'
      endif
      if ((abs(xita)<1.0E-30).or.(abs(xita-180)<1.0E-30))then
       STOP 'The value of xita is not 0 or 180!'
      endif
      !btype
      if (trim(btype(i,1)).eq.'R') then
       iftype=2
      else
       iftype=1
      endif
      call CaluVz(gb,m,n,dx,dy,Xmin,Xmax,Ymin,Ymax,Zmin,Zmax,length
     +,width,h1,h2,alfa,xita,rou,pointx,pointy,rms,igtype,iftype,infft)
       
      if (trim(btype(i,1)).eq.'R') then
       write(*,*) i,btype(i,1)
       call transg(gb,dg,m,n,1)
       gb=dg
      endif
      if (isalldata.eq.1) then
       fbou1=trim(btype(i,2))//'.grd'
	call WriteGRD(fbou1,m,n,gb*facg*fac
     +	,Xmin,Xmax,Ymin,Ymax,Zmin,Zmax)
      endif
      if (isoutbou.eq.1) then
       fbou1=trim(btype(i,2))//'_t.bln'
       fbou2=trim(btype(i,2))//'_b.bln'
       if (trim(btype(i,1)).eq.'L') then
       call outbou(fbou1,width,length,pointx,pointy,1,xita)
       else
       call outbou(fbou1,width,length,pointx,pointy,2,xita)
       endif

       if (abs(alfa-90)>1.0E-30) then
        w1=(h2-h1)/tand(alfa)
        if (trim(btype(i,1)).eq.'L') then
         call outbou2(fbou2,width,length,pointx,pointy,w1,1,xita)
        else
         call outbou2(fbou2,width,length,pointx,pointy,w1,2,xita)  
        endif     
       else
       if (trim(btype(i,1)).eq.'L') then
        call outbou(fbou2,width,length,pointx,pointy,1,xita)
       else
        call outbou(fbou2,width,length,pointx,pointy,2,xita)
       endif
       endif

       
      endif
      gall=gall+gb
      enddo

c      tt2=secnds(tt1)
      call cpu_time(t2)
      write(*,40) t2-t1 !,tt2
40    format(' Total time(seconds) = ', 1f5.2,/)
      write(*,*) 'RMS',rms
c---------------------------output重力数据--------------------------------c
      gall=facg*gall*fac

	call WriteGRD(fgname,m,n,gall
     +	,Xmin+cx0,Xmax+cx0,Ymin+cy0,Ymax+cy0,Zmin,Zmax)


 
c--------------------------------------------------------------- --c
       
10    FORMAT(a26,i10)
20    FORMAT(a26,f10.2)
30    FORMAT(a26,f10.2,f10.2,f10.2,f10.2,f10.2)      
      END PROGRAM
c******************************************************************c
c                                                                  c
c                  读取计算参数                                    c
c------------------------------------------------------------------c
cSUBROUTINE:
c   ReadCmd  
c   ReadG
c   writeG
c   writeresult
c------------------------------------------------------------------c
      SUBROUTINE ReadCmd(inname,fgrid,fbname,fgname,ilen,iunit
     +,infft,igtype,isalldata,isoutbou,m,rou0,cx0,cy0)
      	character*(*) inname,fgrid,fbname,fgname
		integer nunit_in,ilen,iunit,igtype,isalldata,isoutbou,m
	real rou0,cx0,cy0
    
        call search_unit(0,nunit_in)
        open(nunit_in,file=inname,status='old')
        call read_char(nunit_in,fgrid)
        call READ_integer(nunit_in,ilen)
        call READ_integer(nunit_in,iunit)
        call READ_integer(nunit_in,infft)
        call read_char(nunit_in,fbname)
        call READ_integer(nunit_in,igtype)
        call READ_integer(nunit_in,isalldata)
        call READ_integer(nunit_in,isoutbou)
        call read_char(nunit_in,fgname)
        close(nunit_in) 
        open(nunit_in,file=fbname,status='old')
        read(nunit_in,*) m,rou0
        read(nunit_in,*) cx0,cy0
        close(nunit_in) 
      END subroutine

      SUBROUTINE Readbodies(fbname,bodies,btype,m)
      character*(*) fbname
      dimension bodies(m,9),btype(m,2)
      real bodies,r1,r2,r3,r4,r5,r6,r7,r8,r9
      character*20 btype,c1,c2

        call search_unit(0,nunit_in)
        open(nunit_in,file=fbname,status='old')
        read(nunit_in,*)
        read(nunit_in,*) 
        ! 2.61,0.5,2.5,3.0,6.0,L,30,45,GL1
        do i=1,m
         read(nunit_in,*) r1,r2,r3,r4,r5,c1,r6,r7,c2,r8,r9
         bodies(i,1)=r1
         bodies(i,2)=r2
         bodies(i,3)=r3
         bodies(i,4)=r4
         bodies(i,5)=r5
         bodies(i,6)=r6
         bodies(i,7)=r7
         bodies(i,8)=r8
         bodies(i,9)=r9
         btype(i,1)=trim(c1)
         btype(i,2)=trim(c2)
        enddo
        close(nunit_in) 

      end subroutine 
c---------------------------------------------------------------c
      SUBROUTINE transg(g1,g2,m,n,id)
      dimension g1(m,n),g2(m,n) 
      do i=1,m
       do j=1,n
         if (id.eq.1) then    
          g2(i,j)=g1(m-i+1,j)
         else 
          g2(i,j)=g1(i,n-j+1)
         endif
       enddo
      enddo
      end subroutine 
c---------------------------------------------------------------c
      SUBROUTINE outbou(fbou1,width0,length,pointx,pointy,it,xita)
      character*(*) fbou1
      real width0,length,pointx,pointy,xita
         width=width0*sind(xita)
         if (abs(xita-90)>1.0E-30) then
         shift=width0*cosd(xita) !/tand(xita)
         else
         shift=0.0
         endif
         call search_unit(0,nunit_in)
        open(nunit_in,file=fbou1,status='unknown')
        write(nunit_in,*) 5,1 
        if (it.eq.1) then
	  write(nunit_in,*) -length-shift+pointx,width+pointy
	  write(nunit_in,*) length-shift+pointx,width+pointy
	  write(nunit_in,*) length+shift+pointx,-width+pointy
	  write(nunit_in,*) -length+shift+pointx,-width+pointy
	  write(nunit_in,*) -length-shift+pointx,width+pointy  
        else
	  write(nunit_in,*) -(-length-shift+pointx),width+pointy
	  write(nunit_in,*) -(length-shift+pointx),width+pointy
	  write(nunit_in,*) -(length+shift+pointx),-width+pointy
	  write(nunit_in,*) -(-length+shift+pointx),-width+pointy
	  write(nunit_in,*) -(-length-shift+pointx),width+pointy  
        endif
	  close(nunit_in)     
      
      end subroutine   
          
      SUBROUTINE outbou2(fbou1,width0,length,pointx,pointy,w1,it,xita)
      character*(*) fbou1
      real width0,length,pointx,pointy,w1,xita
         if (abs(xita-90)>1.0E-30) then
         st=width0*cosd(xita) !/tand(xita)
         else
         st=0.0
         endif
         width=width0*sind(xita)
         call search_unit(0,nunit_in)
        open(nunit_in,file=fbou1,status='unknown')
        write(nunit_in,*) 5,1 
        if (it.eq.1) then
	  write(nunit_in,*) -length-w1-st+pointx,width+pointy
	  write(nunit_in,*) length-st+pointx,width+pointy
	  write(nunit_in,*) length+st+pointx,-width+pointy
	  write(nunit_in,*) -length-w1+st+pointx,-width+pointy
	  write(nunit_in,*) -length-st-w1+pointx,width+pointy  
        else
	  write(nunit_in,*) -(-length-st-w1+pointx),width+pointy
	  write(nunit_in,*) -(length-st+pointx),width+pointy
	  write(nunit_in,*) -(length+st+pointx),-width+pointy
	  write(nunit_in,*) -(-length-w1+st+pointx),-width+pointy
	  write(nunit_in,*) -(-length-st-w1+pointx),width+pointy   
        endif
	  close(nunit_in)     
      
      end subroutine     
c------------------------------------------------------------------c

c      M----- Number of Points In X-Direction
c	 N----- Number of Points In Y-Direction
c------------------------------------------------------------------c
      SUBROUTINE ReadGRD(grdname,m,n,Xmin,Xmax,Ymin,Ymax,Zmin,Zmax)
      	character*(*) grdname
		integer nunit_in,m,n
		real Xmin,Xmax,Ymin,Ymax,Zmin,Zmax 
        
        call search_unit(0,nunit_in)
        open(nunit_in,file=grdname,status='old')
        read(nunit_in,*) !DSAA 
	  read(nunit_in,*) M,N
	  read(nunit_in,*) Xmin,Xmax
	  read(nunit_in,*) Ymin,Ymax
	  read(nunit_in,*) Zmin,Zmax
        close(nunit_in)  
      END subroutine
c-------------------------------------------------------------------c
      SUBROUTINE ReadGRDG(grdname,g,m,n,Xmin,Xmax,Ymin,Ymax,Zmin,Zmax)
      	character*(*) grdname
		integer nunit_in,m,n
	    dimension g(m,n)
		real Xmin,Xmax,Ymin,Ymax,Zmin,Zmax 
        
        call search_unit(0,nunit_in)
        open(nunit_in,file=grdname,status='old')
        read(nunit_in,*) !DSAA 
	  read(nunit_in,*) M,N
	  read(nunit_in,*) Xmin,Xmax
	  read(nunit_in,*) Ymin,Ymax
	  read(nunit_in,*) Zmin,Zmax
	  do j=1,n 
	    read(nunit_in,*) (g(i,j),i=1,m)
        enddo
        close(nunit_in)  
      END subroutine
c------------------------------------------------------------------c
      SUBROUTINE WriteGRD(grdname,m,n,g,Xmin,Xmax,Ymin,Ymax,Zmin,Zmax)
      	character*(*) grdname
		integer nunit_in,m,n
	    dimension g(m,n)
		real g,Xmin,Xmax,Ymin,Ymax,Zmin,Zmax
        
	zmin=g(1,1)
	zmax=g(1,1)
	  do i=1,m
	    do j=1,n
            if (zmin>g(i,j)) then
               zmin=g(i,j)
	      endif
	      if (zmax<g(i,j)) then
               zmax=g(i,j)
	      endif
	    enddo
	  enddo
        call search_unit(0,nunit_in)
        open(nunit_in,file=grdname,status='unknown')
        write(nunit_in,10) 'DSAA' 
	  write(nunit_in,*) M,N
	  write(nunit_in,*) Xmin,Xmax
	  write(nunit_in,*) Ymin,Ymax
	  write(nunit_in,*) Zmin,Zmax
	  do j=1,n 
	    write(nunit_in,*) (g(i,j),i=1,m)
        enddo
	  close(nunit_in)
10    format(a4)
      END subroutine

c******************************************************************c
c                                                                  c
c        forwar gravity anomaly a 3D model in spectral domain      c
c------------------------------------------------------------------c
c   g＝布格重力异常网格化数据
c   m----x点数 n----y点数
c   dg＝求导后的数据
c   h=向上延拓高度
c   method＝（0＝不求导、1＝Vzz、2＝Vzzz）
c   x,y 方向单位为 m
c   midu 单位为g/cm3
c   重力异常单位gravity_unit为g.u.
c------------------------------------------------------------------c
	SUBROUTINE CaluVz(gb,m,n,dx,dy,Xmin,Xmax,Ymin,Ymax,Zmin,Zmax
     +,a,b,h1,h2,alfa,xita,midu,pointx,pointy,rms,igt,ift,infft)

	integer m,n,MAXNM,me,ne,infft
	dimension gb(m,n)
	complex E,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp0,point,tm10
      real gb,Xmin,Xmax,Ymin,Ymax,zmin,zmax,a,b,h1,h2,alfa,xita,r1,r2,r3
	real,ALLOCATABLE::w(:,:),u(:,:),v(:,:),vzr(:,:),vzi(:,:)   
     +,TREAL(:),TIMAG(:),g(:,:),wc(:,:),uc(:,:),vc(:,:)
	real pi,bigG,midu,pointx,pointy,rms
      real t1,t2,t3,t11,t12,rxt               ! this is for time check  


        real,ALLOCATABLE:: density(:),x1(:),x2(:),
     +   y1(:),y2(:),z1(:),z2(:),dz(:)
        real,ALLOCATABLE:: x_point(:,:),y_point(:,:),
     +    z_point(:,:),field(:,:)
     
               
	pi=3.141592654
	bigG=6.672E-2
	write(*,*) 'Density difference = ',midu

      eigval=1.70141e+038
      EPS=1.E-6
 

      CALL CALCULATE_para(M,N,Me,Ne,infft)
c      Me=1024 !256
c      Ne=Me
      write(*,*) M,N,Me,Ne
	if (me>=ne) then
       MAXNM=me
	endif
	if (me<ne) then
	 MAXNM=ne
	endif
       
      allocate (v(me,ne),u(me,ne),w(me,ne),vc(me,ne),uc(me,ne),wc(me,ne)
     +,vzr(me,ne),vzi(me,ne),TREAL(MAXNM),TIMAG(MAXNM),g(me,ne))
      CALL CALCULATE_xyminmax(M,N,Me,Ne,Xmin,Xmax,Ymin,Ymax,
     +  XXmin,XXmax,YYmin,YYmax)     
     
     
      ! spatial foward for a rectangle prism model

      nbody=1      
      allocate(density(nbody),x1(nbody),x2(nbody),dz(nbody),
     +  y1(nbody),y2(nbody),z1(nbody),z2(nbody),field(me,ne))
      allocate(x_point(me,ne),y_point(me,ne),z_point(me,ne))
      
      X1(1)=-a
      X2(1)=a
      Y1(1)=-b
      Y2(1)=b
      Z1(1)=h1
      Z2(1)=h2
      density(1)=midu

      do i=1,me
	  do j=1,ne
	    x_point(i,j)=XXmin+(i-1)*dx
	    y_point(i,j)=YYmin+(j-1)*dy
	    z_point(i,j)=0.
	  enddo
      enddo
      call cpu_time(t11)
c      tts1=secnds(0.0)
        CALL RECTAN_gravity_sub(me,ne,x_point,y_point,z_point,
     +    field,nbody,X1,X2,Y1,Y2,Z1,Z2,dz,EPS,eigval,density)
      call cpu_time(t12)
c      tts2=secnds(tts1)
      write(*,20) t12-t11 !,tts2
20    format(' spatial forward time(seconds) = ', 1f12.8,/)
      
      call waveuvw2(w,u,v,me,ne,dx,dy)

	    if (abs(alfa-90)>1.0E-30) then
	       r1=1/tand(alfa)
	    else
	       r1=0
	    endif
	    if (abs(xita-90)>1.0E-30) then
             rxt=1/tand(xita)
	    else
	       rxt=0
	    endif
        	
	call cpu_time(t1)
c      tts1=secnds(0.0)	
	do i=1,me
	  do j=1,ne
          tmp4=exp(-w(i,j)*h2)
	    tmp5=exp(-w(i,j)*h1)
	    tmp6=exp(-u(i,j)*a*(0,1))
c          tmp6=cos(-u(i,j)*a)+(0,1)*sin(-u(i,j)*a)
	    tmp7=exp(u(i,j)*a*(0,1))
c         tmp7=cos(u(i,j)*a)+(0,1)*sin(u(i,j)*a)
	    tmp8=exp(u(i,j)*h2*r1*(0,1))
c	    tmp8=cos(u(i,j)*h2*r1)+(0,1)*sin(u(i,j)*h2*r1)
	    tmp9=exp(u(i,j)*h1*r1*(0,1))
c	    tmp9=cos(u(i,j)*h1*r1)+(0,1)*sin(u(i,j)*h1*r1)
      if (abs(w(i,j))>1.0E-30) then

        if (abs(u(i,j))>1.0E-30) then
c            tmp2=(tmp6*(tmp4-tmp5)/w(i,j)+tmp7*(tmp8*tmp4-tmp9*tmp5)
c     +/((0,1)*u(i,j)*r1-w(i,j)))/u(i,j)
             tmp2=(tmp6*(tmp4-tmp5)/w(i,j)+tmp7*(tmp8*tmp4-tmp9*tmp5)
     +/((0,1)*u(i,j)*r1-w(i,j))/tmp9)/u(i,j) !tmp9 new add
	    else
             tmp2=-2*(0,1)*a*(tmp4-tmp5)/w(i,j)
		endif	
        if (abs(v(i,j)-u(i,j)*rxt)>1.0E-30) then
c         write(*,*) v(i,j)
c         tmp3=-2*(0,1)*sin((v(i,j)-u(i,j)*rxt)*b)/(v(i,j)-u(i,j)*rxt) !new add xita
          tmp3=-2*(0,1)*sin((v(i,j)-u(i,j)*rxt)*b*sind(xita))
     +     /(v(i,j)-u(i,j)*rxt) !new add xita        
        else
c         tmp3=-2*(0,1)*b
         tmp3=-2*(0,1)*b*sind(xita)
        endif
          E=tmp2*tmp3		    
      else 
c         E=4*a*b*(h2-h1)+b*r1*(h2*h2-h1*h1)
c         E=4*a*b*(h2-h1)*sind(xita)+b*r1*(h2*h2-h1*h1) !new add xita
         E=4*a*b*(h2-h1)*sind(xita)+b*r1*(h2-h1)*(h2-h1)*sind(xita) !new add xita 20170929
c         write(*,*) E,xita,r1,b*r1*(h2*h2-h1*h1),b*r1*(h2-h1)*(h2-h1)
	endif
	  ! pointx=0.0
	tmp1=exp(((Xmax+Xmin)*u(i,j)/2+(Ymax+Ymin)*v(i,j)/2)*(0,1))
	tmp0=exp(((XXmin-XXmax)*u(i,j)/2+(YYmin-YYmax)*v(i,j)/2)*(0,1))
c	tmp0=exp(((XXmin-XXmax+dx)*u(i,j)/2
c     +	+(YYmin-YYmax+dy)*v(i,j)/2)*(0,1))
	point=exp((0,-1)*(pointx*u(i,j)+pointy*v(i,j)))

          E=E*tmp0*tmp1*point/me/ne/dx/dy
          if (igt.eq.1) then
           E=E*w(i,j)
          else if (igt.eq.2) then
           E=E*u(i,j)*(0,1)
          endif
		vzr(i,j)=real(E)
		vzi(i,j)=imag(E)
	    g(i,j)=sqrt(real(E)*real(E)+imag(E)*imag(E))
          
	  enddo
	enddo
	
      call cpu_time(t2)
c      tts2=secnds(tts1)
      write(*,30) t2-t1 !,tts2
30    format(' spectral forward time(seconds) = ', 1f12.8,/)

	call FFT2(vzr,vzi,me,ne,1,TREAL,TIMAG,MAXNM)
      
      call cpu_time(t3)
      write(*,40) t3-t2
c      write(*,*) minval(vzi),maxval(vzi)
40    format(' iFFT time(seconds) = ', f12.8,/)	
      rmstmp=0.0
	rvalmax=-1.E10
            
	do i=1,m
	 do j=1,n
	  rmstmp1=2*pi*bigG*midu*vzr(i+(me-m)/2,j+(ne-n)/2)
	  rmstmp2=field(i+(me-m)/2,j+(ne-n)/2)
c	  gb(i,j)=2*pi*bigG*midu*vzr(i+(me-m)/2,j+(ne-n)/2)
	  gb(i,j)=rmstmp1 !rmstmp2 !output forward model
	  rmstmp=rmstmp+(rmstmp1-rmstmp2)**2
	   if (rmstmp2>rvalmax) then
	     rvalmax=rmstmp2
	   endif
       enddo
	enddo
	
	rms=sqrt(rmstmp/m/n)
	write(*,*) 'Max,%',rvalmax,rms*100/rvalmax
      END subroutine

c******************************************************************c
c                                                                  c
c             extend boundary using cos funtion                    c
c------------------------------------------------------------------c
c    CALCULATE_para -- get the count for expend (from Mpoint to M) c
c    FFT_expand2-----  cos expend subroutine                       c
c      i.e.
c       (0, 1, ......, n/2-1, n/2, -n/2+1, ......, -1)             c
c------------------------------------------------------------------c	 	
      SUBROUTINE CALCULATE_para(Mpoint,line,M,N,infft)
	integer Mpoint,line,M,N,infft
     
        K_Mpoint=int(LOG(real(Mpoint))/LOG(2.0)+0.1)
	  K_line=int(LOG(real(line))/LOG(2.0)+0.1)
        M=2**(K_Mpoint+infft)       
        N=2**(K_line+infft)

      END SUBROUTINE  
c******************************************************************c
c                Calculate the XXmin,XXmax,YYmin,YYmax  
c------------------------------------------------------------------c
      SUBROUTINE CALCULATE_xyminmax(Mpoint,line,M,N,Xmin,Xmax,Ymin,Ymax,
     +  XXmin,XXmax,YYmin,YYmax) 
      integer Mpoint,line,M,N
	real Xmin,Xmax,Ymin,Ymax,XXmin,XXmax,YYmin,YYmax
	
	dx=(Xmax-Xmin)/(Mpoint-1)
	dy=(Ymax-Ymin)/(line-1)
      M_left=(M-Mpoint)/2+1
	M_right=M_left+Mpoint-1
	N_down=(N-line)/2+1
	N_up=N_down+line-1
      XXmin=Xmin-dx*(M_left-1)
	XXmax=Xmax+dx*(M-M_right)
	YYmin=Ymin-dy*(N_down-1)
	YYmax=Ymax+dy*(N-N_up)
	
     
	END SUBROUTINE
c******************************************************************c
c                                                                  c
c                        wave number only output w                 c
c------------------------------------------------------------------c
c  w----array for output the wave number
c  m----count of x direction  n----count of y direction  
c  dx---x distance  dy--- y distance
c------------------------------------------------------------------c
      SUBROUTINE wave2(w,m,n,dx,dy)
      integer m,n
	dimension w(m,n)
	real w

	pi=3.141592654
	msh=m/2+1
	nsh=n/2+1
	sx=pi*2/float(m)/dx
	sy=pi*2/float(n)/dy
	do j=1,n
	  kky=j
	   if (kky .gt. nsh) then
	     kky=kky-n 
	   endif
         wy2=(sy*(kky-1))*(sy*(kky-1))
	  do i=1,m
	  kkx=i
	   if (kkx .gt. msh) then
	     kkx=kkx-m 
	   endif
	     wx2=(sx*(kkx-1))*(sx*(kkx-1))           
	     w(i,j)=sqrt(wx2+wy2)
	  enddo
	enddo
	 
      END subroutine
c******************************************************************c
c                                                                  c
c                        Wave number with angle                    c
c------------------------------------------------------------------c
c  w,u,v----array for output the wave number
c  m----count of x direction  n----count of y direction  
c  dx---x distance  dy--- y distance
c------------------------------------------------------------------c
      SUBROUTINE wave2lm(w,m,n,dx,dy,a,b)
      integer m,n
	dimension w(m,n)
	real w,a,b

	pi=3.141592654
	msh=m/2+1
	nsh=n/2+1
	sx=pi*2/float(m)/dx
	sy=pi*2/float(n)/dy
	do j=1,n
	  kky=j
	   if (kky .gt. nsh) then
	     kky=kky-n 
	   endif
         wy2=sy*(kky-1)
	  do i=1,m
	  kkx=i
	   if (kkx .gt. msh) then
	     kkx=kkx-m 
	   endif
	     wx2=sx*(kkx-1)
	     w(i,j)=cosd(a)*wx2+cosd(b)*wy2
	  enddo
	enddo
	 
      END subroutine
c******************************************************************c
c                                                                  c
c                        Wave number (w,u,v)                       c
c------------------------------------------------------------------c
c  w,u,v----array for output the wave number
c  m----count of x direction  n----count of y direction  
c  dx---x distance  dy--- y distance
c------------------------------------------------------------------c
      SUBROUTINE waveuvw2(w,u,v,m,n,dx,dy)
      integer m,n
	dimension w(m,n),u(m,n),v(m,n)
	real w,u,v

	pi=3.141592654
	msh=m/2+1
	nsh=n/2+1
	sx=pi*2/float(m)/dx
	sy=pi*2/float(n)/dy
	do j=1,n
	  kky=j
	   if (kky .gt. nsh) then
	     kky=kky-n 
	   endif
         wy2=(sy*(kky-1))*(sy*(kky-1))
	  do i=1,m
	  kkx=i
	   if (kkx .gt. msh) then
	     kkx=kkx-m 
	   endif
	     wx2=(sx*(kkx-1))*(sx*(kkx-1))           
	     w(i,j)=sqrt(wx2+wy2)
		 u(i,j)=sx*(kkx-1)
		 v(i,j)=sy*(kky-1)
	  enddo
	enddo

      END subroutine
!******************************************************************c
!  FUNCTION:                                                       c
!     SEARCH the unopened UNIT                                     c
!                                                                  c
!  INPUT:                                                          c
!    NUNIT_START: the initial UNIT                                 c
!                                                                  c
!  OUTPUT:                                                         c
!    NUNIT_IN: searched UNIT                                       c
!------------------------------------------------------------------c
        SUBROUTINE SEARCH_UNIT(NUNIT_START,NUNIT_IN)
	    logical unit_open
        integer nunit_in,nunit_start

        nunit_in=nunit_start
     	unit_open=.TRUE.
    	do 10 while (unit_open)
        nunit_in=nunit_in+1
        inquire(unit=nunit_in,opened=unit_open)
10      continue

        END SUBROUTINE

!******************************************************************c
!  FUNCTION:                                                       c
!     READ an integer number from NUNIT_IN                         c
!                                                                  c
!  OUTPUT:                                                         c
!   M:  the integer number                                         c
!------------------------------------------------------------------c
        SUBROUTINE READ_INTEGER(NUNIT_IN,Mtemp)
        character*200 ctemp,rr1

        indexposi=0
        DO WHILE (indexposi.eq.0)
        read(nunit_in,'(a)')ctemp
        indexposi=INDEX(ctemp,'=')
        END DO  !while
        read(ctemp,'(<indexposi>X,A)')rr1
        rr1=TRIM(rr1)
        read(rr1,*)Mtemp
        
        RETURN
        END SUBROUTINE

!******************************************************************c
!  FUNCTION:                                                       c
!     READ three real number from NUNIT_IN                         c
!                                                                  c
!  OUTPUT:                                                         c
!   R1,R2,R3:  the real number                                     c
!------------------------------------------------------------------c
        SUBROUTINE READ_REAL(NUNIT_IN,r1)
        character*200 ctemp,rr1
	  real r1
        indexposi=0
        DO WHILE (indexposi.eq.0)
        read(nunit_in,'(a)')ctemp
        indexposi=INDEX(ctemp,'=')
        END DO  !while
        read(ctemp,'(<indexposi>X,A)')rr1
        rr1=TRIM(rr1)
        read(rr1,*)r1
        RETURN
        END SUBROUTINE
!******************************************************************c
!  FUNCTION:                                                       c
!     READ three real number from NUNIT_IN                         c
!                                                                  c
!  OUTPUT:                                                         c
!     R1,R2,R3:  the real number                                   c
!------------------------------------------------------------------c
        SUBROUTINE READ_REAL_2(NUNIT_IN,r1,r2)
        character*200 ctemp,rr1
	  real r1,r2
        indexposi=0
        DO WHILE (indexposi.eq.0)
        read(nunit_in,'(a)')ctemp
        indexposi=INDEX(ctemp,'=')
        END DO  !while
        read(ctemp,'(<indexposi>X,A)')rr1
        rr1=TRIM(rr1)
        read(rr1,*)r1,r2
        RETURN
        END SUBROUTINE


!******************************************************************c
!  FUNCTION:                                                       c
!     READ the character string from NUNIT_IN                      c
!                                                                  c
!  OUTPUT:                                                         c
!     CHAR: the character string                                   c
!------------------------------------------------------------------c
        SUBROUTINE READ_CHAR(NUNIT_IN,CHAR)
        character*(*) char
        character*200 ctemp

        indexposi=0
        DO WHILE (indexposi.eq.0)
        read(nunit_in,'(a)')ctemp
        indexposi=INDEX(ctemp,'=')
        END DO  !while
        read(ctemp,'(<indexposi>X,A)')char
        char=TRIM(char)
        
        END SUBROUTINE
!******************************************************************c
!  FUNCTION:                                                       c
!     read a file name for cmd                                     c
!                                                                  c
!  OUTPUT:                                                         c
!   CMDFILE: filename                                              c
!------------------------------------------------------------------c
        SUBROUTINE readparams(CMDFILE)
        USE MSFLIB
        character*(*) cmdfile
        integer*2 number_arg,number

        if(nargs().ne.2) then
        STOP 'You are illegal User !'
        end if
        number_arg=1
        call getarg(number_arg,cmdfile,number)

        RETURN
        END SUBROUTINE

c******************************************************************c
c                                                                  c
c        2-D  FAST FOURIER TRANSFORM FOR REAL FUNCTION             c
c------------------------------------------------------------------c
c    rreal -- Real part of the function to be transformed          c
c    rimag -- Imaginary part of the function to be transformed     c
c    M ------ The number of points. M must be the mu-th power of 2 c
c    N ------ The number of points. N must be the nu-th power of 2 c
c    MAXNM -- The maximum of M and N                               c
c    NF ----- (=1,reverse transform                                c
c             (=2, normal transform                                c
c                                                                  c
c    TREAL -- Real part of the working array                       c
c    TIMAG -- Imaginary part of the working array                  c
c------------------------------------------------------------------c
	SUBROUTINE FFTA2(rreal,rimag,M,N,NF,TREAL,TIMAG,MAXNM)
	dimension rreal(m,n),rimag(m,n),treal(maxnm),timag(maxnm)

	nmmax=max(m,n)
	if(min(m,n).gt.1) goto 90
	call ffta(rreal,rimag,nmmax,nf)
        RETURN
90      continue
	if(nf.ne.1) then
	do 100 i=1,m
	do 105 j=1,n
	treal(j)=rreal(i,j)
105     timag(j)=rimag(i,j)
	call ffta(treal,timag,n,nf)
	do 100 j=1,n
	rreal(i,j)=treal(j)
100     rimag(i,j)=timag(j)
	do 200 j=1,n/2+1
	jj=n-j+2
	if(j.eq.1) jj=j
	do 110 i=1,m
110     treal(i)=rreal(i,j)
	call ffta(treal,timag,m,nf)
	do 120 i=1,m
	rreal(i,j)=treal(i)
	treal(i)=rimag(i,j)
	rimag(i,j)=timag(i)
	rreal(i,jj)=rreal(i,j)
	rimag(i,jj)=rimag(i,j)
120     continue
	call ffta(treal,timag,m,nf)
	do 130 i=1,m
	rreal(i,j)=rreal(i,j)-timag(i)
	rreal(i,jj)=rreal(i,jj)+timag(i)
	rimag(i,j)=rimag(i,j)+treal(i)
	rimag(i,jj)=rimag(i,jj)-treal(i)
130     continue
200     continue
	else 
	do 300 j=1,n/2+1
	jj=n-j+2
	if(j.eq.1) jj=j
	do 310 i=1,m
	treal(i)=(rreal(i,j)+rreal(i,jj))/2.0
310     timag(i)=(rimag(i,j)+rimag(i,jj))/2.0
	call ffta(treal,timag,m,nf)
	do 320 i=1,m
	dd=treal(i)
	treal(i)=(rimag(i,j)-rimag(i,jj))/2.0
	timag(i)=-(rreal(i,j)-rreal(i,jj))/2.0
320     rreal(i,j)=dd
	call ffta(treal,timag,m,nf)
	do 330 i=1,m
330     rimag(i,j)=treal(i)
300     continue
	do 340 i=1,m
	do 350 j=1,n/2+1
	jj=n-j+2
	if(j.eq.1) jj=j
	treal(j)=rreal(i,j)
	timag(j)=rimag(i,j)
	treal(jj)=treal(j)
350     timag(jj)=-timag(j)
	call ffta(treal,timag,n,nf)
	do 360 j=1,n
	rreal(i,j)=treal(j)
360     rimag(i,j)=timag(j)
340     continue
	end if

        RETURN  
        END SUBROUTINE


c******************************************************************c
c                                                                  c
c         1-D FAST FOURIER TRANSFORM  FOR REAL FUNCTION            c
c------------------------------------------------------------------c
c    XREAL -- Real part of the function to be transformed          c
c    XIMAG -- Imaginary part of the function to be transformed     c
c    N ------ The number of points. N must be the nu-th power of 2 c
c    NF ----- (=1,reverse transform                                c
c             (=2, normal transform                                c
c       (0, 1, ......, n/2-1, n/2, -n/2+1, ......, -1)             c
c------------------------------------------------------------------c
	SUBROUTINE FFTA(XREAL,XIMAG,N,NF)
	dimension xreal(n),ximag(n)

	n2=n/2.0
	arg2=3.141592653589/float(n2)
	if(nf.eq.2) then
	do 10 i=1,n2
	xreal(i)=xreal(2*i-1)
 10     ximag(i)=xreal(2*i)
	else
	arg=0.0
	do 30 i=1,n2
	c=cos(arg)
	s=sin(arg)
	xre=xreal(i)-xreal(n2+i)
	xim=ximag(i)-ximag(n2+i)
	xreal(i)=0.5*(xreal(i)+xreal(n2+i)-xre*s-xim*c)
	ximag(i)=0.5*(ximag(i)+ximag(n2+i)+xre*c-xim*s)
	arg=arg+arg2
 30     continue
	end if
	call fft(xreal,ximag,n2,nf)
	if(nf.eq.2) then
	arg=0.0
	do 50 k=1,n2
	c=cos(arg)
	s=sin(arg)
	kk=n2-k+2
	if(k.eq.1) kk=k
	xrak=xreal(k)+xreal(kk)
	xiak=ximag(k)-ximag(kk)
	xrbk=ximag(k)+ximag(kk)
	xibk=-xreal(k)+xreal(kk)
	xreal(n2+k)=(xrak-xrbk*c-xibk*s)*0.5
	ximag(n2+k)=(xiak-xibk*c+xrbk*s)*0.5
 50     arg=arg+arg2
	xreal(1)=xreal(1)+ximag(1)
	ximag(1)=0.0
	do 55 k=2,n2
	xreal(k)=xreal(n-k+2)
 55     ximag(k)=-ximag(n-k+2)
	else
	do 40 i=n2,1,-1
	xreal(2*i-1)=xreal(i)
	xreal(2*i)=ximag(i)
	ximag(i)=0.0
	ximag(n2+i)=0.0
 40     continue
        end if

        RETURN
        END SUBROUTINE



c******************************************************************c
c                                                                  c
c        2-D  FAST FOURIER TRANSFORM FOR COMPLEX FUNCTION          c
c------------------------------------------------------------------c
c    rreal -- Real part of the function to be transformed          c
c    rimag -- Imaginary part of the function to be transformed     c
c    M ------ The number of points. M must be the mu-th power of 2 c
c    N ------ The number of points. N must be the nu-th power of 2 c
c    MAXNM -- The maximum of M and N                               c
c    NF ----- (=1,reverse transform                                c
c             (=2, normal transform                                c
c                                                                  c
c    TREAL -- Real part of the working array                       c
c    TIMAG -- Imaginary part of the working array                  c
c------------------------------------------------------------------c
	SUBROUTINE FFT2(rreal,rimag,M,N,NF,TREAL,TIMAG,MAXNM)
        real rreal(m,n),rimag(m,n),treal(maxnm),timag(maxnm)
	do 100 i=1,m
	if(n.eq.1) goto 210
	do 80 j=1,n
	treal(j)=rreal(i,j)
  80    timag(j)=rimag(i,j)
	call fft(treal,timag,n,nf)
	do 100 j=1,n
	rreal(i,j)=treal(j)
100     rimag(i,j)=timag(j)
210     do 200 j=1,n
        if(m.eq.1) RETURN
	do 180 i=1,m
	treal(i)=rreal(i,j)
180     timag(i)=rimag(i,j)
	call fft(treal,timag,m,nf)
	do 200 i=1,m
	rreal(i,j)=treal(i)
200     rimag(i,j)=timag(i)

        RETURN  
        END SUBROUTINE



c******************************************************************c
c                                                                  c
c        1-D FAST FOURIER TRANSFORM FOR COMPLEX FUNCTION           c
c------------------------------------------------------------------c
c    XREAL -- Real part of the function to be transformed          c
c    XIMAG -- Imaginary part of the function to be transformed     c
c    N ------ The number of points. N must be the nu-th power of 2 c
c    NF ----- (=1,reverse transform                                c
c             (=2, normal transform                                c
c         (0, 1, ......, n/2-1, n/2, -n/2+1, ......,-1)            c
c------------------------------------------------------------------c
	SUBROUTINE FFT(XREAL,XIMAG,N,NF)
        real xreal(n),ximag(n)
        nu=int(log(float(n))/0.693147+0.001)
	n2=n/2
	nu1=nu-1
	f=float((-1)**nf)
	k=0
	do 100 l=1,nu
102     do 101 i=1,n2
	p=ibitr(k/2**nu1,nu)
	arg=6.2831853*p*f/float(n)
	c=cos(arg)
	s=sin(arg)
	k1=k+1
	k1n2=k1+n2
	treal=xreal(k1n2)*c+ximag(k1n2)*s
	timag=ximag(k1n2)*c-xreal(k1n2)*s
	xreal(k1n2)=xreal(k1)-treal
	ximag(k1n2)=ximag(k1)-timag
	xreal(k1)=xreal(k1)+treal
	ximag(k1)=ximag(k1)+timag
101     k=k+1
	k=k+n2
	if(k.lt.n) goto 102
	k=0
	nu1=nu1-1
100     n2=n2/2
	do 103 k=1,n
	i=ibitr(k-1,nu)+1
	if(i.le.k) goto 103
	treal=xreal(k)
	timag=ximag(k)
	xreal(k)=xreal(i)
	ximag(k)=ximag(i)
	xreal(i)=treal
	ximag(i)=timag
103     continue
        if(nf.eq.1) RETURN
	do 104 i=1,n
	xreal(i)=xreal(i)/float(n)
104     ximag(i)=ximag(i)/float(n)

        RETURN
        END SUBROUTINE
c
        FUNCTION IBITR(J,NU)
	j1=j
	itt=0
	do 20 i=1,nu
	j2=j1/2
	itt=itt*2+(j1-2*j2)
	ibitr=itt
  20    j1=j2
        return
        end function
c------------------------------------------------------------------c
	SUBROUTINE FFT_expand(Mpoint,line,M,N,g_temp,g_expand,
     +  method_expand)
      real g_temp(Mpoint,line),g_expand(M,N)
      integer Mpoint,line,M,N,method_expand
      DOUBLE PRECISION Pi
	 
	M_left=(M-Mpoint)/2+1
	M_right=M_left+Mpoint-1
	N_down=(N-line)/2+1
	N_up=N_down+line-1
      
  	!central
	average=0.0
      do j=1,line
	   do i=1,Mpoint 
	      g_expand(M_left+i-1,N_down+j-1)=g_temp(i,j)
	      average=average+g_temp(i,j)
	   end do
	end do
	average=average/(Mpoint*line)  
	!set value at left and right
      DO j=N_down,N_up
	  g_expand(1,j)=expand_side(method_expand,average,
     +                   g_expand(M_left,j),g_expand(M_right,j))
	  g_expand(M,j)=g_expand(1,j)
	
	!Left
	  do i=2,M_left-1
	    g_expand(i,j)=expand_sub(g_expand(1,j),g_expand(M_left,j),
     +		i,1,M_left)
	  end do
        
	
      !Right 
	  do i=M_right+1,M-1
	   g_expand(i,j)=expand_sub(g_expand(M,j),g_expand(M_right,j),
     +      i,M,M_right)               
        end do
        
      END DO
	
	!set value at up and down
	DO i=1,M
	  g_expand(i,N)=expand_side(method_expand,average,
     +                   g_expand(i,N_down),g_expand(i,N_up))
	  g_expand(i,1)=g_expand(i,N)
      !up
	  do j=N_up+1,N-1
	    g_expand(i,j)=expand_sub(g_expand(i,N),g_expand(i,N_up),
     +       j,N,N_up)               
        end do
	!down 
	  do j=2,N_down-1
	   g_expand(i,j)=expand_sub(g_expand(i,1),g_expand(i,N_down),
     +     j,1,N_down)                
        end do
	END DO

	END SUBROUTINE

c------------------------------------------------------------------c
      FUNCTION expand_sub(left_side,right_side,i,number,M_left)
	real expand_sub
	real left_side,right_side
	integer i,number,M_left
      
	Pii=3.141592654/2.0
      expand_sub=left_side+(right_side-left_side)*cos
     +  (Pii*(i-M_left)/(number-M_left))
	
	END FUNCTION
c------------------------------------------------------------------c
      FUNCTION expand_side(method_expand,average,left_side,right_side)
	real expand_side
	integer method_expand
	real average,left_side,right_side

	if(method_expand.eq.1)then
	  expand_side=0.0
      else if(method_expand.eq.2)then
	  expand_side=average
	else if(method_expand.eq.3)then
	  expand_side=(left_side+right_side)/2
      else
	  write(*,*)'INPUT ERROR!!'
      end if
	
	END FUNCTION


c******************************************************************c
c
c  The subroutine of CALCULATE THE gravity FIELDS OF RECTANGULAR BODIES
c
c******************************************************************c
        SUBROUTINE RECTAN_gravity_sub(mpoint,line,xx,yy,zz,field,
     +    number_body,X1,X2,Y1,Y2,Z1,Z2,dz,EPS,eigval,density)
        REAl dz(number_body)
        REAl x1(number_body),x2(number_body),y1(number_body),
     +    y2(number_body),z1(number_body),z2(number_body)
        REAl density(number_body)
        REAl xx(mpoint,line),yy(mpoint,line),zz(mpoint,line)
        REAL field(mpoint,line)

        REAl u(2),v(2),w(2),ru(2,2,2),rv(2,2,2),rw(2,2,2)
        REAl ra(2,2,2),rb(2,2,2),rc(2,2,2)


          factor=6.672*0.01

          do m=1,number_body
            dz(m)=factor*density(m)
          end do
        
        DO j=1,line
          do i=1,mpoint
            IF (zz(i,j).lt.eigval/2.0) THEN
            field(i,j)=0.0
            do m=1,number_body
              u(1)=x1(m)-xx(i,j)
              u(2)=x2(m)-xx(i,j)
              v(1)=y1(m)-yy(i,j)
              v(2)=y2(m)-yy(i,j)
              w(1)=z1(m)-zz(i,j)
              w(2)=z2(m)-zz(i,j)
              CALL r_rectan(u,v,w,ru,rv,rw,ra,rb,rc,2)
                CALL vz_rectan(u,v,w,ru,vz,dz(m),2,eps)
                field(i,j)=field(i,j)+vz*dz(m)
            end do
            ELSE
            field(i,j)=eigval
            END IF
          end do
        end do

        END subroutine
c******************************************************************c
c                                                                  c
c                         R=SQRT(U*U+V*V+W*W)                      c
c                                                                  c
c******************************************************************c
         SUBROUTINE R_rectan(U,V,W,RU,RV,RW,RA,RB,RC,L)
         DIMENSION u(l),v(l),w(l),ru(l,l,l),rv(l,l,l),rw(l,l,l)
         DIMENSION ra(l,l,l),rb(l,l,l),rc(l,l,l)
         DIMENSION u2(2),v2(2),w2(2)

         do i=1,l
           u2(i)=u(i)*u(i)
           v2(i)=v(i)*v(i)
           w2(i)=w(i)*w(i)
         end do

         do i=1,l
           do j=1,l
             uv2=u2(i)+v2(j)
             do k=1,l
               ru(i,j,k)=sqrt(uv2+w2(k))
               rv(j,i,k)=ru(i,j,k)
               rw(k,j,i)=ru(i,j,k)
               ra(j,k,i)=ru(i,j,k)
               rb(k,i,j)=ru(i,j,k)
               rc(i,k,j)=ru(i,j,k)
             end do
           end do
         end do
         
        END SUBROUTINE



c******************************************************************c
c                                         U*V  |u2|v2|w2           c
c      VZ=-[U*Ln(V+R)+V*Ln(U+R)-W*arctg(-----)]|  |  |             c
c                                         W*R  |u1|v1|w1           c
c------------------------------------------------------------------c
        SUBROUTINE VZ_rectan(U,V,W,R,VZ,CC,L,eps)
        DIMENSION u(l),v(l),w(l),r(l,l,l)

        vz=0.0
        IF(ABS(cc).gt.eps) THEN
        do i=1,l
          do j=1,l
            uv=u(i)*v(j)
            do k=1,l
              ii=(-1)**(i+j+k)
              if(w(k).eq.0.0) then
                ptan=0.0
              else
                ptan=atan(uv/w(k)/r(i,j,k))*w(k)
              end if
              if(v(j)+r(i,j,k).eq.0.0) then
                uvr=0.0
              else
                uvr=u(i)*log(v(j)+r(i,j,k))
              end if
              if(u(i)+r(i,j,k).eq.0.0) then
                vur=0.0
              else
                vur=v(j)*log(u(i)+r(i,j,k))
              end if
              vz=vz+ii*(uvr+vur-ptan)
            end do
          end do
        end do
	vz=-vz
        END IF

        END SUBROUTINE

      program pproc

      implicit none

      integer::ncells,nedge,nvert,i,boundmax
      integer,allocatable,dimension(:,:)::edge,bound,edgept
      integer,allocatable,dimension(:)::nbound,blanking,nmeet,edgeList
      double precision,allocatable,dimension(:,:)::grd
      double precision,allocatable,dimension(:)::rho,rhov,mach,machv,cp,cpv,u,v,    &
     &  uv,vv
      integer,allocatable,dimension(:)::m

      integer::point,nmeetmax,ed,j,nquad,no1,no2,nc,bd,fileno
      integer::possiblenextpt,possiblenexted,nv,unitno,unitno2,np,      &
	 &  nstep,t,ii,jj,nLines,lineLength,lStart,lEnd,ordV,ordE,ordC,     &
	 &  nEdgeBound

      double precision::x

      double precision::rhol,rhor,rhof,ul,ur,uf,machl,machr,machf,cpl,cpr,cpf,vl,vr &
     &  ,vf
	CHARACTER(LEN=15):: formV,formE,formC,strT
      boundmax=100

      nmeetmax=1000

      print *,'Steps:'
      read *,nstep

      open(100,file='griduns',status='unknown')
      open(101,file='flow.dat',status='unknown')
      open(303,file='coords.dat',status='unknown')
      read(100,*) ncells,nedge,nvert
		
      allocate(edge(nedge,4))
      allocate(grd(nvert,2))
      allocate(bound(ncells,boundmax))
      allocate(nbound(ncells))
      allocate(blanking(nedge))
      allocate(edgept(nvert,nmeetmax))
      allocate(nmeet(nvert))
      allocate(cp(ncells))
      allocate(cpv(nvert))
      allocate(mach(ncells))
      allocate(machv(nvert))
      allocate(rho(ncells))
      allocate(rhov(nvert))
      allocate(m(nvert))
      allocate(u(ncells))
      allocate(uv(nvert))
      allocate(v(ncells))
      allocate(vv(nvert))
	  allocate(edgeList(nedge))

      do i=1,nedge
        read(100,*) (edge(i,j),j=1,4)
      enddo

      unitno=222
      open(unitno,file="flowplt.plt",status='unknown')
      unitno2=223
      open(unitno2,file="flowplt_cell.plt",status='unknown')
      write(unitno,*) 'variables = "x" "y" "rho" "u" "v" "M" "cp" '
      write(unitno2,*) 'variables = "x" "y" "rho" "u" "v" "M" "cp" '

      do t=1,nstep

		do nv=1,nvert
		read(303,*) (grd(nv,j),j=1,2)
	  enddo
	  !close(100)

	  do i=1,ncells
		read(101,*) rho(i),u(i),v(i),mach(i),cp(i)
	  enddo
	  !close(101)

	 

	  cpv=0.0
	  uv=0.0
	  vv=0.0
	  machv=0.0
	  rhov=0.0
	  m=0.0

	  do i=1,nedge

		m(edge(i,1))=m(edge(i,1))+1
		m(edge(i,2))=m(edge(i,2))+1

		if(edge(i,3).gt.0)then
		  cpl=cp(edge(i,3))         
		else
		  cpl=cp(edge(i,4))
		endif
		
		if(edge(i,4).gt.0)then
		  cpr=cp(edge(i,4))
		else
		  cpr=cp(edge(i,3))
		endif  
		cpf=0.5*(cpl+cpr)

		if(edge(i,3).gt.0)then
		  machl=mach(edge(i,3))
		else
		  machl=mach(edge(i,4))
		endif
		if(edge(i,4).gt.0)then
		  machr=mach(edge(i,4))
		else
		  machr=mach(edge(i,3))
		endif  
		machf=0.5*(machl+machr)

		if(edge(i,3).gt.0)then
		  rhol=rho(edge(i,3))
		else
		  rhol=rho(edge(i,4))
		endif
		if(edge(i,4).gt.0)then
		  rhor=rho(edge(i,4))
		else
		  rhor=rho(edge(i,3))
		endif  
		rhof=0.5*(rhol+rhor)

		if(edge(i,3).gt.0)then
		  ul=u(edge(i,3))
		else
		  ul=u(edge(i,4))
		endif
		if(edge(i,4).gt.0)then
		  ur=u(edge(i,4))
		else
		  ur=u(edge(i,3))
		endif  
		uf=0.5*(ul+ur)

		if(edge(i,3).gt.0)then
		  vl=v(edge(i,3))
		else
		  vl=v(edge(i,4))
		endif
		if(edge(i,4).gt.0)then
		  vr=v(edge(i,4))
		else
		  vr=v(edge(i,3))
		endif  
		vf=0.5*(vl+vr)


		  cpv(edge(i,1))=cpv(edge(i,1))+cpf
		  cpv(edge(i,2))=cpv(edge(i,2))+cpf

		  machv(edge(i,1))=machv(edge(i,1))+machf
		  machv(edge(i,2))=machv(edge(i,2))+machf

		  rhov(edge(i,1))=rhov(edge(i,1))+rhof
		  rhov(edge(i,2))=rhov(edge(i,2))+rhof

		  uv(edge(i,1))=uv(edge(i,1))+uf
		  uv(edge(i,2))=uv(edge(i,2))+uf

		  vv(edge(i,1))=vv(edge(i,1))+vf
		  vv(edge(i,2))=vv(edge(i,2))+vf
		
	  enddo
	 
	  
	  
    do i=1,nvert
		!if(i.eq.381) then
		!  print *,m(i),rhov(i)
		!endif
		if(m(i).gt.0) cpv(i)=cpv(i)/float(m(i))
		if(m(i).gt.0) uv(i)=uv(i)/float(m(i))
		if(m(i).gt.0) vv(i)=vv(i)/float(m(i))
		if(m(i).gt.0) machv(i)=machv(i)/float(m(i))
		if(m(i).gt.0) rhov(i)=rhov(i)/float(m(i))
      enddo
	    
      if(t==1)then
	      jj=0;
	      do ed=1,nedge
		    if(edge(ed,3).eq.(-1))then
			    jj=jj+1;
			    edgeList(jj)=ed
		    endif
		    if(edge(ed,4).eq.(-1))then
			    jj=jj+1;
			    edgeList(jj)=ed
		    endif
          enddo
      endif
		
		
		
	  write(strT,*) t
	  write(unitno,*) 'ZONE t="1" N=',nvert,'E=',nedge
	  write(unitno,*) 'STRANDID=1 '
	  write(unitno,*) 'SOLUTIONTIME=',trim(adjustl(strT))
	  write(unitno,*) 'DATAPACKING=POINT,'
	  write(unitno,*) 'ZONETYPE=FELINESEG'

	  do np=1,nvert
		write(unitno,*) (grd(np,j),j=1,2),rhov(np),uv(np),vv(np),       &
	 &  machv(np),cpv(np)
	  enddo
	  
	  do ed=1,nedge
		  write(unitno,*) (edge(ed,j),j=1,2)
	  enddo
	  
	  !close(unitno)
	  
	  
	  ordV=ceiling(log10(dble(nvert)))
	  ordE=ceiling(log10(dble(nedge)))
	  ordC=ceiling(log10(dble(ncells)))
	  
	 
	 write(formV,"(I15)") nvert
	 formV = adjustl(formV)
	 write(formE,"(I15)") nedge
	 formE = adjustl(formE)
	 write(formC,"(I15)") ncells
	 formC = adjustl(formC)
	 
	  do i=1,nedge
		if(edge(i,3).gt.0)then
		  edge(i,3)=(edge(i,3))         
		else
		  edge(i,3)=0
		endif
		if(edge(i,4).gt.0)then
		 edge(i,4)=(edge(i,4))
	   else
		 edge(i,4)=0
		endif  
	  enddo
	  
	  !write(unitno2,*) VARLOCATION=([1-2]=NODAL ,[3-5]=CELLCENTERED)
	  write(unitno2,*) 'ZONE '
	  write(unitno2,*) 'STRANDID=1 '
	  write(unitno2,*) 'SOLUTIONTIME=', trim(adjustl(strT))
	  write(unitno2,*) 'VARLOCATION=([1-2]=NODAL,[3-7]=CELLCENTERED)'
	  write(unitno2,"(A,A)") 'NODES=', trim(formV)
	  write(unitno2,"(A,A)") 'ELEMENTS=', trim(formC)
	  write(unitno2,"(A,A)") 'FACES=',trim(formE)
	  write(unitno2,*) 'NUMCONNECTEDBOUNDARYFACES=0'
	  write(unitno2,*) 'TOTALNUMBOUNDARYCONNECTIONS=0'
	  write(unitno2,*) 'DATAPACKING=BLOCK'
	  write(unitno2,*) 'ZONETYPE=FEPOLYGON'
	  
	  ! Calculate Linelengh
	  lineLength=floor(30000.0/(14.0+1.0))
	  nLines=ceiling(dble(nvert)/dble(lineLength))
	  
	  do j=1,2
		lEnd=0
		do ii=1,nLines
			lStart=lEnd+1
			lEnd=min(lStart+lineLength-1,nvert)
			do np=(lStart),(lEnd-1)
				write(unitno2,"(ES14.6E2) ",advance="no") (grd(np,j))
			enddo
			write(unitno2,"(ES14.6E2) ") (grd(lEnd,j))
		enddo
	  enddo
	  lineLength=floor(30000.0/(30.0+1.0))
	  nLines=ceiling(dble(ncells)/dble(lineLength))
	call GENBLOCKVARIABLEREAL(unitno2,lineLength,nLines,ncells,rho)
	call GENBLOCKVARIABLEREAL(unitno2,lineLength,nLines,ncells,u)
	call GENBLOCKVARIABLEREAL(unitno2,lineLength,nLines,ncells,v)
	call GENBLOCKVARIABLEREAL(unitno2,lineLength,nLines,ncells,mach)
	call GENBLOCKVARIABLEREAL(unitno2,lineLength,nLines,ncells,cp)
	  
	  do ed=1,nedge
		  write(unitno2,*) (edge(ed,j),j=1,2)
	  enddo
	  lineLength=floor(30000.0/(30.0+1.0))
	  nLines=ceiling(dble(nedge)/dble(lineLength))
	  
	  do j=3,4
		lEnd=0
		do ii=1,nLines
			lStart=lEnd+1
			lEnd=min(lStart+lineLength-1,nedge)
			
			do np=(lStart),(lEnd-1)
				write(unitno2,"(I14) ",advance="no") (edge(np,j))
			enddo
			write(unitno2,"(I14) ") (edge(lEnd,j))
		enddo
	  enddo
	  
	  write(unitno2,*) 'ZONE t="Surface Data" N=',nvert
	  write(unitno2,*) 'STRANDID=2 '
	  write(unitno2,*) 'SOLUTIONTIME=',trim(adjustl(strT))
	  write(unitno2,*) 'DATAPACKING=POINT'
	  write(unitno2,*) 'ZONETYPE=FELINESEG'
	  do np=1,nvert
		write(unitno2,*) (grd(np,j),j=1,2),rhov(np),uv(np),vv(np),       &
	 &  machv(np),cpv(np)
	  enddo
	  
	  do ed=1,jj
		  write(unitno2,*) (edge(edgeList(ed),j),j=1,2)
	  enddo
	  
	  
	  fileno=5
	  blanking=0
	  call DUMPMESHVTK(nvert,nvert,nedge,nedge,edge,grd,                 &
	 &  blanking,fileno,rhov,uv,vv,machv,cpv)

      enddo
        

      


      end program



	! ********************************************
	subroutine GENBLOCKVARIABLEREAL(unitno2,lineLength,nLines,nvert,array)
	
	implicit none
	double precision, INTENT(IN),DIMENSION(nvert)::array
	INTEGER::unitno2,lineLength,nLines,nvert
	integer::np,ii,lStart,lEnd
	 
	lEnd=0
	do ii=1,nLines
		lStart=lEnd+1
		lEnd=min(lStart+lineLength-1,nvert)
			
		do np=(lStart),(lEnd-1)
			write(unitno2,"(ES14.6E2) ",advance="no") (array(np))
		enddo
		write(unitno2,"(ES14.6E2) ") (array(lEnd))
	enddo
	
	end subroutine
	
		! ********************************************
	
	 
!     **********************************************************      
      subroutine EDGESCELLS(nedgemax,nedge,edge,ncells,boundmax,        &
     &  nbound,bound,blanking)

      implicit none

      integer,intent(in)::nedgemax,ncells,boundmax,nedge
      integer,intent(in),dimension(nedgemax,4)::edge
      integer,intent(in),dimension(nedgemax)::blanking
      integer,intent(out),dimension(ncells,boundmax)::bound
      integer,intent(out),dimension(ncells)::nbound

      integer::i,left,right
      
      nbound=0
      bound=0
      !print *,nedge,ncells
      !pause
      do i=1,nedge

      if(blanking(i).eq.0)then

        left=edge(i,3)
        right=edge(i,4)
        if((left.gt.ncells).or.(right.gt.ncells))then
          print *, 'Too many cells!',i,ncells,left,right
          !stop
        endif
        !print *, 'A',left
        if(left.gt.0) then
          if(nbound(left).eq.boundmax)then
            print *, "Exceeded boundmax"
            stop
          endif
          nbound(left)=nbound(left)+1
          !print *, nbound(left)
          bound(left,nbound(left))=i
        endif
        !print *, 'B'        
        if(right.gt.0) then
          if(nbound(right).eq.boundmax)then
            print *, "Exceeded boundmax"
            stop
          endif 
          nbound(right)=nbound(right)+1
          bound(right,nbound(right))=i
        endif

      endif

      enddo

      end subroutine
!     **********************************************************      


!     **********************************************************
      subroutine DUMPMESHVTK(npt,nptmax,nedge,nedgemax,edge,grd,        &
     &  blanking,fileno,rhov,uv,vv,machv,cpv)

      implicit none

      integer,intent(in)::npt,nptmax,nedge,nedgemax,fileno
      double precision,intent(in),dimension(nptmax,2)::grd
      integer,intent(in),dimension(nedgemax,4)::edge
      integer,intent(in),dimension(nedgemax)::blanking
      double precision,intent(in),dimension(npt)::rhov,uv,vv,machv,cpv
      integer::unitno,ed,j,i,edlimit,np
      character(len=5)::filenochar1
      character(len=5)::filenochar

      unitno=888

      edlimit=0
      do ed=1,nedge
        if(blanking(ed).eq.0) edlimit=edlimit+1
      enddo

      !write(filenochar1,'(i2)') fileno
      !filenochar(1:2)=filenochar1
      !filenochar=adjustl(filenochar)
      !filenochar(2:5)='.vtk'

      open(unitno,file='flowvtk.vtk',status='unknown') 

4008  format('# vtk DataFile Version 2.0')
4007  format('mydata')
4006  format('ASCII')
4009  format('DATASET UNSTRUCTURED_GRID')
      
      write(unitno,4008) !'# vtk DataFile Version 2.0'
      write(unitno,4007) !'CutCell'
      write(unitno,4006) !'ASCII'
      write(unitno,4009) !'DATASET UNSTRUCTURED_GRID'
4005  format('POINTS ',2x,i6,2x,' float')      
      write(unitno,4005) npt
      do np=1,npt
        write(unitno,*) (grd(np,j),j=1,2),' 0.0'
      enddo

5000  format('CELLS',1x,i6,1x,i6) 
      write(unitno,5000) edlimit,3*edlimit           
            
      do ed=1,nedge
        if(blanking(ed).eq.0)then
          write(unitno,*) '2 ',(edge(ed,j)-1,j=1,2)
        endif
      enddo 
      
5001  format('CELL_TYPES',1x,i6)
      write(unitno,5001) edlimit
      do ed=1,nedge
        if(blanking(ed).eq.0)then
          write(unitno,*) '3'
        endif
      enddo
      
6003  format('POINT_DATA',2x,i6)      
6001  format('SCALARS rho float')
6002  format('LOOKUP_TABLE default')
      write(unitno,6003) npt
      write(unitno,6001)
      write(unitno,6002)
      do np=1,npt
        write(unitno,*) rhov(np)
      enddo
      
!6003  format('POINT_DATA',2x,i6)      
6004  format('SCALARS u float')
!6002  format('LOOKUP_TABLE default')
      !write(unitno,6003) npt
      write(unitno,6004)
      write(unitno,6002)
      do np=1,npt
        write(unitno,*) uv(np)
      enddo      
                  
!6003  format('POINT_DATA',2x,i6)      
6005  format('SCALARS v float')
!6002  format('LOOKUP_TABLE default')
      !write(unitno,6003) npt
      write(unitno,6005)
      write(unitno,6002)
      do np=1,npt
        write(unitno,*) vv(np)
      enddo      

!6003  format('POINT_DATA',2x,i6)      
6006  format('SCALARS Mach float')
!6002  format('LOOKUP_TABLE default')
      !write(unitno,6003) npt
      write(unitno,6006)
      write(unitno,6002)
      do np=1,npt
        write(unitno,*) machv(np)
      enddo      
                              
!6003  format('POINT_DATA',2x,i6)      
6007  format('SCALARS Cp float')
!6002  format('LOOKUP_TABLE default')
      !write(unitno,6003) npt
      write(unitno,6007)
      write(unitno,6002)
      do np=1,npt
        write(unitno,*) cpv(np)
      enddo      

      end subroutine
!     **********************************************************



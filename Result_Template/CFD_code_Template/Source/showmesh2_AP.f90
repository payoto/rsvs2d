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


open(100,file='griduns',status='unknown')
read(100,*) ncells,nedge,nvert

allocate(edge(nedge,4))
allocate(grd(nvert,3))
allocate(bound(ncells,boundmax))
allocate(nbound(ncells))
allocate(blanking(nedge))
allocate(edgept(nvert,nmeetmax))
allocate(nmeet(nvert))

allocate(edgeList(nedge))

do i=1,nedge
read(100,*) (edge(i,j),j=1,4)
enddo
do i=1,nvert
read(100,*) (grd(i,j),j=1,3)
enddo
unitno2=223
open(unitno2,file="gridplt_cell.plt",status='unknown')
unitno=222
open(unitno,file="gridplt.plt",status='unknown')
write(unitno2,*) 'variables = "x" "y"'



!close(100)


!close(101)




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

write(unitno,*) 'ZONE t="1" N=',nvert,'E=',nedge
write(unitno,*) 'STRANDID=1 '
write(unitno,*) 'SOLUTIONTIME=1'
write(unitno,*) 'DATAPACKING=POINT,'
write(unitno,*) 'ZONETYPE=FELINESEG'

do np=1,nvert
    write(unitno,*) (grd(np,j),j=2,3)
enddo
	  
do ed=1,nedge
	write(unitno,*) (edge(ed,j),j=1,2)
enddo

!write(unitno2,*) VARLOCATION=([1-2]=NODAL ,[3-5]=CELLCENTERED)
write(unitno2,*) 'ZONE '
write(unitno2,*) 'STRANDID=1 '
write(unitno2,*) 'SOLUTIONTIME=1'
write(unitno2,*) 'VARLOCATION=([1-2]=NODAL,[3-7]=CELLCENTERED)'
write(unitno2,"(A,A)") 'NODES=', trim(formV)
write(unitno2,"(A,A)") 'ELEMENTS=', trim(formC)
write(unitno2,"(A,A)") 'FACES=',trim(formE)
write(unitno2,*) 'NUMCONNECTEDBOUNDARYFACES=0'
write(unitno2,*) 'TOTALNUMBOUNDARYCONNECTIONS=0'
write(unitno2,*) 'DATAPACKING=BLOCK'
write(unitno2,*) 'ZONETYPE=FEPOLYGON'

! Calculate Linelengh
lineLength=floor(30000.0/(30.0+1.0))
nLines=ceiling(dble(nvert)/dble(lineLength))

do j=2,3
lEnd=0
do ii=1,nLines
lStart=lEnd+1
lEnd=min(lStart+lineLength-1,nvert)
do np=(lStart),(lEnd-1)
write(unitno2,"(ES30.20E2) ",advance="no") (grd(np,j))
enddo
write(unitno2,"(ES30.20E2) ") (grd(lEnd,j))
enddo
enddo
lineLength=floor(30000.0/(30.0+1.0))
nLines=ceiling(dble(ncells)/dble(lineLength))

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
write(unitno2,"(ES30.20E2) ",advance="no") (array(np))
enddo
write(unitno2,"(ES30.20E2) ") (array(lEnd))
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



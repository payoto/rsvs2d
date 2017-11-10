!      options(debug)

!     ***********************************************************
!                            FreeCart
!
!              2D Cut-cell cartesian mesh generator 
!                    
!                     Thomas C.S. Rendall
!                  University of Bristol, U.K.
!                         
!                   Version 1.0 of 20/09/09
!
!     ***********************************************************

!     *********************************************************** 

!     The strategy is:
!     1) Setup initial background mesh
!     2) Find intersections, and refine those cells repeatedly, with
!        buffering to neighbour cells to get a nice mesh
!     3) Find intersections
!     4) Sort all intersections along background edges and surface
!        edges
!     5) Introduce edges by stepping along background edges, and then
!        along surface edges
!     6) Final geometry - volumes, then compress lists and dump mesh    

      program cart
     
      implicit none
       
      integer::nedge,doesit,cuts,cnt,imax,jmax,npt,ncellsinit
      integer::boundmax,intsctmax,jj,agged,wallbot,walltop              
      integer::nedgemax,nptmax,nsurfed,nsurfpt
      integer,allocatable,dimension(:,:)::edge,surfed,bound,surfcut,    &
     &  edgenew,blankcon
      integer,allocatable,dimension(:)::dir,nbound,issurfcut,cutof
      real,allocatable,dimension(:,:)::grdout,grd,surfpt,grdnew
      real,allocatable,dimension(:,:)::intsct
      real,dimension(2)::botleft,topright,temp
      integer::cnt2,ct,accumulate
      !integer,allocatable,dimension(:,:)::intsects
      integer,allocatable,dimension(:)::cellpoint,blanking    
      
      integer::i,j,ed,ns,nc,bd,maxedcut,raydir,SF
      integer::IU,BM,DF,inside,np,otherdir,nlevs,offset

      real,dimension(2)::bigpos,pt,ed_a,ed_b 
      real::val,st,nd,test,pi,shiftx,shifty

      integer,allocatable,dimension(:)::inorout,contcell,addcuts,       &
     &  walllist,alias,wallflag
      integer::totalins,totalouts,dired,contained,fileno,tempint

      integer,allocatable,dimension(:,:)::edgecuts,surfcuts
      integer,allocatable,dimension(:)::nedgecuts,nsurfcuts,cellblanking
      integer::nedgecutsmax,nsurfcutsmax,switch,nptnew,nedgenew

      real,dimension(2)::start,next
      real,allocatable,dimension(:)::array,volume,vol
      integer,allocatable,dimension(:)::intrace,outtrace,temparray,     &
     &  ctpos,sfpos
      integer,allocatable,dimension(:,:)::connected

      integer::nextcell,lastcell,nmeetmax,nwallcells,nw,GU
      integer::sumlocal,cloops,splitcells,ncellsmax,ml

      real,allocatable,dimension(:,:)::grdloc
      integer,allocatable,dimension(:,:)::edgecomp,edgecompnew
      real,allocatable,dimension(:)::volnew
      integer,allocatable,dimension(:)::aliascell,buflev,doesita,agg
      integer::nedgenewfinal,ncellfinal,verb

      integer,allocatable,dimension(:)::refine,done,refinecell,nmeet
      real,allocatable,dimension(:,:)::sumgrd
      integer::refines,rf,possiblenexted,possiblenextpt,point,bufcells
      integer::ncellscurrent,ncellscurrent2,no1,no2,k
      integer::centpoint,itopleft,itopright,ibotleft,ibotright
      integer::left,right,outer,firstpt,lastpt,buf
      integer,allocatable,dimension(:,:)::edgept

      integer::other,nomore,nsurf,edsum,edstart,nsurfloc,MEM
      integer::verbfile
      real::volinit,contrib,dx,dy,volrat
      real,dimension(2)::midpt
      integer::nptout
      integer,allocatable,dimension(:)::aliaspoint
      
      print *,'*******************************************************'
      print *,'                        FreeCart                       '
      print *,' '
      print *,'          2D Cut-cell cartesian mesh generator         '
      print *,' '                     
      print *,'                    Thomas C.S. Rendall                '
      print *,'                University of Bristol, U.K.            '
      print *,'                  Version 1.0 of 20/09/09              '
      print *,'*******************************************************'
            
      IU=100
      BM=101
      DF=102 
      SF=103
      MEM=104 

      verb=0
      verbfile=0   

      pi=4.0*atan(1.0) 

      open(MEM,file='memalloc',status='unknown')

!     Max number of edges      
      !nedgemax=500000
      read(MEM,*) nedgemax
      print *,nedgemax
!     max number of points in mesh 
      !nptmax=100000
      read(MEM,*) nptmax
      print *,nptmax
!     Max edges per cell
      !boundmax=100
      read(MEM,*) boundmax
      print *,boundmax
!     Max total intersects      
      !intsctmax=50000
      read(MEM,*) intsctmax
      print *,intsctmax
!     Max number of intersects on a single edge      
      !maxedcut=100
      !read(MEM,*) maxedcut
      !print *,maxedcut
      !nedgecutsmax=200
      read(MEM,*) nedgecutsmax
      print *,nedgecutsmax
      !nsurfcutsmax=100
      read(MEM,*) nsurfcutsmax
      print *,nsurfcutsmax
      !maxsplits=500000
      !read(MEM,*) maxsplits
      !print *,maxsplits
      !ncellsmax=100000
      read(MEM,*) ncellsmax
      print *,ncellsmax
      !nmeetmax=20
      read(MEM,*) nmeetmax
      print *,nmeetmax
      
      close(MEM)

      open(SF,file='cutsettings',status='unknown')
      read(SF,*) imax,jmax,wallbot,walltop    

      nedge=2*(imax-1)*(jmax-1)+(imax-1)+(jmax-1)
      ncellsinit=(imax-1)*(jmax-1) 
      
      npt=imax*jmax
      
      !nsurfed=5
      !nsurfpt=5           
      
      !allocate(intsct(intsctmax,2))
      allocate(nbound(ncellsmax))
      allocate(bound(ncellsmax,boundmax))
      allocate(edge(nedgemax,4))
      allocate(edgenew(nedgemax,4))
      allocate(dir(nedgemax))
      allocate(surfcut(intsctmax,2))
      !allocate(intsects(nedgemax,maxedcut))
      !allocate(nintsects(nedgemax))     
      allocate(grd(nptmax,2))
      allocate(grdnew(nptmax,2))
      allocate(intsct(intsctmax,2))
      allocate(inorout(nptmax))
      allocate(addcuts(nptmax))
      allocate(blanking(nedgemax))
      allocate(edgecuts(nedgemax,nedgecutsmax))
      allocate(nedgecuts(nedgemax))
      allocate(cellblanking(ncellsmax))
      allocate(wallflag(ncellsmax))
      allocate(volume(ncellsmax))
      allocate(vol(ncellsmax))
      !allocate(tempvol(maxsplits))
      allocate(refine(nedgemax))
      allocate(done(ncellsmax))
      allocate(refinecell(nedgemax))
      allocate(sumgrd(ncellsmax,2))
      allocate(edgept(nptmax,nmeetmax))
      allocate(nmeet(nptmax))
      allocate(blankcon(nedgemax,3))
      allocate(doesita(nedgemax))
      allocate(agg(ncellsmax))

      blankcon=0
      
!     Base corners now set in file 'cutsettings'      
      
      read(SF,*) botleft(1),botleft(2)
      read(SF,*) topright(1),topright(2)
      read(SF,*) offset
      read(SF,*) nlevs
      allocate(buflev(nlevs))
      do i=1,nlevs
        read(SF,*) j,buflev(j)
      enddo
      close(SF)
      
      edge=0
      grd=0.0

      open(DF,file='boundary.dat',status='unknown')
      read(DF,*) nsurf
      print *, "Surfaces:",nsurf
      read(DF,*) nsurfed,nsurfpt
      print *, "Total surf edges:",nsurfed
      print *, "Total surf points:",nsurfpt

      allocate(surfpt(nsurfpt,2))
      allocate(surfed(nsurfed,2))
      allocate(issurfcut(nsurfed))
      allocate(cutof(nsurfed))
      allocate(cellpoint(nsurfpt))
      allocate(contcell(nsurfpt))
      allocate(surfcuts(nsurfed,nsurfcutsmax))
      allocate(nsurfcuts(nsurfed))

!     Read boundary.dat
      edsum=0
      do i=1,nsurf
        read(DF,*) nsurfloc
        print *, nsurfloc
        edstart=edsum+1
        do j=1,nsurfloc-1
          edsum=edsum+1
          surfed(edsum,1)=edsum
          surfed(edsum,2)=edsum+1
          read(DF,*) (surfpt(edsum,jj),jj=1,2)         
        enddo
        edsum=edsum+1
        read(DF,*) (surfpt(edsum,jj),jj=1,2)
        surfed(edsum,1)=edsum
        surfed(edsum,2)=edstart
        print *, i,edsum
      enddo
      
      !pause
      
!     Put points in general position
!     This is not a neat solution to degenerate intersections, but
!     works well for the time being.
      if(offset.eq.0)then
        shiftx=0.0
        shifty=0.0
      else
        shiftx=10.0**(-1.0*6.0)*pi/3.14
        shifty=10.0**(-1.0*6.0)*pi/3.14
      endif

      surfpt(:,1)=surfpt(:,1)+shiftx
      surfpt(:,2)=surfpt(:,2)+shifty
      
      open(BM,file='backgroundmesh.plt',status='unknown')
      write(BM,*) 'ZONE t="1" N=',npt,'E=',nedge
      write(BM,*) 'DATAPACKING=POINT,'
      write(BM,*) 'ZONETYPE=FELINESEG'

!     Begin setup of background mesh      
      do j=1,jmax
        do i=1,imax
       
          cnt=i+(j-1)*imax
        grd(cnt,1)=((topright(1)-botleft(1))/(imax-1))*(i-1)+botleft(1)
        grd(cnt,2)=((topright(2)-botleft(2))/(jmax-1))*(j-1)+botleft(2)
          
          write(BM,*) grd(cnt,1),grd(cnt,2)
      
        enddo
      enddo
      
      cnt=0
      ed=0

      volume=((topright(1)-botleft(1))/imax)*                           &
     &  ((topright(2)-botleft(2))/jmax)
      
      do j=1,jmax-1
        do i=1,imax-1
       
          cnt=i+(j-1)*imax
          cnt2=i+(j-1)*(imax-1)
          if(cnt2.gt.ncellsinit) print *, "WARNING 1"
          
          ed=ed+1
          edge(ed,1)=cnt
          edge(ed,2)=cnt+1
          edge(ed,3)=cnt2
          if(j.eq.1)then
            if(wallbot.eq.1)then
              edge(ed,4)=-1
            else
              edge(ed,4)=-2
            endif
          else
            edge(ed,4)=cnt2-(imax-1)
          endif
          
          dir(ed)=2
          write(BM,*) edge(ed,1),edge(ed,2)

          ed=ed+1
          edge(ed,1)=cnt
          edge(ed,2)=cnt+imax
          edge(ed,4)=cnt2
          if(i.eq.1)then
            edge(ed,3)=-2
          else
            edge(ed,3)=cnt2-1
          endif
          
          dir(ed)=1
          write(BM,*) edge(ed,1),edge(ed,2)
              
        enddo
      enddo

      do i=1,imax-1

        j=jmax
        cnt=i+(j-1)*imax
        cnt2=i+(j-2)*(imax-1)
        if(cnt2.gt.ncellsinit) print *, "WARNING 2"
        
        ed=ed+1
        edge(ed,1)=cnt
        edge(ed,2)=cnt+1
        if(walltop.eq.1)then
          edge(ed,3)=-1
        else
          edge(ed,3)=-2
        endif
        edge(ed,4)=cnt2
        dir(ed)=2
        write(BM,*) edge(ed,1),edge(ed,2)

      enddo
      
      do j=1,jmax-1

        i=imax
        cnt=i+(j-1)*imax
        cnt2=i-1+(j-1)*(imax-1)
        if(cnt2.gt.ncellsinit) print *, "WARNING 3"
        
        ed=ed+1
        edge(ed,1)=cnt
        edge(ed,2)=cnt+imax
        edge(ed,3)=cnt2
        edge(ed,4)=-2
        dir(ed)=1
        write(BM,*) edge(ed,1),edge(ed,2)

      enddo

      grdnew=grd
      edgenew=edge

      !fileno=9
      blanking=0
      fileno=500
      open(fileno,file='initialgrid.plt',status='unknown')   
      call DUMPMESHNM(npt,nptmax,nedge,nedgemax,edge,grd,blanking,      &
     &  fileno)
      close(fileno)
!     End background mesh      
      
      if(verbfile.eq.1)then
      open(IU,file='points.plt',status='unknown')
      endif

!     Begin refining split cells      
      ncellscurrent=ncellsinit
      !goto 222
! --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
!     Find all intersected edges XXX
      
! HERE --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------        

       !222 continue
       !stop

       print *, 'Building cell data...'
      bound=0
      nbound=0      
      call EDGESCELLS(nedgemax,nedge,edge,ncellsmax,boundmax,           &
     &  nbound,bound,blanking)
      print *, "...cell data built"
      fileno=7
      if(verbfile.eq.1)then
       call DUMPMESH(npt,nptmax,nedge,nedgemax,edge,grd,blanking,fileno)
      endif
 

!     Work out which cells contain which surface points. This test only
!     works for 4 boundaries, but this is ok since any cell conataining
!     a surf point will always have been refined. By definition, surf
!     points cannot lie on refinement boundaries.
      print *, "Finding containments..."       
      contcell=0
!NEW----
      addcuts=0
      contcell=0

      do np=1,nsurfpt
        
        do nc=1,ncellsinit
          addcuts(np)=0 
          do bd=1,nbound(nc)

          !bigpos(1)=1.0e6
          !bigpos(2)=grd(np,2)
          !raydir=2
          bigpos(1)=surfpt(np,1)!
          bigpos(2)=1.0e6
          raydir=1

          call INTERSECT(grd(edge(bound(nc,bd),1),:),                   &
     & grd(edge(bound(nc,bd),2),:),surfpt(np,:),bigpos,raydir,          &
     & temp,doesit)

          if(doesit.eq.1)then
            addcuts(np)=addcuts(np)+1
          endif

          enddo

          test=float(addcuts(np))/2.0-float(addcuts(np)/2)
          if(abs(test).gt.1.0e-6)then
            contcell(np)=nc
          endif
          
        enddo
      enddo


!     NEW----

      do np=1,nsurfpt
        !print *,contcell(np)
        !write(224,*) "ZONE"
        !write(224,*) (surfpt(np,j),j=1,2)
        nc=contcell(np)
        if(nc.gt.0)then
        do bd=1,nbound(nc)
          ed=bound(nc,bd)
         ! write(224,*) (grd(edge(ed,1),j),j=1,2)
         ! write(224,*) (grd(edge(ed,2),j),j=1,2)
        enddo
        endif
        if(contcell(np).eq.0)then
          print *, "Surface point not contained",np
          !stop
        endif
      enddo

      close(224)
!     Surf point containment complete


      cuts=0
      surfcut=0
      intsct=0.0
      !intsects=0
      !nintsects=0
    
      issurfcut=0   

       print *, 'Building cell data...'            
      call EDGESCELLS(nedgemax,nedge,edge,ncellsmax,boundmax,           &
     &  nbound,bound,blanking)
      print *, "...cell data built"

      cellpoint=0
      do nc=1,ncellsinit
         
        do np=1,nsurfpt
         inside=1  
          do bd=1,nbound(nc)

            ed=bound(nc,bd)
           
            if(dir(ed).eq.1) otherdir=2
            if(dir(ed).eq.2) otherdir=1

            val=surfpt(np,otherdir)
            st=grd(edge(ed,1),otherdir)
            nd=grd(edge(ed,2),otherdir)

            if(((val.le.st).and.(val.ge.nd)).or.((val.ge.st).and.       &
     &        (val.le.nd)))then
            else
              inside=0
            endif

          enddo 

          if(inside.eq.1)then
            cellpoint(np)=nc
          endif

        enddo
      enddo
      !print *, cellpoint


     




!     Find intersections, then work out which background points are
!     in/out. Odd number of intersections = inside, even = outside ('ray
!     tracing', or 'x-ray' algorithm)
      if(verbfile.eq.1)then      
      open(222,file='inside.plt',status='unknown')
      open(223,file='outside.plt',status='unknown')
      open(224,file='contains.plt',status='unknown')
      endif
      addcuts=0

      do np=1,npt
        do ns=1,nsurfed

          !bigpos(1)=1.0e6
          !bigpos(2)=grd(np,2)
          !raydir=2
          bigpos(1)=grd(np,1)!
          bigpos(2)=1.0e6
          raydir=1

          call INTERSECT(surfpt(surfed(ns,1),:),                        &
     &  surfpt(surfed(ns,2),:),grd(np,:),bigpos,raydir,temp,doesit)

          if(doesit.eq.1)then
            addcuts(np)=addcuts(np)+1
          endif

        enddo
      enddo

      inorout=0
      totalins=0
      totalouts=0
      do np=1,npt
        test=float(addcuts(np))/2.0-float(addcuts(np)/2)
        if(abs(test).gt.1.0e-6)then
          inorout(np)=1
          if(verbfile.eq.1)then
          write(222,*) (grd(np,j),j=1,2)
          endif
          totalins=totalins+1
        else
          inorout(np)=0
          if(verbfile.eq.1)then
          write(223,*) (grd(np,j),j=1,2)
          endif
          totalouts=totalouts+1
        endif
      enddo
      print *,"Background points inside:",totalins
      print *,"Background points outside:",totalouts

      !blanking=0
      !do ed=1,nedge
      !  if((inorout(edge(ed,1)).eq.1).and.(inorout(edge(ed,2)).eq.1))then
      !    blanking(ed)=1
      !  endif
      !enddo

!     Background point containment complete      

           

          
      fileno=1
      if(verbfile.eq.1)then    
      call DUMPMESH(npt,nptmax,nedge,nedgemax,edge,grd,blanking,fileno)
      endif


!     Now find all intersections and store them
      if(verbfile.eq.1)then      
      write(IU,*) "ZONE"
      endif
      
      print *, 'Finding intersections...' 

      nedgecuts=0
      nsurfcuts=0
      edgecuts=0
      surfcuts=0

      doesita=0

      do ns=1,nsurfed
        !do nc=1,ncellsinit
        !  do bd=1,nbound(nc)
        do ed=1,nedge
            !ed=bound(nc,bd)

            if(blanking(ed).eq.0)then

            call INTERSECT(surfpt(surfed(ns,1),:),                      &
     &  surfpt(surfed(ns,2),:),grd(edge(ed,1),:),grd(edge(ed,2),:),     &
     &  dir(ed),temp,doesit)
          
            if(doesit.eq.1)then

              doesita(ed)=doesit      
              cuts=cuts+1
              if(cuts.gt.intsctmax)then
                print *, "Intsctmax exceeded",cuts
                stop
              endif
              intsct(cuts,:)=temp(:)

              surfcut(cuts,1)=ns
              surfcut(cuts,2)=ed

!             List of intersections per ed              
              nedgecuts(ed)=nedgecuts(ed)+1
              if(nedgecuts(ed).gt.nedgecutsmax)then
                print *, "Nedgecutsmax exceeded",nedgecuts(ed)
                stop
              endif
              edgecuts(ed,nedgecuts(ed))=cuts

!             List of intersections per ns              
              nsurfcuts(ns)=nsurfcuts(ns)+1
              if(nsurfcuts(ns).gt.nsurfcutsmax)then
                print *, "Nsurfcutsmax exceeded",nsurfcuts(ns)
                stop
              endif
              surfcuts(ns,nsurfcuts(ns))=cuts

              !do i=1,nintsects(ed)
              !  if(ns.eq.intsects(ed,i)) canadd=0
              !enddo

              !if(canadd.eq.1)then
              !  nintsects(ed)=nintsects(ed)+1
              !  intsects(ed,nintsects(ed))=ns
              !endif

            endif

          !enddo

          

         endif

         

        enddo
      enddo

      do ed=1,nedge
        if(doesita(ed).eq.0)then
          if((inorout(edge(ed,1)).eq.1).and.(inorout(edge(ed,2)).eq.1))then
            blanking(ed)=1
          endif
        endif
      enddo


!     Now sort all intersections with distance, first for background
!     edges
      print *, 'Sorting background intersections...'

      if(verbfile.eq.1)then
      open(333,file='multicuts1.plt',status='unknown')
      endif

      do ed=1,nedge

        if(nedgecuts(ed).gt.0)then

          allocate(array(nedgecuts(ed)))
          allocate(intrace(nedgecuts(ed)))
          allocate(outtrace(nedgecuts(ed)))
          allocate(temparray(nedgecuts(ed)))

          start=grd(edge(ed,1),:)
          do j=1,nedgecuts(ed)
            ct=edgecuts(ed,j)
            next=intsct(ct,:)
            array(j)=((next(1)-start(1))**2+(next(2)-start(2))**2)**0.5
          enddo

          do j=1,nedgecuts(ed)
            intrace(j)=j!edgecuts(ed,j)
            temparray(j)=edgecuts(ed,j)
          enddo

          !print *,nedgecuts(ed)
          call INSERTIONSORT(nedgecuts(ed),array,intrace,outtrace)

          !write(777,*) "***"
          do j=1,nedgecuts(ed)
            edgecuts(ed,j)=temparray(outtrace(j))
           ! write(777,*) array(outtrace(j))
          enddo

          if(verbfile.eq.1)then
          write(333,*) "ZONE"
          write(333,*) grd(edge(ed,1),1),grd(edge(ed,1),2)
          do j=1,nedgecuts(ed)
            write(333,*) intsct(edgecuts(ed,j),1),                      &
     &      intsct(edgecuts(ed,j),2)
          enddo
          write(333,*) grd(edge(ed,2),1),grd(edge(ed,2),2)
          endif

          deallocate(array,intrace,outtrace,temparray)

        endif
      enddo


!     Now sort all intersections with distance, second for surface
!     edges
      print *, 'Sorting surface face intersections...'

      if(verbfile.eq.1)then
      open(334,file='multicuts2.plt',status='unknown')
      endif
      !write(777,*) "aaa"

      do ns=1,nsurfed

        if(nsurfcuts(ns).gt.0)then

          allocate(array(nsurfcuts(ns)))
          allocate(intrace(nsurfcuts(ns)))
          allocate(outtrace(nsurfcuts(ns)))
          allocate(temparray(nsurfcuts(ns)))

          start=surfpt(surfed(ns,1),:)
          do j=1,nsurfcuts(ns)
            ct=surfcuts(ns,j)
            next=intsct(ct,:)
            array(j)=((next(1)-start(1))**2+(next(2)-start(2))**2)**0.5
          enddo

          do j=1,nsurfcuts(ns)
            intrace(j)=j!surfcuts(ns,j)
            temparray(j)=surfcuts(ns,j)
          enddo

          !print *,ns,nsurfcuts(ns)
          call INSERTIONSORT(nsurfcuts(ns),array,intrace,outtrace)
          !write(777,*) "***"
          do j=1,nsurfcuts(ns)
            surfcuts(ns,j)=temparray(outtrace(j))
          !  write(777,*) array(outtrace(j))
          enddo

          if(verbfile.eq.1)then
          write(334,*) "ZONE"
          write(334,*) surfpt(surfed(ns,1),1),surfpt(surfed(ns,1),2)
          do j=1,nsurfcuts(ns)
            write(334,*) intsct(surfcuts(ns,j),1),                      &
     &      intsct(surfcuts(ns,j),2)
          enddo
          write(334,*) surfpt(surfed(ns,2),1),surfpt(surfed(ns,2),2)
          endif

          deallocate(array,intrace,outtrace,temparray)

        endif
      enddo





!     Nebulous process of adding edges. First copy unblanked edges and
!     then start adding new edges. Be aware no repeated points allowed
!     so must keep a note of where intersections are added in grd by
!     using ctpos, and where surface points are added by using sfpos.
!     Use inorout to find whether we are inside or outside. Only add
!     edges when inside, which is determined with 'switch'.

      allocate(ctpos(cuts))
      allocate(sfpos(nsurfpt))

      nptnew=npt
      nedgenew=nedge
      do ed=1,nedge

        if(nedgecuts(ed).gt.0)then

          blanking(ed)=1

          if(inorout(edge(ed,1)).eq.0) switch=1
          if(inorout(edge(ed,1)).eq.1) switch=0
                
          do j=1,nedgecuts(ed)

            if(j.eq.1)then
              
                nptnew=nptnew+1
                if(nptnew.gt.nptmax)then
                  print *,"Npt exceeded",nptnew
                  stop
                endif
                ct=edgecuts(ed,j)
                grd(nptnew,:)=intsct(ct,:)
                ctpos(ct)=nptnew

                nedgenew=nedgenew+1
                if(nedgenew.gt.nedgemax)then
                  print *,"Edgemax exceeded",nedgenew
                  stop
                endif
                edge(nedgenew,1)=edge(ed,1)
                edge(nedgenew,2)=nptnew
                edge(nedgenew,3)=edge(ed,3)
                edge(nedgenew,4)=edge(ed,4)

                tempint=switch
                if(switch.eq.0) blanking(nedgenew)=1
                if(tempint.eq.0) switch=1
                if(tempint.eq.1) switch=0

            endif

            if(j.gt.1)then

                nedgenew=nedgenew+1
                if(nedgenew.gt.nedgemax)then
                  print *,"Edgemax exceeded",nedgenew
                  stop
                endif    
                edge(nedgenew,1)=nptnew    

                nptnew=nptnew+1
                if(nptnew.gt.nptmax)then
                  print *,"Npt exceeded",nptnew
                  stop
                endif
                ct=edgecuts(ed,j)
                grd(nptnew,:)=intsct(ct,:)
                ctpos(ct)=nptnew
                
                edge(nedgenew,2)=nptnew
                edge(nedgenew,3)=edge(ed,3)
                edge(nedgenew,4)=edge(ed,4)

                tempint=switch
                if(switch.eq.0) blanking(nedgenew)=1
                if(tempint.eq.0) switch=1
                if(tempint.eq.1) switch=0

            endif

          enddo

          nedgenew=nedgenew+1
          if(nedgenew.gt.nedgemax)then
            print *,"Edgemax exceeded",nedgenew
            stop
          endif
          edge(nedgenew,1)=nptnew
          edge(nedgenew,2)=edge(ed,2)
          edge(nedgenew,3)=edge(ed,3)
          edge(nedgenew,4)=edge(ed,4)

          if(inorout(edge(ed,2)).eq.0) switch=1
          if(inorout(edge(ed,2)).eq.1) switch=0

          if(switch.eq.0) blanking(nedgenew)=1
          !if(tempint.eq.0) switch=1
          !if(tempint.eq.1) switch=0

        endif
      enddo

      fileno=2    
      if(verbfile.eq.1)then
      call DUMPMESH(nptnew,nptmax,nedgenew,nedgemax,edge,grd,blanking,  &
     &  fileno)
      endif





!     Do the same process along surface edges. Instead of finding if we
!     are in or out, the challenge here is to know which cell we are in.
!     Start of with the cell that contains the first edge point, then
!     change it over as we cross the edges. This loop assumes the solid
!     body is ON THE LEFT.
      sfpos=0

      !nptnew=npt
      !nedgenew=nedge
      do ns=1,nsurfed

        !lastcell=contcell(surfed(ns,1))

        
        if(nsurfcuts(ns).gt.0)then

          lastcell=contcell(surfed(ns,1))                
          !print *,ns,nsurfcuts(ns),lastcell

          do j=1,nsurfcuts(ns)

            if(j.eq.1)then

                nedgenew=nedgenew+1
                if(nedgenew.gt.nedgemax)then
                  print *,"Edgemax exceeded",nedgenew
                  stop
                endif

                !nptnew=nptnew+1
                !grd(nptnew,:)=surfpt(surfed(ns,1),:)
                !edge(nedgenew,1)=nptnew

                if(sfpos(surfed(ns,1)).eq.0)then
                  nptnew=nptnew+1
                  if(nptnew.gt.nptmax)then
                    print *,"Npt exceeded",nptnew
                    stop
                  endif
                  grd(nptnew,:)=surfpt(surfed(ns,1),:)
                  edge(nedgenew,1)=nptnew
                  sfpos(surfed(ns,1))=nptnew
                else
                  edge(nedgenew,1)=sfpos(surfed(ns,1))
                endif
                    
                ct=surfcuts(ns,j)

                edge(nedgenew,2)=ctpos(ct)

                edge(nedgenew,3)=-1!edge(ed,3)
                edge(nedgenew,4)=contcell(surfed(ns,1))
                lastcell=contcell(surfed(ns,1))

            endif

            if(j.gt.1)then

                nedgenew=nedgenew+1
                if(nedgenew.gt.nedgemax)then
                  print *,"Edgemax exceeded",nedgenew
                  stop
                endif    
                !edge(nedgenew,1)=nptnew
                edge(nedgenew,1)=ctpos(surfcuts(ns,j-1))

                !nptnew=nptnew+1
                ct=surfcuts(ns,j)
                !grd(nptnew,:)=intsct(ct,:)
                
                edge(nedgenew,2)=ctpos(surfcuts(ns,j))

                ed=surfcut(surfcuts(ns,j-1),2)

                left=edge(ed,3)
                right=edge(ed,4)
                if((lastcell.ne.left).and.(lastcell.ne.right))then
                  print *, "whoa!!!",lastcell,left,right
                  stop
                endif
                if(lastcell.eq.left)then
                  nextcell=right
                endif
                if(lastcell.eq.right)then
                  nextcell=left
                endif
                
                edge(nedgenew,3)=-1
                edge(nedgenew,4)=nextcell
                lastcell=nextcell
                
            endif

          enddo

          nedgenew=nedgenew+1
          if(nedgenew.gt.nedgemax)then
            print *,"Edgemax exceeded",nedgenew
            stop
          endif
          edge(nedgenew,1)=ctpos(surfcuts(ns,nsurfcuts(ns)))

          !nptnew=nptnew+1
          !grd(nptnew,:)=surfpt(surfed(ns,2),:)
          !edge(nedgenew,2)=nptnew
          if(sfpos(surfed(ns,2)).eq.0)then
            nptnew=nptnew+1
            if(nptnew.gt.nptmax)then
              print *,"Npt exceeded",nptnew
              stop
            endif
            grd(nptnew,:)=surfpt(surfed(ns,2),:)
            edge(nedgenew,2)=nptnew
            sfpos(surfed(ns,2))=nptnew
          else
            edge(nedgenew,2)=sfpos(surfed(ns,2))
          endif
          
          edge(nedgenew,3)=-1
          edge(nedgenew,4)=contcell(surfed(ns,2))

          lastcell=contcell(surfed(ns,2))
          

        else  !No cuts exist, so just add the two points and an edge

          nedgenew=nedgenew+1
          if(nedgenew.gt.nedgemax)then
            print *,"Edgemax exceeded",nedgenew
            stop
          endif

          if(sfpos(surfed(ns,1)).eq.0)then
            nptnew=nptnew+1
            if(nptnew.gt.nptmax)then
              print *,"Npt exceeded",nptnew
              stop
            endif
            grd(nptnew,:)=surfpt(surfed(ns,1),:)
            edge(nedgenew,1)=nptnew
            sfpos(surfed(ns,1))=nptnew
          else
            edge(nedgenew,1)=sfpos(surfed(ns,1))
          endif

          if(sfpos(surfed(ns,2)).eq.0)then
            nptnew=nptnew+1
            if(nptnew.gt.nptmax)then
              print *,"Npt exceeded",nptnew
              stop
            endif
            grd(nptnew,:)=surfpt(surfed(ns,2),:)
            edge(nedgenew,2)=nptnew
            sfpos(surfed(ns,2))=nptnew
          else
            edge(nedgenew,2)=sfpos(surfed(ns,2))
          endif

          edge(nedgenew,3)=-1
          edge(nedgenew,4)=contcell(surfed(ns,1))

          lastcell=contcell(surfed(ns,1))

        endif
        !print *,lastcell

      enddo

      fileno=501
      open(501,file='grid.plt',status='unknown')    
      call DUMPMESHNM(nptnew,nptmax,nedgenew,nedgemax,edge,grd,blanking,&
     &  fileno)
      close(501)
      fileno=502
      open(502,file='grid.vtk',status='unknown')
      call DUMPMESHVTKNM(nptnew,nptmax,nedgenew,nedgemax,edge,grd,      &
     &  blanking,fileno)
      close(502)
      !stop

      if(verb.eq.1) print *,nedgenew

      cellblanking=0
!     Get cell blanking. It might be useful, but not used at present      
      do nc=1,ncellsinit
        cellblanking(nc)=0
        do bd=1,nbound(nc)
          if(blanking(bound(nc,bd)).eq.1) cellblanking(nc)=1
        enddo
      enddo      


       print *, 'Building cell data...'            
      call EDGESCELLS(nedgemax,nedgenew,edge,ncellsmax,boundmax,        &
     &  nbound,bound,blanking)
      print *, "...cell data built"

      

      if(verb.eq.1) print *, npt+cuts+nsurfpt,nptnew

!     Find number of wall cells      
      nwallcells=0
      wallflag=0
      do nc=1,ncellsinit
        do i=1,nbound(nc)
       if((edge(bound(nc,i),3).eq.-1).or.(edge(bound(nc,i),4).eq.-1))then
       if(wallflag(nc).eq.0)then
         nwallcells=nwallcells+1
         wallflag(nc)=1
       endif
       endif
        enddo
      enddo

      allocate(walllist(nwallcells))

!     Write list of wall cells      
      nwallcells=0
      wallflag=0
      do nc=1,ncellsinit
        do i=1,nbound(nc)
       if((edge(bound(nc,i),3).eq.-1).or.(edge(bound(nc,i),4).eq.-1))then
       if(wallflag(nc).eq.0)then        
         nwallcells=nwallcells+1
         walllist(nwallcells)=nc
         wallflag(nc)=1
       endif
       endif
        enddo
      enddo


      allocate(alias(nptnew))
      



!     This section checks all wall cells to see if they are split      
      if(verb.eq.1) print *,"wall",nwallcells
      !pause
      open(999,file='splitcells.plt',status='unknown')
      splitcells=0
      accumulate=ncellsinit

      do nw=1,nwallcells

        nc=walllist(nw)
        allocate(connected(nbound(nc),4))

        sumlocal=0
        alias=0
        if(verbfile.eq.1)then
         open(555,file='check.plt',status='unknown')
         open(556,file='checkint.dat',status='unknown')
         write(556,*) nw
        endif
        

        do i=1,nbound(nc)
          if(verbfile.eq.1)then
          write(555,*) "ZONE"
          write(556,*) "ZONE"
          endif
          !if(blanking(bound(nc,i)).eq.0)then

!           INDEGRAPH works in local indices, so make a mini list and
!           keep if for later in 'alias'          
            if(alias(edge(bound(nc,i),1)).eq.0)then
              sumlocal=sumlocal+1
              alias(edge(bound(nc,i),1))=sumlocal
            endif
            if(alias(edge(bound(nc,i),2)).eq.0)then
              sumlocal=sumlocal+1
              alias(edge(bound(nc,i),2))=sumlocal
            endif
            if(verbfile.eq.1)then
      write(555,*) grd(edge(bound(nc,i),1),1),grd(edge(bound(nc,i),1),2)
      write(555,*) grd(edge(bound(nc,i),2),1),grd(edge(bound(nc,i),2),2)
            write(556,*) edge(bound(nc,i),1),edge(bound(nc,i),2)
            endif

          !endif
        enddo

!       For area calc INDEGRAPH needs coordinates, so make a local
!       coordinate array        
        allocate(grdloc(sumlocal,2))
        alias=0
        sumlocal=0
        do i=1,nbound(nc)
            if(alias(edge(bound(nc,i),1)).eq.0)then
              sumlocal=sumlocal+1
              alias(edge(bound(nc,i),1))=sumlocal
              grdloc(sumlocal,:)=grd(edge(bound(nc,i),1),:)
            endif
            if(alias(edge(bound(nc,i),2)).eq.0)then
              sumlocal=sumlocal+1
              alias(edge(bound(nc,i),2))=sumlocal
              grdloc(sumlocal,:)=grd(edge(bound(nc,i),2),:)
            endif
        enddo

        if(verbfile.eq.1)then
        close(555)
        close(556)
        endif

!       Copy the aliased point indices to the local edge connectivity
!       array 'connected'        
        !eflag=0
        do i=1,nbound(nc)
          if(blanking(bound(nc,i)).eq.0)then
            connected(i,1)=alias(edge(bound(nc,i),1))
            connected(i,2)=alias(edge(bound(nc,i),2))
            connected(i,3)=edge(bound(nc,i),3)
            connected(i,4)=edge(bound(nc,i),4)
         !   if(edge(bound(nc,i),1).eq.62639) eflag=1
         !   if(edge(bound(nc,i),2).eq.62639) eflag=1
          endif
        enddo

!       Scan to find all distinct closed loops and their areas
        !print *,verb
        !pause
        call INDEGRAPH(sumlocal,nmeetmax,nbound(nc),connected,          &
     &    cloops,nc,accumulate,verb)

        if(verb.eq.1) print *, nw,ncellsinit,nc,accumulate
        if(accumulate.gt.ncellsmax)then
          print *, "Ncellsmax exceeded with split cells",accumulate
          stop
        endif

!       Copy the areas of the closed loops out again        
!        do clps=1,cloops
!          if(clps.eq.1)then
!            volume(nc)=tempvol(clps)
!          endif
!          if(clps.gt.1)then
!            volume(accumulate)=tempvol(clps)
!          endif
!        enddo

!       Copy the connectivity out again - this will have changed if the
!       cell was split, but will be unchanged if it was not        
        do i=1,nbound(nc)
          if(blanking(bound(nc,i)).eq.0)then
            edge(bound(nc,i),3)=connected(i,3)
            edge(bound(nc,i),4)=connected(i,4)
          endif
        enddo

        deallocate(connected,grdloc)


!       Dump the split cells for viewing
        if(cloops.gt.1)then
        do i=1,nbound(nc)
          write(999,*) "ZONE"
      write(999,*) grd(edge(bound(nc,i),1),1),grd(edge(bound(nc,i),1),2)
      write(999,*) grd(edge(bound(nc,i),2),1),grd(edge(bound(nc,i),2),2)
        enddo
        splitcells=splitcells+1
        endif



      enddo
      
      
      
      
      

      agg=0
      agged=0
      !volrat=100000000
      volrat=10





      vol=0.0
      do i=1,nedgenew
        if(blanking(i).eq.0)then
          midpt=0.5*(grd(edge(i,1),:)+grd(edge(i,2),:))
          dx=grd(edge(i,2),1)-grd(edge(i,1),1)
          dy=grd(edge(i,2),2)-grd(edge(i,1),2)
          contrib=0.5*(midpt(1)*(dy)+midpt(2)*(-dx))                
        if(edge(i,3).gt.0)then
          vol(edge(i,3))=vol(edge(i,3))+contrib
        endif
        if(edge(i,4).gt.0)then
          vol(edge(i,4))=vol(edge(i,4))-contrib
        endif
        endif
      enddo
      vol=abs(vol)

      agg=0

      do i=1,nedgenew

      if(blanking(i).eq.0)then
      
          if((edge(i,3).gt.0).and.(edge(i,4).gt.0))then
            if((agg(edge(i,3)).eq.0).and.(agg(edge(i,4)).eq.0))then
            !print *,vol(edge(i,4)),vol(edge(i,3))

            if((vol(edge(i,4))/vol(edge(i,3))).gt.volrat)then
              agg(edge(i,3))=edge(i,4)
              agg(edge(i,4))=-1
              !vol(edge(i,3))=vol(edge(i,3))+vol(edge(i,4))
              agged=agged+1
              blanking(i)=1
            endif
            if(blanking(i).eq.0)then

            if((vol(edge(i,3))/vol(edge(i,4))).gt.volrat)then
              agg(edge(i,4))=edge(i,3)
              agg(edge(i,3))=-1
              !vol(edge(i,4))=vol(edge(i,4))+vol(edge(i,3))
              agged=agged+1
              blanking(i)=1
            endif

            endif

            endif
          endif
          
      endif
             
      enddo

      do i=1,nedgenew
      
      if(blanking(i).eq.0)then
              
                              
        if((edge(i,3).gt.0))then
          if(agg(edge(i,3)).gt.0)then
            edge(i,3)=agg(edge(i,3))
            !ncellfinal=ncellfinal-1
          endif
        endif
        if((edge(i,4).gt.0))then
          if(agg(edge(i,4)).gt.0)then
            edge(i,4)=agg(edge(i,4))
            !ncellfinal=ncellfinal-1
          endif
        endif

        if(edge(i,3).eq.edge(i,4)) blanking(i)=1
        !endif
       
      endif

      enddo

      fileno=8 
      fileno=503
      open(503,file='grid_merged.plt',status='unknown')   
      call DUMPMESHNM(nptnew,nptmax,nedgenew,nedgemax,edge,grd,blanking,&
     &  fileno)
      close(503)
      fileno=504
      open(504,file='grid_merged.vtk',status='unknown')   
      call DUMPMESHVTKNM(nptnew,nptmax,nedgenew,nedgemax,edge,grd,      &
     &  blanking,fileno)
      close(504)
               

!     Prepare to write the mesh file. This is complicated by the need to
!     compress all the lists first (solver needs the totals)      
      GU=898
      open(GU,file='griduns',status='unknown')

      allocate(edgecomp(nedgenew,4))
      allocate(edgecompnew(nedgenew,4))
      allocate(aliascell(accumulate))
      allocate(aliaspoint(nptnew))
      allocate(grdout(nptnew,2))     

!     Compress edges and points
      nedgenewfinal=0
      nptout=0
      aliaspoint=0
      do i=1,nedgenew
        if(blanking(i).eq.0)then
          nedgenewfinal=nedgenewfinal+1
          edgecomp(nedgenewfinal,:)=edge(i,:)
          edgecompnew(nedgenewfinal,:)=edge(i,:)
          if(aliaspoint(edge(i,1)).eq.0)then
            nptout=nptout+1
            grdout(nptout,:)=grd(edge(i,1),:)
            edgecomp(nedgenewfinal,1)=nptout
            edgecompnew(nedgenewfinal,1)=nptout         
            aliaspoint(edge(i,1))=nptout
          else
!           Cell has an alias, so use that  
            edgecomp(nedgenewfinal,1)=aliaspoint(edge(i,1))
            edgecompnew(nedgenewfinal,1)=aliaspoint(edge(i,1))          
          endif
          if(aliaspoint(edge(i,2)).eq.0)then
            nptout=nptout+1
            grdout(nptout,:)=grd(edge(i,2),:)
            edgecomp(nedgenewfinal,2)=nptout
            edgecompnew(nedgenewfinal,2)=nptout             
            aliaspoint(edge(i,2))=nptout
          else
!           Cell has an alias, so use that  
            edgecomp(nedgenewfinal,2)=aliaspoint(edge(i,2))
            edgecompnew(nedgenewfinal,2)=aliaspoint(edge(i,2))              
          endif
        endif       
      enddo

!     Compress cells. alias cell is the compressed cell index.      
      aliascell=0
      ncellfinal=0

      do i=1,nedgenewfinal

        if(edgecomp(i,3).gt.0)then
!         If cell has no alias                
          if(aliascell(edgecomp(i,3)).eq.0)then
            ncellfinal=ncellfinal+1
            edgecompnew(i,3)=ncellfinal
            !volnew(ncellfinal)=volume(edgecomp(i,3))
            aliascell(edgecomp(i,3))=ncellfinal
          else
!           Cell has an alias, so use that                  
            edgecompnew(i,3)=aliascell(edgecomp(i,3))
            !volnew(aliascell(edgecomp(i,3)))=volume(edgecomp(i,3))
          endif
        endif

        if(edgecomp(i,4).gt.0)then
          if(aliascell(edgecomp(i,4)).eq.0)then
            ncellfinal=ncellfinal+1
            edgecompnew(i,4)=ncellfinal
            !volnew(ncellfinal)=volume(edgecomp(i,4))
            aliascell(edgecomp(i,4))=ncellfinal
          else
            edgecompnew(i,4)=aliascell(edgecomp(i,4))
            !volnew(aliascell(edgecomp(i,4)))=volume(edgecomp(i,4))
          endif
        endif

      enddo

      allocate(volnew(ncellfinal))

      volnew=0.0
      do i=1,nedgenewfinal
        midpt=0.5*(grdout(edgecompnew(i,1),:)+grdout(edgecompnew(i,2),:))
        dx=grdout(edgecompnew(i,2),1)-grdout(edgecompnew(i,1),1)
        dy=grdout(edgecompnew(i,2),2)-grdout(edgecompnew(i,1),2)
        contrib=0.5*(midpt(1)*(dy)+midpt(2)*(-dx))
        if(edgecompnew(i,3).gt.0)then
          volnew(edgecompnew(i,3))=volnew(edgecompnew(i,3))+contrib
        endif
        if(edgecompnew(i,4).gt.0)then
          volnew(edgecompnew(i,4))=volnew(edgecompnew(i,4))-contrib
        endif
      enddo
      volnew=abs(volnew)
      
 
!     Dump the mesh out
      write(GU,*) ncellfinal,nedgenewfinal,nptout
      do i=1,nedgenewfinal
        write(GU,*) (edgecompnew(i,j),j=1,4)
      enddo
      do i=1,nptout
        write(GU,*) i,(grdout(i,j),j=1,2)
      enddo
      do i=1,ncellfinal
        write(GU,*) volnew(i)
      enddo
      close(GU)
    
!     Lets have the news
123   format(1x,'A total of ',i4,' intersections were found.')        
      write(6,123) cuts
124   format(1x,'A total of ',i4,' split cells were found.')          
      write(6,124) splitcells
125   format(1x,'A total of ',i4,' cells were merged.')          
      write(6,125) agged
      !pause      

      end program
!     **********************************************************

                 







subroutine BACGROUNDMESH()

    implicit none
    do j=1,nlevs !ref level loop
        print *, "Refinement level",j
        refine=0
      
        do ns=1,nsurfed
          do ed=1,nedge

              call INTERSECT(surfpt(surfed(ns,1),:),                    &
     &    surfpt(surfed(ns,2),:),grd(edge(ed,1),:),grd(edge(ed,2),:),   &
     &    dir(ed),temp,doesit)
          
              if(doesit.eq.1)then
                refine(ed)=1
              endif

          enddo
        enddo

        print *, 'Building cell data...'            
        call EDGESCELLS(nedgemax,nedge,edge,ncellsmax,boundmax,         &
     &   nbound,bound,blanking)
        print *, "...cell data built"

!       sumgrd is accumulator for cell centres, refinecell is switch to
!       refine any given cell
        sumgrd=0.0
        refinecell=0

!       refines is number of refined cells, done is an array to stop any
!       cell being refined more than once (but child cells may be
!       refined once more, of course!)        
        refines=0   
        done=0
        do nc=1,ncellscurrent
          do i=1,nbound(nc)

!         Don't refine if there are not 4 boundaries - the loops are not
!         set up to do this.
          if(nbound(nc).eq.4)then

            if((refine(bound(nc,i)).eq.1).and.(done(nc).eq.0))then
              refines=refines+1
              refinecell(refines)=nc
              done(nc)=1
            endif

            sumgrd(nc,:)=sumgrd(nc,:)+grd(edge(bound(nc,i),1),:)+       &
     &         grd(edge(bound(nc,i),2),:)

          endif

         enddo
       enddo

       !goto 333
       bufcells=0

!      Need to buffer the refinement a few times to broaden the refined
!      region XXX
       print *, "Buffering..."
       do buf=1,buflev(j)

        do nc=1,ncellscurrent
          nomore=0
          do i=1,nbound(nc)

!           Make sure there are 4 boundaries to the cell. nomore stops
!           any cell being refined as a result of more than one of its
!           edges boardering a refined cell. Obviously, we can only
!           refine once per cell per cycle.           
            if(nbound(nc).eq.4)then
              if(nomore.eq.0)  then      

                if(nc.eq.edge(bound(nc,i),3))then
                  other=edge(bound(nc,i),4)
                 endif
                if(nc.eq.edge(bound(nc,i),4))then
                  other=edge(bound(nc,i),3)
                endif
                if((nc.ne.edge(bound(nc,i),4)).and.                     &
     &             (nc.ne.edge(bound(nc,i),3)))then
                if(verb.eq.1) print *, "hhhhmn"
                  stop
                 endif   

                if(other.gt.0)then
                  if(done(other).eq.buf)then
                    if(done(nc).eq.0)then        
                      refines=refines+1
                      refinecell(refines)=nc
                      done(nc)=buf+1
                      nomore=1
                      bufcells=bufcells+1
                    endif
                  endif
                endif
           
              endif

            endif

         enddo
       enddo

       enddo

       if(verb.eq.1) print *, "buf",bufcells,refines,nomore
       !pause

!333    continue       

       ncellscurrent2=ncellscurrent
!       ncellscurrent1=ncellscurrent

       print *, "Refining..."

!      Start loop to refine all cells tagged for refinement
       do rf=1,refines

         nc=refinecell(rf)
        
         nmeet=0
         edgept=0

         npt=npt+1
         if(npt.gt.nptmax)then
           print *,"Npt exceeded",npt
           stop
         endif
         grd(npt,:)=sumgrd(nc,:)/8.0
         centpoint=npt
     
!        Build the point-edge list to allow traversal of circumference
!        later         
         do i=1,nbound(nc)
      
           no1=edge(bound(nc,i),1)
           no2=edge(bound(nc,i),2)

           nmeet(no1)=nmeet(no1)+1
           if(nmeet(no1).gt.nmeetmax)then
             print *, "Nmeetmax exceeded",nmeet(no1)
             stop
           endif
           edgept(no1,nmeet(no1))=bound(nc,i)

           nmeet(no2)=nmeet(no2)+1
           if(nmeet(no2).gt.nmeetmax)then
             print *, "Nmeetmax exceeded",nmeet(no2)
             stop
           endif
           edgept(no2,nmeet(no2))=bound(nc,i)

         enddo

         do i=1,nbound(nc)
           if(verb.eq.1) print *, nmeet(edge(bound(nc,i),1))
           if(verb.eq.1) print *, nmeet(edge(bound(nc,i),2))
         enddo

         if(verb.eq.1) print *,nc

         if(verb.eq.1) print *, "refining"

!        Important step - we are going to traverse the boundaries
!        anti-clockwise. Remember that edge(ed,3) is the cell on the
!        left. These 2 ifs start loop anti-clockwise         
         ed=bound(nc,1)
         if(nc.eq.edge(ed,3))then
           point=edge(bound(nc,1),1)
         endif
         if(nc.eq.edge(ed,4))then
           point=edge(bound(nc,1),2)
         endif

!        ibotleft,ibotright,itopright,itopleft all point to the four
!        corners of the split cell. Note that the current cell is
!        rearranged to be nc, so ibotleft=nc (i.e. only 3 cells are
!        actually added to the list, the first is just a modification)         
         ibotleft=nc

         ncellscurrent2=ncellscurrent2+1

         ibotright=ncellscurrent2
         ncellscurrent2=ncellscurrent2+1

         itopright=ncellscurrent2
         ncellscurrent2=ncellscurrent2+1

         itopleft=ncellscurrent2
         if(ncellscurrent2.gt.ncellsmax)then
           print *, "Ncells max exceeded at refinement",ncellscurrent2
           stop
         endif

!        Don't forget to split the volumes         
         volinit=volume(nc)
         volume(nc)=volinit/4.0
         volume(ibotright)=volinit/4.0
         volume(itopright)=volinit/4.0
         volume(itopleft)=volinit/4.0

!        Not possible to refine unless 4 boundaries, so check this         
         if(nbound(nc).ne.4) then
           print *, "not 4 boundaries"
           stop
         endif

         do k=1,nbound(nc)

           !write(808,*) "ZONE"

!          Traverse anti-clock to the next point (known to be anti-clock
!          due to pre-loop setup)           
           possiblenextpt=edge(ed,1)
           if(possiblenextpt.ne.point)then
             point=possiblenextpt
           else
             point=edge(ed,2)
           endif

           if(k.gt.4) then
             print *, "EEERRROOO"
             !stop
           endif

!          Need to watch out for the edge being the 'other way around'.
!          We are going aniclockwise and need to preserve this,           
           if(nc.eq.edge(ed,3))then
             outer=edge(ed,4)
             if(verb.eq.1) print *, "onleft"
             firstpt=edge(ed,1)
             lastpt=edge(ed,2)
!             ist=1
!             isp=2
           endif
           if(nc.eq.edge(ed,4))then
             outer=edge(ed,3)
             if(verb.eq.1) print *, "onright"
             firstpt=edge(ed,2)
             lastpt=edge(ed,1)
!             ist=2
!             isp=1
           endif

!          Some last-ditch error catching :-)           
           if((nc.ne.edge(ed,3)).and.(nc.ne.edge(ed,4)))then
             print *, "uh-oh"
             stop
           endif

!                Centre-
!                 post           
!                  |
!                  |
!                  |
!          -----------------           
!          left-t    right-t

!          Increment to do the centrepost edge
!          Blank edges as encountered. blankcon stores the indices for
!          the new midpoint and the first and last new edges for each
!          edge. Need to do this so we avoid duplicate edges and points
!          later on.           
           nedge=nedge+1
           if(nedge.gt.nedgemax)then
             print *,"Edgemax exceeded",nedge
             stop
           endif
           edge(nedge,1)=centpoint
           if(blanking(ed).eq.0)then
             npt=npt+1
             if(npt.gt.nptmax)then
               print *,"Npt exceeded",npt
               stop
             endif
             grd(npt,:)=0.5*(grd(edge(ed,1),:)+grd(edge(ed,2),:))
             edge(nedge,2)=npt
             blankcon(ed,3)=npt
           endif
           if(blanking(ed).eq.1)then
             edge(nedge,2)=blankcon(ed,3)
           endif
        
           if(k.eq.1)then
             edge(nedge,3)=ibotright
             edge(nedge,4)=nc
           endif
           if(k.eq.2)then
             edge(nedge,3)=itopright
             edge(nedge,4)=ibotright
           endif
           if(k.eq.3)then
             edge(nedge,3)=itopleft
             edge(nedge,4)=itopright
           endif
           if(k.eq.4)then
             edge(nedge,3)=ibotleft
             edge(nedge,4)=itopleft
           endif

!          The centrepost is always perpendicular to the edge           
           if(dir(ed).eq.1) otherdir=2
           if(dir(ed).eq.2) otherdir=1
           dir(nedge)=otherdir

           !write(808,*) grd(edge(nedge,1),1),grd(edge(nedge,1),2)
           !write(808,*) grd(edge(nedge,2),1),grd(edge(nedge,2),2)

!          Only add the left t-edge if the edge has not already been
!          added. Check this by looking at blanking.           
           if(blanking(ed).eq.0)then
             nedge=nedge+1
             if(nedge.gt.nedgemax)then
               print *,"Edgemax exceeded",nedge
               stop
             endif
             edge(nedge,1)=firstpt
             edge(nedge,2)=npt
             edge(nedge,4)=outer
             blankcon(ed,1)=nedge
             
             if(k.eq.1)then
             edge(nedge,3)=nc
             endif
             if(k.eq.2)then
             edge(nedge,3)=ibotright
             endif
             if(k.eq.3)then
             edge(nedge,3)=itopright
             endif
             if(k.eq.4)then
             edge(nedge,3)=itopleft
             endif
             dir(nedge)=dir(ed)
           endif

!          The edge has been added already. Get info from blankcon.           
           if(blanking(ed).eq.1)then
             
             if(k.eq.1)then
             edge(blankcon(ed,2),4)=nc!ibotleft
             endif
             if(k.eq.2)then
             edge(blankcon(ed,2),4)=ibotright
             endif
             if(k.eq.3)then
             edge(blankcon(ed,2),4)=itopright
             endif
             if(k.eq.4)then
             edge(blankcon(ed,2),4)=itopleft
             endif

           endif

           !write(808,*) grd(edge(nedge,1),1),grd(edge(nedge,1),2)
           !write(808,*) grd(edge(nedge,2),1),grd(edge(nedge,2),2)

!          Repeat the process for the right t-edge.           
           if(blanking(ed).eq.0)then
           nedge=nedge+1
           if(nedge.gt.nedgemax)then
             print *,"Edgemax exceeded",nedge
             stop
           endif
           edge(nedge,1)=npt
           edge(nedge,2)=lastpt
           edge(nedge,4)=outer
           blankcon(ed,2)=nedge
           if(k.eq.1)then
             edge(nedge,3)=ibotright
           endif
           if(k.eq.2)then
             edge(nedge,3)=itopright
           endif
           if(k.eq.3)then
             edge(nedge,3)=itopleft
           endif
           if(k.eq.4)then
             edge(nedge,3)=nc!ibotleft
           endif
           dir(nedge)=dir(ed)
           endif

           if(blanking(ed).eq.1)then
             
             if(k.eq.1)then
             edge(blankcon(ed,1),4)=ibotright
             endif
             if(k.eq.2)then
             edge(blankcon(ed,1),4)=itopright
             endif
             if(k.eq.3)then
             edge(blankcon(ed,1),4)=itopleft
             endif
             if(k.eq.4)then
             edge(blankcon(ed,1),4)=nc
             endif

           endif
           

           !write(808,*) grd(edge(nedge,1),1),grd(edge(nedge,1),2)
           !write(808,*) grd(edge(nedge,2),1),grd(edge(nedge,2),2)

           !endif

       !right

!          blank the current edge       
           blanking(ed)=1

!          Traverse away from the current point (known to be
!          anti-clockwise due to pre-loop set-up)           
           possiblenexted=edgept(point,1)
           if(possiblenexted.ne.ed)then
             ed=possiblenexted
           else
             ed=edgept(point,2)
           endif

           if(verb.eq.1) print *,"ed",ed,edgept(point,1),edgept(point,2)
          
         enddo

       !pause  

       enddo

       !pause

!       ncellscurrent1=ncellscurrent2
       ncellscurrent=ncellscurrent2

       fileno=6
       if(verbfile.eq.1)then
       call DUMPMESH(npt,nptmax,nedge,nedgemax,edge,grd,blanking,fileno)
       endif
       if(verb.eq.1) print *,ncellsinit
       ncellsinit=ncellscurrent
       if(verb.eq.1) print *, ncellsinit


       enddo  !ref level
      
      
end subroutine





!     **********************************************************
      subroutine INTERSECT(pt_a,pt_b,ed_a,ed_b,dir,intsct,doesit)

      implicit none
      
      real,intent(in),dimension(2)::pt_a,pt_b,ed_a,ed_b
      integer,intent(in)::dir
      real,intent(out),dimension(2)::intsct
      integer,intent(out)::doesit
      integer::otherdir
      
      real::t 
      
!     Check to see if intersection possible. Allow a matching end point
!     as an intersection.                                            
      doesit=0
      if(((ed_a(dir).le.pt_b(dir)).and.(ed_a(dir).ge.pt_a(dir))).or.    &
     &  (ed_a(dir).ge.pt_b(dir)).and.(ed_a(dir).le.pt_a(dir)))then
        doesit=1
      endif
            
      if(doesit.eq.1)then

!       Don't count an edge parallel to a surface as intersection - the
!       end points of the orthogonal edges will take care of this. Need to
!       mark these edges for removal? CHECK THIS.              
        if(abs(pt_b(dir)-pt_a(dir)).gt.0.0)then
          t=(ed_a(dir)-pt_a(dir))/(pt_b(dir)-pt_a(dir))
          intsct=pt_a+t*(pt_b-pt_a)
        else
          doesit=0
        endif

!       Look in other direction                
        if(dir.eq.1) otherdir=2
        if(dir.eq.2) otherdir=1

!       Check intersection lies on edge                        
        if((intsct(otherdir).le.ed_b(otherdir)).and.                    &
     &    (intsct(otherdir).ge.ed_a(otherdir)).or.                      &
     &    (intsct(otherdir).ge.ed_b(otherdir)).and.                     &
     &    (intsct(otherdir).le.ed_a(otherdir)))then
          doesit=1
        else
          doesit=0
        endif
      
      endif      
            
      end subroutine
!     **********************************************************

      
!     **********************************************************
!      subroutine NCELLS(nedge,edge,ncells)

!      implicit none

!      integer,intent(in)::nedge
!      integer,intent(in),dimension(4)::edge
!      integer,intent(out)::ncells

!      integer::i

!      ncells=0

!      do i=1,nedge
!        left=edge(i,3)
!        right=edge(i,4)
!        if(left.gt.0) then
!          ncells=ncells+1
!        endif        
!        if(right.gt.0) then 
!          ncells=ncells+1
!        endif
!      enddo

!      end subroutine
!     **********************************************************
      
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
      subroutine DUMPMESH(npt,nptmax,nedge,nedgemax,edge,grd,blanking,  &
     &  fileno)

      implicit none

      integer,intent(in)::npt,nptmax,nedge,nedgemax,fileno
      real,intent(in),dimension(nptmax,2)::grd
      integer,intent(in),dimension(nedgemax,4)::edge
      integer,intent(in),dimension(nedgemax)::blanking
      integer::unitno,np,ed,j,edlimit
      character(len=5)::filenochar1
      character(len=5)::filenochar

      unitno=888

      edlimit=0
      do ed=1,nedge
        if(blanking(ed).eq.0) edlimit=edlimit+1
      enddo

      write(filenochar1,'(i2)') fileno
      filenochar(1:2)=filenochar1
      filenochar=adjustl(filenochar)
      filenochar(2:5)='.plt'

      open(unitno,file=filenochar,status='unknown')

      write(unitno,*) 'ZONE t="1" N=',npt,'E=',edlimit
      write(unitno,*) 'DATAPACKING=POINT,'
      write(unitno,*) 'ZONETYPE=FELINESEG'

      do np=1,npt
        write(unitno,*) (grd(np,j),j=1,2)
      enddo
      
      do ed=1,nedge
        if(blanking(ed).eq.0)then
          write(unitno,*) (edge(ed,j),j=1,2)
        endif
      enddo

      end subroutine
!     **********************************************************

!     **********************************************************
      subroutine DUMPMESHNM(npt,nptmax,nedge,nedgemax,edge,grd,blanking,&
     &  IO)

      implicit none

      integer,intent(in)::npt,nptmax,nedge,nedgemax,IO
      real,intent(in),dimension(nptmax,2)::grd
      integer,intent(in),dimension(nedgemax,4)::edge
      integer,intent(in),dimension(nedgemax)::blanking
      integer::unitno,np,ed,j,edlimit
      character(len=5)::filenochar1
      character(len=5)::filenochar

      unitno=IO

      edlimit=0
      do ed=1,nedge
        if(blanking(ed).eq.0) edlimit=edlimit+1
      enddo

      !open(unitno,file=filenochar,status='unknown')

      write(unitno,*) 'ZONE t="1" N=',npt,'E=',edlimit
      write(unitno,*) 'DATAPACKING=POINT,'
      write(unitno,*) 'ZONETYPE=FELINESEG'

      do np=1,npt
        write(unitno,*) (grd(np,j),j=1,2)
      enddo
      
      do ed=1,nedge
        if(blanking(ed).eq.0)then
          write(unitno,*) (edge(ed,j),j=1,2)
        endif
      enddo

      end subroutine
!     **********************************************************

!     **********************************************************
      subroutine DUMPMESHVTK(npt,nptmax,nedge,nedgemax,edge,grd,        &
     &  blanking,fileno)

      implicit none

      integer,intent(in)::npt,nptmax,nedge,nedgemax,fileno
      real,intent(in),dimension(nptmax,2)::grd
      integer,intent(in),dimension(nedgemax,4)::edge
      integer,intent(in),dimension(nedgemax)::blanking
      integer::unitno,np,ed,j,edlimit
      character(len=5)::filenochar1
      character(len=5)::filenochar

      unitno=888

      edlimit=0
      do ed=1,nedge
        if(blanking(ed).eq.0) edlimit=edlimit+1
      enddo

      write(filenochar1,'(i2)') fileno
      filenochar(1:2)=filenochar1
      filenochar=adjustl(filenochar)
      filenochar(2:5)='.vtk'

      open(unitno,file=filenochar,status='unknown') 

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


      end subroutine
!     **********************************************************      

!     **********************************************************
      subroutine DUMPMESHVTKNM(npt,nptmax,nedge,nedgemax,edge,grd,      &
     &  blanking,IO)

      implicit none

      integer,intent(in)::npt,nptmax,nedge,nedgemax,IO
      real,intent(in),dimension(nptmax,2)::grd
      integer,intent(in),dimension(nedgemax,4)::edge
      integer,intent(in),dimension(nedgemax)::blanking
      integer::unitno,np,ed,j,edlimit
      character(len=5)::filenochar1
      character(len=5)::filenochar

      unitno=IO

      edlimit=0
      do ed=1,nedge
        if(blanking(ed).eq.0) edlimit=edlimit+1
      enddo

      !open(unitno,file=filenochar,status='unknown') 

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

      end subroutine
!     **********************************************************


!     **********************************************************      
      subroutine INSERTIONSORT(n,array,intrace,outtrace)

      implicit none

      integer,intent(in)::n
      integer,intent(in),dimension(n)::intrace
      real,intent(in),dimension(n)::array
      real,dimension(n)::sorted
      integer,intent(out),dimension(n)::outtrace

      integer::i,j,ivalue
      real::value

      sorted=array
      outtrace=intrace
      
      do i=1,n
        value=sorted(i)
        ivalue=intrace(i)
        j=i-1
        do while((j.ge.1).and.(sorted(j).gt.value))
          sorted(j+1)=sorted(j)
          outtrace(j+1)=outtrace(j)
          j=j-1
        enddo
        sorted(j+1)=value
        outtrace(j+1)=ivalue
      enddo

      end subroutine
!     **********************************************************

!     **********************************************************
      subroutine INDEGRAPH(npt,nmeetmax,ned,connected,cloops,           &
     &  ncell,accumulate,verb)

      implicit none

      integer,intent(in)::npt,nmeetmax,ned,ncell,verb
      integer,intent(inout)::accumulate
      integer,intent(out)::cloops
      integer,intent(inout),dimension(ned,4)::connected
!      real,intent(in),dimension(npt,2)::grdloc
      !real,intent(out),dimension(maxsplits)::tempvol
      integer::i,isloop,no1,no2,ed,point,startpt
      integer::possiblenexted,possiblenext,j,execs,again
      integer,allocatable,dimension(:,:)::connectpt
      integer,allocatable,dimension(:)::nmeet,tochange,colour
      integer::loopchecks

      !real,dimension(2)::base,pta,ptb

      !connect,
      !allocate(nmeet(npt))
      !allocate(colour(npt))

      allocate(connectpt(npt,nmeetmax))
      allocate(nmeet(npt))
      allocate(colour(npt))
      allocate(tochange(npt))

      loopchecks=0
      colour=0

      tochange=0

      nmeet=0
      connectpt=0
      do i=1,ned
      
        no1=connected(i,1)
        no2=connected(i,2)

        nmeet(no1)=nmeet(no1)+1
        connectpt(no1,nmeet(no1))=i

        nmeet(no2)=nmeet(no2)+1
        connectpt(no2,nmeet(no2))=i

      enddo

      if(verb.eq.1)then
        do i=1,npt
          print *, "argh",i,(connectpt(i,j),j=1,nmeet(i))
        enddo
      endif

      do i=1,npt
        if(nmeet(i).ne.2)then
          print *, "There are not two meeting edges at",i,nmeet(i)
         ! print *, (connectpt(i,j),j=1,nmeet(i))
         ! stop
        endif
      enddo

      !stop

      isloop=0
      ed=1
      startpt=connected(ed,1)
      point=connected(ed,1)

      if(verb.eq.1)then
      print *, "Scanning graph..."
        !print *, ed,startpt
        print *, "***"
      endif

      
      cloops=0
!      area=0.0
      again=1

      do while((again.eq.1).and.(loopchecks.lt.500))

        execs=0
        isloop=0
        tochange=0
        

        do while((isloop.ne.1).and.(execs.lt.500))

          !pta=grdloc(point,:)
        
          possiblenext=connected(ed,1)
          if(possiblenext.ne.point)then
            point=possiblenext
          else
            point=connected(ed,2)
          endif

          if(verb.eq.1) print *, ed,point
          colour(point)=1
          tochange(point)=1

!          base=0.0
!          ptb=grdloc(point,:)

!          areainc=0.5*((pta(1)-base(1))*(ptb(2)-base(2))-               &
!     &      (pta(2)-base(2))*(ptb(1)-base(1)))
!          area=area+areainc

          possiblenexted=connectpt(point,1)
          if(possiblenexted.ne.ed)then
            ed=possiblenexted
          else
            ed=connectpt(point,2)
          endif
        
          if(point.eq.startpt)then
            isloop=1
            cloops=cloops+1
            if(cloops.gt.1)then
              accumulate=accumulate+1
            endif
            !if(verb.eq.1) print *, "Loop found",area
            !tempvol(cloops)=abs(area)
          endif

          execs=execs+1

        enddo

        again=0
        do i=1,npt
          if(colour(i).ne.1)then
            again=1
            !print *, "Checking for another loop",i
            startpt=i
            point=i
            ed=connectpt(point,1)
            !area=0.0
          endif
        enddo


        do i=1,ned
          if((tochange(connected(i,1)).eq.1).and.(cloops.gt.1))then
                 ! print *,'adjusting counter'
            if(connected(i,3).eq.ncell)then
              !accumulate=accumulate+1
              connected(i,3)=accumulate
              if(verb.eq.1) print *, "left",accumulate
            endif
            if(connected(i,4).eq.ncell)then
              !accumulate=accumulate+1
              connected(i,4)=accumulate
              if(verb.eq.1) print *, "right",accumulate
            endif
          endif
        enddo

        if(again.eq.1)then
          if(verb.eq.1) print *, "Checking for another loop"
          loopchecks=loopchecks+1
          if(loopchecks.gt.10) then
            !print *, "Max checks exceeded!!!"
            !pause
          endif
        endif
        !pause  

      enddo

      if(verb.eq.1) print *, "Total closed loops",cloops
      if(cloops.gt.1)then
        !pause
      endif

      deallocate(connectpt,nmeet,colour,tochange)

      end subroutine
!     **********************************************************




!      OPTIONS(CHECK,FULL_DEBUG)
!      OPTIONS(CHECK)
!      OPTIONS(OPTIMISE)
!     ***********************************************************
!                           EulerFlow
!
!                 Finite Volume 2D Euler solver 
!                     for unstructured grids
!   
!                     Thomas C.S. Rendall
!                  University of Bristol, U.K.
!                         
!              1st Distributed Version 1.0 of 22/07/10
!
!     ***********************************************************

!     ***********************************************************   
      program E1J
     
      implicit none

!     Local timestepping
!     Characteristic boundary conditions      
!     No multigrid 
!     Not parallel
!     Steady/unsteady  
!     JST unstructured dissipation      
!                                                   
!     ***********************************************************            

!     REAL ARRAYS
      real,allocatable,dimension(:)::rho,u,v,E,p,rhohold
      real,allocatable,dimension(:,:)::grdnp1,grdn,grdnm1,Wnp1,Wn,      &
     &  Wnm1,fvel,cf,grdh,cfh
      double precision,allocatable,dimension(:,:)::cfmat,inv_cfmat
      real,allocatable,dimension(:)::volnp1,voln,volnm1,deltat,mach,    &
     &  lam,k
      
!     REAL VARS
      real::CFL,Runiv,machinf,rhoinf,Tinf,gam,aoa,pinf,cy,cx,cl,cm,cd
      real::res,k1,k2,pi,cp,absspd,fnp1,fn,fnm1,dt,sr,ang,uinf
      real::rfreq,chord,omega,momx,momy,momxh,momyh,xpit,ypit,amp,      &
     &  amp2,disp

!     INTEGER VARS
      integer::t,ncell,nedge,nvert,nc,iunflag,                          &
     &  nreal,nrealstop,tstopi,nsurf,rbf_type,irest,dimnum,npercyc,nv
      integer,allocatable,dimension(:)::dsf,m,tstop,INDX,list
      integer,allocatable,dimension(:,:)::edge
      
      pi=4.0*atan(1.0)
      
      print *,'********************************************************'
      print *,'                        EulerFlow                       '
      print *,'                                                        '
      print *,'               Finite Volume 2D Euler solver            '
      print *,'                  for unstructured grids                '
      print *,'                                                        '
      print *,'                  Thomas C.S. Rendall                   '
      print *,'               University of Bristol, U.K.              '
      print *,'                                                        '
      print *,'          1st Distributed Version 1.0 of 22/07/10       '
      print *,'                                                        '
      print *,'********************************************************'
      
      call READINT(ncell,nedge,nvert)
      
      print *, " "      
      print *, "Sizes read..."
      dimnum=2
      open(256,file='refpoints',status='unknown')
        read(256,*) momxh,momyh
      close(256)
      !momxh=0.25
      !momyh=0.0
      
      allocate(rho(ncell))
      allocate(rhohold(ncell))
      allocate(u(ncell))
      allocate(v(ncell))
      allocate(E(ncell))
      allocate(p(ncell))
      allocate(mach(ncell))
      allocate(lam(ncell))      
                        
      allocate(grdnp1(nvert,2))
      allocate(grdn(nvert,2))
      allocate(grdnm1(nvert,2))
      allocate(grdh(nvert,2))
      allocate(volnp1(ncell))
      allocate(voln(ncell))
      allocate(volnm1(ncell))
      allocate(k(dimnum))
      allocate(deltat(ncell))
      allocate(m(ncell))
                             
      allocate(Wnp1(ncell,4))
      allocate(Wn(ncell,4))
      allocate(Wnm1(ncell,4))
      allocate(fvel(nedge,2))
                  
      allocate(edge(nedge,4))
      
      allocate(dsf(4))
      
      print *, "Arrays allocated..."
      
      call READSETTINGS(ncell,nedge,nvert,                              &
     &             Runiv,machinf,rhoinf,Tinf,gam,aoa,tstopi,CFL,grdnp1  &
     &                  ,edge,volnp1,k1,k2,dsf,m,iunflag,irest,         &
     &                   rfreq,nrealstop,npercyc,chord,rbf_type,sr)

      uinf=machinf*(gam*Runiv*Tinf)**0.5
      
      if(iunflag.eq.1)then
        print *, "Reading forced pitching file..."
        open(700,file='forced',status='unknown')
        read(700,*) xpit,ypit
        read(700,*) amp,amp2
        close(700)
      else
        xpit=0.0
        ypit=0.0
        amp=0.0
      endif    
      print *, "a"
      call CELLN(nedge,nvert,edge,nsurf,grdnp1)
      
      print *, " "
      print *, "Cells:",ncell     
      print *, "Edges:",nedge
      print *, "Vertices:",nvert
      print *, "Surface points:",nsurf
     
      voln=volnp1
      volnm1=volnp1
      grdn=grdnp1
      grdnm1=grdnp1
      grdh=grdnp1
      
      allocate(tstop(nrealstop))
      allocate(cf(nsurf,2))
      allocate(cfh(nsurf,2))
      allocate(cfmat(nsurf,nsurf))
      allocate(inv_cfmat(nsurf,nsurf))
      allocate(INDX(nsurf))
      
      call CELLFN(nedge,nvert,nsurf,grdnp1,edge,cf)
      !(nedge,nvert,nsurf,grd,edge,cf)
      cfh=cf
      k=1.0       
      if(iunflag.eq.1)then
        call WRITEM(nsurf,dimnum,cfh,rbf_type,sr,cfmat)        
        call MIGS(cfmat,nsurf,inv_cfmat,INDX)
        !print *, "a"
      endif
  
      if(iunflag.eq.1)then
      
        tstop(:)=tstopi

        omega=rfreq*2.0*uinf/chord
        dt=2.0*pi/(npercyc*omega)
        
        ! First order 
        !fnp1=1.0/dt
        !fn=-1.0/dt
        !fnm1=0.0
        
        !print *,dt,"dt" 
        ! Second order
        fnp1=1.5/dt
        fn=-2.0/dt
        fnm1=0.5/dt

        !fnp1=0.0
        !fn=0.0
        !fnm1=0.0
              
      else 
      
        tstop(:)=tstopi
        omega=1.0
        nrealstop=1
        dt=1e6
        fnp1=0.0
        fnm1=0.0
        fn=0.0
        ang=0.0
        
      endif
            
      print *, '  '
      print *, 'Flow variables:'
      print *, "-------------------------------------------------------"
      print *, "Mach number:",machinf
      print *, "Density:",rhoinf
      print *, "Temperature:",Tinf
      print *, "Gamma:", gam
      print *, "AoA:",aoa*180.0/pi
      print *, "Maxit:",tstopi
      print *, "-------------------------------------------------------"
      print *, '  '
      
      if(irest.eq.0)then
        call INITIALISE(ncell,Runiv,machinf,rhoinf,Tinf,gam,            &
     &                      aoa,                                        &
     &                      rho,u,v,E,p)
      else
        open(778,file='flow.dump',status='unknown')
        do nc=1,ncell
          read(778,*) rho(nc),u(nc),v(nc),E(nc),p(nc)
        enddo
        close(778)
      endif
                    
      call CON(ncell,rho,u,v,E,Wnp1)
      Wn=Wnp1
      Wnm1=Wnp1
                                              
      open(501,file='res_hist.dat',status='unknown')
      open(401,file='res.plt',status='unknown')  
      open(400,file='cp.dat',status='unknown')
      open(600,file='unsteady_hist.plt',status='unknown')
      open(776,file='coords.dat',status='unknown')
      open(777,file='flow.dat',status='unknown')
      write(501,*) "Iteration   |   Log(Res)  | cl | cm | cd | cx | cy "
      write(501,*) "---------------------------------------------------"
                                            
      print *, "Beginning time-stepping..."
      
200   format("Iteration   |   Log(Res)   |       Cl      |       Cd   ")
300   format("--------------------------------------------------------")
                                                                  
      do nreal=1,nrealstop
        
        if(iunflag.eq.1)then
        
          print *, " "
          print *, "               Real timestep:",nreal
                      
          write(6,200)
          write(6,300)
                                  
          call MOVESURF(nsurf,nreal,omega,xpit,ypit,amp,amp2,dt,cf,     &
     &  cfh,ang,disp,momx,momy,momxh,momyh)
          
          call MOVEGRID(nvert,nsurf,sr,rbf_type,inv_cfmat,grdh,         &
     &  cfh,grdnp1,cf)
          
          call EDGEVEL(nedge,nvert,edge,fnp1,fn,fnm1,grdnp1,grdn,       &
     &  grdnm1,fvel)
          
          call GCL(ncell,nedge,nvert,edge,grdn,fvel,fnp1,fn,fnm1,       &
     &        voln,volnm1,volnp1)
          
          do nc=1,ncell
            if(volnp1(nc).lt.0.0)then
              print*,volnp1(nc),"vol neg"
              stop
            endif
          enddo
          
        else
        
          fvel=0.0
          write(6,200)
          write(6,300)
          momx=momxh
          momy=momyh
          ang=0.0
          disp=0.0
                                    
        endif
      
        do t=1,tstop(nreal)

        rhohold=rho
         
          !print *, 'tstep' 
        call  TSTEP(ncell,nvert,nedge,edge,grdnp1,CFL,gam,rho,u,v,p,    &
     &                 volnp1,deltat,lam,dt)            
         
      !    print *, 'jam'      
        call JAMESON(ncell,nedge,nvert,edge,grdnp1,rho,u,v,E,p,         &
     &                   Runiv,machinf,rhoinf,Tinf,gam,aoa,volnp1,voln, &
     &              volnm1,Wnp1,Wn,Wnm1,fnp1,fn,fnm1,fvel,deltat,       &
     &                   k1,k2,dsf,iunflag,m,lam)      
                                               
          !print *, 'resid'               
        call RESID(ncell,nvert,nedge,grdnp1,edge,rho,rhohold,deltat,    &        
     &            res,machinf,p,gam,rhoinf,Runiv,Tinf,t,tstop(nreal),   &
     &                 momx,momy,cx,cy,cm)

        cl=cy*cos(aoa+ang)-cx*sin(aoa+ang)
        cd=cx*cos(aoa+ang)+cl*sin(aoa+ang)
        
        if(res.gt.0.0)then
          print *,t,log10(res),cl,cd        
          write(501,*) t,log10(res),cl,cm,cd,cx,cy
          write(401,*) t,log10(res)!,cl,cm,cd,cx,cy
        else
          print *,t,res,cl,cd        
          write(501,*) t,res,cl,cm,cd,cx,cy
          write(401,*) t,res!,cl,cm,cd,cx,cy
        endif
     
        enddo

        write(600,*) ang,disp,cl,cm

        do nv=1,nvert
          write(776,*) grdnp1(nv,1),grdnp1(nv,2)
        enddo
        
        !if(iunflag.gt.1)then
        grdnm1=grdn
        grdn=grdnp1
        
        volnm1=voln
        voln=volnp1
        
        Wnm1=Wn
        Wn=Wnp1

        pinf=rhoinf*Runiv*Tinf
        do nc=1,ncell
          absspd=sqrt(u(nc)**2+v(nc)**2)
          mach(nc)=absspd/(sqrt(gam*p(nc)/rho(nc)))
        enddo
        do nc=1,ncell
          cp=(2/(gam*machinf**2))*(p(nc)/pinf-1)
          write(777,*) rho(nc),u(nc),v(nc),mach(nc),cp
        enddo 
        !endif
        !pause
        
      enddo

      close(776)
      close(777)
      
      open(778,file='flow.dump',status='unknown')
      do nc=1,ncell
        write(778,*) rho(nc),u(nc),v(nc),E(nc),p(nc)
      enddo
      close(778)      
      
      close(501)
      close(401)
      close(400)
      close(600)
      
      deallocate(rho)
      deallocate(rhohold)
      deallocate(u)
      deallocate(v)
      deallocate(E)
      deallocate(p)
      deallocate(mach)
      deallocate(lam)       
                        
      deallocate(grdnp1)
      deallocate(grdn)
      deallocate(grdnm1)
      deallocate(grdh)
      deallocate(volnp1)
      deallocate(voln)
      deallocate(volnm1)
      deallocate(k)
      deallocate(deltat)
      deallocate(m)
                             
      deallocate(Wnp1)
      deallocate(Wn)
      deallocate(Wnm1)
      deallocate(fvel)
      
      deallocate(tstop)
      deallocate(cf)
      deallocate(cfh)
      deallocate(cfmat)      
      
      end program E1J
!     **********************************************************        




!     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!     FUNCTIONS
!     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

!     ***********************************************************
      real function RBF(ident,radius,dist)

      implicit none

!     Function for RBF evluation.

      integer, intent (in):: ident 
      real, intent(in):: dist,radius
      real::e,pi!,norm2d,norm3d

      e=2.718281828
      pi=4.0*atan(1.0)
     
      if (dist.lt.0.0)then
        print *, "Negative distance-code error"
      endif

!     Gaussian (G) (Not compact, (Not compact, 
!     'radius' here is the decay parameter)
      if (ident.eq.0)then
        RBF=e**(-radius*dist)
      endif  
        
!     Thin Plate Spline (TPS) (Not compact)
!     Set to RBF=0.0 for dist=0.0, avoid ln 0
!     error
      if (ident.eq.1)then
        if (dist.gt.0.0)then      
          RBF=(dist**2)*log(dist)
        endif
        if (dist.le.0.0)then
          RBF=0.0
        endif
      endif

!     Hardy's Multiquadric (HMQ) (Not compact, 
!     'radius' here is the shape parameter)
 
      if (ident.eq.2)then
        RBF=(radius**2+dist**2)**0.5
      endif        
        
!     Hardy's Inverse Multiquadric (HIMQ) (Not compact, 
!     'radius' here is the shape parameter)
      if (ident.eq.3)then
        RBF=1/((radius**2+(dist)**2)**0.5)
      endif

!     Wendland's C0 (C0, compact)
      if (ident.eq.4)then
        RBF=(1-dist/radius)**2
        if ((dist/radius).gt.1.0)then
          RBF=0.0
        endif
      endif

!     Wendland's C2 (C2, compact)
      if (ident.eq.5)then
        RBF=((1-(dist/radius))**4)*(4*(dist/radius)+1)
        if ((dist/radius).gt.1.0)then
          RBF=0.0
        endif
      endif

!     Wendland's C4 (C4, compact)
      if (ident.eq.6)then
        RBF=((1-(dist/radius))**6)*(35*(dist/radius)**2+                &
     &  18*(dist/radius)+3)
        if ((dist/radius).gt.1.0)then
          RBF=0.0
        endif
      endif

!     Wendland's C6 (C6, compact)
      if (ident.eq.7)then
        RBF=((1-(dist/radius))**8)*(32*(dist/radius)**3+                &
     &  25*(dist/radius)**2+8*(dist/radius)**2+1)
        if ((dist/radius).gt.1.0)then
          RBF=0.0
        endif
      endif
        
!     Euclid's Hat Function (EHF, compact, compact for 2*radius)
      if (ident.eq.8)then
        RBF=pi*((1/12.0)*dist**3-(dist*radius**2)+(4.0/3.0)*radius**3)
        if ((dist/radius).gt.2.0)then
          RBF=0.0
        endif
      endif

      end function RBF
!     ***********************************************************



!     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!     SUBROUTINES
!     &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

!     *****************************************************
      subroutine DISS(ncell,nedge,nvert,W,edge,rho,u,v,p,gam,grd,       &
     &  k1,k2,D,m,lam,nstage,vol)                                           
      
      implicit none
      
      integer,intent(in)::ncell,nedge,nvert,nstage
      real,intent(in),dimension(ncell,4)::W
      real,intent(in),dimension(ncell)::vol
      real,intent(in),dimension(ncell)::rho,u,v,p
      integer,intent(in),dimension(nedge,4)::edge
      real,intent(in)::k1,k2,gam
      real,intent(in),dimension(nvert,2)::grd
      real,intent(out),dimension(ncell,4)::D
      !integer,dimension(ncell)::m
      integer,intent(in),dimension(ncell)::m
      real,intent(inout),dimension(ncell)::lam
      
      integer::ed,nc
      real::alp1,alp2,p3,p4,u01x,u01y,c01,S01,dx,dy,psi1,psi0           &
     &  ,psi01,lam01
      real,dimension(ncell,4)::LAP
      real::s2,s4
      real,dimension(ncell)::pdss,pdssdf,pdsssm
      
      real,dimension(4)::W3,W4
      real::ul,ur,vl,vr,cl,cr
      !real,dimension(ncell)::psi
      
      LAP=0.0                 
      
      do ed=1,nedge      
        if((edge(ed,3).gt.0).and.(edge(ed,4).gt.0)) then
        
          W3=W(edge(ed,3),:)
          W4=W(edge(ed,4),:)

          dx=grd(edge(ed,2),1)-grd(edge(ed,1),1)
          dy=grd(edge(ed,2),2)-grd(edge(ed,1),2)
          S01=sqrt(dx**2+dy**2)
        
          LAP(edge(ed,3),:)=LAP(edge(ed,3),:)+(W4-W3)!*S01
          LAP(edge(ed,4),:)=LAP(edge(ed,4),:)-(W4-W3)!*S01
         
        endif           
      enddo

      !do nc=1,ncell
      !  LAP(nc)=LAP(nc)/vol(nc)
      !enddo
      
      pdss=0.0
      pdssdf=0.0
      pdsssm=0.0
      
      do ed=1,nedge      
        if((edge(ed,3).gt.0).and.(edge(ed,4).gt.0)) then
        
          p3=p(edge(ed,3))
          p4=p(edge(ed,4))
        
          pdssdf(edge(ed,3))=pdssdf(edge(ed,3))+(p4-p3)!/(p4+p3)
          pdssdf(edge(ed,4))=pdssdf(edge(ed,4))-(p4-p3)!/(p4+p3)
          
          pdsssm(edge(ed,3))=pdsssm(edge(ed,3))+(p4+p3)
          pdsssm(edge(ed,4))=pdsssm(edge(ed,4))+(p4+p3)          
          
        endif           
      enddo
      
      do ed=1,nedge      
        if((edge(ed,3).gt.0).and.(edge(ed,4).gt.0)) then
        
        pdss(edge(ed,3))=abs(pdssdf(edge(ed,3)))/(pdsssm(edge(ed,3)))
        pdss(edge(ed,4))=abs(pdssdf(edge(ed,4)))/(pdsssm(edge(ed,4)))
          
        endif           
      enddo

      if(nstage.lt.4)then
              
        lam=0.0

        do ed=1,nedge      
        if(edge(ed,3).gt.0)then
          cl=sqrt((gam*p(edge(ed,3))/rho(edge(ed,3))))
          ul=u(edge(ed,3))
          vl=v(edge(ed,3))
        else
          cl=sqrt((gam*p(edge(ed,4))/rho(edge(ed,4))))
          ul=u(edge(ed,4))
          vl=v(edge(ed,4))
        endif
        
        if(edge(ed,4).gt.0)then
          cr=sqrt((gam*p(edge(ed,4))/rho(edge(ed,4))))
          ur=u(edge(ed,4))
          vr=v(edge(ed,4))
        else
          cr=sqrt((gam*p(edge(ed,3))/rho(edge(ed,3))))
          ur=u(edge(ed,3))
          vr=v(edge(ed,3))
        endif       

        u01x=0.5*(ul+ur)
        u01y=0.5*(vl+vr)
        c01=0.5*(cl+cr)
        
        dx=grd(edge(ed,2),1)-grd(edge(ed,1),1)
        dy=grd(edge(ed,2),2)-grd(edge(ed,1),2)
        S01=sqrt(dx**2+dy**2)
        dx=dx/S01
        dy=dy/S01
            
        if((edge(ed,3).gt.0))then
          lam(edge(ed,3))=lam(edge(ed,3))+                              &
     &    (abs(dy*u01x-dx*u01y)+c01)*S01
        endif
        
        if((edge(ed,4).gt.0))then
          lam(edge(ed,4))=lam(edge(ed,4))+                              &
     &    (abs(-dy*u01x+dx*u01y)+c01)*S01
        endif

        enddo
      endif      
                        
      D=0.0
      
      do ed=1,nedge      
        if((edge(ed,3).gt.0).and.(edge(ed,4).gt.0)) then
        
          W3=W(edge(ed,3),:)
          W4=W(edge(ed,4),:)
          
          c01=0.5*(sqrt(gam*p(edge(ed,3))/rho(edge(ed,3)))+            &
     &   sqrt(gam*p(edge(ed,4))/rho(edge(ed,4))))
     
          u01x=0.5*(u(edge(ed,3))+u(edge(ed,4)))
          u01y=0.5*(v(edge(ed,3))+v(edge(ed,4)))
        
          dx=grd(edge(ed,2),1)-grd(edge(ed,1),1)
          dy=grd(edge(ed,2),2)-grd(edge(ed,1),2)
          S01=sqrt(dx**2+dy**2)
          dx=dx/S01
          dy=dy/S01
          
          lam01=(abs(dy*u01x-dx*u01y)+c01)*S01
          psi0=lam(edge(ed,3))/(4*lam01)
          psi1=lam(edge(ed,4))/(4*lam01)          
          psi01=4*psi1*psi0/(psi0+psi1)
          if(psi01.le.0.0)then
            print *, 'warning no diss edge',ed
            stop
          endif
          !if(psi01.eq.0.0)then
          !  print *, 'warning no diss edge',ed
          !endif 

          s2=3.0*((m(edge(ed,3))+m(edge(ed,4))))/                       &
     &     (m(edge(ed,3))*m(edge(ed,4)))  
          s4=(s2**2)/4.0        
          
          alp1=k1*0.5*(pdss(edge(ed,3))+pdss(edge(ed,4)))*s2      
          alp2=max(0.0,k2-alp1)*s4
          !alp2=k2        
                  
          D(edge(ed,3),:)=D(edge(ed,3),:)+((alp1*(W3(:)-W4(:))-         &
     &    alp2*(LAP(edge(ed,3),:)-LAP(edge(ed,4),:))))*psi01*lam01
     
          D(edge(ed,4),:)=D(edge(ed,4),:)-((alp1*(W3(:)-W4(:))-         &
     &    alp2*(LAP(edge(ed,3),:)-LAP(edge(ed,4),:))))*psi01*lam01
        
        endif        
      enddo
      !d=0.0      
      
      end subroutine DISS
!     *****************************************************

!     *****************************************************            
      subroutine FLUX(ncell,nedge,nvert,edge,grd,rho,u,v,E,p,fvel,R,    &
     &  Runiv,machinf,rhoinf,Tinf,gam,aoa)
      
      implicit none
      
      integer,intent(in)::ncell,nedge,nvert
      real,intent(in),dimension(ncell)::rho,u,v,E,p
      real,intent(in),dimension(nedge,2)::fvel
      integer,intent(in),dimension(nedge,4)::edge  
      real,intent(out),dimension(ncell,4)::R
      real,intent(in)::Runiv,machinf,rhoinf,Tinf,gam,aoa
      real,intent(in),dimension(nvert,2)::grd
      
      real,dimension(4)::F3,F4,G3,G4,FLX,FLXF,FLXG
      real,dimension(ncell)::mach
      real::absspd,dy,dx,S01,vnorm
      integer::ed,nc
       
      !mods=0 
      !angle=0.0
      !print *, 'influx'
      !F=0.0
      !G=0.0
      
      do nc=1,ncell
        absspd=sqrt(u(nc)**2+v(nc)**2)
        if(((p(nc)).le.0.0).or.(rho(nc).le.0.0))then
        print *, "neg sound"
        stop
        endif
        mach(nc)=absspd/(sqrt(gam*p(nc)/rho(nc)))
      enddo
      !print *, 'f1',ncell
      
      R=0.0
      FLX=0.0
      FLXF=0.0
      FLXG=0.0
                          
      do ed=1,nedge
        
                                    
        F3=0.0
        G3=0.0
        F4=0.0
        G4=0.0
        FLX=0.0
        FLXF=0.0
        FLXG=0.0

                
!       aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa              
        if(edge(ed,3).gt.0) then
        
          call GETFG(rho(edge(ed,3)),u(edge(ed,3)),                     &
     &     v(edge(ed,3)),E(edge(ed,3)),p(edge(ed,3)),fvel(ed,:)         &
     &     ,F3,G3)
               
        endif

        if(edge(ed,4).gt.0)then
        
          call GETFG(rho(edge(ed,4)),u(edge(ed,4)),                     &
     &     v(edge(ed,4)),E(edge(ed,4)),p(edge(ed,4)),fvel(ed,:)         &
     &     ,F4,G4)
     
        endif
!       aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa                


!       bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb                
        if(edge(ed,3).eq.-1)then  !wall
          !print *, 'wall'
          call WALL(F3,G3,p(edge(ed,4)),fvel(ed,:))
          !This avoids the averaging at boundaries - 0.5*(a+a)=a etc.
          F4=F3
          G4=G3
        endif
        
        if(edge(ed,3).eq.-2)then  !Farfield
                
          dx=grd(edge(ed,2),1)-grd(edge(ed,1),1)
          dy=grd(edge(ed,2),2)-grd(edge(ed,1),2)
          S01=sqrt(dx**2+dy**2)
          dx=dx/S01
          dy=dy/S01
          vnorm=dy*u(edge(ed,4))-dx*v(edge(ed,4))

          if(vnorm.ge.0.0)then  !inflow          
            if(mach(edge(ed,4)).ge.1.0)then      
              call INFLOWSUP(F3,G3,Runiv,machinf,rhoinf,Tinf,gam,aoa)
            endif
            if(mach(edge(ed,4)).lt.1.0)then
              call INFLOWSUB(F3,G3,rho(edge(ed,4)),u(edge(ed,4)),       &
     &  v(edge(ed,4)),p(edge(ed,4)),Runiv,machinf,                      &
     &  rhoinf,Tinf,gam,aoa)
            endif
          endif

          if(vnorm.lt.0.0)then  !outflow
            if(mach(edge(ed,4)).ge.1.0)then            
              call GETFG(rho(edge(ed,4)),u(edge(ed,4)),                 &
     &     v(edge(ed,4)),E(edge(ed,4)),p(edge(ed,4)),fvel(ed,:),F4,G4)  
              call OUTFLOWSUP(F3,G3,F4,G4) 
            endif
            if(mach(edge(ed,4)).lt.1.0)then
              call OUTFLOWSUB(F3,G3,rho(edge(ed,4)),u(edge(ed,4)),      &
     &  v(edge(ed,4)),p(edge(ed,4)),Runiv,machinf,                      &
     &  rhoinf,Tinf,gam,aoa)
            endif      
          endif
           
        endif  ! end farfield 
!       bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb        


!       cccccccccccccccccccccccccccccccccccccccccccccccccccccccc                
        if(edge(ed,4).eq.-1)then  !wall
          call WALL(F4,G4,p(edge(ed,3)),fvel(ed,:))
          !This avoids the averaging at boundaries - 0.5*(a+a)=a etc.
          F3=F4
          G3=G4
        endif
        
        if(edge(ed,4).eq.-2)then  !Farfield

          dx=grd(edge(ed,2),1)-grd(edge(ed,1),1)
          dy=grd(edge(ed,2),2)-grd(edge(ed,1),2)
          S01=sqrt(dx**2+dy**2)
          dx=dx/S01
          dy=dy/S01
          vnorm=dy*u(edge(ed,3))-dx*v(edge(ed,3))      
          
          if(vnorm.le.0.0)then      
            if(mach(edge(ed,3)).ge.1.0)then
              call INFLOWSUP(F4,G4,Runiv,machinf,rhoinf,Tinf,gam,aoa)
            endif
            if(mach(edge(ed,3)).lt.1.0)then
              call INFLOWSUB(F4,G4,rho(edge(ed,3)),u(edge(ed,3)),       &
     &  v(edge(ed,3)),p(edge(ed,3)),Runiv,machinf,                      &
     &  rhoinf,Tinf,gam,aoa)
            endif
          endif
          
          if(vnorm.gt.0.0)then
            if(mach(edge(ed,3)).ge.1.0)then            
              call GETFG(rho(edge(ed,3)),u(edge(ed,3)),                 &
     &     v(edge(ed,3)),E(edge(ed,3)),p(edge(ed,3)),fvel(ed,:),F3,G3)  
              call OUTFLOWSUP(F4,G4,F3,G3)
            endif
          !F3=F4
          !G3=G4
            if(mach(edge(ed,3)).lt.1.0)then
              call OUTFLOWSUB(F4,G4,rho(edge(ed,3)),u(edge(ed,3)),      &
     &  v(edge(ed,3)),p(edge(ed,3)),Runiv,machinf,                      &
     &  rhoinf,Tinf,gam,aoa)
            endif           
          endif        
          
        endif !end farfield     
!       cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                                

!       dddddddddddddddddddddddddddddddddddddddddddddddddddddddd                                    
        FLXF(:)=0.5*(F3(:)+F4(:))*                                      &
     &    (grd(edge(ed,2),2)-grd(edge(ed,1),2))
        FLXG(:)=0.5*(G3(:)+G4(:))*                                      &
     &   (grd(edge(ed,2),1)-grd(edge(ed,1),1))
             
        FLX(:)=FLXF(:)-FLXG(:)                        
                                        
        if(edge(ed,3).gt.0)then
          R(edge(ed,3),:)=R(edge(ed,3),:)+FLX(:)
        endif
        if(edge(ed,4).gt.0)then
          R(edge(ed,4),:)=R(edge(ed,4),:)-FLX(:)
        endif
!       dddddddddddddddddddddddddddddddddddddddddddddddddddddddd
             
                                                   
      enddo
      
      end subroutine FLUX
!     *****************************************************
      
!     *****************************************************
      subroutine GETFG(rho,u,v,E,p,fvel,F,G)
      
      implicit none
      
      real,intent(in)::rho,u,v,E,p
      real,intent(in),dimension(2)::fvel
      real,intent(out),dimension(4)::F,G
      
      if(rho.le.0)then
        print *, "Negative or zero density"
        stop
      endif
        
      F(1)=rho*(u-fvel(1))
      F(2)=rho*u*(u-fvel(1))+p
      F(3)=rho*v*(u-fvel(1))
      F(4)=rho*(u-fvel(1))*(E)+u*(p)
        
      G(1)=rho*(v-fvel(2))
      G(2)=rho*u*(v-fvel(2))
      G(3)=rho*v*(v-fvel(2))+p
      G(4)=rho*(v-fvel(2))*(E)+v*(p)
      
      end subroutine GETFG
!     **********************************************************

!     **********************************************************
      subroutine WALL(Fx,Gx,p,vel)
      
      implicit none
      
      real,intent(out),dimension(4)::Fx,Gx
      real,intent(in)::p
      real,intent(in),dimension(2)::vel
      
      Fx(1)=0.0
      Fx(2)=1.0*p
      Fx(3)=0.0
      Fx(4)=1.0*p*vel(1)
      
      Gx(1)=0.0
      Gx(2)=0.0
      Gx(3)=1.0*p
      Gx(4)=1.0*p*vel(2)
      
      end subroutine WALL
!     ***********************************************************

!     ***********************************************************
      subroutine INFLOWSUB(Fx,Gx,rho1,u1,v1,p1,Runiv,machinf,           &
     &  rhoinf,Tinf,gam,aoa)
            
      implicit none
      
      real,intent(in)::Runiv,machinf,rhoinf,Tinf,gam,aoa
      real,intent(in)::rho1,u1,v1,p1
      real,intent(out),dimension(4)::Fx,Gx
      real,dimension(4,4)::mat,matinv
      real,dimension(4)::pvveco,pvvecn,cvvec
      real,dimension(2)::fsvel,vel,velc
      !real,dimension(2,2)::rmatinv
      
      real::rhoc,uc,vc,pc
      real::Hc,Ec,signv,signu

      if(u1.gt.0.0)then
        signu=1.0
      else
        signu=-1.0
      endif
      !if(v1.gt.0.0)then
      !  signv=1.0
      !else
      !  signv=-1.0
      !endif

!     Set FF and local velocity vectors.            
      fsvel(1)=(machinf*sqrt(gam*Runiv*Tinf))*cos(aoa)
      fsvel(2)=(machinf*sqrt(gam*Runiv*Tinf))*sin(aoa)
      
      vel(1)=u1
      vel(2)=v1

!     Find rotation and inverse rotation matrices.            
!      call ROTMAT(rmat,rmatinv,angle)
 
!     Rotate.            
!      fsvelrot=matmul(rmat,fsvel)
!      velrot=matmul(rmat,vel)      

!     Get matrix to convert to and from char vars.            
      call CHARM(mat,matinv,rhoinf,rhoinf*Runiv*Tinf)

!     Get deltas.            
      pvveco(1)=rho1-rhoinf
      pvveco(2)=signu*(vel(1)-fsvel(1))
      pvveco(3)=signu*(vel(2)-fsvel(2))
      pvveco(4)=p1-rhoinf*Runiv*Tinf

!     Convert and set char vars appropriately.            
      cvvec=matmul(mat,pvveco)      
      cvvec(1)=0.0
      cvvec(2)=0.0
      cvvec(3)=0.0

!     Convert back to primitives.            
      pvvecn=matmul(matinv,cvvec)
      
      rhoc=pvvecn(1)+rhoinf
      velc(1)=signu*pvvecn(2)+fsvel(1)
      velc(2)=signu*pvvecn(3)+fsvel(2)
      pc=pvvecn(4)+rhoinf*Runiv*Tinf

!     Rotate back velocity.            
!      velc=matmul(rmatinv,velc)

      uc=velc(1)
      vc=velc(2)
      
!     Get fluxes.            
      Ec=pc/(rhoc*(gam-1))+0.5*(uc**2+vc**2)
      Hc=Ec+pc/rhoc
      
      Fx(1)=rhoc*uc
      Fx(2)=rhoc*uc**2+pc
      Fx(3)=rhoc*uc*vc
      Fx(4)=rhoc*uc*Hc
      
      Gx(1)=rhoc*vc
      Gx(2)=rhoc*uc*vc
      Gx(3)=rhoc*vc**2+pc
      Gx(4)=rhoc*vc*Hc
      
      end subroutine INFLOWSUB
!     ***********************************************************

!     ***********************************************************
      subroutine OUTFLOWSUB(Fx,Gx,rho1,u1,v1,p1,Runiv,machinf,          &
     &  rhoinf,Tinf,gam,aoa)
            
      implicit none
      
      real,intent(in)::Runiv,machinf,rhoinf,Tinf,gam,aoa
      real,intent(in)::rho1,u1,v1,p1
      real,intent(out),dimension(4)::Fx,Gx
      real,dimension(4,4)::mat,matinv
      real,dimension(4)::pvveco,pvvecn,cvvec
      real,dimension(2)::fsvel,vel,velc
      !real,dimension(2,2)::rmat,rmatinv
      
      !real::uinf,pinf,vinf
      real::rhoc,uc,vc,pc
      real::Hc,Ec,signv,signu

      if(u1.gt.0.0)then
        signu=1.0
      else
        signu=-1.0
      endif
      !if(v1.gt.0.0)then
      !  signv=1.0
      !else
      !  signv=-1.0
      !endif

!     Set FF and local velocity vectors.            
      fsvel(1)=(machinf*sqrt(gam*Runiv*Tinf))*cos(aoa)
      fsvel(2)=(machinf*sqrt(gam*Runiv*Tinf))*sin(aoa)
      
      vel(1)=u1
      vel(2)=v1     

!     Find rotation and inverse rotation matrices.            
!      call ROTMAT(rmat,rmatinv,angle)
     
!     Rotate.            
!      fsvelrot=matmul(rmat,fsvel)
!      velrot=matmul(rmat,vel)      

!     Get matrix to convert to char vars.            
      call CHARM(mat,matinv,rhoinf,rhoinf*Runiv*Tinf)

!     Get deltas.            
      pvveco(1)=rho1-rhoinf
      pvveco(2)=signu*(vel(1)-fsvel(1))
      pvveco(3)=signu*(vel(2)-fsvel(2))
      pvveco(4)=p1-rhoinf*Runiv*Tinf

!     Convert and set char vars appropriately.            
      cvvec=matmul(mat,pvveco)      
!      cvvec(1)=0.0
!      cvvec(2)=0.0
!      cvvec(3)=0.0
      cvvec(4)=0.0

!     Convert back to primitives.            
      pvvecn=matmul(matinv,cvvec)
      
      rhoc=pvvecn(1)+rhoinf
      velc(1)=signu*pvvecn(2)+fsvel(1)
      velc(2)=signu*pvvecn(3)+fsvel(2)
      pc=pvvecn(4)+rhoinf*Runiv*Tinf

!     Rotate back velocity.            
!      velc=matmul(rmatinv,velc)
      
      uc=velc(1)
      vc=velc(2)

!     Get fluxes.            
      Ec=pc/(rhoc*(gam-1))+0.5*(uc**2+vc**2)
      Hc=Ec+pc/rhoc
      
      Fx(1)=rhoc*uc
      Fx(2)=rhoc*uc**2+pc
      Fx(3)=rhoc*uc*vc
      Fx(4)=rhoc*uc*Hc
      
      Gx(1)=rhoc*vc
      Gx(2)=rhoc*uc*vc
      Gx(3)=rhoc*vc**2+pc
      Gx(4)=rhoc*vc*Hc
      
      end subroutine OUTFLOWSUB
!     ***********************************************************      
      
!     ***********************************************************
      subroutine INFLOWSUP(Fx,Gx,Runiv,machinf,rhoinf,Tinf,gam,aoa)
      
      implicit none
      
      real,intent(in)::Runiv,machinf,rhoinf,Tinf,gam,aoa
      real,intent(out),dimension(4)::Fx,Gx
      
      real::uvel,vvel,absvel,Einf,Hinf,pinf
      
      pinf=rhoinf*Runiv*Tinf             
             
      absvel=machinf*sqrt(gam*Runiv*Tinf)
      
      uvel=absvel*cos(aoa)
      vvel=absvel*sin(aoa)
      
      Einf=pinf/(rhoinf*(gam-1))+0.5*absvel**2
      Hinf=Einf+pinf/rhoinf
      
      Fx(1)=rhoinf*uvel
      Fx(2)=rhoinf*uvel**2+pinf
      Fx(3)=rhoinf*uvel*vvel
      Fx(4)=rhoinf*uvel*Hinf
      
      Gx(1)=rhoinf*vvel
      Gx(2)=rhoinf*uvel*vvel
      Gx(3)=rhoinf*vvel**2+pinf
      Gx(4)=rhoinf*vvel*Hinf
      
      end subroutine INFLOWSUP
!     ***********************************************************

!     ***********************************************************
      subroutine OUTFLOWSUP(Fx,Gx,FF,GG)
      
      implicit none
      
      real,intent(in),dimension(4)::FF,GG
      real,intent(out),dimension(4)::Fx,Gx 
      
      Fx=1.0*FF!-0.5*FFF
      Gx=1.0*GG!-0.5*GGG
            
      end subroutine OUTFLOWSUP
!     ***********************************************************

!     **********************************************************      
      subroutine JAMESON(ncell,nedge,nvert,edge,grd,rho,u,v,E,p,        &
     &                   Runiv,machinf,rhoinf,Tinf,gam,aoa,vol,voln,    &
     &                volnm1,W1,Wn,Wnm1,fnp1,fn,fnm1,fvel,deltat,       &
     &                   k1,k2,dsf,iunflag,m,lam)      
                                    
      implicit none
      
      integer,intent(in)::ncell,nedge,nvert,iunflag            
      real,intent(in),dimension(nvert,2)::grd
      real,intent(in),dimension(ncell)::vol,voln,volnm1,deltat
      real,intent(in),dimension(nedge,2)::fvel
      real,intent(inout),dimension(ncell)::lam
      real,intent(inout),dimension(ncell)::rho,u,v,E,p        
      real,intent(in)::Runiv,machinf,rhoinf,Tinf,gam,aoa,k1,k2      
      integer,intent(in),dimension(4)::dsf
      integer,intent(in),dimension(nedge,4)::edge
      integer,intent(in),dimension(ncell)::m
      !integer,dimension(nedge,4)::edge2 
      
      real,dimension(ncell,4)::R,D,RT
      real,dimension(ncell,4)::W0
      real,intent(in),dimension(ncell,4)::Wn,Wnm1
      real,intent(inout),dimension(ncell,4)::W1
      real,intent(in)::fnp1,fn,fnm1
      
      integer::nc,nstage
            
      call CON(ncell,rho,u,v,E,W0)
       
      W1=W0
      !R=0.0
      D=0.0
                  
      do nstage=4,1,-1        
      
        R=0.0
                                                                   
        call FLUX(ncell,nedge,nvert,edge,grd,rho,u,v,E,p,fvel,R,        &
     &         Runiv,machinf,rhoinf,Tinf,gam,aoa)
                
        if(dsf(nstage).eq.1)then
        D=0.0        
         call DISS(ncell,nedge,nvert,W1,edge,rho,u,v,p,                 &
     &             gam,grd,k1,k2,D,m,lam,nstage,vol)
        endif
        
        if(iunflag.eq.1)then
          !print *, Wn-Wnm1 
          !Wn=W1
          !Wnm1=W1
          call FLUX2(ncell,fnp1,fn,fnm1,W1,Wn,Wnm1,                     &
     &    vol,voln,volnm1,RT)
            !print *, RT
          ! pause
        else
          RT=0.0
        endif
        
        do nc=1,ncell
          W1(nc,:)=W0(nc,:)-(deltat(nc)*                                &
     &      (R(nc,:)+D(nc,:)+RT(nc,:)))/(vol(nc)*nstage)
        enddo
         
        call PRIM(ncell,W1,gam,rho,u,v,E,p)
        
      enddo      
            
      end subroutine JAMESON
!     **********************************************************

!     ***********************************************************
      subroutine CON(ncell,rho,u,v,E,W)
      
      implicit none
      
      integer,intent(in)::ncell            
      real,intent(in),dimension(ncell)::rho,u,v,E            
      real,intent(out),dimension(ncell,4)::W
      
      integer::nc
      
      do nc=1,ncell
                
        W(nc,1)=rho(nc)
        W(nc,2)=rho(nc)*u(nc)
        W(nc,3)=rho(nc)*v(nc)
        W(nc,4)=rho(nc)*E(nc)
        
      enddo
      
      end subroutine CON
!     **********************************************************

!     **********************************************************
      subroutine PRIM(ncell,W,gam,rho,u,v,E,p)
      
      implicit none
      
      integer,intent(in)::ncell
      real,intent(in),dimension(ncell,4)::W
      real,intent(in)::gam
      
      real,intent(out),dimension(ncell)::rho,u,v,E,p
           
      integer::nc
      
      do nc=1,ncell
        
        rho(nc)=W(nc,1)        
        u(nc)=W(nc,2)/W(nc,1)
        v(nc)=W(nc,3)/W(nc,1)
        E(nc)=W(nc,4)/W(nc,1)
        
        p(nc)=(E(nc)-0.5*(u(nc)**2+v(nc)**2))*(gam-1)*rho(nc)
        
      enddo
      
      end subroutine PRIM
!     **********************************************************

!     **********************************************************
      subroutine INITIALISE(ncell,Runiv,machinf,rhoinf,Tinf,gam,        &
     &                      aoa,                                        &
     &                      rho,u,v,E,p)
     
      implicit none
      
      integer,intent(in)::ncell
      real,intent(in)::Runiv,machinf,rhoinf,Tinf,gam,aoa      
      real,intent(out),dimension(ncell)::rho,u,v,E,p
      
      real::SOS,absspd,pinf
      
      integer::nc
      
      SOS=sqrt(gam*Runiv*Tinf)
      absspd=machinf*SOS
      pinf=rhoinf*Runiv*Tinf
      
      do nc=1,ncell
      
        rho(nc)=rhoinf
        u(nc)=absspd*cos(aoa)
        v(nc)=absspd*sin(aoa)
        E(nc)=pinf/((gam-1)*rhoinf)+0.5*(u(nc)**2+v(nc)**2)
        p(nc)=pinf
      
      enddo
      
      end subroutine INITIALISE
!     **********************************************************     
      
!     **********************************************************
      subroutine READINT(ncell,nedge,nvert)
      
      implicit none
      
      integer,intent(out)::ncell,nedge,nvert
      
      open(100,file='griduns',status='unknown')
      
      read(100,*) ncell,nedge,nvert
      
      close(100)
      
      end subroutine READINT
!     **********************************************************         
      
!     **********************************************************   
      subroutine READSETTINGS(ncell,nedge,nvert,                        &
     &                  Runiv,machinf,rhoinf,Tinf,gam,aoa,tstop,CFL,grd &
     &                  ,edge,vol,k1,k2,dsf,m,iunflag,irest,            &
     &                  rfreq,nrealstop,npercyc,chord,rbf_type,sr)
      
      implicit none
      
      integer,intent(in)::ncell,nedge,nvert
      
      real,intent(out),dimension(nvert,2)::grd
      integer,intent(out),dimension(nedge,4)::edge
      real,intent(out),dimension(ncell)::vol
      real,intent(out)::Runiv,machinf,rhoinf,Tinf,gam,aoa,CFL,k1,k2,    &
     &  chord
      integer,intent(out)::tstop,iunflag,nrealstop,irest,npercyc,       &
     &  rbf_type
      real,intent(out)::rfreq,sr
      integer,intent(out),dimension(4)::dsf
      integer,intent(out),dimension(ncell)::m
      real::pi
            
      integer::nc,ed,vt,k,vert
      
      pi=4.0*atan(1.0)
      
      open(100,file='settings',status='unknown')
      
      read(100,*) machinf
      read(100,*) aoa
      aoa=aoa*pi/180.0
      read(100,*) rhoinf
      read(100,*) Tinf
      read(100,*) gam
      read(100,*) Runiv
      read(100,*) tstop
      read(100,*) CFL
      read(100,*) k1
      read(100,*) k2
      read(100,*) dsf(4),dsf(3),dsf(2),dsf(1)
      read(100,*) iunflag,irest
      if(iunflag.eq.1)then
        read(100,*) rfreq,nrealstop,npercyc,chord
        read(100,*) rbf_type,sr
      else
        !dt=1e6
        nrealstop=1
      endif
      
      close(100)
      
      open(101,file='griduns',status='unknown')
      
      read(101,*)
            
      do ed=1,nedge
        read(101,*) (edge(ed,k),k=1,4)
        !print *, (edge(ed,k),k=1,4)
      enddo
      !pause
      do vt=1,nvert
        read(101,*) vert,(grd(vert,k),k=1,2)
      enddo      
      !print *, "c"
      do nc=1,ncell
        read(101,*) vol(nc)
      enddo
      
      m=0
      do ed=1,nedge      
        if((edge(ed,3).gt.0).and.(edge(ed,4).gt.0)) then
        
          m(edge(ed,3))=m(edge(ed,3))+1
          m(edge(ed,4))=m(edge(ed,4))+1
         
        endif           
      enddo      
           
      close(101)
      
      end subroutine READSETTINGS                                              
!     **********************************************************

!     **********************************************************            
      subroutine TSTEP(ncell,nvert,nedge,edge,grd,CFL,gam,rho,u,v,p,    &
     &                 vol,deltat,lam,dt)      
      
      implicit none
                                                
      integer,intent(in)::ncell,nvert,nedge
      integer,intent(in),dimension(nedge,4)::edge
      real,intent(in)::CFL,gam,dt
      real,intent(in),dimension(ncell)::rho,u,v,p,vol
      real,intent(out),dimension(ncell)::deltat
      real,intent(out),dimension(ncell)::lam
      real,intent(in),dimension(nvert,2)::grd
            
      integer::nc,ed      
      real::u01x,u01y,dx,dy,S01,c01
      real,dimension(ncell)::deltatnew
      real::minedt,ul,ur,vl,vr,cl,cr
      
      deltat=12345.0
      
      lam=0.0
      
      do ed=1,nedge  

        if(edge(ed,3).gt.0)then
          cl=sqrt((gam*p(edge(ed,3))/rho(edge(ed,3))))
          ul=u(edge(ed,3))
          vl=v(edge(ed,3))
        else
          cl=sqrt((gam*p(edge(ed,4))/rho(edge(ed,4))))
          ul=u(edge(ed,4))
          vl=v(edge(ed,4))
        endif
        
        if(edge(ed,4).gt.0)then
          cr=sqrt((gam*p(edge(ed,4))/rho(edge(ed,4))))
          ur=u(edge(ed,4))
          vr=v(edge(ed,4))
        else
          cr=sqrt((gam*p(edge(ed,3))/rho(edge(ed,3))))
          ur=u(edge(ed,3))
          vr=v(edge(ed,3))
        endif       

        u01x=0.5*(ul+ur)
        u01y=0.5*(vl+vr)
        c01=0.5*(cl+cr)
        
        dx=grd(edge(ed,2),1)-grd(edge(ed,1),1)
        dy=grd(edge(ed,2),2)-grd(edge(ed,1),2)
        S01=sqrt(dx**2+dy**2)
        dx=dx/S01
        dy=dy/S01
            
        if((edge(ed,3).gt.0))then
          lam(edge(ed,3))=lam(edge(ed,3))+                              &
     &    (abs(dy*u01x-dx*u01y)+c01)*S01
        endif
        
        if((edge(ed,4).gt.0))then
          lam(edge(ed,4))=lam(edge(ed,4))+                              &
     &    (abs(-dy*u01x+dx*u01y)+c01)*S01
        endif
        !if(edge(ed,4).eq.50)print *, lam(edge(ed,4))        
        !if(edge(ed,3).eq.50)print *, lam(edge(ed,3))
            
      enddo
      
      !print *,lam(50)
      do nc=1,ncell
        if(lam(nc).gt.0.0)then
        deltat(nc)=CFL*vol(nc)/lam(nc)
        else
        print *, 'zero lam',nc,rho(nc),u(nc),v(nc)
        stop
        endif
      enddo
     
      !if(deltatmin.eq.0.0)then
      !  print *, "Zero timestep at edge",ed,length,cfl
      !  stop
      !endif
      
      !for any cell, set time step to minimum for that cell 
      !and any of its neighbours
      deltatnew=deltat
      do ed=1,nedge
        if((edge(ed,3).gt.0).and.(edge(ed,4).gt.0))then
          minedt=min(deltat(edge(ed,3)),deltat(edge(ed,4)))
          deltatnew(edge(ed,3))=min(deltatnew(edge(ed,3)),minedt)
          deltatnew(edge(ed,4))=min(deltatnew(edge(ed,4)),minedt)
        endif
      enddo      
      !copy back
      deltat=deltatnew                                                  
      
      ! safety check
      do nc=1,ncell
        !if(deltat(nc).eq.12345.0)then
        !  print *, "No timestep for cell",nc
        !  stop
        !endif
        !if(deltat(nc).lt.0.0)then
        !  print *, "Neg timestep for cell",nc
        !  stop
        !endif
        if(deltat(nc).gt.0.5*dt)then
          deltat(nc)=0.5*dt
          !print *, "tstep check",nc
        endif
      enddo

      !do nc=1,ncell
      !  print *, nc,deltat(nc)
      !enddo
      !pause
      
      end subroutine TSTEP       
!     **********************************************************

!     **********************************************************      
      subroutine RESID(ncell,nvert,nedge,grd,edge,rho,rhohold,deltat,   &
     &                 res,machinf,p,gam,rhoinf,Runiv,Tinf,t,tstop,     &
     &                 momx,momy,cx,cy,cm)

      implicit none

      integer,intent(in)::ncell,nvert,nedge,t,tstop
      integer,intent(in),dimension(nedge,4)::edge
      real,intent(in),dimension(ncell)::rho,rhohold,deltat,p
      real,intent(in)::machinf,gam,rhoinf,Runiv,Tinf,momx,momy
      real,intent(in),dimension(nvert,2)::grd
      real,intent(out)::res,cx,cy,cm

      real::cp,pinf,px,py,sgn,dx,dy,xc,yc
      
      integer::nc,ed,usedat

      res=0.0
      pinf=rhoinf*Runiv*Tinf
      
      do nc=1,ncell
        res=((rho(nc)-rhohold(nc))/deltat(nc))**2+res
      enddo       
      res=sqrt(res/float(ncell))
      
      !cx=0.0
      !cy=0.0
      
      !write(400,*) '***'
      !if(t.eq.tstop)then
      cy=0.0
      cx=0.0
      cm=0.0

      !print *, momx,momy
      
      do ed=1,nedge                                
        
        if(edge(ed,3).eq.-1)then
          usedat=edge(ed,4)
          sgn=-1.0
          cp=(2/(gam*machinf**2))*(p(usedat)/pinf-1)
          dx=(grd(edge(ed,2),1)-grd(edge(ed,1),1))
          dy=(grd(edge(ed,2),2)-grd(edge(ed,1),2))
          xc=grd(edge(ed,1),1)+0.5*dx-momx
          yc=grd(edge(ed,1),2)+0.5*dy-momy
          cy=cy+(-cp*dx)*sgn
          cx=cx+(cp*dy)*sgn
          cm=cm+(cp*dx*xc)*sgn+(cp*dy*yc)*sgn        
        endif
        if(edge(ed,4).eq.-1)then
          usedat=edge(ed,3)
          sgn=1.0
          cp=(2/(gam*machinf**2))*(p(usedat)/pinf-1)
          dx=(grd(edge(ed,2),1)-grd(edge(ed,1),1))
          dy=(grd(edge(ed,2),2)-grd(edge(ed,1),2))
          xc=grd(edge(ed,1),1)+0.5*dx-momx
          yc=grd(edge(ed,1),2)+0.5*dy-momy
          cy=cy+(-cp*dx)*sgn
          cx=cx+(cp*dy)*sgn
          cm=cm+(cp*dx*xc)*sgn+(cp*dy*yc)*sgn       
        endif
        
      enddo
      !endif      
      
      if(t.eq.tstop)then
      write(400,*) "ZONE"
      do ed=1,nedge                                
        
        if(edge(ed,3).eq.-1)then
          usedat=edge(ed,4)
          sgn=-1.0
          cp=(2/(gam*machinf**2))*(p(usedat)/pinf-1)
          px=0.5*(grd(edge(ed,1),1)+grd(edge(ed,2),1))
          py=0.5*(grd(edge(ed,1),2)+grd(edge(ed,2),2))        
          write(400,*) px,py,-cp
        endif
        if(edge(ed,4).eq.-1)then
          usedat=edge(ed,3)
          sgn=1.0
          cp=(2/(gam*machinf**2))*(p(usedat)/pinf-1)
          px=0.5*(grd(edge(ed,1),1)+grd(edge(ed,2),1))
          py=0.5*(grd(edge(ed,1),2)+grd(edge(ed,2),2))        
          write(400,*) px,py,-cp
        endif

        !cp=(2/(gam*machinf**2))*(p(usedat)/pinf-1)
        !px=0.5*(grd(edge(ed,1),1)+grd(edge(ed,2),1))
        !py=0.5*(grd(edge(ed,2),2)+grd(edge(ed,2),2))
        
        !write(400,*) px,py,cp

      enddo
      endif

      end subroutine RESID
!     *********************************************************
      
!     ***********************************************************
      subroutine CHARM(mat,matinv,rho,p)
      
      implicit none
      
      real,intent(out),dimension(4,4)::mat,matinv
      real,intent(in)::rho,p
      real::c,gam
      
      mat=0.0
      matinv=0.0
      gam=1.4
      
      c=sqrt(gam*p/rho)
      
      mat(1,1)=-c**2
      mat(1,4)=1
      mat(2,3)=rho*c
      mat(3,2)=rho*c
      mat(3,4)=1
      mat(4,2)=-rho*c
      mat(4,4)=1
      
      matinv(1,1)=-1/(c**2)
      matinv(1,3)=1/(2*c**2)
      matinv(1,4)=1/(2*c**2)
      matinv(2,3)=1/(2*rho*c)
      matinv(2,4)=-1/(2*rho*c)
      matinv(3,2)=1/(rho*C)
      matinv(4,3)=0.5
      matinv(4,4)=0.5
      
      end subroutine CHARM     
!     ***********************************************************
      
!     **********************************************************
      subroutine GCL(ncell,nedge,nvert,edge,grd,fvel,fnp1,fn,fnm1,      &
     &        voln,volnm1,volnp1)
      
      implicit none
      
      integer,intent(in)::ncell,nedge,nvert
      real,intent(in),dimension(ncell)::voln,volnm1
      real,intent(in),dimension(nedge,2)::fvel
      integer,intent(in),dimension(nedge,4)::edge
      real,intent(in),dimension(nvert,2)::grd
      real,intent(out),dimension(ncell)::volnp1
      real,intent(in)::fnp1,fn,fnm1
            
      integer::nc,ed
      real,dimension(ncell)::R
      real::flux,dx,dy
      
      R=0.0
      do ed=1,nedge
      
        dx=grd(edge(ed,2),1)-grd(edge(ed,1),1)
        dy=grd(edge(ed,2),2)-grd(edge(ed,1),2)
        !S01=(dx**2+dy**2)**0.5
        !dx=dx/S01
        !dy=dy/S01
        flux=(fvel(ed,1)*dy-fvel(ed,2)*dx) 
        !print *, flux
        if(edge(ed,3).gt.0)then
          R(edge(ed,3))=R(edge(ed,3))+flux
        endif
        if(edge(ed,4).gt.0)then
          R(edge(ed,4))=R(edge(ed,4))-flux
        endif
        !print *,flux
        
      enddo
      !pause

      do nc=1,ncell
        volnp1(nc)=(R(nc)-fn*voln(nc)-fnm1*volnm1(nc))/fnp1
      enddo
      
      end subroutine GCL
!     **********************************************************
             
!     **********************************************************
      subroutine EDGEVEL(nedge,nvert,edge,fnp1,fn,fnm1,grdnp1,grdn,     &
     &  grdnm1,fvel)
      
      implicit none
      
      integer,intent(in)::nedge,nvert
      real,intent(out),dimension(nedge,2)::fvel
      integer,intent(in),dimension(nedge,4)::edge
      real,intent(in),dimension(nvert,2)::grdnp1,grdn,grdnm1
      real,intent(in)::fnp1,fn,fnm1
            
      integer::ed
      real::vel1x,vel1y,vel2x,vel2y,velfx,velfy
      
      do ed=1,nedge 
      
        vel1x=fnp1*grdnp1(edge(ed,1),1)+                                &
     &    fn*grdn(edge(ed,1),1)+fnm1*grdnm1(edge(ed,1),1)
        vel1y=fnp1*grdnp1(edge(ed,1),2)+                                &
     &    fn*grdn(edge(ed,1),2)+fnm1*grdnm1(edge(ed,1),2)
     
        vel2x=fnp1*grdnp1(edge(ed,2),1)+                                &
     &    fn*grdn(edge(ed,2),1)+fnm1*grdnm1(edge(ed,2),1)
        vel2y=fnp1*grdnp1(edge(ed,2),2)+                                &
     &    fn*grdn(edge(ed,2),2)+fnm1*grdnm1(edge(ed,2),2)

        !dx=grd(edge(ed,2),1)-grd(edge(ed,1),1)
        !dy=grd(edge(ed,2),2)-grd(edge(ed,1),2)
        !S01=(dx**2+dy**2)**0.5
        !dx=dx/S01
        !dy=dy/S01
               
        velfx=0.5*(vel1x+vel2x)
        velfy=0.5*(vel1y+vel2y)

        fvel(ed,1)=velfx
        fvel(ed,2)=velfy

      enddo
      
      end subroutine EDGEVEL
!     **********************************************************
             
!     **********************************************************
      subroutine FLUX2(ncell,fnp1,fn,fnm1,Wnp1,Wn,Wnm1,                 &
     &  volnp1,voln,volnm1,RT)
      
      implicit none
      
      integer,intent(in)::ncell
      real,intent(in),dimension(ncell)::volnp1,voln,volnm1
      real,intent(in),dimension(ncell,4)::Wnp1,Wn,Wnm1
      real,intent(out),dimension(ncell,4)::RT
      real,intent(in)::fnp1,fn,fnm1
      
      integer::nc
      
      do nc=1,ncell
        RT(nc,:)=fnp1*(volnp1(nc)*Wnp1(nc,:))+fn*(voln(nc)*Wn(nc,:))    &
     &     +fnm1*(volnm1(nc)*Wnm1(nc,:))
      enddo
     
      end subroutine FLUX2
!     **********************************************************      

!     **********************************************************      
      subroutine CELLN(nedge,nvert,edge,nsurf,grd)
      
      implicit none
      
      integer,intent(in)::nedge,nvert
      integer,intent(in),dimension(nedge,4)::edge
      integer,intent(out)::nsurf
      real,intent(in),dimension(nvert,2)::grd
      
      integer::ed,points,index1,index2,hit1,hit2,p

      integer,dimension(10000)::list
      real::dist1,dist2
       
      nsurf=0
      points=0
      list=-999
      
      do ed=1,nedge
        if((edge(ed,3).eq.-1).or.(edge(ed,4).eq.-1)) then
          !nsurf=nsurf+1
          index1=edge(ed,1)
          index2=edge(ed,2)
          hit1=0
          hit2=0
          dist1=1.0e6
          dist2=1.0e6
          do p=1,points
            if(ed.gt.1)then
            dist1=((grd(index1,1)-grd(list(p),1))**2+                   &
     &  (grd(index1,2)-grd(list(p),2))**2)**0.5
            dist2=((grd(index2,1)-grd(list(p),1))**2+                   &
     &  (grd(index2,2)-grd(list(p),2))**2)**0.5
            endif
            if(points.gt.1)then
            if(index1.eq.list(p))then 
              hit1=1
            endif
            if(dist1.lt.1e-9)then
              hit1=1
              !print *, '*** Double addressed vertex found ***'
            endif
            if(index2.eq.list(p))then
              hit2=1
            endif
            if(dist2.lt.1e-9)then
              hit2=1
              !print *, '*** Double addressed vertex found ***'
            endif
            endif
          enddo
          if(hit1.eq.0)then 
            points=points+1
            list(points)=index1            
          endif
          if(hit2.eq.0)then 
            points=points+1
            list(points)=index2            
          endif
        endif         
      enddo

      nsurf=points
      
      end subroutine CELLN
!     **********************************************************           
      
!     **********************************************************      
      subroutine CELLFN(nedge,nvert,nsurf,grd,edge,cf)
      !(nedge,nvert,nsurf,grdnp1,edge,cf)
      
      implicit none
      
      integer,intent(in)::nedge,nvert
      integer,intent(in),dimension(nedge,4)::edge
      integer,intent(in)::nsurf
      integer,dimension(nsurf)::list
      real,intent(in),dimension(nvert,2)::grd
      real,intent(out),dimension(nsurf,2)::cf
      
      integer::ed,points,index1,index2,hit1,hit2,p
      real::dist1,dist2
       
      !nsurf=0
      points=0
      list=-999
      
      do ed=1,nedge
        if((edge(ed,3).eq.-1).or.(edge(ed,4).eq.-1)) then
          !nsurf=nsurf+1
          index1=edge(ed,1)
          index2=edge(ed,2)
          hit1=0
          hit2=0
          dist1=1.0e6
          dist2=1.0e6
          do p=1,points
            if(ed.gt.1)then
            dist1=((grd(index1,1)-grd(list(p),1))**2+                   &
     &  (grd(index1,2)-grd(list(p),2))**2)**0.5
            dist2=((grd(index2,1)-grd(list(p),1))**2+                   &
     &  (grd(index2,2)-grd(list(p),2))**2)**0.5
            endif
            if(points.gt.1)then
            if(index1.eq.list(p))then 
              hit1=1
            endif
            if(dist1.lt.1e-9)then
              hit1=1
            endif
            if(index2.eq.list(p))then
              hit2=1
            endif
            if(dist2.lt.1e-9)then
              hit2=1
            endif
            endif
          enddo
          if(hit1.eq.0)then
            points=points+1 
            list(points)=index1            
          endif
          if(hit2.eq.0)then
            points=points+1 
            list(points)=index2            
          endif
        endif         
      enddo
 
      open(600,file='surf_vertices.plt',status='unknown')
      do p=1,nsurf
        cf(p,:)=grd(list(p),:)
        write(600,*) cf(p,1),cf(p,2)
      enddo
      close(600)
      
      end subroutine CELLFN
!     **********************************************************            
            
!     ***********************************************************
      subroutine WRITEM(n_struct,dimnum,struct_grid,rbf_type,sr,M)

      implicit none
      
      integer,intent(in)::n_struct,dimnum,rbf_type
      real,intent(in),dimension(n_struct,dimnum)::                      &
     &  struct_grid
      double precision,intent(out),dimension(n_struct,n_struct)::M      
      real,intent(in)::sr
      !real,intent(in),dimension(dimnum)::k
      real::RBF,e_dist
      integer::i,j 

!     Write matrix M 

      do i=1,n_struct
        do j=1,n_struct
          e_dist=((struct_grid(j,1)-struct_grid(i,1))**2+                &
     &      (struct_grid(j,2)-struct_grid(i,2))**2)**0.5
          M(i,j)=RBF(rbf_type,sr,e_dist)            
        enddo
      enddo

      end subroutine WRITEM
!     ***********************************************************

!     ***********************************************************                                    
      subroutine MOVEGRID(nvert,nsurf,sr,rbf_type,inv_cfmat,grdh,       &
     &  cfh,grdnp1,cf)
      
      implicit none
      
      integer,intent(in)::nvert,nsurf,rbf_type
      double precision,intent(in),dimension(nsurf,nsurf)::inv_cfmat
      real,intent(in),dimension(nvert,2)::grdh                          
      real,intent(in),dimension(nsurf,2)::cf,cfh
      real,intent(out),dimension(nvert,2)::grdnp1
      real,intent(in)::sr
      
      integer::i,j,d
      real::RBF,summ,e_dist
      real,dimension(nsurf,2)::x
      !real,dimension(2)::k
      
      !delta=0.00000001
      !maxit=10000
      !k=1.0
      
      x=0.1
      
      do d=1,2
        x(:,d)=matmul(inv_cfmat,cf(:,d)-cfh(:,d))
      enddo
      
      !maxx=0.0
      do d=1,2
        do i=1,nvert
          summ=0.0
          do j=1,nsurf
            !e_dist=norm(2,grdh(i,:),cfh(j,:),k)
            e_dist=((grdh(i,1)-cfh(j,1))**2+                            &
     &            (grdh(i,2)-cfh(j,2))**2)**0.5
            summ=summ+RBF(rbf_type,sr,e_dist)*x(j,d)
          enddo
          grdnp1(i,d)=grdh(i,d)+summ
          !if(summ**2.gt.maxx**2)then
          !  maxx=summ
          !endif
        enddo
        !print *,maxx,"maxx"
        !pause
      enddo
       
      !open(999,file='test')
      !do i=1,nvert
      !  write(999,*) grdnp1(i,1),grdnp1(i,2)
      !enddo
      
      end subroutine MOVEGRID
!     **********************************************************      

!     **********************************************************      
      subroutine MOVESURF(nsurf,nreal,omega,xpit,ypit,amp,amp2,dt,cf,   &
     &  cfh,aoa,disp,momx,momy,momxh,momyh)

      implicit none

      integer,intent(in)::nsurf,nreal
      real,intent(in)::dt,omega,xpit,ypit,amp,amp2
      real,intent(out),dimension(nsurf,2)::cf
      real,intent(in),dimension(nsurf,2)::cfh
      real,intent(out)::aoa,disp
      real,intent(out)::momx,momy
      real,intent(in)::momxh,momyh
      
      real::time,pi
      integer::i

      pi=4.0*atan(1.0)

      !omega=2*0.0814*uinf/1.0
      time=float(nreal-1)*dt

      !cf(:,2)=cfh(:,2)+0.25*sin(2*pi*f*time)
      aoa=amp*sin(omega*time)
      print *, "AoA:",aoa
      aoa=aoa*pi/180.0

      do i=1,nsurf
        cf(i,1)=(cfh(i,1)-xpit)*cos(aoa)+(cfh(i,2)-ypit)*sin(aoa)+xpit
        cf(i,2)=-(cfh(i,1)-xpit)*sin(aoa)+(cfh(i,2)-ypit)*cos(aoa)+ypit
      enddo
      !print *,amp,amp2
      
      disp=amp2*sin(omega*time)
      cf(:,2)=cf(:,2)+disp
      momx=momxh
      momy=momyh+disp
        
      end subroutine MOVESURF
!     **********************************************************      

!     ***********************************************************
      SUBROUTINE MIGS (A,N,X,INDX)

!     Subroutine to invert matrix A(N,N) with the inverse stored
!     in X(N,N) in the output.  Copyright (c) Tao Pang 2001.

      IMPLICIT NONE
      INTEGER, INTENT (IN) :: N
      INTEGER :: I,J,K
      INTEGER, INTENT (OUT), DIMENSION (N) :: INDX
      double precision, INTENT (INOUT), DIMENSION (N,N):: A
      double precision, INTENT (OUT), DIMENSION (N,N):: X
      double precision, DIMENSION (N,N) :: B

      DO I = 1, N
        DO J = 1, N
          B(I,J) = 0.0
        END DO
      END DO
      DO I = 1, N
        B(I,I) = 1.0
      END DO

      CALL ELGS (A,N,INDX)

      DO I = 1, N-1
        DO J = I+1, N
          DO K = 1, N
            B(INDX(J),K) = B(INDX(J),K)-A(INDX(J),I)*B(INDX(I),K)
          END DO
        END DO
      END DO

      DO I = 1, N
        X(N,I) = B(INDX(N),I)/A(INDX(N),N)
        DO J = N-1, 1, -1
          X(J,I) = B(INDX(J),I)
            DO K = J+1, N
              X(J,I) = X(J,I)-A(INDX(J),K)*X(K,I)
            END DO
            X(J,I) =  X(J,I)/A(INDX(J),J)
          END DO
        END DO
      END SUBROUTINE MIGS
!     ***********************************************************

!     ***********************************************************
      SUBROUTINE ELGS (A,N,INDX)

!     Subroutine to perform the partial-pivoting Gaussian elimination.
!     A(N,N) is the original matrix in the input and transformed matrix
!     plus the pivoting element ratios below the diagonal in the output.
!     INDX(N) records the pivoting order.  Copyright (c) Tao Pang 2001.

      IMPLICIT NONE
      INTEGER, INTENT (IN) :: N
      INTEGER :: I,J,K,ITMP
      INTEGER, INTENT (OUT), DIMENSION (N) :: INDX
      double precision :: C1,PI,PI1,PJ
      double precision, INTENT (INOUT), DIMENSION (N,N) :: A
      double precision, DIMENSION (N) :: C

!     Initialize the index

      DO I = 1, N
        INDX(I) = I
      END DO

!     Find the rescaling factors, one from each row

      DO I = 1, N
        C1= 0.0
        DO J = 1, N
          C1 = DMAX1(C1,ABS(A(I,J)))
        END DO
        C(I) = C1
      END DO

!     Search the pivoting (largest) element from each column

      DO J = 1, N-1
        PI1 = 0.0
        DO I = J, N
          PI = ABS(A(INDX(I),J))/C(INDX(I))
          IF (PI.GT.PI1) THEN
            PI1 = PI
            K   = I
          ENDIF
        END DO

!     Interchange the rows via INDX(N) to record pivoting order

      ITMP    = INDX(J)
      INDX(J) = INDX(K)
      INDX(K) = ITMP
      DO I = J+1, N
      PJ  = A(INDX(I),J)/A(INDX(J),J)

!     Record pivoting ratios below the diagonal

      A(INDX(I),J) = PJ

!     Modify other elements accordingly

      DO K = J+1, N
        A(INDX(I),K) = A(INDX(I),K)-PJ*A(INDX(J),K)
      END DO
      END DO
      END DO

      END SUBROUTINE ELGS
!     ***********************************************************


      

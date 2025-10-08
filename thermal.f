c August 31/2025
c***********************************************************************
c THERMAL sets up the principal A(1,1) and B(1) matrices for solving
c the linear heat balance equation set. It also solves for the new
c temp T and melt BMELT by application of the tridiagonal matrix
c algorithm.  For nodes within the meltzone, T is determined from
c BMELT within subroutine FTEMP.
c Modified to accept a Neumann bottom boundary condition. 10/12/2023
c 2023 changes involve the heat balance in node 1 and which tridiag to
c access Neumann. Search on Neumann (as well as RJ)
c***********************************************************************

      subroutine thermal(dtmin,dzmin,recross,*,errtallowd,Neumann,
     &jday,ihour,prnt_hour,icalcstep,iunit,thinnode,TsurfEst) 
     &  !RJ 10/12/23 & RJ12/17!REJ_2025/08/11 
c      

c 
c Called from MAIN
c
c Calls TRIDIAG,TRIDIAGN, FTEMP function FLIQUID
      implicit none
      include 'const'
      include 'arrays'
c Arguments
c
c dtmin : minimum allowed time step [s]
c dzmin : minimum nodal thickness (except for precipitation cases) [m]
c recross : Used in trapping overshot phase boundary.
c errtallowd: maximum allowed linearization error per time-step in heat balance
c             eq, expressed in units of temperature [K]
c
      logical recross,prnt_hour  !RJ12/17
      integer jday,ihour,icalcstep !REJ_2025/08/11;08/28
      double precision dtmin,errtallowd,dzmin,TsurfEst
c
c Local
c
c blmax: Bulk density of liquid water corresponding to upper melt
c         zone limit [kg/m^3]
c blmin: Bulk density of liquid water corresponding to lower melt
c         zone limit [kg/m^3]
c i: looping index.
c j: looping index.
c m: m=k(i)
c ok: Used in trapping overshot phase boundary.
c redo: Logical flag denoting that it is ok to redo the iteration
c       = dt .gt.1.1*dtmin .and. istart.eq.0.and.iskip(n).eq.0
c
c Passed through common (incomplete) 
c
c a(nd,3) : Coefficient array in linear equation matrix
c b(nd) : Constant array in linear equation matrix
c bbo(nd) : .5*(Old nodal conducted and convected heat plus absorbed solar)
c dsol(nd) : Nodal absorbed solar radiation [W/m^2]
c gk(nd) : Constant for calculating temperature from melt [K]
c gv(nd) : Coefficient for calculating temperature from melt
c iskip(nd) : 1 signifies that convergence criteria are skipped for
c              this CV
c n: Number of CVs in model
c nsoil: Number of soil CVs
c melt(nd): Signifies if CV in meltzone if = 1
c qf(nd) : .5*specific heat of water*CV mass water flux 
c qk(nd) : .5*(thermal conductivity at upper nodal boundary/nodal thickness), 
c qs(nd) : Coefficient on CV stored heat [W/m^2 K]
c To(nd) : Past T(nd)  ! 11/10/23  RJ
c unbar(nd) : Average nodal convective mass flux of water  [K/m^2 s]
c convect : = -ci(n)*us(n)*dt Fix later for rainfall
c prnt_hour : if .true., code prints tridiag/tridiagN input & other. RJ12/17
      integer i,m,j,iunit
      double precision blmin,blmax,precip_coef !10/12/2023
     &    !+REJ_2025/08/07
      logical ok,redo,Neumann,thinnode !10/12/2023 added option for a bottom Neumann boundary
c function
c
c fliquid: liquid water as function of temperature and liquid water content.
c
      double precision fliquid
c
11    format(10g13.5)
c
c interior CVs  RJ note: Neumann BC only impacts bottom CV
      do 10 i = 2,n-1
         a(i,1)=-qk(i)
         if(i .gt. nsoil)then
           a(i,2)=qs(i)+qf(i)+qk(i+1)+qk(i)
           a(i,3)=-(qf(i+1)+qk(i+1))
         else
c Temp  Temporary creation of sink at snow/soil interface
           a(i,2)=qs(i)+qk(i+1)+qk(i)
           a(i,3)=-qk(i+1)
         endif
c Temp
         b(i)=qs(i)*to(i)+unbar(i)*ho(i)+.5*dsol(i)+bbo(i)-a(i,1)*
     &   gk(i-1)-a(i,2)*gk(i)-a(i,3)*gk(i+1)
         a(i,1)=a(i,1)*gv(i-1)
c        if(melt(i).gt.0)
         if(melt(i).gt.0 .and. bw(i) .ge. 1d0)
     &                a(i,2)=dlm*dz(i)/dt+a(i,2)*gv(i)
         a(i,3)=a(i,3)*gv(i+1)
 10   continue
c
c top node
c         if(icalcstep .ge. 1192) stop 'THERMAL step 1192'
c        if(icalcstep .gt. 500)write(iunit,*)'THERMAL:mass,limit,n',
c     &   n,dmass(n),0.8d2*dzmin,iskip(n)
c     Compare with earlier version sent to ADI.
c      if(dmass(n) .gt. 0.8d2*dzmin)then  !Next redundant.  Add in rainfall
CR      if(icalcstep .gt. 700)write(*,*)'TH: thinnode check',thinnode(n),
CR     &    icalcstep
      IF(.not. thinnode)THEN
c         if (icalcstep .ge. 550)write(120,*)'TH:COMPUTE BRANCH qs(n)',
c     &    icalcstep, qs(n),'cx',qs(n)*dt
         if(us(n) .lt. 0d0)then  !REJ_2025/08/07 snowfall.  
          convect  = - ci(n)*us(n)*dt  !Add rainfall and wet ssnow
c                     Added *dt  !REJ_2025/08/28
         else
          convect  =  0d0
         endif
         a(n,1) =-qk(n)
         a(n,2) = qs(n) + qf(n) + qk(n) - topfluxv !REJ2025/08/07 convect omitted .5 omitted
         a(n,3) = 0d0
         b(n) = qs(n)*to(n)+unbar(n)*ho(n)+.5d0*dsol(n) + topfluxk !REJ2025/08/07 08/11 ditto
     &     + bbo(n) - a(n,1)*gk(n-1)-a(n,2)*gk(n) !REJ_2025/08/07 
         a(n,1)=a(n,1)*gv(n-1)
c        b(n)=qs(n)*to(n)+convect*273.15d0 !REJ2025/08/07 08/11 omit convect term
c 1      gk(n-1)-a(n,2)*gk(n)+(topfluxk-ci(n)*us(n)*(tprecip-to(n)))/2d0 REJ_2025/08/07

c        if (icalcstep .ge. 550)write(iunit,*)'Thermal b(n)',b(n)

c        if(melt(n) .gt. 0)then
         if(melt(i).gt.0 .and. bw(i) .ge. 1d0)then
            a(n,2)=dlm*dz(n)/dt+a(n,2)*gv(n)
         end if
      ELSE
c     circumvent thermal balance for top elements with minimal mass
         iskip(n)=1 !REJ:2025/08/30 Moved up
         melt(n)=0
         a(n,3)=0d0
         a(n,2)=1d0
         a(n,1)=0d0
         if(prcp .le. 0d0)then !no precip
           b(n)=To(n)
         else  !snowfall
c          b(n)=tprecip
           b(n)=(tprecip+to(n-1))/2d0
c           if(ltype(k(n)) .le. 1)b(n)=dmin1(th(k(n)),b(n))!REJ copied
c                                 from sntherm89 recheck before using.
         endif
c        iskip(n)=1 !REJ:2025/08/30 Moved up                         
      ENDIF
c
c     Bottom node,qk(2),T(2)-T(1),botflux
      If(Neumann)then
c Bottom specified constant flux boundary. Added by RJ 10/12-10/30/23
         a(1,1)= 0d0
         if(1 .ge. nsoil)then 
c          This loop has potential water flow, but is not used 
           a(1,3)=-(qf(2)+qk(2))
           a(1,2)=qs(1)+qf(1)+qk(2)+qk(1)
           b(1)=qs(1)*to(1)+unbar(1)*ho(1)+.5*dsol(1)+bbo(1)
     &      -a(1,2)*gk(1)-a(1,3)*gk(2)
         else
c          This loop does not have water flow. Melt left in,
c          But still needs to be tested.      
           a(1,3)=-qk(2)
           a(1,2)=qs(1)+qk(2)
           b(1)=qs(1)*to(1)+.5*dsol(1)+bbo(1)
     &      -a(1,2)*gk(1)-a(1,3)*gk(2)
         endif
          a(1,3)=a(1,3)*gv(2)
          if(melt(1).gt.0 .and. bw(1) .ge. 1d0)
     &    a(1,2)=dlm*dz(1)/dt+a(1,2)*gv(1) ! I need to check this.  RJ
      Else
c Bottom specified constant temperature boundary 
         a(1,3)=0d0
         a(1,2)=1d0
         a(1,1)=0d0
         b(1)=to(1)
      Endif
c
c solve the tridiagonal matrix       
c
cprint------------------------------------------------------------------
c      All print blocks are new.  I cleaned this up alot on 12/17/23 and 
c      suggest replacing the whole print block.
c       write(*,*)jday, ihour,prnt_hour
       IF(prnt_hour)THEN  
       write(*,*)jday,ihour
c        stop  !restore   FIX
        open(69,file='thermal_out',status='unknown')
        write(69,*) 'thermal_out','jday=',jday,'ihour=',ihour
        write(69,*) 'This print-out picked a date-time where the 4 bottom
     &  snow nodes are in'
        write(69,*)'the temperature regime and the 4 top nodes are 
     &  melting.'  
        write(69,*)'Otherwise, the routine will need some adjustments'
        write(69,*)'Correct coeffs for the top node not printed yet'
        write(69,*)  
        write(69,*) 'INPUT MATRICES TO TRIDIAG/TRIDIAGN IN T() REGIME'
        write(69,*) 'a(i,3)=-qf(i+1)-qk(i+1))'
        write(69,*) 'a(i,2)=qs(i)+qf(i)+qk(i+1)+qk(i)'
        write(69,*) 'a(i,1)=-qk(i)'
        write(69,*) 'b(i)=qs(i)*to(i)+unbar(i)*ho(i)+.5*dsol(i)+bbo(i)'
        write(69,*)   '-a(i,1)*gk(i-1)-a(i,2)*gk(i)-a(i,3)*gk(i+1)'
        write(69,*) 
16      format(a5,4i15)
13      format(a7,4f16.6)
14      format(a7,4i16)
        write(69,16)'node=',nsoil+1,nsoil+2,nsoil+3,nsoil+4
        write(69,14) 'melt',(melt(i),i=nsoil+1,nsoil+4)
        write(69,13)'to',     (to(i),i=nsoil+1,nsoil+4)
        write(69,13)' b',      (b(i),i=nsoil+1,nsoil+4)
        write(69,13)'a3',    (a(i,3),i=nsoil+1,nsoil+4)
        write(69,13)'a2',    (a(i,2),i=nsoil+1,nsoil+4)
        write(69,13)'a1',    (a(i,1),i=nsoil+1,nsoil+4)
        write(69,*)
        write(69,*)'INPUT MATRICES TO TRIDIAG/TRIDIAGN IN Pmelt() REGME'
        write(69,*) ' a(i,3)=a(i,3)*gv(i+1)'
        write(69,*) ' a(i,2)=dlm*dz(i)/dt+a(i,2)*gv(i)' 
        write(69,*) ' a(i,1)=a(i,1)*gv(i-1)'
        write(69,16)'node=', n-3,n-2,n-3,n
        write(69,14)'melt', (melt(i),i=n-3,n)
        write(69,13)'to',     (to(i),i=n-3,n)
        write(69,13)' b',      (b(i),i=n-3,n)
        write(69,13)'a3',    (a(i,3),i=n-3,n)
        write(69,13)'a2',    (a(i,2),i=n-3,n)
        write(69,13)'a1',    (a(i,1),i=n-3,n)
        write(69,*)
        write(69,*)'OTHER PARAMETERS: Definitions'
        write(69,*) 'qk(nd)=.5*(thermal conductivity at upper nodal'
        write(69,*)           ' boundary/dz)' 
        write(69,*) 'qf(nd)=.5*specific heat of water*'
        write(69,*)            'nodal mass water flux'
        write(69,*) 'T( ) REGIME'
        write(69,*) 'qs(i)=(bt(i)*ct(i)+f(i)*dlm+dbvdt(i)*dls)'
        write(69,*) '      *dz/dt'
        write(69,16)'node=',1,2,3,4
        write(69,13) 'gv',      (gv(i),i=1,4)
        write(69,13) 'gk',      (gk(i),i=1,4)
        write(69,13) 'qs',      (qs(i),i=1,4)
        write(69,13) 'qk',      (qk(i),i=1,4)
        write(69,13) 'qf',      (qf(i),i=1,4)
        write(69,13) 'bbo',     (bbo(i),i=1,4)
        write(69,*)  'Pmelt() REGIME'
        write(69,*) 'qs(i)=(bt(i)*ct(i)+dbvdt(i)*dls)*dz(i)/dt'
        write(69,16)'node=', n-3,n-2,n-3,n
        write(69,13) 'gv',      (gv(i),i=n-3,n)
        write(69,13) 'gk',      (gk(i),i=n-3,n)
        write(69,13) 'qs',      (qs(i),i=n-3,n)
        write(69,13) 'qk',      (qk(i),i=n-3,n)
        write(69,13) 'qf',      (qf(i),i=n-3,n)
        write(69,13) 'bbo',     (bbo(i),i=n-3,n)
        write(69,*)
        close(69)
c       endif   !RJ12/17
      ENDIF 
cprintend---------------------------------------------------------------
      if(.not. Neumann)then
        call tridiag(n,a,t,b,Neumann)
      else
        call tridiagN(n,a,t,b,Neumann)
      endif 
c-----------------------------------------------------------------------
c      This block prints to converge_out, iunit = 120
          if(icalcstep .gt. 550)then
            write(iunit,*)'2. THERMAL'
            write(iunit,6)n,icalcstep,dt,dz(i),qs(i)*dt,bt(i),ct(i),
     &            f(i), dbvdt(i)*dls,T(i),To(i)
6           format(2i5,f10.3,f10.6,5f12.6,2f7.2)
            write(iunit,*)'Surface balance in TridMatrix:',
     &      a(n,2)*T(n) + a(n,1)*T(n-1)-b(n),'T(n)',T(n)
          endif
c-----------------------------------------------------------------------
c
c Next first checks to see if melt-zone has been skipped.  If so,
c reset variables and redo this time-step.
c For elements which were within the meltzone, call FTEMP to solve for
c new temp as a function of bmelt.  FTEMP also checks to see if 
c the node has over-shot the meltzone. If the phase boundary was 
c over-shot by more than 5%, the time step is shortened and the itera- 
c tion is repeated.  See routine FTEMP for details.
c     redo=dt .gt. 1.1*dtmin .and. istart .eq.0 .and. iskip(n) .eq. 0
      redo=dt .gt. 1.1*dtmin .and. istart .eq.0
      ok=.true.
      do 40 j=1,n  ! 11/7/23  RJ
c      do 40 j=1,n-1  !11/7/23  RJ
      i=n+1-j  !11/7/23  Moved up RJ
      m=k(i)   !11/7/23  Moved up RJ
      if(bw(i) .le. 1d0)goto 40 !11.10.03
      if(melt(i) .le. 0 .or. bw(i) .lt. 1d0)then  !RJ
c Check to see if melt-zone was totally skipped 
        if(to(i) .lt. tl(m).and. t(i).gt.th(m) .or.
     1   to(i) .gt. th(m) .and. t(i) .lt. tl(m))then
           ok=.false.
           dt=dmin1(0.1d0,dtmin)
           recross=.true.
           call reset
           return 1
         endif
       else
         blmin=fliquid(bw(i),tdl(m,1),bdjp(m),a243(m),tdl13(m,1)
     &          ,flglim(m,1))
         blmax=fliquid(bw(i),tdl(m,2),bdjp(m),a243(m),tdl13(m,2)
     &          ,flglim(m,2))
c Routine FTEMP estimates liquid water and temperature from melt  RJ 12/15/23
         call ftemp(t(i),dt,bl(i),bw(i),a1(m),
     1      a243(m),unbar(i),dz(i),dzo(i),bdjp(m),td13(i),tl(m),
     2      th(m),blo(i),bmelt(i),
     3      flgo(i),blmin,blmax,dtmin,redo,ok,recross,f(i),
     4      errtallowd,*50) !Apologies for this goto! Will fix later. RJ
      endif
cREJ      if(i .eq. n .and. iskip(n) .eq.1)write(*,*)'THERM:iskip,To,T',
cREJ     & iskip(n),To(n),T(n)
      goto 40
c  Next is situation where iteration is being repeated (which was
c  determined in FTEMP)
50    return 1
40    continue
      return
8888      end

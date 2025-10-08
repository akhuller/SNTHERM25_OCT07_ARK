c***********************************************************************
c THPARAM computes thermal parameters
c***********************************************************************
      subroutine thparam(Neumann,jday,ihour,prnthour,iunit,icalcstep) !REJ_2025/08/10

      implicit none
      include 'const'
      include 'arrays'
c
c Called from MAIN
c
c Calls function dblcurve 
c
c
c Local
c
c i: Nodal looping index
c m: Layer index = k(i)
c cl12: cl12=cl/2
c dum0: dum0=a223*td13*td13.
c dum1: dum1=dum0+1d0.
c dum2: dum2=1.4142136d0*a213*td13.
c dzm:  dzm=dz(i-1)
c thkm: thkm=thk(i-1)
c
      integer i,m,loopstart,jday,ihour,iunit,icalcstep !REJ_2025/08/10 
      double precision dum1,dum2,dum0
      double precision dzm,thkm,cl12
      logical Neumann,prnt,prnthour  !  10/26/22  RJ RJ12/17 line28

c Passed through common
c
c to(nd) : Old nodal temperature [K]
c ci(nd) : Specific heat of ice        [J/kg-K]
c ct(nd) : Nodal combined specific heat of all constituents [J/kg-K]
c bw(nd) : Nodal water constituents bulk density [kg/m^3]
c td(nd) : Nodal depression temperature, t-273.15.  [K]
c td13(nd) : Nodal temperature depression from 0 C raised to (1/3) power
c flo(nd) : Nodal fraction of liquid (unfrozen) water due to both capillary
c          and adsorbed potential, blo/bwo
c gk(nd) : Constant for calculating temperature from melt [K]
c gv(nd) : Coefficient for calculating temperature from melt
c ho(nd) : Water flux enthapy adjustment
c f(nd) : Slope of freezing curve for node [kg/m^3 K] (Note:f is
c   defined here as the change in bl over dT, and in the report as
c   the change in the fraction of liquid water over dT).
c bdjp(ld) : 0.75*Dry soil bulk density*plasticity index (0.75*bd*djp)
c a1(ld) : Constant in unfrozen water formula
c a2(ld) : Constant in unfrozen water formula
c a243(ld) : a2^^(4/3)
c ltype(ld) : Layer type 1=snow >=90 user supplied soil layer
c dz(nd) :  Elemental thickness [m]
c thk(nd) : Current Nodal thermal conductivity [W/m-K]
c bt(nd) : Nodal total bulk density [kg/m^3]
c unbar(nd) : Average nodal convective mass flux of water  [K/m^2 s]
c dt :  Current time step [s]
c melt(nd) : Signifies if node in meltzone if = 1
c cl : Specific heat of water at 273.15 K (4217.7) [J/kg-K]
c dlm : Latent heat of fusion for ice(3.335E5) [J/kg]
c qk(nd) :.5*(thermal conductivity at upper nodal boundary/dz) 
c qs(nd) : Coefficient on nodal stored heat
c bdcd(ld) : Dry soil bulk density*specific heat (bd*cd)
c flgo(nd) : Old fractional unfrozen water due to capillary potential
c qf(nd) : .5*specific heat of water*nodal mass water flux 
c u(nd) : Nodal convective mass flux of water [kg/m^2-s]
c a213(ld) : a2^^(1/3)
c a223(ld) : a2^^(2/3)
c dbvdt(nd) :1)Change in saturation vapor density per degree K [kg/m^3 K]
c           :2)Change in bulk vapor density per degree K [kg/m^3 K]
c dls : Latent heat of sublimation of ice (2.838E6) [J/kg] at 273.15
c hs(nd) : Water flux enthapy adjustment.
c uvapor(nd) : Net diffusive mass flux of water vapor for element [K/m^2 s]
c n : Number of nodes in model
c prnt : optional print-out for testing if .true.  User sets this RJ12/17 line 73
c prnthour : print time for prnt = .true.  code sets prnthour RJ12/17 line 74
c jday :  print day    ! RJ12/17 line 74
c ihour :  print hour  ! RJ12/17 line 75
c
c Functions
c
      double precision dblcurve
c
      data prnt/.false./ !RJ12/17 Line 82 Make .true. to option print

c-----------------Optional printout-------------------------------------
c RJ 10/03/23 Optional print block that activates for prnt = .true.
c Also set a print-out date/time. 'prnthour' turns positive when time is
c reached.  Program stops after print-out completes.
      if(prnt .and.jday .eq. 7 .and. ihour .eq. 16)
     &   prnthour = .true.  
      if(prnthour)open(69, file= 'thparam_out', status = 'unknown')
      if(prnthour)write(69,*)'thparam_out'  
      if(prnthour)write(69,*)'jday =',jday,'  ihour =',ihour,'dt = ',dt
      if(prnthour)write(69,*)
     &'  node      qs         bl          ct          f           dz'!RJ12/17   
cRJ12/16 bw to bl in above
c-----------------EndOptional print-------------------------------------


      cl12 = 0.5*cl


c RJ 10/26/23. For a Neumann lower BC, the Jacobian matrix requires
c storage term parameters for node 1. Influx parameters are already  
c computed. Outflux is constant = botflux  
           
      if(Neumann)      loopstart = 1   !RJ 10/26/23
      if(.not. Neumann)loopstart = 2   !RJ 10/26/23
c      do 100 i=2,n                    !RJ 10/26/23
       DO 100 i= loopstart,n           !RJ 10/26/23
         m= k(i)  ! m is the layer code. 1 
         if(i .gt.1)then               !RJ 10/26/23
           dzm = dz(i-1) 
           thkm= thk(i-1)
         endif                         !RJ 10/26/23
         IF(to(i) .le. 273.15d0)THEN
c Node is below freezing
c REJ 1/17/24             if(i.eq.1 .and. unbar(1).gt.1d-8 .and. Neumann)STOP 
c REJ 1/17/24      &     'Neumann with water flow in bottom node not treated yet'
            ci(i)=-13.3d0+7.8d0*to(i)
            if(dabs(unbar(i)) .gt. 0d0)then ! Water flow occurring
c REJ 1/17/24     stop 'wrong water flow loop in thparam'!RJ rechking this 
               hs(i)=(cl-ci(i))*273.15d0+ci(i)*to(i)
              if(ltype(m) .gt. 1)then ! soil
                  dum0=a223(m)*td13(i)*td13(i)
                  dum1=dum0+1d0
                  dum2=1.4142136d0*a213(m)*td13(i)
                  hs(i)=hs(i)-((cl-ci(i))*(1d0/bw(i)))*
     &                 ((bw(i)-bdjp(m))*datan(a1(m)*td(i))/a1(m)+
     1                 .75d0*bdjp(m)*(dlog((dum1-dum2)/(dum1+dum2))
     2                 +2d0*datan(dum2/(1d0-dum0)))/(1.414213d0*a2(m)))
              endif
               ho(i)=hs(i)-(1.-flgo(i))*dlm
            else ! No water flow
c           It isn't necessary to calculate ho when there is no flow
               ho(i)=0d0
               hs(i)=0d0
            end if
            f(i)=dblcurve(td(i),bw(i),bdjp(m),a1(m),a243(m),td13(i),
     &           flgo(i),melt(i))
            ct(i)=(bdcd(m)+bw(i)*(cl*flo(i)+ci(i)*(1.-flo(i))))/bt(i)
c           if(melt(i) .eq. 0 )then
            if(melt(i) .eq. 0 .or. bw(i) .lt. 1d0)then !temp zone RJ12/17

cRJ12/16  above comment
               qs(i)=(bt(i)*ct(i)+f(i)*dlm+dbvdt(i)*dls)*dz(i)/dt

c-----------------------------------------------------------------------
            goto 9999
c Prints to converge_out. Compare these surface arrays with those in 
c converge_out (iunit = 120).  They should be equal !REJ_2025/08/28
            if(icalcstep.gt.550 .and. i .eq.n)then
               write(iunit,*) '1. THPARAM'        
               write(iunit,6)
     &         '     n   step    dt           dz      ',
     &         'qs*dt          bt        ct f        dbvdt*Ls'
               write(iunit,6)n,icalcstep,dt,dz(i),qs(i)*dt,bt(i),ct(i),
     &            f(i), dbvdt(i)*dls
c6              format(2i5,7f12.6)
             endif
9999  continue
c-----------------------------------------------------------------------
c REJ/08/28/2025 Next print superceded by above block
c               if(prnthour)write(69,5)i,qs(i),bl(i),ct(i),f(i),dz(i),
c    &              To(i)
c5              format(i5,7f12.6)
               gk(i)=0.0
               gv(i)=1.0
            else                                       !melt zone RJ12/17
cRJ12/16  above comment                                        
               gv(i)=1/f(i)
               gk(i)=to(i)
               qs(i)=(bt(i)*ct(i)+dbvdt(i)*dls)*dz(i)/dt
c REJ/08/28/2025 Next print superceded by above block
c               if(prnthour)write(69,5)i,qs(i),bl(i),ct(i),f(i),dz(i),
c     &              To(i)
            endif
         ELSE
c Node is above freezing
            ho(i)=cl*to(i)
            hs(i)=ho(i)
            f(i)=0.0
            ct(i)=(bdcd(m)+cl*bw(i))/bt(i)
            qs(i)=bt(i)*ct(i)*dz(i)/dt
            gk(i)=0.0
            gv(i)=1.0
         ENDIF
c-----------------------------------------------------------------------
c Prints to converge_out. Compare these surface arrays with those in 
c converge_out (iunit = 120).  They should be equal !REJ_2025/08/28
            if(icalcstep.gt.550 .and. i .eq.n)then
              write(iunit,*) '1. THPARAM'        
              write(iunit,*)
     &        '  n   step    dt           dz       qs*dt       bt',
     &        '      ct            f      dbvdt*Ls   Test      To '
               write(iunit,6)n,icalcstep,dt,dz(i),qs(i)*dt,bt(i),ct(i),
     &            f(i),dbvdt(i)*dls,T(i),To(i)
6              format(2i5,f10.3,f10.6,5f12.6,2f7.2)
             endif
c-----------------------------------------------------------------------
c        if(i.eq.n)qs(i)=qs(i)+ct(i)*uvapor(i)/2d0
         if(i.gt.1)qk(i)=thk(i)*thkm/(thk(i)*dzm+thkm*dz(i))
         qf(i)=-cl12*u(i)
 100  CONTINUE
      if(prnthour)then
         write(69,*)
         close (69)
      endif

      return
      end


c**********************************************************************  
c  NEWSNOW adds new nodes. New precip is added in elements of thickness 
c  dzinc (default 4cm for snow/1cm for rain). Precip is combined with
c  existing mass within an element until thickness dzinc is reached. 
c  Excess precip is placed in a new top element.  After precip stops, 
c  if top element exceeds dzn or next to top element exceeds dznm, 
c  they are subdivided within SUBDIVIDE. Maximum timestep 'snowfalldtmax' 
c  is set to preclude more than one new node being added per iteration.
c  REJ_2025.  MASS FLUXES ARE POSITIVE UPWARDS.  To be changed later.
c***********************************************************************
      subroutine newsnow(rate,rateo,dzinc,bwfall,ssminsnow,dsnowfall,
c     &   repeat,ssisnow,blfall,fullnode,a1snow) REJ_2025/05/28
     &   repeat,ssisnow,blfall,fullnode,a1snow,
     &  thinnode)!REJ_2025/08/31 
     

c IMPORT ANT. THIS ROUTINE WAS REVISED IN NOV, 1995 TO ACCEPT
c SWE RATHER THAN SNOW ACCUMULATION AS THE PRECIP RATE.
c
c Called from MAIN.
c
      implicit none
      include 'const'
      include 'arrays'
      include 'l_arrays.txt'
c
c Arguments
c rate, rateo: Is snowrate(o) or rainrate(o) when called from MAIN 
c                                                       (m-SWE/s)
c dzinc: Thickness of element of new precip addition (default 4cm/1cm)
c bifall: Bulk density of newly fallen snow (dry + wet)
c         (0d2 if rainfall occurs)
c ssminsnow : -ssisnow/(1d0-ssisnow)
c ssisnow : irreducible water saturation for snow [.04]
c dsnowfall:  Effective grain diameter of falling snow particle (m)
c repeat : logical flag denoting a repeat of this iteration (if true)
c blfall: Bulk density of liquid precipitation
c fullnode: Flag used to indicate whether top node is full (dz exceeds dzinc)
c a1snow : unfrozen water constant
c thinnode : .true. if dz .lt. dzmin .or. dmass(i) .lt. 80d0*dzmin
c      dzmin set to 2d-3 in MAIN data
c
      double precision rate,rateo,dzinc,bifall,ssminsnow,dsnowfall
cNov13,1995      double precision ssisnow,blfall,a1snow,bwfall,binew(nd)
      double precision ssisnow,blfall,a1snow,bwfall
      logical fullnode,repeat,thinnode !scalar thinnode in newsnow.f

c Local
c 
c avrate: average of old and current precipitation rate (m-swe/s) 
c April 4. Now using current value for avrate.
c         (avrate=(rate+rateo)/2d0).
c dum: dum=a1snow*(273.15d0-to(n))
c dzfall: dzfall=avrate*dt.
c
c Passed through common (incomplete list)
c
c bbo(n)
c bwo(nd)
c d(nd),do(nd) :
c ddztp(nd) : Change in CV thickness of snow due to added 
c             precip [m snowdepth/s] = prcp/(3.6d0*bifall)                        
c dt: Current time step [s]
c dzo(nd): Old nodal thickness [m]
c H2Olayer(nd) : water based CV if .true.
c ipond(nd): For case of rain falling on impermeable layer,set ipond to 1
c        Note: Ponding disallowed in this version
c impermeable(nd): 1 = Node impermeable to water, 0 = Not.  See DENSITY
c        Note: Impermeable node disallowed in this version
c istart: Snow just started flag 1=yes 0=no
c k(nd) :
c n: Number of nodes in model
c newelm: Number of new nodes to be added to model
c newCV: same as above, but change terminology from node to CV
c nosnowcover: no snow cover flag 1= no snow cover 0=snow cover
c nn(ln)
c nold: Old number of nodes in model
c porosity(nd) :
c ss(nd): Effective water saturation = (s - ssi)/(1. - ssi)
c sso(nd): Old effective water saturation
c to(nd) : Old nodal temperature [K]
c u(nd) :
c unbar(nd) :
c us(nd) : Mass flux of snow [kg/m^2*s] us(n) = -rate*1000d0
c                                          or = -prcp/3.6d0
c 
c
      double precision dzfall,avrate
      integer m
c
c  Calculate rate of change in element thickness due to snow falling
c  during this iteration
c
cREJ_2025/04/29S      stop 'should not be hitting newsnow for Khuller project'
      nold=n
c March28, 1995   avrate=(rate+rateo)/2d0
      avrate=rate
c    Precip has just started or previous top node is full.  Initiate a
c    new node.
      if((istart .eq.1 .and. .not. repeat).or. fullnode)then
         newelm=1
         thinnode = .true.    !REJ_2025/08/30
         H2Olayer(n) = .true. !REJ_2025/06/21
         nosnowcover=0
         n=n+newelm
         ddzdtp(n) = 0d0 !REJ_2025/05/09 Initialized next 3 lines
         us(n) = 0d0
         dz(n) = 0d0
        if(n .gt.nd)stop '**Number of nodes has exceeded array size**'
         us(n)=-1d3*avrate  !REJ note. us() is positive upwards
         bifall=bwfall-blfall !REJ_CK  compute this from tair?
         ddzdtp(n)=avrate*1d3/bifall
         u(n+1)=-blfall*rate  ! liquid downward flux?
c         istart=0  ! REJ_2025 removed on 6/11. don't reset until old.f
         fullnode=.false.
         binew(n)=bifall ! Not really needed? !REJ_CK 2025/05/25
         bwo(n) = bwfall

c    Precip is continuing.  Add precip to top node.
      else
         newelm=0
         us(n)=-1d3*avrate
cS         write(*,*)'NEW:prcp,avrate,us(n),prcp/3.6',prcp,avrate,us(n),
cS     &     prcp/3.6d0
         bifall=bwfall-blfall
         ddzdtp(n)=avrate*1d3/bifall
         dzfall=ddzdtp(n)*dt
         
      end if
c
c  If top node is full (dz exceeds dzinc), initiate a new top node
c  for the NEXT iteration.
      if(dzfall+dzo(n) .ge. dzinc)fullnode=.true.
c
c  If there is a new node, initiate necessary parameters
c
      if(newelm .gt. 0)then
         k(n)=ln
         nn(ln)=nn(ln)+newelm
c      Next will cause new nodes to be retained in case of repeat iteration
c         nold=n
c
         write(80,*) 'New element ',n,' has been added'
cremove       Note that the new node is initialized with 0. bulk density and is
cremove       augmented by an amount u or us in sec. 13.
c        REJ_2025/05/25. Folowing block rewritten.  Partially based on newsnow
c        from SNTHERMP_2013
         to(n)=tprecip ! ARK_2025/10/7 to match OCT_07_ARK file
         bwo(n)=bwfall
         bbo(n)=0d0
         bifall=bwfall-blfall
         flo(n)=flfall
         if(blfall .lt. 9d2)then !REJ_2025/05/25 Don't follow this?? For ponding?
c        snow/firn is being added
           ss(n)=ssminsnow
           sso(n)=ssminsnow
           ssi(n)=ssisnow
           porosity(n)=1d0-bifall/dice
           unbar(n)=0d0
           u(n)=0d0
           u(n+1)=0d0
           do(n)=dsnowfall
           d(n)=dsnowfall
           dmass(n)= -us(n)*dt
           ci(n) =-13.3d0+7.8d0*Tprecip
         else
c        rainwater is ponding on top
            ss(n)=1d0
            sso(n)=1d0
            ssi(n)=0d0
            porosity(n)=1d0
            do(n)=dtol1
            d(n)=do(n)
c            u(n+1)=-1d3*rate  !REJ_2025/05/25 U is water flux??
            dmass(n+1)= -us(n+1)*dt
c          Next will cause unold in FILTRATE to be 0.0
            uo(n+1)=-1d3*rateo
            u(n)=0d0
            bw(n)=1d+3
         endif
c         write(*,*)to(n),bwo(n),td(n),td13(n),flo(n),blo(n),bi(n),
c     1        bt(n),dmass(n),bdjp(ln),a243(ln),bd(ln),a1(ln),
c     2        0,dzo(n),flgo(n),sso(n),dice,ssi(n),porosity(n),
c     &        dmvol(ln),ltype(ln),impermeable(n),idelete(n),
c     &        solidporosity(n),ipond(n),dicevol(n),dliqvol(n), 
c     &        rhowater,H2Olayer(n)
c         write(*,*)'END CALL'
c         goto 50  !REJ-CK  skip for now
         call density(to(n),bwo(n),td(n),td13(n),flo(n),blo(n),bi(n),
     1        bt(n),dmass(n),bdjp(ln),a243(ln),bd(ln),a1(ln),
     2        0,dzo(n),flgo(n),sso(n),dice,ssi(n),porosity(n),
     3        dmvol(ln),ltype(ln),impermeable(n),idelete(n),
     4        solidporosity(n),ipond(n),dicevol(n),dliqvol(n), 
     5        rhowater,H2Olayer(n)) !REJ_2025/5/29
cS        write(*,*)to(n),bwo(n),td(n),td13(n),flo(n),blo(n),bi(n),
cS     1        bt(n),dmass(n),bdjp(ln),a243(ln),bd(ln),a1(ln),
cS     2        0,dzo(n),flgo(n),sso(n),dice,ssi(n),porosity(n),
cS     3        dmvol(ln),ltype(ln),impermeable(n),idelete(n),
cS     4        solidporosity(n),ipond(n),dicevol(n),dliqvol(n), 
cS     5        rhowater,H2Olayer(n)
cS              stop 'newsnow'
50       continue
         melt(n)=0
         layertype(n)=soiltype(ltype(k(n)))
c        Next is a bandaid. For some reason, call density reset
c        this to F.  Then, thrk did not get called
         H2Olayer(n) = .true. !REJ_2025/06/21  Bandaid
c     For case of rain falling on impermeable layer, set ipond to 1
        if(blfall .gt. 9d2 .and. (ss(n-1) .gt. 1d0-1d-14 .or.
     &    impermeable(n-1).eq. 1 .or. sso(n-1).gt. 1d0-1d-14))then
          ipond(n)=1
        else
          ipond(n)=0
        endif
      end if
c      write(*,*)'End newsnow: us(n)',us(n)
c      write(*,*)'End newsnow:ddzdtp(n)', ddzdtp(n)
c      write(*,*)'End newsnow: n,bwfall,bwo(n)',n,bwfall,bwo(n),
c     & 'array(44,4)',array1(n,4)
      return
      end

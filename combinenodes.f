c***********************************************************************
c  COMBINENODES checks for elements which are below prescribed minima
c  for thickness or mass. It then determines which neighboring element
c  to best combine with, and executes the combination in routine COMBO.
c***********************************************************************
      subroutine combinenodes(dzmin,a1snow,thsnow,tlsnow,ssisnow,
     &  icalcstep,iunit) !REJ_2025/08/04 Added iunit
c
c  If snow element thickness is less than a prescribed minimum or
c  element has totally melted, combine with neighboring element.
c  Only one combination of elements is permitted per time step.
c  Idelete=1 if element has totally melted, which is determined
c  within DENSITY.
c  REJ_2025/01/19/ Changes for sntherm24 and beyond.  Layertypes
c  designated as "snow" in sntherm89 are designated as "H2Olayer" in
c  sntherm24, which contains the subtypes: snow, firn, and ice.  The   
c  best neighboring CV/element to combine with is now determined by CV
c  mass rather than thickness, without regard to the subtype.

c  Called from: MAIN
c
c  Calls: COMBO,DENSITY, function NMELT
c
      include 'const'
      include 'arrays'
      include 'l_arrays.txt' !REJ 6/24/24
c

c  Arguments
c
c dzmin : minimum nodal thickness (except for precipitation cases) [m]
c  Currently dzmin = 2d-3m.  Read in routine getinput.f
c a1snow : unfrozen water constant
c thsnow : meltzone upper limit temperatures [k]
c tlsnow : meltzone lower limit temperatures [k}
c ssisnow : irreducible water saturation for snow [.04]
c
      double precision dzmin,a1snow,thsnow,tlsnow,ssisnow
      integer icalcstep  !REJ remove after testing
      integer iunit !REJ_2025/08/04 Added iunit
c
c Local
c
c i: Nodal looping index.
c iar: Nodal index in common arrays.
c j: Node index.
c l: Node index.
c m: i+1 nodal array element.
c neighbor: Adjacent node selected for combination.
c
c Passed through common (incomplete list)
c
c dzo(nd) : Old elemental thickness [m]
c dmass(nd) : Elemental mass [kg/m^2]
c idelete(nd) : 1 signifies node to be removed from model 
c iskip(nd) : 1 signifies that convergence criteria are skipped for
c              this node
c ipond(nd) : Ponding flag.  Disallowed for this version
c n : Number of nodes in model
c nn : Number of nodes in each layer type, lowest to surface
c nsoil : Number of soil nodes
c nosnowcover : no snow cover flag 1= no snow cover 0=snow cover
c ltype(ld) : Layer type 1=snow; >99=user supplied soil layer
c         2-3=soil types contained in Block Data file soil
c         90-94 = firn/Ice.  see getinput for more details
c flo(nd) : Nodal fraction of liquid (unfrozen) water due to both
c       capillary and adsorbed potential
c k(nd) : Layer type assoicated with node
c nold : Old number of nodes in model
c ln : Number of different layers
c dt : Current time step [s]
c
c Functions
c
c nmelt: Function NMELT determines if element is in melt state.
c
      integer i,j,l,m,neighbor,nmelt,iar
      logical thin_or_lowmass_CV, do_combinenodes

c      do 50 i=nsoil+1,n,-1 !REJ_2025_06/25
      do 50 i = n,nsoil+1,-1 !REJ_2025_06/25
        do_combinenodes = .false.
        thin_or_lowmass_CV = .false.

        if(i .ge. n-1)then
c        Skip during snowfall
         if(i .eq. n .and. snowrate .le. 0d0)then !REJ_2025/08/24
c       Use usual combination criteria for CVs n and n-1
          if(dzo(i) .lt. dzmin .or. dmass(i) .lt. 80d0*dzmin)
     &    thin_or_lowmass_CV = .true.
CR          if(thin_or_lowmass_CV)write(iunit,*)'COMBO',
CR     &    'dzo(i),dzmin,dmass(i),80d0*dzmin,i',
CR     &     dzo(i),dzmin,dmass(i),80d0*dzmin,i  
         endif          
c          write(*,*)thin_or_lowmass_CV
c        if(icalstep .gt. 100)stop 'combine step 100'

        elseif(bi(i) .gt.  800d0)then  ! Recheck ice bi cutoff
c       Exclude ice lens from combination
          thin_or_lowmass_CV = .false.
        else ! i .lt. n-1
c       For Summit run increase combine threshhold from dzmin (2mm) to 1cm
         if(dzo(i) .lt. 0.01d0 .or. dmass(i) .lt.80d0*0.01d0)
     &   thin_or_lowmass_CV = .true. !REJ_2025/08/20     
      
        endif
c REJ_2025/08/23  next 3 out for now
c         if(icalcstep .ge.  1 .and. thin_or_lowmass_CV)
c     &   write(iunit,*)'COMBO1,i',i,thin_or_lowmass_CV,dzo(i),
c     &   dmass(i)
 
c REJ_2025/06/02 if( ((dzo(i) .lt. dzmin .or. dmass(i) .lt. 80d0*dzmin )
c REJ_2025/06/02  &  .and.((snowrate.le.0d0.and.ipond(i).eq.0).or.i.ne.n)).or.
c REJ_2025/06/02  &   idelete(i) .eq. 1)then

cREJ_2025/06/02 Replaced above 3 lines with the following IF block

         if(idelete(i) .eq. 1)then
            do_combinenodes = .true.
         elseif(i.ne. n)then
            if(thin_or_lowmass_CV)do_combinenodes = .true.
cREJ         if(icalcstep .ge. 6655)stop 'combo 6655'  'REMOVE'
cREJ     &     (snowrate.le.0d0.and.ipond(i).eq.0) )then
cREJ         elseif(i.eq. n .and.snowrate.ge.0d0 )then    
         elseif(i.eq. n .and.snowrate.lt.0d0 )then       
            if(thin_or_lowmass_CV)do_combinenodes = .true.
         endif      
c        if(icalcstep .ge. 6627 .and.do_combinenodes.and.i.eq.n)  !REMOVE
c     &    write(*,*)'snowrate',snowrate  !REMOVE
c         if(icalcstep .ge. 6627 .and. do_combinenodes)then
c           write(*,*)'got do_combinenodes',do_combinenodes,'step',
c     &      icalcstep, 'thin_or_lowmass_CV',thin_or_lowmass_CV,
c     &    'i',i,'del',idelete(i)  
c          if(icalcstep .ge. 6655)stop 'COMBO'  !REMOVE
c        endif
c         if(icalcstep .ge. 6655)stop 'COMBO'  !REMOVE
 
c      write(*,*)'do_combinenodes',do_combinenodes       
      IF(do_combinenodes)THEN     

c       stop 'Hit combinenodes. H2Olayer not programmed for this routine'  !reinstate c'd out
            dt=1d0
c     First determine which neighbor to combine with.  This is
c     generally the thinnest neighbor, except for the cases noted below. REJ_01/19
c     generally the lightest neighbor, except for the cases noted below.
c
c      If top node is removed, combine with bottom neighbor.
            if(i.eq. n)then
               if(nn(ln) .gt. 1)then
                  neighbor=n-1
               else
c             All H2Olayers are gone REJ_2025/01/19
                  n=n-1
                  nn(ln)=0
                  nosnowcover=1
c Nov 5, 1966
                  istart=1
        write(80,'(23h All H2Olayers are gone)') !REJ_2025/01/19
c March28, 1995   return
                  goto 99
               end if
c       If the bottom neighbor is not snow, combine with the top neighbor REJ_2025/01/19
c       If the bottom neighbor is soil, combine with the top neighbor
c            else if(ltype(k(i-1)) .gt. 1)then  !REJ_2025/01/19
             else if(H2Olayer(i-1) .neqv. .true.)then
               neighbor=i+1
c       If element is entirely melted, combine with bottom neighbor
            else if (flo(i) .ge. 1d0)then
               neighbor=i-1
c       If none of the above special cases apply, combine with the
c       thinnest neighbor
            else
               neighbor=i-1
c               if(dzo(i+1) .lt. dzo(i-1))neighbor=i+1 REJ_2025/01/19
               if(dmass(i+1) .lt. dmass(i-1))neighbor=i+1
            end if
            goto 60
         end if
 50   continue
c    No thin element was found.  return.   REJ_2025/01/19
c    No light CV was found.  return.
      return
c
c  node l and j are combined and stored as node j.
c
 60   if(neighbor .lt. i)then
         j=neighbor
         l=i
      else
         j=i
         l=neighbor
      end if
      iskip(j)=1
      nold=n
      n=n-1
c
c  Now combine
c
      if(idelete(i) .le.0)then
         write(80,*)icalcstep
         write(80,'(7h Nodes ,i4,4h AND,i4,14h combined into,i4)') l,j,j
      else
         write(80,'(6h Node ,i4,36h is totally melted and combined into,
     15h node,i4)')i,j
      endif
c
c     REJ_2025/08/23 Changed dmass to dmasso and bl to blo in    
c     following calls to combo and density. Also chg bi and bt??     
      call combo(dzo(j),dzo(l),to(j),to(l),bwo(j),bwo(l),dmasso(j)
     &   ,dmasso(l),bbo(j),bbo(l),dlm,ci(j),tlsnow,dice, a1snow,
     &   dsolo(j),dsolo(l),flo(j),do(l),do(j), H2Olayer, icalcstep,
     &   iunit)!REJ_2025/06/29+08/04
      nn(ln)=nn(ln)-1
      melt(j)=nmelt(to(j),thsnow,tlsnow)
      call density(to(j),bwo(j),td(j),0d0,flo(j),blo(j),bi(j),bt(j),
     1  dmasso(j),0d0,0d0,0d0,a1snow,0,dzo(j),flgo(j),sso(j),dice,
     2  ssisnow,porosity(j),0d0,ltype(k(j)),impermeable(j),0,
     3  solidporosity(j),ipond(j),dicevol(j),dliqvol(j),rhowater,
     4  H2Olayer(j))  ! REJ 6/24/24
 
c

c REJ_2025/08/07. It makes sense the most of the ICVPs are returned 
c as their past, or 'old', values--since routine old.f was called
c just prior to combinenodes.f. Dmass, however,is a puzzling exception,
c as the Dmasso omission causes erroneous bulk density values in the next
c call to the new massbal.f routine.  Therefore, I now return both dmass 
c and dmasso. REJ_2025/08/07.  Reason probably was that dmasso was not
c included in common in sntherm89. Added in this 2024 version.
c Next 3 lines no longer needed.
c      do 145 i = j,n  !New block setting dmasso = dmass for i .ge. j
c       dmasso(i) = dmass(i)
c 145   continue

c Now shift all elements above this down one.
      if(n .gt. j)then

      do 150 i=j+1,n
      m=i+1
      do 24 iar=1,narray1
         array1(i,iar)=array1(m,iar)
24    continue
      do 25 iar=1,narray2
         array2(i,iar)=array2(m,iar)
25    continue
      do 26 iar=1,niarray
         iarray(i,iar)=iarray(m,iar)
26    continue
150   continue
      end if
c  zero out arrays for top node
c March 28, 1995      do 20 iar=1,narray1
99    do 20 iar=1,narray1
      array1(nold,iar)=0d0
20    continue
      do 21 iar=1,narray2
      array2(nold,iar)=0d0
21    continue
      do 22 iar=1,niarray
      iarray(nold,iar)=0
22    continue
      uo(n+1)=uo(nold+1)
      u(nold+1)=0d0
      uo(nold+1)=0d0
cREJ_2025/04/29      stop 'combo'
      return
      end

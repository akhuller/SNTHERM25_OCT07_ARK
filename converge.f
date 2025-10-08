c**********************************************************************
c CONVERGE checks the convergence of the solution scheme with
c regards to linearization error in estimating temperature and
c nodal change in saturation
c**********************************************************************
      subroutine converge(repeat,igood,ngoodmin,dtsum,dssallowed,
     &   errtallowd,dtmin,dtmax,dtsmin,dtssmax,print,Neumann,!11/7/23 RJ
     &   icalcstep,iunit) !REJ_2025/06/25

      implicit none

      include 'const'
      include 'arrays'
      include 'l_arrays.txt'

c
c calling routine : main
c
c calls routines  : reset
c
c Argument list
c
c   repeat : logical flag denoting a repeat of this iteration (if true)
c   igood      : consecutive good calculations counter
c   ngoodmin   : minimum number of consecutive good calculations required
c              before increasing time step
c   iter       : if = 0 denotes new basic time period, not = 1
c   dtsum      : running total time limited to a maximum of dtbase
c   dssallowed : maximum allowed nodal change in saturation per
c               calculation time step
c   errtallowd: maximum allowed linearization error per time-step in heat balance
c               eq, expressed in units of temperature [K]
c   dtmin      : minimum allowed time step
c   dtsmin     : minimum allowed time step when water flow is present
c   dtssmax    : maximum time step for repeat calculation for
c               saturation criteria
c   dtmax      : maximum allowed time step
c   print      : print (if true) flag, usually at the hour
c   Neumann    : use a Neumann bottom boundary condition (if true) 11/05/23 RJ
c
c
c
c      integer igood,ngoodmin,nbot  ! 11/05/23  RJ
      integer igood,ngoodmin,nbot,icalcstep,idebug,l  !REJ_2025/06/23 + REJ_2025/07/25
      double precision dtsum,dtsmin
      double precision dssallowed,errtallowd,dtmin,dtssmax,dtmax,dum !RJ

      logical print,repeat,Neumann,debug ! 11/05/23  RJ + REJ_2025/07/25
c
c Local
c
c   convchk   : Flag used to test whether convergence test is needed 
c             : (convchk=idelete(i).eq.0 .and. iskip(i) .eq. 0).
c   dss       : nodal change in saturation.
c   dssmax    : maximum nodal change in saturation
c   deltamassl : time change in mass of liquid water, used in
c                temperature criterion.
c   cxmass    : dummy variable used in temperature criterion.
c   flux      : dummy variable used in temperature criterion.
c   errtempmax: maximum error in temperature estimation
c   errtemp   : error in temperature estimation
c   hconverge : logical flag denoting temperaure convergence met (if true)
c   h1rstwarn : logical flag denoting first temperaure convergence
c             : warning for this basic time step (if true)
c   i         : loop counter
c   nodeh     : node not meeting change in stored heat convergence
c   nodes     : node not meeting change in saturation convergence
c   notmet    : =true means convergence criterion not met.
c   ok        : logical flag denoting temperature and saturation convergence 
c             : met (if true).
c   ratio     : temperature error/ allowed temperature error
c   sconverge : logical flag denoting saturation convergence met (if true)
c   s1rstwarn : logical flag denoting first saturation convergence
c             : warning for this basic time step (if true)
c   tolerance : used in testing convergence criterion.
c   tolerance2: used in assuring that timestep hits dtbase
c
c Passed through common
c
c   dt        : Current time step [s]
c   dto       : Old time step [s]
c   idelete(nd) : 1 signifies node to be removed from model
c   iskip(nd) : 1 signifies that convergence criteria are skipped for
c               this node      
c   iskipo(nd): Past value for iskip
c   n         : Number of nodes in model
c   ct(nd)    : Nodal combined specific heat of all constituents [J/kg-K]
c   dbvdt(nd) :1)Change in saturation vapor density per degree K [kg/m^3 K]
c             :2)Change in bulk vapor density per degree K [kg/m^3 K]
c   dls       : Latent heat of sublimation of ice (2.838E6) [J/kg] at 273.15
c   dlm       : Latent heat of fusion for ice(3.335E5) [J/kg]
c   dz(nd)    : Elemental thickness [m]
c   dzo(nd)   : Old elemental thickness [m]
c   bl(nd)    : Nodal liquid water bulk density [kg/m^3]
c   blo(nd)   : Old Bulk density of liquid water [kg/m^3]
c   bt(nd)    : Nodal total bulk density [kg/m^3]
c   bb(nd)    : Nodal conducted and convected heat plus absorbed solar
c   bbo(nd)   : Old nodal conducted and convected heat plus absorbed solar
c   dmass(nd) : Nodal mass [kg/m^2]
c   us(nd)    : Nodal convective mass flux of snow [kg/m^2-s] 
c   ci        : Specific heat of ice  [J/kg-K]
c   tprecip   : Precipitation temperature [K]
c   to(nd)    : Old nodal temperature [K]
c   flfall    : Fraction of liquid water within falling precip.
c   flo(nd)   : Nodal fraction of liquid (unfrozen) water due to both capillary
c               and adsorbed potential
c   ss(nd)    : Effective water saturation = (s - ssi)/(1. - ssi)
c   sso(nd)   : Old effective water saturation
c a(nd,3) : Coefficient array in linear equation matrix  RJ
c b(nd) : Constant array in linear equation matrix  RJ
c qs(nd)  RJ
c
      logical hconverge,sconverge,ok,convchk,notmet,h1rstwarn,
     &s1rstwarn
      integer i,nodeh,nodes,iunit
      double precision tolerance,dss,deltamassl,tolerance2,cxmass,flux
      double precision dssmax,errtempmax,errtemp,ratio  
      double precision errtemp2 !REJ_2025/06/23    
      parameter (tolerance=0.99d-3,tolerance2=0.99d-6)
c
c      if(dtsum .gt. 2500) !These writes generate "confirm _dtbase_hit'
c     &   write(*,*)'start: dtsum,dt,dto',dtsum,dt,dto
      dto=dt
      dssmax=0.0
      errtempmax=0.0
      if(.not. Neumann)then  ! added block on 11/05/23
        nbot = 2
      else
        nbot = 1
      endif  
c      nodeh=2  11/05/23
c      nodes=2
      nodeh=nbot
      nodes=nbot
      ratio=0.0
      print=.false.
      if(iter .le. 1) then
        h1rstwarn = .true.
        s1rstwarn = .true.
      endif
c
c determine nodes causing worst convergence problems.
c
cc   do 10 i=2,n  !11/05/23
      do 10  i = nbot,n
c
c if current node is to be deleted (as if from melt) or combined with
c another, then convergence test skipped. if converge criteria are not
c met, converge will cause code to reduce time step until minimum step
c is reached.
c
         convchk = .false.  !REJ_2025/09/01        
         convchk=idelete(i).eq.0.and.iskip(i).eq.0.and.iskipo(i).eq.0!REJ_2025/09/01

c         if(icalcstep .gt. 559 .and. i .eq. n)
c     &      write(*,*)'convchk = ',convchk
c         if(i .eq. n.and. icalcstep .gt. 570)write(iunit,*)'conv.f',
c     &   convchk,idelete(i),iskip(i),iskipo(i)
c        REJ_2025/06/24  Next is a debugging block.  Remove when done.
c         if(.not. convchk)then
c          write(*,*)'convchk is F,i,iskip(i)', i,iskip(i)
c          stop 'converge.f'
c         endif
c          if(i .eq. n .and. icalcstep .gt. 570)write(iunit,*)'convchk',
c     &      convchk,icalcstep
c         if(i .eq. n .and. icalcstep .gt. 560)write(iunit,*)'convchk',
c     &     convchk

         IF( convchk ) THEN
c          Next checks temperature criteria
c redid next for clarity cxmass= (ct(i)+dbvdt(i)*dls/bt(i))*dmass(i)
c accurate to about 1d-13 !RJ
            cxmass= (ct(i)*bt(i)+dbvdt(i)*dls)*dz(i)

            flux=bb(i)+bbo(i)
cM            if(us(i) .lt. 0d0)flux=flux-us(i)*ci(i)*(tprecip-to(i))
cM     &         -us(i)*dlm*(flfall-flo(i))
            deltamassl=bl(i)*dz(i)-blo(i)*dzo(i)
            errtemp=dabs(dto*(flux-dlm*(deltamassl/dt+unbar(i))
     &        +unbar(i)*hs(i))/cxmass-(t(i)-to(i)))
c           errtemp2 is for "dry" snow
            errtemp2=dabs(dto*flux/cxmass-(t(i)-to(i)))
c**********************************************************************
c----------------------------------------------------------------------
c Print block for comparing surface arrays after thparam and within this
c converge routine. Written to converge_out. Arrays should be equal.
cT_RESTORE            if(i.eq.n.and. icalcstep .gt. 550) then 
C           if(i.eq.n) then  !c_REMOVE
C            write(iunit,*)'4. CONVERGE                       cxmass'
C            write(iunit,6)n,icalcstep,dt,dz(i),cxmass,bt(i),ct(i),
C    &            f(i),dbvdt(i)*dls,T(i),To(i)
C6           format(2i5,f10.3,f10.6,5f12.6,2f7.2)
C            write(iunit,*)'bb(n) = topfluxk+topfluxv*T(n)+dsol(n)',
C    &       'topflux',topfluxk+topfluxv*T(n),'dsol',dsol(n),
C    &       'bb(n)',bb(n),'bbo(n)', bbo(n)
C            write(iunit,*)'dt,dto',dt,dto
C            write(iunit,*)'errtemp2 = ',
C    &           'dabs(dto*flux/cxmass-(T(i)-To(i)))'
C            write(120,*)
C    &       'errtemp,errtemp2,T(n)',errtemp,errtemp2,n,T(n),icalcstep
C             write(120,*)
C           endif
c---------------------------------------------------------------------------
            if(errtemp.gt.errtempmax) then
               nodeh=i
               errtempmax=errtemp
            endif
c          Next checks saturation criteria
            dss= dabs(ss(i)-sso(i))
            if(dss.gt.dssmax .and.ss(i).gt.0d0) then
               nodes=i
               dssmax=dss
            endif
         ENDIF
 10   continue

      debug = .false.; idebug = 0

c compare actual error/change to allowed
c
      ratio=errtempmax/errtallowd
      hconverge= errtempmax .le. errtallowd
      sconverge= dssmax .le. dabs(dssallowed)
c
      ok = hconverge .and. sconverge
      notmet=.false.
c
      if( ok ) then
C     if(icalcstep .gt. 550)write(iunit,*)'ok iunit',iunit
c  ok convergence : increment the good calculation counter, estimate
c  the next time step, and keep going.
c
         repeat=.false.
         igood=igood+1
         iter=iter+1
         if(dabs(unbar(nodes)) .gt. 1d-7)then
c        This is a period of water infiltration
            dt=dmin1(dto,porosity(nodes)/(4d0*dabs(unbar(nodes))))
            dt=dmax1(dt,1d-2)
            if(dt .lt. dto)igood=0
          endif
c
c       Try increasing the time step if the number of consecutive good
c       calculations exceeds the minimum required.
         if(igood .ge. ngoodmin)then
c            if(icalcstep .gt. 550)write(iunit,*)'igood',igood
c March22, 1995            dt = dmin1( dabs(dto/ratio),1.5*dto,dtmax)
            dum=dt !dum
c         write(*,*)'BEF: igood,ngoodmin,dt',igood,ngoodmin,dt  !remove
c           REJ_2025 note: Next usually increases dt by 50%
            dt = dmin1( dabs(dto/ratio),1.5*dto,(1-iwet)*dtmax+
     &      iwet*dtssmax)
C           if(icalcstep.gt.550)write(iunit,*)'dabs(dto/ratio),1.5*dto'
C    &         ,dabs(dto/ratio),1.5*dto,'ratio',ratio
CR     &      ,(1-iwet)*dtmax+iwet*dtssmax,dto  ! remove
            dt=dmax1(dt,1d-2)
c November 13, 1996. Set minimum on qs for top node.  Increases accuracy
c of temperature prediction-important for stability correction.
            dt=dmin1(qs(n)*dto/4d0,dt)
            igood=0
         endif
      else
c Convergence criteria not met:
         igood=0
         repeat=.true.
         if(dto .le. dtmin+tolerance) notmet=.true.
      end if
c
c Convergence criteria not met at minimum time step. notify
c user and continue.

c***********************************************************************
c REJ_2025/07/25
c     optional block to check input variables for debug.
       l = nodeh 
       if(icalcstep .ge. 999999)debug = .true.         
       IF(debug)THEN
c        write(iunit,*)icalcstep,l,prcp,errtempmax,T(l),
        write(iunit,*)'CNV.f:',icalcstep,l,prcp,errtempmax,T(l),
     &  ct(l),dbvdt(l),bt(l),dmass(l),bb(l),bbo(l),us(l),ci(l),tprecip
c     &  ,flfall,flo(l)
CR         write(*,*)'step,dt,ertempmax,l',icalcstep,dt,errtempmax,
CR     &      errtallowd,l,ok
       ENDIf
c      if(icalcstep .ge. 700)stop  'converge 700'  
c***********************************************************************
      if(notmet)then
c         stop 'hit notmet'
         repeat=.false.
cREJ_2025/09/04         dt=dtmin
         write(*,*)'dt',dt
         dt=dmax1(dt,1d-2)
         iter=iter+1
CNOTE!!: Fix next statements later so that the max occurs for that basic
c time period are printed out (instead of first occurrance)
         if(h1rstwarn)then
           write(80,100)dto
         if(.not.hconverge) write(80,101) nodeh,errtempmax,errtallowd !REJ_2025/08/10
         h1rstwarn = .false.
         endif
         if(s1rstwarn)then
           write(80,100)dto
         if(.not.sconverge) write(80,102) nodes,dssmax,dssallowed
         s1rstwarn = .false.
        end if
      end if
c
      if(.not. repeat)then
          dtsum= dtsum + dto
c
c       insure that time step + running time hit dtbase
c
c      if(dtsum .gt. 2500) ! do not delete
c     & write(*,*)'midway: dtsum,dt,dto,print',dtsum,dt,dto,print

      if( dabs(dtsum-dtbase) .le. tolerance2 .or.dtsum.gt.dtbase) then
c          if(icalcstep .gt. 615)stop 'dtbase.  08/14  check hit'
          dtsum=0.0
          iter=0
          igood=0
          dt=dto
          dt=dmax1(dt,1d-2)
          print=.true.
      else
           dt = dmin1(dt, dabs(dtbase - dtsum + tolerance2) )
           dt=dmax1(dt,1d-2)
      endif
c      write(*,*)'end: dtsum,dt,dto,print',dtsum,dt,dto,print  !do not delete
c      write(*,*)
c      if(dabs(dtsum).lt. 1d-8)stop 'BJ converge'
c
      else
c
c convergence criteria not met.
c reduce time step and repeat calculation.
c
      if(.not.hconverge) then
         dt=dmin1(dabs(0.67*dto*1.0/ratio),0.67*dto)
c         write(iunit,*)'convergence not met. new dt',dt ! ARK_2025/09/26
      else
         if(dabs(unbar(nodes)) .gt. 1d-9)then
           dt=dmin1(dtssmax,0.67*dto,porosity(nodes)/(4d0*dabs
     &            (unbar(nodes))))
         else
           dt=dmin1(dtssmax,0.67*dto)
         end if
      endif
      dt=dmin1(dt,dtmax)
      dt=dmax1(dt,1d-2)
c
c reset appropriate conditions
c
c       write(*,*)'WARNING: Reset is being called in CONVERGE' !REJ_2025/09/04 ! ARK_2025/09/26
c       write(80,*)'WARNING: Reset is being called in CONVERGE' !REJ_2025/09/04 ! ARK_2025/09/26
cR      stop 'calling reset in CONVERGE'
       call reset

      endif
c
c       write(*,*)'dt & dtsum bottom of converge',dt,dtsum
      return
 100  format(/,'***convergence not obtained for time step of ',
     &1pe11.4,'sec in CONVERGE ***',/)
 101  format('   at node=',i4,' error in temperature est= ',1pe11.4,' >
     & ',' allowed= ',1pe11.4,'icalstep=  ',i8) !REJ_2025/08/10
 102  format('   at node=',i4,' change in saturation= ',1pe11.4,' > ',
     &' allowed= ',1pe11.4)
      end

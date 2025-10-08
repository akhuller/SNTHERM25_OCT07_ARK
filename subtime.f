c***********************************************************************
c Subroutine SUBTIME interpolates linearly between consecutive 
c meteorological data points, hourlymet_old and hourlymet
c REJ_2025/06/10  Checking that Summit halfhour step works OK.
c dtbase passing correctly as 1800
c Not called when this is the initial met step in the 
c***********************************************************************
      subroutine subtime(repeat,timesum,tkairsave,icalcstep)!REJ_2025/01/21+06/10

      implicit none
      include 'const'
      include 'arrays'
      double precision timesum,tkairsave  !REJ_2025/01/21
      integer icalcstep  !REJ_2025/06/10
c
c called from main
c
c arguments
c
c repeat : logical flag denoting a repeat of this iteration (if true)
c 
c delarations : !REJ_2025/01/22
c           Met variables declared in arrays, Common block PARAM
c
c REJ_2025/01/22.  Note on testing this routine. 
c Summed dt over the hourly met step.  It summed to 3600s, as it should.
c Saved value of newly read hourlymet value as hourlymetsave. Checked 
c that hourlymet = hourlymetsave when sum(dt) = 3600 seconds. The preceding is
c true, but only to a 0.1d0 accuracy.
      logical repeat 
      if(iter .le. 1 .and. .not. repeat)then
c       Next forces step function for tkair at start of precip.
         if(istart .eq. 1)tkairo=tkair
         da=(tkair-tkairo)/dtbase
         dr=(dirdown-dirdowno)/dtbase
         drh=(rh-rho)/dtbase
         dis=(solar-solaro)/dtbase
         dw=(wsp-wspo)/dtbase
         if(istart .eq. 1)tprecipo=tprecip
         dtprecip=(tprecip-tprecipo)/dtbase
c        Recheck the best way to interpolate precip
c        Past read-in precip value is stored , since it is
c        not hit bt interpoaltion routine due to built-in tolerance
         prcpo=prcpstore
         timesum = 0d0  !REJ_2025/01/21
         Tkairsave = Tkair
c  Next added on March 27,1995.  Eliminates interpolation of precip
         prcpo=prcp
         dprecip=(prcp-prcpo)/dtbase
         prcpstore=prcp
      end if
c      write(*,*)'timesum2',timesum
       timesum = timesum + dt
c      write(*,*)'timesum2',timesum
c      write(*,*)'iter,dt,da*dt',iter,dt,da*dt,tkair,timesum
cREJ_2025/06/03      if(timesum .ge. 3600)then
cREJ_2025/06/03  Passed test for dtbase = 1800
c      if(timesum .ge. dtbase)then !REJ_2025/06/23
c        write(*,*)timesum, tkair,tkairsave  !RESTORE THIS WRITE LATER
c        stop 'end subtime test'
c      endif
      tkair=tkairo+da*dt
      tprecip=tprecipo+dtprecip*dt
      wsp=wspo+dw*dt
      rh=rho+drh*dt
      dirdown=dirdowno+dr*dt
      solar=solaro+dis*dt
      solar=dmax1(0d0,solar)
      prcp=prcpo+dprecip*dt  !REJ_2025/06/26 Check this
      if(prcp .gt. 0d0)then
c        write(*,*)'FROM SUBTIME:prcp,dprecip',prcp,dprecip
      endif
      return
      end

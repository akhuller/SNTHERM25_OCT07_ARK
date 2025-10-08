***********************************************************************
c  WRITEFLUXMATLAB creates output file containing summary of surface
c  fluxes and meteorological parameters
c***********************************************************************
      subroutine writefluxmatlab(iy,jday,ihour,ibasestep,iptype,height
     &  ,clearness,imin,isolarcalc,effceiling,islope,mu_not,ffmt)
c %W% %G%
c
c Called from MAIN 
c
      implicit none
      include 'const'
      include 'arrays'
c arguments
c
c jday : Julian day of simulation, (read from met time-hack)
c ihour : Hour of simulation, (read from met time-hack)
c ibasestep: Current number of basestep from start of problem
c initial : Flag =1 for first calculation then set to 0
c iptype : precipitation type (from met data) 1= rain 2 =snow
c height : height above ground of measured met values
c iy : Year of simulation: last 2 digits, (read from met time-hack)
c clearness: approximate clearness factor for the sky
c imin: Minute of simulation, (read from met time-hack)
c isolarcalc: Estimate solar radiation flag 2=estimate missing values only 
c             1=yes 0=no
c islope : Adjust solar radiation for slope  1=yes  0=no
c effceiling: Effective cloud ceiling, computed in GETMETc 
c 
c July 10 1996 - added bert's output format  (jcm)
c ffmt: format type 1=default -1=bert's
c
      integer jday,ihour,ibasestep,iptype,iy,imin,isolarcalc,islope
      double precision height(3),clearness,effceiling,mu_not
      integer ffmt
c
c local
c
c dqest: net change in heat content of snowcover(dqest=solar+dlong+hg+
c        sensheat+dlatheat).
c i: index of array height.
c type: Precipitation type
c     :type = 'R' means rain.
c     :type = 'S' means snow.
c     :type = ' ' otherwise.
c doheading : Logical flag. T=print heading   F=do not print heading
c doformat945 : T=use format 945 for first print of hour=0 or 12 
c       F=use format 928.
c
c Passed through common
c
c solar: Net solar radiative flux [W/m^2]
c dlong: Net long wave radiative flux
c hg: Geothermal heat flux (across snow/soil interface) [W/m^2]
c sensheat: Turbulent flux of sensible heat [W/m^2]
c dlatheat: Turbulent flux of latent heat [W/m^2]
c prcp: Precipitation Value (m/hour)
c tkair: Air temperature [K]
c wsp: Wind speed [m/s]
c rh: Relative humidity
c
      double precision dqest
      
      

      
     
      write(49,691) iy,jday,ihour,imin,dqest,solar,dlong,hg,sensheat,
     &    dlatheat,prcp,tkair,wsp,rh,mu_not
 

691   format(i2,1x,i3,1x,i2,1x,i3,1x,f8.2,f8.2,f8.2,f8.2,f8.2,f8.2,1x,
     &       f8.2,f8.2,f8.2,f8.2,f8.2)

      
 
      return
      end

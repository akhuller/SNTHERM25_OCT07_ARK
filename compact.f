c***********************************************************************
c COMPACT calculates the  natural compaction rate of snow cover
c***********************************************************************
      subroutine compact(bi,t,bl,overburden,pdzdt,ss,dice,
     1 bwo,flo,floo,dto,unbar,dzo,dzoo,melt,dmlimit,eta0,binew, 
c    10/22/2023 RJ   2 i,n)
c    2 i,n,Khuller,icalcstep) !10/22/2023 RJ !remove icalcstep c REJ_2025/01/19
     2 i,n,ICVPfixd,icalcstep) !10/22/2023 RJ !remove icalcstep
c  SNTHERM89.REV4
cxx Routine adapted from e. anderson , p. 82-83.  Units have been changed
cxx to si system.  

c Routine changed extensively in Oct-Nov, 1995. The recommended coeffi-
c cients on the 'overburden' components are those used by Brun 
c et al. (1989), J. of Glac., 35,121.  The cut-off for compaction due
c to destructive metamorphism has been made an input.  The cut-off is
c taken as the smaller of this value and the density of the newly
c fallen snow compacted by 15%.   The cut-off for melt compaction 
c has been lowered to 350 kg/m3, except for the top node, where
c ablation alwalys occurs through a loss in snowdepth.   

c May 12, 1992.  Changes made to add compaction due to snow melt.
c
c CALLED from MAIN
c
c No CALLS to subroutines or functions.
c
c arguments
c
c bi : Nodal bulk density of ice [kg/m^3]
c t : Nodal temperature [K]
c bl : Nodal liquid water bulk density [kg/m^3]
c overburden : pressure of overlying snow [kg/m^2]
c pdzdt : Nodal rate of change in fractional-thickness due to 
c             compaction [fraction/s] Note. Passed back as pdcdtc
c ss : Effective water saturation = (s - ssi)/(1. - ssi)
c dice : density of ice (917) [kg/m^3]
c bwo : Old nodal bulk water density [kg/m^3]
c flo : Old nodal mass fraction of liquid water
c floo : Nodal mass fraction of liquid water from time period before last
c dto : Old time step [s]
c unbar : Average nodal convective mass flux of water [kg/m^2-s]
c dzo : Old elemental thickness [m]
c dzoo : Elemental thickness from time period before last [m]
c melt : Signifies if node in meltzone if = 1
c ICVPfixd : if .true., the inherent CV properties are fixed
c
      implicit none

      double precision bi,t,bl,overburden,pdzdt,ss,dice,dmlimit
      double precision bwo,flo,floo,dto,unbar,dzo,dzoo,eta0,binew
      integer melt,i,n,icalcstep
c      logical Khuller !REJ 10/22/2023
      logical ICVPfixd !REJ_2025/01/19
c
c local
c
c Now Variable c1: = 2.777d-7 m2/(kg s).
c c2: = 21d-3 m3/kg.
c changed above t0 value used by Brun et al. 23d-3 m3/kg
c c3: = 2.777d-6 1/s.
c c4: = 0.04 1/k.
c c5: = 2.0. 
c c6: = 5.15d-7.
c c7: = 4d0.
c ddz1: Rate of settling of snowpack due to destructive metamorphism.
c ddz2: Rate of compaction of snowpack due to overburden.
c dexpf: dexpf=dexp(-c4*(273.15-t)).
c dm: Critical density.  Input. Above this, settling function is damped
c nodalmelt : Water flux resulting from nodal melt  [kg/m^2-s]
c ddz3 : Rate of compaction of snowpack due to melt [1/s]
c
      double precision c2,c3,c4,c5,dm,ddz1,ddz2,c6,c7,dexpf
      double precision nodalmelt,ddz3
      data c2,c3,c4,c5/23d-3,2.777d-6,0.04,2.0/
c     data c1,c2,c3,c4,c5,dm/2.777d-7,21d-3,2.777d-6,0.04,2.0,0.15d3/
      data c6/5.15d-7/,c7/4d0/

c Disallow compaction for ice and for water saturated node.
      if(bi  .ge. dice .or. ss .ge. 1.)return
c REJ_2025/01/19      IF(.not. Khuller)THEN  !11/11/23  RJ
      IF(.not. ICVPfixd)THEN  !11/11/23  RJ
        dexpf=dexp(-c4*(273.15-t))
c Settling as a result of destructive metamorphism
        ddz1=-c3*dexpf 
        dm=dmin1(dmlimit,1.15d0*binew) ! REMOVE for now !REJ_2025/08/26
c       REJ_2025/06/19. Note that the Summit KM run
c       is intialized with Bw > dm
c REJ_2025/06/09. Check out the next damping (dm) term with
c published data, other models.
        if(bi .gt. dm) ddz1=ddz1*dexp(-46.0d-3*(bi-dm))

c Liquid water term
        if(bl .gt. 0.01)ddz1=ddz1*c5
c Compaction due to overburden
c       ddz2=-overburden*c1*dexp(-0.08*(273.15-t)-c2*bi)
        ddz2=-overburden*dexp(-0.08*(273.15-t)-c2*bi)/eta0
      ENDIF
c Compaction occurring during melt
      if(melt .gt. 0)then  ! does this exclude small liqvol?
       nodalmelt=bwo*(flo-floo)*dzo/dto+unbar*(1d0-floo)
c      Changed next on Nov 30, 1992 RJ
       if(nodalmelt .gt. 0d0 .and. bi*dzoo .gt. 0d0
     &    .and.(bi .lt. 250d0 .or. i .eq.n))then
         ddz3=-nodalmelt/(bi*dzoo)
CR         write(46,*)'nodalmelt',nodalmelt,bi
       else
         ddz3=0d0
       endif
      else
       ddz3=0d0
      endif
c Time rate of fractional change in dz (units of s-1)
c       if(icalcstep .ge. 567020)
c REJ_2025/01/19      if(.not. Khuller)then  ! Change in next 5 statement
      if(.not. ICVPfixd)then  
       pdzdt=ddz1+ddz2+ddz3  ! 11/11/23  RJ
      else
       pdzdt = ddz3
      endif
c Next used by Khuller to disallow compaction for snow density < 550 kg/m3 REJ_2025/01/19
c Next used by ICVPfixd to disallow compaction for snow density > 550 kg/m3
c Ask if he used bi or bw? Or are we omitting compaction altogether? Yes
c      if(Khuller .and. bi .gt. 550d0)pdzdt = 0d0  ! 10/22/2023 RJ
       if(ICVPfixd .and. bi .gt. 550d0)pdzdt = 0d0  ! 10/22/2023 & REJ_2025/01/19
     &     !added dz change for melt alone.
      return
      end

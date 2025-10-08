c***********************************************************************
      Subroutine massbal(icalcstep,ICVPfixd,dtmin,dzmin,bwfall,
     &     sum_vapor_dt,totaltime,sum_delta_dz,thinnode)!REJ_2025/07/22
c      Last updated 08/30/2025
c***********************************************************************

c
c Called from MAIN.  Section 13.
c Calls None
c
      implicit none
      include 'const'  
      include 'arrays'
      include 'l_arrays.txt'

c Arguments
      integer icalcstep
      double precision dtmin,dzmin,bwfall,sum_vapor_dt  
      logical ICVPfixd,thinnode
      
c Passed Through ARRAYS
c      bw(nd)
c      bwo(nd)
c      bl*nd)
c      dz(nd)
c      dzo(nd)
c      k(nd)
c      porosity(nd)
c      us(nd)
c      unbar(nd)
c      uvapor(nd)
c      ss(nd)
c      idelete(nd)
c      iskip(nd)
c      snow_age(nd)  !REJ_2025/07/24
c      n,nold,nsoil,dt,dto,prcp,istart
      
c Local 
      double precision Ddzprcp_rate,Ddzvapor_rate,Ddzcompact_rate 
      double precision bwvflux,Ddz_rate(nd)
      double precision mass_changes,totaltime  !REJ_2025
      integer i,m,icalcstep_print
      double precision sum_delta_dz  !REJ_REMOVEW

c Functions
      double precision fliquid, fCV_thickness
c
c IMPORTANT SIGN CONVENTION.  CONTRARY to the usual SNTHERM 
c POSITIVE DOWNWARDS THERMAL convention, the MASS flux convention
c IS POSITIVE UPWARDS. This will be reversed as time permits.
      
      DO 395 i=n,nsoil+1,-1 !CV LOOP. No soil treatment, for now

c     INITIALIZE
 
        Ddzprcp_rate = 0d0; Ddzvapor_rate = 0d0; Ddzcompact_rate = 0d0
        Ddzvapor_rate = 0d0; bwvflux = 0d0; Ddz_rate(nd) = 0d0
        mass_changes = 0d0
        if(icalcstep .eq.1)sum_delta_dz = 0d0 !REJ_2025/08/27 Testing


      IF( i .eq. n)THEN

c     PRELIMINARY COMPUTATIONS FOR SURFACE CVs !REJ_2025

c       MASS CHANGES  (kg/m2)
c       Mass changes due to snowfall, rainfall, and vapor exchange 
c       at the surface are accommodated by changes in CV thickness,
c       with bulk density specified according to water phase and
c       flux direction as listed below.
c         
c                               Bulk density     [kg/m3]
c         Snowfall:             bwfall (read in getmet.f)
c         Rainfall              1000 
c         Vapor_down, no prcp   hoarfrost 50?, else bwfall
c         Vapor_down, precip    bwfall
c         Vapor_up              bwo(n)

         If(uvapor(n) .le. 0d0)Then             !Vapor Downward
c          Cutoff 1d-10 is arbitrary            !Condensing
           if(prcp .lt. 1d-10)then !no prcp             
              prcp    = 0d0; us(n) = 0d0         
              bwvflux = bwfall !Otherwise
c             Other options can be added here
           else     
              us(n)   = -prcp/3.6d0
              bwvflux = bwfall
           endif
         Else                                   !Vapor Upward
              bwvflux = bwo(n)                  !Evaporating 
         Endif

c        SURFACE CV THICKNESS CHANGE (m)
        
           if(us(n) .lt. 0d0)then               !Snowfall
            Ddzprcp_rate = us(n)/bwfall 
           else      
            Ddzprcp_rate = 0d0              
           endif                                
          
           Ddzvapor_rate = uvapor(n)/bwvflux   !Vapor phase chg                                         

           Ddzcompact_rate = pdzdtc(n)*dzo(i)  !Compaction (-)

c         Original exponential solution. 
c         Above solution is first order Taylor's expansion.
c         dz(n) = fCV_thickness(dzo(n),pdzdtc(n),Ddzprcp_rate,
c     &          Ddzvapor_rate,i,dt)

c        FINAL COMPUTATIONS FOR SURFACE CV (New CVs treated differently)

c         NEW SURFACE CV INITIATED WHEN SNOWFALL STARTS
           If(istart .eq. 1)Then
c          Mass Balance
             Mass_changes = -(us(n)+unbar(n)+uvapor(n))  
cREJ_2025/07/25             dmass(n)     =  dmasso(n) + Mass_changes * dt
             dmass(n)     =  0d0 + Mass_changes * dt

c          CV thickness change.  
c          Compaction not implemented for istart = 1

             dzo(n) = 0d0  ! Case for istart = 1
             dz(n) = dzo(n) - (Ddzprcp_rate +Ddzvapor_rate)*dt
             bw(n) = dmass(n)/dz(n) !REJ considering a harmonic mean

         Else

c         REGULAR SURFACE CV with or without continuing precip.
c         Mass Balance
           Mass_changes = -( us(n)+unbar(n)+uvapor(n))
           dmass(n)  = dmasso(n) + Mass_changes * dt

c         Snow/Firn Compaction rate                             
            Ddzcompact_rate  = pdzdtc(n)*dzo(n)

c         Total CV thickness change rate
            Ddz_rate(n) = Ddzcompact_rate - Ddzprcp_rate - Ddzvapor_rate
            dz(n) = dzo(n) +  Ddz_rate(n) * dt

c         Update CV bulk density

           bw(n) = dmass(n)/dz(n)   !REJ added June 15,2025

        Endif

c        New CVs where evaporation exceeds precipitation are flagged
c        for deletion and subtracted from the adjacent lower CV

            newelm = 0  !REJ: track down where this is used!
c  REJ.  Will include next section later after testing
c           If(istart .and. dmass(n) .lt. dtol2)Then
c           If(newelm .eq.1 .and. dmass(n).lt. dtol2)Then
c             write(80,*)'New'// 'snow'// 'is'// 'evaporating'// 'or',
c     &      'melting faster'// 'than it'// 'is'// 'accumulating. '//
c     &      'Flag for deletion and combine with lower neighbor' 
c             dz(n)=dzo(i)+ddzdtp(n)*dt MB
c             idelete(n)=1
c             iskip(n)=1
c             goto 61 
c           Endif
 
       ELSE

c      INTERIOR CVs        

c        Mass balance. Adjust water mass for liquid water and vapor flow

          us(i) = 0d0 !REJ note. While this is an array, it doesn't need to be.
cREJ_2025/10/06          unbar(i) = 0d0 !water flow needs more testing before inclusion

         Mass_changes = -( unbar(i)+uvapor(i)) !REJ_2025/10/06
c          Mass_changes = - uvapor(i)         
          dmass(i)  = dmasso(i) + Mass_changes * dt

c         Change in CV thickess, depending only on compaction.

          Ddzvapor_rate  = 0d0; Ddzprcp_rate = 0d0
          Ddzcompact_rate = pdzdtc(i)*dzo(i)
          Ddz_rate(i) = Ddzcompact_rate
          dz(i) = dzo(i) + Ddz_rate(i)*dt 
          if(i .eq. 42)then   !REJ_REMOVE this testing block
           goto 555
           write(*,*)'dz,dzo,Ddz_rate(42)',dz(42),dzo(42),dt,
     &         Ddz_rate(42)*dt,icalcstep !CR
           write(*,*)'sum_delta_dz_old',sum_delta_dz
           sum_delta_dz = sum_delta_dz + dz(i)-dzo(i)
           write(*,*)'sum_delta_dz_new,dz(i),dzold',i,sum_delta_dz,
     &     dz(i)-dzo(i),Ddz_rate(i)*dt
           write(*,*)'start - final',0.02d0-dz(42)
           write(*,*)'Ddz',dz(i)-dzo(i),'rate',Ddz_rate(i)*dt,dt
555        continue
          endif
c           if(icalcstep .eq. 676)stop 'MB 676'
c          if(icalcstep .eq. 123)stop 'MB 123'

c         Original exponential solution. 
c         Above solution is first order Taylor's expansion.
c         dz(n ) = fCV_thickness(dzo(i),pdzdtc(i),Ddzprcp_rate,
c     &          Ddzvapor_rate,i,dt)

          bw(i) = dmass(i)/dz(i)   !REJ added June 15,2025 
        ENDIf
                 
c WRITE OUTPUT. OUTPUT FILE 115 IS MASSBAL_OUT
      icalcstep_print = 999999999 !Large integer for no output 
c     Print Header
      if(icalcstep .eq. icalcstep_print .and. i .eq. n)then
        write(115,*)
     &  'i icalcstep    dz(i)        dmass(i)       dmasso(i)     diff   
     &    uvapor(i)*dt   us(i)*dt Ddzcompact_rate*dt', !REJ_2025/08/24 
     & ' Ddz_rate(i)*dt  bw(i)'  !REJ_2025/08/24 
      endif
        if(icalcstep .ge. icalcstep_print)then
         if(i .eq. n)then         
          write(115,'(i4,i7,e15.7,2f13.6,5e14.6,e14.6)' ) 
     &    i,icalcstep, dz(i), dmass(i), dmasso(i), dmass(i)-
     &    dmasso(i),uvapor(i)*dt,us(i)*dt,Ddzcompact_rate*dt,    
     &    Ddz_rate(i)*dt,bw(i) !REJ_2025/08/24 added '*dt to above'
         else
         If(i .gt. 41)
     &    write(115,'(i4,i7,e15.7,4f13.6,14x,2e14.6,e14.6)' ) 
     &    i,icalcstep, dz(i), dmass(i), dmasso(i), dmass(i)-
     &    dmasso(i),uvapor(i)*dt,Ddzcompact_rate*dt,Ddz_rate(i)*dt,
     &    bw(i)
         endif
       endif
c      if(icalcstep .gt. 9)stop 'MB'
c Next section does a totaled mass balance.  added by REJ_2025/07/20
        if (i .eq. 42)then
          if(icalcstep .eq.1)sum_vapor_dt = 0d0
          sum_vapor_dt = sum_vapor_dt + uvapor(i)*dt  !REJ_2025
          write(125,'(4g28.20,f7.2,g16.9,i6)')uvapor(i),-sum_vapor_dt,
c     &    bw(i)*dz(i)-5d0,-sum_vapor_dt-bw(i)*dz(i)+5d0,bw(i),   !CV41
     &    bw(i)*dz(i)-3d0,-sum_vapor_dt-bw(i)*dz(i)+3d0,bw(i),   !CV42
c     &    bw(i)*dz(i)-1d0,-sum_vapor_dt-bw(i)*dz(i)+1d0,bw(i),   !CV43 
     &    totaltime/3600d0,icalcstep 
        endif

 395  Continue  !End CV loop
      if(icalcstep .ge. icalcstep_print)
     & write(115,*)!Skips a line for readability'

c***********************************************************************
cc****************************NOT ACTIVATED*****************************
c   Next block computes changes in dz and bulk water density due
c   to melt and water flow. Not yet activated   
c  
c!REJ_2025/10/06 goto 999  !Temporarily skips liquid water computation.  
c      Disallow for water saturated  or totally solid CVs   
c          If((ss(i).lt.1d0.and.porosity(i).gt.0d0))Then
c            If(i .eq. n)Then !For Surface CV
c           Change due to snowfall   
c               if(dabs(us(n)).gt. 0d0)dz(n)=dzo(n)+ddzdtp(n)*dt 
c               if(dabs(us(n)).gt. 0d0)write(*,*)'dz(n),bifallin,bw(n)'
c     &          ,dzo(n),dz(n),ddzdtp(n)*dt,bifallin,bw(n)
c              For top snow node, adjustments to mass because of
c              vapor flux are made by changing thickness dz rather
c              than bulk density. 
c              dz(n)=(wmass(n)-dt*uvapor(n))/bw(n)
c            Endif
 
836   format(i6,i9,5e18.7)   

        
c      Compute liquid water fraction  
         m = k(i) !REJ_2025/10/06
         bl(i)=fliquid(bw(i),td(i),bdjp(m),a243(m),td13(i),
     &    flgo(i)) !REJ_2025/10/06
c      Compute total mass and density, bd, for soil

          if(.not. H2Olayer(i))then !Treatment for soil
            bt(i)=bw(i)+bd(m)
            dmass(i)=bt(i)*dz(i)
          endif
999     Continue
c*************************END NOT ACTIVATED***********************************
     
61    if(prcp .le. 0d0.and.dmass(n) .lt. 80d0*dzmin)iskip(n)=1  !REJ_2025/06/28
c      REJ_2025 added the no precip condition to above on 06/28.  No longer
c      triggered iskip(n) = 1, although the surface CV is still set to thinnode.
c     REJ_2025.  I need to recheck the following section.
       goto 9999  !REJ_2025/07/23  Taken out for now
      if(iskip(n) .eq. 1)then
        write(80,*)'WARNING: skip found'  !REJ_2025/07/23
        write(80,*)'MB:step',icalcstep,iskip(n),'dzmin',dzmin,'80*',
     &      80*dzmin,'dmass(n)',dmass(n),'prcp',prcp
      endif
9999   continue
      if(dz(n) .lt. 0d0)idelete(n)=1 !REJ_2025/06/24 Added 61
        thinnode=.false.
        if(dmass(n) .lt. 80d0*dzmin)thinnode=.true.
c1.27.03  Without next, wrong gv and gk switches were selected
        if(iskip(n) .eq. 1)melt(n)=0 !1.27.03  REJ_2025  Recheck this.
            
      RETURN
      END
c***********************************************************************
      DOUBLE PRECISION FUNCTION 
     &        FCV_THICKNESS(dzo,pdzdtc,Ddzprcp_rate,Ddzvapor_rate,i,dt)
c     Function solves dz/dt = (Az + B) by substitution of Z = Az + B,
c     z = dz (CV thickness), A = pdzdtc, B = -us/bwfall - uvapor/bwvflux
c     dz = {exp(A*dt)*[A*dzo + B] - B}/A 
c***********************************************************************
      implicit none
C Arguments
c
      double precision dzo,pdzdtc,Ddzprcp_rate,Ddzvapor_rate,dt
      integer i
c Local
      double precision A,B

            if(dabs(pdzdtc) .lt. 1d-20)then
              write(115,*)'pdzdtc = ',pdzdtc,i  
              stop 'near zero pdzdtc found in fCV_thickness'
            endif
            A  = pdzdtc
            B  = -Ddzprcp_rate  -Ddzvapor_rate 
      RETURN
      END



c***********************************************************************
c  DIFFUSION computes mass diffusive vapor flux within snowcover.
c  pp.8-9 CRREL SR 91-16
c***********************************************************************
c      subroutine diffusion(de0,bp,bvi0,dlsdrw,frh,dlvdrw,bvw0,Khuller) !REJ 1/5/2024
      subroutine diffusion(de0,bp,bvi0,dlsdrw,frh,dlvdrw,bvw0,ICVPfixd,
     &    H2Olayer,icalcstep) !REJ_2025/05/18
c
c Called from MAIN
c
      implicit none
      include 'const'
      include 'arrays'
c
c Arguments
c
c de0   : effective diffusion coefficient for water vapor in snow  
c          [m^2/s] at 1000mb and 273.15k  !Set to 0.9d-4 in MAIN
c de    : effective diffusion coefficient for water vapor in snow
c bp    : approximate average barometric pressure [mb]
c Next 4 constants calculated in calconstant !REJ_2025/07/21
c bvi0  : constant for vapor diffusion, with respect to ice
cREJ_2025/07/27        =100.*e0*dexp(dlsdrw/273.15)/rw
c                      =4.85d-3*dexp(dlsdrw/273.15)/rw = 4.5234e3
c bvw0  : constant for vapor diffusion, with respect to water
cREJ_2025/07/27        =100.*e0*dexp(dlsdrw/273.15)/rw
c                      =4.85d-3*dexp(dlsdrw/273.15)/rw = 6.3597e4

c dlsdrw: constant for vapor diffusion dls/rw   ! = dls/rw
c dlvdrw: constant for vapor diffusion dlv/rw   ! = dlv/rw 
c frh   : Fractional saturation of water vapor  ! Set in MAIN
cc logical Khuller REJ_2025/01/11
c H2Olayer  : Water material  REJ_2025/05/18
      double precision de0,bp,dlsdrw,bvi0,frh(ld),vaporvol,dlvdrw
      double precision bvw0
      logical ICVPfixd,H2Olayer(nd) !REJ_2025/05/18
c
c Local
c
c de: diffusion coefficient
c i: looping index.
c
      integer i,icalcstep  !REJ_2025/06/04
      double precision de,dum,twodz
c
c Passed through common (incomplete)
c
c df(nd)     : Coefficient for vapor diffision
c              = de*change in sat. vapor pressure with temp
c uvapor(nd) : 1) diffusive mass flux of water vapor for CV interfaces [Kg/m^2 s]
c            : 2) Net diffusive mass flux of water vapor for CV [Kg/m^2 s]
c ufvapor(nd : "Average" diffusive flux at nodal centroid [(upper+
c              lower)/2]. Used in computation of graingrowth.[Kg/m^2 s]               
c dbvdT(nd)  : 1)Variation in saturated vapor density with temp 2) Variation
c        in saturated bulk vapor density with temp=vaporvolr*(1)  [kg/m^3 s]
c vaporvol(nd): volume fraction available for vapor
c dum         : = qlat when icalcstep = 1, else = qlato
c dum         : temporarily store the lower CV interfacial vapor flux

c  Vapor exchange at Surface of H2O Material or Bare Soil
c    (Vapor flow within soil disallowed at this point)
c
c

      if(icalcstep .eq. 1)qlato = qlat !REJ_2025/06/04 Otherwise, qlato = 0d0
      dum = 8.4879831653432862d-9 ! ARK_2025/09/10
c     REJ 06/-4/2025  recheck next
c   Vapor flow across Top Surface Boundary
c REJ note July 26,2025:  does not look like "dum" is passed in. Check this out

      if(dum .eq. 0)then
c               
        uvapor(n+1)=0d0
        write(80,*)'qlat(o) = 0d0 found in Diffusion' 
        stop 'diffusion: qlat(o) = zero'
      else
        if(dliqvol(n) .gt. 0.02. or. t(n) .gt. 273.15)then !wrt water
         uvapor(n+1)=-(eao-frh(k(n))*eso)*qlato/dlv
        else                                               !wrt ice
         uvapor(n+1)=-(eao-frh(k(n))*eso)*qlato/dls 
         endif
      endif
c      write(*,*)'qlato,uvapor(n+1)',qlato,uvapor(n+1),frh(k(n))
c       write(*,*)'de0,bp,bvi0,dlsdrw,frh,dlvdrw,bvw0',
c     &          de0,bp,bvi0,dlsdrw,frh,dlvdrw,bvw0
c       stop 'diffusion'

c  Vapor exchange at Interior Interfaces of H2O CVs 
c  REJ_2025/07/27 Replaced Saturation Vapor Pressure with
c                          Saturation Vapor Density

c
c      do 10 i=1,n !REJ_2025/05/18 !REJ_2025/05/18
c     REJ note 07/21/25.  Tried i = 1,n again. Gave identical results
      do 10 i=n,1,-1  !REJ_2025/05/18
         if(ss(i) .lt. 1d0 .and. porosity(i) .gt. 0d0)then 
c              Include a minimum porosity in the above
      
c       Compute diffusion coefficient
c        if(ltype(k(i)) .eq. 1)then !REJ_2025/05/18
         if(H2Olayer(i))then !REJ_2025/05/18
            de=de0*(1000./bp)*(to(i)/273.16)**6d0
         else
          return
c        For soil.  From Farouki, p. 55 Saving this for future reference
c        Note:  Soil option is not included in this version
c           de=0.161d-5*porosity(i)*(1000./bp)*(to(i)/273.)**2.3d0
         endif
         if(H2Olayer(i))then  !REJ_May20
c        if(ltype(k(i)) .eq. 1)then !REJ_May20
          if(dliqvol(i) .gt. 0.02)then
           dbvdT(i)=bvw0*dexp(-dlvdrw/to(i))*(dlvdrw/to(i)-1.0)/to(i)**2
          else
           dbvdT(i)=bvi0*dexp(-dlsdrw/to(i))*(dlsdrw/to(i)-1.0)/to(i)**2
          endif
           df(i)=de*dbvdT(i)
           vaporvol=solidporosity(i)-dliqvol(i)
           dbvdT(i)=vaporvol*dbvdT(i)
         else
           dbvdT(i)=0d0
           df(i)=0d0
           vaporvol=solidporosity(i)-dliqvol(i)
c           write(*,*)'vaporvol,solidporosity',i,
c     &        vaporvol,solidporosity(i)
           dbvdT(i)=vaporvol*dbvdT(i)
         end if
      else
c  vapor diffusion cannot occur within a water saturated element
        df(i)=0d0
        dbvdT(i)=0d0
      end if
10    continue
c
c  Vapor flux across lower CV boundary
c      do 20 i=nsoil+1,n
c      do 20 i=n,1,-1
       do 20 i=n,2,-1 !ARK_2025/09/25 Changed to 2 so we avoid df(0)
      if(df(i).le.0d0 .or. df(i-1).le.0d0 .or. ltype(k(i)).gt.1)then
     &    !RECHECK only for snow?
      uvapor(i)=0d0
      else
c REJ_2025/07/27 out      twodz = (dz(i)+dz(i-1)) !REJ_2025_May 20
c REJ_2025/07/27 OUT returned to ealrier form: Next is a significant 
c REJ_2025/07/27      uvapor(i)=twodz *  !REJ_2025_May 20 !Harmonic mean
c     & df(i-1)*df(i)*(to(i)-to(i-1))/(dzo(i)*df(i-1)+dzo(i-1)*df(i))
c REJ_2025/07/27 Restored original form   
       uvapor(i)=
     &2d0*df(i-1)*df(i)*(to(i)-to(i-1))/(dzo(i)*df(i-1)+dzo(i-1)*df(i))
c      if(i.ge.40 .and. i .le. 43)
c     &   write(*,*)'icalcstep,uvapor(i)',i,icalcstep,uvapor(i)
      end if
20    continue
c
c  Vapor flux across upper nodal boundary, net flux and "average" flux
c  at centroid (used in grain growth algorithm)
      do 30 i=nsoil+1,n-1
c      if(df(i).gt.0d0 .and. df(i+1).gt.0d0 .and. ltype(k(i)).eq.1)then
      if(df(i).gt.0d0 .and. df(i+1).gt.0d0 .and. ltype(k(i)).eq.1)then
      dum=
c REJ_2025/07/27     &twodz*df(i)*df(i+1)*(to(i+1)-to(i))/
     & 2d0*df(i)*df(i+1)*(to(i+1)-to(i))/ ! REL_2025/07/27
     & (dzo(i+1)*df(i)+dzo(i)*df(i+1))
       ufvapor(i)=0.5d0*(dabs(uvapor(i))+dabs(dum))
       uvapor(i)=uvapor(i)-dum
c Testing       if(i .eq. 42)write(*,*)'uvapor(i),dum',uvapor(i),dum 
c Testing       if(i .eq. 42 .and. icalcstep .gt. 2)stop      
      end if
30    continue
c Note:mass vapor flux across air interface is (eao-frh*eso)*qlato/dls
      if(dliqvol(i) .le. 0.02)then
        dum=(eao-frh(k(n))*eso)*qlato/dls
      else
        dum=(eao-frh(k(n))*eso)*qlato/dlv
      endif
      ufvapor(n)=0.5d0*(dabs(uvapor(n))+dabs(dum))
      uvapor(n)=uvapor(n)-dum
c      write(*,*)'uvapor(n)',uvapor(n)  !REJ_testing
c      stop 'diffusion' !REJ_testing
c        do i = n,20,-1
c              write(*,*)i,uvapor(i),df(i),T(i)
c        enddo

      return
      end


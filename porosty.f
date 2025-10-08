c***********************************************************************
c  Subroutine POROSTY calculates porosity of snow/soil medium
c***********************************************************************
      subroutine porosty(dicevol,dmvol,H2Olayer,porosity,
     &  solidporosity)!REJ6/24/24 !Passed as scalars
c
c Called from density
c
c arguments
c
c dicevol:  Volume fraction of ice
c dmvol : Fractional volume of dry soil material
c ltype : Layer type 1=snow >=90=user supplied soil layer
c                2-3=soil types contained in Block Data file soil
c porosity : fractional volume of voids: between ice matrix in
c                snow and between dry matrix in soil
c solidporosity : Volume of voids between dry soil and ice
c
      implicit none
      double precision dicevol,dmvol,porosity,solidporosity
      logical H2Olayer
c
cREJ 6/24/24      if (ltype .le. 1)then
      if (H2Olayer)then
        porosity=1d0-dicevol
        solidporosity=porosity
      else
        porosity=1d0-dmvol
        solidporosity=porosity-dicevol
      endif
      if(porosity .gt. 1d0)porosity=1d0
      if(porosity .lt. 0d0)porosity=0d0
      if(solidporosity .gt. 1d0)solidporosity=1d0
      if(solidporosity .lt. 0d0)solidporosity=0d0
      return
      end

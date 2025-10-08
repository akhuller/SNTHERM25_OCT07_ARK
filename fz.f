c***********************************************************************
c  Subroutine FZ computes snowdepth and nodal mid-point positions
c  relative to snow/soil interface. REJ changed on 1/7/24 to accommodate
c  case where the whole layer is snow
c***********************************************************************
      subroutine fz(z,dz,snowdepth,nsoil,n)
c %W% %G%
      implicit none
      include 'const'
c
c Called by MAIN
c
c arguments
c
c z : Nodal Cumulative distance from snow/soil interface to
c         centroid of node [m]
c dz :  Nodal thickness [m]
c snowdepth : snow depth [m]
c nsoil : Number of soil nodes
c n : Number of nodes in model
c
      integer nsoil,n
      double precision z(nd),dz(nd),snowdepth
c local
c
c i: looping index.
c
      integer i
      snowdepth=0d0
c      z(1)=dz(1)/2.  REJ 1/7/24
      do 10 i=1,n ! Replaced whole loop.  REJ 1/7/24
         if(i.ge. nsoil+1)snowdepth=snowdepth+dz(i)
         if(i .gt.1)then
          z(i)=z(i-1)+(dz(i-1)+dz(i))/2.
         else
          z(i)=dz(i)/2d0
         endif
 10   continue ! End replace loop REJ
      return
      end

c***********************************************************************
c  GETINPUT reads input data from file
c***********************************************************************
      subroutine getinput(ifluxout,isolarcalc,ircalc,islope,itracks,
     &   itm,ioutfiltrate,imetcalc,ngoodmin,itimezone,iy2,jday2,
cRJ 12/10/23& ihour2,pinv,bp,bext,albsnow,height,dtmin,dtsmin,dtmax,dtssmax,
     &   ihour2,pinv,bp,frvis,albsnow,height,dtmin,dtsmin,dtmax,dtssmax,
     &   dssallowed,errtallowd,emsnow,dlatt,elev,dlongt, azslope,dzmin,
     &   dzn,dznm,fnm,ssisnow,frh,bifallin,dmlimit,istboff,iqturb,eta0,
     &   exp1)

c RJ note 12/11/23. dfrvis,fvis and exp1 are similiar!! See getmet. 
c REJ  IMPORTANT input change.  Bext replaced with fvis 
      implicit none
      include 'const'
      include 'arrays'  
      
c**Arguments (no nomenclature yet)

      integer ifluxout,isolarcalc,ircalc,islope,itracks,itm
      integer ioutfiltrate,imetcalc,ngoodmin,itimezone,iy2
      integer jday2,ihour2,istboff,iqturb
c      double precision pinv,bp,bext,albsnow,height(3) RJ 12/10/23
      double precision pinv,bp,frvis,albsnow,height(3)! RJ 12/10/23
      double precision dtmin,dtsmin,dtmax,dtssmax,dssallowed,errtallowd
      double precision emsnow,dlatt,elev,ssisnow,frh(ld),eta0
      double precision dlongt,azslope,dmlimit
      double precision dzmin,dzn,dznm,bifallin
c  July 8 1996 - added experemental parameter  (jcm)
      double precision exp1
      character*160 fnm(nfiles)

c     Need to fill in definitions!!  Can see list in Documentation in meantime
c     itm : read measured temperature flag 1=yes 0=no
 
c** Local
      integer i

      read(90,*)ln,pinv,ifluxout,bp,isolarcalc,ircalc,islope,itracks  
     1,frvis,itm,ioutfiltrate,dtbase,imetcalc,albsnow,ssisnow,bifallin ! REJ
     2,dmlimit,istboff,iqturb,eta0

      if(itm.eq.1) open(89,file=fnm(3),status='old')

c     If met generation is optioned, radiation must be computed
      if(imetcalc .eq. 1)then
        isolarcalc=1
        ircalc=1
      endif

c  Measurement heights for airtemp, wsp and rh(or dewpoint)
      Read(90,*)height(1),height(2),height(3)
c  Layer parameters
      read(90,*)(nn(i),ltype(i),qtz(i),znaught(i),cd(i),rce(i),rch(i),
     1           ck(i),csk(i),frh(i),i=1,ln)
c  Convergence related input
      read(90,*)ngoodmin,dtmin,dtsmin,dtmax,dtssmax,dssallowed,
     & errtallowd
      if(dint(dtmax) .gt. dint(dtbase/2d0))dtmax=real(dint(dtbase/2d0))
      if(dint(dtssmax) .gt. dint(dtbase/2d0))dtssmax=real(dint(dtbase/
     &  2d0))
      if(dtsmin .gt. 10.001)stop '**Max allowable DTSMIN is 10 sec**'
c  now call look-up chart to get soil parameters for this soil type.
c  if ltype >= 90, then the subroutine reads in user-suppllied
c  parameters from layer.in.
      do 24 i=1,ln
         n=n+nn(i)
         if (ltype(i) .eq. 1.or. ltype(i) .eq. 90)then   !snow/firn
c     This is a snow layer or firn layer  !REJ_2025/01/31
            if(nn(i) .gt. 0)then
               nosnowcover=0
            else
               nosnowcover=1
            endif
            bd(i)=0.0
            djp(i)=0.0
            alb(i)=albsnow
            em(i)=emsnow
c  Changes below to add pure ice, blue ice, and sea ice REJ 1/5/2024
         elseif(ltype(i) .ge. 2 .and. ltype(i) .le. 4)then     !soil
c  soilchart read from main input file for material codes >=90
c  else gets info from BLOCK Data file
            call soilchart(ltype(i),i,90)
            nsoil= nsoil+nn(i)
c REJ removed firn if statement, which is now in the snow branch.!2025/03/10
         elseif(ltype(i) .eq. 91)then     !blue ice
            bd(i)=0.0
            djp(i)=0.0
            alb(i)=albsnow
            em(i)=emsnow
         elseif(ltype(i) .eq. 92)then     !pure ice
            bd(i)=0.0
            djp(i)=0.0
            alb(i)=albsnow
            em(i)=emsnow
         elseif(ltype(i) .eq. 95)then     !sea ice
            bd(i)=0.0
            djp(i)=0.0
            alb(i)=albsnow
            em(i)=emsnow
         elseif(ltype(i) .gt.99)then          !User-supplied
c           set up to read from the layer_in  file
            stop 'usersupplied characteristics not implemented'
         end if
c  REJ 1/5/2024   End changes
 24   continue

 
      if(n .ge.nd)stop'Array size (nd)is less than node number (n)+1'
     
c  insolation and slope parameters
      if(islope + isolarcalc .gt. 0)then
        read(90,*)dlatt,dlongt,elev,azslope,itimezone
        if(dlatt .le. 0)stop 'Slope opt. not implemented for S. Hemisph'
      endif
c  July 8 1996 - changed track parameters to experemental parameters  (jcm)
      if(itracks.eq.2) then
      read(90,*)exp1
      else
      exp1 = 999.0
      endif
c  initial node values
      do 666 i=1,n
         read(90,*,end=39)to(i),dzo(i),bwo(i),do(i)
         t(i) = to(i); dz(i) = dzo(i); bw(i) = bwo(i); d(i) = do(i) !2025_03/03
         concen(i) = 0d0  !REJ added 4/25/24
         if(dzo(i) .lt. dzmin)print *,
     &   'WARNING: Initial thickness of node ',i,
     &   'is less than prescribed minimum of ',sngl(dzmin)
         go to 666
 39      stop 'ERROR: End of data list found in layer.in file.'
666   continue
c      close(90)  !REJ 12/29/23  Moved to MAIN  
      if(dzo(n) .gt. dzn)write(80,*)'WARNING: Initial thickness of node'
     &,n,'exceeds prescribed maximum of ',sngl(dzn)
      if(dzo(n-1) .gt. dznm)write(80,*)'WARNING: Initial thickness of no
     &de',n-1,'exceeds prescribed maximum of ',sngl(dznm)
      return
      end

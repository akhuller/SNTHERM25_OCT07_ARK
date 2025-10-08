c***********************************************************************
c  SUBDIVIDE subdivides top or next to top CV if they exceed 
c  respectively dzn or dznm. 
c  REJ_2025/08/21 Also added lower CVs. But this may not be needed.
c***********************************************************************
      subroutine subdivide(dzn,dznm,a1snow,dzinc,floo,iunit,icalcstep)
     &    !REJ_2025/08/20. iunit=45. Currently 'summit_2007_precip'
      implicit none
      include 'const'
      include 'arrays'
      include 'l_arrays.txt'
c 
c Called from MAIN
c
c Calls DENSITY
c
c Arguments
c
c dzn : surface node maximum thickness[m] after snowfall has stopped (1.67cm)
c dznm : node under surface node maximum thickness[m] after snowfall (3.33cm)
c        has stopped
c a1snow : unfrozen water constant
c dzinc : maximum thickness for new snow nodes [m]
c
c April 4, 1995  double precision dzn,dznm,a1snow,dzinc
      double precision dzn,dznm,a1snow,dzinc,floo(nd)
c
c Local
c
c bold: Old BB value of node being subdivided (bold=bbo(ndiv)).
c dzold: Old thickness of node being subdivided (dzold=dzo(ndiv)).
c propor: mass/thickness proportion.
c
c Passed through common (nomenclature incomplete)
c
c dzo(nd): Old elemental thickness [m]
c i: nodal index.
c iar: array element index.
c m: array element.
c n: Number of nodes in model
c narray1: Vectors in common block array1 =30 (in ARRAYS)
c narray2: Vectors in common block array2 =29 (in ARRAYS). Should be 
c          increased if vectors are added.
c niarray: Vectors in common block iarray =8 (in ARRAYS). Should be
c          increased if vectors are added.
c ndiv: index denoting node to be subdivided
c nbot: temporary lowest CV checked for subdivide.  Currently = 30
      
c H2Olayer(nd) In l_arrays_txt

c 
      double precision propor,dzold,bold
      integer i,iar,m,ndiv,iunit,nbot,icalcstep !REJ_2025/08/20
      logical thickCV_found

cREJ_2025/01/12  stop 'Hit subdivide. H2Olayer not programmed for this routine'
cREJ_2025/08/06  Now programmed for H2O
      nbot = 30  ! Just for Summit run ! ARK_2025/09/22 Check with Rachel
      thickCV_found = .false.
      
      DO 100 i = n, max(n-nbot,nsoil-1), -1 !REJ_2025/08/20
       IF(i .eq. n)THEN
        if(dzo(n) .gt. dzn)then !REJ_2025
c  Top CV will be subdivided
          write(80,*)'Subdivided Top CV : ',n
          write(iunit,*)'Subdivided CV : ',n,'dzo(n),dzn',dzo(n),dzn
          ndiv=n
          goto 30
        endif
       ELSEIF(i .eq. n-1)THEN !REJ_2025
        if(dzo(n-1) .gt. dznm .or. dzo(n-1) .gt. .9*dzinc)then !REJ_2025
c  Next to top CV will be subdivided. Move top node up.
           write(80,*)'Subdivided CV : ',n-1
           write(iunit,*)'Subdiv n-1 CV :',n-1,'dzo(n-1),dzn',dzo(n-1),
     &                                        dznm
           ndiv=n-1
           m=n+1
           do 24 iar=1,narray1
              array1(m,iar)=array1(n,iar)
24         continue
           do 25 iar=1,narray2
              array2(m,iar)=array2(n,iar)
25         continue
           do 26 iar=1,niarray
              iarray(m,iar)=iarray(n,iar)
26         continue
           goto 30
         endif
c         u(m+1)=u(n+1)   !REJ_2025/08/20 Recheck these
c         uo(m+1)=uo(n+1)
c         layertype(m)=soiltype(ltype(k(n)))
c         write(*,*)'SUBD:m,layertype',m,layertype(m)
       ELSE ! Check remaining CVs REJ_2025/08/21
         if(dzo(i) .gt. 4d-2)then 
          ndiv = i
          write(80,*)'Subdivided CV : ',i
          write(iunit,*)'Subdivided CV :',i,'dzo(i),dzn',dzo(n),.04d0
          do 27 iar=1,narray1
            array1(i+1,iar)=array1(i,iar)
27        continue
          do 28 iar=1,narray2
            array2(i+1,iar)=array2(i,iar)
28        continue
          do 29 iar=1,niarray
            iarray(i+1,iar)=iarray(i,iar)
29        continue
         endif
         goto 30
       ENDIF
100   CONTINUE
         u(m+1)=u(n+1)
         uo(m+1)=uo(n+1)
c         layertype(m)=soiltype(ltype(k(n)))  ???

c         No thick node was found
          return
         
c
c  Now start subdivision process.  First set all new nodal array values
c  to that of node being subdivided.
c
30    m=ndiv+1
      uvapor(ndiv)=0d0
      do 34 iar=1,narray1
         array1(m,iar)=array1(ndiv,iar)
34    continue
      do 35 iar=1,narray2
         array2(m,iar)=array2(ndiv,iar)
35    continue
      do 36 iar=1,niarray
         iarray(m,iar)=iarray(ndiv,iar)
36    continue
      n=n+1
      nn(k(ndiv))=nn(k(ndiv))+1
      layertype(m)=soiltype(ltype(k(ndiv)))
c REJ 6/24/24      if(ltype(k(ndiv)).gt.1)nsoil=nsoil+1
      if(.not.H2Olayer(ndiv))nsoil=nsoil+1 !case where a soil layer is subdived
c      if(H2Olayer(ndiv) .eqv. .true.)
c     &   H2Olayer(ndiv+1) = .true. !REJ_2025/01/20  ESSENTIAL!!
      if(H2Olayer(ndiv))H2Olayer(ndiv+1) = .true. !REJ_2025/01/20  ESSENTIAL!!
c
c  If ndiv = n or n-1, assign 1/3 of mass/thickness to upper node and 2/3 
c  to lower node. Otherwise assign 1/2 to each. !REJ_2025/08/24
c
      if(ndiv .ge. n-1)then !REJ_2025/08/24 Entire if/then block
        propor=2d0/3d0
      else
        propor= 0.5d0
      endif
      dzold=dzo(ndiv)
      bold=bbo(ndiv)
      binew(ndiv+1)=binew(ndiv)
      do 10 i=ndiv,ndiv+1
         dzo(i)=dzold*propor
         bbo(i)=bold*propor
c        First melt=0 forces calculation of bl in density
         m=k(ndiv)
c 6/24/24         if(ltype(m).le. 1)then
         if(H2Olayer(i))then
         call density(to(i),bwo(i),td(i),0d0,flo(i),blo(i),bi(i),
c     1      bt(i),dmass(i),0d0,0d0,0d0,a1snow,0,dzo(i),!REJ_2025/08/22
     1    bt(i),dmasso(i),0d0,0d0,0d0,a1snow,0,dzo(i),!REJ_2025/08/22 
     2    flgo(i),ss(i),dice,ssi(i),porosity(i),
     3    dmvol(m), ltype(m),impermeable(i),idelete(i),
     4    solidporosity(m),ipond(i),dicevol(i),dliqvol(i),rhowater,
     5    H2Olayer(i)) ! ARK_2025/09/22
         else
         call density(to(i),bwo(i),td(i),td13(i),flo(i),blo(i),bi(i),
     1      bt(i),dmass(i),bdjp(m),a243(m),bd(m),a1(m),0,dzo(i),
     2      flgo(i),ss(i),dice,ssi(i),porosity(i),dmvol(m),
     3      ltype(m),impermeable(i),idelete(i),
     4      solidporosity(m),ipond(i),dicevol(i),dliqvol(i),rhowater,
     5    H2Olayer(i)) ! ARK_2025/09/22
         endif
         melt(i)=melt(ndiv)
c         propor=1d0/3d0 !REJ_2025/08/24
         propor = 1d0-propor  !REJ_2025/08/24
c Next added on April 4, 1995
         floo(i)=flo(i)
 10   continue
      write(iunit,*)'dzo(ndiv),dzo(ndiv+1)',ndiv,dzo(ndiv),
     &  ndiv+1,dzo(ndiv+1) !REJ_2025/08/24

      return
      end

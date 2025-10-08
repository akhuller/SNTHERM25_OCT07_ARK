c***********************************************************************
c  Subroutine FBB returns the sum of one-half the nodal net conductive
c  and convective energy fluxes and the nodal absorbed solar energy,
c  defined as bb.  This is done for efficiency purposes, to preclude
c  repeated re-summation of these terms, which appear of the
c  right-hand-side of the heat equation.
c  Because a Runge-Kutta temporal scheme is used, the total energy
c  flux is one half the past value, bbo, and one half the current value,
c  bb, which is based on updated temperatures computed in thermal. 
c  Fbb updates bb after thermal. It becomes bbo in the next step.  
c***********************************************************************
      subroutine fbb(qk,t,dsol,qf,topflux,botflux,bb,n,nsoil,
c     &      jday,ihour,Neumann,prnt_hour)  !10/30/23 RJ & RJ12/17
     &     jday,ihour,Neumann,prnt_hour,H2Olayer,heatfluxbtop,icalcstep,
     &     bbo,iunit) !10/30/23 RJ & RJ12/17;REJ_2025/02/09+08/17+09/04 
c     &   
c
c  Called by MAIN. Sec 15 (Initialize), Sec 16 (After thermal balance)
c  
c
c  arguments
c
c qk : .5*(thermal conductivity at upper nodal boundary/nodal thickness)
c Note the factor 0.5 that is included in qk and qf 
c t : Nodal temperature [K]
c dsol : Nodal absorbed solar radiation
c qf : .5*specific heat of water*nodal mass water flux 
c topflux: 0.5*Energy flux across air/media interface, excluding solar 
c          flux [W/m^2]
c botflux: Geothermal heat flux across bottom interface. Positive
c          downward [W/m^2]  ! RJ November 2023. 
c          -0.051 W/m2 Burton-Johnson for general Liston location
c to : Old nodal temperature [K]
c bb : .5*(Nodal conducted and convected heat plus absorbed solar)
c n : Number of nodes in model
c Neumann :  Flux bottom BC if .true.
c iter  :
c prnt_hour : Set in THPARAM.f and passed to FBB.f  !RJ17/23
      implicit none 
      include 'const'
      integer n,nsoil,ibasestep,jday,ihour,iunit !REJ_2025/09/04
      double precision qk(nd),t(nd),dsol(nd),qf(nd),bb(nd) !10/22/23  RJ
      double precision topflux,botflux  !10/22/23  RJ
      logical Neumann,prnt_hour,H2Olayer(nd)! RJ12/17

c  local
c
c heatfluxb : .5*Heatflux across lower boundary of node [W/m^2]
c heatfluxbtop : .5*Heatflux across lower boundary of top node [W/m^2] !REJ_2025/02/09
c heatfluxt : .5*Heatflux across upper boundary of node [W/m^2]
c heatfluxbsave  :  Saved for print-out
c i         : looping index.
c nsoil     : Number of soil nodes

      double precision heatfluxb,heatfluxt,heatfluxbtop !REJ_2025/02/09 
      double precision bbo(nd)
      integer i,iprnt,icalcstep !REJ_2025/08/17 +icalstep
      logical CallTestFbb, iprnthour

      CallTestFbb = .true. !generates a testing print-out for a specified
      ! jday and ihour.

      CallTestFbb = .false.; iprntHour = .false.  ! RJ12/17
      if(CallTestFbb .and.jday .eq. 160 .and. ihour .eq. 13) ! RJ12/17_REJ_2025/08/17
     &   iprnthour = .false.! RJ12/17 was prnt_hour  
      if(iprnthour)open(20,file='fbb_out',status='unknown')
      if(iprnthour)write(20,*)'fbb_out'
      if(iprnthour)write(20,*)'jday =',jday, '  ihour =',ihour ,
     & 'icalcstep =',icalcstep !RJ12/17 was prnt_hour,REJ_2025/08/17+icalcstep
      if(iprnthour)write(20,*)'Followng variabls, excpt for T, are *1/2'  ! RJ12/17
      if(iprnthour)write(20,*)  ! RJ12/17
     & '  node      bbo       hfluxt      hfluxb     dsol     qk      T'   

c  Top node
      heatfluxb = qk(n)*(t(n)-t(n-1))+qf(n)*t(n)
      heatfluxbtop = heatfluxb   !REJ_2025/02/09 
c      write(*,*)'heatfluxbtop',heatfluxbtop,'qk(n),t(n),t(n-1)',
c     &      qk(n),t(n),t(n-1),'topflux',topflux
      heatfluxt = topflux !REJ_2025/08/18
      bb(n) = topflux-heatfluxb +.5*dsol(n)
C     if(icalcstep .ge. 550)then !REJ_2025/09/04.  block added
C      write(iunit,*)
C    &  '3. FBB topflux  heatfluxbtop  0.5*dsol(n)   bb(n)   bbo(n)'
C      write(iunit,5)topflux,heatfluxbtop,.5d0*dsol(n),bb(n),bbo(n)         
C5      format(5x,5f11.4)
C     endif
      if(iprnthour)call TestFbb(n,bb,heatfluxt,heatfluxb,  ! RJ12/17 was prnt_hour
     &   topflux,botflux,dsol,qk,T,1)

c  Interior nodes
      do 20 i=n-1,2,-1
cREJ 6/24/24        if(i .gt. nsoil)then !Snow/ice
         if(H2Olayer(i))then  !Snow/Firn/Ice
          heatfluxt = heatfluxb
          heatfluxb = qk(i)*(t(i)-t(i-1))+qf(i)*t(i)
         elseif(i .eq. nsoil)then !Top Soil Node. No water flow for now
          heatfluxt = qk(i)*(t(i+1)-t(i))
          heatfluxb = qk(i)*(t(i)-t(i-1))
        else  ! Remaining soil nodes.  No water flow
          heatfluxt = heatfluxb
          heatfluxb = qk(i)*(t(i)-t(i-1))
        end if
        bb(i) = heatfluxt-heatfluxb + 0.5*dsol(i)
        if(iprnthour)call TestFbb(i,bb,heatfluxt,heatfluxb,  ! RJ12/17 was prnt_hour
     &   topflux,botflux,dsol,qk,T,2)
20    continue

c   Bottom node
      heatfluxt = heatfluxb
      if(.not. Neumann)then 
        heatfluxb=heatfluxt  !REJ3 1/21/24  !necessary dor bottom Dirichlet BC 
      else
        heatfluxb = 0.5d0*botflux
        heatfluxb = 0.5d0*botflux +qf(1)*T(1)
      endif
      bb(1) = heatfluxt-heatfluxb + 0.5*dsol(1)
      if(iprnthour)call TestFbb(1,bb,heatfluxt,heatfluxb,  ! RJ12/17
     &   topflux,botflux,dsol,qk,T,3)

      return
      end
c***********************************************************************
      Subroutine TestFbb(j,bb,heatfluxt,heatfluxb,topflux,
     &   botflux,dsol,qk,T,ipassin)
c***********************************************************************
      include 'const'
      integer j,ipassin
      double precision bb(nd),heatfluxt,heatfluxb,topflux,botflux,
     &  dsol(nd),qk(nd),T(nd)
5      format(i5,7f12.6) 
      if(ipassin .eq.1)then !Top node
        write(20,5)j,bb(j),
     &  Topflux,heatfluxb,.5d0*dsol(j),qk(j),T(j)
      elseif(ipassin .eq.3)then !Bottom node
        write(20,5)j,bb(j),
     &  heatfluxt,0.5d0*Botflux,0.5d0*dsol(j),qk(j),T(j)
c         stop 'End of fbb print-out'  !RJ12/17
      else  ! Interior nodes
        write(20,5)j,bb(j),
     &  heatfluxt,heatfluxb,.5d0*dsol(j),qk(j),T(j)
      endif

      end


c***********************************************************************
c Subroutine SDSOL computes absorption of solar radiation within snow
c cover.  Equation for extinction coefficient is eq. 4.26 from
c e. anderson, p.58, where the constant cv is 3.795d-3 (si units), p.80.
c This version was last used in November,2023
c***********************************************************************
c RJ note on 11/3/2023.  Lots of clean-ups to code.  Not flagged since
c there were many changes, the routine is condensed, but the
c results are the same.
      subroutine sdsol(dsol,d,dmass,fext,botnode,nsoil,n,solar,bextnir)
c
      include 'const'
c
c Called from MAIN
c
c Arguments
c
c dsol= shortwave radiation absorbed within element
c d=grain diameter
c dmass=bw*dz
c fext(nd) : Nodal transmission coefficient for the solar flux, 
c       exp(-cv*dmass/sqt(d)) or exp(-cv*dmass/sqt(d))*exp(-bext*.002)
c botnode : nodal solar penetration limit. 
c solar : Net solar radiative flux [W/m^2]
c bextnir : Snow extinction coefficient for near-IR
c nsoil : Number of soil nodes
c n : Number of nodes in model

      integer botnode,nsoil,n
      double precision dsol(nd),d(nd),dmass(nd),fext(nd)
      double precision solar,bextnir
c
c Local
c
c cv: empirical constant (3.795d-3 (si units))
c j: array index (Top down) 
c botdsol  : solar absorption is taken as 0d0 beyond this minimum 
c solflux(nd) Descending solar flux in snow or ice

      integer j
      double precision cv,solflux(nd),botdsol
c
c 10/26  RJ replaced "depth" with lowest node, "botnode" 
      cv=3.795d-3
      botdsol =1.0d-6
c
c Initialize dsol and solflux to zero 
c
      do 5 j=n,nsoil+1,-1
        dsol(j) = 0d0
        solflux(j) = 0d0
5     continue

c  Next loop computes solar flux at bottom of descending nodes. Dsol(j) is
c  then computed as solflux(j)-solflux(j-1). The flux at the surface, solfux(n), 
c  is the net solar from met file, solar.
 
      solflux(n) = solar
      botnode = n + 1     
      do 10 j=n,nsoil+1,-1
c  A minimum dsol for solar penetration (botdsol) is assumed. dsol is not
c  computed below this. 
         IF(botnode .eq. n + 1)THEN
          if(d(j).le.0.0) then
            stop '** ERROR SDSOL: Grain Diameter <= 0.0. Execution Stopp
     &ed**'
          end if
c         Note: cv*dmass(j)/dsqrt(d(j)) = bext*dz = fext  
          fext(j)=dexp(-cv*dmass(j)/dsqrt(d(j)))
c          write(*,*)j,fext(j),dexp(-bextnir*2d-3)           
          if(j .eq. n)then
             fext(n)=dexp(-bextnir*2d-3)*fext(n)
          end if
c          write(*,*)j,'fext(j)',fext(j)
          solflux(j-1)=solflux(j)*fext(j)
          dsol(j)=solflux(j)-solflux(j-1) 
c          write(*,*)solflux(j),solflux(j-1),'dsol',dsol(j)
c          if(j .eq. 33)stop
          if(dsol(j) .lt. botdsol)botnode = j                
        ENDIF
 10   continue

c  If solar flux penetrates to soil/rock base, remaining solflux is 
c  absorbed in top node of the base. Also need to consider case where 
c  base is ice.
c
c      if(botnode .eq. nsoil)dsol(nsoil) = solflux(nsoil)!Need to test this.
c            
c       write(*,*)'solar',solar,botnode,botdsol
c        write(*,*)'dsol', (dsol(j),j=botnode,n)
c        pause
      return
      end

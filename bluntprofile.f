c-----------------------------------------------------------------------
      subroutine bluntprofile_setup(wdin,trptin)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      real wdin(lx1,ly1,lz1,*)

      integer iel,ifc,i,j,k,i0,i1,j0,j1,k0,k1,ia
      real fact,aa,term,glsum,glmax
      real factor,wdmax,trpt,trptin
      common /bvel/ factor,wdmax,trpt

      character*3 cb3
      character cb1(3)

      equivalence (cb3,cb1)

      wdmax=-1.0
      do 15 iel = 1,nelv
      do 15 ifc = 1,2*ldim
        cb3=cbc(ifc,iel,1)
        if(cb1(1).eq.'v') then
          call facind(i0,i1,j0,j1,k0,k1,lx1,ly1,lz1,ifc)
          do 25 k=k0,k1
          do 25 j=j0,j1
          do 25 i=i0,i1
            wdmax=max(wdmax,wdin(i,j,k,iel))
  25      continue
        endif
  15  continue
      wdmax=glmax(wdmax,1)

      trpt = trptin

      factor = 0.0
      aa = 0.0
      do 10 iel = 1,nelv
      do 10 ifc = 1,2*ldim
        cb3=cbc(ifc,iel,1)
        if(cb1(1).eq.'v')then
          call facind(i0,i1,j0,j1,k0,k1,lx1,ly1,lz1,ifc)
          ia = 0
          do 20 k=k0,k1
          do 20 j=j0,j1
          do 20 i=i0,i1
            ia = ia + 1
            term=wdin(i,j,k,iel)/wdmax
            if(term.lt.trpt) then
              fact=1.0-(1.0-term/trpt)**2
            else
              fact=1.0
            endif
            factor=factor+area(ia,1,ifc,iel)*fact
            aa=aa+area(ia,1,ifc,iel)
 20       continue
        endif
 10   continue
      factor=glsum(factor,1)
      aa=glsum(aa,1)
      factor=aa/factor

      return
      end
c-----------------------------------------------------------------------
      real function blunt_vel(vel,dist)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      real vel,dist

      real factor,wdmax,trpt,term
      common /bvel/ factor,wdmax,trpt

      term=dist/wdmax
      if(term.lt.trpt) then
        blunt_vel=vel*(1.0-(1.0-term/trpt)**2)*factor
      else
        blunt_vel=vel*factor
      endif

      return
      end
c-----------------------------------------------------------------------

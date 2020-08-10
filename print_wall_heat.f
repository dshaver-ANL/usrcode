      subroutine print_wall_heat
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer iel,ifc,i0,i1,i,j0,j1,j,k0,k1,k,ipt
      logical ifdograd
 
      real glsum
      real dtdx(lx1,ly1,lz1)
      real dtdz(lx1,ly1,lz1)
      real dtdy(lx1,ly1,lz1)
      real qflux,qtot,a1
      real norm(3)

      qtot = 0.0
      do iel=1,nelv
      ifdograd=.true.
        do ifc=1,2*ndim
        if(cbc(ifc,iel,1).eq.'W  ') then
          if(ifdograd) then !only evaluate gradient once per element
            call gradm1(dtdx,dtdy,dtdz,t,iel)
            ifdograd=.false.
          endif
          ipt=1
          do 20 k=k0,k1
          do 20 j=j0,j1
          do 20 i=i0,i1
            call getSnormal(norm,i,j,k,ifc,iel)
            qflux=(dtdx(i,j,k)*norm(1)+
     &             dtdy(i,j,k)*norm(2)+
     &             dtdz(i,j,k)*norm(3))*vdiff(i,j,k,iel)
            qtot=qtot+qflux*area(ipt,1,ifc,iel)
            ipt=ipt+1
 20       continue
        endif
        enddo
      enddo
      qtot=glsum(qtot,1)      

      if(nio.eq.0) write(*,*) 'total heat through walls = ',qtot    

      return
      end

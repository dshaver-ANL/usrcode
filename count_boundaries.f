C-----------------------------------------------------------------------
      subroutine count_boundaries
      include 'SIZE'
      include 'TOTAL'

      integer lxyz,ielem,iside,n
      parameter(lxyz=lx1*ly1*lz1)
      character*3 uid(ldimt1)
      integer wcnt(ldimt1),symcnt(ldimt1),ocnt(ldimt1)
      integer tcnt(ldimt1),fcnt(ldimt1),axicnt(ldimt1)
      integer inscnt(ldimt1),pcnt(ldimt1),othcnt(ldimt1)
      integer vcnt(ldimt1),vlcnt(ldimt1),trcnt(ldimt1),ukncnt(ldimt1)
      integer mtrcnt(ldimt1),prcnt(ldimt1),intcnt(ldimt1)
      integer vreacnt(ldimt1),treacnt(ldimt1),convcnt(ldimt1)

      call izero(wcnt,ldimt1)
      call izero(trcnt,ldimt1)
      call izero(mtrcnt,ldimt1)
      call izero(vcnt,ldimt1)
      call izero(vlcnt,ldimt1)
      call izero(vreacnt,ldimt1)
      call izero(symcnt,ldimt1)
      call izero(ocnt,ldimt1)
      call izero(tcnt,ldimt1)
      call izero(treacnt,ldimt1)
      call izero(axicnt,ldimt1)
      call izero(inscnt,ldimt1)
      call izero(intcnt,ldimt1)
      call izero(fcnt,ldimt1)
      call izero(pcnt,ldimt1)
      call izero(prcnt,ldimt1)
      call izero(convcnt,ldimt1)
      call izero(othcnt,ldimt1)
      call izero(ukncnt,ldimt1)

      do ifld=1,nfield
        n=nelv
        if(iftmsh(ifld)) n=nelt
        do ielem=1,n
        do iside=1,2*ldim
          if(cbc(iside,ielem,ifld).eq.'W  ')then
            wcnt(ifld)=wcnt(ifld)+1
          elseif(cbc(iside,ielem,ifld).eq.'shl')then
            trcnt(ifld)=trcnt(ifld)+1
          elseif(cbc(iside,ielem,ifld).eq.'sml')then
            mtrcnt(ifld)=mtrcnt(ifld)+1
          elseif(cbc(iside,ielem,ifld).eq.'v  ')then
            vcnt(ifld)=vcnt(ifld)+1
          elseif(cbc(iside,ielem,ifld).eq.'vl ')then
            vlcnt(ifld)=vlcnt(ifld)+1
          elseif(cbc(iside,ielem,ifld).eq.'V  ')then
            vreacnt(ifld)=vreacnt(ifld)+1
          elseif(cbc(iside,ielem,ifld).eq.'t  ')then
            tcnt(ifld)=tcnt(ifld)+1
          elseif(cbc(iside,ielem,ifld).eq.'T  ')then
            treacnt(ifld)=treacnt(ifld)+1
          elseif(cbc(iside,ielem,ifld).eq.'O  ')then
            ocnt(ifld)=ocnt(ifld)+1
          elseif(cbc(iside,ielem,ifld).eq.'o  ')then
            prcnt(ifld)=prcnt(ifld)+1
          elseif(cbc(iside,ielem,ifld).eq.'P  ')then
            pcnt(ifld)=pcnt(ifld)+1
          elseif(cbc(iside,ielem,ifld).eq.'f  ')then 
            fcnt(ifld)=fcnt(ifld)+1
          elseif(cbc(iside,ielem,ifld).eq.'c  ')then 
            convcnt(ifld)=convcnt(ifld)+1
          elseif(cbc(iside,ielem,ifld).eq.'I  ')then
            inscnt(ifld)=inscnt(ifld)+1
          elseif(cbc(iside,ielem,ifld).eq.'SYM')then
            symcnt(ifld)=symcnt(ifld)+1
          elseif(cbc(iside,ielem,ifld).eq.'A  ')then
            axicnt(ifld)=axicnt(ifld)+1
          elseif(cbc(iside,ielem,ifld).eq.'int')then
            intcnt(ifld)=intcnt(ifld)+1
          elseif(cbc(iside,ielem,ifld).eq.'   ')then
            othcnt(ifld)=othcnt(ifld)+1
          elseif(cbc(iside,ielem,ifld).ne.'E  ')then
c           if(ukncnt(ifld).eq.0) then  !handle one unknown BC
c             uid(ifld)=cbc(iside,ielem,ifld) !can't synch uid across domains
c             ukncnt(ifld)=1
c           elseif(cbc(iside,ielem,ifld).eq.uid(ifld)) then
              ukncnt(ifld)=ukncnt(ifld)+1
c           endif
          endif
        enddo
        enddo
        wcnt(ifld)=iglsum(wcnt(ifld),1)
        trcnt(ifld)=iglsum(trcnt(ifld),1)
        mtrcnt(ifld)=iglsum(mtrcnt(ifld),1)
        vcnt(ifld)=iglsum(vcnt(ifld),1)
        vlcnt(ifld)=iglsum(vlcnt(ifld),1)
        vreacnt(ifld)=iglsum(vreacnt(ifld),1)
        tcnt(ifld)=iglsum(tcnt(ifld),1)
        treacnt(ifld)=iglsum(treacnt(ifld),1)
        ocnt(ifld)=iglsum(ocnt(ifld),1)
        prcnt(ifld)=iglsum(prcnt(ifld),1)
        pcnt(ifld)=iglsum(pcnt(ifld),1)
        fcnt(ifld)=iglsum(fcnt(ifld),1)
        convcnt(ifld)=iglsum(convcnt(ifld),1)
        inscnt(ifld)=iglsum(inscnt(ifld),1)
        symcnt(ifld)=iglsum(symcnt(ifld),1)
        axicnt(ifld)=iglsum(axicnt(ifld),1)
        intcnt(ifld)=iglsum(intcnt(ifld),1)
        othcnt(ifld)=iglsum(othcnt(ifld),1)
        ukncnt(ifld)=iglsum(ukncnt(ifld),1)
      enddo

      if(nid.eq.0) then
        write(*,*)
        write(*,255) 'Found the following Boundary Conditions'
        write(*,*)
        do ifld=1,nfield
          write(*,254) 'for field',ifld,':'
          if(wcnt(ifld).gt.0)write(*,256)'Wall',wcnt(ifld)
          if(trcnt(ifld).gt.0)write(*,256)'Traction',trcnt(ifld)
          if(mtrcnt(ifld).gt.0)write(*,256)'Mixed Traction'
     &                                                   ,mtrcnt(ifld)
          if(vcnt(ifld).gt.0)write(*,256)'Velocity',vcnt(ifld)
          if(vlcnt(ifld).gt.0)write(*,256)'Velocity (local)',vlcnt(ifld)
          if(vreacnt(ifld).gt.0)write(*,256)'Velocity (REA)',
     &                                                   vreacnt(ifld)
          if(tcnt(ifld).gt.0)write(*,256)'Dirichlet',tcnt(ifld)
          if(treacnt(ifld).gt.0)write(*,256)'Dirichlet (REA)',
     &                                                   treacnt(ifld)
          if(pcnt(ifld).gt.0)write(*,256)'Periodic',pcnt(ifld)
          if(fcnt(ifld).gt.0)write(*,256)'Flux',fcnt(ifld)
          if(convcnt(ifld).gt.0)write(*,256)'Convection',convcnt(ifld)
          if(ocnt(ifld).gt.0)write(*,256)'Outlet',ocnt(ifld)
          if(prcnt(ifld).gt.0)write(*,256)'Pressure',prcnt(ifld)
          if(inscnt(ifld).gt.0)write(*,256)'Insulated',inscnt(ifld)
          if(intcnt(ifld).gt.0)write(*,256)'Interpolated',intcnt(ifld)
          if(symcnt(ifld).gt.0)write(*,256)'Symmetry',symcnt(ifld)
          if(axicnt(ifld).gt.0)write(*,256)
     &                                     'Axisymmetric',axicnt(ifld)
          if(othcnt(ifld).gt.0)write(*,256)'"   "',othcnt(ifld)
          if(ukncnt(ifld).gt.0)write(*,256)'Unknown',ukncnt(ifld)
          write(*,*)
        enddo
      endif

 254  format(5x,a,i2,a)
 255  format(2x,a)
 256  format(2x,a16,1x,i12)
c257  format(2x,'Unknown boundary of type: "',a,'" ',i9)

      return
      end
C-----------------------------------------------------------------------

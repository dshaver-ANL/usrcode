      subroutine load_inlet(fname)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      character*(*) fname

      integer i,n

      real uin,vin,win,tin
      common /INLBCs/ uin(lx1*ly1*lz1*lelv)
     &               ,vin(lx1*ly1*lz1*lelv)
     &               ,win(lx1*ly1*lz1*lelv)
     &               ,tin(lx1*ly1*lz1*lelv,ldimt)

      n=lx1*ly1*lz1*nelv

      if(nio.eq.0) write(*,*)
     &                     'loading inlet data from file ','"',fname,'"'

      call store_solution

      call gfldr(fname)

      call copy(uin,vx,n)
      call copy(vin,vy,n)
      if(if3d) call copy(win,vz,n)
      do i=1,ldimt
        call copy(tin(1,i),t(1,1,1,1,i),n)
      enddo

      call reload_solution

      return
      end
c-----------------------------------------------------------------------
      subroutine load_inlet_region(fname,dir,xx0,xx1)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      character*(*) fname
      character dc(3)

      integer i,j,k,n1,n2,dir
      real eps,fact,x0,xx0,x1,xx1,xx

      real uin,vin,win,tin,prin
      common /INLBCs/ uin(lx1*ly1*lz1,lelv)
     &               ,vin(lx1*ly1*lz1,lelv)
     &               ,win(lx1*ly1*lz1,lelv)
     &               ,tin(lx1*ly1*lz1,lelv,ldimt)
     &               ,prin(lx2*ly2*lz2,lelv)

      data dc /'X','Y','Z'/

      n1=lx1*ly1*lz1
      n2=lx2*ly2*lz2
      x0=xx0
      x1=xx1

      if(nio.eq.0) write(*,'(a,f5.3,3a,f5.3,4a)') 'loading data for '
     &              ,x0,'<=',dc(dir),'<=',x1,' from file ','"',fname,'"'

      call store_solution

      call gfldr_rescale(fname,dir,x0,x1)

      call copy(uin,vx,n1*nelv)
      call copy(vin,vy,n1*nelv)
      if(if3d) call copy(win,vz,n1*nelv)
      call copy(prin,pr,n2*nelv)
      do i=1,ldimt
        call copy(tin(1,1,i),t(1,1,1,1,i),n1*nelv)
      enddo

      call reload_solution

      eps=0.005*(x1-x0)
      x1=x0+0.975*(x1-x0)
      do j=1,nelv
        do i=1,n1
          if(dir.eq.1)xx=xm1(i,1,1,j)
          if(dir.eq.2)xx=ym1(i,1,1,j)
          if(dir.eq.3)xx=zm1(i,1,1,j)
          fact=0.5*(1.0+tanh((xx-x1)/eps))
          vx(i,1,1,j)=uin(i,j)*(1.0-fact)+vx(i,1,1,j)*fact
          vy(i,1,1,j)=vin(i,j)*(1.0-fact)+vy(i,1,1,j)*fact
          vz(i,1,1,j)=win(i,j)*(1.0-fact)+vz(i,1,1,j)*fact
          do k=1,ldimt
            t(i,1,1,j,k)=tin(i,j,k)*(1.0-fact)+t(i,1,1,j,k)*fact
          enddo
        enddo
        do i=1,n2
          if(dir.eq.1)xx=xm2(i,1,1,j)
          if(dir.eq.2)xx=ym2(i,1,1,j)
          if(dir.eq.3)xx=zm2(i,1,1,j)
          fact=0.5*(1.0+tanh((xx-x1)/eps))
          pr(i,1,1,j)=prin(i,j)*(1.0-fact)+pr(i,1,1,j)*fact
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine gfldr_rescale(sourcefld,dir,x0,x1)
c
c     generic field file reader
c     reads sourcefld and interpolates all avaiable fields
c     onto current mesh
c
c     memory requirement: 
c     nelgs*nxs**ldim < np*(4*lelt*lx1**ldim)
c
      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'
      include 'GFLDR'

      character sourcefld*(*)

      common /scrcg/  pm1(lx1*ly1*lz1,lelv)
      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      character*1   hdr(iHeaderSize)

      integer*8 dtmp8
      integer dir

      logical if_byte_swap_test
      real*4 bytetest

      etime_t = dnekclock_sync()
      if(nio.eq.0) write(6,*) 'call gfldr ',trim(sourcefld) 

      ! open source field file
      ierr = 0
      if(nid.eq.0) then
        open (90,file=sourcefld,status='old',err=100)
        close(90)
        goto 101
 100    ierr = 1
 101  endif
      call err_chk(ierr,' Cannot open source fld file!$')
      call byte_open_mpi(sourcefld,fldh_gfldr,.true.,ierr)

      ! read and parse header
      call byte_read_mpi(hdr,iHeaderSize/4,0,fldh_gfldr,ierr)
      call byte_read_mpi(bytetest,1,0,fldh_gfldr,ierr)

      call mfi_parse_hdr(hdr,ierr)
      call err_chk(ierr,' Invalid header!$')
      ifbswp = if_byte_swap_test(bytetest,ierr)
      call err_chk(ierr,' Invalid endian tag!$')

      nelgs   = nelgr
      nxs     = nxr
      nys     = nyr
      nzs     = nzr
      if(nzs.gt.1) then 
        ldims = 3
      else
        ldims = 2
      endif

      ! distribute elements across all ranks
      nels = nelgs/np
      do i = 0,mod(nelgs,np)-1
         if(i.eq.nid) nels = nels + 1
      enddo
      nxyzs      = nxs*nys*nzs
      dtmp8      = nels
      ntots_b    = dtmp8*nxyzs*wdsizr
      rankoff_b  = igl_running_sum(nels) - dtmp8
      rankoff_b  = rankoff_b*nxyzs*wdsizr  
      dtmp8      = nelgs
      nSizeFld_b = dtmp8*nxyzs*wdsizr
      noff0_b    = iHeaderSize + iSize + iSize*dtmp8

      ! do some checks
      if(ldims.ne.ldim) 
     $ call exitti('ldim of source does not match target!$',0)
      if(ntots_b/wdsize .gt. ltots) then
        dtmp8 = nelgs
        lelt_req = dtmp8*nxs*nys*nzs / (np*ltots/lelt)
        lelt_req = lelt_req + 1
        if(nio.eq.0) write(6,*)
     $   'ABORT: buffer too small, increase lelt > ', lelt_req
        call exitt
      endif

      ifldpos = 0
      if(ifgetxr) then
        ! read source mesh coordinates
        call gfldr_getxyz(xm1s,ym1s,zm1s)
        ifldpos = ldim
      else
        call exitti('source does not contain a mesh!$',0)
      endif

      ns=nels*nxyzs
      if(nio.eq.0) write(*,*) dir,ns,ltots
      if(dir.eq.1) then
        xmin=glmin(xm1s,nels*nxyzs)
        xmax=glmax(xm1s,nels*nxyzs)
        scale=(x1-x0)/(xmax-xmin)
        do i=1,nels*nxyzs
          xm1s(i)=x0+scale*(xm1s(i)-xmin)
        enddo
      elseif(dir.eq.2) then
        xmin=glmin(ym1s,nels*nxyzs)
        xmax=glmax(ym1s,nels*nxyzs)
        scale=(x1-x0)/(xmax-xmin)
        do i=1,nels*nxyzs
          ym1s(i)=x0+scale*(ym1s(i)-xmin)
        enddo
      elseif(dir.eq.3) then
        xmin=glmin(zm1s,nels*nxyzs)
        xmax=glmax(zm1s,nels*nxyzs)
        scale=(x1-x0)/(xmax-xmin)
        do i=1,nels*nxyzs
          zm1s(i)=x0+scale*(zm1s(i)-xmin)
        enddo
      endif
      if(nio.eq.0) write(*,*) 2,dir,ns,ltots

      if(if_full_pres) then
        call exitti('no support for if_full_pres!$',0)
      endif

      ! initialize interpolation tool using source mesh
      nxf   = 2*nxs
      nyf   = 2*nys
      nzf   = 2*nzs
      nhash = nels*nxs*nys*nzs 
      nmax  = 128

      call fgslib_findpts_setup(inth_gfldr,nekcomm,np,ldim,
     &                          xm1s,ym1s,zm1s,nxs,nys,nzs,
     &                          nels,nxf,nyf,nzf,bb_t,
     &                          nhash,nhash,nmax,tol)

      ! read source fields and interpolate
      if(ifgetur) then
        if(nid.eq.0 .and. loglevel.gt.2) write(6,*) 'reading vel'
        ntot = nx1*ny1*nz1*nelv 
        call gfldr_getfld(vx,vy,vz,ntot,ldim,ifldpos+1)
        ifldpos = ifldpos + ldim
      endif
      if(ifgetpr) then
        if(nid.eq.0 .and. loglevel.gt.2) write(6,*) 'reading pr'
        ntot = nx1*ny1*nz1*nelv 
        call gfldr_getfld(pm1,dum,dum,ntot,1,ifldpos+1)
        ifldpos = ifldpos + 1
c       if (ifaxis) call axis_interp_ic(pm1)
c       call map_pm1_to_pr(pm1,1)
      endif
      if(ifgettr .and. ifheat) then
        if(nid.eq.0 .and. loglevel.gt.2) write(6,*) 'reading temp'
        ntot = nx1*ny1*nz1*nelfld(2) 
        call gfldr_getfld(t(1,1,1,1,1),dum,dum,ntot,1,ifldpos+1)
        ifldpos = ifldpos + 1
      endif
      do i = 1,ldimt-1
         if(ifgtpsr(i)) then
           if(nid.eq.0 .and. loglevel.gt.2) 
     $       write(6,*) 'reading scalar',i
           ntot = nx1*ny1*nz1*nelfld(i+2) 
           call gfldr_getfld(t(1,1,1,1,i+1),dum,dum,ntot,1,ifldpos+1) 
           ifldpos = ifldpos + 1
         endif
      enddo

      call byte_close_mpi(fldh_gfldr,ierr)
      etime_t = dnekclock_sync() - etime_t
      call fgslib_findpts_free(inth_gfldr)
      if(nio.eq.0) write(6,'(A,1(1g9.2),A)')
     &                   ' done :: gfldr  ', etime_t, ' sec'

      return
      end
c-----------------------------------------------------------------------

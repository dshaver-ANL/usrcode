c-----------------------------------------------------------------------
      subroutine lineplot(pt1,pt2,lpts)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer nptmax
      parameter(nptmax=1001)

      real pt1(ldim),pt2(ldim)
      integer npts,lpts,iplot

      character*32 fname
      character*14 afmt
      character*10 rfmt
      integer intp_h,i,j,nt,nfld
      save intp_h
      logical ifset,ifdo
      real dx,pts(nptmax,ldim)
      real fwrk(lx1*ly1*lz1*lelt,ldim+1+ldimt)
      real fpts(nptmax*(ldim+1+ldimt))
      real uout(nptmax),vout(nptmax),wout(nptmax)
      real prout(nptmax),tout(nptmax,ldimt)
      character*4 outname(ldim+1+ldimt)

      real rwrk(nptmax,ldim+1)
      integer iwrk(nptmax,3)
      save rwrk,iwrk

      save ifdo,ifset
      data ifdo /.true./
      data ifset /.true./

      save iplot
      data iplot /1/

      if(.not.ifdo) return

      nt=lx1*ly1*lz1*nelt

      npts=max(lpts,2)
      if(npts.gt.nptmax) then
        if(nio.eq.0) write(*,*)
     &               "Error in lineplot, recompile with nptmax >= ",npts
        ifdo=.false.
        return
      endif

      call rzero(pts,npts*ndim)
      do j=1,ndim 
        pts(1,j)=pt1(j)
        dx=(pt2(j)-pt1(j))/(real(npts-1))
        do i=2,npts
          pts(i,j)=pts(i-1,j)+dx
        enddo
      enddo

      if(ifset)then
        ifset=.false.
        call interp_setup(intp_h,0.0,0,nelt)
      endif

      nfld=0
      if(ifvo) then
        write(outname(1),'(a4)')"VELX"
        write(outname(2),'(a4)')"VELY"
        call copy(fwrk(1,1),vx,nt)
        call copy(fwrk(1,2),vy,nt)
        nfld=2
      endif
      if(if3d.and.ifvo)then
        nfld=nfld+1
        write(outname(nfld),'(a4)')"VELZ"
        call copy(fwrk(1,nfld),vz,nt)
      endif
      if(ifpo) then
        nfld=nfld+1
        write(outname(nfld),'(a4)')"PRES"
        call copy(fwrk(1,nfld),pr,nt)
      endif
      if(ifheat) then
        if(ifto) then
          nfld=nfld+1
          write(outname(nfld),'(a4)')"TEMP"
          call copy(fwrk(1,nfld),t,nt)
        endif
        do i=1,ldimt-1
          if(ifpsco(i)) then
            nfld=nfld+1
            write(outname(nfld),'(a2,i2)')"PS",i
            call copy(fwrk(1,nfld),t(1,1,1,1,i+1),nt)
          endif
        enddo
      endif

      if(nfld.gt.0) then
        call blank(fname,32)
        if(iplot.lt.10) then
          write(fname,'(a,i1,a)') "plot",iplot,".dat"
        elseif(iplot.lt.100) then
          write(fname,'(a,i2,a)') "plot",iplot,".dat"
        else
          write(fname,'(a,i3,a)') "plot",iplot,".dat"
        endif

        if(nio.eq.0) then
          write(*,*)'   Writing line plot data to file ',fname
          if(if3d)then
            write(*,'(7x,3es15.6)')pt1(1),pt1(2),pt1(3)
            write(*,'(7x,3es15.6)')pt2(1),pt2(2),pt2(3)
          else
            write(*,'(7x,2es15.6)')pt1(1),pt1(2)
            write(*,'(7x,2es15.6)')pt2(1),pt2(2)
          endif
          write(*,*)
        endif
   
        call interp_nfld(fpts,fwrk,nfld,pts(1,1),pts(1,2),pts(1,3),npts
     &                                  ,iwrk,rwrk,nptmax,.true.,intp_h)

        call blank(afmt,14)
        call blank(rfmt,10)
        if(if3d) then
          write(afmt,'(a1,i2,a11)')"(",nfld+3,"a16,es16.8)"
          write(rfmt,'(a1,i2,a7)')"(",nfld+3,"es16.8)"
        else
          write(afmt,'(a1,i2,a11)')"(",nfld+2,"a16,es16.8)"
          write(rfmt,'(a1,i2,a7)')"(",nfld+2,"es16.8)"
        endif
  
        if(nio.eq.0) then
          open(unit=10,file=fname,status='unknown',form='formatted')
          if(if3d) then
            write(10,afmt)"X","Y","Z",(outname(i),i=1,nfld),time
          else
            write(10,afmt)"X","Y",(outname(i),i=1,nfld),time
          endif
          do i=1,npts
            if(if3d) then
              write(10,rfmt)pts(i,1),pts(i,2),pts(i,3)
     &                               ,(fpts(i+j),j=0,(npts*nfld-1),npts)
            else
              write(10,rfmt)pts(i,1),pts(i,2)
     &                               ,(fpts(i+j),j=0,(npts*nfld-1),npts)
            endif
          enddo
        endif
  
        close(10)
  
        iplot=iplot+1
      endif

      return
      end
c-----------------------------------------------------------------------

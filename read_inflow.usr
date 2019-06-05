      subroutine read_inflow(uin,tin)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer lxyz,nmax,iglsum
      parameter(lxyz=lx1*ly1*lz1)
      parameter(nmax=max(lelx*max(lely,lelz),lely*lelz)*lx1*lx1)

      real uin(1),tin(lx1*ly1*lz1*lelt,1)

      character sdir
      integer i,j,npts,aflag,iset,ntot,npsf
      real scrd(nmax),tcrd(nmax),uf(nmax),tf(nmax,ldimt),wk(nmax)

      aflag=0
      npts=0
      iset=0
      ntot=lxyz*nelv
 
      npsf=0
      if(ifheat) then
        do i=1,ldimt
          if(idpss(i).ne.-1) npsf=npsf+1
        enddo
      endif

      call rzero(scrd,nmax)
      call rzero(tcrd,nmax)
      call rzero(uf,nmax)
      call rzero(tf,nmax*ldimt)

      if(nio.eq.0) then
        if(npsf.eq.0) then
          write(*,'(5x,2a)')"Setting inflow boundary and "
     &          ,"initial conditions for velocity from inflow.dat"
        else
          write(*,'(5x,2a,i1,1a)')"Setting inflow boundary and initial "
     &    ,"conditions for velocty and ",npsf," scalars from inflow.dat"
        endif
      endif

      if(nid.eq.0) then
        open(unit=25,file='inflow.dat',status='old'
     &                                         ,form='formatted',err=25)

        read(25,*) npts,sdir
 26     if(npts.gt.nmax) then !(26) returning from open file error
          write(*,'(7x,a)') "Warning: Too many points in inflow file"
          aflag=1
          write(*,'(9x,a)')"inflow data NOT set!"
          write(*,*)
        else
          write(*,'(7x,a,i5,a)')'reading ',npts,' points'
          write(*,*)
          do i=1,npts
            if(if3d) then
              read(25,*) scrd(i),tcrd(i),uf(i),(tf(i,j),j=1,npsf)
              if(nid.eq.0.and.npts.le.300) write(*,'(i3,25es19.12)')
     &                        i,scrd(i),tcrd(i),uf(i),(tf(i,j),j=1,npsf)
            else
              read(25,*) scrd(i),uf(i),(tf(i,j),j=1,npsf)
              if(nid.eq.0.and.npts.le.300) write(*,'(i3,25es19.12)')
     &                        i,scrd(i),uf(i),(tf(i,j),j=1,npsf)
            endif
          enddo
        endif
        close(25)
      endif

      aflag=iglsum(aflag,1)
      if(aflag.gt.0) return

      npts=iglsum(npts,1)
      sdir='Y'

      call gop(scrd,wk,'+  ',npts)
      if(if3d) call gop(tcrd,wk,'+  ',npts)
      call gop(uf,wk,'+  ',npts)
      do i=1,npsf
        call gop(tf(1,i),wk,'+  ',npts)
      enddo

      if(if3d) then
        if(sdir.eq.'Y')then !uf=ux, compare y,z coordinates
          call set_bc_arrays(ym1,zm1,uin,tin,ntot
     &                                  ,scrd,tcrd,uf,tf,npts,npsf,iset)
        elseif(sdir.eq.'Z')then !uf=uy, compare z,x coordinates
          call set_bc_arrays(zm1,xm1,uin,tin,ntot
     &                                  ,scrd,tcrd,uf,tf,npts,npsf,iset)
        else !uf=uz, compare x,y coordinates
          call set_bc_arrays(xm1,ym1,uin,tin,ntot
     &                                  ,scrd,tcrd,uf,tf,npts,npsf,iset)
        endif
      else
        if(sdir.eq.'Y')then !uf=ux, compare y coordinates
          call set_bc_arrays2(ym1,uin,tin,ntot
     &                                  ,scrd,uf,tf,npts,npsf,iset)
        else !uf=uz, compare x coordinates
          call set_bc_arrays2(xm1,uin,tin,ntot
     &                                  ,scrd,uf,tf,npts,npsf,iset)
        endif
      endif

      iset=iglsum(iset,1)
      if(nio.eq.0)then
        write(*,'(7x,a,1x,i9,1x,a)') "set",iset,"points"
        write(*,*)
      endif

      return

 25   if(nio.eq.0) write(*,'(7x,a)')"Unable to open inflow.dat"
      aflag=1
      goto 26

      end
C-----------------------------------------------------------------------
      subroutine set_bc_arrays(x1,y1,u1,t1,n1,x2,y2,u2,t2,n2,n3,iset)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer nmax
      parameter(nmax=max(lelx*max(lely,lelz),lely*lelz)*lx1*lx1)
      real x1(1),y1(1),u1(1),t1(lx1*ly1*lz1*lelt,1)
      real x2(1),y2(1),u2(1),t2(nmax,1)
      integer n1,n2,n3,i,j,k,iset

      do i=1,n1
        do j=1,n2
          if((abs(x1(i)-x2(j)).lt.1.0d-6).and.
     &                                (abs(y1(i)-y2(j)).lt.1.0d-6)) then
            u1(i)=u2(j)
            do k=1,n3
              t1(i,k)=t2(j,k)
            enddo
            iset=iset+1
            goto 27
          endif
        enddo
 27     continue
      enddo

      return
      end
C-----------------------------------------------------------------------
      subroutine set_bc_arrays2(x1,u1,t1,n1,x2,u2,t2,n2,n3,iset)
      implicit none
      include 'SIZE'
      include 'TOTAL'

      integer nmax
      parameter(nmax=max(lelx*max(lely,lelz),lely*lelz)*lx1*lx1)
      real x1(1),u1(1),t1(lx1*ly1*lz1*lelt,1)
      real x2(1),u2(1),t2(nmax,1)
      integer n1,n2,n3,i,j,k,iset

      do i=1,n1
        do j=1,n2
          if(abs(x1(i)-x2(j)).lt.1.0d-6) then
            u1(i)=u2(j)
            do k=1,n3
              t1(i,k)=t2(j,k)
            enddo
            iset=iset+1
            goto 28
          endif
        enddo
 28     continue
      enddo

      return
      end

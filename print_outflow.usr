      subroutine print_outflow(rchar,rloc)
      include 'mpif.h'
      include 'SIZE'
      include 'TOTAL'

      integer lxyz,nmax,lpmax
      parameter(lxyz=lx1*ly1*lz1, lpmax = 128)
      parameter(nmax=max(lelx*max(lely,lelz),lely*lelz)*lx1*lx1)

      character rchar
      real rloc,tmploc
      real rcrd(lx1*ly1*lz1*lelv)
      real scrd(lx1*ly1*lz1*lelv)
      real tcrd(lx1*ly1*lz1*lelv)

      character fname*32,schar,tchar
      character*4 tname(ldimt)
      integer ipoint,npts,aflag,lpts(nmax),gpts(lpmax),gpos(lpmax),ierr
      integer i,i0,i1,j,j0,j1,k,k0,k1,n,npsout
      real norm(3),fnorm(3),dp,vel
      real lscrd(nmax),ltcrd(nmax),lvel(nmax),lt(nmax,ldimt)
      real gscrd(nmax),gtcrd(nmax),gvel(nmax),gt(nmax,ldimt)

      data aflag /0/
      save aflag
      
      if(np.gt.lpmax) then
        if(nio.eq.0) write(6,*)
     &                      "Not configured for more than 128 processes"
        return
      endif
      if(aflag.gt.0) return

      npts=0
      n=nx1*ny1*nz1*nelv

      call cfill(norm,0.0,3)
      if(rchar.eq.'x'.or.rchar.eq.'X')then
        norm(1)=1.0
        schar="Y"
        tchar="Z"
        call copy(rcrd,xm1,n)
        call copy(scrd,ym1,n)
        if(if3d) call copy(tcrd,zm1,n)
      elseif(rchar.eq.'y'.or.rchar.eq.'Y')then
        norm(2)=1.0
        call copy(rcrd,ym1,n)
        if(if3d)then 
          schar="Z"
          tchar="X"
          call copy(scrd,zm1,n)
          call copy(tcrd,xm1,n)
        else
          schar="X"
          call copy(scrd,xm1,n)
        endif
      elseif(rchar.eq.'z'.or.rchar.eq.'Z')then
        norm(3)=1.0
        schar="X"
        tchar="Y"
        call copy(rcrd,zm1,n)
        call copy(scrd,xm1,n)
        call copy(tcrd,ym1,n)
      endif

      npsout=0
      if(ifheat) then
        j=0
        if(ifto) then
          write(tname(1),'(a4)')"temp"
          j=j+1
        endif
        do i=1,npscal
          if(ifpsco(i))then
            write(tname(j+1),'(a2,i2)')"PS",i
            j=j+1
          endif
        enddo
        npsout=j
      endif

      tmploc=get_nearest_face(rloc,rcrd,norm)

      do ielem=1,nelv
      do iside=1,2*ldim
        call facind(i0,i1,j0,j1,k0,k1,lx1,ly1,lz1,iside)
        i=(i0+i1)/2
        j=(j0+j1)/2
        k=(k0+k1)/2
        call getSnormal(fnorm,i,j,k,iside,ielem)
        dp=fnorm(1)*norm(1)+fnorm(2)*norm(2)+fnorm(3)*norm(3)
        if(abs(1.0-dp).lt.1.0d-8) then
          ipoint=i+lx1*(j-1)+lx1*ly1*(k-1)+lxyz*(ielem-1)
          if(abs(rcrd(ipoint)-tmploc).lt.1.0d-8) then
            do i=i0,i1
            do j=j0,j1
            do k=k0,k1
              ipoint=i+lx1*(j-1)+lx1*ly1*(k-1)+lxyz*(ielem-1)
              npts=npts+1
              if(npts.gt.nmax)then
                aflag=1
                goto 25
              endif
              lpts(npts)=ipoint
            enddo
            enddo
            enddo
          endif
        endif
      enddo
      enddo

 25   aflag=iglsum(aflag,1)
      if(aflag.gt.0) then
        if(nid.eq.0) then
          write(*,'(2x,a)')"Warning: Too many points in outflow plane!"
          write(*,'(5x,a)')"no outflow data will be printed!"
        endif
        return
      endif

      do i=1,npts
        ipoint=lpts(i)
        lscrd(i)=scrd(ipoint)
        if(if3d) ltcrd(i)=tcrd(ipoint)
        vel=vx(ipoint,1,1,1)*norm(1)+vy(ipoint,1,1,1)*norm(2)
        if(if3d) vel=vel+vz(ipoint,1,1,1)*norm(3)
        lvel(i) =vel
        if(ifheat) then
          k=1
          if(ifto) then
            lt(i,1)=t(ipoint,1,1,1,1)
            k=k+1
          endif
          do j=1,npscal
            if(ifpsco(j))then
              lt(i,k)=t(ipoint,1,1,1,j+1)
              k=k+1
            endif
          enddo
        endif
      enddo

      call MPI_Gather(npts,1,MPI_INTEGER,gpts,1,MPI_INTEGER,0
     &                                             ,MPI_COMM_WORLD,ierr)
      gpos(1)=0
      do i=2,np
        gpos(i)=gpos(i-1)+gpts(i-1)
      enddo

      call MPI_GatherV(lscrd,npts,MPI_DOUBLE,gscrd,gpts,gpos,MPI_DOUBLE
     &                                           ,0,MPI_COMM_WORLD,ierr)
      if(if3d)
     &call MPI_GatherV(ltcrd,npts,MPI_DOUBLE,gtcrd,gpts,gpos,MPI_DOUBLE
     &                                           ,0,MPI_COMM_WORLD,ierr)
      call MPI_GatherV(lvel ,npts,MPI_DOUBLE,gvel ,gpts,gpos,MPI_DOUBLE
     &                                           ,0,MPI_COMM_WORLD,ierr)
      if(ifheat) then
        do i=1,npsout
          call MPI_GatherV(lt(1,i),npts,MPI_DOUBLE,gt(1,i),gpts,gpos
     &                                ,MPI_DOUBLE,0,MPI_COMM_WORLD,ierr)
        enddo
      endif
      npts=iglsum(npts,1)

      if(nid.eq.0) then
        write(fname,'(a,i3,a)')"outflow.dat"
        open(unit=10,file=fname,status='unknown',form='formatted')
        if(if3d) then
          write(10,'(i6,a14,21a20)')npts,schar,tchar,"velocity"
     &                                            ,(tname(i),i=1,npsout)
          do i=1,npts 
            write(10,'(23(ES20.12))') gscrd(i),gtcrd(i),gvel(i)
     &                                             ,(gt(i,j),j=1,npsout)
          enddo
        else
          write(10,'(i6,a14,21a20)')npts,schar,"velocity"
     &                                            ,(tname(i),i=1,npsout)
          do i=1,npts 
            write(10,'(23(ES20.12))') gscrd(i),gvel(i)
     &                                             ,(gt(i,j),j=1,npsout)
          enddo
        endif
        close(10)
        write(*,'(5x,i6,1x,a,a)')
     &                     npts,"GLL points written to file ",fname
        write(*,*)
      endif

      return
      end

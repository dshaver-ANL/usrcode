c----------------------------------------------------------------------
      subroutine my_mv_mesh(umeshx,umeshy,umeshz)

c     This subroutine solves for and applies the overall mesh 
c     displacement when provided with the displacment vector on boundaries
c     with the 'mv ' BC in field 0.

      include 'SIZE'
      include 'TOTAL'
      real umeshx(lx1,ly1,lz1,lelt)
      real umeshy(lx1,ly1,lz1,lelt)
      real umeshz(lx1,ly1,lz1,lelt)
      parameter (lt = lx1*ly1*lz1*lelt)
      common /mrthoi/ napprx(2),nappry(2),napprz(2)
      common /mrthov/ apprx(lt,0:mxprev)
     $              , appry(lt,0:mxprev)
     $              , apprz(lt,0:mxprev)
      common /mstuff/ d(lt),h1(lt),h2(lt),mask(lt)
      real mask,pmax,pmin
      real srfbl,volbl,delta,deltap1,deltap2,arg1,arg2
      real zero,one
      integer e,f
      integer icalld
      save    icalld
      data    icalld /0/

c     use local arrays to avoid relying on lx1m
      real delx(lx1,ly1,lz1,lelt)
      real dely(lx1,ly1,lz1,lelt)
      real delz(lx1,ly1,lz1,lelt)

      n = nx1*ny1*nz1*nelv
      nface = 2*ndim
      zero = 0.
      one  = 1.

      if (icalld.eq.0) then
        icalld=1
        napprx(1)=0
        nappry(1)=0
        napprz(1)=0
        nxz   = nx1*nz1
        nxyz  = nx1*ny1*nz1
        srfbl = 0.   ! Surface area of elements in b.l.
        volbl = 0.   ! Volume of elements in boundary layer
        do e=1,nelv
        do f=1,nface
          if (cbc(f,e,0).eq.'mv ') then
            srfbl = srfbl + vlsum(area(1,1,f,e),nxz )
            volbl = volbl + vlsum(bm1 (1,1,1,e),nxyz)
          endif
        enddo
        enddo
        srfbl = glsum(srfbl,1)  ! Sum over all processors
        volbl = glsum(volbl,1)
        delta = volbl / srfbl   ! Avg thickness of b.l. elements
c       delta = 0.02            ! 1/2 separation of cylinders
        call rone (h1,n)
        call rzero(h2,n)
     
        call cheap_dist(d,0,'mv ')

        if (nid.eq.0) write(6,*) "delta: ",delta
        deltap1 = 1.0*delta  ! Protected b.l. thickness
        deltap2 = 2.0*delta

c       magic distribution - it really does a better job of preseving BLs 
        do i=1,n
          arg1   = -(d(i)/deltap1)**2
          arg2   = -(d(i)/deltap2)**2
          h1(i)  = h1(i) + 1000.0*exp(arg1) + 10.0*exp(arg2)
        enddo

        call rone(mask,n)
        do e=1,nelv
        do f=1,nface
          if(cbc(f,e,0).eq.'W  ')call facev(mask,e,f,zero,nx1,ny1,nz1)
          if(cbc(f,e,0).eq.'W1 ')call facev(mask,e,f,one ,nx1,ny1,nz1) !! for sides to be moved.
          if(cbc(f,e,0).eq.'v  ')call facev(mask,e,f,zero,nx1,ny1,nz1)
          if(cbc(f,e,0).eq.'mv ')call facev(mask,e,f,zero,nx1,ny1,nz1)
          if(cbc(f,e,0).eq.'O  ')call facev(mask,e,f,zero,nx1,ny1,nz1)
        enddo
        enddo
        call dsop(mask,'*  ',nx1,ny1,nz1)    ! dsop mask
        call opzero(delx,dely,delz)
      endif
  
      do e=1,nelv
      do f=1,nface
        if (cbc(f,e,0).eq.'mv ') then
         call facec(delx,umeshx,e,f,nx1,ny1,nz1,nelv)
         call facec(dely,umeshy,e,f,nx1,ny1,nz1,nelv)
         call facec(delz,umeshz,e,f,nx1,ny1,nz1,nelv)
        endif
      enddo
      enddo
      tol = -1.e-3

      pmax=glamax(delx,n)  
      if (nid.eq.0) write(6,*) "delx max: ",pmax 
 
      utx_usr=glamax(umeshx,n)
      uty_usr=glamax(umeshy,n)  
      utz_usr=glamax(umeshz,n)

      if (nid.eq.0) write(6,*) "utx_usr: ",utx_usr
      if (nid.eq.0) write(6,*) "uty_usr: ",uty_usr
      if (nid.eq.0) write(6,*) "utz_usr: ",utz_usr
  
      if (utx_usr.gt.1e-8)
     & call laplaceh('mshx',delx,h1,h2,mask,vmult,1,tol,
     & 500,apprx,napprx)
      if (uty_usr.gt.1e-8) 
     & call laplaceh('mshy',dely,h1,h2,mask,vmult,1,tol,
     & 500,appry,nappry)
      if (utz_usr.gt.1e-8)
     & call laplaceh('mshz',delz,h1,h2,mask,vmult,1,tol,
     & 500,apprz,napprz)

      call add2(xm1,delx,n)
      call add2(ym1,dely,n)
      call add2(zm1,delz,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine laplaceh
     $     (name,u,h1,h2,mask,mult,ifld,tli,maxi,approx,napprox)
c
c     Solve Laplace's equation, with projection onto previous solutions.
c
c     Boundary condition strategy:
c
c     u = u0 + ub
c
c        u0 = 0 on Dirichlet boundaries
c        ub = u on Dirichlet boundaries
c
c        _
c        A ( u0 + ub ) = 0
c
c        _            _
c        A  u0  =   - A ub
c
c        _             _
c       MAM u0  =   -M A ub,    M is the mask
c
c                      _
c        A  u0  =   -M A ub ,  Helmholtz solve with SPD matrix A
c
c        u = u0+ub
c
      include 'SIZE'
      include 'TOTAL'
      include 'CTIMER'
c
      character*4 name
      real u(1),h1(1),h2(1),mask(1),mult(1),approx (1)
      integer   napprox(1)

      parameter (lt=lx1*ly1*lz1*lelt)
      common /scruz/ r (lt),ub(lt)

      logical ifstdh
      character*4  cname
      character*6  name6

      logical ifwt,ifvec

      call chcopy(cname,name,4)
      call capit (cname,4)

      call blank (name6,6)
      call chcopy(name6,name,4)
      ifwt  = .true.
      ifvec = .false.
      isd   = 1
      imsh  = 1
      nel   = nelfld(ifld)

      n = nx1*ny1*nz1*nel

      call copy (ub,u,n)             ! ub = u on boundary
      call dsavg(ub)                 ! Make certain ub is in H1
                                     !     _
      call axhelm (r,ub,h1,h2,1,1)   ! r = A*ub

      do i=1,n                       !        _
         r(i)=-r(i)*mask(i)          ! r = -M*A*ub
      enddo

      call dssum  (r,nx1,ny1,nz1)    ! dssum rhs

      call project1
     $    (r,n,approx,napprox,h1,h2,mask,mult,ifwt,ifvec,name6)

      tol = abs(tli)
      p22=param(22)
      param(22)=abs(tol)
      if (nel.eq.nelv) then
        call hmhzpf (name,u,r,h1,h2,mask,mult,imsh,tol,maxi,isd,binvm1)
      else
        call hmhzpf (name,u,r,h1,h2,mask,mult,imsh,tol,maxi,isd,bintm1)
      endif
      param(22)=p22

      call project2
     $     (u,n,approx,napprox,h1,h2,mask,mult,ifwt,ifvec,name6)

      call add2(u,ub,n)

      return
      end
C-----------------------------------------------------------------------

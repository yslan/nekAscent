      program eddy2d
      implicit none
      include 'mpif.h'

      integer ldim,lx1,ly1,lz1,lelt,ldimt,lx2,ly2,lz2
      parameter (ldim=2,lelt=1000,lx1=8,ldimt=2)
      parameter (ly1=lx1,lz1=1+(ldim-2)*(lx1-1))

      parameter (lx2=lx1,ly2=lx2,lz2=1+(ldim-2)*(lx2-1))

      logical mpi_is_initialized
      integer nid, np, nekcomm, nekgroup, nekreal
      common /nekmpi/ nid,np,nekcomm,nekgroup,nekreal

      real xm1(lx1,ly1,lz1,lelt)
      real ym1(lx1,ly1,lz1,lelt)
      real zm1(lx1,ly1,lz1,lelt)
      real xm2(lx2,ly2,lz2,lelt)
      real ym2(lx2,ly2,lz2,lelt)
      real zm2(lx2,ly2,lz2,lelt)

      real vx (lx1,ly1,lz1,lelt)
      real vy (lx1,ly1,lz1,lelt)
      real vz (lx1,ly1,lz1,lelt)
      real pr (lx2,ly2,lz2,lelt)
      real pm1(lx1,ly1,lz1,lelt)
      real t  (lx1,ly1,lz1,lelt,ldimt)

      ! common block is needed to have fixed address space
      common /nelsol/ xm1,ym1,zm1,xm2,ym2,zm2,vx,vy,vz,pr,pm1,t 


      integer istep, nsteps, nelx, nely, nelz
     $      , nelgv, nelv, nelt, nfldt, ismvbd
      real time, dt

      integer n, n2, ierr, nel0, nel1
      real visc, u0, v0

      integer e,ix,iy ! dbg

      call mpi_initialized(mpi_is_initialized, ierr)
      if (.not.mpi_is_initialized) call mpi_init(ierr)

      call mpi_comm_dup(MPI_COMM_WORLD,nekcomm,ierr)
      call mpi_comm_size(nekcomm,np,ierr)
      call mpi_comm_rank(nekcomm,nid,ierr)
     

      nelx = 10 !10 ! 10 x 10 box
      nely = nelx
      nelz = nelx

      nelgv = nelx*nely
      if (ldim.eq.3) nelgv=nelgv*nelz

      call set_nelv(nelv,nel0,nel1,nelgv,nid,np)
      nelt = nelv

      nfldt = 1
      ismvbd = 0 ! moving mesh

      nsteps = 1000
      dt = 1.e-2

      visc = 1./20.
      u0 = 1.0
      v0 = 0.3

      n    = lx1*ly1*lz1*nelv
      n2   = lx2*ly2*lz2*nelv

      call gen_mesh(xm1,ym1,zm1,lx1,ly1,lz1,nelx,nely,nelz
     $             ,ldim,nel0,nel1)
      call copy(xm2,xm1,n)
      call copy(ym2,ym1,n)
      call copy(zm2,zm1,n)

      call fnekascent_setup(nekcomm,
     $                      ldim, lx1, ly1, lz1, nelt, lelt, nfldt,
     $                      xm1, ym1, zm1,
     $                      vx, vy, vz, pm1, t)

      time = 0.0
      do istep=1,nsteps
        call MPI_Barrier(nekcomm, ierr)
        time = time + dt
        if (nid.eq.0) write(*,*) istep,time

        call exact  (vx,vy,xm1,ym1,n,time,visc,u0,v0)
        call exactp (pr,xm2,ym2,n2,time,visc,u0,v0)
        call copy(pm1,pr,n)

        call fnekascent_update(time, istep, ismvbd);
      enddo

      call fnekascent_finalize()

      call mpi_finalize(ierr)

      end
c-----------------------------------------------------------------------
      subroutine set_nelv(nelv,nel0,nel1,nelgv,nid,np)
      implicit none
      integer nelv,nel0,nel1,nelgv,nid,np,ntmp

      nelv = nelgv / np
      ntmp = nelgv - nelv*np

      if (nid.lt.ntmp) then
         nel0 = nid * (nelv+1)      ! total #nel before this rank
         nelv = nelv + 1
         nel1 = nel0 + nelv         ! total #nel including this rank
      else 
         nel0 = ntmp * (nelv+1) + (nid-ntmp) * nelv
         nel1 = nel0 + nelv
      endif

      write(*,20)'nel chk, np',np,' nelg',nelgv
     $          ,', nid',nid,' nel',nelv,nel0,nel1

   20 format(ai3ai6ai3a,3(i6))

      return
      end
c-----------------------------------------------------------------------
      subroutine gen_mesh(x,y,z,lx1,ly1,lz1,nelx,nely,nelz
     $                   ,ndim,nel0,nel1)
c     Generate uniform mesh in [0,2,pi] with chebyshev nodes (with endpt)
      implicit none
      integer lx1,ly1,lz1,nelx,nely,nelz
     $      , ndim,nel0,nel1,ix,iy,iz,ex,ey,ez,eg,e
      real xmin,xmax,xlen,xdel,ydel,zdel,x0,x1,y0,y1,z0,z1
      real x(lx1,ly1,lz1,1)
      real y(lx1,ly1,lz1,1)
      real z(lx1,ly1,lz1,1)
      real pi
      data pi /3.14159265358979323846264338/

      real zcheb(lx1)

      xmin = 0.0
      xmax = 2*pi
      xlen = xmax - xmin
      xdel = xlen / nelx

      ydel = xdel
      xdel = xdel

      do ix=1,lx1 ! zcheb in [0,2]
        zcheb(ix) = cos( (lx1-ix)*pi / real(lx1-1.0) ) + 1.0
      enddo
      

      if (ndim.eq.3) then

        do ez=1,nelz
        do ey=1,nely
        do ex=1,nelx
          eg = (ez-1)*nely*nelx + (ey-1)*nelx + ex
          if (eg.gt.nel0.AND.eg.le.nel1) then
            e = eg - nel0

            z0 = xmin + (ez-1)*zdel
            z1 = xmin + (ez)*zdel
            y0 = xmin + (ey-1)*ydel
            y1 = xmin + (ey)*ydel
            x0 = xmin + (ex-1)*xdel
            x1 = xmin + (ex)*xdel
  
            do iz=1,lz1
            do iy=1,ly1
            do ix=1,lx1
              x(ix,iy,iz,e) = x0 + zcheb(ix)*(x1-x0)*0.5
              y(ix,iy,iz,e) = y0 + zcheb(iy)*(y1-y0)*0.5
              z(ix,iy,iz,e) = z0 + zcheb(iz)*(z1-z0)*0.5
            enddo
            enddo
            enddo
          endif
        enddo
        enddo
        enddo

      else

        do ey=1,nely
        do ex=1,nelx
          eg = (ey-1)*nelx + ex
          if (eg.gt.nel0.AND.eg.le.nel1) then
            e = eg - nel0

            y0 = xmin + (ey-1)*ydel
            y1 = xmin + (ey)*ydel
            x0 = xmin + (ex-1)*xdel
            x1 = xmin + (ex)*xdel
  
            do iy=1,ly1
            do ix=1,lx1
              x(ix,iy,1,e) = x0 + zcheb(ix)*(x1-x0)*0.5
              y(ix,iy,1,e) = y0 + zcheb(iy)*(y1-y0)*0.5
            enddo
            enddo
          endif
        enddo
        enddo

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine copy(a,b,n)
      implicit none
      real a(1), b(1)
      integer i, n
      do i=1,n
        a(i) = b(i)
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine exact(uu,vv,xx,yy,n,time,visc,u0,v0)
      implicit none
      real uu(n),vv(n),xx(n),yy(n)
      real cpsi(2,5), a(2,5)
      save cpsi     , a
      data a / -.2,-.2, .25,0.,   0,0  ,  0,0  ,  0,0  /
      data cpsi / 0, 5 ,  3, 4 ,  0,0  ,  0,0  ,  0,0  /
      real time, visc, u0, v0, one, pi, aa, arg, e
      real x, y, sx, sy, cx, cy, u, v, c1, c2
     $    , c1x, c1y, c2x, c2y, s1x, s1y, s2x, s2y
      integer i, k, n

      one   = 1.
      pi    = 4.*atan(one)

      aa    = cpsi(2,1)**2
      arg   = -visc*time*aa  ! domain is [0:2pi]
      e     = exp(arg)

      do i=1,n
         x = xx(i) - u0*time
         y = yy(i) - v0*time

         sx = sin(cpsi(2,1)*x)
         cx = cos(cpsi(2,1)*x)
         sy = sin(cpsi(2,1)*y)
         cy = cos(cpsi(2,1)*y)
         u  =  a(1,1)*cpsi(2,1)*cy 
         v  =  a(2,1)*cpsi(2,1)*sx

         do k=2,5
            s1x = sin(cpsi(1,k)*x)
            c1x = cos(cpsi(1,k)*x)
            s2x = sin(cpsi(2,k)*x)
            c2x = cos(cpsi(2,k)*x)

            s1y = sin(cpsi(1,k)*y)
            c1y = cos(cpsi(1,k)*y)
            s2y = sin(cpsi(2,k)*y)
            c2y = cos(cpsi(2,k)*y)
            
            c1  = cpsi(1,k)
            c2  = cpsi(2,k)

            if (k.eq.2) u = u + a(1,k)*s1x*c2y*c2
            if (k.eq.2) v = v - a(1,k)*c1x*s2y*c1
            if (k.eq.2) u = u - a(2,k)*s2x*c1y*c1
            if (k.eq.2) v = v + a(2,k)*c2x*s1y*c2

            if (k.eq.3) u = u - a(1,k)*s1x*c2y*c2
            if (k.eq.3) v = v + a(1,k)*c1x*s2y*c1
            if (k.eq.3) u = u - a(2,k)*c2x*c1y*c1
            if (k.eq.3) v = v - a(2,k)*s2x*s1y*c2

            if (k.eq.4) u = u + a(1,k)*c1x*c2y*c2
            if (k.eq.4) v = v + a(1,k)*s1x*s2y*c1
            if (k.eq.4) u = u + a(2,k)*c2x*c1y*c1
            if (k.eq.4) v = v + a(2,k)*s2x*s1y*c2

            if (k.eq.5) u = u - a(1,k)*s1x*c2y*c2
            if (k.eq.5) v = v + a(1,k)*c1x*s2y*c1
            if (k.eq.5) u = u - a(2,k)*s2x*c1y*c1
            if (k.eq.5) v = v + a(2,k)*c2x*s1y*c2
         enddo
         uu(i) = u*e + u0
         vv(i) = v*e + v0
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine exactp(pe,x2,y2,n2,time,visc,u0,v0)
      implicit none
      integer n2
      real x2(n2),y2(n2),pe(n2)
      real time, visc, u0, v0

      real e, x, y
      integer i

      e = exp(-50*time*visc)

      do i=1,n2
         x = x2(i) - u0*time
         y = y2(i) - v0*time

         pe(i) = (1.0/64.0)*e*(16*cos(6*x) + 8*cos(8*x-4*y)
     $         - 32*cos(2*(x-2*y)) + 9*cos(8*y) - 8*cos(4*(2*x+y))
     $         + 32*cos(2*(x+2*y)) - 4*sin(3*(x-3*y)) + 32*sin(5*(x-y))
     $         + 36*sin(3*x-y) - 32*sin(5*(x+y)) + 36*sin(3*x+y)
     $         - 4*sin(3*(x+3*y)))
      enddo

c     call ortho(pe)

      return
      end


!     3d semi-analytical streamline tracking code including diffusion based on
!     Pollock, D.W., Semianalytical Computation of Path Lines for
!     Finite-Difference Models. Ground Water, 1988. 26(6): p.743-750.
!
!     Daniel W. Meyer, February 2017
!     Christop Schultheiss, April 2017
!     Institute of Fluid Dynamics
!     ETH Zurich
!
      program streamlinesnt3d
        use pollock_lib, only: pollock, getg, getv
        use streamlinesnt3d_lib
        use omp_lib
        implicit none
!       simulation parameters
        integer, parameter :: d = 3 ! spatial dimensionality
        integer :: mode ! 0 = generate one s.l., 1,2 = transport sampling from inflow or void cells
        character*(256) :: Q1file, Q2file, Q3file ! name of files with velocity data
        character*(256) :: Gfile ! name of file with geometry data
        character*(256) :: Ofile ! name of output file
        real :: Uref, Tref ! reference velocity [L/T] and time [T]
        real :: Dm ! molecular diffusion coefficient [L^2/T]
        integer :: slsmax ! maximum number of streamline sections (0 if unlimited)
        real :: tmax ! maximum streamline travel time (0 if unlimited) [T]
        integer :: csf ! streamline cell storage frequency, to reduce file size in s.l. mode
        real :: dx(d) ! grid cell size [L]
        logical :: rflag, oflag ! rotate domain when sampling streamline sections? Mask output?
        integer :: n(d) ! number of grid cells
        ! variables transport sampling mode
        integer :: np, nt ! number of particles and time snapshots
        real, allocatable :: tl(:) ! snapshot time list
!       other variables
        integer, parameter :: bnd(d,2) = reshape((/3, 1, 1, 3, 1, 1/), shape(bnd)) ! bnd's for pollock
        real, allocatable :: v(:,:) ! velocity [L/T]
        logical, allocatable :: g(:) ! geometry [.true. = solid, .false. = void]
        integer, allocatable, dimension(:) :: gl ! geometry data for reading [1 = solid, 0 = void]
        real :: ts, xs(d) ! cumulative time [T] and position [L] along one streamline section
        real :: xsi(d) ! s.l. section start point within start cell [L]
        integer:: ic(d) ! starting point cell index
        real, pointer :: dt(:,:) ! streamline time increments [T]
        real, pointer, dimension(:,:) :: xi, u, a ! s.l. incr. [L], veloc. [L/T], & accel. [L/T^2]
        real, pointer :: ns(:,:), trace(:,:) ! counter of steps in travelcell1d, indices of cells visited
        real :: xe(d) ! s.l. position within last cell of section [L]
        integer :: m ! number of streamline cells in section
        logical, allocatable, dimension(:,:,:) :: visited ! true if visited by a s.l.
        logical :: pflag ! outcome pollock algorithm
        integer :: i1, i2, i3, i, j, k ! loop counters
        logical:: dflag ! direction flag for s.l. section interfaces
        real :: vi(d) ! velocity to find next s.l. section start point [L]
        real :: R(d,d), alpha ! rotation matrix and angle
        double precision :: xc(d), tc, lc ! cumulative s.l. point [L], time [T], and length [T]
        integer :: idum ! seed for random number generator
        integer :: date_time(8) ! used for generation of random seed
        integer :: tid ! thread id
        ! variables transport sampling mode
        double precision, allocatable :: xp(:,:,:), l(:,:) ! particle point, s.l. length [L]
        integer :: it, caseno ! snapshot time index, caseno for plume files
        ! AH 2019
        integer :: bnd_shared(d,2) = bnd
        integer :: d_shared = d
        real    :: t_start, t_finish, t_load_1, t_load_2
        call cpu_time(t_start)
        ! AH 2019
!
        print '(a)', 'running streamlinesnt3d (dwm&cs,ETHZ,2017)'
!       get user input
        print '(a,/,a,/a)', 'enter mode,', ' (0) streamline generation from void cell,', &
          ' (1,2) transport sampling from inflow or void cells:'
        read(unit=*,fmt=*) mode
        print '(a)', 'enter names of files with velocity distribution:'
        read(unit=*,fmt='(a)') Q1file
        read(unit=*,fmt='(a)') Q2file
        read(unit=*,fmt='(a)') Q3file
        print '(a)', 'enter name of file with geometry data:'
        read(unit=*,fmt='(a)') Gfile
        print '(a)', 'enter molecular diffusion coefficient [L^2/T]:'
        read(unit=*,fmt=*) Dm
        print '(a)', 'enter reference velocity [L/T] and time [T]:'
        read(unit=*,fmt=*) Uref, Tref
        print '(a)', 'enter number of grid cells in each direction:'
        read(unit=*,fmt=*) n
        print '(a)', 'enter grid cell size [L]:'
        read(unit=*,fmt=*) dx(1); dx = dx(1)
        print '(a)', 'rotate domain when tracking streamline (1=yes, 0=no):'
        read(unit=*,fmt=*) i1; rflag = (i1 == 1)
        if (mode == 0) then ! streamline generation
          print '(a)', 'enter name of streamline output file:'
          read(unit=*,fmt='(a)') Ofile
          np = 1 ! just track one streamline
          print '(a)', 'enter maximum number of downstream s.l. sections (0 for no limit):'
          read(unit=*,fmt=*) slsmax
          print '(a)', 'enter maximum streamline travel time (0 for no limit) [T]:'
          read(unit=*,fmt=*) tmax
          print '(a)', 'enter s.l. cell storage frequency (e.g. 2 for every 2nd cell):'
          read(unit=*,fmt=*) csf
        else ! transport sampling
          print '(a)', 'enter number of particles:'
          read(unit=*,fmt=*) np
          print '(a)', 'enter number of time snapshots:'
          read(unit=*,fmt=*) nt
          allocate(tl(nt))
          print '(a)', 'enter snapshot times:'
          read(unit=*,fmt=*) tl
          slsmax = 0; tmax = tl(nt)
        end if
        !AH 2019
        print '(a)', 'Mask stdio output ? (1=yes, 0=no):'
        read(unit=*,fmt=*) i2; oflag = (i2 == 1)     
        !AH 2019 
!       echo parameters
        print '(/,a)', 'input parameters:'
        print '("velocity data files",/,a,/,a,/,a)', &
          trim(Q1file), trim(Q2file), trim(Q3file)
        print '("geometry data file",/,a)', trim(Gfile)
        print '("molecular diffusion coefficient [L^2/T]",es10.3)', Dm
        print '("Uref [L/T] ",es10.4,", Tref [T] ",es10.4)', Uref, Tref
        if (Dm > 0.0) print '("Pe = Uref^2*Tref/Dm = ",es10.4)', Uref**2*Tref/Dm
        print '("number of grid cells ",2(i0,"x"),i0)', n
        print '("grid cell size [L] ",es10.3)', dx(1)
        print '("domain rotation during streamline tracking ",l1)', rflag
        if (mode == 0) then ! streamline generation
          print '(a)', 'streamline generation mode'
          print '("streamline output file """,a,"""")', trim(Ofile)
          print '("maximum number of downstream s.l. sections ",i0)', slsmax
          print '("maximum streamline travel time [T] ",es10.3)', tmax
          print '("streamline cell storage frequency ",i0)', csf
        else ! transport sampling
          print '(a)', 'transport sampling mode'
          if (mode == 2) print '(a)', 'launching streamlines from void cells'
          print '("np = ",i0,", nt = ",i0)', np, nt
          print '("snapshot times:")'
          do j = 1, nt; print '(es9.3)', tl(j); end do
          call openandwriteplumefiles(caseno)
          print '("case number of plume files ",i0)', caseno
        end if
        print * ! line break
!       read geometry data from file
        call cpu_time(t_load_1) !AH 2019 measuring loading time
        print '("reading geometry data...")'
        allocate(g(product(n)), gl(n(1)))
        open(unit=4,file=Gfile,form='formatted',status='old',action='read')
        do i3 = 1, n(3)
          do i2 = 1, n(2)
            read(unit=4,fmt=*) (gl(i1), i1=1,n(1))
            do i1 = 1, n(1)
              g(n(1)*n(2)*(i3-1)+n(1)*(i2-1)+i1) = (gl(i1) == 1)
            end do
          end do
        end do
        close(unit=4)
        deallocate(gl)
        call cpu_time(t_load_2) !AH 2019 measuring loading time
        print '("Loading time = ",f6.3," seconds.")',t_load_2-t_load_1 !AH 2019
!       read velocity data from file
        allocate(v(product(n+1),d)); v = 0.0
        !
        print '("reading u3 velocity data...")'
        call cpu_time(t_load_1) !AH 2019 measuring loading time
        open(unit=1,file=Q3file,form='formatted',status='old',action='read')
        do i3 = 1, n(3)
          do i2 = 1, n(2)
            read(unit=1,fmt=*) (v((n(1)+1)*(n(2)+1)*(i3-1)+(n(1)+1)*(i2-1)+i1,3), i1=1,n(1))
          end do
        end do
        close(unit=1)
        call cpu_time(t_load_2) !AH 2019 measuring loading time
        print '("Loading time = ",f6.3," seconds.")',t_load_2-t_load_1 !AH 2019
        !
        print '("reading u2 velocity data...")'
        call cpu_time(t_load_1) !AH 2019 measuring loading time
        open(unit=2,file=Q2file,form='formatted',status='old',action='read')
        do i3 = 1, n(3)
          do i2 = 1, n(2)
            read(unit=2,fmt=*) (v((n(1)+1)*(n(2)+1)*(i3-1)+(n(1)+1)*(i2-1)+i1,2), i1=1,n(1))
          end do
        end do
        close(unit=2)
        call cpu_time(t_load_2) !AH 2019 measuring loading time
        print '("Loading time = ",f6.3," seconds.")',t_load_2-t_load_1 !AH 2019
        !
        print '("reading u1 velocity data...")'
        call cpu_time(t_load_1) !AH 2019 measuring loading time
        open(unit=3,file=Q1file,form='formatted',status='old',action='read')
        do i3 = 1, n(3)
          do i2 = 1, n(2)
            read(unit=3,fmt=*) (v((n(1)+1)*(n(2)+1)*(i3-1)+(n(1)+1)*(i2-1)+i1,1), i1=1,n(1))
            v((n(1)+1)*(n(2)+1)*(i3-1)+(n(1)+1)*(i2-1)+n(1)+1,1) &
              = getv(d ,n+1, v, (/n(1),i2,i3/),1) &
              + getv(d ,n+1, v, (/n(1),i2,i3/),2) - getv(d ,n+1, v, (/n(1),i2+1,i3/),2) &
              + getv(d ,n+1, v, (/n(1),i2,i3/),3) - getv(d ,n+1, v, (/n(1),i2,i3+1/),3)
          end do
        end do
        close(unit=3)
        call cpu_time(t_load_2) !AH 2019 measuring loading time
        print '("Loading time = ",f6.3," seconds.")',t_load_2-t_load_1 !AH 2019
!       make things non-dimensional
        print '("making things non-dimensional...")'
        dx = dx / (Uref*Tref); v = v/Uref; Dm = Dm/(Uref**2*Tref)
!       inspect dead-end cells (flow sinks)
        j = 0; k = 0 ! counters
        do i3 = 1, n(3)
          do i2 = 1, n(2)
            do i1 = 1, n(1)
              if (.not.getg(d, n, g, (/i1,i2,i3/))) then ! identify dead-end cells in void space
                k = k+1 ! count void cells
                if ((getv(d,n+1,v,(/i1,i2,i3/),1)>=0.0).and.(getv(d,n+1,v,(/i1+1,i2,i3/),1)<=0.0).and. &
                    (getv(d,n+1,v,(/i1,i2,i3/),2)>=0.0).and.(getv(d,n+1,v,(/i1,i2+1,i3/),2)<=0.0).and. &
                    (getv(d,n+1,v,(/i1,i2,i3/),3)>=0.0).and.(getv(d,n+1,v,(/i1,i2,i3+1/),3)<=0.0)) then
                  j = j+1 ! count dead-end cells
                end if
              end if
            end do
          end do
        end do
        print '("found ",i0," dead-end cells in ",i0," void cells")', j, k
!       track streamline
        allocate(visited(n(1),n(2),n(3)))
        if (mode /= 0) then
          allocate(xp(np,nt,d), l(np,nt))
          l = -1.0 ! set to -1 to determine which particles have reached snapshot time
          xp = 0.0
        end if
!$omp   parallel if(mode /= 0) default(private) &
!$omp   shared(mode, Ofile, Dm, tmax, slsmax, csf, rflag, oflag, np, nt, tl, &
!$omp   n, dx, v, g, bnd_shared, d_shared, visited, xp, l, caseno)
!       initialize things
        tid = omp_get_thread_num()
        call date_and_time (values = date_time)
        idum = date_time(8) + 1000*(date_time(7) + 60*(date_time(6) &
            + 60*date_time(5))) + 1000*tid
        !idum = 80893252 + 1000*tid ! for debugging
        idum = -idum
        print '("random seed ",i0," (thread ",i0,")")', idum, tid
        i = int(ran1(idum))
        visited = .false.
        if (mode == 0) open(unit=2,file=Ofile,form='formatted',action='write')
!       transport sampling or tracking of one streamline
!$omp   do schedule(dynamic) reduction(.or.: visited)
        do i = 1, np ! loop over particles
!         initialize sample angle for current s.l. section
          alpha = 0.0 ! in the case of no rotation
          if ((mode /= 0).and.rflag) then ! random initial angle for sampling mode
            alpha = 2.0 * 4.0*atan(1.0) * ran1(idum) ! 2*pi* random number
          end if
!         pick streamline origin
          if (mode == 1) then
            i1 = 1
            call pickrandomoriginatinflow(v, g, n, idum, ic)
          else
            call pickrandomorigininvoid(g, n, idum, ic)
          end if
          print '("s.l. ",i0," from cell (",i4,", ",i4,", ",i4,a,f4.2," (thread ",i0,")")', &
            i, ic, ') at sample angle ', alpha, tid
!         track streamline sections in (rotated) domain
          k = 1 ! streamline section counter
          xc = 0.0d0; tc = 0.0d0; lc = 0.0d0 ! initialize current streamline point
          it = 1 ! current snapshot time index for transport sampling
          do ! loop over streamline sections
            ts = 0.0 ! initialize section time
            vi = (/(getv(d,n+1,v,(/ic(1),ic(2),ic(3)/),j), j = 1,d)/) ! initial velocity
!           add random noise to starting point
            if (mode == 1) then ! to avoid recurring streamline patterns
              xsi = dx*(/0.0,ran1(idum),ran1(idum)/)
            else ! to arrive at uniformly distributed intial positions
              xsi = dx*(/ran1(idum),ran1(idum),ran1(idum)/)
            end if
            xs = dx*(ic-(/1,1,1/)) + xsi
!           streamline tracking
            if (mode == 0.and.oflag) then !AH 2019 (trying to reduce IO)
              write(*,fmt='("(",es10.3,2(",",es10.3),") ",3(es10.3,1x),i0,1x)',advance="no") &
                xc, sqrt(vi(1)**2 + vi(2)**2 + vi(3)**2), tc, lc, k
            end if
            pflag = pollock(d_shared, n, dx, v, g, bnd_shared, ic, xsi, Dm, 0, real(tmax-tc), &
                            m, dt, xi, u, a, ns, trace, xe, idum)
            if (mode == 0.and.oflag) print '(i0)', m !AH 2019 (trying to reduce IO)
            ic = int(trace(:,m)) ! index of last cell in s.l. section
            dflag = (ic(1) == n(1)) ! .true. if particle exits at outflow boundary
!           store streamline
            ! setup rotation such that velocity over consecutive streamline sections stays smooth
            R(1,:) = (/1.0, 0.0, 0.0/)
            R(2,:) = (/0.0, cos(alpha), sin(alpha)/)
            R(3,:) = (/0.0, -sin(alpha), cos(alpha)/)
            ! rotate start point
            xs = matmul(R, xs); vi = matmul(R, vi)
            ! write streamline section header
            if (mode == 0) then
              write(unit=2,fmt='(i0,2(1x,i0))') 2, 1+3*d, 2 + int((m+1 - 2)/csf) ! array size
              write(unit=2,fmt='(4(es15.8,1x),6(es12.5,1x))') 0.0, xs, vi, (/(0.0,j=1,d)/) ! start
            end if
            ! store streamline section or snapshots depending on mode
            do j = 1, m ! loop over streamline cells
              ! rotate streamline section
              xi(:,j) = matmul(R, xi(:,j))
              u(:,j) = matmul(R, u(:,j)) 
              a(:,j) = matmul(R, a(:,j)) 
              ! update current cumulative time & distance
              tc = tc + dt(1,j)
              ts = ts + dt(1,j)
              lc = lc + norm2(xi(:,j))
              xc = xc + xi(:,j)
              xs = xs + xi(:,j)
              if (mode == 0) then ! write streamline section to file
                if ((mod(j,csf) == 0).or.(j == m)) &
                  write(unit=2,fmt='(4(es15.8,1x),6(es12.5,1x))') ts, xs, u(:,j), a(:,j)
              else ! extract particle positions at snapshot times
                if ((.not.pflag).and.(j == m)) then ! trapped particle
                  do while (it <= nt)
                    xp(i,it,:) = xc; l(i,it) = lc
!$omp               critical
                    call openandwriteplumefiles(caseno, it, (/xp(i,it,:), l(i,it)/))
!$omp               end critical  
                    it = it+1 ! increment snapshot time counter
                  end do    
                else ! normal case
                  do while(real(tl(it)-tc) < 0.0)
                    xp(i,it,:) = xc - xi(:,j)*(tc - tl(it))/dt(1,j)
                    l(i,it) = lc - (norm2(xi(:,j)))*(tc- tl(it))/dt(1,j)
!$omp               critical
                    call openandwriteplumefiles(caseno, it, (/xp(i,it,:), l(i,it)/))
!$omp               end critical
                    it = it + 1; if (it > nt) exit ! increment snapshot time counter
                  end do
                end if
                ! finish when largest snapshot time has been reached
                if (it > nt) then
                  print '("s.l. ",i0," completed after ",i0," sections (thread ",i0,")")', i, k, tid
                  exit
                end if
              end if ! mode
            end do ! loop over streamline cells to store streamline
!           end current s.l. section
            vi = u(:,m) ! initial velocity of next streamline section
            do j = 1, m
!$omp         critical
              visited(int(trace(1,j)),int(trace(2,j)),int(trace(3,j))) = .true.
!$omp         end critical
            end do
            deallocate(dt, xi, u, a, ns, trace) ! clean up
            if (mode /= 0) then ! finish s.l. tracking if largest snapshot time has been reached
              if (it > nt) exit
            else ! number of downstream s.l. sections completed
              if ((k >= slsmax).and.(slsmax /= 0.0)) exit
            end if
            ! stop s.l. in case of pollock failure
            if (.not.pflag) then
              print '("warning: pollock failed, ending s.l. ",i0," (thread ",i0,")")', i, tid
              exit
            end if
            ! tmax exceeded
            if ((tc > tmax).and.(tmax /= 0.0)) then
              print '("s.l. ",i0," ended due to max. time (",es10.3,") thread ",i0)', i, tc, tid
              exit
            end if
!           prepare next s.l. section
            k = k + merge(1,-1,dflag) ! update streamline section counter
            call startpointnextsection(v, g, n, vi, rflag, ic, alpha, dflag)
          end do ! loop over s.l. sections
        end do ! loop over particles
!$omp	end do
!       report key streamline quantities
        if ((mode == 0).and.pflag) then
          print '(a,2(1x,es10.3),1x,i0)', &
            's.l. length, time, and downstream sect. count are ', lc, tc, k
          print '(a,3(1x,es10.3))', 's.l. end point is', xc
          print '(a,es11.4)', 'resulting normalized mean x1-velocity ', xc(1)/tc
          print '(a,es11.4)', 'tortuosity is ', lc / xc(1)
        end if
!$omp   end parallel
!       report how many void cells have been visited
        j = 0; k = 0
        do i3 = 1, n(3)
          do i2 = 1, n(2)
            do i1 = 1, n(1)
              if (visited(i1,i2,i3)) k = k+1 ! count visited cells
              if (.not.getg(d, n, g, (/i1,i2,i3/))) j = j+1 ! count void cells
            end do
          end do
        end do
        print '(f5.2,"% of ",i0," void cells have been visited by s.l. sections")', &
          100.0 * k/j, j
!       clean up
        if (mode == 0) close(unit=2)
        deallocate(v, g, visited)
        if (mode /= 0) deallocate(tl, xp, l)
        print '(a)', 'streamlinesnt3d completed'
        !AH 2019
        call cpu_time(t_finish)
        print '("CPU Time = ",f6.3," seconds.")',t_finish-t_start
        !AH 2019
      end program streamlinesnt3d

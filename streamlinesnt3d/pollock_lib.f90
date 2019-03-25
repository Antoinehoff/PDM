!     pollock calculates a single streamline from the starting point xs and molecular
!     diffusion Dm in the velocity field v given on a grid with n cells of
!     size Lc and dimensionality d based on a semi-analytical algorithm inspired by
!     Pollock, D.W., Semianalytical Computation of Path Lines for Finite-Difference
!     Models. Ground Water, 1988. 26(6): p.743-750.
!
!     Daniel W. Meyer, February 2017
!     Institute of Fluid Dynamics
!     ETH Zurich
!
      module pollock_lib
      
        implicit none
        private O, L, interfacev, increasemem, travelcell, travelcell1d, newtonsearch, &
          tea, tdg, xn, un, dawson, gasdev, ran1
        real, parameter :: pi = 3.14159265358979323846
        real, parameter :: infty = huge(infty)
        real, parameter :: eps = epsilon(eps) ! small epsilon such that 1.0 + eps > 1.0
        real, parameter :: fail = -1.0 ! return value if travelcell1d failed [-]
        real, parameter :: O = 0.0, L = 1.0 ! non-dim. cell origin & size [-]

      contains

      logical function pollock(d, n, Lc, v, g, bnd, is, xs, Dm, cmax, tmax, m, dt, dx, u, a, ns, trace, xe, idum)
        implicit none
!       interface variables
        integer, intent(in) :: d ! spatial dimensionality
        integer, intent(in) :: n(d) ! number of grid cells
        real, intent(in) :: Lc(d) ! grid cell size [L]
        real, intent(in) :: v(product(n+1),d) ! cell interface velocities [L/T]
        logical, intent(in) :: g(product(n)) ! geometry [.true. = solid, .false. = void]
        integer, intent(in) :: bnd(d,2) ! boundaries, 1 = solid, 2 = periodic, 3 = in/out
        integer, intent(in) :: is(d) ! streamline starting point, cell index
        real, intent(in) :: xs(d) ! streamline starting point, position within cell is [L]
        real, intent(in) :: Dm ! molecular diffusion coefficient [L^2/T]
        integer, intent(in) :: cmax ! maximum number of streamline cells (<= 0 for no limit)
        real, intent(in) :: tmax ! maximum travel time of streamline (<= 0 for no limit) [T]
        integer, intent(out) :: m ! number of streamline cells
        real, pointer :: dt(:,:), dx(:,:) ! streamline increments in time [T] & space [L]
        real, pointer :: u(:,:), a(:,:) ! streamline velocities [L/T] & accelerations [L/T^2]
        real, pointer :: ns(:,:), trace(:,:) ! counter of steps in travelcell1d, indices of cells visited
        real, intent(out) :: xe(d) ! relative position within last cell [L]
        integer, intent(inout) :: idum ! seed for random number generator
!       other variables
        integer :: i(d), ei ! current cell index, cell counter, exit interface & direction
        integer :: mem ! current number of steps in streamline arrays
        real :: xi(d) ! particle position before/after travelcell [L]
        real :: zeta(d) ! diffusion step [-]
        real :: dtz, dtze ! diffusion step time and elapsed time [T]
        real :: vl(d), vr(d) ! grid cell velocities [L/T]
        double precision :: t ! current streamline travel time [T]
        integer :: j, k ! loop counters
!       initialize
        dtz = 0.0; dtze = dtz
        ! cell index and initial position
        i = is; xi = xs
        ei = 0 ! for simplicity, stay away from interfaces initially
        where (xi < eps) xi = eps
        where (Lc*(1.0-eps) < xi) xi = Lc*(1.0-eps)
!       march in time
        mem = 0; m = 0; t = 0.0
        do ! loop over streamline cells
!         extend memory as streamline gets longer
          if (mem <= m) then
            j = ceiling(product(n)**(1.0/d))
            call increasemem(mem, 1, j, dt)
            call increasemem(mem, d, j, dx)
            call increasemem(mem, d, j, u)
            call increasemem(mem, d, j, a)
            call increasemem(mem, 1, j, ns)
            call increasemem(mem, d, j, trace)
            mem = mem+j
          end if
!         determine next cell (if at cell interface)
 1        if ((Dm > 0.0).and.(dtze >= dtz)) then ! initialize diffusion
            ! set zeta vector (non-zero consistent with next step)
            do
              zeta = (/(gasdev(idum),k=1,d)/)
              if (.not.any(zeta == 0.0)) exit
            end do
            ! set diffusion step time
            call interfacev(d, n, v, i, vl, vr)
            dtz = diffusiontime(d, vl, vr, Lc, Dm); dtze = 0.0
            ! next cell and interface indices based on zeta
            if (ei /= 0) then
              i(abs(ei)) = i(abs(ei)) + (ei/abs(ei) + int(sign(1.0,zeta(abs(ei))))) / 2
              ei = abs(ei)*int(sign(1.0,-zeta(abs(ei))))
            end if
          else ! pure advection or diffusion with dtze /= dtz, proceed along flow path
            if (ei /= 0) then
              i(abs(ei)) = i(abs(ei)) + ei/(abs(ei)); ei = -ei
            end if
          end if
!         check boundaries
          do k = 1,d
            if ((i(k) < 1).or.(i(k) > n(k))) then ! particle exits at left or right boundary
              select case (bnd(k,merge(1,2,i(k) < 1)))
              case (1) ! solid
                if (Dm == 0.0) then ! pure advection
                  goto 3 ! dead end
                else ! advection & diffusion
                  zeta(k) = -zeta(k) ! flip zeta
                  i(k) = i(k) + merge(1,-1,i(k) < 1)
                  ei = -ei
                end if
              case (2) ! periodic
                i(k) = merge(n(k),1,i(k) < 1)
              case (3) ! in/out
                goto 2 ! streamline end
              end select
            end if
          end do
!         check geometry
          if (getg(d, n, g, i)) then ! moving into solid?
            if (Dm == 0.0) then ! pure advection
              goto 3 ! dead end
            else ! advection & diffusion
              ! flip zeta if at solid interface
              zeta(abs(ei)) = -zeta(abs(ei))
              i(abs(ei)) = i(abs(ei)) + merge(1,-1,ei > 0)
              ei = -ei
            end if
          end if
!         move particle
          m = m+1 ! streamline cell counter
          trace(:,m) = i ! keep track of cells visited
          ! extract interface velocities from cell i
          call interfacev(d, n, v, i, vl, vr)
          ! put particle on interface
          if (ei /= 0) then
            xi(abs(ei)) = merge(O, Lc(abs(ei)), ei < 0)
          end if
          ! set travel time limit
          dt(1,m) = infty; if (tmax > 0.0) dt(1,m) = real(tmax-t);
          ! move particle in cell i
          call travelcell(d, vl, vr, Dm, zeta, dtz, dtze, Lc, xi, xe, u(:,m), a(:,m), dt(1,m), ei, k)
          ns(1,m) = real(k)
          ! process output of travelcell
          if (dt(1,m) == fail) then ! failure in case with Dm > 0
            m = m-1; dtze = dtz
            goto 1 ! try again with different zeta
          end if
          dx(:,m) = xe-xi
          if (dt(1,m) == infty) goto 3 ! dead end
          t = t+dt(1,m) ! current streamline travel time
!         check for streamline end
          if ((m >= cmax).and.(cmax > 0)) exit ! maximum number of cells
          if ((t >= tmax).and.(tmax > 0.0)) exit ! maximum travel time
          ! prepare next move
          xi = xe ! current particle position
        end do
 2      pollock = .true. ! success
        return
 3      print '("warning: dead end")'
        pollock = .false. ! failure
      end function pollock


!     extract interface velocities vl & vr from cell i
      subroutine interfacev(d, n, v, i, vl, vr)
        implicit none
        integer, intent(in) :: d
        integer, intent(in) :: n(d)
        real, intent(in) :: v(product(n+1),d)
        real, dimension(d), intent(out) :: vl, vr
        integer, intent(in) :: i(d)
        integer :: j, k
        do k = 1,d
          vl(k) = getv(d, n+1, v, i, k)
          vr(k) = getv(d, n+1, v, i + &
            merge((/(1,j=1,d)/),(/(0,j=1,d)/),(/(j,j=1,d)/) == k), k)
        end do
      end subroutine interfacev


!     allocate new larger field, copy old data to new field, deallocate old field
      subroutine increasemem(n, d, dn, field)
        implicit none
        integer, intent(inout) :: n ! current number of entries in field
        integer, intent(in) :: d, dn ! dimensionality of field, numb. of entries to add
        real, pointer :: field(:,:) ! new and old field
        real, pointer :: newfield(:,:) ! dummy for new field
        allocate(newfield(d,n+dn)); newfield = 0.0 ! increase numb. of elements
        if (n /= 0) then ! copy content from old to new field if any
          newfield(1:d,1:n) = field(1:d,1:n) ! copy
          deallocate(field) ! free memory of old field
        end if
        field => newfield ! make increase effective
      end subroutine increasemem


!     access function for d-dimensional velocity field v of size n
      real function getv(d, n, v, i, j)
        implicit none
        integer, intent(in) :: d, n(d) ! dimensionality, size of v
        real, intent(in) :: v(product(n),d) ! field
        integer, intent(in) :: i(d), j ! index vector, component
        integer :: m, k
        m = i(1)
        do k = 2,d
          m = (i(k)-1)*product(n(1:k-1)) + m
        end do
        getv = v(m,j)
      end function getv


!     access function for d-dimensional geometry field g of size n
      logical function getg(d, n, g, i)
        implicit none
        integer, intent(in) :: d, n(d) ! dimensionality, size of g
        logical, intent(in) :: g(product(n)) ! field
        integer, intent(in) :: i(d) ! index vector
        integer :: m, k
        m = i(1)
        do k = 2,d
          m = (i(k)-1)*product(n(1:k-1)) + m
        end do
        getg = g(m)
      end function getg


!     solves dimensional advective/diffusive particle motion in a cell
!     set dt = infty for unlimited tracking in time within one cell
!     return status: dt = infty part. trapped, = fail no solution,
!       otherwise dt is equal to exit time or input dt if exit time > input dt
!     xe, ve, ae are end position, advection-only velocity and acceleration at dt
!     ei = index of exit cell interface, e.g., = -3 for left face on dimension d = 3
!        = 0 if no interface was hit within input dt
!     ns is number of initial guess and Newton steps
      subroutine travelcell(d, vl, vr, Dmi, zetai, dtz, dtze, dx, xi, xe, ve, ae, dt, ei, ns)
        implicit none
        ! interface variables
        integer, intent(in) :: d ! spatial dimensionality
        real, intent(in) :: vl(d), vr(d) ! velocities at left and right cell interfaces [L/T]
        real, intent(in) :: Dmi ! molecular diffusion coefficient [L^2/T]
        real, intent(in) :: zetai(d) ! diffusion vector [-]
        real, intent(in) :: dtz ! diffusion step time [T]
        real, intent(inout) :: dtze ! elapsed time of diffusion step dtz [T]
        real, intent(in) :: xi(d) ! start point at cell border or inside cell in [0,dx] [L]
        real, intent(in) :: dx(d) ! grid cell size [L]
        real, intent(out) :: xe(d) ! trajectory ending or cell exit position [L]
        real, intent(inout) :: dt ! trajectory ending time (in) or cell exit time (out) [T]
        real, intent(out) :: ve(d), ae(d) ! exit velocity [L/T] and acceleration [L/T^2]
        integer, intent(inout) :: ei ! exit face, -k or k for left or right face on dimension k, 0 for in cell
        integer, intent(out) :: ns ! counter of steps in travelcell1d
        ! other variables
        real :: Dm, zeta(d) ! local copies of input arguments Dmi and zetai
        real :: A(d) ! velocity gradient [1/T]
        real :: x1(d), v1(d), Pe(d), te(d), tzeta(d) ! variables used by travelcell1d [-]
        real :: t, v ! time and initial velocity used by un [-]
        integer :: k, eit ! loop counter, exit interface
        ! in case of diffusion, determine section ending time [T]
        if (Dmi > 0.0) dt = min(dt, dtz - dtze)
        ! handle case with no diffusion (keep normalization applicable)
        Dm = Dmi; zeta = zetai; if (Dmi <= 0.0) then; Dm = 1.0; zeta = 0.0; end if
        ! determine shortest dimensional exit time
        eit = 0
        do k = 1,d
          A(k) = (vr(k)-vl(k)) / dx(k)
          if (xi(k) == O) then; v1(k) = vl(k) ! at left interface
          elseif (xi(k) == dx(k)) then; v1(k) = vr(k) ! at right interface
          else; v1(k) = vl(k) + A(k)*xi(k) ! inside cell
          end if
          ! non-dimensionalize
          x1(k) = xi(k)/dx(k); v1(k) = v1(k)*dx(k)/Dm; Pe(k) = A(k)*dx(k)**2/Dm
          tzeta(k) = -dtze / (dx(k)**2/Dm)
          ! calculate non-dimensional exit time along spatial direction k
          call travelcell1d(x1(k), v1(k), Pe(k), zeta(k), tzeta(k), xe(k), te(k), ns)
          ! to be comparable among different k, make time dimensional
          if ((te(k) /= infty).and.(te(k) /= fail)) te(k) = te(k)*dx(k)**2/Dm ! [-] -> [T]
          ! find shortest exit time dt [T] and set exit face
          if ((dt > te(k)).and.(te(k) /= fail)) then
            dt = te(k)
            eit = k*(2*int(xe(k)/L) - 1)
            t = 0.0; v = merge(vl(k),vr(k),eit<0)*dx(k)/Dm
          end if
        end do
        ! abort if travelcell1d failed
        if (any(te == fail)) then; dt = fail; return; end if
        ! set exit interface
        ei = eit
        ! abort if particle is trapped
        if (dt == infty) then; xe = xe*dx; ve = 0.0; ae = 0.0; return; end if
        ! dt = 0 in case where advective & diffusive motion cancel
        if ((dt == 0.0).and.(Dmi > 0.0)) dt = dtz*1.0E-3 ! make tiny time step to avoid infinite loop
        if (dt == 0.0) then; print '("error: dt = 0")'; dt = fail; return; end if ! pure advection
        ! in case of diffusion, update diffusion step time [T]
        if (Dmi > 0.0) dtze = dtze + dt
        ! re-normalize exit time (possibly different among dimensions)
        te = dt * Dm/dx**2 ! [T] -> [-]
        ! determine dimensional exit pos. and advection-only veloc. & accel.
        do k = 1,d
          ! position and velocity
          if (abs(ei) == k) then ! exit face, no need to calculate
            xe(k) = dx(k); if (ei < 0) xe(k) = O
            ve(k) = vr(k); if (ei < 0) ve(k) = vl(k)
          else ! calculate exit position
            ! position
            xe(k) = xn(te(k), x1(k), v1(k), Pe(k), zeta(k), tzeta(k)) ! [-]
            ! clipping, use eps to avoid corners (cross only one interface at a time)
            if ((xe(k) < eps).or.(L*(1.0-eps) < xe(k))) then
              if (max(xe(k)-1,-xe(k)) > 0.01) &
                print '("warning: clipping ",es12.5)', xe(k)
              xe(k) = max(eps,min(L*(1.0-eps),xe(k)))
            end if
            xe(k) = xe(k)*dx(k) ! [-] -> [L]
            ! velocity
            ve(k) = vl(k) + A(k)*xe(k) ! [L/T]
          end if
          ! acceleration
          if (A(k)*dt > log(infty)/2) then
            ae(k) = infty
          elseif (A(k)*dt < -log(infty)/2) then
            ae(k) = 0.0
          else
            v1(k) = vl(k) + A(k)*xi(k) ! start point velocity [L/T]
            ae(k) = A(k) * v1(k) * exp(A(k)*dt) ! from eq.(12a) in Pollock (1988)
          end if
        end do
      end subroutine travelcell


!     solves the non-dimensional advective/diffusive particle motion
!     in a 1d cell [0,1] with the semi-analytical algorithm derived in
!     pollock_with_diffusion.nb
!     the particle starts at point x1 with the velocity v1
!     and a diffusion increment zeta, it travels for the time te
!     until it hits the cell interface xe at either 0 or 1
!     return status: te = infty part. trapped, = fail no solution, otherwise success
      subroutine travelcell1d(x1, v1, Pe, zeta, tzeta, xe, te, ns)
        implicit none
!       interface variables
        real, intent(in) :: x1 ! starting point inside [0,1] in unit cell size [dx]
        real, intent(in) :: v1 ! starting velocity [Dm/dx]
        real, intent(in) :: Pe ! non-dimensional velocity gradient = dv*dx/Dm [-]
        real, intent(in) :: zeta ! standard normal random number [-]
        real, intent(in) :: tzeta ! start time of diffusion step zeta [-] tzeta < 0
        real, intent(out) :: xe, te ! exit position and time [dx, dx^2/Dm]
        integer, intent(out) :: ns ! counter of initial guess and Newton steps
!       other variables
        real :: xa, xs ! exit interface, stationary point
        real :: ta, ta1, ta2 ! advection exit time or guess for te
        integer :: k ! loop counter
        integer, parameter :: nsdmax = 50, nsamax = 1000
!       implementation
        !print '("x1=",es10.3,",v1=",es10.3,",Pe=",es10.3,",zeta=",es10.3,",tzeta=",es10.3)', &
        !   x1,v1,Pe,zeta,tzeta
        ns = 0 ! init step counter
        xs = -L; if (Pe < 0.0) xs = x1 - v1/Pe ! stationary point with u = 0
        !print '("stationary point at ",es10.3," (-1 for no point) ")', xs
        ! pure advection / advection and/or diffusion
        if (zeta == 0.0) then
!         pure advection
          !print '("pure advection")'
          if (v1 == 0.0) then ! trapped
            te = infty; xe = x1
          else ! determine exit position
            xa = O; if (v1 > 0.0) xa = L ! exit interface
            if (x1 == xa) then ! already at exit interface
              te = 0.0; xe = xa
            else ! not at exit position yet
              if ((xs <= O).or.(L <= xs)) then
                ! stationary point outside of [0,1] or Pe >= 0
                ! intersecting with interface, use Pollock solution
                te = tea(x1, v1, Pe, xa); xe = xa
              else ! trapped at stationary point
                te = infty; xe = xs
              end if
            end if
          end if
        else
!         advection and/or diffusion
          xa = O; te = 0.0; if (((tzeta == 0.0).and.(zeta > 0.0)).or. &
                                ((tzeta < 0.0).and.(un(te,v1,Pe,zeta,tzeta) > 0.0))) xa = L
          ! exit interface, diffusion dominates at small t
          if (x1 == xa) then ! already at exit position
            te = 0.0; xe = xa
          else ! not at exit position yet
            if ((v1 == 0.0).and.(Pe == 0.0)) then
!             pure diffusion
              !print '("pure diffusion")'
              te = (x1-xa)*(x1-xa-2*sqrt(-2*tzeta)*zeta) / (2*zeta**2); xe = xa
            else ! advection & diffusion
              if (((zeta*v1 <= 0.0).or.(Pe <= 0.0)).and. &
                (((tzeta < 0.0).and.((xa-x1)*(sqrt(-tzeta)*v1+zeta/sqrt(2.0)) > 0) & ! tdg
                               .and.(abs(2*sqrt(-tzeta)*v1 + sqrt(2.0)*zeta) > eps)).or. & ! tdg
                 ((tzeta == 0.0).and.(zeta**2 >= 2*v1*(x1-xa))))) then ! tdg
!               diffusion at small-time probably dominates
                !print '("diffusion-dominated regime")'
                ! initial guess from root of polynomial approximation
                ta = tdg(x1, v1, zeta, tzeta, xa)
                ! find root with Newton search
                if (newtonsearch(x1, v1, Pe, zeta, tzeta, xa, ta, te, k, nsdmax)) then ! success
                  xe = xa
                else ! failed
                  !print '("failed")'
                  te = infty
                end if
                ns = ns+k ! Newton steps
              else ! proceed with advection-dominated large-time regime
                te = infty
              end if
              if (te >= infty) then
!               advection-dominated large-time regime
                !print '("advection-dominated regime")'
                if ((xs < O).or.(L < xs)) then ! stationary point outside of [0,1] or Pe >= 0?
                  ! intersecting with interface, use Newton search
!                 find good initial guess for exit time ta first
                  xa = O; if ((v1 > 0.0).or.((v1 == 0.0).and.(zeta > 0.0))) xa = L ! exit interface
                  ta = tea(x1, v1, Pe, xa) ! advection exit time
                  if ((ta <= 0.0).or.(ta >= infty)) then ! in case of useless initial guess from tea
                    ta = 1 / sqrt(v1**2 + Pe*(v1 - 2*v1*x1 + Pe*(1.0/3 - x1 + x1**2))) ! 1/(rms of u(x))
                  end if
                  ! find point outside of domain as initial guess for Newton search
                  ! first, increase ta to find point xe outside [0,1]
                  !print '("initial guess...")'
                  k = 0
                  do
                    ta1 = exp(real(k))*ta; ta2 = ta1
                    xe = xn(ta2, x1, v1, Pe, zeta, tzeta)
                    !print '(i3,1x,es10.3)', k, xe
                    if ((xe <= O).or.(L <= xe)) exit ! xe outside [0,1]
                    if ((ta1 /= ta2).or.(ta2 >= infty)) exit ! avoid overflow in xn(t) with t too large
                    k = k+1
                  end do
                  ns = ns+k ! initial guess steps
                  if ((ta1 == ta2).and.(ta2 < infty)) then ! check if guessed t is reasonably small
                    ! second, reduce ta if xn(ta) is far outside [0,1]
                    if (k == 0) then
                      xa = xe
                      k = -1
                      do
                        xe = xa ! xa from previous iteration is outside [0,1]
                        ta1 = exp(real(k))*ta
                        xa = xn(ta1, x1, v1, Pe, zeta, tzeta)
                        !print '(i3,1x,es10.3)', k, xa
                        if (((O <= xa).and.(xa <= L)).or.(k < -23)) exit
                        k = k-1
                      end do
                    end if
                    ns = ns-k ! initial guess steps
                    ! update exit interface xa based on point xe outside
                    xa = O; if (xe >= L) xa = L
!                   Newton search
                    if (newtonsearch(x1, v1, Pe, zeta, tzeta, xa, exp(real(k))*ta, te, k, nsamax)) then ! success
                      xe = xa
                    else ! failed
                      !print '("failed")'
                      te = fail; xe = fail
                      print '("warning: failed (x1,v1,Pe,zeta,tzeta",/,5("  ",es10.3),")")', &
                        x1, v1, Pe, zeta, tzeta
                    end if
                    ns = ns+k ! Newton steps
                  else ! guess for te > tmax
                    te = infty; xe = xa
                  end if
                else ! trapped at stationary point
                  te = infty; xe = xs
                end if
              end if ! small-time solution, te >= infty
            end if ! pure diffusion
          end if ! already at exit position, x1 == xa
        end if ! advection and/or diffusion / pure advection, zeta == 0
        !print '("te = ",es10.3,", xe = ",f8.5)', te, xe
      end subroutine travelcell1d


!     Newton-Raphson/bisection algorithm
!     find t from eq. x(t) = xe, using underrelaxation between iterations and
!     bisectioning if Newton-iteration steps out of brackets, if there are any
      logical function newtonsearch(x1, v1, Pe, zeta, tzeta, xe, tg, te, k, kmax)
        implicit none
        real, intent(in) :: x1, v1 ! start point and velocity
        real, intent(in) :: Pe, zeta, tzeta
        real, intent(in) :: xe, tg ! exit interface, guess for exit time te
        integer, intent(out) :: k ! iteration counter
        integer, intent(in) :: kmax ! max. number iterations
        real, intent(out) :: te ! exit time
        real :: t0, t1, xs0, xs1, u  ! times, positions, velocity
        real, parameter :: rtol = 1.0E-5, atol = 1.0E-15, urf = 1.0/100 ! tolerances, underrelax. fact.
        real :: tb1, tb2, xb1, xb2 ! brackets
        logical :: bf ! bracketing state
        ! initialize
        !print '("Newton search...")'
        newtonsearch = .true.; t0 = tg; xs0 = xn(t0, x1, v1, Pe, zeta, tzeta)
        bf = .false.; tb1 = -1.0; tb2 = -1.0
        ! Newton or bisection search
        do k = 1,kmax
          if ((.not.bf).or.(mod(k,2) == 0)) then ! brackets do not exist -> Newton step
            u = un(t0,v1,Pe,zeta,tzeta)
            if (u /= 0.0) then
              t1 = t0 - (xs0 - xe) / u
            else
              t1 = t0 - sign(infty/2,xs0 - xe)
            end if
            t1 = max(0.0,t1) ! time must be positive
            t1 = (1 - urf)*t1 + urf*t0 ! underrelaxation
          end if
          if (bf) then ! bisection step
            if((mod(k,2) == 1).or.((t1<tb1).or.(tb2<t1))) t1 = (tb1+tb2)/2
          end if
          ! update x
          xs1 = xn(t1, x1, v1, Pe, zeta, tzeta)
          ! maintain brackets b1, b2 for root at xe
          if (bf) then ! tighten existing brackets
            if (sign(1.0,xb1-xe) /= sign(1.0,xs1-xe)) then
              tb2=t1; xb2=xs1 ! tighten right bracket on t-axis
            else
              tb1=t1; xb1=xs1 ! tighten left bracket on t-axis
            end if
          else ! establish brackets if possible
            if (sign(1.0,xs0-xe) /= sign(1.0,xs1-xe)) then
              if (t0 < t1) then
                tb1=t0; xb1=xs0; tb2=t1; xb2=xs1
              else
                tb1=t1; xb1=xs1; tb2=t0; xb2=xs0
              end if
              bf = .true. ! found brackets
              !print '("found brackets for root")'
            end if
          end if
          ! check convergence
          if (abs(t1-t0) <= abs(t1)*rtol + atol) exit
          ! complete step
          t0 = t1; xs0 = xs1
          !print '(es10.3)', t0
        end do
        newtonsearch = (.not.(k >= kmax))
        te = t1
      end function newtonsearch


!     diffusion step time [T]
      real function diffusiontime(d, vl, vr, Lc, Dm)
        implicit none
        integer, intent(in) :: d ! spatial dimensionality
        real, intent(in) :: vl(d), vr(d) ! grid cell velocities [L/T]
        real, intent(in) :: Lc(d) ! grid cell size [L]
        real, intent(in) :: Dm ! molecular diffusion coefficient [L^2/T]
        real :: A(d), v1(d), Pe(d), xe, ta(d)
        integer :: k ! loop counter
        real, parameter :: cells2travel = 2.0 ! mean numb. cells travelled in diffusiontime
        A = (vr-vl)/Lc; v1 = (vl+vr)/2;
        v1 = v1*Lc/Dm; Pe = A*Lc**2/Dm
        do k = 1,d
          xe = O; if (v1(k) > 0.0) xe = L; ! exit cell interface
          ta(k) = tea(0.5, v1(k), Pe(k), xe) ! travel time from cell center
          xe = 2*cells2travel * Lc(k)**2/Dm ! advection travel time factor
          if (ta(k) < infty / max(1.0,xe)) & ! ta may be infinity
            ta(k) = ta(k) * xe ! [-] -> [T]
        end do
        diffusiontime = min(minval(ta),minval(cells2travel * Lc**2/Dm))
      end function diffusiontime


!     advection-only exit time [-]
      real function tea(x1, v1, Pe, xe)
        implicit none
        real, intent(in) :: x1, v1, Pe, xe
        double precision :: a
        if (v1 == 0.0) then
          tea = infty
        elseif (Pe == 0.0) then
          tea = (xe - x1) / v1 ! limit of next expression for Pe -> 0 (zero veloc. grad.)
        else
           a = Pe*(xe - x1) / v1
          if (a <= -1.0) then ! handle stationary point
            tea = infty ! stationary point inside [0,1] with u = 0
          else
            tea = real(log(1.0 + a) / Pe) ! advection only (Pollock solution)
          end if
        end if
      end function tea


!     diffusion-dominated guess for exit time [-]
!     input condition: zeta not = 0
      real function tdg(x1, v1, zeta, tzeta, xe)
        implicit none
        real, intent(in) :: x1, v1, zeta, tzeta, xe
        double precision :: a, b, c
        ! guess from second-order polynomial approximation
        if (tzeta == 0) then
          if (v1**2 == 0.0) then
            tdg = (x1 - xe)**2 / (2*zeta**2)
          else
            a = sqrt(2*v1 * (xe - x1) + zeta**2)
            tdg = real(min((-a - zeta)**2, (a - zeta)**2))
            ! in case of truncation error due to zeta ~= to a, use alternative form
            if (tdg == 0.0) then
              b = exp(-zeta)*exp(-a); c = exp(a)*exp(-zeta)
              tdg = real(min(log(b)**2, log(c)**2))
            end if
            tdg = 2.0 * tdg / (4 * v1**2)
          end if
        else ! tzeta < 0
          tdg = (2*sqrt(-tzeta)*(xe-x1)) / (2*sqrt(-tzeta)*v1 + sqrt(2.0)*zeta)
        end if
      end function tdg


!     position evolution equation [-]
      real function xn(t, x1, v1, Pe, zeta, tzeta)
        implicit none
        real, intent(in) :: x1, v1, Pe, zeta, tzeta
        real, intent(inout) :: t ! limited t
        real :: erfttzeta, erftzeta
        if (Pe == 0.0) then
          xn = x1 + t*v1 + sqrt(2.0)*(sqrt(t-tzeta) - sqrt(-tzeta))*zeta
        elseif (Pe > 0.0) then
          ! if needed, limit t to avoid overflow in xn = ... exp(Pe*t)
          if (t >= 1) then
            if (Pe > log(infty)/2/t) t = log(infty)/2/Pe
          elseif (Pe >= 1) then
            if (t > log(infty)/2/Pe) t = log(infty)/2/Pe
          end if
          xn = x1 + v1/Pe*(exp(Pe*t)-1)
          ! for large negative tzeta, erf term drops out
          erfttzeta = erf(sqrt(Pe*(t-tzeta))); erftzeta = erf(sqrt(-Pe*tzeta))
          if (erfttzeta > erftzeta) then ! erf term is non-zero
            xn = xn + sqrt(pi/(2*Pe)) * exp(Pe*t)*zeta * &
              exp(-Pe*tzeta + log(erfttzeta - erftzeta))
          end if
        else ! Pe < 0
          ! if needed, limit t to avoid underflow in xn = ... exp(Pe*t) ...
          if (t >= 1) then
            if (-Pe > log(infty)/2/t) t = log(infty)/2/(-Pe)
          elseif (-Pe >= 1) then
            if (t > log(infty)/2/(-Pe)) t = log(infty)/2/(-Pe)
          end if
          xn = x1 + v1/Pe*(exp(Pe*t)-1) + sqrt(2/(-Pe))*zeta* &
            (dawson(sqrt(-Pe*(t-tzeta))) - exp(Pe*t)*dawson(sqrt(Pe*tzeta)))
        end if
      end function xn


!     velocity evolution equation [-]
      real function un(t, v1, Pe, zeta, tzeta)
        implicit none
        real, intent(in) :: v1, Pe, zeta, tzeta
        real, intent(inout) :: t ! limited t
        real :: erfttzeta, erftzeta
        if ((t == 0.0).and.(tzeta == 0.0)) then
          un = sign(infty,zeta)
        elseif (Pe == 0.0) then
          un = v1 + zeta/sqrt(2*(t-tzeta))
        elseif (Pe > 0.0) then
          if (log(infty)/2 < Pe*t) t = log(infty)/2/Pe ! avoid overflow in exp(Pe*t)
          un = zeta/sqrt(2*(t-tzeta)) + v1*exp(Pe*t)
          ! for large negative tzeta, erf term drops out
          erfttzeta = erf(sqrt(Pe*(t-tzeta))); erftzeta = erf(sqrt(-Pe*tzeta))
          if (erfttzeta > erftzeta) then ! erf term is non-zero
            un = un + sqrt(pi*Pe/2) * exp(Pe*t)*zeta * &
              exp(-Pe*tzeta + log(erfttzeta - erftzeta))
          end if
        else ! Pe < 0
          if (log(infty)/2 < -Pe*t) t = log(infty)/2/(-Pe) ! avoid overflow in exp(-Pe*t)
          un = v1*exp(Pe*t) + zeta*(1.0/sqrt(2*(t-tzeta)) + &
            sqrt(-2*Pe)*(exp(Pe*t)*dawson(sqrt(Pe*tzeta)) - dawson(sqrt(-Pe*(t-tzeta)))))
        end if
      end function un


!     Dawson integral function
!     copied from "Numerical Recipes in Fortran 77 Version 2.10"
      real function dawson(x)
        implicit none
        real, intent(in) :: x
        integer nmax
        real h,a1,a2,a3
        parameter (nmax=6,h=0.4,a1=2./3.,a2=0.4,a3=2./7.)
        integer i,n0
        real d1,d2,e1,e2,sum,x2,xp,xx
        integer, save :: init = 0
        real, save :: c(nmax)
!$omp   threadprivate(init, c)
        if (init.eq.0) then
          init=1
          do 11 i=1,nmax
            c(i)=exp(-((2.*float(i)-1.)*h)**2)
 11       continue
        endif
        if (abs(x).lt.0.2) then
          x2=x**2
          dawson=x*(1.-a1*x2*(1.-a2*x2*(1.-a3*x2)))
        else
          xx=abs(x)
          n0=2*nint(0.5*xx/h)
          xp=xx-float(n0)*h
          e1=exp(2.*xp*h)
          e2=e1**2
          d1=float(n0+1)
          d2=d1-2.
          sum=0.
          do 12 i=1,nmax
            sum=sum+c(i)*(e1/d1+1./(d2*e1))
            d1=d1+2.
            d2=d2-2.
            e1=e2*e1
 12       continue
          dawson=0.5641895835*sign(exp(-xp**2),x)*sum
        endif
      end function dawson


!     standard normally distributed random deviates
!     copied from "Numerical Recipes in Fortran 77 Version 2.10"
      real function gasdev(idum)
        implicit none
        integer, intent(inout) :: idum
        real fac,rsq,v1,v2
        integer, save :: iset = 0
        real, save :: gset
!$omp   threadprivate(iset, gset)
        if (idum.lt.0) iset=0
        if (iset.eq.0) then
 1        v1=2.*ran1(idum)-1.
          v2=2.*ran1(idum)-1.
          rsq=v1**2+v2**2
          if(rsq.ge.1..or.rsq.eq.0.)goto 1
          fac=sqrt(-2.*log(rsq)/rsq)
          gset=v1*fac
          gasdev=v2*fac
          iset=1
        else
          gasdev=gset
          iset=0
        endif
      end function gasdev


!     random deviate, minimal standard plus shuffle
!     returns a uniform random deviate between 0 and 1 (exclusive of the endpoint values)
!     copied from "Numerical Recipes in Fortran 77 Version 2.10"
      real function ran1(idum)
        implicit none
        integer, intent(inout) :: idum
        integer ia,im,iq,ir,ntab,ndiv
        real am,epsr,rnmx
        parameter (ia=16807,im=2147483647,am=1./im,iq=127773,ir=2836, &
                   ntab=32,ndiv=1+(im-1)/ntab,epsr=1.2e-7,rnmx=1.-epsr)
        integer j,k
        integer, save :: iv(ntab) = 0, iy = 0
!$omp   threadprivate(iv, iy)
        if (idum.le.0.or.iy.eq.0) then
          idum=max(-idum,1)
          do 11 j=ntab+8,1,-1
            k=idum/iq
            idum=ia*(idum-k*iq)-ir*k
            if (idum.lt.0) idum=idum+im
            if (j.le.ntab) iv(j)=idum
 11       continue
          iy=iv(1)
        endif
        k=idum/iq
        idum=ia*(idum-k*iq)-ir*k
        if (idum.lt.0) idum=idum+im
        j=1+iy/ndiv
        iy=iv(j)
        iv(j)=idum
        ran1=min(am*iy,rnmx)
      end function ran1

      end module pollock_lib


!     Library with the following functions and subroutines:
!      pickrandomorigininvoid: randomly pick cell position in void space
!      pickrandomoriginatinflow: randomly pick flux-weighted initial cell on inflow plane
!      startpointnextsection: determine starting point on plane for next s.l. section
!      openandwriteplumefiles: open new and append existing plume files
!      rani: generate random integer deviates in the range 1..n
!      ran1: generate random numbers that are uniformly distributed in the interval between 0 and 1
!
!     Daniel W. Meyer, February 2017
!     Christoph Schultheiss, April 2017
!     Institute of Fluid Dynamics
!     ETH Zurich
!
      module streamlinesnt3d_lib
     
        implicit none

      contains

!     randomly pick initial cell in void space
      subroutine pickrandomorigininvoid(g, n, idum, ic)
        use pollock_lib, only: getg, getv
        implicit none
!       interface variables
        integer, intent(in) :: n(3) ! number of grid cells
        logical, intent(in) :: g(product(n)) ! geometry data [true = solid, false = void]
        integer, intent(inout) :: idum ! seed for random number generator
        integer, intent(out) :: ic(3) ! streamline start cell index
!       implementation
        do ! loop until a void cell was randomly chosen
          ic = (/rani(n(1),idum), rani(n(2),idum), rani(n(3),idum)/)
          if (.not.getg(3, n, g, ic)) exit
        end do
      end subroutine pickrandomorigininvoid


!     randomly pick flux-weighted initial cell on inflow plane
      subroutine pickrandomoriginatinflow(v, g, n, idum, ic)
        use pollock_lib, only: getg, getv
        implicit none
!       interface variables
        integer, intent(in) :: n(3) ! number of grid cells
        logical, intent(in) :: g(product(n)) ! geometry data [true = solid, false = void]
        real, intent(in) :: v(product(n+1)) ! velocity [L/T]
        integer, intent(inout) :: idum ! seed for random number generator
        integer, intent(out) :: ic(3) ! streamline start cell index
!       other variables
        real :: vimax ! maximum v1 velocity on inflow plane
        integer :: i2, i3 ! cell indices
!       implementation
        ! maximum velocity on inflow plane
        vimax = 0.0
        do i2 = 1, n(2)
          do i3 = 1, n(3)
            if (.not.getg(3, n, g, (/1,i2,i3/))) then ! only void or inflow cells
              vimax = max(vimax, getv(3, n+1, v, (/1,i2,i3/), 1))
            end if
          end do
        end do
        ! pick inflow plane cell via acceptance/rejection sampling
        do ! probability of cell is proportional to local inflow velocity
          ic = (/1, rani(n(2),idum), rani(n(3),idum)/)
          if (.not.getg(3, n, g, ic)) then ! only void or inflow cells
            if (ran1(idum) <= getv(3, n+1, v, ic, 1)/vimax) exit
          end if
        end do
      end subroutine pickrandomoriginatinflow


!     determine starting point on plane for next s.l. section
      subroutine startpointnextsection(v, g, n, vi, rflag, ic, alpha, dflag)
        use pollock_lib, only: getg, getv
        implicit none
!       interface variables
        integer, intent(in) :: n(3) ! number of grid cells
        logical, intent(in) :: g(product(n)) ! geometry data [true = solid, false = void]
        real, intent(in) :: v(product(n+1)) ! velocity [L/T]
        real, intent(in) :: vi(3) ! target velocity [L/T]
        logical, intent(in) :: rflag ! rotate domain when sampling streamline sections?
        integer, intent(out) :: ic(3) ! streamline start cell index
        real, intent(inout) :: alpha ! rotation angle
        logical, intent(in):: dflag ! direction flag (true/false = down/upstream)
!       other variables
        real :: dv, dvmin, vr(3) ! velocity differences, current velocity [L/T]
        integer :: i2, i3 ! loop counters
        integer :: i1 ! index of plane for next s.l. section
!       implementation
        dvmin = -1.0
        i1 = merge(1,n(1),dflag)
        ! find velocity that is closest to vi on x2-x3-plane at i1
        do i3 = 1, n(3)
          do i2 = 1, n(2)
            if (.not.getg(3, n, g, (/i1,i2,i3/))) then
              vr = (/getv(3, n+1, v, (/i1,i2,i3/), 1), &
                     getv(3, n+1, v, (/i1,i2,i3/), 2), &
                     getv(3, n+1, v, (/i1,i2,i3/), 3)/) ! velocity in cell (i1,i2,i3)
              ! calculate distance between vr and vi
              if (rflag) then ! vr is rotated around x1-axis into vi - x1-axis - plane
                dv = sqrt((vr(1)-vi(1))**2 + (sqrt(vr(2)**2+vr(3)**2) - sqrt(vi(2)**2+vi(3)**2))**2)
              else ! no rotation of vr and domain
                dv = sqrt((vr(1)-vi(1))**2 + (vr(2)-vi(2))**2 + (vr(3)-vi(3))**2)
              end if
              ! check if current vr is closer to vi
              if ((dv < dvmin).or.(dvmin < 0.0)) then
                dvmin = dv ! keep track of shortest distance
                ic = (/i1,i2,i3/) ! store grid cell index
                if (rflag) alpha = atan2(vr(3),vr(2)) - atan2(vi(3),vi(2)) ! rotation angle
              end if
            end if
          end do
        end do
      end subroutine startpointnextsection


!     in case plumeno (and d) is present, write data d to plume file identified by case and plume no.
!     otherwise, by searching for an inexistent file plume_x_1.txt, where x is the caseno,
!     identify a new caseno
      subroutine openandwriteplumefiles(caseno, plumeno, d)
        implicit none
!       interface variables
        integer, intent(inout) :: caseno ! new or previously identified case number
        integer, optional :: plumeno ! number of plume
        double precision, optional :: d(:) ! data
!       other variables
        character*(256) :: Ofile ! name of plume file
        logical :: existflag ! used to check if plume file exists already
!       implementation
        if (present(plumeno).and.present(d)) then ! append data d to existing plume file
          write(Ofile,fmt='("plume_",i0,"_",i0,".txt")') caseno, plumeno
          open(unit=8,file=Ofile,form='formatted',action='write',position='append')
          write(unit=8,fmt='(999(es13.6,1x))') d
          close(unit=8)
        else ! identify new case number
          caseno = 1
          do
            write(Ofile,fmt='("plume_",i0,"_",i0,".txt")') caseno, 1
            inquire(file=Ofile, exist=existflag) ! check if plume file exists
            if (.not.existflag) exit
            caseno = caseno + 1 ! next case
          end do
          ! reserve case number by creating first plume file
          open(unit=8,file=Ofile,form='formatted',action='write'); close(unit=8)
        end if
      end subroutine openandwriteplumefiles


!     generate random integer deviates in the range 1..n
      integer function rani(n, idum)
        implicit none
        integer, intent(in) :: n
        integer, intent(inout) :: idum
        rani = ceiling(real(n)*ran1(idum))
        if (rani .eq. 0) rani = 1 ! this is very, very unlikely
      end function rani


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

      end module streamlinesnt3d_lib

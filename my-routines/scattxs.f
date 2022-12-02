      subroutine scattxs(nmaxp, wp, itrans, theta, !inp
     1                   skn)                      !out
c
c     Compute the exact Compton cross section taking into account
c     the Klein-Nishina corrections at high energies, and the 
c     relativistic corrections at high temperatures (T>1.e8 K)
c
!$    use omp_lib
      use constants
      implicit none
      integer nmaxp, np, itrans, iz
      real*8 xloc
      real*8 skn(nmaxp,itrans), ixloc, wp(nmaxp), theta(itrans)
c
c
      do iz = 1, itrans
!$omp parallel
!$omp& shared(iz,nmaxp,wp,theta,skn)
!$omp& private(np,xloc,ixloc)
!$omp do
         do np=1, nmaxp
            xloc = 3.913894d-6*wp(np)
            ixloc = 1.d0/xloc
            if (theta(iz).lt.0.0169d0) then    ! T < 1.e8 K
               if (xloc.lt.0.002d0) then
                  skn(np,iz) = sigma_t
               else
                  skn(np,iz) = 4.9875d-25 * ((1.d0-4.d0 * ixloc-8.d0 *
     1                     ixloc**2.0)* dlog(1.d0 + xloc) + 0.5d0 +
     2                     8.d0*ixloc - 0.5d0/((1.d0+xloc)**2.0))*ixloc
               endif
            else 
               call crsexact(theta(iz), wp(np)/511.d3, skn(np,iz))
            endif
         enddo
!$omp end do
!$omp end parallel
      enddo
c
      return
      end

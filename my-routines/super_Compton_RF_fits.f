!-------------------------------------------------------------------------------------------------
      subroutine super_Compton_RF_fits(itrans, temps, theta, 
     &                                 nmaxp, wp, df, skn, 
     &                                 nksmps, limit)
c
c     This routine writes a file with the super redistribution function (SRF) for Compton
c     scatterting. For a given gas temperature T and final photon energy Ef, the SRF is
c     defined for a set of initial photon energies Ei as:
c
c         SRF(T,Ef,Ei) = IRF(Ef,Ei)/N(Ei)*skn(Ei)*dEi/Ei
c
c     where IRF(Ef,Ei) is in fact the inverse redistribution function for the Compton scattering
c     of a photon from initial energy Ei to final energy Ef; N(Ei) is the normalization to
c     ensure photon number conservation; and skn(Ei) is the Klein-Nishina cross section.
c     This routine implements the exact Compton RF from Madej et al. (2017).
c     The SRF contains all the information needed for the convolution of a given spectrum
c     to account for the Compton scattering at the given temperature.
c     Only significant values of the SRF are actually written, i.e., when RF > limit.
c
c     Input arguments:
c         itrans: Total number of temperatures
c         temps: Array (itrans) of temperatures in Kelvin
c         theta: Array (itrans) of temperatures in kT/mec2
c         nmaxp: Total number of energy points
c         wp: Array (nmaxp) of energies in eV
c         df: Array (nmaxp) of delta energies in eV
c         skn: Array (nmaxp) of Kein-Nishina cross sections
c         nksmps: Total number of angles for averaging
c
c         limit: Fractional limit below which the SRF is set to 0
c
c     Output arguments: 
c         None
c
c     Internal aruments:
c         smit: Array (nksmps) Legendre ordinates (angles)
c         agt: Array (nksmps) weights
c
c     Requires:
c         probab.f: Routine for the RF calculation
c
!$    use omp_lib
      use constants
      use fits_writing
      implicit none
      integer itrans, nmaxp, nksmps, iz, np, jj
      integer j_pmax, jlow, jhigh
      logical GOUP, GODOWN
      real*8 problim
      real*8 theta(itrans), wp(nmaxp), df(nmaxp), skn(nmaxp,itrans)
      real*8 smit(nksmps), agt(nksmps)
      real*8 check
      real*8, allocatable :: prob(:,:) !real*8 prob(nmaxp,nmaxp)
      real*8 ecen, temps(itrans)
      real*8 limit
      integer n, curr_ind_s, curr_ind_e
      integer, allocatable :: fSInd(:), fLen(:)

c     These arrays need to be allocatable to avoid segfaults at large nmaxp
      allocate( prob(nmaxp,nmaxp) )   ! Prob(initial, final)
      allocate( fSInd(nmaxp*itrans), fLen(nmaxp*itrans) )

c     Check1 is to ensure photon number is conserved in scatterings

      call set_filename('table.fits') !name of the fits file
      n = 1 ! column number of the fits file

c     Get the Gaussian quadratures for angular integration
      call gaulegf(-1.d0, 1.d0, smit, agt, nksmps)

c     Set up a new fits file
      call setup_new_file(nmaxp, itrans, 0,
     &                    wp, df, temps, smit, agt,
     &                    skn, limit, .TRUE., nksmps) !create the fits file

      do iz = 1, itrans
         prob = 0.d0
c         temp = theta(iz)*ikbol*mec2          ! temperature in K
!$omp parallel 
!$omp& shared(iz, wp, prob, nmaxp, nksmps,
!$omp&         smit, agt, theta)
!$omp& private(np, ecen, jj, problim, j_pmax,
!$omp&         jlow, jhigh, GOUP, GODOWN)
!$omp do
         ! Initial energy
         do np = 1, nmaxp
            ecen = wp(np)/mec2
            call probab(theta(iz),ecen,ecen,nksmps,smit,
     &                  agt,prob(np,np))
            problim = prob(np,np)*limit
            j_pmax = np
            jlow = np-1
            jhigh = np+1
            GOUP = .TRUE.
            GODOWN = .TRUE.
            ! Final energy
            do while (GOUP .or. GODOWN)
               if(jlow.LT.1)       GODOWN = .FALSE.
               if(jhigh.GT.nmaxp)  GOUP = .FALSE.
               if(GODOWN)then
                  call probab(theta(iz),wp(jlow)/mec2,ecen,nksmps,smit,
     &                     agt,prob(np,jlow))
                  if( prob(np,jlow) .LT. problim )then
                     GODOWN = .FALSE.
                  elseif ( prob(np,jlow) .GT. prob(np, j_pmax) ) then
                     j_pmax = jlow
                     problim = prob(np,jlow)*limit
                  endif
                  jlow = jlow -1
               endif
               if(GOUP)then
                  call probab(theta(iz),wp(jhigh)/mec2,ecen,nksmps,smit,
     &                     agt,prob(np,jhigh))
                  if( prob(np,jhigh) .LT. problim )then
                     GOUP = .FALSE.
                  elseif ( prob(np,jhigh) .GT. prob(np, j_pmax) ) then
                     j_pmax = jhigh
                     problim = prob(np,jhigh)*limit
                  endif
                  jhigh = jhigh + 1
               endif
            enddo
         enddo
!$omp end do
!$omp end parallel


         curr_ind_s = (iz-1)*nmaxp
         curr_ind_e = iz*nmaxp
         
         call srf_nonlimit( prob(:,:), nmaxp, limit, 
     &                      fSInd(curr_ind_s+1 : curr_ind_e), 
     &                      fLen(curr_ind_s+1 : curr_ind_e) )
         
         do jj=1, nmaxp
             prob(1:fSInd(curr_ind_s+jj)-1,jj)=0.0
             prob(fSInd(curr_ind_s+jj)+fLen(curr_ind_s+jj):nmaxp,jj)=0.0
         enddo


         do np=1,nmaxp ! Initial
            check = 0.0
            do jj=1,nmaxp ! Final 
               check = check + df(jj) * prob(np, jj)
            enddo
            prob(np,:) = prob(np,:) *skn(np,iz)*df(np)/wp(np)/check
         enddo
         
c     Write the iSRF of this temperature to the fits file
         do jj = 1, nmaxp
            call write_SRFs(n, 1, prob(:,jj), 
     &                      fSInd(curr_ind_s+jj), 
     &                      fLen(curr_ind_s+jj) )
            n = n + 1  
         enddo
         print *,iz
      enddo   
 

      call write_SRF_pointers(1, fSInd, fLen)

c     The FITS file must always be closed before exiting the program. 
      call close_and_save_fits()
      deallocate( prob, fSInd, fLen )
      return
      end subroutine

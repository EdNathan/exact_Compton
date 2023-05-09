!-------------------------------------------------------------------------------------------------      
      subroutine super_Compton_RF_fits_angle(itrans, temps, theta,
     &                                       nmaxp, wp, df, skn, mgi,
     &                                       smit, agt, nksmps, limit) 
c     This routine writes a file with the super redistribution function (SRF) for Compton
c     scatterting. For a given gas temperature T and final photon energy Ef, the SRF is
c     defined for a set of initial photon energies Ei as:
c
c         SRF(T,f,i) = IRF(f,i)/N(i)*skn(Ei)*dEi/Ei
c
c     where IRF(f,i) is in fact the inverse redistribution function for the Compton scattering
c     of a photon from initial energy Ei and initial angle mui, to final energy Ef and  
c     final angle muf.  N(i) a normalization for an initial energy and angle to
c     ensure photon number conservation; and skn(Ei) is the Klein-Nishina cross section.
c     This routine implements the exact Compton RF from Madej et al. (2017).
c     The SRF contains all the information needed for the convolution of a given spectrum
c     to account for the Compton scattering at the given temperature.
c     Only significant values of the SRF are actually written, i.e., when RF > limit.
c
c     Input arguments: 
c         itrans: Total number of temperatures
c         theta: Array (itrans) of temperatures in kT/mec2
c         nmaxp: Total number of energy points
c         wp: Array (nmaxp) of energies in eV
c         df: Array (nmaxp) of delta energies in eV
c         skn: Array (nmaxp) of Kein-Nishina cross sections
c         mgi: Total number of angles
c         smit: Array (mgi) Legendre ordinates (angles)
c         agt: Array (mgi) weights 
c         nksmps: Number of sample points for azimuthal integration
c         limit: Fractional limit below which the SRF is set to 0
c
c     Output arguments: 
c         None
c         
c     Requires:
c         probab.f: Routine for the RF calculation
c
c
!$    use omp_lib
      use constants
      use fits_writing
      implicit none
      integer itrans, nmaxp, mgi, iz, np, jj, qq
      integer nksmps, inang, outang
      real*8 A, B
      real*8 theta(itrans), wp(nmaxp), df(nmaxp), skn(nmaxp,itrans),x
      real*8 smit(mgi), agt(mgi)
      real*8 kords(nksmps), kweights(nksmps)
      real*8 check
      real*8 profil, temps(itrans)
      real*8 limit
      real*8, target :: srf(mgi,nmaxp,mgi,nmaxp)
      real*8, pointer :: flatsrf(:,:)
      integer, target :: fSInd(mgi*nmaxp*itrans), fLen(mgi*nmaxp*itrans)
      integer, pointer :: flatfs(:,:), flatfl(:,:)
      integer n, curr_ind_s, curr_ind_e

c     flatsrf is a pointer to allow accessing to the response function as a matrix
c     flatsrf( init, final )
      flatsrf(1:mgi*nmaxp, 1:mgi*nmaxp) => srf

c     flatfs and flatfl are pointers to fSInd and fLen to make it easier to deal with the running temperature
      flatfs(1:mgi*nmaxp, 1:itrans) => fSInd
      flatfl(1:mgi*nmaxp, 1:itrans) => fLen

      call set_filename('angle.fits') !name of the fits file
      n = 1 ! column number of the fits file
      
c     Set up a new fits file
      call setup_new_file(nmaxp, itrans, mgi,
     &                    wp, temps, smit, agt,
     &                    skn, limit, .FALSE., nksmps) !create the fits file

      do iz = 1, itrans
         x=1/theta(iz)
        
         srf = 0.d0

c        srf( init_ang, init_en, final_ang, final_en  )
         do inang=1, mgi
!$omp parallel
!$omp& shared(iz,wp,nmaxp,mgi,smit,x,srf,kords,kweights)
!$omp& private(np,jj,qq,A,B)
!$omp do
            do np = 1, nmaxp ! initial energy
               do outang=inang, mgi
                  A = sqrt((1-smit(inang)**2)*(1-smit(outang)**2))
                  B = smit(inang) * smit(outang)
                  do jj=1,nmaxp ! final energy
c                    Do azimuthal integration
                     do qq=1,nksmps
                        srf(inang,np,outang,jj) = 
     1                     srf(inang,np,outang,jj) 
     2                     + profil(1,wp(np)/mec2,wp(jj)/mec2, 
     3                              A*kords(qq)+B,x)
     4                     / sqrt(1 - kords(qq)**2)
     5                     * kweights(qq)
                     enddo ! qq
                     srf(inang,np,outang,jj)=srf(inang,np,outang,jj)/pi
                  enddo ! jj
               enddo ! outang
            enddo ! np
!$omp end do
!$omp end parallel
         enddo ! inang
         do inang=2, mgi
            do outang=1,inang-1
               do np=1, nmaxp
                  do jj=1, nmaxp
c                    Symmetry:
                     srf(inang,np,outang,jj)=srf(outang,np,inang,jj)
                  enddo
               enddo
            enddo
         enddo


c        srf( init_ang, init_en, final_ang, final_en  )
c        flatsrf( init_pair, final_pair)
         call srf_nonlimit( flatsrf(:,:), nmaxp*mgi, limit, 
     &                      flatfs(:,iz), 
     &                      flatfl(:,iz)) 

c        Go through each 'final' energy/angle row, and set the areas to 0 when below limit
         do jj=1, nmaxp*mgi
            flatsrf( 1 :flatfs(jj,iz)-1 , jj) = 0.0
            flatsrf( flatfs(jj,iz) + flatfl(jj,iz) : nmaxp*mgi, jj) = 0.0
         enddo

c        Check to ensure photon number is conserved in scatterings.  
c        0.5 * Int_{-1}^{+1} Int_0^\inf R[init_mu, init_en, fin_mu, fin_en] d{fin_en} d{fin_my} = 1

c        srf( init_ang, init_en, final_ang, final_en  )
         do inang=1,mgi ! inital angle
            do np=1,nmaxp ! initial energy
               check = 0.0
               do outang=1,mgi ! final angle
                  do jj=1,nmaxp ! final energy
                     check = check  
     1                 + df(jj) * srf(inang,np,outang,jj) * agt(outang)
                  enddo
               enddo
               check = check / 2
               srf(inang,np,:,:)=srf(inang,np,:,:)
     1                           * skn(np,iz)*df(np)/wp(np)/check
            enddo
         enddo
      
         do jj=1,nmaxp ! final energy
            do outang=1, mgi
               do np = 1, nmaxp ! initial energy
                   do inang=1, mgi
                     write(70,*) srf(inang,np,outang,jj) / skn(np,iz)
     &                             / df(np)  * wp(np) 
                  enddo
               enddo
            enddo
         enddo


c     Write the iSRF of this temperature to the fits file
         do jj = 1, nmaxp*mgi
            call write_SRFs(n, flatsrf(:,jj), 
     &                      flatfs(jj,iz), 
     &                      flatfl(jj,iz) )
            n = n + 1  
         enddo
         print *,iz
      enddo      

      call write_SRF_pointers(fSInd, fLen)
      
c     The FITS file must always be closed before exiting the program. 
      call close_and_save_fits()

      return
      end subroutine

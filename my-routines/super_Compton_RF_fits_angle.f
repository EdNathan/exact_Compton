!-------------------------------------------------------------------------------------------------      
      subroutine super_Compton_RF_fits_angle(itrans, temps, theta,
     &                                       nmaxp, wp, df, skn, mgi,
     &                                       smit, agt, limit) 
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
c         theta: Array (itrans) of temperatures in kT/mec2
c         nmaxp: Total number of energy points
c         wp: Array (nmaxp) of energies in eV
c         df: Array (nmaxp) of delta energies in eV
c         skn: Array (nmaxp) of Kein-Nishina cross sections
c         mgi: Total number of angles 
c         smit: Array (mgi) Legendre ordinates (angles)
c         agt: Array (mgi) weights 
c
c     Output arguments: 
c         None
c         
c     Requires:
c         probab.f: Routine for the RF calculation
c
!$    use omp_lib
      use constants
      use fits_writing
      implicit none
      integer itrans, nmaxp, mgi, iz, np, jj, qq
      real*8 theta(itrans), wp(nmaxp), df(nmaxp), skn(nmaxp,itrans),x
      real*8 smit(mgi), agt(mgi)
      real*8 check1(nmaxp)
      real*8 prob(nmaxp,nmaxp), srf(mgi,nmaxp,nmaxp)
      real*8 ecen, profil, temps(itrans)
      real*8 limit
      integer fSInd(mgi, nmaxp*itrans), fLen(mgi, nmaxp*itrans)
      integer n, curr_ind_s, curr_ind_e

c     Check1 is to ensure photon number is conserved in scatterings

      call set_filename('angle.fits') !name of the fits file
      n = 1 ! column number of the fits file

c     Set up a new fits file
      call setup_new_file(nmaxp, itrans, mgi,
     &                    wp, temps, smit, agt,
     &                    skn, limit, .FALSE.) !create the fits file

      do iz = 1, itrans
         x=1/theta(iz)

         check1 = 0.d0
        ! prob = 0.d0
         srf = 0.d0
c
        ! temp = theta(iz)*ikbol*mec2          ! temperature in K
         do np = 1, nmaxp
            ecen = wp(np)
!$omp parallel
!$omp& shared(iz,wp,prob,nmaxp,np,mgi,smit,
!$omp& agt)
!$omp& private(jj)
!$omp do
            do jj=1,nmaxp
               do qq=1,mgi
                  srf(qq,jj,np)=profil(1,wp(np)/mec2,wp(jj)/mec2,
     &             smit(qq),x)
               !   prob(jj,np)=prob(jj,np)+agt(qq)*srf(qq,jj,np)
               enddo
            enddo
!$omp end do
!$omp end parallel
         enddo
         
         curr_ind_s = (iz-1)*nmaxp
         curr_ind_e = iz*nmaxp
         do qq = 1, mgi
             call srf_nonlimit( srf(qq,:,:)*agt(qq), df, nmaxp, limit, 
     &                          fSInd(qq, curr_ind_s+1 : curr_ind_e), 
     &                          fLen(qq, curr_ind_s+1 : curr_ind_e), 
     &                          check1 )
         enddo

         do np=1,nmaxp
            srf(:,:,np)=srf(:,:,np)*skn(np,iz)*df(np)/wp(np)/check1(np)
         enddo
         
c     Write the iSRF of this temperature to the fits file
         do jj = 1, nmaxp
            call write_SRFs(n, srf(:,jj,:), 
     &                      fSInd(:,curr_ind_s+jj), 
     &                      fLen(:, curr_ind_s+jj) )
            n = n + 1  
         enddo
         print *,iz
      enddo      

      call write_SRF_pointers(fSInd, fLen)
      
c     The FITS file must always be closed before exiting the program. 
      call close_and_save_fits()

      return
      end subroutine

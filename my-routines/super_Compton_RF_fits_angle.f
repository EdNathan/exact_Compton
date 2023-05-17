!-------------------------------------------------------------------------------------------------      
      subroutine super_Compton_RF_fits_angle(itrans, temps, theta,
     &                                       nmaxp, wp, df, skn, mgi,
     &                                       smit, agt, nksmps, limit) 
c     This routine writes a file with the super redistribution function (SRF) for Compton
c     scatterting. For a given gas temperature T and final photon energy Ef, the SRF is
c     defined for a set of initial photon energies Ei as:
c
c         SRF(T,f,i) = IRF(f,i)/N(i)*skn(Ei)*dmui*dEi/Ei
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
      integer nksmps, inang, outang, angflip
      real*8 A, B
      real*8 theta(itrans), wp(nmaxp), df(nmaxp), skn(nmaxp,itrans),x
      real*8 smit(mgi), agt(mgi)
      real*8 kords(nksmps), kweights(nksmps)
      real*8 check, factor
      real*8 profil, temps(itrans)
      real*8 limit
      real*8 srfval
      real*8, target :: srfe(mgi,nmaxp,mgi,nmaxp), 
     &                  srfo(mgi,nmaxp,mgi,nmaxp)
      real*8, pointer :: flatsrfe(:,:), flatsrfo(:,:)
      integer, target :: fSInde(mgi*nmaxp*itrans), 
     &                   fLene(mgi*nmaxp*itrans)
      integer, target :: fSIndo(mgi*nmaxp*itrans), 
     &                   fLeno(mgi*nmaxp*itrans)
      integer, pointer :: flatfse(:,:), flatfle(:,:)
      integer, pointer :: flatfso(:,:), flatflo(:,:)
      integer n, curr_ind_s, curr_ind_e

c     flatsrf is a pointer to allow accessing to the response function as a matrix
c     flatsrf( init, final )
      flatsrfe(1:mgi*nmaxp, 1:mgi*nmaxp) => srfe
      flatsrfo(1:mgi*nmaxp, 1:mgi*nmaxp) => srfo

c     flatfs and flatfl are pointers to fSInd and fLen to make it easier to deal with the running temperature
      flatfse(1:mgi*nmaxp, 1:itrans) => fSInde
      flatfle(1:mgi*nmaxp, 1:itrans) => fLene

      flatfso(1:mgi*nmaxp, 1:itrans) => fSIndo
      flatflo(1:mgi*nmaxp, 1:itrans) => fLeno

      call set_filename('angle.fits') !name of the fits file
      n = 1 ! column number of the fits file
 
c     Set up a new fits file
      call setup_new_file(nmaxp, itrans, mgi,
     &                    wp, df, temps, smit, agt,
     &                    skn, limit, .FALSE., nksmps) !create the fits file

      call gaulegf(-1.d0, 1.d0, kords, kweights, nksmps)

      do iz = 1, itrans
         x=1/theta(iz)
        
         srfe = 0.d0
         srfo = 0.d0

c        srf( init_ang, init_en, final_ang, final_en  )
         do inang=1, mgi
            do angflip=-1,1,2
!$omp parallel
!$omp& shared(iz,wp,nmaxp,mgi,smit,x,srfe,srfo,kords,kweights,angflip)
!$omp& private(np,jj,qq,A,B,srfval)
!$omp do
               do np = 1, nmaxp ! initial energy
                  do outang=inang, mgi
                     A = sqrt((1-smit(inang)**2)*(1-smit(outang)**2))
                     B = smit(inang) * smit(outang) * angflip
                     do jj=1,nmaxp ! final energy
c                      Do azimuthal integration
                        srfval = 0.0
                        do qq=1,nksmps
                           srfval = srfval 
     1                        + profil(1,wp(np)/mec2,wp(jj)/mec2, 
     2                                 A*kords(qq)+B,x)
     3                        / sqrt(1 - kords(qq)**2)
     4                        * kweights(qq)
                        enddo ! qq
                        srfval = srfval / pi
                        srfe(inang,np,outang,jj) = 
     1                     srfe(inang,np,outang,jj) + srfval
                        srfo(inang,np,outang,jj) = 
     1                     srfo(inang,np,outang,jj) + angflip*srfval
                     enddo ! jj
                  enddo ! outang
               enddo ! np
!$omp end do
!$omp end parallel
            enddo ! angflip
         enddo ! inang
         do inang=2, mgi
            do outang=1,inang-1
               do np=1, nmaxp
                  do jj=1, nmaxp
c                    Symmetry:
                     srfe(inang,np,outang,jj)=srfe(outang,np,inang,jj)
                     srfo(inang,np,outang,jj)=srfo(outang,np,inang,jj)
                  enddo
               enddo
            enddo
         enddo


c        srf( init_ang, init_en, final_ang, final_en  )
c        flatsrf( init_pair, final_pair)
         call srf_nonlimit( flatsrfe(:,:), nmaxp*mgi, limit, 
     &                      flatfse(:,iz), 
     &                      flatfle(:,iz)) 

         call srf_nonlimit( flatsrfo(:,:), nmaxp*mgi, limit, 
     &                      flatfso(:,iz), 
     &                      flatflo(:,iz)) 


c        Go through each 'final' energy/angle row, and set the areas to 0 when below limit
         do jj=1, nmaxp*mgi
            flatsrfe(1:flatfse(jj,iz)-1,jj) = 0.0
            flatsrfe(flatfse(jj,iz)+flatfle(jj,iz):nmaxp*mgi,jj) = 0.0

            flatsrfo(1:flatfso(jj,iz)-1,jj) = 0.0
            flatsrfo(flatfso(jj,iz)+flatflo(jj,iz):nmaxp*mgi,jj) = 0.0
         enddo

c        Check to ensure photon number is conserved in scatterings.  
c        0.5 * Int_{-1}^{+1} Int_0^\inf R[init_mu, init_en, fin_mu, fin_en] d{fin_en} d{fin_my} = 1
c        Which becomes:
c        Int_{0}^{+1} Int_0^\inf Re[init_mu, init_en, fin_mu, fin_en] d{fin_en} d{fin_my} = 1

c        srf( init_ang, init_en, final_ang, final_en  )
         do inang=1,mgi ! inital angle
            do np=1,nmaxp ! initial energy
               check = 0.0
               do outang=1,mgi ! final angle
                  do jj=1,nmaxp ! final energy
                     check = check  
     1                  + srfe(inang,np,outang,jj)*agt(outang)*df(jj)
                  enddo
               enddo
               check = check / 2
               factor = skn(np,iz)*df(np)*agt(inang)/wp(np)/check
               srfe(inang,np,:,:)=srfe(inang,np,:,:) * factor
               srfo(inang,np,:,:)=srfo(inang,np,:,:) * factor
            enddo
         enddo

c     Write the iSRF of this temperature to the fits file
         do jj = 1, nmaxp*mgi
            call write_SRFs(n, 1, flatsrfe(:,jj), 
     &                            flatfse(jj,iz), 
     &                            flatfle(jj,iz) )

            call write_SRFs(n, 2, flatsrfo(:,jj), 
     &                            flatfso(jj,iz), 
     &                            flatflo(jj,iz) )
            n = n + 1  
         enddo
         print *,iz
      enddo      

      call write_SRF_pointers(1, fSInde, fLene)
      call write_SRF_pointers(2, fSIndo, fLeno)
      
c     The FITS file must always be closed before exiting the program. 
      call close_and_save_fits()

      return
      end subroutine

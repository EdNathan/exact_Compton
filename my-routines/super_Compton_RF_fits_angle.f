!-------------------------------------------------------------------------------------------------      
      subroutine super_Compton_RF_fits_angle(itrans, theta, nmaxp, wp, 
     & df,skn,  mgi, smit, agt) 
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
      implicit none
      include 'omp_lib.h'
      integer itrans, nmaxp, mgi, iz, np,qq, jj, kk
      real*8 theta(itrans), wp(nmaxp), df(nmaxp), skn(nmaxp,itrans),x
      real*8 smit(mgi), agt(mgi)
      real*8 check1(nmaxp), limit, pmax(nmaxp)
      real*8 prob(nmaxp,nmaxp), srf(mgi,nmaxp,nmaxp)
      real*8 mec2, ecen, isp, ikbol,profil, temp(itrans)
      integer point(nmaxp), indices(nmaxp,nmaxp)
      integer n ,indmax(nmaxp)
      character (len=200) filename
! Added by Gullo
      integer unit, status
      double precision hhh,skkk,ccc,smhy,smel,esu,eee,sigma,pi,eve,evf
      parameter (hhh= 6.62620d-27, skkk= 1.38062d-16, ccc= 2.99793d+10,
     &           smhy= 1.67333d-24, smel= 9.10956d-28, esu= 4.80325d-10,
     &           eee= 4.80325d-10, sigma=5.66961d-05, pi=3.141592654d0)
c Initialize status
      status=0
c
      ikbol   = 1.16d4      ! inverse of kbol (K * ev-1)
      mec2  = 5.11d5        ! m_e c^2 (eV)
      isp = 0.5641895835d0  ! 1/sqrt(pi)
      limit = 1.d-3         ! Limit for the redistribution function
c
c$$$100   format(2i8)
c$$$101   format(i8,ES16.8E3)
c$$$102   format(i8,ES16.8E3,i8,ES16.8E3)
c$$$103   format(2i8,2ES16.8E3)
c
c     Check1 is to ensure photon number is conserved in scatterings

      filename = 'angle.fits' !name of the fits file
      n = 1 ! column number of the fits file

      call create_fits(filename) !create the fits file

      do iz=1, itrans
            temp(iz) = theta(iz)*ikbol*mec2
      enddo
!crate and fill the first extension with temperature and energy
      call write_param_angle(itrans, nmaxp, temp, wp, filename,mgi,
     & smit,agt)
!append and other extension (3rd) to store the SRF
!IMPORTANT: this routine leaves the fits file opened      
      call add_HDU(itrans, nmaxp,filename, unit)
      do iz = 1, itrans
         x=1/theta(iz)
         do kk=1,nmaxp
            point(kk)=0
            pmax(kk)=0.d0
            check1(kk) = 0.d0
            indmax(kk)=0
            do np = 1, nmaxp
               prob(np,kk) = 0.d0
               indices(np,kk)=0
               do qq=1,mgi
                  srf(qq,kk,np)=0.d0
               enddo
            enddo
         enddo
c
c         temp = theta(iz)*ikbol*mec2          ! temperature in K
         do np = 1, nmaxp
            ecen = wp(np)
c$omp parallel num_threads(30)
c$omp& shared(iz,wp,prob,point,pmax,nmaxp,np,ecen,temp,mgi,smit,
c$omp& agt,mec2,indmax)
c$omp& private(jj)
c$omp do
            do jj=1,nmaxp
               do qq=1,mgi
                  srf(qq,jj,np)=profil(1,wp(np)/mec2,wp(jj)/mec2,
     &             smit(qq),x)
                  prob(jj,np)=prob(jj,np)+agt(qq)*srf(qq,jj,np)
               enddo
               if(prob(jj,np).gt.pmax(np))then
                  pmax(np)=prob(jj,np)
                  indmax(np)=jj
               endif
            enddo
c$omp end do
c$omp end parallel
         enddo
         do np=1,nmaxp
            do jj=indmax(np),1,-1
                  kk=point(jj)+1
                  if(prob(jj,np).ge.(pmax(np)*limit))then
                        indices(jj,kk)=np
                        check1(np)=check1(np)+df(jj)*prob(jj,np)
                        point(jj)=kk
                  else
                        exit
                  endif
            enddo

            do jj=indmax(np)+1,nmaxp
                  kk=point(jj)+1
                  if(prob(jj,np).ge.(pmax(np)*limit))then
                        indices(jj,kk)=np
                        check1(np)=check1(np)+df(jj)*prob(jj,np)
                        point(jj)=kk
                  else
                        exit
                  endif

            enddo
            srf(:,:,np)=srf(:,:,np)*skn(np,iz)*df(np)/wp(np)/check1(np)
         enddo

         
!     From here crated by Gullo Dec 2020
!     Once it calculates the srf it writes directly the fits file
!     It needs to differentiate the first call, where it writes the extension from all the other calls 
      do jj = 1, nmaxp
            if(point(jj).gt.0)then
            call add_row_HDU_angle(n,nmaxp,point(jj),indices(jj,1),mgi,
     &       skn(jj,iz),
     &        srf(:,jj,indices(jj,1):indices(jj,point(jj))),unit)
            else
            call add_row_HDU_angle(n,nmaxp,1,jj,mgi,
     &       skn(jj,iz),
     &        srf(:,jj,jj),unit)
            endif
            n = n + 1  
      enddo
               
      print *,iz
      enddo      
c
     
c The FITS file must always be closed before exiting the program. 
c Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      write(*,*) 'final status', status
      call ftclos(unit, status)
      write(*,*) 'final status', status
      call ftfiou(unit, status)
      write(*,*) 'final status', status
      if (status .gt. 0) then
         call printerror(status)
c$$$         write(*,*) 'end file'
      endif
      return
      end subroutine

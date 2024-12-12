!-------------------------------------------------------------------------------------------------      
      subroutine super_Compton_RF_fits_angle(itrans, temps, theta,
     &                                       nmaxp, wp, df, skn,
     &                                       mgi, nksmps, limit) 
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
      integer itrans, nmaxp, mgi, iz, np, jj, qq, kk, mu
      integer nksmps, ang1, ang2, inang, outang
      integer j_pmax, jlow, jhigh
      logical GOUP, GODOWN
      integer, parameter :: MUintnum = 1000
      real*8 MUintp(MUintnum), MUintw(MUintnum)
      real*8 A, B, ecen
      real*8 theta(itrans), wp(nmaxp), df(nmaxp), skn(nmaxp,itrans),x
      real*8 smit(mgi), agt(mgi), halfsmit(mgi+1)
      real*8 thetak(nksmps+1), muk(nksmps+1), dthetak
      real*8 phik, phik1, phiint
      real*8 check, factor
      real*8 profil, temps(itrans)
      real*8 limit, problim, total, maxtot
      real*8 srfval, kord
      real*8, target, allocatable :: srfe(:,:,:,:), srfo(:,:,:,:)                    
      real*8, pointer :: flatsrfe(:,:), flatsrfo(:,:)
      integer, target, allocatable :: fSInde(:), fLene(:) 
      integer, target, allocatable :: fSIndo(:), fLeno(:) 
      real*8, allocatable :: we(:, :, :), wo(:, :, :)
      real*8 wval
      integer, pointer :: flatfse(:,:), flatfle(:,:)
      integer, pointer :: flatfso(:,:), flatflo(:,:)
      integer n

      INTERFACE
         PURE SUBROUTINE kloop(e1, e2, x, nksmps, muk, we, wo, mgi,
     1                         total, se, so)
            integer, intent(IN) :: nksmps, mgi
            real*8, intent(IN) :: e1, e2, x
            real*8, intent(IN) :: muk(nksmps)
            real*8, intent(IN) :: we(mgi, mgi, nksmps)
            real*8, intent(IN) :: wo(mgi, mgi, nksmps)
            real*8, intent(OUT) :: total
            real*8, intent(OUT) :: se(mgi, mgi), so(mgi, mgi)
         END SUBROUTINE kloop
         PURE DOUBLE PRECISION FUNCTION getphi(v)
            real*8, INTENT(IN) :: v
         END FUNCTION getphi
      END INTERFACE

      call set_filename('angle.fits') !name of the fits file
      n = 1 ! column number of the fits file
 
c     Get the angles - we only need positive ones
      call gaulegf(0.d0, 1.d0, smit, agt, mgi)

c     Set up a new fits file
      call setup_new_file(nmaxp, itrans, mgi,
     &                    wp, df, temps, smit, agt,
     &                    skn, limit, .FALSE., nksmps) !create the fits file


      

      allocate(srfe(mgi,nmaxp,mgi,nmaxp), srfo(mgi,nmaxp,mgi,nmaxp))
      allocate(fSInde(mgi*nmaxp*itrans), fLene(mgi*nmaxp*itrans))
      allocate(fSIndo(mgi*nmaxp*itrans), fLeno(mgi*nmaxp*itrans))
      allocate(we(mgi, mgi, nksmps), wo(mgi, mgi, nksmps) )
c     flatsrf is a pointer to allow accessing to the response function as a matrix
c     flatsrf( init, final )
      flatsrfe(1:mgi*nmaxp, 1:mgi*nmaxp) => srfe
      flatsrfo(1:mgi*nmaxp, 1:mgi*nmaxp) => srfo

c     flatfs and flatfl are pointers to fSInd and fLen to make it easier to deal with the running temperature
      flatfse(1:mgi*nmaxp, 1:itrans) => fSInde
      flatfle(1:mgi*nmaxp, 1:itrans) => fLene

      flatfso(1:mgi*nmaxp, 1:itrans) => fSIndo
      flatflo(1:mgi*nmaxp, 1:itrans) => fLeno

!     Get half steps in angular grid
      halfsmit(1) = 0.d0
      do qq = 1, mgi
         halfsmit(qq+1) = halfsmit(qq) + agt(qq)
      enddo

!     halfsmit should end at 1. This is a small correction for numerical precision
      halfsmit = halfsmit / halfsmit(mgi+1)     


!     0 < thetak < pi      - note, this theta is angle, not temperature
!     Not allowed to be 0, as cos(theta)=1 causes problems
      do qq = 1, nksmps+1
          thetak(qq) = pi * (DBLE(qq) / DBLE(nksmps+2) )
          muk(qq) = COS(thetak(qq))
      enddo
      dthetak = pi / DBLE(nksmps)
      
      we = 0.d0
      wo = 0.d0
      ! These loops produce the w(j', j, k) from eqs A6/A7.
      ! Split into symmetric and anti-symmetric parts
      do ang1=1, mgi
         call gaulegf( halfsmit(ang1), halfsmit(ang1+1),  
     1                 MUintp, MUintw, MUintnum )
         do ang2=ang1,mgi
!$omp parallel
!$omp& shared(ang1, ang2, muk, smit, MUintp, MUintw,
!$omp&        we, wo, nksmps)  
!$omp& private(mu, kk, A, B, phik, phik1, phiint)
!$omp do         
            do mu=1,MUintnum
               A = 1/sqrt((1-MUintp(mu)**2)*(1-smit(ang2)**2))
               B = MUintp(mu)*smit(ang2)
               do kk=1, nksmps
                  ! When mu' and muj have same sign
                  phik  = getphi( A*(muk(kk) - B) )
                  phik1 = getphi( A*(muk(kk+1) - B) )
                  phiint = (phik1 - phik) * MUintw(mu)
                  we(ang2, ang1, kk) = we(ang2, ang1, kk)
     1                                 + phiint
                  wo(ang2, ang1, kk) = wo(ang2, ang1, kk)
     1                                 + phiint
                  ! When mu' and muj have same sign
                  phik  = getphi( A*(muk(kk) + B) )
                  phik1 = getphi( A*(muk(kk+1) + B) )
                  phiint = (phik1 - phik) * MUintw(mu)
                  we(ang2, ang1, kk) = we(ang2, ang1, kk)
     1                                 + phiint
                  wo(ang2, ang1, kk) = wo(ang2, ang1, kk)
     1                                 - phiint
               enddo
            enddo
!$omp end do
!$omp end parallel
         enddo
      enddo
      print *,"Weights computed"

   !   ! Fill in the other triangle of the matrix
   !   ! Could possibly be done with pointers, but it's a fairly small matrix
   !   do ang1=1, mgi
   !      do ang2=1,ang1-1
   !         do kk=1, nksmps
   !            we(kk, ang2, ang1) = we(kk, ang2, ang1)
   !            wo(kk, ang2, ang1) = wo(kk, ang1, ang2)
   !         enddo
   !      enddo
   !   enddo


      ! Now, loop through and compute Fk(x',x)
      do iz = 1, itrans
         x=1/theta(iz)
         srfe = 0.d0
         srfo = 0.d0
!$omp parallel
!$omp& shared(wp, nmaxp, x, muk, limit)
!$omp& private(np, jlow, jhigh, ecen, problim,
!$omp&         GOUP, GODOWN, total, maxtot)
!$omp do
         !!!
         ! Prepare Fk(x',x) as in Eq. A5.
         ! As stated in 2nd paragraph of page '257',
         ! "The Fk's are computed simply by evaluating R at the midpoint of each interval" 
         !!!
         do np=1, nmaxp ! initial energy
            ecen = wp(np)/mec2
            call kloop(ecen, ecen, x, nksmps, muk, we, wo, mgi,
     1                 maxtot, srfe(:,np,:,np), srfo(:,np,:,np))
            problim = limit * maxtot
            jlow = np-1
            jhigh = np+1
            GOUP = .TRUE.
            GODOWN = .TRUE.
            ! final energy
            do while (GOUP .or. GODOWN)
               if(jlow.LT.1)       GODOWN = .FALSE.
               if(jhigh.GT.nmaxp)  GOUP = .FALSE.
               if(GODOWN)then
                  call kloop(ecen, wp(jlow)/mec2, x, 
     1                       nksmps, muk, we, wo, mgi, 
     2                       total, 
     3                       srfe(:,np,:,jlow), srfo(:,np,:,jlow))
                  if(total.LT.problim)then
                     GODOWN = .FALSE.
                  elseif(total .GT. maxtot)then
                     maxtot = total
                     problim = total*limit
                  endif
                  jlow = jlow-1
               endif
               if(GOUP)then
                  call kloop(ecen, wp(jhigh)/mec2, x, 
     1                       nksmps, muk, we, wo, mgi, 
     2                       total, 
     3                       srfe(:,np,:,jhigh), srfo(:,np,:,jhigh))
                  if(total.LT.problim)then
                     GOUP = .FALSE.
                  elseif(total .GT. maxtot)then
                     maxtot = total
                     problim = total*limit
                  endif
                  jhigh = jhigh+1
               endif
            enddo ! GOUP/GODOWN while
         enddo ! np
!$omp end do
!$omp end parallel

        
c        srf( init_ang, init_en, final_ang, final_en  )
c        flatsrf( init_pair, final_pair)
         call srf_nonlimit( flatsrfe(:,:), nmaxp*mgi, limit, 
     &                      flatfse(:,iz), 
     &                      flatfle(:,iz)) 
c
         call srf_nonlimit( flatsrfo(:,:), nmaxp*mgi, limit, 
     &                      flatfso(:,iz), 
     &                      flatflo(:,iz)) 
         

c        Go through each 'final' energy/angle row, and set the areas to 0 when below limit
         do jj=1, nmaxp*mgi
            flatsrfe(1:flatfse(jj,iz)-1,jj) = 0.d0
            flatsrfe(flatfse(jj,iz)+flatfle(jj,iz):nmaxp*mgi,jj) = 0.d0

            flatsrfo(1:flatfso(jj,iz)-1,jj) = 0.d0
            flatsrfo(flatfso(jj,iz)+flatflo(jj,iz):nmaxp*mgi,jj) = 0.d0
         enddo

c        Check to ensure photon number is conserved in scatterings.  
c        0.5 * Int_{-1}^{+1} Int_0^\inf R[init_mu, init_en, fin_mu, fin_en] d{fin_en} d{fin_my} = 1
c        Which becomes:
c        Int_{0}^{+1} Int_0^\inf Re[init_mu, init_en, fin_mu, fin_en] d{fin_en} d{fin_my} = 1

c        srf( init_ang, init_en, final_ang, final_en  )
         do inang=1,mgi ! inital angle
            do np=1,nmaxp ! initial energy
               check = 0.d0
               do outang=1,mgi ! final angle
                  do jj=1,nmaxp ! final energy
                     check = check  
     1                  + srfe(inang,np,outang,jj)*agt(outang)*df(jj)
                  enddo
               enddo
               ! write(*,*) check
               factor = skn(np,iz)*df(np)  ! *smit(inang)
     1                  *agt(inang)/wp(np)/(check)
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
      deallocate(srfe, srfo, fSInde, fLene, fSIndo, fLeno, we, wo)
      return
      end subroutine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      PURE DOUBLE PRECISION FUNCTION getphi(v)
         use constants, only: pi
         implicit none
         real*8, INTENT(IN) :: v
         if (v.LT.-1)then
           getphi = pi
         elseif (v.GT.1)then
           getphi = 0.d0
         else
           getphi = ACOS(v)
         endif
      END FUNCTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      PURE SUBROUTINE kloop(e1, e2, x, nksmps, muk, we, wo, mgi,
     1                      total, se, so)
         implicit none
         INTERFACE
            PURE double precision function profil(nr,eps,eps1,costh,x)
               integer, intent(IN) :: nr
               double precision, intent(IN) :: eps, eps1, costh, x
            end function profil
         END INTERFACE
         integer, intent(IN) :: nksmps, mgi
         real*8, intent(IN) :: e1, e2, x
         real*8, intent(IN) :: muk(nksmps)
         real*8, intent(IN) :: we(mgi, mgi, nksmps)
         real*8, intent(IN) :: wo(mgi, mgi, nksmps)
         real*8, intent(OUT) :: total
         real*8, intent(OUT) :: se(mgi, mgi), so(mgi, mgi)
         real*8 Fk
         integer kk, inang, outang

         total = 0
         se = 0
         so = 0
         do kk=1, nksmps
            Fk =  profil(1, e1, e2, muk(kk), x)
            total = total + Fk
            do outang=1, mgi
               ! Backfill angles from symmetry
               do inang=1, outang-1
                  se(inang,outang) = se(inang,outang)
               enddo

               do inang=outang, mgi
                  se(inang,outang) = se(inang,outang) + 
     1                                 Fk * we(inang, outang, kk)
                  so(inang,outang) = so(inang,outang) + 
     1                                 Fk * wo(inang, outang, kk)
                  
               enddo
            enddo
         enddo
         ! Integration is over evenly space theta
         ! mu variable is just a shortcut

         ! As total is just used for a prob limit, any scaling factor is irrelevent
         return
      END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

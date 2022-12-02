      program drive_SRF 
c     Driving program for the writting of the Super Redistribution Function SRF
c     due to Compton scattering (see more details in super_Compton_RF.f).
c     The output of this code is used by the XILLVER code. 
c
c     It writes to an ascii file for the moment.
c
c     To-Do:
c    - The minimum number of angles (mgi) needed for accurate angular integration
c      depends on temp and initial photon energy. 3000 seems to be the worse case  
c      scenario, but we could make mgi to change according to the conditions for speed up.
c    - Write the SRF's into a FITS file.
c
c      Done in v0.2.1:
c    - Fixed the NaNs printed in the output by setting df(1)=df(2) in enegrd.f
c
c      Done in v0.2a:
c    - Replaced the kleinni.f SR for scattxs.f, which computes the exact Compton
c      scatteting cross section including Klein-Nishina corrections at high energies
c      and relativistic corrections at high temperatures
c    - The total cross section skn is now written in the output file for each
c      temperature and finel energy.
c
c     Version: 0.2.1 - Wed Aug 21 19:13:39 PDT 2019
c     Version: 0.2a -Wed Apr 24 19:01:50 PDT 2019
c     Version: 0.1a -Tue Apr 23 16:02:54 PDT 2019
c
c     Authors: Javier Garcia (javier@caltech.edu)
c              Ekaterina Sokolova-Lapa (ekaterina.sokolova-lapa@fau.de)
c
c........    
      use omp_lib
      implicit none
      integer nmaxp, itrans, mgi
      parameter (nmaxp=500, itrans=70, mgi=8)
      logical avangle
      parameter (avangle=.false.)
      integer ii
      double precision pemin, pemax, pemax2
      double precision theta(itrans), wp(nmaxp), df(nmaxp)
      double precision skn(nmaxp,itrans)
      double precision smit(mgi), agt(mgi)
      double precision tini, tfin, tcpu, temp
     

C     This line is only ran if the compiler can handle parallisation
!$    write(*,*)"Parallised over ",OMP_get_max_threads()

c     Get current time
      call cpu_time(tini)
c
c

c
c     Array of temperatures
      temp = 1.d4                ! Gas temp in K
      do ii=1,itrans
         call theta_from_temp(temp, theta(ii))
         temp = temp*(1.d10/1.d4)**(1.d0/dfloat(itrans-1))
      enddo
c
c     Photon energy grid
      pemin = 0.1d0
      pemax = 9.9d5
      pemax2 = 1.d6
      call enegrd(nmaxp, pemin, pemax, pemax2, wp, df)
c
c     Calculate the Compton Cross Section
      call scattxs(nmaxp, wp, itrans, theta, skn)
c
c     Get the Gaussian quadratures for angular integration
      call gaulegf(-1.d0, 1.d0, smit, agt, mgi)
c     Produce file with all SRF's
      if (avangle) then
            call super_Compton_RF_fits(itrans, theta, nmaxp, wp, df,
     &                                 skn, mgi, smit, agt)
      else
            call super_Compton_RF_fits_angle(itrans, theta, nmaxp, wp,
     &                                       df, skn, mgi, smit, agt)
      endif
c     Get current time
      call cpu_time(tfin)
      tcpu=tfin-tini
c
      print *, ' '
      print *, 'CPU time (s) =',real(tcpu)
c
      end program drive_SRF

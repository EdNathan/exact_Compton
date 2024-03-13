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
      integer nmaxp, itrans, mgi, knsamps
      ! parameter (nmaxp=500, itrans=70, mgi=8)
      logical avangle
      ! parameter (avangle=.false.)
      character*(*) input_file
      character (len=255) file_line
      parameter (input_file = 'drive_params.dat') 
      logical exist
      integer ii, file_unit, file_stat, convert_stat
      double precision pemin, pemax, pemax2
      double precision, allocatable :: theta(:), temps(:), wp(:), df(:)
      double precision, allocatable ::  skn(:,:)
      double precision, allocatable ::  smit(:), agt(:)
      double precision tini, tfin, tcpu, temp
      double precision limit

c     Default values
      limit = 1.d-3
      nmaxp = 500
      itrans = 70
      mgi = 3000
      knsamps = 0
      avangle = .true.

      exist = .false.
      file_unit = 101
      file_stat = 1

c     Handle parameter file
      inquire( file=input_file, exist=exist )
      if (exist) then
        open(newunit=file_unit, file=input_file, 
     &       action='READ', iostat=file_stat )
c       If file exists, loop through lines            
        do while (file_stat .eq. 0)
          read(file_unit, '(A255)', iostat=file_stat) file_line
c         Exit if EOF
          if (file_stat .ne. 0) exit
c         Read possible parameters
          file_line = trim(file_line)
          convert_stat = 1
          if (file_line(:6) .eq. 'nmaxp') then
            read( file_line(6:), *, iostat=convert_stat) nmaxp
          elseif (file_line(:6) .eq. 'itrans') then
            read( file_line(7:), *, iostat=convert_stat) itrans
          elseif (file_line(:3) .eq. 'mgi') then
            read( file_line(4:), *, iostat=convert_stat) mgi
          elseif (file_line(:7) .eq. 'knsamps') then
            read( file_line(8:), *, iostat=convert_stat) knsamps
          elseif (file_line(:5) .eq. 'limit') then
            read( file_line(6:), *, iostat=convert_stat) limit
          endif
          if (convert_stat .ne. 0) then
            write(*,*)"Failed to understand ",trim(file_line)
          endif
        end do
        close(file_unit)
      endif

      avangle = (knsamps .le. 0)

      write(*,*)'Using nmaxp:   ', nmaxp
      write(*,*)'Using itrans:  ', itrans
      write(*,*)'Using mgi:     ', mgi
      write(*,*)'Using knsamps: ', knsamps
      write(*,*)'Using avangle: ', avangle
      write(*,*)'Using limit:   ', limit

c     Allocate arrays
      allocate( theta(itrans), temps(nmaxp), wp(nmaxp), df(nmaxp) )
      allocate( skn(nmaxp,itrans) )
      

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
         temps(ii) = temp
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

c     Produce file with all SRF's
      if (avangle) then
            allocate ( smit(mgi), agt(mgi) )
c           Get the Gaussian quadratures for angular integration
            call gaulegf(-1.d0, 1.d0, smit, agt, mgi)
            call super_Compton_RF_fits(itrans, temps, theta, 
     1                                 nmaxp, wp, df, skn, mgi,
     2                                 smit, agt, limit)
      else
c           We only need postive angles
            allocate ( smit(mgi), agt(mgi) )
            call gaulegf(0.d0, 1.d0, smit, agt, mgi)
            call super_Compton_RF_fits_angle(itrans, temps, theta, 
     1                                       nmaxp, wp,df, skn, mgi,
     2                                       smit, 
     3                                       agt,
     4                                       knsamps, limit)
      endif
c     Get current time
      call cpu_time(tfin)
      tcpu=tfin-tini
c
      print *, ' '
      print *, 'CPU time (s) =',real(tcpu)
c
      end program drive_SRF

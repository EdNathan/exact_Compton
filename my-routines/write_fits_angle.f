c-----------------------------------------------------------------------
c     
c      Subroutines to write the fits file
c
c-----------------------------------------------------------------------
c
c      Angle specific routines
c
c-----------------------------------------------------------------------


c-----------------------------------------------------------------------
      subroutine write_param_angle(temp_dim, en_dim, temp, en,
     &                             filename, mgi, smit, agt)
c     Create fits file with one extension and  columns 
      implicit none
      integer          :: temp_dim, en_dim, mgi
      double precision :: temp(temp_dim), en(en_dim),smit(mgi),agt(mgi)
      character* (*) filename
c     Internal variables
      integer status,unit,blocksize
      integer rownum,readwrite,colnum,tfields
      character (len=16) extname, name
      integer,parameter :: coldim = 4
      character (LEN=16) ttype(coldim),tform(coldim),tunit(coldim)
      integer nrows, varidat
c     inverse of kbol (K * ev-1)
!     double precision, parameter :: ikbol   = 1.16d4
c     m_e c^2 (eV)
!     double precision, parameter :: mec2  = 5.11d5
      
c     Initialize status
      status=0
      blocksize = 1

c     Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

c     Open the FITS file, with write access.
      readwrite=1
      call ftopen(unit,filename,readwrite,blocksize,status)   

c     Append/create a new empty HDU onto the end of the file and move to it.
      call ftcrhd(unit,status)

      ttype(1) = 'TEMP'
      ttype(2) = 'ENERGIES'
      ttype(3) = 'POINTS'
      ttype(4) = 'WEIGHTS'

      write(name, '(I4,A1)') temp_dim, 'D'
      tform(1) = trim(name)
      
      write(name, '(I4,A1)') en_dim, 'D'
      tform(2) = trim(name)

      write(name, '(I4,A1)') mgi, 'D'
      tform(3) = trim(name)

      write(name, '(I4,A1)') mgi, 'D'
      tform(4) = trim(name)

      write(*,*) tform(1) , tform(2) , tform(3) , tform(4)
!     tform(1) = '70D'
!     tform(2) = '500D'

      tunit(1) = 'kT/mec2'
      tunit(2) = 'eV'
      tunit(3) = ''
      tunit(4) = ''

c     Define parameters for the binary table (see the above data statements)
      tfields  = coldim
      nrows    = 1
      extname  = 'PARAMETERS'
      varidat  = 0

c     FTPHBN writes all the required header keywords which define the
c     structure of the binary table. NROWS and TFIELDS gives the number of
c     rows and columns in the table, and the TTYPE, TFORM, and TUNIT arrays
c     give the column name, format, and units, respectively of each column.
      call ftphbn(unit, nrows, tfields, ttype, tform, tunit, 
     &            extname, varidat, status)

c     Filling the columns 
      rownum = 1
      colnum = 1
      call ftpcld(unit, colnum, rownum, 1, temp_dim, temp, status)
      colnum = 2
      call ftpcld(unit, colnum, rownum, 1, en_dim, en, status)      
      colnum = 3
      call ftpcld(unit, colnum, rownum, 1, mgi, smit, status)     
      colnum = 4
      call ftpcld(unit, colnum, rownum, 1, mgi, agt, status)
       

c     Write keywords to this extension
c     ftpkyj to write integer
c     tpkye to write real (4)
c     ftpkyd to write real (8)

c     The FITS file must always be closed before exiting the program. 
c     ny unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(unit, status)
      call ftfiou(unit, status)

c     Check for any error, and if so print out error messages.
      if (status .gt. 0) call printerror(status)
      return
      end
c-----------------------------------------------------------------------

      

c-----------------------------------------------------------------------
      subroutine add_row_HDU_angle(n, nmaxp, out_en_dim, out_en_ind, 
     &                             mgi, kn_cross, srf, unit)
      implicit none
      integer          :: unit, n, nmaxp, out_en_dim, out_en_ind,mgi
      double precision :: srf(mgi,nmaxp)
      double precision :: kn_cross(1)
      integer status, hdutype
!     integer readwrite, blocksize
      integer colnum,rownum
      
c     Initialize status
      status=0
!     blocksize = 1

c     Get an unused Logical Unit Number to use to open the FITS file.
!     call ftgiou(unit,status)

c     Open the FITS file, with write access.
!     readwrite=1
!     call ftopen(unit,filename,readwrite,blocksize,status)

c     Move to the last (3nd) HDU in the file (the paameter values table).
      call ftmahd(unit,3,hdutype,status)

c     Filling the columns 
      rownum = n
      colnum = 1
      call ftpcld(unit, colnum, rownum, 1, 1, kn_cross(1), status)    

      colnum = 2
      call ftpcli(unit, colnum, rownum, 1, 1, out_en_ind, status) 

      colnum = 3
      call ftpcli(unit, colnum, rownum, 1, 1, out_en_dim, status)   

      colnum = 4
c     Fill the column with arrays of different size
      call ftpcld(unit, colnum, rownum, 1, out_en_dim*mgi, srf, status)         

c     The FITS file must always be closed before exiting the program. 
c     Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
!     call ftclos(unit, status)
!     call ftfiou(unit, status)

c     Check for any error, and if so print out error messages.
            if (status .gt. 0) then
               call printerror(status)
            endif
      return
      end
c-----------------------------------------------------------------------
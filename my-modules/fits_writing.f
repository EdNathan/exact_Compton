      module fits_writing
         implicit none
c        Everything in here is, by default, private
         private
c        Explicitly export public entities
         public set_filename, setup_new_file
         public write_SRF_pointers, write_SRFs   
         public close_and_save_fits

         integer file_en_dim, file_temp_dim, file_mgi
         integer :: unit=-1, status=0, curr_tab=-1
         character*(100) :: filename='SRF.fits'



      contains
      
c-----------------------------------------------------------------------       
         subroutine set_filename(new_name)
            implicit none
            character*(*)new_name

            if(unit.LE.0)then
               filename = TRIM(new_name)
            endif
         end subroutine set_filename
c-----------------------------------------------------------------------



c-----------------------------------------------------------------------
         subroutine print_any_errors()
            implicit none
            character errtext*30,errmessage*8
c           This subroutine prints out the descriptive text corresponding to the
c           error status value and prints out the contents of the internal
c           error message stack generated by FITSIO whenever an error occurs.
c           character errtext*30,errmessage*80

c           Check if status is OK (no error); if so, simply return
            if (status .le. 0)return

c           The FTGERR subroutine returns a descriptive 30-character text string that
c           corresponds to the integer error status number.  A complete list of all
c           the error numbers can be found in the back of the FITSIO User's Guide.
            call ftgerr(status,errtext)
            print *,'FITSIO Error Status =',status,': ',errtext

c           FITSIO usually generates an internal stack of error messages whenever
c           an error occurs.  These messages provide much more information on the
c           cause of the problem than can be provided by the single integer error
c           status value.  The FTGMSG subroutine retrieves the oldest message from
c           the stack and shifts any remaining messages on the stack down one
c           position.  FTGMSG is called repeatedly until a blank message is
c           returned, which indicates that the stack is empty.  Each error message
c           may be up to 80 characters in length.  Another subroutine, called
c           FTCMSG, is available to simply clear the whole error message stack in
c           cases where one is not interested in the contents.
            call ftgmsg(errmessage)
            do while (errmessage .ne. ' ')
                  print *,errmessage
                  call ftgmsg(errmessage)
            end do
         end subroutine print_any_errors
c-----------------------------------------------------------------------



c-----------------------------------------------------------------------
         subroutine setup_new_file(en_dim, temp_dim, mgi,
     &                             en, den, temp, smit, agt,
     &                             skn, limit, avangle, nksamps)
            implicit none
c           Input
            integer en_dim, temp_dim, mgi, nksamps
            double precision en(en_dim), den(en_dim), temp(temp_dim)
            double precision smit(mgi), agt(mgi)
            double precision skn(en_dim, temp_dim)
            double precision limit
            logical avangle
c           Internal variables
            logical simple, extnd
            integer bitpix,naxis,naxes(1)
            integer blocksize, readwrite
            integer coldim 
            character (LEN=16),allocatable:: ttype(:),tform(:),tunit(:)
            character (len=16) extname
            integer colnum, rownum
            integer tfields, nrows, varidat
            integer i

            blocksize = 1
            if (unit .GT. 0) then
               return
            endif
            
            file_en_dim = en_dim 
            file_temp_dim = temp_dim
            if( avangle )then
               file_mgi = 1
            else
               file_mgi = mgi
            endif

            coldim = 5 ! MAX(5, file_mgi)

c           If file already exists, delete it!
            call deletefile()


c           Get unused unit number
            call ftgiou(unit,status)

c           Open new fits file for reading 
            readwrite=1
            call ftinit(unit,TRIM(filename),blocksize,status)
            
            
            ! call ftopen(unit,filename,readwrite,blocksize,status)
            
c           Initialize parameters about the FITS primary image - need image for structure
            simple=.true.
            bitpix=16
            naxis=1
            naxes(1)=0
            extnd=.true.
c           Write the required header keywords to the file
            call ftphpr(unit,simple,bitpix,naxis,
     &                   naxes,0,1,extnd,status)
            
c           Allocate arrays
            ALLOCATE( ttype(coldim), tform(coldim), tunit(coldim) )

c          Now, we add the tables we will want!

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c           First table / 2nd extension - Parameters
c           Add new table, move to it      
            call ftcrhd(unit,status)
            if ( avangle ) then
              tfields  = 3
            else
              tfields  = 5
            endif
            nrows    = 1
            extname  = 'Parameters'
            varidat  = 0

c           Write header
c           - Column names
            ttype(1) = 'TEMP'
            ttype(2) = 'ENERGIES'
            ttype(3) = 'del ENERGY'
            ttype(4) = 'MU'
            ttype(5) = 'WEIGHTS'
c           - Column data formats 
            write(tform(1), '(I15,A1)') temp_dim, 'D'
            tform(1) = trim(tform(1))
            write(tform(2), '(I15,A1)') en_dim, 'D'
            tform(2) = trim(tform(2))
            write(tform(3), '(I15,A1)') en_dim, 'D'
            tform(3) = trim(tform(3))
            write(tform(4), '(I15,A1)') mgi, 'D'
            tform(4) = trim(tform(4))
            write(tform(5), '(I15,A1)') mgi, 'D'
            tform(5) = trim(tform(5))
c           - Column units
            tunit(1) = 'K'! 'kT/mec2'
            tunit(2) = 'eV'
            tunit(3) = 'eV'
            tunit(4) = ''
            tunit(5) = ''

            call ftphbn(unit, nrows, tfields, 
     &            ttype(1:tfields), tform(1:tfields), tunit(1:tfields),
     &            extname, varidat, status )
c           Write useful headers
            call ftpkyj(unit, 'EN_DIM', en_dim, 
     &             "Number of energy points in grid", 
     &              status)
            call ftpkyj(unit, 'TEMP_DIM', temp_dim, 
     &             "Number of temperatures in grid",
     &              status)
            call ftpkyj(unit, 'ANGLES', mgi, 
     &             "Number of angles used", 
     &              status)
            if(.not.AVANGLE)then
               call ftpkyj(unit, 'AZIMUTHAL', nksamps, 
     &               "Number of points used in azimuthal integration", 
     &               status)
            endif
            call ftpkyl(unit, 'AVANGLE', avangle, 
     &             "Whether the SRF has been integrated over angles", 
     &              status)
            call ftpkyd(unit, 'LIMIT', limit, 10, 
     &             "Fractional limit of maximum
     &              below which SRF is set to zero", 
     &              status)    
c           Insert data 
            rownum = 1
            colnum = 1
            call ftpcld(unit, colnum, rownum, 1, temp_dim, temp, 
     &                  status)
            colnum = 2
            call ftpcld(unit, colnum, rownum, 1, en_dim, en, status)
            colnum=3
            call ftpcld(unit, colnum, rownum, 1, en_dim, den, status)
            if ( .not. avangle ) then
               colnum = 4
               call ftpcld(unit, colnum, rownum, 1, mgi, smit, status)
               colnum = 5
               call ftpcld(unit, colnum, rownum, 1, mgi, agt, status)
            endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c           Second table / 3rd extension - Cross section

c           Add new table, move to it      
            call ftcrhd(unit,status)

            tfields  = 1
            nrows    = 1
            extname  = 'SKN'
            varidat  = 0

c           Write header 
c           - Column names
            ttype(1) = 'KN_cross'
c           - Column data formats 
            write(tform(1), '(I15,A1)') en_dim * temp_dim, 'D'
            tform(1) = trim(tform(1))
c           - Column units
            tunit(1) = 'cm^2'

            call ftphbn(unit, nrows, tfields, 
     &            ttype(1:tfields), tform(1:tfields), tunit(1:tfields),
     &            extname, varidat, status )

c           Insert data 
            rownum = 1
            colnum = 1
            call ftpcld(unit, colnum, rownum, 1, 
     &                  (en_dim * temp_dim), skn, status)
            
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c           Third table  / 4th extension - SRF pointers   

c           Add new table, move to it      
            call ftcrhd(unit,status)
            
            tfields  = 2
            if(avangle)then
               nrows    = 1
            else
               nrows    = 2
            endif
            extname  = 'iSRF_pointers'
            varidat  = 0
         
c           Write header 
c           - Column names
            ttype(1) = 'IND'
            ttype(2) = 'LEN'
c           - Column data formats      
            write(tform(1),'(I15,A1)') (en_dim*temp_dim*file_mgi),'I'  
            tform(1:2) = trim(tform(1))
c           - Column units
            tunit(1:2) = ''   

            call ftphbn(unit, nrows, tfields, 
     &            ttype(1:tfields), tform(1:tfields), tunit(1:tfields),
     &            extname, varidat, status )
            
c           Not inserting data yet
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c           Fourth table / 5th extension - SRF

c           Add new table, move to it      
            call ftcrhd(unit,status)
            if (avangle)then
               tfields  = 1
            else
               tfields = 2
            endif
            nrows    = file_mgi*en_dim * temp_dim
            extname  = 'iSRF'
            varidat  = 0
c           Write header 

c           - Column name           
            if (avangle)then
               ttype(1) = 'SRF'
            else
               ttype(1) = 'SRFe'
               ttype(2) = 'SRFo' 
            endif
c           - Column data type
            tform(1:tfields) = '1PD'
c           - Column unit
            tunit(1:tfields) = ''                        
               
            call ftphbn(unit, nrows, tfields, 
     &            ttype(1:tfields), tform(1:tfields), tunit(1:tfields),
     &            extname, varidat, status )
            
c           Not inserting data yet
cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
            
c           Deallocate arrays
            DEALLOCATE( ttype, tform, tunit )

            call update_curr_tab()
            call print_any_errors()
         end subroutine setup_new_file
c-----------------------------------------------------------------------



c-----------------------------------------------------------------------
         subroutine write_SRF_pointers(rownum, SRF_ind, SRF_len)
            implicit none
            integer hdutype, i
            integer rownum, colnum
            integer writeLen
            integer SRF_ind(file_mgi*file_en_dim*file_temp_dim)
            integer SRF_len(file_mgi*file_en_dim*file_temp_dim)

c           If file isn't open, then don't proceed
            if(unit .LE. 0)then
               return
            endif
c           Move to the fourth tab, if that's not current
            call ensure_tab(4)
            
            writeLen = file_mgi * file_en_dim * file_temp_dim

            colnum = 1
c              Write with option "j" instead of "i" as we are using
c              regular integers, not short.  Perhaps we can swap to
c              short integers save memory, but it isn't significant
c              with the current values expected.
            call ftpclj(unit, colnum, rownum, 1, writeLen, 
     &                  SRF_ind, status)
            
            colnum = 2
            call ftpclj(unit, colnum, rownum, 1, writeLen,
     &                  SRF_len, status)
      

            call print_any_errors()
         end subroutine write_SRF_pointers
c-----------------------------------------------------------------------



c-----------------------------------------------------------------------
         subroutine write_SRFs(rownum, colnum, SRF, ind, num)
            implicit none
            integer hdutype, endind, i
            integer rownum, colnum
            integer ind, num
            double precision SRF(file_mgi*file_en_dim)
       
c           If file isn't open, then don't proceed
            if(unit .LE. 0)then
               return
            endif
c           Move to the fifth tab, if that's not current
            call ensure_tab(5)

c           Write the iSRF
            endind = ind + num
            call ftpcld(unit, colnum, rownum, 1, num, 
     &                  SRF(ind:endind), status)
            call print_any_errors()
         end subroutine write_SRFs
c-----------------------------------------------------------------------



c-----------------------------------------------------------------------
         subroutine deletefile()
c           A simple little routine to delete a FITS file
            implicit none
            integer blocksize

c           Simply return if status is greater than zero
            if (status .gt. 0)return

c         If file is currently open, don't do anything
            if(unit.GT.0)then
               write(*,*)'Cannot delete currently open file'
               return
            endif

c           Get an unused Logical Unit Number to use to open the FITS file
            call ftgiou(unit,status)

c           Try to open the file, to see if it exists
            call ftopen(unit,TRIM(filename),1,blocksize,status)
            if (status .eq. 0)then
c              File was opened;  so now delete it
               call ftdelt(unit,status)
            else if (status .eq. 103)then
c              File doesn't exist, so just reset status to zero and clear errors
               status=0
               call ftcmsg
            else
c              There was some other error opening the file; delete the file anyway
               status=0
               call ftcmsg
               call ftdelt(unit,status)
            end if
      
c           Free the unit number for later reuse
            call ftfiou(unit, status)

            call print_any_errors()

c           Reset unit and status

            unit = -1
            status = 0
            curr_tab = -1
         end subroutine deletefile
c-----------------------------------------------------------------------



c-----------------------------------------------------------------------
      subroutine update_curr_tab()
         implicit none
         if(unit.LE.0)then
            curr_tab = -1
         else
            call FTGHDN(unit, curr_tab)
            curr_tab = curr_tab
         endif
      end subroutine update_curr_tab
c-----------------------------------------------------------------------



c-----------------------------------------------------------------------
      subroutine ensure_tab(tabnum)
         implicit none
         integer tabnum, hdutype
         call update_curr_tab()
         if(curr_tab.NE.tabnum)then
            call ftmahd(unit, tabnum, hdutype, status)
         endif
      end subroutine ensure_tab
c-----------------------------------------------------------------------
      


c-----------------------------------------------------------------------
         subroutine close_and_save_fits()
            implicit none
            if(unit.LE.0)then
               return
            endif

c           Close file
            call ftclos(unit, status)
c           Free up unit number
            call ftfiou(unit, status)

            call print_any_errors()

            unit = -1
            curr_tab = -1
         end subroutine close_and_save_fits
c-----------------------------------------------------------------------
    


      end module fits_writing

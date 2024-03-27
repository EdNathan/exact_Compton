      subroutine srf_nonlimit(prob, nmaxp, limit, 
     &                        fSInd, fLen) ! wp , check
c     This finds the region of the SRF which is with the fractional limit of the maxmimum
c     probability.
c     ! It then integrates the SRF for each initial energy, to produce the check.
c
c     Input arguments: 
c         prob:   The SRF probabilities, given in coordinated (initial energy, final energy)
c         ! df:     Array (nmaxp) of delta energies in eV
c         nmaxp:  Total number of energy points
c         limit:  The fraction of the maximum value, below which the SRF is considered 0
c
c     Output arguments: 
c         fSInd:  For each final energy, the first energy bin where the SRF above the limit
c         fEInd:  For each final energy, the last energy bin where the SRF above the limit
c         ! check:  For each initial energy, the integrated SRF.
c

      implicit none
      integer nmaxp
      integer eini, efin
      real*8 problimit(nmaxp)
      real*8 prob(nmaxp,nmaxp) ! , df(nmaxp)
      !real*8 check(nmaxp)
      real*8 limit
      integer maxind(nmaxp)
      integer iSInd(nmaxp), iEInd(nmaxp)
      integer fSInd(nmaxp), fEInd(nmaxp)
      integer fLen(nmaxp)
      
      maxind(:) = 1
      
      fSInd(:) = nmaxp+1
      fEInd(:) = 0

c     For each initial energy, find the index of the most likely 
c     final energy.   (If fortran 95, could use MAXLOC)
      do eini=1, nmaxp
         do efin=2, nmaxp  
            if( ABS(prob(eini, efin))
     1          .GT. ABS(prob(eini, maxind(eini)))
     2         )then
               maxind(eini) = efin
            endif
         enddo

         iSInd(eini) = maxind(eini)
         iEInd(eini) = maxind(eini)
         problimit(eini) = ABS(prob(eini, maxind(eini)) * limit)
         
c        Decrement the value of the "start" index until reaches 1 or
c        the probability would be below the limit
         do efin = 1, maxind(eini)
            if(ABS(prob(eini, efin)) .GT. problimit(eini))then
               iSInd(eini) = efin
               exit
            endif
         enddo
 !        do while ( (iSInd(eini).GT.1) .AND.
 !    &              (prob(iSInd(eini)-1, eini) .GT. problimit(eini)) )
  !          iSInd(eini) =  iSInd(eini) - 1
   !      enddo
c        Increment the value of the "end" index until reaches nmaxp or
c        the probability would be below the limit
         do efin = nmaxp, maxind(eini),-1
            if(ABS(prob(eini, efin)) .GT. problimit(eini))then
               iEInd(eini) = efin
               exit
            endif
         enddo

 !        do while ( (iEInd(eini) .LT. nmaxp) .AND.
  !   &              (prob(iEInd(eini)+1, eini) .GT. problimit(eini)) )
  !          iEInd(eini) =  iEInd(eini) + 1
  !       enddo

c       !  Integrate over final energies which will be used
c        Also, mark the bounds of the final energies
         do efin=iSInd(eini), iEInd(eini)
            ! check(eini) = check(eini) + df(efin) * prob(efin, eini)
            fSInd(efin) = MIN( fSInd(efin), eini  )
            fEInd(efin) = MAX( fEInd(efin), eini  )
         enddo
         fLen = 1 + fEInd - fSInd
         do efin=1, nmaxp
            fLen(efin) = MAX(fLen(efin), 0)
         enddo
      enddo



      end subroutine
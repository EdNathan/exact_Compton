      module constants  
c
c     This module contains a number of physical and mathematical constants.
c     This is included so they don't need repeated definitions.
c      
      implicit none 
      double precision ergsev, ikbol, pi, mec2, isp, limit, sigma_t
      double precision hhh, skkk, ccc, smhy, smel, esu, eee, sigma
      parameter (ergsev  = 1.602197d-12     ) ! Convert eV to ergs
      parameter (ikbol   = 1.16d4           ) ! inverse of kbol (K * ev-1)
      parameter (pi      = 4.d0*datan(1.d0) ) ! pi number
      parameter (mec2    = 5.11d5           ) ! m_e c^2 (eV)
      parameter (isp     = 0.5641895835d0   ) ! 1/sqrt(pi)
      parameter (limit   = 1.d-3            ) ! Limit for the redistribution function
      parameter (sigma_t = 6.65d-25         ) ! Thomson cross section
      
      parameter (hhh     = 6.62620d-27      )
      parameter (skkk    = 1.38062d-16      )
      parameter (ccc     = 2.99793d+10      )
      parameter (smhy    = 1.67333d-24      )
      parameter (smel    = 9.10956d-28      )
      parameter (esu     = 4.80325d-10      )     
      parameter (eee     = 4.80325d-10      )
      parameter (sigma   = 5.66961d-05      )
      ! parameter (pi      = 3.141592654d0    ) 
      end module constants

      subroutine theta_from_temp(temp, theta)
            use constants
            implicit none
            double precision temp, theta
            theta = temp/ikbol/mec2
      end subroutine theta_from_temp
      
      subroutine temp_from_theta(theta, temp)
            use constants
            implicit none
            double precision temp, theta
            temp = theta*ikbol*mec2
      end subroutine temp_from_theta
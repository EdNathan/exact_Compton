      module constants  
c
c     This module contains a number of physical and mathematical constants.
c     This is included so they don't need repeated definitions.
c      
      implicit none 
      double precision ergsev, ikbol, pi, mec2, isp, limit, sigma_t
      parameter (ergsev  = 1.602197d-12     ) ! Convert eV to ergs
      parameter (ikbol   = 1.16d4           ) ! inverse of kbol (K * ev-1)
      parameter (pi      = 4.d0*datan(1.d0) ) ! pi number
      parameter (mec2    = 5.11d5           ) ! m_e c^2 (eV)
      parameter (isp     = 0.5641895835d0   ) ! 1/sqrt(pi)
      parameter (limit   = 1.d-3            ) ! Limit for the redistribution function
      parameter (sigma_t = 6.65d-25         ) ! Thomson cross section
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
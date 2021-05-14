! boundaries.F90
!
!=============================================================================

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine Christoffel_Boundaries( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT ierr
  CCTK_INT, parameter :: one = 1
  CCTK_INT bndsize

  if (compute_every .ge. 0 .and. MOD(cctk_iteration, compute_every) .eq. 0) then

     if (derivs_order == 6) then
        bndsize = 5
     else if (derivs_order == 4) then
        bndsize = 3
     else if (derivs_order == 2) then
        bndsize = 1
     else
        call CCTK_ERROR("derivs_order not yet implemented.")
     end if

     ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, bndsize, -one, &
          "Christoffel::Gammas", "flat")
     if (ierr < 0)                                                           &
          call CCTK_ERROR("Failed to register BC for Christoffel:Gammas!")

     if (save_dgab) then
        ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, bndsize, -one, &
             "Christoffel::derivatives_gab", "flat")
        if (ierr < 0)                                                           &
             call CCTK_ERROR("Failed to register BC for Christoffel:derivatives_gab!")
     end if

  end if
end subroutine Christoffel_Boundaries

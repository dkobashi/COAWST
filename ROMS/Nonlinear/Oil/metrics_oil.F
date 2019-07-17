#include "cppdefs.h"
      MODULE metrics_oil_mod

#if defined FLOAT_OIL
!
!=======================================================================
!                                                                      !
!========================================== Dmitry S Dukhovskoy      ===
!                                                                      !
!     Computes diagnsotics of oil 3D fields and reports                !
!                                                                      !
!=======================================================================
!
      implicit none

      PRIVATE
      PUBLIC  :: metrics_oil

      CONTAINS
!
!***********************************************************************
      SUBROUTINE metrics_oil (ng, tile, model, Lstr, Lend)
!***********************************************************************
!
      USE mod_param
      USE mod_grid
      USE mod_floats
      USE mod_stepping
      USE mod_oil_Eulvar
!
!  Imported variable declarations
!
      integer, intent(in) :: ng, tile, model
      integer, intent(in) :: Lstr, Lend
!
!  Local variable declaration
!
#include "tile.h"
!
!      print*,'metrics_oil, tile=',tile

      CALL metrics_oil_tile (ng, tile,                                  &
     &                       LBi, UBi, LBj, UBj,                        &
     &                       IminS, ImaxS, JminS, JmaxS,                &
     &                       Lstr, Lend, model,                         &
     &                       OIL3D(ng) % Coil,                          &
     &                       OIL3D(ng) % Doil,                          &
     &                       OIL3D(ng) % NFLT3D,                        &
     &                       nf(ng),                                    &
     &                       DRIFTER(ng) % track,                       &
     &                       DRIFTER(ng) % bounded)

      RETURN
      END SUBROUTINE metrics_oil
!
!***********************************************************************
      SUBROUTINE metrics_oil_tile (ng, tile,                            &
     &                             LBi, UBi, LBj, UBj,                  &
     &                             IminS, ImaxS, JminS, JmaxS,          &
     &                             Lstr, Lend, model,                   & 
     &                             Coil, Doil, NFLT3D, nf,              &
     &                             track, bounded)                   
!***********************************************************************
!
      USE mod_param
      USE mod_floats
      USE mod_scalars
      USE mod_grid
      USE mod_parallel
      USE mod_iounits
!
#ifdef DISTRIBUTE
      USE distribute_mod, ONLY : mp_reduce
#endif

      implicit none

!
!  Imported variable declarations
!
      integer, intent(in) :: ng, tile
      integer, intent(in) :: LBi, UBi, LBj, UBj
      integer, intent(in) :: IminS, ImaxS, JminS, JmaxS
      integer, intent(in) :: Lstr, Lend, model, nf
#ifdef ASSUMED_SHAPE
      real(r8), intent(in) :: track(:,0:,:)

      logical, intent(in) :: bounded(:)

      integer, intent(inout) :: NFLT3D(LBi:,LBj:,:)
      real(r8), intent(inout) :: Coil(LBi:,LBj:,:,:)
      real(r8), intent(inout) :: Doil(LBi:,LBj:,:)
#else
      real(r8), intent(in) :: track(NFV(ng),0:NFT,Nfloats(ng))

      logical, intent(in) :: bounded(Nfloats(ng))

      integer, intent(inout) :: NFLT3D(LBi:UBi,LBj:UBj,N(ng))
      real(r8), intent(inout) :: Coil(LBi:UBi,LBj:UBj,N(ng),Nocmp)
      real(r8), intent(inout) :: Doil(LBi:UBi,LBj:UBj,N(ng))
#endif
! 
! Local Variables
!
      integer :: nfltLagr, nfltEulr, nff
      integer :: i, j, k, ic, l
#ifdef DISTRIBUTE
      real(r8), dimension(2) :: buffer
      character (len=3), dimension(2) :: op_handle
#endif

#include "set_bounds.h"

#ifdef DISTRIBUTE !Added by DJ
      buffer = 0.0
!  Calculate total number of active floats in Lagr. and Eulerian spaces
      DO l=Lstr,Lend
        IF (bounded(l)) THEN
          buffer(1)=buffer(1)+1.0
          op_handle(1)='SUM'
        ENDIF
      ENDDO

      DO k=1,N(ng)
        DO j=Jstr-1,Jstr+1
          DO i=Istr-1,Istr+1
            nff=NFLT3D(i,j,k)
            buffer(2)=buffer(2)+REAL(nff)
            op_handle(2)='SUM'
          ENDDO
        ENDDO
      ENDDO
#endif

#ifdef DISTRIBUTE
      CALL mp_reduce (ng, model, 2, buffer, op_handle) ! note 2 should = dim of buffer & op_handle
      nfltLagr=buffer(1)
      nfltEulr=buffer(2)
#endif      

      IF (Master) THEN
        write(stdout,10) ng, nfltLagr, nfltEulr
      ENDIF

 10   FORMAT(/,' Lagrangian -> Eulerian Mapping (ng=',i2,') :',/,       &
     &         '================================',/,                    &
     &         'Total N released floats: Lagr=',i6,                     &
     &         '  N floats Mapped to Eulerian=',i6,/)


      END SUBROUTINE metrics_oil_tile

#endif
      END MODULE metrics_oil_mod


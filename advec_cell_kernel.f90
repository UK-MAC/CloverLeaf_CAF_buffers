!Crown Copyright 2012 AWE.
!
! This file is part of CloverLeaf.
!
! CloverLeaf is free software: you can redistribute it and/or modify it under 
! the terms of the GNU General Public License as published by the 
! Free Software Foundation, either version 3 of the License, or (at your option) 
! any later version.
!
! CloverLeaf is distributed in the hope that it will be useful, but 
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
! details.
!
! You should have received a copy of the GNU General Public License along with 
! CloverLeaf. If not, see http://www.gnu.org/licenses/.

!>  @brief Fortran cell advection kernel.
!>  @author Wayne Gaudin
!>  @details Performs a second order advective remap using van-Leer limiting
!>  with directional splitting.

MODULE advec_cell_kernel_module

USE clover_module

CONTAINS

SUBROUTINE advec_cell_kernel(x_min,       &
                             x_max,       &
                             y_min,       &
                             y_max,       &
                             dir,         &
                             sweep_number,&
                             vertexdx,    &
                             vertexdy,    &
                             volume,      &
                             density1,    &
                             energy1,     &
                             mass_flux_x, &
                             vol_flux_x,  &
                             mass_flux_y, &
                             vol_flux_y,  &
                             pre_vol,     &
                             post_vol,    &
                             pre_mass,    &
                             post_mass,   &
                             advec_vol,   &
                             post_ener,   &
                             ener_flux,   &
                             fields, depth, exchange)

  IMPLICIT NONE

  INTEGER :: x_min,x_max,y_min,y_max
  INTEGER :: sweep_number,dir
  INTEGER :: g_xdir=1,g_ydir=2

  INTEGER :: fields(:), depth
  LOGICAL :: exchange

  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: volume
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density1
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: energy1
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2) :: vol_flux_x
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3) :: vol_flux_y
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2) :: mass_flux_x
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3) :: mass_flux_y
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: pre_vol
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: post_vol
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: pre_mass
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: post_mass
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: advec_vol
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: post_ener
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: ener_flux

  REAL(KIND=8), DIMENSION(x_min-2:x_max+3) :: vertexdx
  REAL(KIND=8), DIMENSION(y_min-2:y_max+3) :: vertexdy


  IF(dir.EQ.g_xdir) THEN

        !compute first loop all boundaries 
        CALL xdir_loopblock1(x_min-2      , x_min+depth  , y_min-2      , y_max+2      , x_min, x_max, y_min, y_max, sweep_number, volume, vol_flux_x, vol_flux_y, pre_vol, post_vol)
        CALL xdir_loopblock1(x_max-depth  , x_max+2      , y_min-2      , y_max+2      , x_min, x_max, y_min, y_max, sweep_number, volume, vol_flux_x, vol_flux_y, pre_vol, post_vol)
        CALL xdir_loopblock1(x_min+depth+1, x_max-depth-1, y_min-2      , y_min+depth-1, x_min, x_max, y_min, y_max, sweep_number, volume, vol_flux_x, vol_flux_y, pre_vol, post_vol)
        CALL xdir_loopblock1(x_min+depth+1, x_max-depth-1, y_max-depth+1, y_max+2      , x_min, x_max, y_min, y_max, sweep_number, volume, vol_flux_x, vol_flux_y, pre_vol, post_vol)


        !compute second loop all boundaries
        CALL xdir_loopblock2(x_min        , x_min+depth, y_min        , y_max        , x_min, x_max, y_min, y_max, density1, energy1, vol_flux_x, mass_flux_x, pre_vol, ener_flux, vertexdx)
        CALL xdir_loopblock2(x_max-depth+1, x_max+2    , y_min        , y_max        , x_min, x_max, y_min, y_max, density1, energy1, vol_flux_x, mass_flux_x, pre_vol, ener_flux, vertexdx)
        CALL xdir_loopblock2(x_min+depth+1, x_max-depth, y_min        , y_min+depth-1, x_min, x_max, y_min, y_max, density1, energy1, vol_flux_x, mass_flux_x, pre_vol, ener_flux, vertexdx)
        CALL xdir_loopblock2(x_min+depth+1, x_max-depth, y_max-depth+1, y_max        , x_min, x_max, y_min, y_max, density1, energy1, vol_flux_x, mass_flux_x, pre_vol, ener_flux, vertexdx)


        !compute thirs loop all boundaries
        CALL xdir_loopblock3(x_min        , x_min+depth-1, y_min        , y_max        , x_min, x_max, y_min, y_max,                    &
                             density1, energy1, vol_flux_x, mass_flux_x, pre_vol, pre_mass, post_mass, advec_vol, post_ener, ener_flux)
        CALL xdir_loopblock3(x_max-depth+1, x_max        , y_min        , y_max        , x_min, x_max, y_min, y_max,                    &
                             density1, energy1, vol_flux_x, mass_flux_x, pre_vol, pre_mass, post_mass, advec_vol, post_ener, ener_flux)
        CALL xdir_loopblock3(x_min+depth  , x_max-depth  , y_min        , y_min+depth-1, x_min, x_max, y_min, y_max,                    &
                             density1, energy1, vol_flux_x, mass_flux_x, pre_vol, pre_mass, post_mass, advec_vol, post_ener, ener_flux)
        CALL xdir_loopblock3(x_min+depth  , x_max-depth  , y_max-depth+1, y_max        , x_min, x_max, y_min, y_max,                    &
                             density1, energy1, vol_flux_x, mass_flux_x, pre_vol, pre_mass, post_mass, advec_vol, post_ener, ener_flux)

        !execute send and recvs
        CALL clover_exchange_send_async(depth, fields)

        !compute remain centre
        CALL xdir_loopblock1(x_min+depth+1, x_max-depth-1, y_min+depth, y_max-depth, x_min, x_max, y_min, y_max,                        &
                             sweep_number, volume, vol_flux_x, vol_flux_y, pre_vol, post_vol)
        CALL xdir_loopblock2(x_min+depth+1, x_max-depth  , y_min+depth, y_max-depth, x_min, x_max, y_min, y_max,                        &
                             density1, energy1, vol_flux_x, mass_flux_x, pre_vol, ener_flux, vertexdx)
        CALL xdir_loopblock3(x_min+depth  , x_max-depth  , y_min+depth, y_max-depth, x_min, x_max, y_min, y_max,                        &
                             density1, energy1, vol_flux_x, mass_flux_x, pre_vol, pre_mass, post_mass, advec_vol, post_ener, ener_flux)

#ifdef LOCAL_SYNC
        sync images( chunks(1)%imageNeighbours )
#else
        sync all
#endif

        !execite WAITALL and unpack receives
        CALL clover_exchange_receive_async(depth, fields)

  ELSEIF(dir.EQ.g_ydir) THEN

        !compute 1st loop all boundaries
        CALL ydir_loopblock1(x_min-2      , x_min+depth-1, y_min-2    , y_max+2    , x_min, x_max, y_min, y_max, sweep_number, volume, vol_flux_x, vol_flux_y, pre_vol, post_vol)
        CALL ydir_loopblock1(x_max-depth+1, x_max+2      , y_min-2    , y_max+2    , x_min, x_max, y_min, y_max, sweep_number, volume, vol_flux_x, vol_flux_y, pre_vol, post_vol)
        CALL ydir_loopblock1(x_min+depth  , x_max-depth  , y_min-2    , y_min+depth, x_min, x_max, y_min, y_max, sweep_number, volume, vol_flux_x, vol_flux_y, pre_vol, post_vol)
        CALL ydir_loopblock1(x_min+depth  , x_max-depth  , y_max-depth, y_max+2    , x_min, x_max, y_min, y_max, sweep_number, volume, vol_flux_x, vol_flux_y, pre_vol, post_vol)


        !compute second loop all boundaries
        CALL ydir_loopblock2(x_min        , x_min+depth-1, y_min        , y_max+2    , x_min, x_max, y_min, y_max, density1, energy1, vol_flux_y, mass_flux_y, pre_vol, ener_flux, vertexdy)
        CALL ydir_loopblock2(x_max-depth+1, x_max        , y_min        , y_max+2    , x_min, x_max, y_min, y_max, density1, energy1, vol_flux_y, mass_flux_y, pre_vol, ener_flux, vertexdy)
        CALL ydir_loopblock2(x_min+depth  , x_max-depth  , y_min        , y_min+depth, x_min, x_max, y_min, y_max, density1, energy1, vol_flux_y, mass_flux_y, pre_vol, ener_flux, vertexdy)
        CALL ydir_loopblock2(x_min+depth  , x_max-depth  , y_max-depth+1, y_max+2    , x_min, x_max, y_min, y_max, density1, energy1, vol_flux_y, mass_flux_y, pre_vol, ener_flux, vertexdy)


        !compute third loop all boundaries
        CALL ydir_loopblock3(x_min        , x_min+depth-1, y_min        , y_max        , x_min, x_max, y_min, y_max,                    &
                             density1, energy1, vol_flux_y, mass_flux_y, pre_vol, post_vol, pre_mass, post_mass, advec_vol, post_ener, ener_flux)
        CALL ydir_loopblock3(x_max-depth+1, x_max        , y_min        , y_max        , x_min, x_max, y_min, y_max,                    &
                             density1, energy1, vol_flux_y, mass_flux_y, pre_vol, post_vol, pre_mass, post_mass, advec_vol, post_ener, ener_flux)
        CALL ydir_loopblock3(x_min+depth  , x_max-depth  , y_min        , y_min+depth-1, x_min, x_max, y_min, y_max,                    &
                             density1, energy1, vol_flux_y, mass_flux_y, pre_vol, post_vol, pre_mass, post_mass, advec_vol, post_ener, ener_flux)
        CALL ydir_loopblock3(x_min+depth  , x_max-depth  , y_max-depth+1, y_max        , x_min, x_max, y_min, y_max,                    &
                             density1, energy1, vol_flux_y, mass_flux_y, pre_vol, post_vol, pre_mass, post_mass, advec_vol, post_ener, ener_flux)

        !execute send and recvs
        CALL clover_exchange_send_async(depth, fields)

        !compute remaining centre
        CALL ydir_loopblock1(x_min+depth, x_max-depth, y_min+depth+1, y_max-depth-1, x_min, x_max, y_min, y_max,                        &
                             sweep_number, volume, vol_flux_x, vol_flux_y, pre_vol, post_vol)
        CALL ydir_loopblock2(x_min+depth, x_max-depth, y_min+depth+1, y_max-depth  , x_min, x_max, y_min, y_max,                        &
                             density1, energy1, vol_flux_y, mass_flux_y, pre_vol, ener_flux, vertexdy)
        CALL ydir_loopblock3(x_min+depth, x_max-depth, y_min+depth  , y_max-depth  , x_min, x_max, y_min, y_max,                        &
                             density1, energy1, vol_flux_y, mass_flux_y, pre_vol, post_vol, pre_mass, post_mass, advec_vol, post_ener, ener_flux)

#ifdef LOCAL_SYNC
        sync images( chunks(1)%imageNeighbours )
#else
        sync all
#endif

        !execute WAITALL and unpack receives
        CALL clover_exchange_receive_async(depth, fields)

  ENDIF

END SUBROUTINE advec_cell_kernel

SUBROUTINE xdir_loopblock1(j_start, j_end, k_start, k_end, x_min, x_max, y_min, y_max, sweep_number, volume, vol_flux_x, vol_flux_y, pre_vol, post_vol)

    IMPLICIT NONE

    INTEGER :: x_min,x_max,y_min,y_max
    INTEGER :: j_start, j_end, k_start, k_end
    INTEGER :: sweep_number

    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: volume
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2) :: vol_flux_x
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3) :: vol_flux_y
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: pre_vol
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: post_vol

    INTEGER :: j,k

    IF(sweep_number.EQ.1)THEN
      DO k=k_start,k_end
        DO j=j_start,j_end
          pre_vol(j,k)=volume(j,k)+(vol_flux_x(j+1,k  )-vol_flux_x(j,k)+vol_flux_y(j  ,k+1)-vol_flux_y(j,k))
          post_vol(j,k)=pre_vol(j,k)-(vol_flux_x(j+1,k  )-vol_flux_x(j,k))
        ENDDO
      ENDDO 
    ELSE
      DO k=k_start,k_end
        DO j=j_start,j_end
          pre_vol(j,k)=volume(j,k)+vol_flux_x(j+1,k)-vol_flux_x(j,k)
          post_vol(j,k)=volume(j,k)
        ENDDO
      ENDDO 
    ENDIF

END SUBROUTINE xdir_loopblock1

SUBROUTINE xdir_loopblock2(j_start, j_end, k_start, k_end, x_min, x_max, y_min, y_max, density1, energy1, vol_flux_x, mass_flux_x, pre_vol, ener_flux, vertexdx)

    IMPLICIT NONE

    INTEGER :: x_min,x_max,y_min,y_max
    INTEGER :: j_start, j_end, k_start, k_end

    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density1
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: energy1
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2) :: vol_flux_x
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2) :: mass_flux_x
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: pre_vol
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: ener_flux

    REAL(KIND=8), DIMENSION(x_min-2:x_max+3) :: vertexdx

    INTEGER :: j,k,upwind,donor,downwind,dif

    REAL(KIND=8) :: wind,sigma,sigmat,sigmav,sigmam,sigma3,sigma4
    REAL(KIND=8) :: diffuw,diffdw,limiter
    REAL(KIND=8) :: one_by_six=1.0_8/6.0_8

    DO k=k_start,k_end
      DO j=j_start,j_end

        IF(vol_flux_x(j,k).GT.0.0)THEN
          upwind   =j-2
          donor    =j-1
          downwind =j
          dif      =donor
        ELSE
          upwind   =MIN(j+1,x_max+2)
          donor    =j
          downwind =j-1
          dif      =upwind
        ENDIF

        sigmat=ABS(vol_flux_x(j,k))/pre_vol(donor,k)
        sigma3=(1.0_8+sigmat)*(vertexdx(j)/vertexdx(dif))
        sigma4=2.0_8-sigmat

        sigma=sigmat
        sigmav=sigmat

        diffuw=density1(donor,k)-density1(upwind,k)
        diffdw=density1(downwind,k)-density1(donor,k)
        wind=1.0_8
        IF(diffdw.LE.0.0) wind=-1.0_8
        IF(diffuw*diffdw.GT.0.0)THEN
          limiter=(1.0_8-sigmav)*wind*MIN(ABS(diffuw),ABS(diffdw)&
              ,one_by_six*(sigma3*ABS(diffuw)+sigma4*ABS(diffdw)))
        ELSE
          limiter=0.0
        ENDIF
        mass_flux_x(j,k)=vol_flux_x(j,k)*(density1(donor,k)+limiter)

        sigmam=ABS(mass_flux_x(j,k))/(density1(donor,k)*pre_vol(donor,k))
        diffuw=energy1(donor,k)-energy1(upwind,k)
        diffdw=energy1(downwind,k)-energy1(donor,k)
        wind=1.0_8
        IF(diffdw.LE.0.0) wind=-1.0_8
        IF(diffuw*diffdw.GT.0.0)THEN
          limiter=(1.0_8-sigmam)*wind*MIN(ABS(diffuw),ABS(diffdw)&
              ,one_by_six*(sigma3*ABS(diffuw)+sigma4*ABS(diffdw)))
        ELSE
          limiter=0.0
        ENDIF

        ener_flux(j,k)=mass_flux_x(j,k)*(energy1(donor,k)+limiter)

      ENDDO
    ENDDO

END SUBROUTINE xdir_loopblock2

SUBROUTINE xdir_loopblock3(j_start, j_end, k_start, k_end, x_min, x_max, y_min, y_max, density1, energy1, vol_flux_x, mass_flux_x, pre_vol, pre_mass, post_mass, advec_vol, post_ener, ener_flux)

    IMPLICIT NONE

    INTEGER :: x_min,x_max,y_min,y_max
    INTEGER :: j_start, j_end, k_start, k_end

    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density1
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: energy1
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2) :: vol_flux_x
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2) :: mass_flux_x
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: pre_vol
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: pre_mass
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: post_mass
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: advec_vol
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: post_ener
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: ener_flux

    INTEGER :: j,k

    DO k=k_start,k_end
      DO j=j_start,j_end
        pre_mass(j,k)=density1(j,k)*pre_vol(j,k)
        post_mass(j,k)=pre_mass(j,k)+mass_flux_x(j,k)-mass_flux_x(j+1,k)
        post_ener(j,k)=(energy1(j,k)*pre_mass(j,k)+ener_flux(j,k)-ener_flux(j+1,k))/post_mass(j,k)
        advec_vol(j,k)=pre_vol(j,k)+vol_flux_x(j,k)-vol_flux_x(j+1,k)
        density1(j,k)=post_mass(j,k)/advec_vol(j,k)
        energy1(j,k)=post_ener(j,k)
      ENDDO
    ENDDO

END SUBROUTINE xdir_loopblock3


SUBROUTINE ydir_loopblock1(j_start, j_end, k_start, k_end, x_min, x_max, y_min, y_max, sweep_number, volume, vol_flux_x, vol_flux_y, pre_vol, post_vol)

    IMPLICIT NONE

    INTEGER :: x_min,x_max,y_min,y_max
    INTEGER :: j_start, j_end, k_start, k_end
    INTEGER :: sweep_number

    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: volume
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2) :: vol_flux_x
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3) :: vol_flux_y
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: pre_vol
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: post_vol

    INTEGER :: j,k

    IF(sweep_number.EQ.1)THEN
      DO k=k_start,k_end
        DO j=j_start,j_end
          pre_vol(j,k)=volume(j,k)+(vol_flux_y(j  ,k+1)-vol_flux_y(j,k)+vol_flux_x(j+1,k  )-vol_flux_x(j,k))
          post_vol(j,k)=pre_vol(j,k)-(vol_flux_y(j  ,k+1)-vol_flux_y(j,k))
        ENDDO
      ENDDO
    ELSE
      DO k=k_start,k_end
        DO j=j_start,j_end
          pre_vol(j,k)=volume(j,k)+vol_flux_y(j  ,k+1)-vol_flux_y(j,k)
          post_vol(j,k)=volume(j,k)
        ENDDO
      ENDDO
    ENDIF

END SUBROUTINE ydir_loopblock1

SUBROUTINE ydir_loopblock2(j_start, j_end, k_start, k_end, x_min, x_max, y_min, y_max, density1, energy1, vol_flux_y, mass_flux_y, pre_vol, ener_flux, vertexdy)

    IMPLICIT NONE

    INTEGER :: x_min,x_max,y_min,y_max
    INTEGER :: j_start, j_end, k_start, k_end

    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density1
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: energy1
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3) :: vol_flux_y
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3) :: mass_flux_y
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: pre_vol
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: ener_flux

    REAL(KIND=8), DIMENSION(y_min-2:y_max+3) :: vertexdy

    INTEGER :: j,k,upwind,donor,downwind,dif

    REAL(KIND=8) :: wind,sigma,sigmat,sigmav,sigmam,sigma3,sigma4
    REAL(KIND=8) :: diffuw,diffdw,limiter
    REAL(KIND=8) :: one_by_six=1.0_8/6.0_8

    DO k=k_start,k_end
      DO j=j_start,j_end

        IF(vol_flux_y(j,k).GT.0.0)THEN
          upwind   =k-2
          donor    =k-1
          downwind =k
          dif      =donor
        ELSE
          upwind   =MIN(k+1,y_max+2)
          donor    =k
          downwind =k-1
          dif      =upwind
        ENDIF

        sigmat=ABS(vol_flux_y(j,k))/pre_vol(j,donor)
        sigma3=(1.0_8+sigmat)*(vertexdy(k)/vertexdy(dif))
        sigma4=2.0_8-sigmat

        sigma=sigmat
        sigmav=sigmat

        diffuw=density1(j,donor)-density1(j,upwind)
        diffdw=density1(j,downwind)-density1(j,donor)
        wind=1.0_8
        IF(diffdw.LE.0.0) wind=-1.0_8
        IF(diffuw*diffdw.GT.0.0)THEN
          limiter=(1.0_8-sigmav)*wind*MIN(ABS(diffuw),ABS(diffdw)&
              ,one_by_six*(sigma3*ABS(diffuw)+sigma4*ABS(diffdw)))
        ELSE
          limiter=0.0
        ENDIF
        mass_flux_y(j,k)=vol_flux_y(j,k)*(density1(j,donor)+limiter)

        sigmam=ABS(mass_flux_y(j,k))/(density1(j,donor)*pre_vol(j,donor))
        diffuw=energy1(j,donor)-energy1(j,upwind)
        diffdw=energy1(j,downwind)-energy1(j,donor)
        wind=1.0_8
        IF(diffdw.LE.0.0) wind=-1.0_8
        IF(diffuw*diffdw.GT.0.0)THEN
          limiter=(1.0_8-sigmam)*wind*MIN(ABS(diffuw),ABS(diffdw)&
              ,one_by_six*(sigma3*ABS(diffuw)+sigma4*ABS(diffdw)))
        ELSE
          limiter=0.0
        ENDIF
        ener_flux(j,k)=mass_flux_y(j,k)*(energy1(j,donor)+limiter)

      ENDDO
    ENDDO

END SUBROUTINE ydir_loopblock2

SUBROUTINE ydir_loopblock3(j_start, j_end, k_start, k_end, x_min, x_max, y_min, y_max,          &
                           density1, energy1, vol_flux_y, mass_flux_y, pre_vol, post_vol,       &
                           pre_mass, post_mass, advec_vol, post_ener, ener_flux)

    IMPLICIT NONE

    INTEGER :: x_min,x_max,y_min,y_max
    INTEGER :: j_start, j_end, k_start, k_end

    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density1
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: energy1
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3) :: vol_flux_y
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3) :: mass_flux_y
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: pre_vol
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: post_vol
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: pre_mass
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: post_mass
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: advec_vol
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: post_ener
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: ener_flux

    INTEGER :: j,k

    DO k=k_start,k_end
      DO j=j_start,j_end
        pre_mass(j,k)=density1(j,k)*pre_vol(j,k)
        post_mass(j,k)=pre_mass(j,k)+mass_flux_y(j,k)-mass_flux_y(j,k+1)
        post_ener(j,k)=(energy1(j,k)*pre_mass(j,k)+ener_flux(j,k)-ener_flux(j,k+1))/post_mass(j,k)
        advec_vol(j,k)=pre_vol(j,k)+vol_flux_y(j,k)-vol_flux_y(j,k+1)
        density1(j,k)=post_mass(j,k)/advec_vol(j,k)
        energy1(j,k)=post_ener(j,k)
      ENDDO
    ENDDO

END SUBROUTINE ydir_loopblock3

END MODULE advec_cell_kernel_module


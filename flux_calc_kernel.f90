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

!>  @brief Fortran flux kernel.
!>  @author Wayne Gaudin
!>  @details The edge volume fluxes are calculated based on the velocity fields.

MODULE flux_calc_kernel_module

USE clover_module

CONTAINS

SUBROUTINE flux_calc_kernel(x_min,x_max,y_min,y_max,dt,              &
                            xarea,                           &
                            yarea,                           &
                            xvel0,                           &
                            yvel0,                           &
                            xvel1,                           &
                            yvel1,                           &
                            vol_flux_x,                      &
                            vol_flux_y,                      &
                            fields, depth, exchange          )

    IMPLICIT NONE

    INTEGER       :: x_min, x_max, y_min, y_max
    REAL(KIND=8) :: dt
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2) :: xarea
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3) :: yarea
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xvel0,yvel0
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xvel1,yvel1
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2) :: vol_flux_x
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3) :: vol_flux_y

    INTEGER :: fields(:), depth
    LOGICAL :: exchange

    !compute left boundary
    CALL real_vol_flux_x(x_min, x_min+depth-1, y_min, y_max  , x_min, x_max, y_min, y_max, dt, vol_flux_x, xarea, xvel0, xvel1)
    CALL real_vol_flux_y(x_min, x_min+depth-1, y_min, y_max+1, x_min, x_max, y_min, y_max, dt, vol_flux_y, yarea, yvel0, yvel1)

    !compute right boundary
    CALL real_vol_flux_x(x_max-depth+1, x_max+1,y_min, y_max  , x_min, x_max, y_min, y_max, dt, vol_flux_x, xarea, xvel0, xvel1)
    CALL real_vol_flux_y(x_max-depth+1, x_max  ,y_min, y_max+1, x_min, x_max, y_min, y_max, dt, vol_flux_y, yarea, yvel0, yvel1)

    !compute bottom boundary
    CALL real_vol_flux_x(x_min+depth, x_max-depth, y_min, y_min+depth-1, x_min, x_max, y_min, y_max, dt, vol_flux_x, xarea, xvel0, xvel1)
    CALL real_vol_flux_y(x_min+depth, x_max-depth, y_min, y_min+depth-1, x_min, x_max, y_min, y_max, dt, vol_flux_y, yarea, yvel0, yvel1)

    !compute top boundary 
    CALL real_vol_flux_x(x_min+depth, x_max-depth, y_max-depth+1, y_max  ,x_min, x_max, y_min, y_max, dt, vol_flux_x, xarea, xvel0, xvel1)
    CALL real_vol_flux_y(x_min+depth, x_max-depth, y_max-depth+1, y_max+1,x_min, x_max, y_min, y_max, dt, vol_flux_y, yarea, yvel0, yvel1)

    !send all the fields and post receives before doing the remainder of the compute
    CALL clover_exchange_send_async(depth, fields)

    CALL real_vol_flux_x(x_min+depth, x_max-depth, y_min+depth, y_max-depth, x_min, x_max, y_min, y_max, dt, vol_flux_x, xarea, xvel0, xvel1)
    CALL real_vol_flux_y(x_min+depth, x_max-depth, y_min+depth, y_max-depth, x_min, x_max, y_min, y_max, dt, vol_flux_y, yarea, yvel0, yvel1)

#ifdef LOCAL_SYNC
    sync images( chunks(parallel%task+1)%imageNeighbours )
#else
    sync all
#endif

    !call method in clover to do to waitall and unpack the buffers, could prepost the buffers at some point
    CALL clover_exchange_receive_async(depth, fields)

END SUBROUTINE flux_calc_kernel

SUBROUTINE real_vol_flux_x(j_start,j_end,k_start,k_end,x_min,x_max,y_min,y_max,dt,vol_flux_x,xarea,xvel0,xvel1)

    IMPLICIT NONE

    INTEGER       :: x_min, x_max, y_min, y_max
    INTEGER       :: j_start, j_end, k_start, k_end
    REAL(KIND=8) :: dt
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2) :: xarea
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xvel0
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xvel1
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+2) :: vol_flux_x

    INTEGER :: j,k

    DO k=k_start,k_end
      DO j=j_start,j_end
        vol_flux_x(j,k)=0.25_8*dt*xarea(j,k)                  &
                       *(xvel0(j,k)+xvel0(j,k+1)+xvel1(j,k)+xvel1(j,k+1))
      ENDDO
    ENDDO

END SUBROUTINE real_vol_flux_x

SUBROUTINE real_vol_flux_y(j_start,j_end,k_start,k_end,x_min,x_max,y_min,y_max,dt,vol_flux_y,yarea,yvel0,yvel1)

    IMPLICIT NONE

    INTEGER       :: x_min, x_max, y_min, y_max
    INTEGER       :: j_start, j_end, k_start, k_end
    REAL(KIND=8) :: dt
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3) :: yarea
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: yvel0
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: yvel1
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+3) :: vol_flux_y

    INTEGER :: j,k

    DO k=k_start,k_end
      DO j=j_start,j_end
        vol_flux_y(j,k)=0.25_8*dt*yarea(j,k)                  &
                       *(yvel0(j,k)+yvel0(j+1,k)+yvel1(j,k)+yvel1(j+1,k))
      ENDDO
    ENDDO

END SUBROUTINE real_vol_flux_y

END MODULE flux_calc_kernel_module

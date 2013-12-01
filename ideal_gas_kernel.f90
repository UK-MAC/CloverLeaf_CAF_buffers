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

!>  @brief Fortran ideal gas kernel.
!>  @author Wayne Gaudin
!>  @details Calculates the pressure and sound speed for the mesh chunk using
!>  the ideal gas equation of state, with a fixed gamma of 1.4.

MODULE ideal_gas_kernel_module

USE clover_module

CONTAINS

SUBROUTINE ideal_gas_kernel(x_min,x_max,y_min,y_max,                &
                            density,                                &
                            energy,                                 &
                            pressure,                               &
                            soundspeed,                             &                              
                            fields, depth, exchange                 )

    IMPLICIT NONE

    INTEGER :: x_min,x_max,y_min,y_max
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: energy
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: pressure
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: soundspeed

    INTEGER :: j,k, depth
    INTEGER :: fields(:)

    REAL(KIND=8) :: sound_speed_squared,v,pressurebyenergy,pressurebyvolume

    LOGICAL :: exchange

    IF (exchange) THEN
        !compute left boundary
        CALL ideal_gas_kernel_real(x_min,x_min+depth-1,y_min,y_max,x_min,x_max,y_min,y_max,density,energy,pressure,soundspeed)
        !compute right boundary
        CALL ideal_gas_kernel_real(x_max-depth+1,x_max,y_min,y_max,x_min,x_max,y_min,y_max,density,energy,pressure,soundspeed)
        !compute bottom boundary
        CALL ideal_gas_kernel_real(x_min+depth,x_max-depth,y_min,y_min+depth-1,x_min,x_max,y_min,y_max,density,energy,pressure,soundspeed)
        !compute top boundary
        CALL ideal_gas_kernel_real(x_min+depth,x_max-depth,y_max-depth+1,y_max,x_min,x_max,y_min,y_max,density,energy,pressure,soundspeed)

        !send all the fields and post receives before doing the remainder of the compute
        CALL clover_exchange_send_async(depth, fields)

        ! compute remaining cells in the centre, overlapping with comms
        CALL ideal_gas_kernel_real(x_min+depth,x_max-depth,y_min+depth,y_max-depth,x_min,x_max,y_min,y_max,density,energy,pressure,soundspeed)

        !call method in clover to do to waitall and unpack the buffers, could prepost the buffers at some point
        CALL clover_exchange_receive_async(depth, fields)

    ELSE
        ! compute the entire ideal gas domain as no comms needed
        CALL ideal_gas_kernel_real(x_min,x_max,y_min,y_max,x_min,x_max,y_min,y_max,density,energy,pressure,soundspeed)
    ENDIF

END SUBROUTINE ideal_gas_kernel

SUBROUTINE ideal_gas_kernel_real(j_start,j_end,k_start,k_end,x_min,x_max,y_min,y_max,density,energy,pressure,soundspeed)

    IMPLICIT NONE 

    INTEGER :: x_min,x_max,y_min,y_max
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: energy
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: pressure
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: soundspeed

    REAL(KIND=8) :: sound_speed_squared,v,pressurebyenergy,pressurebyvolume

    INTEGER :: j, k, j_start, j_end, k_start, k_end

    DO k=k_start, k_end
        DO j=j_start, j_end
            v=1.0_8/density(j,k)  
            pressure(j,k)=(1.4_8-1.0_8)*density(j,k)*energy(j,k)
            pressurebyenergy=(1.4_8-1.0_8)*density(j,k)
            pressurebyvolume=-density(j,k)*pressure(j,k)
            sound_speed_squared=v*v*(pressure(j,k)*pressurebyenergy-pressurebyvolume)
            soundspeed(j,k)=SQRT(sound_speed_squared)
        ENDDO
    ENDDO

END SUBROUTINE ideal_gas_kernel_real

END MODULE ideal_gas_kernel_module

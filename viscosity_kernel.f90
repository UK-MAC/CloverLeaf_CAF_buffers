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

!>  @brief Fortran viscosity kernel.
!>  @author Wayne Gaudin
!>  @details Calculates an artificial viscosity using the Wilkin's method to
!>  smooth out shock front and prevent oscillations around discontinuities.
!>  Only cells in compression will have a non-zero value.

MODULE viscosity_kernel_module

USE clover_module

CONTAINS

SUBROUTINE viscosity_kernel(x_min,x_max,y_min,y_max,    &
                            celldx,celldy,              &
                            density0,                   &
                            pressure,                   &
                            viscosity,                  &
                            xvel0,                      &
                            yvel0,                      &
                            fields, depth, exchange     )


    IMPLICIT NONE

    INTEGER     :: x_min,x_max,y_min,y_max
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2)                     :: celldx
    REAL(KIND=8), DIMENSION(y_min-2:y_max+2)                     :: celldy
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2)     :: density0
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2)     :: pressure
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2)     :: viscosity
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3)     :: xvel0,yvel0

    INTEGER       :: depth, fields(:)
    LOGICAL       :: exchange

    IF (exchange) THEN

        CALL viscosity_kernel_real(x_min        ,x_min+depth-1,y_min        ,y_max        ,x_min,x_max,y_min,y_max,viscosity,xvel0,yvel0,celldx,celldy,pressure,density0)
        CALL viscosity_kernel_real(x_max-depth+1,x_max        ,y_min        ,y_max        ,x_min,x_max,y_min,y_max,viscosity,xvel0,yvel0,celldx,celldy,pressure,density0)
        CALL viscosity_kernel_real(x_min+depth  ,x_max-depth  ,y_min        ,y_min+depth-1,x_min,x_max,y_min,y_max,viscosity,xvel0,yvel0,celldx,celldy,pressure,density0)
        CALL viscosity_kernel_real(x_min+depth  ,x_max-depth  ,y_max-depth+1,y_max        ,x_min,x_max,y_min,y_max,viscosity,xvel0,yvel0,celldx,celldy,pressure,density0)

        !send all the fields and post receives before doing the remainder of the compute
        CALL clover_exchange_send_async(depth, fields)

        ! compute remaining cells in the centre, overlapping with comms
        CALL viscosity_kernel_real(x_min+depth,x_max-depth,y_min+depth,y_max-depth,x_min,x_max,y_min,y_max,viscosity,xvel0,yvel0,celldx,celldy,pressure,density0)

#ifdef LOCAL_SYNC
        sync images( chunks(1)%imageNeighbours )
#else
        sync all
#endif

        !call method in clover to do to waitall and unpack the buffers, could prepost the buffers at some point
        CALL clover_exchange_receive_async(depth, fields)

    ELSE
        ! compute entire viscosity domain as no comms required 
        CALL viscosity_kernel_real(x_min,x_max,y_min,y_max,x_min,x_max,y_min,y_max,viscosity,xvel0,yvel0,celldx,celldy,pressure,density0)
    ENDIF

END SUBROUTINE viscosity_kernel




SUBROUTINE viscosity_kernel_real(j_start,j_end,k_start,k_end,x_min,x_max,y_min,y_max,viscosity,xvel0,yvel0,celldx,celldy,pressure,density0)

    IMPLICIT NONE

    INTEGER     :: x_min,x_max,y_min,y_max
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2)                     :: celldx
    REAL(KIND=8), DIMENSION(y_min-2:y_max+2)                     :: celldy
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2)     :: density0
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2)     :: pressure
    REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2)     :: viscosity
    REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3)     :: xvel0,yvel0

    INTEGER       :: j,k
    INTEGER       :: j_start,j_end,k_start,k_end
    REAL(KIND=8)  :: ugrad,vgrad,grad2,pgradx,pgrady,pgradx2,pgrady2,grad     &
                    ,ygrad,pgrad,xgrad,div,strain2,limiter,dirx,diry

    DO k=k_start,k_end
      DO j=j_start,j_end
        ugrad=(xvel0(j+1,k  )+xvel0(j+1,k+1))-(xvel0(j  ,k  )+xvel0(j  ,k+1))

        vgrad=(yvel0(j  ,k+1)+yvel0(j+1,k+1))-(yvel0(j  ,k  )+yvel0(j+1,k  ))

        div = (celldx(j)*(ugrad)+  celldy(k)*(vgrad))

        strain2 = 0.5_8*(xvel0(j,  k+1) + xvel0(j+1,k+1)-xvel0(j  ,k  )-xvel0(j+1,k  ))/celldy(k) &
                + 0.5_8*(yvel0(j+1,k  ) + yvel0(j+1,k+1)-yvel0(j  ,k  )-yvel0(j  ,k+1))/celldx(j)

        pgradx=(pressure(j+1,k)-pressure(j-1,k))/(celldx(j)+celldx(j+1))
        pgrady=(pressure(j,k+1)-pressure(j,k-1))/(celldy(k)+celldy(k+1))

        pgradx2 = pgradx*pgradx
        pgrady2 = pgrady*pgrady

        limiter = ((0.5_8*(ugrad)/celldx(j))*pgradx2+(0.5_8*(vgrad)/celldy(k))*pgrady2+strain2*pgradx*pgrady)  &
                /MAX(pgradx2+pgrady2,1.0e-16_8)

        IF ((limiter.GT.0.0).OR.(div.GE.0.0))THEN
          viscosity(j,k) = 0.0
        ELSE
          dirx=1.0_8
          IF(pgradx.LT.0.0) dirx=-1.0_8
          pgradx = dirx*MAX(1.0e-16_8,ABS(pgradx))
          diry=1.0_8
          IF(pgradx.LT.0.0) diry=-1.0_8
          pgrady = diry*MAX(1.0e-16_8,ABS(pgrady))
          pgrad = SQRT(pgradx**2+pgrady**2)
          xgrad = ABS(celldx(j)*pgrad/pgradx)
          ygrad = ABS(celldy(k)*pgrad/pgrady)
          grad  = MIN(xgrad,ygrad)
          grad2 = grad*grad

          viscosity(j,k)=2.0_8*density0(j,k)*grad2*limiter*limiter
        ENDIF

      ENDDO
    ENDDO

END SUBROUTINE viscosity_kernel_real

END MODULE viscosity_kernel_module

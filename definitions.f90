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

!>  @brief Holds the high level Fortran data types
!>  @author Wayne Gaudin
!>  @details The high level data types used to store the mesh and field data
!>  are defined here.
!>
!>  Also the global variables used for defining the input and controlling the
!>  scheme are defined here.

MODULE definitions_module

   USE data_module
   
   IMPLICIT NONE

   TYPE state_type
      LOGICAL            :: defined

      REAL(KIND=8)       :: density          &
                           ,energy           &
                           ,xvel             &
                           ,yvel

      INTEGER            :: geometry

      REAL(KIND=8)       :: xmin               &
                           ,xmax               &
                           ,ymin               &
                           ,ymax               &
                           ,radius
   END TYPE state_type

   TYPE(state_type), ALLOCATABLE             :: states(:)
   INTEGER                                   :: number_of_states

   TYPE grid_type
     REAL(KIND=8)       :: xmin            &
                          ,ymin            &
                          ,xmax            &
                          ,ymax
                     
     INTEGER            :: x_cells              &
                          ,y_cells
   END TYPE grid_type

   INTEGER      :: step

   LOGICAL      :: advect_x

   INTEGER      :: error_condition

   INTEGER      :: test_problem
   LOGICAL      :: complete

   LOGICAL      :: use_fortran_kernels
   LOGICAL      :: use_C_kernels
   LOGICAL      :: use_OA_kernels

   LOGICAL      :: use_vector_loops ! Some loops work better in serial depending on the hardware

   REAL(KIND=8) :: end_time

   INTEGER      :: end_step

   REAL(KIND=8) :: dtold          &
                  ,dt             &
                  ,time           &
                  ,dtinit         &
                  ,dtmin          &
                  ,dtmax          &
                  ,dtrise         &
                  ,dtu_safe       &
                  ,dtv_safe       &
                  ,dtc_safe       &
                  ,dtdiv_safe     &
                  ,dtc            &
                  ,dtu            &
                  ,dtv            &
                  ,dtdiv

   INTEGER      :: visit_frequency   &
                  ,summary_frequency

   INTEGER         :: jdt,kdt

   TYPE field_type
     REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: density0,density1
     REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: energy0,energy1
     REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: pressure
     REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: viscosity
     REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: soundspeed
     REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: xvel0,xvel1
     REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: yvel0,yvel1
     REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: vol_flux_x,mass_flux_x
     REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: vol_flux_y,mass_flux_y
     REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: work_array1 !node_flux, stepbymass, volume_change, pre_vol
     REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: work_array2 !node_mass_post, post_vol
     REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: work_array3 !node_mass_pre,pre_mass
     REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: work_array4 !advec_vel, post_mass
     REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: work_array5 !mom_flux, advec_vol
     REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: work_array6 !pre_vol, post_ener
     REAL(KIND=8),    DIMENSION(:,:), ALLOCATABLE :: work_array7 !post_vol, ener_flux

     INTEGER         :: left            &
                       ,right           &
                       ,bottom          &
                       ,top             &
                       ,left_boundary   &
                       ,right_boundary  &
                       ,bottom_boundary &
                       ,top_boundary

     INTEGER         :: x_min  &
                       ,y_min  &
                       ,x_max  &
                       ,y_max

     REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: cellx    &
                                                 ,celly    &
                                                 ,vertexx  &
                                                 ,vertexy  &
                                                 ,celldx   &
                                                 ,celldy   &
                                                 ,vertexdx &
                                                 ,vertexdy

     REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: volume  &
                                                 ,xarea   &
                                                 ,yarea

   END TYPE field_type
   
   TYPE chunk_type

     INTEGER         :: task   !caf task minus 1

     INTEGER         :: chunk_neighbours(8) ! Chunks, not tasks, so we can overload in the future

     ! Idealy, create an array to hold the buffers for each field so a commuincation only needs
     !  one send and one receive per face, rather than per field.
     ! If chunks are overloaded, i.e. more chunks than tasks, might need to pack for a task to task comm 
     !  rather than a chunk to chunk comm. See how performance is at high core counts before deciding
     REAL(KIND=8),ALLOCATABLE:: density0_left_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: density0_left_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: density0_right_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: density0_right_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: density0_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: density0_bottom_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: density0_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: density0_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: density0_left_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: density0_left_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: density0_right_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: density0_right_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: density0_right_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: density0_right_bottom_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: density0_left_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: density0_left_bottom_rcv_buffer(:)

     REAL(KIND=8),ALLOCATABLE:: density1_left_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: density1_left_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: density1_right_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: density1_right_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: density1_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: density1_bottom_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: density1_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: density1_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: density1_left_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: density1_left_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: density1_right_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: density1_right_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: density1_right_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: density1_right_bottom_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: density1_left_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: density1_left_bottom_rcv_buffer(:)

     REAL(KIND=8),ALLOCATABLE:: energy0_left_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: energy0_left_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: energy0_right_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: energy0_right_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: energy0_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: energy0_bottom_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: energy0_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: energy0_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: energy0_left_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: energy0_left_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: energy0_right_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: energy0_right_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: energy0_right_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: energy0_right_bottom_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: energy0_left_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: energy0_left_bottom_rcv_buffer(:)

     REAL(KIND=8),ALLOCATABLE:: energy1_left_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: energy1_left_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: energy1_right_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: energy1_right_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: energy1_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: energy1_bottom_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: energy1_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: energy1_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: energy1_left_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: energy1_left_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: energy1_right_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: energy1_right_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: energy1_right_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: energy1_right_bottom_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: energy1_left_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: energy1_left_bottom_rcv_buffer(:)

     REAL(KIND=8),ALLOCATABLE:: pressure_left_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: pressure_left_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: pressure_right_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: pressure_right_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: pressure_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: pressure_bottom_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: pressure_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: pressure_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: pressure_left_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: pressure_left_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: pressure_right_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: pressure_right_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: pressure_right_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: pressure_right_bottom_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: pressure_left_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: pressure_left_bottom_rcv_buffer(:)

     REAL(KIND=8),ALLOCATABLE:: viscosity_left_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: viscosity_left_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: viscosity_right_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: viscosity_right_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: viscosity_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: viscosity_bottom_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: viscosity_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: viscosity_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: viscosity_left_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: viscosity_left_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: viscosity_right_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: viscosity_right_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: viscosity_right_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: viscosity_right_bottom_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: viscosity_left_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: viscosity_left_bottom_rcv_buffer(:)

     REAL(KIND=8),ALLOCATABLE:: soundspeed_left_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: soundspeed_left_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: soundspeed_right_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: soundspeed_right_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: soundspeed_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: soundspeed_bottom_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: soundspeed_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: soundspeed_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: soundspeed_left_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: soundspeed_left_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: soundspeed_right_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: soundspeed_right_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: soundspeed_right_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: soundspeed_right_bottom_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: soundspeed_left_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: soundspeed_left_bottom_rcv_buffer(:)

     REAL(KIND=8),ALLOCATABLE:: xvel0_left_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: xvel0_left_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: xvel0_right_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: xvel0_right_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: xvel0_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: xvel0_bottom_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: xvel0_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: xvel0_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: xvel0_left_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: xvel0_left_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: xvel0_right_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: xvel0_right_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: xvel0_right_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: xvel0_right_bottom_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: xvel0_left_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: xvel0_left_bottom_rcv_buffer(:)

     REAL(KIND=8),ALLOCATABLE:: xvel1_left_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: xvel1_left_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: xvel1_right_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: xvel1_right_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: xvel1_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: xvel1_bottom_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: xvel1_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: xvel1_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: xvel1_left_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: xvel1_left_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: xvel1_right_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: xvel1_right_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: xvel1_right_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: xvel1_right_bottom_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: xvel1_left_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: xvel1_left_bottom_rcv_buffer(:)

     REAL(KIND=8),ALLOCATABLE:: yvel0_left_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: yvel0_left_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: yvel0_right_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: yvel0_right_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: yvel0_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: yvel0_bottom_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: yvel0_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: yvel0_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: yvel0_left_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: yvel0_left_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: yvel0_right_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: yvel0_right_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: yvel0_right_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: yvel0_right_bottom_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: yvel0_left_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: yvel0_left_bottom_rcv_buffer(:)

     REAL(KIND=8),ALLOCATABLE:: yvel1_left_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: yvel1_left_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: yvel1_right_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: yvel1_right_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: yvel1_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: yvel1_bottom_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: yvel1_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: yvel1_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: yvel1_left_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: yvel1_left_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: yvel1_right_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: yvel1_right_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: yvel1_right_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: yvel1_right_bottom_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: yvel1_left_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: yvel1_left_bottom_rcv_buffer(:)

     REAL(KIND=8),ALLOCATABLE:: volflux_x_left_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: volflux_x_left_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: volflux_x_right_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: volflux_x_right_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: volflux_x_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: volflux_x_bottom_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: volflux_x_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: volflux_x_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: volflux_x_left_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: volflux_x_left_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: volflux_x_right_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: volflux_x_right_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: volflux_x_right_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: volflux_x_right_bottom_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: volflux_x_left_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: volflux_x_left_bottom_rcv_buffer(:)

     REAL(KIND=8),ALLOCATABLE:: volflux_y_left_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: volflux_y_left_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: volflux_y_right_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: volflux_y_right_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: volflux_y_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: volflux_y_bottom_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: volflux_y_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: volflux_y_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: volflux_y_left_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: volflux_y_left_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: volflux_y_right_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: volflux_y_right_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: volflux_y_right_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: volflux_y_right_bottom_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: volflux_y_left_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: volflux_y_left_bottom_rcv_buffer(:)

     REAL(KIND=8),ALLOCATABLE:: massflux_x_left_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: massflux_x_left_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: massflux_x_right_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: massflux_x_right_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: massflux_x_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: massflux_x_bottom_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: massflux_x_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: massflux_x_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: massflux_x_left_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: massflux_x_left_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: massflux_x_right_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: massflux_x_right_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: massflux_x_right_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: massflux_x_right_bottom_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: massflux_x_left_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: massflux_x_left_bottom_rcv_buffer(:)

     REAL(KIND=8),ALLOCATABLE:: massflux_y_left_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: massflux_y_left_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: massflux_y_right_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: massflux_y_right_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: massflux_y_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: massflux_y_bottom_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: massflux_y_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: massflux_y_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: massflux_y_left_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: massflux_y_left_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: massflux_y_right_top_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: massflux_y_right_top_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: massflux_y_right_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: massflux_y_right_bottom_rcv_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: massflux_y_left_bottom_snd_buffer(:)
     REAL(KIND=8),ALLOCATABLE:: massflux_y_left_bottom_rcv_buffer(:)


#ifdef LOCAL_SYNC
     INTEGER, ALLOCATABLE:: imageNeighbours(:)
#endif

     TYPE(field_type):: field

  END TYPE chunk_type


  TYPE(chunk_type),  ALLOCATABLE       :: chunks(:)[:]
  INTEGER                              :: number_of_chunks

  TYPE(grid_type)                      :: grid

END MODULE definitions_module

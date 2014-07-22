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

!>  @brief Communication Utilities
!>  @author Wayne Gaudin
!>  @details Contains all utilities required to run CloverLeaf in a distributed
!>  environment, including initialisation, mesh decompostion, reductions and
!>  halo exchange using explicit buffers.
!>
!>  Note the halo exchange is currently coded as simply as possible and no 
!>  optimisations have been implemented, such as post receives before sends or packing
!>  buffers with multiple data fields. This is intentional so the effect of these
!>  optimisations can be measured on large systems, as and when they are added.
!>
!>  Even without these modifications CloverLeaf weak scales well on moderately sized
!>  systems of the order of 10K cores.

MODULE clover_module

  USE data_module
  USE definitions_module

  IMPLICIT NONE

CONTAINS

SUBROUTINE clover_barrier

    sync all

END SUBROUTINE clover_barrier

SUBROUTINE clover_abort

    ERROR STOP

END SUBROUTINE clover_abort

SUBROUTINE clover_finalize

  INTEGER :: err

  CLOSE(g_out)
  CALL FLUSH(0)
  CALL FLUSH(6)
  CALL FLUSH(g_out)

END SUBROUTINE clover_finalize

SUBROUTINE clover_init_comms

  IMPLICIT NONE

  parallel%parallel=.TRUE.

  parallel%image = this_image()
  parallel%max_image = num_images()

  parallel%task = parallel%image - 1

  IF(parallel%image.EQ.1) THEN
    parallel%boss=.TRUE.
  ENDIF

  parallel%boss_task=0
  parallel%max_task=parallel%max_image

END SUBROUTINE clover_init_comms

SUBROUTINE clover_get_num_chunks(count)

  IMPLICIT NONE

  INTEGER :: count

! Should be changed so there can be more than one chunk per mpi task

  count=parallel%max_task

END SUBROUTINE clover_get_num_chunks

SUBROUTINE clover_decompose(x_cells,y_cells,left,right,bottom,top)

  ! This decomposes the mesh into a number of chunks.
  ! The number of chunks may be a multiple of the number of mpi tasks
  ! Doesn't always return the best split if there are few factors
  ! All factors need to be stored and the best picked. But its ok for now

  IMPLICIT NONE

  INTEGER :: x_cells,y_cells,left(:),right(:),top(:),bottom(:)
  INTEGER :: c,delta_x,delta_y

  REAL(KIND=8) :: mesh_ratio,factor_x,factor_y
  INTEGER  :: chunk_x,chunk_y,mod_x,mod_y,split_found

  INTEGER  :: cx,cy,chunk,add_x,add_y,add_x_prev,add_y_prev
#ifdef LOCAL_SYNC
  INTEGER :: numNeighbours,n
#endif

  ! 2D Decomposition of the mesh

  mesh_ratio=real(x_cells)/real(y_cells)

  chunk_x=number_of_chunks
  chunk_y=1

  split_found=0 ! Used to detect 1D decomposition
  DO c=1,number_of_chunks
    IF (MOD(number_of_chunks,c).EQ.0) THEN
      factor_x=number_of_chunks/real(c)
      factor_y=c
      !Compare the factor ratio with the mesh ratio
      IF(factor_x/factor_y.LE.mesh_ratio) THEN
        chunk_y=c
        chunk_x=number_of_chunks/c
        split_found=1
        EXIT
      ENDIF
    ENDIF
  ENDDO

  IF(split_found.EQ.0.OR.chunk_y.EQ.number_of_chunks) THEN ! Prime number or 1D decomp detected
    IF(mesh_ratio.GE.1.0) THEN
      chunk_x=number_of_chunks
      chunk_y=1
    ELSE
      chunk_x=1
      chunk_y=number_of_chunks
    ENDIF
  ENDIF

  delta_x=x_cells/chunk_x
  delta_y=y_cells/chunk_y
  mod_x=MOD(x_cells,chunk_x)
  mod_y=MOD(y_cells,chunk_y)

  ! Set up chunk mesh ranges and chunk connectivity

    add_x_prev=0
    add_y_prev=0
    chunk=1
    DO cy=1,chunk_y
        DO cx=1,chunk_x
            add_x=0
            add_y=0
            IF(cx.LE.mod_x)add_x=1
            IF(cy.LE.mod_y)add_y=1

            IF (chunk .EQ. parallel%task+1) THEN

                left(1)   = (cx-1)*delta_x+1+add_x_prev
                right(1)  = left(1)+delta_x-1+add_x
                bottom(1) = (cy-1)*delta_y+1+add_y_prev
                top(1)    = bottom(1)+delta_y-1+add_y

                chunks(1)%chunk_neighbours(chunk_left)=chunk_x*(cy-1)+cx-1
                chunks(1)%chunk_neighbours(chunk_right)=chunk_x*(cy-1)+cx+1
                chunks(1)%chunk_neighbours(chunk_bottom)=chunk_x*(cy-2)+cx
                chunks(1)%chunk_neighbours(chunk_top)=chunk_x*(cy)+cx

                IF(cx.EQ.1)       chunks(1)%chunk_neighbours(chunk_left)=external_face
                IF(cx.EQ.chunk_x) chunks(1)%chunk_neighbours(chunk_right)=external_face
                IF(cy.EQ.1)       chunks(1)%chunk_neighbours(chunk_bottom)=external_face
                IF(cy.EQ.chunk_y) chunks(1)%chunk_neighbours(chunk_top)=external_face

#ifdef LOCAL_SYNC
                numNeighbours=0
                IF (chunks(1)%chunk_neighbours(chunk_left).NE.external_face) THEN
                    numNeighbours = numNeighbours +1
                ENDIF
                IF (chunks(1)%chunk_neighbours(chunk_right).NE.external_face) THEN
                   numNeighbours = numNeighbours +1
                ENDIF
                IF (chunks(1)%chunk_neighbours(chunk_top).NE.external_face) THEN
                   numNeighbours = numNeighbours +1
                ENDIF
                IF (chunks(1)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
                   numNeighbours = numNeighbours +1
                ENDIF
                ALLOCATE(chunks(1)%imageNeighbours(numNeighbours))

                !caf:may need to update this when multiple chunks per image so that the image is recorded correctly 
                IF (numNeighbours > 0) THEN
                   n=1
                   IF (chunks(1)%chunk_neighbours(chunk_left).NE.external_face) THEN
                      chunks(1)%imageNeighbours(n) = chunks(1)%chunk_neighbours(chunk_left)
                      n=n+1
                   ENDIF
                   IF (chunks(1)%chunk_neighbours(chunk_right).NE.external_face) THEN
                      chunks(1)%imageNeighbours(n) = chunks(1)%chunk_neighbours(chunk_right)
                      n=n+1
                   ENDIF
                   IF (chunks(1)%chunk_neighbours(chunk_top).NE.external_face) THEN
                      chunks(1)%imageNeighbours(n) = chunks(1)%chunk_neighbours(chunk_top)
                      n=n+1
                   ENDIF
                   IF (chunks(1)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
                      chunks(1)%imageNeighbours(n) = chunks(1)%chunk_neighbours(chunk_bottom)
                      n=n+1
                   ENDIF
                ENDIF
#endif

            ENDIF

            IF(cx.LE.mod_x)add_x_prev=add_x_prev+1
            chunk=chunk+1
        ENDDO
        add_x_prev=0
        IF(cy.LE.mod_y)add_y_prev=add_y_prev+1
    ENDDO

  IF(parallel%boss)THEN
    WRITE(g_out,*)
    WRITE(g_out,*)"Mesh ratio of ",mesh_ratio
    WRITE(g_out,*)"Decomposing the mesh into ",chunk_x," by ",chunk_y," chunks"
    WRITE(g_out,*)
  ENDIF

END SUBROUTINE clover_decompose

SUBROUTINE clover_allocate_buffers(chunk)

  IMPLICIT NONE

  INTEGER      :: chunk
  
  ! Unallocated buffers for external boundaries caused issues on some systems so they are now
  !  all allocated
  IF(parallel%task.EQ.chunks(chunk)%task)THEN
    !IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
      ALLOCATE(left_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(left_rcv_buffer(2*(chunks(chunk)%field%y_max+5))[*])
    !ENDIF
    !IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
      ALLOCATE(right_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(right_rcv_buffer(2*(chunks(chunk)%field%y_max+5))[*])
    !ENDIF
    !IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
      ALLOCATE(bottom_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(bottom_rcv_buffer(2*(chunks(chunk)%field%x_max+5))[*])
    !ENDIF
    !IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
      ALLOCATE(top_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(top_rcv_buffer(2*(chunks(chunk)%field%x_max+5))[*])
    !ENDIF
  ENDIF

END SUBROUTINE clover_allocate_buffers

SUBROUTINE clover_exchange(fields,depth)

  IMPLICIT NONE

  INTEGER      :: fields(:),depth,location_of_tasks_chunk

  ! Assuming 1 patch per task, this will be changed
  ! Also, not packing all fields for each communication, doing one at a time

    location_of_tasks_chunk = 1

  IF(fields(FIELD_DENSITY0).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,                                               &
                                 chunks(location_of_tasks_chunk)%field%density0,      &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_DENSITY1).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,                                               &
                                 chunks(location_of_tasks_chunk)%field%density1,      &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_ENERGY0).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,                                               &
                                 chunks(location_of_tasks_chunk)%field%energy0,       &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_ENERGY1).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,                                               &
                                 chunks(location_of_tasks_chunk)%field%energy1,       &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_PRESSURE).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,                                               &
                                 chunks(location_of_tasks_chunk)%field%pressure,      &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_VISCOSITY).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,                                               &
                                 chunks(location_of_tasks_chunk)%field%viscosity,     &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,                                               &
                                 chunks(location_of_tasks_chunk)%field%soundspeed,    &
                                 depth,CELL_DATA)
  ENDIF

  IF(fields(FIELD_XVEL0).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,                                               &
                                 chunks(location_of_tasks_chunk)%field%xvel0,         &
                                 depth,VERTEX_DATA)
  ENDIF

  IF(fields(FIELD_XVEL1).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,                                               &
                                 chunks(location_of_tasks_chunk)%field%xvel1,         &
                                 depth,VERTEX_DATA)
  ENDIF

  IF(fields(FIELD_YVEL0).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,                                               &
                                 chunks(location_of_tasks_chunk)%field%yvel0,         &
                                 depth,VERTEX_DATA)
  ENDIF

  IF(fields(FIELD_YVEL1).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,                                               &
                                 chunks(location_of_tasks_chunk)%field%yvel1,         &
                                 depth,VERTEX_DATA)
  ENDIF

  IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,                                               &
                                 chunks(location_of_tasks_chunk)%field%vol_flux_x,    &
                                 depth,X_FACE_DATA)
  ENDIF

  IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,                                               &
                                 chunks(location_of_tasks_chunk)%field%vol_flux_y,    &
                                 depth,Y_FACE_DATA)
  ENDIF

  IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,                                               &
                                 chunks(location_of_tasks_chunk)%field%mass_flux_x,   &
                                 depth,X_FACE_DATA)
  ENDIF

  IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
    CALL clover_exchange_message(location_of_tasks_chunk,                                               &
                                 chunks(location_of_tasks_chunk)%field%mass_flux_y,   &
                                 depth,Y_FACE_DATA)
  ENDIF


END SUBROUTINE clover_exchange

SUBROUTINE clover_exchange_message(chunk,field,                            &
                                   depth,field_type)

    USE pack_kernel_module

    IMPLICIT NONE

    REAL(KIND=8) :: field(-1:,-1:) ! This seems to work for any type of mesh data

    INTEGER      :: chunk,depth,field_type

    INTEGER      :: size,err,x_inc,y_inc
    INTEGER      :: receiver,sender
    INTEGER      :: left_neighbour_chunk, right_neighbour_chunk, bottom_neighbour_chunk, top_neighbour_chunk

    ! Field type will either be cell, vertex, x_face or y_face to get the message limits correct

    ! I am packing my own buffers. I am sure this could be improved with MPI data types
    !  but this will do for now

    ! I am also sending buffers to chunks with the same task id for now.
    ! This can be improved in the future but at the moment there is just 1 chunk per task anyway

    ! The tag will be a function of the sending chunk and the face it is coming from
    !  like chunk 6 sending the left face

    ! No open mp in here either. May be beneficial will packing and unpacking in the future, though I am not sure.

    ! Change this so it will allow more than 1 chunk per task


    ! Pack and send

    ! These array modifications still need to be added on, plus the donor data location changes as in update_halo
    IF(field_type.EQ.CELL_DATA) THEN
      x_inc=0
      y_inc=0
    ENDIF
    IF(field_type.EQ.VERTEX_DATA) THEN
      x_inc=1
      y_inc=1
    ENDIF
    IF(field_type.EQ.X_FACE_DATA) THEN
      x_inc=1
      y_inc=0
    ENDIF
    IF(field_type.EQ.Y_FACE_DATA) THEN
      x_inc=0
      y_inc=1
    ENDIF

    ! Pack real data into buffers
    IF(parallel%task.EQ.chunks(chunk)%task) THEN
        size=(1+(chunks(chunk)%field%y_max+y_inc+depth)-(chunks(chunk)%field%y_min-depth))*depth

        IF(use_fortran_kernels) THEN
          CALL pack_left_right_buffers(chunks(chunk)%field%x_min,chunks(chunk)%field%x_max, &
                                       chunks(chunk)%field%y_min,chunks(chunk)%field%y_max, &
                                       chunks(chunk)%chunk_neighbours(chunk_left),          &
                                       chunks(chunk)%chunk_neighbours(chunk_right),         &
                                       external_face,                                       &
                                       x_inc,y_inc,depth,size,                              &
                                       field,left_snd_buffer,right_snd_buffer)
        ELSEIF(use_C_kernels)THEN
          CALL pack_left_right_buffers_c(chunks(chunk)%field%x_min,chunks(chunk)%field%x_max, &
                                         chunks(chunk)%field%y_min,chunks(chunk)%field%y_max, &
                                         chunks(chunk)%chunk_neighbours(chunk_left),          &
                                         chunks(chunk)%chunk_neighbours(chunk_right),         &
                                         external_face,                                       &
                                         x_inc,y_inc,depth,size,                              &
                                         field,left_snd_buffer,right_snd_buffer)
        ENDIF

        IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
          receiver=chunks(chunk)%chunk_neighbours(chunk_left)
          right_rcv_buffer(1:size)[receiver] = left_snd_buffer(1:size)
        ENDIF

        IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
          receiver = chunks(chunk)%chunk_neighbours(chunk_right)
          left_rcv_buffer(1:size)[receiver] = right_snd_buffer(1:size)
        ENDIF
    ENDIF

  ! Wait for the messages
#ifdef LOCAL_SYNC
    sync images( chunks(chunk)%imageNeighbours )
#else
    sync all
#endif


    ! Unpack buffers in halo cells
    IF(parallel%task.EQ.chunks(chunk)%task) THEN

        IF(use_fortran_kernels) THEN
          CALL unpack_left_right_buffers(chunks(chunk)%field%x_min,chunks(chunk)%field%x_max, &
                                         chunks(chunk)%field%y_min,chunks(chunk)%field%y_max, &
                                         chunks(chunk)%chunk_neighbours(chunk_left),          &
                                         chunks(chunk)%chunk_neighbours(chunk_right),         &
                                         external_face,                                       &
                                         x_inc,y_inc,depth,size,                              &
                                         field,left_rcv_buffer,right_rcv_buffer)
        ELSEIF(use_C_kernels)THEN
          CALL unpack_left_right_buffers_c(chunks(chunk)%field%x_min,chunks(chunk)%field%x_max, &
                                           chunks(chunk)%field%y_min,chunks(chunk)%field%y_max, &
                                           chunks(chunk)%chunk_neighbours(chunk_left),          &
                                           chunks(chunk)%chunk_neighbours(chunk_right),         &
                                           external_face,                                       &
                                           x_inc,y_inc,depth,size,                              &
                                           field,left_rcv_buffer,right_rcv_buffer)
        ENDIF
    ENDIF


    IF(parallel%task.EQ.chunks(chunk)%task) THEN
        size=(1+(chunks(chunk)%field%x_max+x_inc+depth)-(chunks(chunk)%field%x_min-depth))*depth
        IF(use_fortran_kernels) THEN
          CALL pack_top_bottom_buffers(chunks(chunk)%field%x_min,chunks(chunk)%field%x_max, &
                                       chunks(chunk)%field%y_min,chunks(chunk)%field%y_max, &
                                       chunks(chunk)%chunk_neighbours(chunk_bottom),        &
                                       chunks(chunk)%chunk_neighbours(chunk_top),           &
                                       external_face,                                       &
                                       x_inc,y_inc,depth,size,                              &
                                       field,bottom_snd_buffer,top_snd_buffer)
        ELSEIF(use_C_kernels)THEN
          CALL pack_top_bottom_buffers_c(chunks(chunk)%field%x_min,chunks(chunk)%field%x_max, &
                                         chunks(chunk)%field%y_min,chunks(chunk)%field%y_max, &
                                         chunks(chunk)%chunk_neighbours(chunk_bottom),        &
                                         chunks(chunk)%chunk_neighbours(chunk_top),           &
                                         external_face,                                       &
                                         x_inc,y_inc,depth,size,                              &
                                         field,bottom_snd_buffer,top_snd_buffer)
        ENDIF

        IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
            receiver=chunks(chunk)%chunk_neighbours(chunk_bottom)
            top_rcv_buffer(1:size)[receiver] = bottom_snd_buffer(1:size)
        ENDIF

        IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
            receiver=chunks(chunk)%chunk_neighbours(chunk_top)
            bottom_rcv_buffer(1:size)[receiver] = top_snd_buffer(1:size)
        ENDIF
    ENDIF


  ! Wait for the messages
#ifdef LOCAL_SYNC
    sync images( chunks(chunk)%imageNeighbours )
#else
    sync all
#endif


    ! Unpack buffers in halo cells
    IF(parallel%task.EQ.chunks(chunk)%task) THEN

        IF(parallel%task.EQ.chunks(chunk)%task) THEN
            IF(use_fortran_kernels) THEN
              CALL unpack_top_bottom_buffers(chunks(chunk)%field%x_min,chunks(chunk)%field%x_max, &
                                             chunks(chunk)%field%y_min,chunks(chunk)%field%y_max, &
                                             chunks(chunk)%chunk_neighbours(chunk_bottom),        &
                                             chunks(chunk)%chunk_neighbours(chunk_top),           &
                                             external_face,                                       &
                                             x_inc,y_inc,depth,size,                              &
                                             field,bottom_rcv_buffer,top_rcv_buffer)
            ELSEIF(use_C_kernels)THEN
              CALL unpack_top_bottom_buffers_c(chunks(chunk)%field%x_min,chunks(chunk)%field%x_max, &
                                               chunks(chunk)%field%y_min,chunks(chunk)%field%y_max, &
                                               chunks(chunk)%chunk_neighbours(chunk_bottom),        &
                                               chunks(chunk)%chunk_neighbours(chunk_top),           &
                                               external_face,                                       &
                                               x_inc,y_inc,depth,size,                              &
                                               field,bottom_rcv_buffer,top_rcv_buffer)
            ENDIF
        ENDIF
    ENDIF

END SUBROUTINE clover_exchange_message

SUBROUTINE clover_sum(value)

  ! Only sums to the master

  IMPLICIT NONE

  REAL(KIND=8) :: value

  REAL(KIND=8) :: total

  total=value

  CALL CO_SUM(value, total, result_image=1)

  value=total

END SUBROUTINE clover_sum

SUBROUTINE clover_min(value)

  IMPLICIT NONE

  REAL(KIND=8) :: value

  REAL(KIND=8) :: minimum

  minimum=value

  CALL CO_MIN(value, minimum)

  value=minimum

END SUBROUTINE clover_min

SUBROUTINE clover_max(value)

  IMPLICIT NONE

  REAL(KIND=8) :: value

  REAL(KIND=8) :: maximum

  INTEGER :: err

  maximum=value

  !CALL MPI_ALLREDUCE(value,maximum,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,err)
  CALL CO_MAX(value, maximum)

  value=maximum

END SUBROUTINE clover_max

SUBROUTINE clover_allgather(value,values)

    IMPLICIT NONE

    REAL(KIND=8) :: value

    REAL(KIND=8) :: values(parallel%max_task)

    INTEGER :: err

    values(1)=value ! Just to ensure it will work in serial

    totals(parallel%image)[1] = value

END SUBROUTINE clover_allgather

SUBROUTINE clover_check_error(error)

  IMPLICIT NONE

  INTEGER :: error

  INTEGER :: maximum

  maximum=error

  CALL CO_MAX(error, maximum)

  error=maximum

END SUBROUTINE clover_check_error


END MODULE clover_module

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
  USE pack_kernel_module

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
      left(chunk)=(cx-1)*delta_x+1+add_x_prev
      right(chunk)=left(chunk)+delta_x-1+add_x
      bottom(chunk)=(cy-1)*delta_y+1+add_y_prev
      top(chunk)=bottom(chunk)+delta_y-1+add_y

      chunks(chunk)%chunk_neighbours(chunk_left)=chunk_x*(cy-1)+cx-1
      chunks(chunk)%chunk_neighbours(chunk_right)=chunk_x*(cy-1)+cx+1
      chunks(chunk)%chunk_neighbours(chunk_bottom)=chunk_x*(cy-2)+cx
      chunks(chunk)%chunk_neighbours(chunk_top)=chunk_x*(cy)+cx

      chunks(chunk)%chunk_neighbours(CHUNK_LEFT_TOP) = chunk_x*cy+cx-1
      chunks(chunk)%chunk_neighbours(CHUNK_LEFT_BOTTOM) = chunk_x*(cy-2)+cx-1
      chunks(chunk)%chunk_neighbours(CHUNK_RIGHT_TOP) = chunk_x*cy+cx+1
      chunks(chunk)%chunk_neighbours(CHUNK_RIGHT_BOTTOM) = chunk_x*(cy-2)+cx+1

      IF(cx.EQ.1) THEN 
        chunks(chunk)%chunk_neighbours(chunk_left)=external_face
        chunks(chunk)%chunk_neighbours(CHUNK_LEFT_TOP)=external_face
        chunks(chunk)%chunk_neighbours(CHUNK_LEFT_BOTTOM)=external_face
      ENDIF
      IF(cx.EQ.chunk_x) THEN 
        chunks(chunk)%chunk_neighbours(chunk_right)=external_face
        chunks(chunk)%chunk_neighbours(CHUNK_RIGHT_TOP)=external_face
        chunks(chunk)%chunk_neighbours(CHUNK_RIGHT_BOTTOM)=external_face
      ENDIF
      IF(cy.EQ.1) THEN 
       chunks(chunk)%chunk_neighbours(chunk_bottom)=external_face
        chunks(chunk)%chunk_neighbours(CHUNK_LEFT_BOTTOM)=external_face
        chunks(chunk)%chunk_neighbours(CHUNK_RIGHT_BOTTOM)=external_face
      ENDIF
      IF(cy.EQ.chunk_y) THEN 
        chunks(chunk)%chunk_neighbours(chunk_top)=external_face
        chunks(chunk)%chunk_neighbours(CHUNK_LEFT_TOP)=external_face
        chunks(chunk)%chunk_neighbours(CHUNK_RIGHT_TOP)=external_face
      ENDIF


#ifdef LOCAL_SYNC
      numNeighbours=0
      IF (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
          numNeighbours = numNeighbours +1
      ENDIF
      IF (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
         numNeighbours = numNeighbours +1
      ENDIF
      IF (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
         numNeighbours = numNeighbours +1
      ENDIF
      IF (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
         numNeighbours = numNeighbours +1
      ENDIF
      IF (chunks(chunk)%chunk_neighbours(chunk_left_top).NE.external_face) THEN
          numNeighbours = numNeighbours +1
      ENDIF
      IF (chunks(chunk)%chunk_neighbours(chunk_right_top).NE.external_face) THEN
         numNeighbours = numNeighbours +1
      ENDIF
      IF (chunks(chunk)%chunk_neighbours(chunk_right_bottom).NE.external_face) THEN
         numNeighbours = numNeighbours +1
      ENDIF
      IF (chunks(chunk)%chunk_neighbours(chunk_left_bottom).NE.external_face) THEN
         numNeighbours = numNeighbours +1
      ENDIF
      ALLOCATE(chunks(chunk)%imageNeighbours(numNeighbours))

      !caf:may need to update this when multiple chunks per image so that the image is recorded correctly 
      IF (numNeighbours > 0) THEN
         n=1
         IF (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
            chunks(chunk)%imageNeighbours(n) = chunks(chunk)%chunk_neighbours(chunk_left)
            n=n+1
         ENDIF
         IF (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
            chunks(chunk)%imageNeighbours(n) = chunks(chunk)%chunk_neighbours(chunk_right)
            n=n+1
         ENDIF
         IF (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
            chunks(chunk)%imageNeighbours(n) = chunks(chunk)%chunk_neighbours(chunk_top)
            n=n+1
         ENDIF
         IF (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
            chunks(chunk)%imageNeighbours(n) = chunks(chunk)%chunk_neighbours(chunk_bottom)
            n=n+1
         ENDIF
         IF (chunks(chunk)%chunk_neighbours(chunk_left_top).NE.external_face) THEN
            chunks(chunk)%imageNeighbours(n) = chunks(chunk)%chunk_neighbours(chunk_left_top)
            n=n+1
         ENDIF
         IF (chunks(chunk)%chunk_neighbours(chunk_right_top).NE.external_face) THEN
            chunks(chunk)%imageNeighbours(n) = chunks(chunk)%chunk_neighbours(chunk_right_top)
            n=n+1
         ENDIF
         IF (chunks(chunk)%chunk_neighbours(chunk_right_bottom).NE.external_face) THEN
            chunks(chunk)%imageNeighbours(n) = chunks(chunk)%chunk_neighbours(chunk_right_bottom)
            n=n+1
         ENDIF
         IF (chunks(chunk)%chunk_neighbours(chunk_left_bottom).NE.external_face) THEN
            chunks(chunk)%imageNeighbours(n) = chunks(chunk)%chunk_neighbours(chunk_left_bottom)
            n=n+1
         ENDIF
      ENDIF
#endif

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
      ALLOCATE(chunks(chunk)%density0_left_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%density0_left_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%density0_right_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%density0_right_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%density0_bottom_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%density0_bottom_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%density0_top_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%density0_top_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%density0_left_top_snd_buffer(4))
      ALLOCATE(chunks(chunk)%density0_left_top_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%density0_right_top_snd_buffer(4))
      ALLOCATE(chunks(chunk)%density0_right_top_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%density0_right_bottom_snd_buffer(4))
      ALLOCATE(chunks(chunk)%density0_right_bottom_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%density0_left_bottom_snd_buffer(4))
      ALLOCATE(chunks(chunk)%density0_left_bottom_rcv_buffer(4))

      ALLOCATE(chunks(chunk)%density1_left_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%density1_left_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%density1_right_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%density1_right_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%density1_bottom_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%density1_bottom_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%density1_top_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%density1_top_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%density1_left_top_snd_buffer(4))
      ALLOCATE(chunks(chunk)%density1_left_top_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%density1_right_top_snd_buffer(4))
      ALLOCATE(chunks(chunk)%density1_right_top_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%density1_right_bottom_snd_buffer(4))
      ALLOCATE(chunks(chunk)%density1_right_bottom_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%density1_left_bottom_snd_buffer(4))
      ALLOCATE(chunks(chunk)%density1_left_bottom_rcv_buffer(4))

      ALLOCATE(chunks(chunk)%energy0_left_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%energy0_left_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%energy0_right_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%energy0_right_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%energy0_bottom_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%energy0_bottom_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%energy0_top_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%energy0_top_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%energy0_left_top_snd_buffer(4))
      ALLOCATE(chunks(chunk)%energy0_left_top_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%energy0_right_top_snd_buffer(4))
      ALLOCATE(chunks(chunk)%energy0_right_top_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%energy0_right_bottom_snd_buffer(4))
      ALLOCATE(chunks(chunk)%energy0_right_bottom_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%energy0_left_bottom_snd_buffer(4))
      ALLOCATE(chunks(chunk)%energy0_left_bottom_rcv_buffer(4))

      ALLOCATE(chunks(chunk)%energy1_left_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%energy1_left_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%energy1_right_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%energy1_right_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%energy1_bottom_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%energy1_bottom_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%energy1_top_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%energy1_top_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%energy1_left_top_snd_buffer(4))
      ALLOCATE(chunks(chunk)%energy1_left_top_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%energy1_right_top_snd_buffer(4))
      ALLOCATE(chunks(chunk)%energy1_right_top_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%energy1_right_bottom_snd_buffer(4))
      ALLOCATE(chunks(chunk)%energy1_right_bottom_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%energy1_left_bottom_snd_buffer(4))
      ALLOCATE(chunks(chunk)%energy1_left_bottom_rcv_buffer(4))

      ALLOCATE(chunks(chunk)%pressure_left_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%pressure_left_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%pressure_right_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%pressure_right_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%pressure_bottom_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%pressure_bottom_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%pressure_top_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%pressure_top_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%pressure_left_top_snd_buffer(4))
      ALLOCATE(chunks(chunk)%pressure_left_top_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%pressure_right_top_snd_buffer(4))
      ALLOCATE(chunks(chunk)%pressure_right_top_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%pressure_right_bottom_snd_buffer(4))
      ALLOCATE(chunks(chunk)%pressure_right_bottom_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%pressure_left_bottom_snd_buffer(4))
      ALLOCATE(chunks(chunk)%pressure_left_bottom_rcv_buffer(4))

      ALLOCATE(chunks(chunk)%viscosity_left_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%viscosity_left_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%viscosity_right_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%viscosity_right_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%viscosity_bottom_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%viscosity_bottom_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%viscosity_top_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%viscosity_top_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%viscosity_left_top_snd_buffer(4))
      ALLOCATE(chunks(chunk)%viscosity_left_top_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%viscosity_right_top_snd_buffer(4))
      ALLOCATE(chunks(chunk)%viscosity_right_top_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%viscosity_right_bottom_snd_buffer(4))
      ALLOCATE(chunks(chunk)%viscosity_right_bottom_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%viscosity_left_bottom_snd_buffer(4))
      ALLOCATE(chunks(chunk)%viscosity_left_bottom_rcv_buffer(4))

      ALLOCATE(chunks(chunk)%soundspeed_left_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%soundspeed_left_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%soundspeed_right_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%soundspeed_right_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%soundspeed_bottom_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%soundspeed_bottom_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%soundspeed_top_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%soundspeed_top_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%soundspeed_left_top_snd_buffer(4))
      ALLOCATE(chunks(chunk)%soundspeed_left_top_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%soundspeed_right_top_snd_buffer(4))
      ALLOCATE(chunks(chunk)%soundspeed_right_top_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%soundspeed_right_bottom_snd_buffer(4))
      ALLOCATE(chunks(chunk)%soundspeed_right_bottom_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%soundspeed_left_bottom_snd_buffer(4))
      ALLOCATE(chunks(chunk)%soundspeed_left_bottom_rcv_buffer(4))

      ALLOCATE(chunks(chunk)%xvel0_left_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%xvel0_left_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%xvel0_right_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%xvel0_right_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%xvel0_bottom_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%xvel0_bottom_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%xvel0_top_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%xvel0_top_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%xvel0_left_top_snd_buffer(4))
      ALLOCATE(chunks(chunk)%xvel0_left_top_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%xvel0_right_top_snd_buffer(4))
      ALLOCATE(chunks(chunk)%xvel0_right_top_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%xvel0_right_bottom_snd_buffer(4))
      ALLOCATE(chunks(chunk)%xvel0_right_bottom_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%xvel0_left_bottom_snd_buffer(4))
      ALLOCATE(chunks(chunk)%xvel0_left_bottom_rcv_buffer(4))

      ALLOCATE(chunks(chunk)%xvel1_left_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%xvel1_left_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%xvel1_right_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%xvel1_right_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%xvel1_bottom_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%xvel1_bottom_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%xvel1_top_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%xvel1_top_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%xvel1_left_top_snd_buffer(4))
      ALLOCATE(chunks(chunk)%xvel1_left_top_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%xvel1_right_top_snd_buffer(4))
      ALLOCATE(chunks(chunk)%xvel1_right_top_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%xvel1_right_bottom_snd_buffer(4))
      ALLOCATE(chunks(chunk)%xvel1_right_bottom_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%xvel1_left_bottom_snd_buffer(4))
      ALLOCATE(chunks(chunk)%xvel1_left_bottom_rcv_buffer(4))

      ALLOCATE(chunks(chunk)%yvel0_left_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%yvel0_left_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%yvel0_right_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%yvel0_right_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%yvel0_bottom_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%yvel0_bottom_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%yvel0_top_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%yvel0_top_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%yvel0_left_top_snd_buffer(4))
      ALLOCATE(chunks(chunk)%yvel0_left_top_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%yvel0_right_top_snd_buffer(4))
      ALLOCATE(chunks(chunk)%yvel0_right_top_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%yvel0_right_bottom_snd_buffer(4))
      ALLOCATE(chunks(chunk)%yvel0_right_bottom_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%yvel0_left_bottom_snd_buffer(4))
      ALLOCATE(chunks(chunk)%yvel0_left_bottom_rcv_buffer(4))

      ALLOCATE(chunks(chunk)%yvel1_left_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%yvel1_left_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%yvel1_right_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%yvel1_right_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%yvel1_bottom_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%yvel1_bottom_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%yvel1_top_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%yvel1_top_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%yvel1_left_top_snd_buffer(4))
      ALLOCATE(chunks(chunk)%yvel1_left_top_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%yvel1_right_top_snd_buffer(4))
      ALLOCATE(chunks(chunk)%yvel1_right_top_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%yvel1_right_bottom_snd_buffer(4))
      ALLOCATE(chunks(chunk)%yvel1_right_bottom_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%yvel1_left_bottom_snd_buffer(4))
      ALLOCATE(chunks(chunk)%yvel1_left_bottom_rcv_buffer(4))

      ALLOCATE(chunks(chunk)%volflux_x_left_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%volflux_x_left_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%volflux_x_right_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%volflux_x_right_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%volflux_x_bottom_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%volflux_x_bottom_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%volflux_x_top_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%volflux_x_top_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%volflux_x_left_top_snd_buffer(4))
      ALLOCATE(chunks(chunk)%volflux_x_left_top_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%volflux_x_right_top_snd_buffer(4))
      ALLOCATE(chunks(chunk)%volflux_x_right_top_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%volflux_x_right_bottom_snd_buffer(4))
      ALLOCATE(chunks(chunk)%volflux_x_right_bottom_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%volflux_x_left_bottom_snd_buffer(4))
      ALLOCATE(chunks(chunk)%volflux_x_left_bottom_rcv_buffer(4))

      ALLOCATE(chunks(chunk)%volflux_y_left_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%volflux_y_left_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%volflux_y_right_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%volflux_y_right_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%volflux_y_bottom_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%volflux_y_bottom_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%volflux_y_top_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%volflux_y_top_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%volflux_y_left_top_snd_buffer(4))
      ALLOCATE(chunks(chunk)%volflux_y_left_top_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%volflux_y_right_top_snd_buffer(4))
      ALLOCATE(chunks(chunk)%volflux_y_right_top_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%volflux_y_right_bottom_snd_buffer(4))
      ALLOCATE(chunks(chunk)%volflux_y_right_bottom_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%volflux_y_left_bottom_snd_buffer(4))
      ALLOCATE(chunks(chunk)%volflux_y_left_bottom_rcv_buffer(4))

      ALLOCATE(chunks(chunk)%massflux_x_left_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%massflux_x_left_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%massflux_x_right_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%massflux_x_right_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%massflux_x_bottom_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%massflux_x_bottom_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%massflux_x_top_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%massflux_x_top_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%massflux_x_left_top_snd_buffer(4))
      ALLOCATE(chunks(chunk)%massflux_x_left_top_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%massflux_x_right_top_snd_buffer(4))
      ALLOCATE(chunks(chunk)%massflux_x_right_top_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%massflux_x_right_bottom_snd_buffer(4))
      ALLOCATE(chunks(chunk)%massflux_x_right_bottom_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%massflux_x_left_bottom_snd_buffer(4))
      ALLOCATE(chunks(chunk)%massflux_x_left_bottom_rcv_buffer(4))

      ALLOCATE(chunks(chunk)%massflux_y_left_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%massflux_y_left_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%massflux_y_right_snd_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%massflux_y_right_rcv_buffer(2*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%massflux_y_bottom_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%massflux_y_bottom_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%massflux_y_top_snd_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%massflux_y_top_rcv_buffer(2*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%massflux_y_left_top_snd_buffer(4))
      ALLOCATE(chunks(chunk)%massflux_y_left_top_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%massflux_y_right_top_snd_buffer(4))
      ALLOCATE(chunks(chunk)%massflux_y_right_top_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%massflux_y_right_bottom_snd_buffer(4))
      ALLOCATE(chunks(chunk)%massflux_y_right_bottom_rcv_buffer(4))
      ALLOCATE(chunks(chunk)%massflux_y_left_bottom_snd_buffer(4))
      ALLOCATE(chunks(chunk)%massflux_y_left_bottom_rcv_buffer(4))
  ENDIF

END SUBROUTINE clover_allocate_buffers


SUBROUTINE clover_exchange(fields,depth)

    IMPLICIT NONE

    INTEGER      :: fields(:),depth

    CALL clover_exchange_send_async(parallel%task+1, depth, fields)

    CALL clover_exchange_receive_async(parallel%task+1, depth, fields)

END SUBROUTINE clover_exchange


SUBROUTINE clover_exchange_send_async(chunk, depth, fields)

    IMPLICIT NONE

    INTEGER :: chunk, depth, fields(NUM_FIELDS), receiver 

    IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
        !call method to write all buffers in this direction 
        CALL clover_exchange_write_all_buffers_left(chunk, depth, fields)
    ENDIF
    IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
        !call method to write all buffers in this direction 
        CALL clover_exchange_write_all_buffers_right(chunk, depth, fields)
    ENDIF
    IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
        !call method to write all buffers in this direction 
        CALL clover_exchange_write_all_buffers_bottom(chunk, depth, fields)
    ENDIF
    IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
        !call method to write all buffers in this direction 
        CALL clover_exchange_write_all_buffers_top(chunk, depth, fields)
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
        !call method to write all buffers in this direction 
        CALL clover_exchange_write_all_buffers_left_top(chunk, depth, fields)
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
        !call method to write all buffers in this direction 
        CALL clover_exchange_write_all_buffers_right_top(chunk, depth, fields)
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
        !call method to write all buffers in this direction 
        CALL clover_exchange_write_all_buffers_right_bottom(chunk, depth, fields)
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
        !call method to write all buffers in this direction 
        CALL clover_exchange_write_all_buffers_left_bottom(chunk, depth, fields)
    ENDIF

    ! Wait for the messages
#ifdef LOCAL_SYNC
    sync images( chunks(chunk)%imageNeighbours )
#else
    sync all
#endif

END SUBROUTINE clover_exchange_send_async


SUBROUTINE clover_exchange_receive_async(chunk, depth, fields)

    IMPLICIT NONE

    INTEGER :: chunk, depth, fields(NUM_FIELDS), receiver 

    IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
        !unpack the buffer
        CALL clover_exchange_unpack_all_buffers_left(chunk, depth, fields)
    ENDIF
    IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
        !unpack the buffer
        CALL clover_exchange_unpack_all_buffers_right(chunk, depth, fields)
    ENDIF
    IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
        !unpack the buffer
        CALL clover_exchange_unpack_all_buffers_bottom(chunk, depth, fields)
    ENDIF
    IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
        !unpack the buffer
        CALL clover_exchange_unpack_all_buffers_top(chunk, depth, fields)
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
        !unpack the buffer
        CALL clover_exchange_unpack_all_buffers_left_top(chunk, depth, fields)
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
        !unpack the buffer
        CALL clover_exchange_unpack_all_buffers_right_top(chunk, depth, fields)
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
        !unpack the buffer
        CALL clover_exchange_unpack_all_buffers_right_bottom(chunk, depth, fields)
    ENDIF
    IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
        !unpack the buffer
        CALL clover_exchange_unpack_all_buffers_left_bottom(chunk, depth, fields)
    ENDIF

    ! Wait for the messages
#ifdef LOCAL_SYNC
    sync images( chunks(chunk)%imageNeighbours )
#else
    sync all
#endif

END SUBROUTINE clover_exchange_receive_async






        
SUBROUTINE clover_exchange_write_all_buffers_left(chunk, depth, fields)

    IMPLICIT NONE 

    INTEGER :: chunk, depth, receiver, topedge, bottomedge, fields(NUM_FIELDS), left_neighbour_chunk, x_inc, y_inc, size
    
    left_neighbour_chunk = chunks(chunk)%chunk_neighbours(chunk_left)
    receiver=chunks(left_neighbour_chunk)%task + 1

    topedge = 0
    bottomedge = 0
    IF (chunks(chunk)%chunk_neighbours(chunk_top).EQ.external_face) THEN
        topedge = depth
    ENDIF
    IF (chunks(chunk)%chunk_neighbours(chunk_bottom).EQ.external_face) THEN
        bottomedge = depth
    ENDIF

    x_inc=0
    y_inc=0

    size=(1+(chunks(chunk)%field%y_max+y_inc+topedge)-(chunks(chunk)%field%y_min-bottomedge))*depth

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
        !CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
        !                                        chunks(chunk)%density0_left_snd_buffer, chunks(chunk)%density0_right_rcv_buffer, chunks(chunk)%field%density0)

        CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%density0_left_snd_buffer, chunks(chunk)%field%density0)

        chunks(left_neighbour_chunk)[receiver]%density0_right_rcv_buffer(1:size) = chunks(chunk)%density0_left_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        !CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
        !                                        chunks(chunk)%density1_left_snd_buffer, chunks(chunk)%density1_right_rcv_buffer, chunks(chunk)%field%density1)

        CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%density1_left_snd_buffer, chunks(chunk)%field%density1)

        chunks(left_neighbour_chunk)[receiver]%density1_right_rcv_buffer(1:size) = chunks(chunk)%density1_left_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        !CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
        !                                        chunks(chunk)%energy0_left_snd_buffer, chunks(chunk)%energy0_right_rcv_buffer, chunks(chunk)%field%energy0)
        CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%energy0_left_snd_buffer, chunks(chunk)%field%energy0)

        chunks(left_neighbour_chunk)[receiver]%energy0_right_rcv_buffer(1:size) = chunks(chunk)%energy0_left_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        !CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
        !                                        chunks(chunk)%energy1_left_snd_buffer, chunks(chunk)%energy1_right_rcv_buffer, chunks(chunk)%field%energy1)
        CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%energy1_left_snd_buffer, chunks(chunk)%field%energy1)

        chunks(left_neighbour_chunk)[receiver]%energy1_right_rcv_buffer(1:size) = chunks(chunk)%energy1_left_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        !CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
        !                                        chunks(chunk)%pressure_left_snd_buffer, chunks(chunk)%pressure_right_rcv_buffer, chunks(chunk)%field%pressure)
        CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%pressure_left_snd_buffer, chunks(chunk)%field%pressure)

        chunks(left_neighbour_chunk)[receiver]%pressure_right_rcv_buffer(1:size) = chunks(chunk)%pressure_left_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        !CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
        !                                        chunks(chunk)%viscosity_left_snd_buffer, chunks(chunk)%viscosity_right_rcv_buffer, chunks(chunk)%field%viscosity)
        CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%viscosity_left_snd_buffer, chunks(chunk)%field%viscosity)

        chunks(left_neighbour_chunk)[receiver]%viscosity_right_rcv_buffer(1:size) = chunks(chunk)%viscosity_left_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        !CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
        !                                        chunks(chunk)%soundspeed_left_snd_buffer, chunks(chunk)%soundspeed_right_rcv_buffer, chunks(chunk)%field%soundspeed)
        CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%soundspeed_left_snd_buffer, chunks(chunk)%field%soundspeed)

        chunks(left_neighbour_chunk)[receiver]%soundspeed_right_rcv_buffer(1:size) = chunks(chunk)%soundspeed_left_snd_buffer(1:size)
    ENDIF

    x_inc=1
    y_inc=1
    size=(1+(chunks(chunk)%field%y_max+y_inc+topedge)-(chunks(chunk)%field%y_min-bottomedge))*depth

    IF(fields(FIELD_XVEL0).EQ.1) THEN
        !CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, VERTEX_DATA, & 
        !                                        chunks(chunk)%xvel0_left_snd_buffer, chunks(chunk)%xvel0_right_rcv_buffer, chunks(chunk)%field%xvel0)
        CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%xvel0_left_snd_buffer, chunks(chunk)%field%xvel0)

        chunks(left_neighbour_chunk)[receiver]%xvel0_right_rcv_buffer(1:size) = chunks(chunk)%xvel0_left_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        !CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, VERTEX_DATA, & 
        !                                        chunks(chunk)%xvel1_left_snd_buffer, chunks(chunk)%xvel1_right_rcv_buffer, chunks(chunk)%field%xvel1)
        CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%xvel1_left_snd_buffer, chunks(chunk)%field%xvel1)

        chunks(left_neighbour_chunk)[receiver]%xvel1_right_rcv_buffer(1:size) = chunks(chunk)%xvel1_left_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        !CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, VERTEX_DATA, & 
        !                                        chunks(chunk)%yvel0_left_snd_buffer, chunks(chunk)%yvel0_right_rcv_buffer, chunks(chunk)%field%yvel0)
        CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%yvel0_left_snd_buffer, chunks(chunk)%field%yvel0)

        chunks(left_neighbour_chunk)[receiver]%yvel0_right_rcv_buffer(1:size) = chunks(chunk)%yvel0_left_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        !CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, VERTEX_DATA, & 
        !                                        chunks(chunk)%yvel1_left_snd_buffer, chunks(chunk)%yvel1_right_rcv_buffer, chunks(chunk)%field%yvel1)
        CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%yvel1_left_snd_buffer, chunks(chunk)%field%yvel1)

        chunks(left_neighbour_chunk)[receiver]%yvel1_right_rcv_buffer(1:size) = chunks(chunk)%yvel1_left_snd_buffer(1:size)
    ENDIF

    x_inc=1
    y_inc=0
    size=(1+(chunks(chunk)%field%y_max+y_inc+topedge)-(chunks(chunk)%field%y_min-bottomedge))*depth

    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        !CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, X_FACE_DATA, & 
        !                                        chunks(chunk)%volflux_x_left_snd_buffer, chunks(chunk)%volflux_x_right_rcv_buffer, chunks(chunk)%field%vol_flux_x)
        CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%volflux_x_left_snd_buffer, chunks(chunk)%field%vol_flux_x)

        chunks(left_neighbour_chunk)[receiver]%volflux_x_right_rcv_buffer(1:size) = chunks(chunk)%volflux_x_left_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        !CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, X_FACE_DATA, & 
        !                                        chunks(chunk)%massflux_x_left_snd_buffer, chunks(chunk)%massflux_x_right_rcv_buffer, chunks(chunk)%field%mass_flux_x)
        CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%massflux_x_left_snd_buffer, chunks(chunk)%field%mass_flux_x)

        chunks(left_neighbour_chunk)[receiver]%massflux_x_right_rcv_buffer(1:size) = chunks(chunk)%massflux_x_left_snd_buffer(1:size)
    ENDIF

    x_inc=0
    y_inc=1
    size=(1+(chunks(chunk)%field%y_max+y_inc+topedge)-(chunks(chunk)%field%y_min-bottomedge))*depth

    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        !CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, Y_FACE_DATA, & 
        !                                        chunks(chunk)%volflux_y_left_snd_buffer, chunks(chunk)%volflux_y_right_rcv_buffer, chunks(chunk)%field%vol_flux_y)
        CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%volflux_y_left_snd_buffer, chunks(chunk)%field%vol_flux_y)

        chunks(left_neighbour_chunk)[receiver]%volflux_y_right_rcv_buffer(1:size) = chunks(chunk)%volflux_y_left_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        !CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, Y_FACE_DATA, & 
        !                                        chunks(chunk)%massflux_y_left_snd_buffer, chunks(chunk)%massflux_y_right_rcv_buffer, chunks(chunk)%field%mass_flux_y)
        CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%massflux_y_left_snd_buffer, chunks(chunk)%field%mass_flux_y)

        chunks(left_neighbour_chunk)[receiver]%massflux_y_right_rcv_buffer(1:size) = chunks(chunk)%massflux_y_left_snd_buffer(1:size)
    ENDIF

END SUBROUTINE clover_exchange_write_all_buffers_left

SUBROUTINE clover_exchange_write_all_buffers_right(chunk, depth, fields)

    IMPLICIT NONE 

    INTEGER :: chunk, depth, receiver, topedge, bottomedge, fields(NUM_FIELDS), right_neighbour_chunk, x_inc, y_inc, size

    right_neighbour_chunk = chunks(chunk)%chunk_neighbours(chunk_right)
    receiver=chunks(right_neighbour_chunk)%task + 1

    topedge = 0
    bottomedge = 0
    IF (chunks(chunk)%chunk_neighbours(chunk_top).EQ.external_face) THEN
        topedge = depth
    ENDIF
    IF (chunks(chunk)%chunk_neighbours(chunk_bottom).EQ.external_face) THEN
        bottomedge = depth
    ENDIF

    x_inc=0
    y_inc=0
    size=(1+(chunks(chunk)%field%y_max+y_inc+topedge)-(chunks(chunk)%field%y_min-bottomedge))*depth

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
        !CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
        !                                         chunks(chunk)%density0_right_snd_buffer, chunks(chunk)%density0_left_rcv_buffer, chunks(chunk)%field%density0)

        CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%density0_right_snd_buffer, chunks(chunk)%field%density0)

        chunks(right_neighbour_chunk)[receiver]%density0_left_rcv_buffer(1:size) = chunks(chunk)%density0_right_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        !CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
        !                                         chunks(chunk)%density1_right_snd_buffer, chunks(chunk)%density1_left_rcv_buffer, chunks(chunk)%field%density1)
        CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%density1_right_snd_buffer, chunks(chunk)%field%density1)

        chunks(right_neighbour_chunk)[receiver]%density1_left_rcv_buffer(1:size) = chunks(chunk)%density1_right_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        !CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
        !                                         chunks(chunk)%energy0_right_snd_buffer, chunks(chunk)%energy0_left_rcv_buffer, chunks(chunk)%field%energy0)
        CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%energy0_right_snd_buffer, chunks(chunk)%field%energy0)

        chunks(right_neighbour_chunk)[receiver]%energy0_left_rcv_buffer(1:size) = chunks(chunk)%energy0_right_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        !CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
        !                                         chunks(chunk)%energy1_right_snd_buffer, chunks(chunk)%energy1_left_rcv_buffer, chunks(chunk)%field%energy1)
        CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%energy1_right_snd_buffer, chunks(chunk)%field%energy1)

        chunks(right_neighbour_chunk)[receiver]%energy1_left_rcv_buffer(1:size) = chunks(chunk)%energy1_right_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        !CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
        !                                         chunks(chunk)%pressure_right_snd_buffer, chunks(chunk)%pressure_left_rcv_buffer, chunks(chunk)%field%pressure)
        CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%pressure_right_snd_buffer, chunks(chunk)%field%pressure)

        chunks(right_neighbour_chunk)[receiver]%pressure_left_rcv_buffer(1:size) = chunks(chunk)%pressure_right_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        !CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
        !                                         chunks(chunk)%viscosity_right_snd_buffer, chunks(chunk)%viscosity_left_rcv_buffer, chunks(chunk)%field%viscosity)
        CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%viscosity_right_snd_buffer, chunks(chunk)%field%viscosity)

        chunks(right_neighbour_chunk)[receiver]%viscosity_left_rcv_buffer(1:size) = chunks(chunk)%viscosity_right_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        !CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
        !                                         chunks(chunk)%soundspeed_right_snd_buffer, chunks(chunk)%soundspeed_left_rcv_buffer, chunks(chunk)%field%soundspeed)
        CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%soundspeed_right_snd_buffer, chunks(chunk)%field%soundspeed)

        chunks(right_neighbour_chunk)[receiver]%soundspeed_left_rcv_buffer(1:size) = chunks(chunk)%soundspeed_right_snd_buffer(1:size)
    ENDIF

    x_inc=1
    y_inc=1
    size=(1+(chunks(chunk)%field%y_max+y_inc+topedge)-(chunks(chunk)%field%y_min-bottomedge))*depth

    IF(fields(FIELD_XVEL0).EQ.1) THEN
        !CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, VERTEX_DATA, & 
        !                                         chunks(chunk)%xvel0_right_snd_buffer, chunks(chunk)%xvel0_left_rcv_buffer, chunks(chunk)%field%xvel0)
        CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%xvel0_right_snd_buffer, chunks(chunk)%field%xvel0)

        chunks(right_neighbour_chunk)[receiver]%xvel0_left_rcv_buffer(1:size) = chunks(chunk)%xvel0_right_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        !CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, VERTEX_DATA, & 
        !                                         chunks(chunk)%xvel1_right_snd_buffer, chunks(chunk)%xvel1_left_rcv_buffer, chunks(chunk)%field%xvel1)
        CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%xvel1_right_snd_buffer, chunks(chunk)%field%xvel1)

        chunks(right_neighbour_chunk)[receiver]%xvel1_left_rcv_buffer(1:size) = chunks(chunk)%xvel1_right_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        !CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, VERTEX_DATA, & 
        !                                         chunks(chunk)%yvel0_right_snd_buffer, chunks(chunk)%yvel0_left_rcv_buffer, chunks(chunk)%field%yvel0)
        CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%yvel0_right_snd_buffer, chunks(chunk)%field%yvel0)

        chunks(right_neighbour_chunk)[receiver]%yvel0_left_rcv_buffer(1:size) = chunks(chunk)%yvel0_right_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        !CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, VERTEX_DATA, & 
        !                                         chunks(chunk)%yvel1_right_snd_buffer, chunks(chunk)%yvel1_left_rcv_buffer, chunks(chunk)%field%yvel1)
        CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%yvel1_right_snd_buffer, chunks(chunk)%field%yvel1)

        chunks(right_neighbour_chunk)[receiver]%yvel1_left_rcv_buffer(1:size) = chunks(chunk)%yvel1_right_snd_buffer(1:size)
    ENDIF

    x_inc=1
    y_inc=0
    size=(1+(chunks(chunk)%field%y_max+y_inc+topedge)-(chunks(chunk)%field%y_min-bottomedge))*depth

    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        !CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, X_FACE_DATA, & 
        !                                         chunks(chunk)%volflux_x_right_snd_buffer, chunks(chunk)%volflux_x_left_rcv_buffer, chunks(chunk)%field%vol_flux_x)
        CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%volflux_x_right_snd_buffer, chunks(chunk)%field%vol_flux_x)

        chunks(right_neighbour_chunk)[receiver]%volflux_x_left_rcv_buffer(1:size) = chunks(chunk)%volflux_x_right_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        !CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, X_FACE_DATA, & 
        !                                         chunks(chunk)%massflux_x_right_snd_buffer, chunks(chunk)%massflux_x_left_rcv_buffer, chunks(chunk)%field%mass_flux_x)
        CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%massflux_x_right_snd_buffer, chunks(chunk)%field%mass_flux_x)

        chunks(right_neighbour_chunk)[receiver]%massflux_x_left_rcv_buffer(1:size) = chunks(chunk)%massflux_x_right_snd_buffer(1:size)
    ENDIF

    x_inc=0
    y_inc=1
    size=(1+(chunks(chunk)%field%y_max+y_inc+topedge)-(chunks(chunk)%field%y_min-bottomedge))*depth
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        !CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, Y_FACE_DATA, & 
        !                                         chunks(chunk)%volflux_y_right_snd_buffer, chunks(chunk)%volflux_y_left_rcv_buffer, chunks(chunk)%field%vol_flux_y)
        CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%volflux_y_right_snd_buffer, chunks(chunk)%field%vol_flux_y)

        chunks(right_neighbour_chunk)[receiver]%volflux_y_left_rcv_buffer(1:size) = chunks(chunk)%volflux_y_right_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        !CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, Y_FACE_DATA, & 
        !                                         chunks(chunk)%massflux_y_right_snd_buffer, chunks(chunk)%massflux_y_left_rcv_buffer, chunks(chunk)%field%mass_flux_y)
        CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%massflux_y_right_snd_buffer, chunks(chunk)%field%mass_flux_y)

        chunks(right_neighbour_chunk)[receiver]%massflux_y_left_rcv_buffer(1:size) = chunks(chunk)%massflux_y_right_snd_buffer(1:size)
    ENDIF

END SUBROUTINE clover_exchange_write_all_buffers_right

SUBROUTINE clover_exchange_write_all_buffers_bottom(chunk, depth, fields)

    IMPLICIT NONE 

    INTEGER :: chunk, depth, receiver, topedge, bottomedge, fields(NUM_FIELDS), leftedge, rightedge, bottom_neighbour_chunk, x_inc, y_inc, size

    bottom_neighbour_chunk = chunks(chunk)%chunk_neighbours(chunk_bottom)
    receiver=chunks(bottom_neighbour_chunk)%task + 1

    leftedge= 0
    rightedge= 0
    IF (chunks(chunk)%chunk_neighbours(chunk_left).EQ.external_face) THEN
        leftedge = depth
    ENDIF
    IF (chunks(chunk)%chunk_neighbours(chunk_right).EQ.external_face) THEN
        rightedge = depth
    ENDIF

    x_inc=0
    y_inc=0
    size=(1+(chunks(chunk)%field%x_max+x_inc+rightedge)-(chunks(chunk)%field%x_min-leftedge))*depth

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
        !CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
        !                                          chunks(chunk)%density0_bottom_snd_buffer, chunks(chunk)%density0_top_rcv_buffer, chunks(chunk)%field%density0)
        CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%density0_bottom_snd_buffer, chunks(chunk)%field%density0)

        chunks(bottom_neighbour_chunk)[receiver]%density0_top_rcv_buffer(1:size) = chunks(chunk)%density0_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        !CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
        !                                          chunks(chunk)%density1_bottom_snd_buffer, chunks(chunk)%density1_top_rcv_buffer, chunks(chunk)%field%density1)
        CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%density1_bottom_snd_buffer, chunks(chunk)%field%density1)

        chunks(bottom_neighbour_chunk)[receiver]%density1_top_rcv_buffer(1:size) = chunks(chunk)%density1_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        !CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
        !                                          chunks(chunk)%energy0_bottom_snd_buffer, chunks(chunk)%energy0_top_rcv_buffer, chunks(chunk)%field%energy0)
        CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%energy0_bottom_snd_buffer, chunks(chunk)%field%energy0)

        chunks(bottom_neighbour_chunk)[receiver]%energy0_top_rcv_buffer(1:size) = chunks(chunk)%energy0_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        !CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
        !                                          chunks(chunk)%energy1_bottom_snd_buffer, chunks(chunk)%energy1_top_rcv_buffer, chunks(chunk)%field%energy1)
        CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%energy1_bottom_snd_buffer, chunks(chunk)%field%energy1)

        chunks(bottom_neighbour_chunk)[receiver]%energy1_top_rcv_buffer(1:size) = chunks(chunk)%energy1_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        !CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
        !                                          chunks(chunk)%pressure_bottom_snd_buffer, chunks(chunk)%pressure_top_rcv_buffer, chunks(chunk)%field%pressure)
        CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%pressure_bottom_snd_buffer, chunks(chunk)%field%pressure)

        chunks(bottom_neighbour_chunk)[receiver]%pressure_top_rcv_buffer(1:size) = chunks(chunk)%pressure_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        !CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
        !                                          chunks(chunk)%viscosity_bottom_snd_buffer, chunks(chunk)%viscosity_top_rcv_buffer, chunks(chunk)%field%viscosity)
        CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%viscosity_bottom_snd_buffer, chunks(chunk)%field%viscosity)

        chunks(bottom_neighbour_chunk)[receiver]%viscosity_top_rcv_buffer(1:size) = chunks(chunk)%viscosity_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        !CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
        !                                          chunks(chunk)%soundspeed_bottom_snd_buffer, chunks(chunk)%soundspeed_top_rcv_buffer, chunks(chunk)%field%soundspeed)
        CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%soundspeed_bottom_snd_buffer, chunks(chunk)%field%soundspeed)

        chunks(bottom_neighbour_chunk)[receiver]%soundspeed_top_rcv_buffer(1:size) = chunks(chunk)%soundspeed_bottom_snd_buffer(1:size)
    ENDIF

    x_inc=1
    y_inc=1
    size=(1+(chunks(chunk)%field%x_max+x_inc+rightedge)-(chunks(chunk)%field%x_min-leftedge))*depth

    IF(fields(FIELD_XVEL0).EQ.1) THEN
        !CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, VERTEX_DATA, & 
        !                                          chunks(chunk)%xvel0_bottom_snd_buffer, chunks(chunk)%xvel0_top_rcv_buffer, chunks(chunk)%field%xvel0)
        CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%xvel0_bottom_snd_buffer, chunks(chunk)%field%xvel0)

        chunks(bottom_neighbour_chunk)[receiver]%xvel0_top_rcv_buffer(1:size) = chunks(chunk)%xvel0_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        !CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, VERTEX_DATA, & 
        !                                          chunks(chunk)%xvel1_bottom_snd_buffer, chunks(chunk)%xvel1_top_rcv_buffer, chunks(chunk)%field%xvel1)
        CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%xvel1_bottom_snd_buffer, chunks(chunk)%field%xvel1)

        chunks(bottom_neighbour_chunk)[receiver]%xvel1_top_rcv_buffer(1:size) = chunks(chunk)%xvel1_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        !CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, VERTEX_DATA, & 
        !                                          chunks(chunk)%yvel0_bottom_snd_buffer, chunks(chunk)%yvel0_top_rcv_buffer, chunks(chunk)%field%yvel0)
        CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%yvel0_bottom_snd_buffer, chunks(chunk)%field%yvel0)

        chunks(bottom_neighbour_chunk)[receiver]%yvel0_top_rcv_buffer(1:size) = chunks(chunk)%yvel0_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        !CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, VERTEX_DATA, & 
        !                                          chunks(chunk)%yvel1_bottom_snd_buffer, chunks(chunk)%yvel1_top_rcv_buffer, chunks(chunk)%field%yvel1)
        CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%yvel1_bottom_snd_buffer, chunks(chunk)%field%yvel1)

        chunks(bottom_neighbour_chunk)[receiver]%yvel1_top_rcv_buffer(1:size) = chunks(chunk)%yvel1_bottom_snd_buffer(1:size)
    ENDIF

    x_inc=1
    y_inc=0
    size=(1+(chunks(chunk)%field%x_max+x_inc+rightedge)-(chunks(chunk)%field%x_min-leftedge))*depth

    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        !CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, X_FACE_DATA, & 
        !                                          chunks(chunk)%volflux_x_bottom_snd_buffer, chunks(chunk)%volflux_x_top_rcv_buffer, chunks(chunk)%field%vol_flux_x)
        CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%volflux_x_bottom_snd_buffer, chunks(chunk)%field%vol_flux_x)

        chunks(bottom_neighbour_chunk)[receiver]%volflux_x_top_rcv_buffer(1:size) = chunks(chunk)%volflux_x_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        !CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, X_FACE_DATA, & 
        !                                          chunks(chunk)%massflux_x_bottom_snd_buffer, chunks(chunk)%massflux_x_top_rcv_buffer, chunks(chunk)%field%mass_flux_x)
        CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%massflux_x_bottom_snd_buffer, chunks(chunk)%field%mass_flux_x)

        chunks(bottom_neighbour_chunk)[receiver]%massflux_x_top_rcv_buffer(1:size) = chunks(chunk)%massflux_x_bottom_snd_buffer(1:size)
    ENDIF

    x_inc=0
    y_inc=1
    size=(1+(chunks(chunk)%field%x_max+x_inc+rightedge)-(chunks(chunk)%field%x_min-leftedge))*depth

    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        !CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, Y_FACE_DATA, & 
        !                                          chunks(chunk)%volflux_y_bottom_snd_buffer, chunks(chunk)%volflux_y_top_rcv_buffer, chunks(chunk)%field%vol_flux_y)
        CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%volflux_y_bottom_snd_buffer, chunks(chunk)%field%vol_flux_y)

        chunks(bottom_neighbour_chunk)[receiver]%volflux_y_top_rcv_buffer(1:size) = chunks(chunk)%volflux_y_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        !CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, Y_FACE_DATA, & 
        !                                          chunks(chunk)%massflux_y_bottom_snd_buffer, chunks(chunk)%massflux_y_top_rcv_buffer, chunks(chunk)%field%mass_flux_y)
        CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%massflux_y_bottom_snd_buffer, chunks(chunk)%field%mass_flux_y)

        chunks(bottom_neighbour_chunk)[receiver]%massflux_y_top_rcv_buffer(1:size) = chunks(chunk)%massflux_y_bottom_snd_buffer(1:size)
    ENDIF

END SUBROUTINE clover_exchange_write_all_buffers_bottom


SUBROUTINE clover_exchange_write_all_buffers_top(chunk, depth, fields)

    IMPLICIT NONE 

    INTEGER :: chunk, depth, receiver, topedge, bottomedge, fields(NUM_FIELDS), leftedge, rightedge, top_neighbour_chunk, x_inc, y_inc, size

    top_neighbour_chunk = chunks(chunk)%chunk_neighbours(chunk_top)
    receiver=chunks(top_neighbour_chunk)%task + 1

    leftedge= 0
    rightedge= 0
    IF (chunks(chunk)%chunk_neighbours(chunk_left).EQ.external_face) THEN
        leftedge = depth
    ENDIF
    IF (chunks(chunk)%chunk_neighbours(chunk_right).EQ.external_face) THEN
        rightedge = depth
    ENDIF

    x_inc=0
    y_inc=0
    size=(1+(chunks(chunk)%field%x_max+x_inc+rightedge)-(chunks(chunk)%field%x_min-leftedge))*depth

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
        !CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
        !                                       chunks(chunk)%density0_top_snd_buffer, chunks(chunk)%density0_bottom_rcv_buffer, chunks(chunk)%field%density0)
        CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%density0_top_snd_buffer, chunks(chunk)%field%density0)

        chunks(top_neighbour_chunk)[receiver]%density0_bottom_rcv_buffer(1:size) = chunks(chunk)%density0_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        !CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
        !                                       chunks(chunk)%density1_top_snd_buffer, chunks(chunk)%density1_bottom_rcv_buffer, chunks(chunk)%field%density1)
        CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%density1_top_snd_buffer, chunks(chunk)%field%density1)

        chunks(top_neighbour_chunk)[receiver]%density1_bottom_rcv_buffer(1:size) = chunks(chunk)%density1_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        !CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
        !                                       chunks(chunk)%energy0_top_snd_buffer, chunks(chunk)%energy0_bottom_rcv_buffer, chunks(chunk)%field%energy0)
        CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%energy0_top_snd_buffer, chunks(chunk)%field%energy0)

        chunks(top_neighbour_chunk)[receiver]%energy0_bottom_rcv_buffer(1:size) = chunks(chunk)%energy0_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        !CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
        !                                       chunks(chunk)%energy1_top_snd_buffer, chunks(chunk)%energy1_bottom_rcv_buffer, chunks(chunk)%field%energy1)
        CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%energy1_top_snd_buffer, chunks(chunk)%field%energy1)

        chunks(top_neighbour_chunk)[receiver]%energy1_bottom_rcv_buffer(1:size) = chunks(chunk)%energy1_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        !CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
        !                                       chunks(chunk)%pressure_top_snd_buffer, chunks(chunk)%pressure_bottom_rcv_buffer, chunks(chunk)%field%pressure)
        CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%pressure_top_snd_buffer, chunks(chunk)%field%pressure)

        chunks(top_neighbour_chunk)[receiver]%pressure_bottom_rcv_buffer(1:size) = chunks(chunk)%pressure_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        !CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
        !                                       chunks(chunk)%viscosity_top_snd_buffer, chunks(chunk)%viscosity_bottom_rcv_buffer, chunks(chunk)%field%viscosity)
        CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%viscosity_top_snd_buffer, chunks(chunk)%field%viscosity)

        chunks(top_neighbour_chunk)[receiver]%viscosity_bottom_rcv_buffer(1:size) = chunks(chunk)%viscosity_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        !CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
        !                                       chunks(chunk)%soundspeed_top_snd_buffer, chunks(chunk)%soundspeed_bottom_rcv_buffer, chunks(chunk)%field%soundspeed)
        CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%soundspeed_top_snd_buffer, chunks(chunk)%field%soundspeed)

        chunks(top_neighbour_chunk)[receiver]%soundspeed_bottom_rcv_buffer(1:size) = chunks(chunk)%soundspeed_top_snd_buffer(1:size)
    ENDIF

    x_inc=1
    y_inc=1
    size=(1+(chunks(chunk)%field%x_max+x_inc+rightedge)-(chunks(chunk)%field%x_min-leftedge))*depth

    IF(fields(FIELD_XVEL0).EQ.1) THEN
        !CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, VERTEX_DATA, & 
        !                                       chunks(chunk)%xvel0_top_snd_buffer, chunks(chunk)%xvel0_bottom_rcv_buffer, chunks(chunk)%field%xvel0)
        CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%xvel0_top_snd_buffer, chunks(chunk)%field%xvel0)

        chunks(top_neighbour_chunk)[receiver]%xvel0_bottom_rcv_buffer(1:size) = chunks(chunk)%xvel0_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        !CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, VERTEX_DATA, & 
        !                                       chunks(chunk)%xvel1_top_snd_buffer, chunks(chunk)%xvel1_bottom_rcv_buffer, chunks(chunk)%field%xvel1)
        CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%xvel1_top_snd_buffer, chunks(chunk)%field%xvel1)

        chunks(top_neighbour_chunk)[receiver]%xvel1_bottom_rcv_buffer(1:size) = chunks(chunk)%xvel1_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        !CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, VERTEX_DATA, & 
        !                                       chunks(chunk)%yvel0_top_snd_buffer, chunks(chunk)%yvel0_bottom_rcv_buffer, chunks(chunk)%field%yvel0)
        CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%yvel0_top_snd_buffer, chunks(chunk)%field%yvel0)

        chunks(top_neighbour_chunk)[receiver]%yvel0_bottom_rcv_buffer(1:size) = chunks(chunk)%yvel0_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        !CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, VERTEX_DATA, & 
        !                                       chunks(chunk)%yvel1_top_snd_buffer, chunks(chunk)%yvel1_bottom_rcv_buffer, chunks(chunk)%field%yvel1)
        CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%yvel1_top_snd_buffer, chunks(chunk)%field%yvel1)

        chunks(top_neighbour_chunk)[receiver]%yvel1_bottom_rcv_buffer(1:size) = chunks(chunk)%yvel1_top_snd_buffer(1:size)
    ENDIF

    x_inc=1
    y_inc=0
    size=(1+(chunks(chunk)%field%x_max+x_inc+rightedge)-(chunks(chunk)%field%x_min-leftedge))*depth
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        !CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, X_FACE_DATA, & 
        !                                       chunks(chunk)%volflux_x_top_snd_buffer, chunks(chunk)%volflux_x_bottom_rcv_buffer, chunks(chunk)%field%vol_flux_x)
        CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%volflux_x_top_snd_buffer, chunks(chunk)%field%vol_flux_x)

        chunks(top_neighbour_chunk)[receiver]%volflux_x_bottom_rcv_buffer(1:size) = chunks(chunk)%volflux_x_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        !CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, X_FACE_DATA, & 
        !                                       chunks(chunk)%massflux_x_top_snd_buffer, chunks(chunk)%massflux_x_bottom_rcv_buffer, chunks(chunk)%field%mass_flux_x)
        CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%massflux_x_top_snd_buffer, chunks(chunk)%field%mass_flux_x)

        chunks(top_neighbour_chunk)[receiver]%massflux_x_bottom_rcv_buffer(1:size) = chunks(chunk)%massflux_x_top_snd_buffer(1:size)
    ENDIF

    x_inc=0
    y_inc=1
    size=(1+(chunks(chunk)%field%x_max+x_inc+rightedge)-(chunks(chunk)%field%x_min-leftedge))*depth
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        !CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, Y_FACE_DATA, & 
        !                                       chunks(chunk)%volflux_y_top_snd_buffer, chunks(chunk)%volflux_y_bottom_rcv_buffer, chunks(chunk)%field%vol_flux_y)
        CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%volflux_y_top_snd_buffer, chunks(chunk)%field%vol_flux_y)

        chunks(top_neighbour_chunk)[receiver]%volflux_y_bottom_rcv_buffer(1:size) = chunks(chunk)%volflux_y_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        !CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, Y_FACE_DATA, & 
        !                                       chunks(chunk)%massflux_y_top_snd_buffer, chunks(chunk)%massflux_y_bottom_rcv_buffer, chunks(chunk)%field%mass_flux_y)
        CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%massflux_y_top_snd_buffer, chunks(chunk)%field%mass_flux_y)

        chunks(top_neighbour_chunk)[receiver]%massflux_y_bottom_rcv_buffer(1:size) = chunks(chunk)%massflux_y_top_snd_buffer(1:size)
    ENDIF

END SUBROUTINE clover_exchange_write_all_buffers_top

SUBROUTINE clover_exchange_write_all_buffers_left_top(chunk, depth, fields)

    IMPLICIT NONE 

    INTEGER :: size, chunk, depth, receiver, topedge, bottomedge, fields(NUM_FIELDS), left_top_neighbour_chunk, x_inc, y_inc

    left_top_neighbour_chunk = chunks(chunk)%chunk_neighbours(chunk_left_top)
    receiver = chunks(left_top_neighbour_chunk)%task + 1

    size = depth*depth

    x_inc=0
    y_inc=0

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
        !CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, CELL_DATA, & 
        !                                            chunks(chunk)%density0_left_top_snd_buffer, chunks(chunk)%density0_right_bottom_rcv_buffer, chunks(chunk)%field%density0)
        CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%density0_left_top_snd_buffer, chunks(chunk)%field%density0)

        chunks(left_top_neighbour_chunk)[receiver]%density0_right_bottom_rcv_buffer(1:size) = chunks(chunk)%density0_left_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        !CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, CELL_DATA, & 
        !                                            chunks(chunk)%density1_left_top_snd_buffer, chunks(chunk)%density1_right_bottom_rcv_buffer, chunks(chunk)%field%density1)
        CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%density1_left_top_snd_buffer, chunks(chunk)%field%density1)

        chunks(left_top_neighbour_chunk)[receiver]%density1_right_bottom_rcv_buffer(1:size) = chunks(chunk)%density1_left_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        !CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, CELL_DATA, & 
        !                                            chunks(chunk)%energy0_left_top_snd_buffer, chunks(chunk)%energy0_right_bottom_rcv_buffer, chunks(chunk)%field%energy0)
        CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%energy0_left_top_snd_buffer, chunks(chunk)%field%energy0)

        chunks(left_top_neighbour_chunk)[receiver]%energy0_right_bottom_rcv_buffer(1:size) = chunks(chunk)%energy0_left_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        !CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, CELL_DATA, & 
        !                                            chunks(chunk)%energy1_left_top_snd_buffer, chunks(chunk)%energy1_right_bottom_rcv_buffer, chunks(chunk)%field%energy1)
        CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%energy1_left_top_snd_buffer, chunks(chunk)%field%energy1)

        chunks(left_top_neighbour_chunk)[receiver]%energy1_right_bottom_rcv_buffer(1:size) = chunks(chunk)%energy1_left_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        !CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, CELL_DATA, & 
        !                                            chunks(chunk)%pressure_left_top_snd_buffer, chunks(chunk)%pressure_right_bottom_rcv_buffer, chunks(chunk)%field%pressure)
        CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%pressure_left_top_snd_buffer, chunks(chunk)%field%pressure)

        chunks(left_top_neighbour_chunk)[receiver]%pressure_right_bottom_rcv_buffer(1:size) = chunks(chunk)%pressure_left_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        !CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, CELL_DATA, & 
        !                                            chunks(chunk)%viscosity_left_top_snd_buffer, chunks(chunk)%viscosity_right_bottom_rcv_buffer, chunks(chunk)%field%viscosity)
        CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%viscosity_left_top_snd_buffer, chunks(chunk)%field%viscosity)

        chunks(left_top_neighbour_chunk)[receiver]%viscosity_right_bottom_rcv_buffer(1:size) = chunks(chunk)%viscosity_left_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        !CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, CELL_DATA, & 
        !                                            chunks(chunk)%soundspeed_left_top_snd_buffer, chunks(chunk)%soundspeed_right_bottom_rcv_buffer, chunks(chunk)%field%soundspeed)
        CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%soundspeed_left_top_snd_buffer, chunks(chunk)%field%soundspeed)

        chunks(left_top_neighbour_chunk)[receiver]%soundspeed_right_bottom_rcv_buffer(1:size) = chunks(chunk)%soundspeed_left_top_snd_buffer(1:size)
    ENDIF

    x_inc=1
    y_inc=1
    IF(fields(FIELD_XVEL0).EQ.1) THEN
        !CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, VERTEX_DATA, & 
        !                                            chunks(chunk)%xvel0_left_top_snd_buffer, chunks(chunk)%xvel0_right_bottom_rcv_buffer, chunks(chunk)%field%xvel0)
        CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%xvel0_left_top_snd_buffer, chunks(chunk)%field%xvel0)

        chunks(left_top_neighbour_chunk)[receiver]%xvel0_right_bottom_rcv_buffer(1:size) = chunks(chunk)%xvel0_left_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        !CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, VERTEX_DATA, & 
        !                                            chunks(chunk)%xvel1_left_top_snd_buffer, chunks(chunk)%xvel1_right_bottom_rcv_buffer, chunks(chunk)%field%xvel1)
        CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%xvel1_left_top_snd_buffer, chunks(chunk)%field%xvel1)

        chunks(left_top_neighbour_chunk)[receiver]%xvel1_right_bottom_rcv_buffer(1:size) = chunks(chunk)%xvel1_left_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        !CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, VERTEX_DATA, & 
        !                                            chunks(chunk)%yvel0_left_top_snd_buffer, chunks(chunk)%yvel0_right_bottom_rcv_buffer, chunks(chunk)%field%yvel0)
        CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%yvel0_left_top_snd_buffer, chunks(chunk)%field%yvel0)

        chunks(left_top_neighbour_chunk)[receiver]%yvel0_right_bottom_rcv_buffer(1:size) = chunks(chunk)%yvel0_left_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        !CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, VERTEX_DATA, & 
        !                                            chunks(chunk)%yvel1_left_top_snd_buffer, chunks(chunk)%yvel1_right_bottom_rcv_buffer, chunks(chunk)%field%yvel1)
        CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%yvel1_left_top_snd_buffer, chunks(chunk)%field%yvel1)

        chunks(left_top_neighbour_chunk)[receiver]%yvel1_right_bottom_rcv_buffer(1:size) = chunks(chunk)%yvel1_left_top_snd_buffer(1:size)
    ENDIF

    x_inc=1
    y_inc=0
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        !CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, X_FACE_DATA, & 
        !                                            chunks(chunk)%volflux_x_left_top_snd_buffer, chunks(chunk)%volflux_x_right_bottom_rcv_buffer, chunks(chunk)%field%vol_flux_x)
        CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%volflux_x_left_top_snd_buffer, chunks(chunk)%field%vol_flux_x)

        chunks(left_top_neighbour_chunk)[receiver]%volflux_x_right_bottom_rcv_buffer(1:size) = chunks(chunk)%volflux_x_left_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        !CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, X_FACE_DATA, & 
        !                                            chunks(chunk)%massflux_x_left_top_snd_buffer, chunks(chunk)%massflux_x_right_bottom_rcv_buffer, chunks(chunk)%field%mass_flux_x)
        CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%massflux_x_left_top_snd_buffer, chunks(chunk)%field%mass_flux_x)

        chunks(left_top_neighbour_chunk)[receiver]%massflux_x_right_bottom_rcv_buffer(1:size) = chunks(chunk)%massflux_x_left_top_snd_buffer(1:size)
    ENDIF

    x_inc=0
    y_inc=1
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        !CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, Y_FACE_DATA, & 
        !                                            chunks(chunk)%volflux_y_left_top_snd_buffer, chunks(chunk)%volflux_y_right_bottom_rcv_buffer, chunks(chunk)%field%vol_flux_y)
        CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%volflux_y_left_top_snd_buffer, chunks(chunk)%field%vol_flux_y)

        chunks(left_top_neighbour_chunk)[receiver]%volflux_y_right_bottom_rcv_buffer(1:size) = chunks(chunk)%volflux_y_left_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        !CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, Y_FACE_DATA, & 
        !                                            chunks(chunk)%massflux_y_left_top_snd_buffer, chunks(chunk)%massflux_y_right_bottom_rcv_buffer, chunks(chunk)%field%mass_flux_y)
        CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%massflux_y_left_top_snd_buffer, chunks(chunk)%field%mass_flux_y)

        chunks(left_top_neighbour_chunk)[receiver]%massflux_y_right_bottom_rcv_buffer(1:size) = chunks(chunk)%massflux_y_left_top_snd_buffer(1:size)
    ENDIF

END SUBROUTINE clover_exchange_write_all_buffers_left_top

SUBROUTINE clover_exchange_write_all_buffers_right_top(chunk, depth, fields)

    IMPLICIT NONE 

    INTEGER :: size, chunk, depth, receiver, topedge, bottomedge, fields(NUM_FIELDS), right_top_neighbour_chunk, x_inc, y_inc

    right_top_neighbour_chunk = chunks(chunk)%chunk_neighbours(chunk_right_top)
    receiver = chunks(right_top_neighbour_chunk)%task + 1

    size = depth*depth

    x_inc=0
    y_inc=0

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
        !CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, CELL_DATA, & 
        !                                             chunks(chunk)%density0_right_top_snd_buffer, chunks(chunk)%density0_left_bottom_rcv_buffer, chunks(chunk)%field%density0)
        CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%density0_right_top_snd_buffer, chunks(chunk)%field%density0)

        chunks(right_top_neighbour_chunk)[receiver]%density0_left_bottom_rcv_buffer(1:size) = chunks(chunk)%density0_right_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        !CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, CELL_DATA, & 
        !                                             chunks(chunk)%density1_right_top_snd_buffer, chunks(chunk)%density1_left_bottom_rcv_buffer, chunks(chunk)%field%density1)
        CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%density1_right_top_snd_buffer, chunks(chunk)%field%density1)

        chunks(right_top_neighbour_chunk)[receiver]%density1_left_bottom_rcv_buffer(1:size) = chunks(chunk)%density1_right_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        !CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, CELL_DATA, & 
        !                                             chunks(chunk)%energy0_right_top_snd_buffer, chunks(chunk)%energy0_left_bottom_rcv_buffer, chunks(chunk)%field%energy0)
        CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%energy0_right_top_snd_buffer, chunks(chunk)%field%energy0)

        chunks(right_top_neighbour_chunk)[receiver]%energy0_left_bottom_rcv_buffer(1:size) = chunks(chunk)%energy0_right_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        !CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, CELL_DATA, & 
        !                                             chunks(chunk)%energy1_right_top_snd_buffer, chunks(chunk)%energy1_left_bottom_rcv_buffer, chunks(chunk)%field%energy1)
        CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%energy1_right_top_snd_buffer, chunks(chunk)%field%energy1)

        chunks(right_top_neighbour_chunk)[receiver]%energy1_left_bottom_rcv_buffer(1:size) = chunks(chunk)%energy1_right_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        !CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, CELL_DATA, & 
        !                                             chunks(chunk)%pressure_right_top_snd_buffer, chunks(chunk)%pressure_left_bottom_rcv_buffer, chunks(chunk)%field%pressure)
        CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%pressure_right_top_snd_buffer, chunks(chunk)%field%pressure)

        chunks(right_top_neighbour_chunk)[receiver]%pressure_left_bottom_rcv_buffer(1:size) = chunks(chunk)%pressure_right_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        !CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, CELL_DATA, & 
        !                                             chunks(chunk)%viscosity_right_top_snd_buffer, chunks(chunk)%viscosity_left_bottom_rcv_buffer, chunks(chunk)%field%viscosity)
        CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%viscosity_right_top_snd_buffer, chunks(chunk)%field%viscosity)

        chunks(right_top_neighbour_chunk)[receiver]%viscosity_left_bottom_rcv_buffer(1:size) = chunks(chunk)%viscosity_right_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        !CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, CELL_DATA, & 
        !                                             chunks(chunk)%soundspeed_right_top_snd_buffer, chunks(chunk)%soundspeed_left_bottom_rcv_buffer, chunks(chunk)%field%soundspeed)
        CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%soundspeed_right_top_snd_buffer, chunks(chunk)%field%soundspeed)

        chunks(right_top_neighbour_chunk)[receiver]%soundspeed_left_bottom_rcv_buffer(1:size) = chunks(chunk)%soundspeed_right_top_snd_buffer(1:size)
    ENDIF

    x_inc=1
    y_inc=1
    IF(fields(FIELD_XVEL0).EQ.1) THEN
        !CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, VERTEX_DATA, & 
        !                                             chunks(chunk)%xvel0_right_top_snd_buffer, chunks(chunk)%xvel0_left_bottom_rcv_buffer, chunks(chunk)%field%xvel0)
        CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%xvel0_right_top_snd_buffer, chunks(chunk)%field%xvel0)

        chunks(right_top_neighbour_chunk)[receiver]%xvel0_left_bottom_rcv_buffer(1:size) = chunks(chunk)%xvel0_right_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        !CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, VERTEX_DATA, & 
        !                                             chunks(chunk)%xvel1_right_top_snd_buffer, chunks(chunk)%xvel1_left_bottom_rcv_buffer, chunks(chunk)%field%xvel1)
        CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%xvel1_right_top_snd_buffer, chunks(chunk)%field%xvel1)

        chunks(right_top_neighbour_chunk)[receiver]%xvel1_left_bottom_rcv_buffer(1:size) = chunks(chunk)%xvel1_right_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        !CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, VERTEX_DATA, & 
        !                                             chunks(chunk)%yvel0_right_top_snd_buffer, chunks(chunk)%yvel0_left_bottom_rcv_buffer, chunks(chunk)%field%yvel0)
        CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%yvel0_right_top_snd_buffer, chunks(chunk)%field%yvel0)

        chunks(right_top_neighbour_chunk)[receiver]%yvel0_left_bottom_rcv_buffer(1:size) = chunks(chunk)%yvel0_right_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        !CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, VERTEX_DATA, & 
        !                                             chunks(chunk)%yvel1_right_top_snd_buffer, chunks(chunk)%yvel1_left_bottom_rcv_buffer, chunks(chunk)%field%yvel1)
        CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%yvel1_right_top_snd_buffer, chunks(chunk)%field%yvel1)

        chunks(right_top_neighbour_chunk)[receiver]%yvel1_left_bottom_rcv_buffer(1:size) = chunks(chunk)%yvel1_right_top_snd_buffer(1:size)
    ENDIF

    x_inc=1
    y_inc=0
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        !CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, X_FACE_DATA, & 
        !                                             chunks(chunk)%volflux_x_right_top_snd_buffer, chunks(chunk)%volflux_x_left_bottom_rcv_buffer, chunks(chunk)%field%vol_flux_x)
        CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%volflux_x_right_top_snd_buffer, chunks(chunk)%field%vol_flux_x)

        chunks(right_top_neighbour_chunk)[receiver]%volflux_x_left_bottom_rcv_buffer(1:size) = chunks(chunk)%volflux_x_right_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        !CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, X_FACE_DATA, & 
        !                                             chunks(chunk)%massflux_x_right_top_snd_buffer, chunks(chunk)%massflux_x_left_bottom_rcv_buffer, chunks(chunk)%field%mass_flux_x)
        CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%massflux_x_right_top_snd_buffer, chunks(chunk)%field%mass_flux_x)

        chunks(right_top_neighbour_chunk)[receiver]%massflux_x_left_bottom_rcv_buffer(1:size) = chunks(chunk)%massflux_x_right_top_snd_buffer(1:size)
    ENDIF

    x_inc=0
    y_inc=1
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        !CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, Y_FACE_DATA, & 
        !                                             chunks(chunk)%volflux_y_right_top_snd_buffer, chunks(chunk)%volflux_y_left_bottom_rcv_buffer, chunks(chunk)%field%vol_flux_y)
        CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%volflux_y_right_top_snd_buffer, chunks(chunk)%field%vol_flux_y)

        chunks(right_top_neighbour_chunk)[receiver]%volflux_y_left_bottom_rcv_buffer(1:size) = chunks(chunk)%volflux_y_right_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        !CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, Y_FACE_DATA, & 
        !                                             chunks(chunk)%massflux_y_right_top_snd_buffer, chunks(chunk)%massflux_y_left_bottom_rcv_buffer, chunks(chunk)%field%mass_flux_y)
        CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%massflux_y_right_top_snd_buffer, chunks(chunk)%field%mass_flux_y)

        chunks(right_top_neighbour_chunk)[receiver]%massflux_y_left_bottom_rcv_buffer(1:size) = chunks(chunk)%massflux_y_right_top_snd_buffer(1:size)
    ENDIF

END SUBROUTINE clover_exchange_write_all_buffers_right_top

SUBROUTINE clover_exchange_write_all_buffers_right_bottom(chunk, depth, fields)

    IMPLICIT NONE 

    INTEGER :: size, chunk, depth, receiver, topedge, bottomedge, fields(NUM_FIELDS), right_bottom_neighbour_chunk, x_inc, y_inc

    right_bottom_neighbour_chunk = chunks(chunk)%chunk_neighbours(chunk_right_bottom)
    receiver = chunks(right_bottom_neighbour_chunk)%task + 1

    size = depth*depth

    x_inc=0
    y_inc=0

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
        !CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, CELL_DATA, & 
        !                                                chunks(chunk)%density0_right_bottom_snd_buffer, chunks(chunk)%density0_left_top_rcv_buffer, chunks(chunk)%field%density0)
        CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%density0_right_bottom_snd_buffer, chunks(chunk)%field%density0)

        chunks(right_bottom_neighbour_chunk)[receiver]%density0_left_top_rcv_buffer(1:size) = chunks(chunk)%density0_right_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        !CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, CELL_DATA, & 
        !                                                chunks(chunk)%density1_right_bottom_snd_buffer, chunks(chunk)%density1_left_top_rcv_buffer, chunks(chunk)%field%density1)
        CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%density1_right_bottom_snd_buffer, chunks(chunk)%field%density1)

        chunks(right_bottom_neighbour_chunk)[receiver]%density1_left_top_rcv_buffer(1:size) = chunks(chunk)%density1_right_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        !CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, CELL_DATA, & 
        !                                                chunks(chunk)%energy0_right_bottom_snd_buffer, chunks(chunk)%energy0_left_top_rcv_buffer, chunks(chunk)%field%energy0)
        CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%energy0_right_bottom_snd_buffer, chunks(chunk)%field%energy0)

        chunks(right_bottom_neighbour_chunk)[receiver]%energy0_left_top_rcv_buffer(1:size) = chunks(chunk)%energy0_right_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        !CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, CELL_DATA, & 
        !                                                chunks(chunk)%energy1_right_bottom_snd_buffer, chunks(chunk)%energy1_left_top_rcv_buffer, chunks(chunk)%field%energy1)
        CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%energy1_right_bottom_snd_buffer, chunks(chunk)%field%energy1)

        chunks(right_bottom_neighbour_chunk)[receiver]%energy1_left_top_rcv_buffer(1:size) = chunks(chunk)%energy1_right_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        !CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, CELL_DATA, & 
        !                                                chunks(chunk)%pressure_right_bottom_snd_buffer, chunks(chunk)%pressure_left_top_rcv_buffer, chunks(chunk)%field%pressure)
        CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%pressure_right_bottom_snd_buffer, chunks(chunk)%field%pressure)

        chunks(right_bottom_neighbour_chunk)[receiver]%pressure_left_top_rcv_buffer(1:size) = chunks(chunk)%pressure_right_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        !CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, CELL_DATA, & 
        !                                                chunks(chunk)%viscosity_right_bottom_snd_buffer, chunks(chunk)%viscosity_left_top_rcv_buffer, chunks(chunk)%field%viscosity)
        CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%viscosity_right_bottom_snd_buffer, chunks(chunk)%field%viscosity)

        chunks(right_bottom_neighbour_chunk)[receiver]%viscosity_left_top_rcv_buffer(1:size) = chunks(chunk)%viscosity_right_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        !CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, CELL_DATA, & 
        !                                                chunks(chunk)%soundspeed_right_bottom_snd_buffer, chunks(chunk)%soundspeed_left_top_rcv_buffer, chunks(chunk)%field%soundspeed)
        CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%soundspeed_right_bottom_snd_buffer, chunks(chunk)%field%soundspeed)

        chunks(right_bottom_neighbour_chunk)[receiver]%soundspeed_left_top_rcv_buffer(1:size) = chunks(chunk)%soundspeed_right_bottom_snd_buffer(1:size)
    ENDIF

    x_inc=1
    y_inc=1
    IF(fields(FIELD_XVEL0).EQ.1) THEN
        !CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, VERTEX_DATA, & 
        !                                                chunks(chunk)%xvel0_right_bottom_snd_buffer, chunks(chunk)%xvel0_left_top_rcv_buffer, chunks(chunk)%field%xvel0)
        CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%xvel0_right_bottom_snd_buffer, chunks(chunk)%field%xvel0)

        chunks(right_bottom_neighbour_chunk)[receiver]%xvel0_left_top_rcv_buffer(1:size) = chunks(chunk)%xvel0_right_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        !CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, VERTEX_DATA, & 
        !                                                chunks(chunk)%xvel1_right_bottom_snd_buffer, chunks(chunk)%xvel1_left_top_rcv_buffer, chunks(chunk)%field%xvel1)
        CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%xvel1_right_bottom_snd_buffer, chunks(chunk)%field%xvel1)

        chunks(right_bottom_neighbour_chunk)[receiver]%xvel1_left_top_rcv_buffer(1:size) = chunks(chunk)%xvel1_right_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        !CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, VERTEX_DATA, & 
        !                                                chunks(chunk)%yvel0_right_bottom_snd_buffer, chunks(chunk)%yvel0_left_top_rcv_buffer, chunks(chunk)%field%yvel0)
        CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%yvel0_right_bottom_snd_buffer, chunks(chunk)%field%yvel0)

        chunks(right_bottom_neighbour_chunk)[receiver]%yvel0_left_top_rcv_buffer(1:size) = chunks(chunk)%yvel0_right_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        !CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, VERTEX_DATA, & 
        !                                                chunks(chunk)%yvel1_right_bottom_snd_buffer, chunks(chunk)%yvel1_left_top_rcv_buffer, chunks(chunk)%field%yvel1)
        CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%yvel1_right_bottom_snd_buffer, chunks(chunk)%field%yvel1)

        chunks(right_bottom_neighbour_chunk)[receiver]%yvel1_left_top_rcv_buffer(1:size) = chunks(chunk)%yvel1_right_bottom_snd_buffer(1:size)
    ENDIF

    x_inc=1
    y_inc=0
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        !CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, X_FACE_DATA, & 
        !                                                chunks(chunk)%volflux_x_right_bottom_snd_buffer, chunks(chunk)%volflux_x_left_top_rcv_buffer, chunks(chunk)%field%vol_flux_x)
        CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%volflux_x_right_bottom_snd_buffer, chunks(chunk)%field%vol_flux_x)

        chunks(right_bottom_neighbour_chunk)[receiver]%volflux_x_left_top_rcv_buffer(1:size) = chunks(chunk)%volflux_x_right_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        !CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, X_FACE_DATA, & 
        !                                                chunks(chunk)%massflux_x_right_bottom_snd_buffer, chunks(chunk)%massflux_x_left_top_rcv_buffer, chunks(chunk)%field%mass_flux_x)
        CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%massflux_x_right_bottom_snd_buffer, chunks(chunk)%field%mass_flux_x)

        chunks(right_bottom_neighbour_chunk)[receiver]%massflux_x_left_top_rcv_buffer(1:size) = chunks(chunk)%massflux_x_right_bottom_snd_buffer(1:size)
    ENDIF

    x_inc=0
    y_inc=1
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        !CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, Y_FACE_DATA, & 
        !                                                chunks(chunk)%volflux_y_right_bottom_snd_buffer, chunks(chunk)%volflux_y_left_top_rcv_buffer, chunks(chunk)%field%vol_flux_y)
        CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%volflux_y_right_bottom_snd_buffer, chunks(chunk)%field%vol_flux_y)

        chunks(right_bottom_neighbour_chunk)[receiver]%volflux_y_left_top_rcv_buffer(1:size) = chunks(chunk)%volflux_y_right_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        !CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, Y_FACE_DATA, & 
        !                                                chunks(chunk)%massflux_y_right_bottom_snd_buffer, chunks(chunk)%massflux_y_left_top_rcv_buffer, chunks(chunk)%field%mass_flux_y)
        CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%massflux_y_right_bottom_snd_buffer, chunks(chunk)%field%mass_flux_y)

        chunks(right_bottom_neighbour_chunk)[receiver]%massflux_y_left_top_rcv_buffer(1:size) = chunks(chunk)%massflux_y_right_bottom_snd_buffer(1:size)
    ENDIF

END SUBROUTINE clover_exchange_write_all_buffers_right_bottom

SUBROUTINE clover_exchange_write_all_buffers_left_bottom(chunk, depth, fields)

    IMPLICIT NONE 

    INTEGER :: size, chunk, depth, receiver, topedge, bottomedge, fields(NUM_FIELDS), left_bottom_neighbour_chunk, x_inc, y_inc

    left_bottom_neighbour_chunk = chunks(chunk)%chunk_neighbours(chunk_left_bottom)
    receiver = chunks(left_bottom_neighbour_chunk)%task + 1

    x_inc=0
    y_inc=0
    size = depth*depth

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
        !CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, CELL_DATA, & 
        !                                               chunks(chunk)%density0_left_bottom_snd_buffer, chunks(chunk)%density0_right_top_rcv_buffer, chunks(chunk)%field%density0)
        CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%density0_left_bottom_snd_buffer, chunks(chunk)%field%density0)

        chunks(left_bottom_neighbour_chunk)[receiver]%density0_right_top_rcv_buffer(1:size) = chunks(chunk)%density0_left_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        !CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, CELL_DATA, & 
        !                                               chunks(chunk)%density1_left_bottom_snd_buffer, chunks(chunk)%density1_right_top_rcv_buffer, chunks(chunk)%field%density1)
        CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%density1_left_bottom_snd_buffer, chunks(chunk)%field%density1)

        chunks(left_bottom_neighbour_chunk)[receiver]%density1_right_top_rcv_buffer(1:size) = chunks(chunk)%density1_left_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        !CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, CELL_DATA, & 
        !                                               chunks(chunk)%energy0_left_bottom_snd_buffer, chunks(chunk)%energy0_right_top_rcv_buffer, chunks(chunk)%field%energy0)
        CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%energy0_left_bottom_snd_buffer, chunks(chunk)%field%energy0)

        chunks(left_bottom_neighbour_chunk)[receiver]%energy0_right_top_rcv_buffer(1:size) = chunks(chunk)%energy0_left_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        !CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, CELL_DATA, & 
        !                                               chunks(chunk)%energy1_left_bottom_snd_buffer, chunks(chunk)%energy1_right_top_rcv_buffer, chunks(chunk)%field%energy1)
        CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%energy1_left_bottom_snd_buffer, chunks(chunk)%field%energy1)

        chunks(left_bottom_neighbour_chunk)[receiver]%energy1_right_top_rcv_buffer(1:size) = chunks(chunk)%energy1_left_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        !CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, CELL_DATA, & 
        !                                               chunks(chunk)%pressure_left_bottom_snd_buffer, chunks(chunk)%pressure_right_top_rcv_buffer, chunks(chunk)%field%pressure)
        CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%pressure_left_bottom_snd_buffer, chunks(chunk)%field%pressure)

        chunks(left_bottom_neighbour_chunk)[receiver]%pressure_right_top_rcv_buffer(1:size) = chunks(chunk)%pressure_left_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        !CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, CELL_DATA, & 
        !                                               chunks(chunk)%viscosity_left_bottom_snd_buffer, chunks(chunk)%viscosity_right_top_rcv_buffer, chunks(chunk)%field%viscosity)
        CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%viscosity_left_bottom_snd_buffer, chunks(chunk)%field%viscosity)

        chunks(left_bottom_neighbour_chunk)[receiver]%viscosity_right_top_rcv_buffer(1:size) = chunks(chunk)%viscosity_left_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        !CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, CELL_DATA, & 
        !                                               chunks(chunk)%soundspeed_left_bottom_snd_buffer, chunks(chunk)%soundspeed_right_top_rcv_buffer, chunks(chunk)%field%soundspeed)
        CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%soundspeed_left_bottom_snd_buffer, chunks(chunk)%field%soundspeed)

        chunks(left_bottom_neighbour_chunk)[receiver]%soundspeed_right_top_rcv_buffer(1:size) = chunks(chunk)%soundspeed_left_bottom_snd_buffer(1:size)
    ENDIF

    x_inc=1
    y_inc=1
    IF(fields(FIELD_XVEL0).EQ.1) THEN
        !CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, VERTEX_DATA, & 
        !                                               chunks(chunk)%xvel0_left_bottom_snd_buffer, chunks(chunk)%xvel0_right_top_rcv_buffer, chunks(chunk)%field%xvel0)
        CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%xvel0_left_bottom_snd_buffer, chunks(chunk)%field%xvel0)

        chunks(left_bottom_neighbour_chunk)[receiver]%xvel0_right_top_rcv_buffer(1:size) = chunks(chunk)%xvel0_left_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        !CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, VERTEX_DATA, & 
        !                                               chunks(chunk)%xvel1_left_bottom_snd_buffer, chunks(chunk)%xvel1_right_top_rcv_buffer, chunks(chunk)%field%xvel1)
        CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%xvel1_left_bottom_snd_buffer, chunks(chunk)%field%xvel1)

        chunks(left_bottom_neighbour_chunk)[receiver]%xvel1_right_top_rcv_buffer(1:size) = chunks(chunk)%xvel1_left_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        !CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, VERTEX_DATA, & 
        !                                               chunks(chunk)%yvel0_left_bottom_snd_buffer, chunks(chunk)%yvel0_right_top_rcv_buffer, chunks(chunk)%field%yvel0)
        CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%yvel0_left_bottom_snd_buffer, chunks(chunk)%field%yvel0)

        chunks(left_bottom_neighbour_chunk)[receiver]%yvel0_right_top_rcv_buffer(1:size) = chunks(chunk)%yvel0_left_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        !CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, VERTEX_DATA, & 
        !                                               chunks(chunk)%yvel1_left_bottom_snd_buffer, chunks(chunk)%yvel1_right_top_rcv_buffer, chunks(chunk)%field%yvel1)
        CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%yvel1_left_bottom_snd_buffer, chunks(chunk)%field%yvel1)

        chunks(left_bottom_neighbour_chunk)[receiver]%yvel1_right_top_rcv_buffer(1:size) = chunks(chunk)%yvel1_left_bottom_snd_buffer(1:size)
    ENDIF

    x_inc=1
    y_inc=0
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        !CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, X_FACE_DATA, & 
        !                                               chunks(chunk)%volflux_x_left_bottom_snd_buffer, chunks(chunk)%volflux_x_right_top_rcv_buffer, chunks(chunk)%field%vol_flux_x)
        CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%volflux_x_left_bottom_snd_buffer, chunks(chunk)%field%vol_flux_x)

        chunks(left_bottom_neighbour_chunk)[receiver]%volflux_x_right_top_rcv_buffer(1:size) = chunks(chunk)%volflux_x_left_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        !CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, X_FACE_DATA, & 
        !                                               chunks(chunk)%massflux_x_left_bottom_snd_buffer, chunks(chunk)%massflux_x_right_top_rcv_buffer, chunks(chunk)%field%mass_flux_x)
        CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%massflux_x_left_bottom_snd_buffer, chunks(chunk)%field%mass_flux_x)

        chunks(left_bottom_neighbour_chunk)[receiver]%massflux_x_right_top_rcv_buffer(1:size) = chunks(chunk)%massflux_x_left_bottom_snd_buffer(1:size)
    ENDIF

    x_inc=0
    y_inc=1
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        !CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, Y_FACE_DATA, & 
        !                                               chunks(chunk)%volflux_y_left_bottom_snd_buffer, chunks(chunk)%volflux_y_right_top_rcv_buffer, chunks(chunk)%field%vol_flux_y)
        CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%volflux_y_left_bottom_snd_buffer, chunks(chunk)%field%vol_flux_y)

        chunks(left_bottom_neighbour_chunk)[receiver]%volflux_y_right_top_rcv_buffer(1:size) = chunks(chunk)%volflux_y_left_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        !CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, Y_FACE_DATA, & 
        !                                               chunks(chunk)%massflux_y_left_bottom_snd_buffer, chunks(chunk)%massflux_y_right_top_rcv_buffer, chunks(chunk)%field%mass_flux_y)
        CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%massflux_y_left_bottom_snd_buffer, chunks(chunk)%field%mass_flux_y)

        chunks(left_bottom_neighbour_chunk)[receiver]%massflux_y_right_top_rcv_buffer(1:size) = chunks(chunk)%massflux_y_left_bottom_snd_buffer(1:size)
    ENDIF

END SUBROUTINE clover_exchange_write_all_buffers_left_bottom










!SUBROUTINE clover_exchange(fields,depth)
!
!  IMPLICIT NONE
!
!  INTEGER      :: fields(:),depth
!
!  ! Assuming 1 patch per task, this will be changed
!  ! Also, not packing all fields for each communication, doing one at a time
!
!  IF(fields(FIELD_DENSITY0).EQ.1) THEN
!    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%density0,      &
!                                 chunks(parallel%task+1)%left_snd_buffer,                      &
!                                 chunks(parallel%task+1)%left_rcv_buffer,                      &
!                                 chunks(parallel%task+1)%right_snd_buffer,                     &
!                                 chunks(parallel%task+1)%right_rcv_buffer,                     &
!                                 chunks(parallel%task+1)%bottom_snd_buffer,                    &
!                                 chunks(parallel%task+1)%bottom_rcv_buffer,                    &
!                                 chunks(parallel%task+1)%top_snd_buffer,                       &
!                                 chunks(parallel%task+1)%top_rcv_buffer,                       &
!                                 chunks(parallel%task+1)%left_top_snd_buffer,                  &
!                                 chunks(parallel%task+1)%left_top_rcv_buffer,                  &
!                                 chunks(parallel%task+1)%right_top_snd_buffer,                 &
!                                 chunks(parallel%task+1)%right_top_rcv_buffer,                 &
!                                 chunks(parallel%task+1)%right_bottom_snd_buffer,              &
!                                 chunks(parallel%task+1)%right_bottom_rcv_buffer,              &
!                                 chunks(parallel%task+1)%left_bottom_snd_buffer,               &
!                                 chunks(parallel%task+1)%left_bottom_rcv_buffer,               &
!                                 depth,CELL_DATA)
!  ENDIF
!
!  IF(fields(FIELD_DENSITY1).EQ.1) THEN
!    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%density1,      &
!                                 chunks(parallel%task+1)%left_snd_buffer,                      &
!                                 chunks(parallel%task+1)%left_rcv_buffer,                      &
!                                 chunks(parallel%task+1)%right_snd_buffer,                     &
!                                 chunks(parallel%task+1)%right_rcv_buffer,                     &
!                                 chunks(parallel%task+1)%bottom_snd_buffer,                    &
!                                 chunks(parallel%task+1)%bottom_rcv_buffer,                    &
!                                 chunks(parallel%task+1)%top_snd_buffer,                       &
!                                 chunks(parallel%task+1)%top_rcv_buffer,                       &
!                                 chunks(parallel%task+1)%left_top_snd_buffer,                  &
!                                 chunks(parallel%task+1)%left_top_rcv_buffer,                  &
!                                 chunks(parallel%task+1)%right_top_snd_buffer,                 &
!                                 chunks(parallel%task+1)%right_top_rcv_buffer,                 &
!                                 chunks(parallel%task+1)%right_bottom_snd_buffer,              &
!                                 chunks(parallel%task+1)%right_bottom_rcv_buffer,              &
!                                 chunks(parallel%task+1)%left_bottom_snd_buffer,               &
!                                 chunks(parallel%task+1)%left_bottom_rcv_buffer,               &
!                                 depth,CELL_DATA)
!  ENDIF
!
!  IF(fields(FIELD_ENERGY0).EQ.1) THEN
!    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%energy0,       &
!                                 chunks(parallel%task+1)%left_snd_buffer,                      &
!                                 chunks(parallel%task+1)%left_rcv_buffer,                      &
!                                 chunks(parallel%task+1)%right_snd_buffer,                     &
!                                 chunks(parallel%task+1)%right_rcv_buffer,                     &
!                                 chunks(parallel%task+1)%bottom_snd_buffer,                    &
!                                 chunks(parallel%task+1)%bottom_rcv_buffer,                    &
!                                 chunks(parallel%task+1)%top_snd_buffer,                       &
!                                 chunks(parallel%task+1)%top_rcv_buffer,                       &
!                                 chunks(parallel%task+1)%left_top_snd_buffer,                  &
!                                 chunks(parallel%task+1)%left_top_rcv_buffer,                  &
!                                 chunks(parallel%task+1)%right_top_snd_buffer,                 &
!                                 chunks(parallel%task+1)%right_top_rcv_buffer,                 &
!                                 chunks(parallel%task+1)%right_bottom_snd_buffer,              &
!                                 chunks(parallel%task+1)%right_bottom_rcv_buffer,              &
!                                 chunks(parallel%task+1)%left_bottom_snd_buffer,               &
!                                 chunks(parallel%task+1)%left_bottom_rcv_buffer,               &
!                                 depth,CELL_DATA)
!  ENDIF
!
!  IF(fields(FIELD_ENERGY1).EQ.1) THEN
!    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%energy1,       &
!                                 chunks(parallel%task+1)%left_snd_buffer,                      &
!                                 chunks(parallel%task+1)%left_rcv_buffer,                      &
!                                 chunks(parallel%task+1)%right_snd_buffer,                     &
!                                 chunks(parallel%task+1)%right_rcv_buffer,                     &
!                                 chunks(parallel%task+1)%bottom_snd_buffer,                    &
!                                 chunks(parallel%task+1)%bottom_rcv_buffer,                    &
!                                 chunks(parallel%task+1)%top_snd_buffer,                       &
!                                 chunks(parallel%task+1)%top_rcv_buffer,                       &
!                                 chunks(parallel%task+1)%left_top_snd_buffer,                  &
!                                 chunks(parallel%task+1)%left_top_rcv_buffer,                  &
!                                 chunks(parallel%task+1)%right_top_snd_buffer,                 &
!                                 chunks(parallel%task+1)%right_top_rcv_buffer,                 &
!                                 chunks(parallel%task+1)%right_bottom_snd_buffer,              &
!                                 chunks(parallel%task+1)%right_bottom_rcv_buffer,              &
!                                 chunks(parallel%task+1)%left_bottom_snd_buffer,               &
!                                 chunks(parallel%task+1)%left_bottom_rcv_buffer,               &
!                                 depth,CELL_DATA)
!  ENDIF
!
!  IF(fields(FIELD_PRESSURE).EQ.1) THEN
!    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%pressure,      &
!                                 chunks(parallel%task+1)%left_snd_buffer,                      &
!                                 chunks(parallel%task+1)%left_rcv_buffer,                      &
!                                 chunks(parallel%task+1)%right_snd_buffer,                     &
!                                 chunks(parallel%task+1)%right_rcv_buffer,                     &
!                                 chunks(parallel%task+1)%bottom_snd_buffer,                    &
!                                 chunks(parallel%task+1)%bottom_rcv_buffer,                    &
!                                 chunks(parallel%task+1)%top_snd_buffer,                       &
!                                 chunks(parallel%task+1)%top_rcv_buffer,                       &
!                                 chunks(parallel%task+1)%left_top_snd_buffer,                  &
!                                 chunks(parallel%task+1)%left_top_rcv_buffer,                  &
!                                 chunks(parallel%task+1)%right_top_snd_buffer,                 &
!                                 chunks(parallel%task+1)%right_top_rcv_buffer,                 &
!                                 chunks(parallel%task+1)%right_bottom_snd_buffer,              &
!                                 chunks(parallel%task+1)%right_bottom_rcv_buffer,              &
!                                 chunks(parallel%task+1)%left_bottom_snd_buffer,               &
!                                 chunks(parallel%task+1)%left_bottom_rcv_buffer,               &
!                                 depth,CELL_DATA)
!  ENDIF
!
!  IF(fields(FIELD_VISCOSITY).EQ.1) THEN
!    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%viscosity,     &
!                                 chunks(parallel%task+1)%left_snd_buffer,                      &
!                                 chunks(parallel%task+1)%left_rcv_buffer,                      &
!                                 chunks(parallel%task+1)%right_snd_buffer,                     &
!                                 chunks(parallel%task+1)%right_rcv_buffer,                     &
!                                 chunks(parallel%task+1)%bottom_snd_buffer,                    &
!                                 chunks(parallel%task+1)%bottom_rcv_buffer,                    &
!                                 chunks(parallel%task+1)%top_snd_buffer,                       &
!                                 chunks(parallel%task+1)%top_rcv_buffer,                       &
!                                 chunks(parallel%task+1)%left_top_snd_buffer,                  &
!                                 chunks(parallel%task+1)%left_top_rcv_buffer,                  &
!                                 chunks(parallel%task+1)%right_top_snd_buffer,                 &
!                                 chunks(parallel%task+1)%right_top_rcv_buffer,                 &
!                                 chunks(parallel%task+1)%right_bottom_snd_buffer,              &
!                                 chunks(parallel%task+1)%right_bottom_rcv_buffer,              &
!                                 chunks(parallel%task+1)%left_bottom_snd_buffer,               &
!                                 chunks(parallel%task+1)%left_bottom_rcv_buffer,               &
!                                 depth,CELL_DATA)
!  ENDIF
!
!  IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
!    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%soundspeed,    &
!                                 chunks(parallel%task+1)%left_snd_buffer,                      &
!                                 chunks(parallel%task+1)%left_rcv_buffer,                      &
!                                 chunks(parallel%task+1)%right_snd_buffer,                     &
!                                 chunks(parallel%task+1)%right_rcv_buffer,                     &
!                                 chunks(parallel%task+1)%bottom_snd_buffer,                    &
!                                 chunks(parallel%task+1)%bottom_rcv_buffer,                    &
!                                 chunks(parallel%task+1)%top_snd_buffer,                       &
!                                 chunks(parallel%task+1)%top_rcv_buffer,                       &
!                                 chunks(parallel%task+1)%left_top_snd_buffer,                  &
!                                 chunks(parallel%task+1)%left_top_rcv_buffer,                  &
!                                 chunks(parallel%task+1)%right_top_snd_buffer,                 &
!                                 chunks(parallel%task+1)%right_top_rcv_buffer,                 &
!                                 chunks(parallel%task+1)%right_bottom_snd_buffer,              &
!                                 chunks(parallel%task+1)%right_bottom_rcv_buffer,              &
!                                 chunks(parallel%task+1)%left_bottom_snd_buffer,               &
!                                 chunks(parallel%task+1)%left_bottom_rcv_buffer,               &
!                                 depth,CELL_DATA)
!  ENDIF
!
!  IF(fields(FIELD_XVEL0).EQ.1) THEN
!    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%xvel0,         &
!                                 chunks(parallel%task+1)%left_snd_buffer,                      &
!                                 chunks(parallel%task+1)%left_rcv_buffer,                      &
!                                 chunks(parallel%task+1)%right_snd_buffer,                     &
!                                 chunks(parallel%task+1)%right_rcv_buffer,                     &
!                                 chunks(parallel%task+1)%bottom_snd_buffer,                    &
!                                 chunks(parallel%task+1)%bottom_rcv_buffer,                    &
!                                 chunks(parallel%task+1)%top_snd_buffer,                       &
!                                 chunks(parallel%task+1)%top_rcv_buffer,                       &
!                                 chunks(parallel%task+1)%left_top_snd_buffer,                  &
!                                 chunks(parallel%task+1)%left_top_rcv_buffer,                  &
!                                 chunks(parallel%task+1)%right_top_snd_buffer,                 &
!                                 chunks(parallel%task+1)%right_top_rcv_buffer,                 &
!                                 chunks(parallel%task+1)%right_bottom_snd_buffer,              &
!                                 chunks(parallel%task+1)%right_bottom_rcv_buffer,              &
!                                 chunks(parallel%task+1)%left_bottom_snd_buffer,               &
!                                 chunks(parallel%task+1)%left_bottom_rcv_buffer,               &
!                                 depth,VERTEX_DATA)
!  ENDIF
!
!  IF(fields(FIELD_XVEL1).EQ.1) THEN
!    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%xvel1,         &
!                                 chunks(parallel%task+1)%left_snd_buffer,                      &
!                                 chunks(parallel%task+1)%left_rcv_buffer,                      &
!                                 chunks(parallel%task+1)%right_snd_buffer,                     &
!                                 chunks(parallel%task+1)%right_rcv_buffer,                     &
!                                 chunks(parallel%task+1)%bottom_snd_buffer,                    &
!                                 chunks(parallel%task+1)%bottom_rcv_buffer,                    &
!                                 chunks(parallel%task+1)%top_snd_buffer,                       &
!                                 chunks(parallel%task+1)%top_rcv_buffer,                       &
!                                 chunks(parallel%task+1)%left_top_snd_buffer,                  &
!                                 chunks(parallel%task+1)%left_top_rcv_buffer,                  &
!                                 chunks(parallel%task+1)%right_top_snd_buffer,                 &
!                                 chunks(parallel%task+1)%right_top_rcv_buffer,                 &
!                                 chunks(parallel%task+1)%right_bottom_snd_buffer,              &
!                                 chunks(parallel%task+1)%right_bottom_rcv_buffer,              &
!                                 chunks(parallel%task+1)%left_bottom_snd_buffer,               &
!                                 chunks(parallel%task+1)%left_bottom_rcv_buffer,               &
!                                 depth,VERTEX_DATA)
!  ENDIF
!
!  IF(fields(FIELD_YVEL0).EQ.1) THEN
!    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%yvel0,         &
!                                 chunks(parallel%task+1)%left_snd_buffer,                      &
!                                 chunks(parallel%task+1)%left_rcv_buffer,                      &
!                                 chunks(parallel%task+1)%right_snd_buffer,                     &
!                                 chunks(parallel%task+1)%right_rcv_buffer,                     &
!                                 chunks(parallel%task+1)%bottom_snd_buffer,                    &
!                                 chunks(parallel%task+1)%bottom_rcv_buffer,                    &
!                                 chunks(parallel%task+1)%top_snd_buffer,                       &
!                                 chunks(parallel%task+1)%top_rcv_buffer,                       &
!                                 chunks(parallel%task+1)%left_top_snd_buffer,                  &
!                                 chunks(parallel%task+1)%left_top_rcv_buffer,                  &
!                                 chunks(parallel%task+1)%right_top_snd_buffer,                 &
!                                 chunks(parallel%task+1)%right_top_rcv_buffer,                 &
!                                 chunks(parallel%task+1)%right_bottom_snd_buffer,              &
!                                 chunks(parallel%task+1)%right_bottom_rcv_buffer,              &
!                                 chunks(parallel%task+1)%left_bottom_snd_buffer,               &
!                                 chunks(parallel%task+1)%left_bottom_rcv_buffer,               &
!                                 depth,VERTEX_DATA)
!  ENDIF
!
!  IF(fields(FIELD_YVEL1).EQ.1) THEN
!    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%yvel1,         &
!                                 chunks(parallel%task+1)%left_snd_buffer,                      &
!                                 chunks(parallel%task+1)%left_rcv_buffer,                      &
!                                 chunks(parallel%task+1)%right_snd_buffer,                     &
!                                 chunks(parallel%task+1)%right_rcv_buffer,                     &
!                                 chunks(parallel%task+1)%bottom_snd_buffer,                    &
!                                 chunks(parallel%task+1)%bottom_rcv_buffer,                    &
!                                 chunks(parallel%task+1)%top_snd_buffer,                       &
!                                 chunks(parallel%task+1)%top_rcv_buffer,                       &
!                                 chunks(parallel%task+1)%left_top_snd_buffer,                  &
!                                 chunks(parallel%task+1)%left_top_rcv_buffer,                  &
!                                 chunks(parallel%task+1)%right_top_snd_buffer,                 &
!                                 chunks(parallel%task+1)%right_top_rcv_buffer,                 &
!                                 chunks(parallel%task+1)%right_bottom_snd_buffer,              &
!                                 chunks(parallel%task+1)%right_bottom_rcv_buffer,              &
!                                 chunks(parallel%task+1)%left_bottom_snd_buffer,               &
!                                 chunks(parallel%task+1)%left_bottom_rcv_buffer,               &
!                                 depth,VERTEX_DATA)
!  ENDIF
!
!  IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
!    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%vol_flux_x,    &
!                                 chunks(parallel%task+1)%left_snd_buffer,                      &
!                                 chunks(parallel%task+1)%left_rcv_buffer,                      &
!                                 chunks(parallel%task+1)%right_snd_buffer,                     &
!                                 chunks(parallel%task+1)%right_rcv_buffer,                     &
!                                 chunks(parallel%task+1)%bottom_snd_buffer,                    &
!                                 chunks(parallel%task+1)%bottom_rcv_buffer,                    &
!                                 chunks(parallel%task+1)%top_snd_buffer,                       &
!                                 chunks(parallel%task+1)%top_rcv_buffer,                       &
!                                 chunks(parallel%task+1)%left_top_snd_buffer,                  &
!                                 chunks(parallel%task+1)%left_top_rcv_buffer,                  &
!                                 chunks(parallel%task+1)%right_top_snd_buffer,                 &
!                                 chunks(parallel%task+1)%right_top_rcv_buffer,                 &
!                                 chunks(parallel%task+1)%right_bottom_snd_buffer,              &
!                                 chunks(parallel%task+1)%right_bottom_rcv_buffer,              &
!                                 chunks(parallel%task+1)%left_bottom_snd_buffer,               &
!                                 chunks(parallel%task+1)%left_bottom_rcv_buffer,               &
!                                 depth,X_FACE_DATA)
!  ENDIF
!
!  IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
!    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%vol_flux_y,    &
!                                 chunks(parallel%task+1)%left_snd_buffer,                      &
!                                 chunks(parallel%task+1)%left_rcv_buffer,                      &
!                                 chunks(parallel%task+1)%right_snd_buffer,                     &
!                                 chunks(parallel%task+1)%right_rcv_buffer,                     &
!                                 chunks(parallel%task+1)%bottom_snd_buffer,                    &
!                                 chunks(parallel%task+1)%bottom_rcv_buffer,                    &
!                                 chunks(parallel%task+1)%top_snd_buffer,                       &
!                                 chunks(parallel%task+1)%top_rcv_buffer,                       &
!                                 chunks(parallel%task+1)%left_top_snd_buffer,                  &
!                                 chunks(parallel%task+1)%left_top_rcv_buffer,                  &
!                                 chunks(parallel%task+1)%right_top_snd_buffer,                 &
!                                 chunks(parallel%task+1)%right_top_rcv_buffer,                 &
!                                 chunks(parallel%task+1)%right_bottom_snd_buffer,              &
!                                 chunks(parallel%task+1)%right_bottom_rcv_buffer,              &
!                                 chunks(parallel%task+1)%left_bottom_snd_buffer,               &
!                                 chunks(parallel%task+1)%left_bottom_rcv_buffer,               &
!                                 depth,Y_FACE_DATA)
!  ENDIF
!
!  IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
!    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%mass_flux_x,   &
!                                 chunks(parallel%task+1)%left_snd_buffer,                      &
!                                 chunks(parallel%task+1)%left_rcv_buffer,                      &
!                                 chunks(parallel%task+1)%right_snd_buffer,                     &
!                                 chunks(parallel%task+1)%right_rcv_buffer,                     &
!                                 chunks(parallel%task+1)%bottom_snd_buffer,                    &
!                                 chunks(parallel%task+1)%bottom_rcv_buffer,                    &
!                                 chunks(parallel%task+1)%top_snd_buffer,                       &
!                                 chunks(parallel%task+1)%top_rcv_buffer,                       &
!                                 chunks(parallel%task+1)%left_top_snd_buffer,                  &
!                                 chunks(parallel%task+1)%left_top_rcv_buffer,                  &
!                                 chunks(parallel%task+1)%right_top_snd_buffer,                 &
!                                 chunks(parallel%task+1)%right_top_rcv_buffer,                 &
!                                 chunks(parallel%task+1)%right_bottom_snd_buffer,              &
!                                 chunks(parallel%task+1)%right_bottom_rcv_buffer,              &
!                                 chunks(parallel%task+1)%left_bottom_snd_buffer,               &
!                                 chunks(parallel%task+1)%left_bottom_rcv_buffer,               &
!                                 depth,X_FACE_DATA)
!  ENDIF
!
!  IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
!    CALL clover_exchange_message(parallel%task+1,chunks(parallel%task+1)%field%mass_flux_y,   &
!                                 chunks(parallel%task+1)%left_snd_buffer,                      &
!                                 chunks(parallel%task+1)%left_rcv_buffer,                      &
!                                 chunks(parallel%task+1)%right_snd_buffer,                     &
!                                 chunks(parallel%task+1)%right_rcv_buffer,                     &
!                                 chunks(parallel%task+1)%bottom_snd_buffer,                    &
!                                 chunks(parallel%task+1)%bottom_rcv_buffer,                    &
!                                 chunks(parallel%task+1)%top_snd_buffer,                       &
!                                 chunks(parallel%task+1)%top_rcv_buffer,                       &
!                                 chunks(parallel%task+1)%left_top_snd_buffer,                  &
!                                 chunks(parallel%task+1)%left_top_rcv_buffer,                  &
!                                 chunks(parallel%task+1)%right_top_snd_buffer,                 &
!                                 chunks(parallel%task+1)%right_top_rcv_buffer,                 &
!                                 chunks(parallel%task+1)%right_bottom_snd_buffer,              &
!                                 chunks(parallel%task+1)%right_bottom_rcv_buffer,              &
!                                 chunks(parallel%task+1)%left_bottom_snd_buffer,               &
!                                 chunks(parallel%task+1)%left_bottom_rcv_buffer,               &
!                                 depth,Y_FACE_DATA)
!  ENDIF
!
!
!END SUBROUTINE clover_exchange

!SUBROUTINE clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, field_type, left_snd_buffer, right_rcv_buffer, field)
!
!    IMPLICIT NONE
!
!    INTEGER :: chunk, depth, receiver, size, field_type, x_inc, y_inc, topedge, bottomedge, left_neighbour_chunk
!    REAL(KIND=8) :: left_snd_buffer(:), right_rcv_buffer(:)
!    REAL(KIND=8) :: field(-1:,-1:)
!
!    IF(field_type.EQ.CELL_DATA) THEN
!        x_inc=0
!        y_inc=0
!    ENDIF
!    IF(field_type.EQ.VERTEX_DATA) THEN
!        x_inc=1
!        y_inc=1
!    ENDIF
!    IF(field_type.EQ.X_FACE_DATA) THEN
!        x_inc=1
!        y_inc=0
!    ENDIF
!    IF(field_type.EQ.Y_FACE_DATA) THEN
!        x_inc=0
!        y_inc=1
!    ENDIF
!
!
!    size=(1+(chunks(chunk)%field%y_max+y_inc+topedge)-(chunks(chunk)%field%y_min-bottomedge))*depth
!
!    CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, left_snd_buffer, field)
!
!    left_neighbour_chunk = chunks(chunk)%chunk_neighbours(chunk_left)
!    !receiver=chunks(left_neighbour_chunk)%task
!    chunks(left_neighbour_chunk)[receiver]%right_rcv_buffer(1:size) = left_snd_buffer(1:size)
!
!END SUBROUTINE clover_exchange_write_message_left
!
!SUBROUTINE clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, field_type, right_snd_buffer, left_rcv_buffer, field)
!
!    IMPLICIT NONE
!
!    INTEGER :: chunk, depth, receiver, size, field_type, x_inc, y_inc, topedge, bottomedge, right_neighbour_chunk
!    REAL(KIND=8) :: right_snd_buffer(:), left_rcv_buffer(:)
!    REAL(KIND=8) :: field(-1:,-1:)
!
!    IF(field_type.EQ.CELL_DATA) THEN
!        x_inc=0
!        y_inc=0
!    ENDIF
!    IF(field_type.EQ.VERTEX_DATA) THEN
!        x_inc=1
!        y_inc=1
!    ENDIF
!    IF(field_type.EQ.X_FACE_DATA) THEN
!        x_inc=1
!        y_inc=0
!    ENDIF
!    IF(field_type.EQ.Y_FACE_DATA) THEN
!        x_inc=0
!        y_inc=1
!    ENDIF
!
!    size=(1+(chunks(chunk)%field%y_max+y_inc+topedge)-(chunks(chunk)%field%y_min-bottomedge))*depth
!
!    CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, right_snd_buffer, field)
!
!    right_neighbour_chunk = chunks(chunk)%chunk_neighbours(chunk_right)
!    !receiver = chunks(right_neighbour_chunk)%task
!    chunks(right_neighbour_chunk)[receiver]%left_rcv_buffer(1:size) = right_snd_buffer(1:size)
!
!END SUBROUTINE clover_exchange_write_message_right
!
!
!SUBROUTINE clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, field_type, bottom_snd_buffer, top_rcv_buffer, field)
!
!    IMPLICIT NONE
!
!    INTEGER :: chunk, depth, receiver, size, field_type, x_inc, y_inc, leftedge, rightedge, bottom_neighbour_chunk
!    REAL(KIND=8) :: bottom_snd_buffer(:), top_rcv_buffer(:)
!    REAL(KIND=8) :: field(-1:,-1:)
!
!    IF(field_type.EQ.CELL_DATA) THEN
!        x_inc=0
!        y_inc=0
!    ENDIF
!    IF(field_type.EQ.VERTEX_DATA) THEN
!        x_inc=1
!        y_inc=1
!    ENDIF
!    IF(field_type.EQ.X_FACE_DATA) THEN
!        x_inc=1
!        y_inc=0
!    ENDIF
!    IF(field_type.EQ.Y_FACE_DATA) THEN
!        x_inc=0
!        y_inc=1
!    ENDIF
!
!    size=(1+(chunks(chunk)%field%x_max+x_inc+rightedge)-(chunks(chunk)%field%x_min-leftedge))*depth
!
!    CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, bottom_snd_buffer, field)
!
!    bottom_neighbour_chunk = chunks(chunk)%chunk_neighbours(chunk_bottom)
!    !receiver=chunks(bottom_neighbour_chunk)%task
!    chunks(bottom_neighbour_chunk)[receiver]%top_rcv_buffer(1:size) = bottom_snd_buffer(1:size)
!
!END SUBROUTINE clover_exchange_write_message_bottom
!
!SUBROUTINE clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, field_type, top_snd_buffer, bottom_rcv_buffer, field)
!
!    IMPLICIT NONE
!
!    INTEGER :: chunk, depth, receiver, size, field_type, x_inc, y_inc, leftedge, rightedge, top_neighbour_chunk
!    REAL(KIND=8) :: top_snd_buffer(:), bottom_rcv_buffer(:)
!    REAL(KIND=8) :: field(-1:,-1:)
!
!    IF(field_type.EQ.CELL_DATA) THEN
!        x_inc=0
!        y_inc=0
!    ENDIF
!    IF(field_type.EQ.VERTEX_DATA) THEN
!        x_inc=1
!        y_inc=1
!    ENDIF
!    IF(field_type.EQ.X_FACE_DATA) THEN
!        x_inc=1
!        y_inc=0
!    ENDIF
!    IF(field_type.EQ.Y_FACE_DATA) THEN
!        x_inc=0
!        y_inc=1
!    ENDIF
!
!    size=(1+(chunks(chunk)%field%x_max+x_inc+rightedge)-(chunks(chunk)%field%x_min-leftedge))*depth
!
!    CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, top_snd_buffer, field)
!
!    top_neighbour_chunk = chunks(chunk)%chunk_neighbours(chunk_top)
!    !receiver=chunks(top_neighbour_chunk)%task
!    chunks(top_neighbour_chunk)[receiver]%bottom_rcv_buffer(1:size) = top_snd_buffer(1:size)
!
!END SUBROUTINE clover_exchange_write_message_top

!SUBROUTINE clover_exchange_write_message_left_top(chunk, depth, receiver, size, field_type, left_top_snd_buffer, right_bottom_rcv_buffer, field)
!
!    IMPLICIT NONE
!
!    INTEGER :: chunk, depth, receiver, size, field_type, x_inc, y_inc, left_top_neighbour_chunk
!    REAL(KIND=8) :: left_top_snd_buffer(:), right_bottom_rcv_buffer(:)
!    REAL(KIND=8) :: field(-1:,-1:)
!
!    IF(field_type.EQ.CELL_DATA) THEN
!        x_inc=0
!        y_inc=0
!    ENDIF
!    IF(field_type.EQ.VERTEX_DATA) THEN
!        x_inc=1
!        y_inc=1
!    ENDIF
!    IF(field_type.EQ.X_FACE_DATA) THEN
!        x_inc=1
!        y_inc=0
!    ENDIF
!    IF(field_type.EQ.Y_FACE_DATA) THEN
!        x_inc=0
!        y_inc=1
!    ENDIF
!
!    CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, left_top_snd_buffer, field)
!
!    left_top_neighbour_chunk = chunks(chunk)%chunk_neighbours(chunk_left_top)
!    !receiver = chunks(left_top_neighbour_chunk)%task + 1
!    chunks(left_top_neighbour_chunk)[receiver]%right_bottom_rcv_buffer(1:size) = left_top_snd_buffer(1:size)
!
!END SUBROUTINE clover_exchange_write_message_left_top
!
!SUBROUTINE clover_exchange_write_message_right_top(chunk, depth, receiver, size, field_type, right_top_snd_buffer, left_bottom_rcv_buffer, field)
!
!    IMPLICIT NONE
!
!    INTEGER :: chunk, depth, receiver, size, field_type, x_inc, y_inc, right_top_neighbour_chunk
!    REAL(KIND=8) :: right_top_snd_buffer(:), left_bottom_rcv_buffer(:)
!    REAL(KIND=8) :: field(-1:,-1:)
!
!    IF(field_type.EQ.CELL_DATA) THEN
!        x_inc=0
!        y_inc=0
!    ENDIF
!    IF(field_type.EQ.VERTEX_DATA) THEN
!        x_inc=1
!        y_inc=1
!    ENDIF
!    IF(field_type.EQ.X_FACE_DATA) THEN
!        x_inc=1
!        y_inc=0
!    ENDIF
!    IF(field_type.EQ.Y_FACE_DATA) THEN
!        x_inc=0
!        y_inc=1
!    ENDIF
!
!    CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, right_top_snd_buffer, field)
!
!    right_top_neighbour_chunk = chunks(chunk)%chunk_neighbours(chunk_right_top)
!    !receiver = chunks(right_top_neighbour_chunk)%task + 1
!    chunks(right_top_neighbour_chunk)[receiver]%left_bottom_rcv_buffer(1:size) = right_top_snd_buffer(1:size)
!
!END SUBROUTINE clover_exchange_write_message_right_top
!
!SUBROUTINE clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, field_type, right_bottom_snd_buffer, left_top_rcv_buffer, field)
!
!    IMPLICIT NONE
!
!    INTEGER :: chunk, depth, receiver, size, field_type, x_inc, y_inc, right_bottom_neighbour_chunk
!    REAL(KIND=8) :: right_bottom_snd_buffer(:), left_top_rcv_buffer(:)
!    REAL(KIND=8) :: field(-1:,-1:)
!
!    IF(field_type.EQ.CELL_DATA) THEN
!        x_inc=0
!        y_inc=0
!    ENDIF
!    IF(field_type.EQ.VERTEX_DATA) THEN
!        x_inc=1
!        y_inc=1
!    ENDIF
!    IF(field_type.EQ.X_FACE_DATA) THEN
!        x_inc=1
!        y_inc=0
!    ENDIF
!    IF(field_type.EQ.Y_FACE_DATA) THEN
!        x_inc=0
!        y_inc=1
!    ENDIF
!
!    CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, right_bottom_snd_buffer, field)
!
!    right_bottom_neighbour_chunk = chunks(chunk)%chunk_neighbours(chunk_right_bottom)
!    !receiver = chunks(right_bottom_neighbour_chunk)%task + 1
!    chunks(right_bottom_neighbour_chunk)[receiver]%left_top_rcv_buffer(1:size) = right_bottom_snd_buffer(1:size)
!
!END SUBROUTINE clover_exchange_write_message_right_bottom
!
!SUBROUTINE clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, field_type, left_bottom_snd_buffer, right_top_rcv_buffer, field)
!
!    IMPLICIT NONE
!
!    INTEGER :: chunk, depth, receiver, size, field_type, x_inc, y_inc, left_bottom_neighbour_chunk
!    REAL(KIND=8) :: left_bottom_snd_buffer(:), right_top_rcv_buffer(:)
!    REAL(KIND=8) :: field(-1:,-1:)
!
!    IF(field_type.EQ.CELL_DATA) THEN
!        x_inc=0
!        y_inc=0
!    ENDIF
!    IF(field_type.EQ.VERTEX_DATA) THEN
!        x_inc=1
!        y_inc=1
!    ENDIF
!    IF(field_type.EQ.X_FACE_DATA) THEN
!        x_inc=1
!        y_inc=0
!    ENDIF
!    IF(field_type.EQ.Y_FACE_DATA) THEN
!        x_inc=0
!        y_inc=1
!    ENDIF
!
!    CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, left_bottom_snd_buffer, field)
!
!    left_bottom_neighbour_chunk = chunks(chunk)%chunk_neighbours(chunk_left_bottom)
!    !receiver = chunks(left_bottom_neighbour_chunk)%task + 1
!    chunks(left_bottom_neighbour_chunk)[receiver]%right_top_rcv_buffer(1:size) = left_bottom_snd_buffer(1:size)
!
!END SUBROUTINE clover_exchange_write_message_left_bottom



SUBROUTINE clover_exchange_unpack_all_buffers_left(chunk, depth, fields)

    IMPLICIT NONE 

    INTEGER :: chunk, depth, fields(NUM_FIELDS)

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
        CALL unpack_left_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%density0_left_rcv_buffer, chunks(chunk)%field%density0)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        CALL unpack_left_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%density1_left_rcv_buffer, chunks(chunk)%field%density1)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        CALL unpack_left_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%energy0_left_rcv_buffer, chunks(chunk)%field%energy0)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        CALL unpack_left_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%energy1_left_rcv_buffer, chunks(chunk)%field%energy1)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        CALL unpack_left_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%pressure_left_rcv_buffer, chunks(chunk)%field%pressure)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        CALL unpack_left_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%viscosity_left_rcv_buffer, chunks(chunk)%field%viscosity)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        CALL unpack_left_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%soundspeed_left_rcv_buffer, chunks(chunk)%field%soundspeed)
    ENDIF
    IF(fields(FIELD_XVEL0).EQ.1) THEN
        CALL unpack_left_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%xvel0_left_rcv_buffer, chunks(chunk)%field%xvel0)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        CALL unpack_left_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%xvel1_left_rcv_buffer, chunks(chunk)%field%xvel1)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        CALL unpack_left_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%xvel0_left_rcv_buffer, chunks(chunk)%field%yvel0)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        CALL unpack_left_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%yvel1_left_rcv_buffer, chunks(chunk)%field%yvel1)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        CALL unpack_left_buffer_seq(chunk, depth, X_FACE_DATA, chunks(chunk)%volflux_x_left_rcv_buffer, chunks(chunk)%field%vol_flux_x)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        CALL unpack_left_buffer_seq(chunk, depth, Y_FACE_DATA, chunks(chunk)%volflux_y_left_rcv_buffer, chunks(chunk)%field%vol_flux_y)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        CALL unpack_left_buffer_seq(chunk, depth, X_FACE_DATA, chunks(chunk)%massflux_x_left_rcv_buffer, chunks(chunk)%field%mass_flux_x)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        CALL unpack_left_buffer_seq(chunk, depth, Y_FACE_DATA, chunks(chunk)%massflux_y_left_rcv_buffer, chunks(chunk)%field%mass_flux_y)
    ENDIF

END SUBROUTINE clover_exchange_unpack_all_buffers_left

SUBROUTINE clover_exchange_unpack_all_buffers_right(chunk, depth, fields)

    IMPLICIT NONE 

    INTEGER :: chunk, depth, fields(NUM_FIELDS)

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
        CALL unpack_right_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%density0_right_rcv_buffer, chunks(chunk)%field%density0)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        CALL unpack_right_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%density1_right_rcv_buffer, chunks(chunk)%field%density1)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        CALL unpack_right_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%energy0_right_rcv_buffer, chunks(chunk)%field%energy0)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        CALL unpack_right_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%energy1_right_rcv_buffer, chunks(chunk)%field%energy1)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        CALL unpack_right_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%pressure_right_rcv_buffer, chunks(chunk)%field%pressure)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        CALL unpack_right_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%viscosity_right_rcv_buffer, chunks(chunk)%field%viscosity)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        CALL unpack_right_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%soundspeed_right_rcv_buffer, chunks(chunk)%field%soundspeed)
    ENDIF
    IF(fields(FIELD_XVEL0).EQ.1) THEN
        CALL unpack_right_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%xvel0_right_rcv_buffer, chunks(chunk)%field%xvel0)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        CALL unpack_right_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%xvel1_right_rcv_buffer, chunks(chunk)%field%xvel1)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        CALL unpack_right_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%xvel0_right_rcv_buffer, chunks(chunk)%field%yvel0)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        CALL unpack_right_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%yvel1_right_rcv_buffer, chunks(chunk)%field%yvel1)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        CALL unpack_right_buffer_seq(chunk, depth, X_FACE_DATA, chunks(chunk)%volflux_x_right_rcv_buffer, chunks(chunk)%field%vol_flux_x)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        CALL unpack_right_buffer_seq(chunk, depth, Y_FACE_DATA, chunks(chunk)%volflux_y_right_rcv_buffer, chunks(chunk)%field%vol_flux_y)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        CALL unpack_right_buffer_seq(chunk, depth, X_FACE_DATA, chunks(chunk)%massflux_x_right_rcv_buffer, chunks(chunk)%field%mass_flux_x)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        CALL unpack_right_buffer_seq(chunk, depth, Y_FACE_DATA, chunks(chunk)%massflux_y_right_rcv_buffer, chunks(chunk)%field%mass_flux_y)
    ENDIF

END SUBROUTINE clover_exchange_unpack_all_buffers_right


SUBROUTINE clover_exchange_unpack_all_buffers_bottom(chunk, depth, fields)

    IMPLICIT NONE 

    INTEGER :: chunk, depth, fields(NUM_FIELDS)

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
        CALL unpack_bottom_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%density0_bottom_rcv_buffer, chunks(chunk)%field%density0)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        CALL unpack_bottom_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%density1_bottom_rcv_buffer, chunks(chunk)%field%density1)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        CALL unpack_bottom_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%energy0_bottom_rcv_buffer, chunks(chunk)%field%energy0)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        CALL unpack_bottom_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%energy1_bottom_rcv_buffer, chunks(chunk)%field%energy1)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        CALL unpack_bottom_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%pressure_bottom_rcv_buffer, chunks(chunk)%field%pressure)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        CALL unpack_bottom_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%viscosity_bottom_rcv_buffer, chunks(chunk)%field%viscosity)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        CALL unpack_bottom_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%soundspeed_bottom_rcv_buffer, chunks(chunk)%field%soundspeed)
    ENDIF
    IF(fields(FIELD_XVEL0).EQ.1) THEN
        CALL unpack_bottom_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%xvel0_bottom_rcv_buffer, chunks(chunk)%field%xvel0)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        CALL unpack_bottom_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%xvel1_bottom_rcv_buffer, chunks(chunk)%field%xvel1)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        CALL unpack_bottom_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%xvel0_bottom_rcv_buffer, chunks(chunk)%field%yvel0)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        CALL unpack_bottom_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%yvel1_bottom_rcv_buffer, chunks(chunk)%field%yvel1)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        CALL unpack_bottom_buffer_seq(chunk, depth, X_FACE_DATA, chunks(chunk)%volflux_x_bottom_rcv_buffer, chunks(chunk)%field%vol_flux_x)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        CALL unpack_bottom_buffer_seq(chunk, depth, Y_FACE_DATA, chunks(chunk)%volflux_y_bottom_rcv_buffer, chunks(chunk)%field%vol_flux_y)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        CALL unpack_bottom_buffer_seq(chunk, depth, X_FACE_DATA, chunks(chunk)%massflux_x_bottom_rcv_buffer, chunks(chunk)%field%mass_flux_x)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        CALL unpack_bottom_buffer_seq(chunk, depth, Y_FACE_DATA, chunks(chunk)%massflux_y_bottom_rcv_buffer, chunks(chunk)%field%mass_flux_y)
    ENDIF

END SUBROUTINE clover_exchange_unpack_all_buffers_bottom

SUBROUTINE clover_exchange_unpack_all_buffers_top(chunk, depth, fields)

    IMPLICIT NONE 

    INTEGER :: chunk, depth, fields(NUM_FIELDS)

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
        CALL unpack_top_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%density0_top_rcv_buffer, chunks(chunk)%field%density0)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        CALL unpack_top_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%density1_top_rcv_buffer, chunks(chunk)%field%density1)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        CALL unpack_top_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%energy0_top_rcv_buffer, chunks(chunk)%field%energy0)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        CALL unpack_top_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%energy1_top_rcv_buffer, chunks(chunk)%field%energy1)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        CALL unpack_top_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%pressure_top_rcv_buffer, chunks(chunk)%field%pressure)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        CALL unpack_top_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%viscosity_top_rcv_buffer, chunks(chunk)%field%viscosity)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        CALL unpack_top_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%soundspeed_top_rcv_buffer, chunks(chunk)%field%soundspeed)
    ENDIF
    IF(fields(FIELD_XVEL0).EQ.1) THEN
        CALL unpack_top_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%xvel0_top_rcv_buffer, chunks(chunk)%field%xvel0)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        CALL unpack_top_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%xvel1_top_rcv_buffer, chunks(chunk)%field%xvel1)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        CALL unpack_top_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%xvel0_top_rcv_buffer, chunks(chunk)%field%yvel0)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        CALL unpack_top_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%yvel1_top_rcv_buffer, chunks(chunk)%field%yvel1)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        CALL unpack_top_buffer_seq(chunk, depth, X_FACE_DATA, chunks(chunk)%volflux_x_top_rcv_buffer, chunks(chunk)%field%vol_flux_x)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        CALL unpack_top_buffer_seq(chunk, depth, Y_FACE_DATA, chunks(chunk)%volflux_y_top_rcv_buffer, chunks(chunk)%field%vol_flux_y)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        CALL unpack_top_buffer_seq(chunk, depth, X_FACE_DATA, chunks(chunk)%massflux_x_top_rcv_buffer, chunks(chunk)%field%mass_flux_x)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        CALL unpack_top_buffer_seq(chunk, depth, Y_FACE_DATA, chunks(chunk)%massflux_y_top_rcv_buffer, chunks(chunk)%field%mass_flux_y)
    ENDIF

END SUBROUTINE clover_exchange_unpack_all_buffers_top

SUBROUTINE clover_exchange_unpack_all_buffers_left_top(chunk, depth, fields)

    IMPLICIT NONE 

    INTEGER :: chunk, depth, fields(NUM_FIELDS)

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
        CALL unpack_left_top_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%density0_left_top_rcv_buffer, chunks(chunk)%field%density0)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        CALL unpack_left_top_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%density1_left_top_rcv_buffer, chunks(chunk)%field%density1)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        CALL unpack_left_top_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%energy0_left_top_rcv_buffer, chunks(chunk)%field%energy0)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        CALL unpack_left_top_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%energy1_left_top_rcv_buffer, chunks(chunk)%field%energy1)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        CALL unpack_left_top_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%pressure_left_top_rcv_buffer, chunks(chunk)%field%pressure)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        CALL unpack_left_top_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%viscosity_left_top_rcv_buffer, chunks(chunk)%field%viscosity)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        CALL unpack_left_top_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%soundspeed_left_top_rcv_buffer, chunks(chunk)%field%soundspeed)
    ENDIF
    IF(fields(FIELD_XVEL0).EQ.1) THEN
        CALL unpack_left_top_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%xvel0_left_top_rcv_buffer, chunks(chunk)%field%xvel0)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        CALL unpack_left_top_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%xvel1_left_top_rcv_buffer, chunks(chunk)%field%xvel1)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        CALL unpack_left_top_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%xvel0_left_top_rcv_buffer, chunks(chunk)%field%yvel0)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        CALL unpack_left_top_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%yvel1_left_top_rcv_buffer, chunks(chunk)%field%yvel1)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        CALL unpack_left_top_buffer_seq(chunk, depth, X_FACE_DATA, chunks(chunk)%volflux_x_left_top_rcv_buffer, chunks(chunk)%field%vol_flux_x)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        CALL unpack_left_top_buffer_seq(chunk, depth, Y_FACE_DATA, chunks(chunk)%volflux_y_left_top_rcv_buffer, chunks(chunk)%field%vol_flux_y)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        CALL unpack_left_top_buffer_seq(chunk, depth, X_FACE_DATA, chunks(chunk)%massflux_x_left_top_rcv_buffer, chunks(chunk)%field%mass_flux_x)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        CALL unpack_left_top_buffer_seq(chunk, depth, Y_FACE_DATA, chunks(chunk)%massflux_y_left_top_rcv_buffer, chunks(chunk)%field%mass_flux_y)
    ENDIF

END SUBROUTINE clover_exchange_unpack_all_buffers_left_top

SUBROUTINE clover_exchange_unpack_all_buffers_right_top(chunk, depth, fields)

    IMPLICIT NONE 

    INTEGER :: chunk, depth, fields(NUM_FIELDS)

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
        CALL unpack_right_top_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%density0_right_top_rcv_buffer, chunks(chunk)%field%density0)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        CALL unpack_right_top_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%density1_right_top_rcv_buffer, chunks(chunk)%field%density1)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        CALL unpack_right_top_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%energy0_right_top_rcv_buffer, chunks(chunk)%field%energy0)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        CALL unpack_right_top_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%energy1_right_top_rcv_buffer, chunks(chunk)%field%energy1)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        CALL unpack_right_top_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%pressure_right_top_rcv_buffer, chunks(chunk)%field%pressure)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        CALL unpack_right_top_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%viscosity_right_top_rcv_buffer, chunks(chunk)%field%viscosity)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        CALL unpack_right_top_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%soundspeed_right_top_rcv_buffer, chunks(chunk)%field%soundspeed)
    ENDIF
    IF(fields(FIELD_XVEL0).EQ.1) THEN
        CALL unpack_right_top_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%xvel0_right_top_rcv_buffer, chunks(chunk)%field%xvel0)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        CALL unpack_right_top_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%xvel1_right_top_rcv_buffer, chunks(chunk)%field%xvel1)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        CALL unpack_right_top_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%xvel0_right_top_rcv_buffer, chunks(chunk)%field%yvel0)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        CALL unpack_right_top_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%yvel1_right_top_rcv_buffer, chunks(chunk)%field%yvel1)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        CALL unpack_right_top_buffer_seq(chunk, depth, X_FACE_DATA, chunks(chunk)%volflux_x_right_top_rcv_buffer, chunks(chunk)%field%vol_flux_x)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        CALL unpack_right_top_buffer_seq(chunk, depth, Y_FACE_DATA, chunks(chunk)%volflux_y_right_top_rcv_buffer, chunks(chunk)%field%vol_flux_y)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        CALL unpack_right_top_buffer_seq(chunk, depth, X_FACE_DATA, chunks(chunk)%massflux_x_right_top_rcv_buffer, chunks(chunk)%field%mass_flux_x)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        CALL unpack_right_top_buffer_seq(chunk, depth, Y_FACE_DATA, chunks(chunk)%massflux_y_right_top_rcv_buffer, chunks(chunk)%field%mass_flux_y)
    ENDIF

END SUBROUTINE clover_exchange_unpack_all_buffers_right_top

SUBROUTINE clover_exchange_unpack_all_buffers_right_bottom(chunk, depth, fields)

    IMPLICIT NONE 

    INTEGER :: chunk, depth, fields(NUM_FIELDS)

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
        CALL unpack_right_bottom_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%density0_right_bottom_rcv_buffer, chunks(chunk)%field%density0)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        CALL unpack_right_bottom_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%density1_right_bottom_rcv_buffer, chunks(chunk)%field%density1)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        CALL unpack_right_bottom_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%energy0_right_bottom_rcv_buffer, chunks(chunk)%field%energy0)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        CALL unpack_right_bottom_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%energy1_right_bottom_rcv_buffer, chunks(chunk)%field%energy1)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        CALL unpack_right_bottom_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%pressure_right_bottom_rcv_buffer, chunks(chunk)%field%pressure)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        CALL unpack_right_bottom_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%viscosity_right_bottom_rcv_buffer, chunks(chunk)%field%viscosity)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        CALL unpack_right_bottom_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%soundspeed_right_bottom_rcv_buffer, chunks(chunk)%field%soundspeed)
    ENDIF
    IF(fields(FIELD_XVEL0).EQ.1) THEN
        CALL unpack_right_bottom_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%xvel0_right_bottom_rcv_buffer, chunks(chunk)%field%xvel0)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        CALL unpack_right_bottom_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%xvel1_right_bottom_rcv_buffer, chunks(chunk)%field%xvel1)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        CALL unpack_right_bottom_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%xvel0_right_bottom_rcv_buffer, chunks(chunk)%field%yvel0)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        CALL unpack_right_bottom_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%yvel1_right_bottom_rcv_buffer, chunks(chunk)%field%yvel1)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        CALL unpack_right_bottom_buffer_seq(chunk, depth, X_FACE_DATA, chunks(chunk)%volflux_x_right_bottom_rcv_buffer, chunks(chunk)%field%vol_flux_x)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        CALL unpack_right_bottom_buffer_seq(chunk, depth, Y_FACE_DATA, chunks(chunk)%volflux_y_right_bottom_rcv_buffer, chunks(chunk)%field%vol_flux_y)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        CALL unpack_right_bottom_buffer_seq(chunk, depth, X_FACE_DATA, chunks(chunk)%massflux_x_right_bottom_rcv_buffer, chunks(chunk)%field%mass_flux_x)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        CALL unpack_right_bottom_buffer_seq(chunk, depth, Y_FACE_DATA, chunks(chunk)%massflux_y_right_bottom_rcv_buffer, chunks(chunk)%field%mass_flux_y)
    ENDIF

END SUBROUTINE clover_exchange_unpack_all_buffers_right_bottom

SUBROUTINE clover_exchange_unpack_all_buffers_left_bottom(chunk, depth, fields)

    IMPLICIT NONE 

    INTEGER :: chunk, depth, fields(NUM_FIELDS)

    IF(fields(FIELD_DENSITY0).EQ.1) THEN
        CALL unpack_left_bottom_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%density0_left_bottom_rcv_buffer, chunks(chunk)%field%density0)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        CALL unpack_left_bottom_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%density1_left_bottom_rcv_buffer, chunks(chunk)%field%density1)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        CALL unpack_left_bottom_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%energy0_left_bottom_rcv_buffer, chunks(chunk)%field%energy0)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        CALL unpack_left_bottom_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%energy1_left_bottom_rcv_buffer, chunks(chunk)%field%energy1)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        CALL unpack_left_bottom_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%pressure_left_bottom_rcv_buffer, chunks(chunk)%field%pressure)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        CALL unpack_left_bottom_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%viscosity_left_bottom_rcv_buffer, chunks(chunk)%field%viscosity)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        CALL unpack_left_bottom_buffer_seq(chunk, depth, CELL_DATA, chunks(chunk)%soundspeed_left_bottom_rcv_buffer, chunks(chunk)%field%soundspeed)
    ENDIF
    IF(fields(FIELD_XVEL0).EQ.1) THEN
        CALL unpack_left_bottom_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%xvel0_left_bottom_rcv_buffer, chunks(chunk)%field%xvel0)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        CALL unpack_left_bottom_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%xvel1_left_bottom_rcv_buffer, chunks(chunk)%field%xvel1)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        CALL unpack_left_bottom_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%xvel0_left_bottom_rcv_buffer, chunks(chunk)%field%yvel0)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        CALL unpack_left_bottom_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%yvel1_left_bottom_rcv_buffer, chunks(chunk)%field%yvel1)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        CALL unpack_left_bottom_buffer_seq(chunk, depth, X_FACE_DATA, chunks(chunk)%volflux_x_left_bottom_rcv_buffer, chunks(chunk)%field%vol_flux_x)
    ENDIF
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        CALL unpack_left_bottom_buffer_seq(chunk, depth, Y_FACE_DATA, chunks(chunk)%volflux_y_left_bottom_rcv_buffer, chunks(chunk)%field%vol_flux_y)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        CALL unpack_left_bottom_buffer_seq(chunk, depth, X_FACE_DATA, chunks(chunk)%massflux_x_left_bottom_rcv_buffer, chunks(chunk)%field%mass_flux_x)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        CALL unpack_left_bottom_buffer_seq(chunk, depth, Y_FACE_DATA, chunks(chunk)%massflux_y_left_bottom_rcv_buffer, chunks(chunk)%field%mass_flux_y)
    ENDIF

END SUBROUTINE clover_exchange_unpack_all_buffers_left_bottom







!SUBROUTINE clover_exchange_message(chunk,field,                            &
!                                   left_snd_buffer,                        &
!                                   left_rcv_buffer,                        &
!                                   right_snd_buffer,                       &
!                                   right_rcv_buffer,                       &
!                                   bottom_snd_buffer,                      &
!                                   bottom_rcv_buffer,                      &
!                                   top_snd_buffer,                         &
!                                   top_rcv_buffer,                         &
!                                   left_top_snd_buffer,                    &
!                                   left_top_rcv_buffer,                    &
!                                   right_top_snd_buffer,                   &
!                                   right_top_rcv_buffer,                   &
!                                   right_bottom_snd_buffer,                &
!                                   right_bottom_rcv_buffer,                &
!                                   left_bottom_snd_buffer,                 &
!                                   left_bottom_rcv_buffer,                 &
!                                   depth,field_type)
!
!
!  IMPLICIT NONE
!
!  REAL(KIND=8) :: field(-1:,-1:) ! This seems to work for any type of mesh data
!  REAL(KIND=8) :: left_snd_buffer(:),left_rcv_buffer(:),right_snd_buffer(:),right_rcv_buffer(:)
!  REAL(KIND=8) :: bottom_snd_buffer(:),bottom_rcv_buffer(:),top_snd_buffer(:),top_rcv_buffer(:)
!
!  REAL(KIND=8) :: left_top_snd_buffer(:),left_top_rcv_buffer(:),right_top_snd_buffer(:),right_top_rcv_buffer(:)
!  REAL(KIND=8) :: right_bottom_snd_buffer(:),right_bottom_rcv_buffer(:),left_bottom_snd_buffer(:),left_bottom_rcv_buffer(:)
!
!  INTEGER      :: chunk,depth,field_type
!
!  INTEGER      :: size,err,tag,j,k,x_inc,y_inc,index
!  INTEGER      :: receiver,sender
!  INTEGER      :: left_neighbour_chunk, right_neighbour_chunk, bottom_neighbour_chunk, top_neighbour_chunk
!  INTEGER      :: left_top_neighbour_chunk, right_top_neighbour_chunk, right_bottom_neighbour_chunk, left_bottom_neighbour_chunk
!  INTEGER      :: leftedge, rightedge, bottomedge, topedge
!
!  ! Field type will either be cell, vertex, x_face or y_face to get the message limits correct
!
!  ! I am packing my own buffers. I am sure this could be improved with MPI data types
!  !  but this will do for now
!
!  ! I am also sending buffers to chunks with the same task id for now.
!  ! This can be improved in the future but at the moment there is just 1 chunk per task anyway
!
!  ! The tag will be a function of the sending chunk and the face it is coming from
!  !  like chunk 6 sending the left face
!
!  ! No open mp in here either. May be beneficial will packing and unpacking in the future, though I am not sure.
!
!  ! Change this so it will allow more than 1 chunk per task
!
!
!  ! Pack and send
!
!  ! These array modifications still need to be added on, plus the donor data location changes as in update_halo
!  IF(field_type.EQ.CELL_DATA) THEN
!    x_inc=0
!    y_inc=0
!  ENDIF
!  IF(field_type.EQ.VERTEX_DATA) THEN
!    x_inc=1
!    y_inc=1
!  ENDIF
!  IF(field_type.EQ.X_FACE_DATA) THEN
!    x_inc=1
!    y_inc=0
!  ENDIF
!  IF(field_type.EQ.Y_FACE_DATA) THEN
!    x_inc=0
!    y_inc=1
!  ENDIF
!
!    ! Pack real data into buffers
!    IF(parallel%task.EQ.chunks(chunk)%task) THEN
!
!        !pack all the data buffers required 
!        IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
!            CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, left_snd_buffer, field)
!        ENDIF
!        IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
!            CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, right_snd_buffer, field)
!        ENDIF
!        IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
!            CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, bottom_snd_buffer, field)
!        ENDIF
!        IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
!            CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, top_snd_buffer, field)
!        ENDIF
!        IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
!            CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, left_top_snd_buffer, field)
!        ENDIF
!        IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
!            CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, right_top_snd_buffer, field)
!        ENDIF
!        IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
!            CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, right_bottom_snd_buffer, field)
!        ENDIF
!        IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
!            CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, left_bottom_snd_buffer, field)
!        ENDIF
!
!        topedge = 0
!        bottomedge = 0
!        IF (chunks(chunk)%chunk_neighbours(chunk_top).EQ.external_face) THEN
!          topedge = depth
!        ENDIF
!        IF (chunks(chunk)%chunk_neighbours(chunk_bottom).EQ.external_face) THEN
!          bottomedge = depth
!        ENDIF
!
!        size=(1+(chunks(chunk)%field%y_max+y_inc+topedge)-(chunks(chunk)%field%y_min-bottomedge))*depth
!
!        !CALL pack_left_right_buffers(chunks(chunk)%field%x_min,chunks(chunk)%field%x_max, &
!        !                             chunks(chunk)%field%y_min,chunks(chunk)%field%y_max, &
!        !                             chunks(chunk)%chunk_neighbours(chunk_left),          &
!        !                             chunks(chunk)%chunk_neighbours(chunk_right),         &
!        !                             external_face,                                       &
!        !                             x_inc,y_inc,depth,size,                              &
!        !                             field,left_snd_buffer,right_snd_buffer)
!
!        IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
!            left_neighbour_chunk = chunks(chunk)%chunk_neighbours(chunk_left)
!            receiver=chunks(left_neighbour_chunk)%task
!            chunks(left_neighbour_chunk)[receiver+1]%right_rcv_buffer(1:size) = left_snd_buffer(1:size)
!        ENDIF
!
!        IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
!            right_neighbour_chunk = chunks(chunk)%chunk_neighbours(chunk_right)
!            receiver = chunks(right_neighbour_chunk)%task
!            chunks(right_neighbour_chunk)[receiver+1]%left_rcv_buffer(1:size) = right_snd_buffer(1:size)
!        ENDIF
!
!
!        leftedge= 0
!        rightedge= 0
!        IF (chunks(chunk)%chunk_neighbours(chunk_left).EQ.external_face) THEN
!          leftedge = depth
!        ENDIF
!        IF (chunks(chunk)%chunk_neighbours(chunk_right).EQ.external_face) THEN
!          rightedge = depth
!        ENDIF
!
!        size=(1+(chunks(chunk)%field%x_max+x_inc+rightedge)-(chunks(chunk)%field%x_min-leftedge))*depth
!
!        !CALL pack_top_bottom_buffers(chunks(chunk)%field%x_min,chunks(chunk)%field%x_max, &
!        !                             chunks(chunk)%field%y_min,chunks(chunk)%field%y_max, &
!        !                             chunks(chunk)%chunk_neighbours(chunk_bottom),        &
!        !                             chunks(chunk)%chunk_neighbours(chunk_top),           &
!        !                             external_face,                                       &
!        !                             x_inc,y_inc,depth,size,                              &
!        !                             field,bottom_snd_buffer,top_snd_buffer)
!
!        IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
!            bottom_neighbour_chunk = chunks(chunk)%chunk_neighbours(chunk_bottom)
!            receiver=chunks(bottom_neighbour_chunk)%task
!            chunks(bottom_neighbour_chunk)[receiver+1]%top_rcv_buffer(1:size) = bottom_snd_buffer(1:size)
!        ENDIF
!
!        IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
!            top_neighbour_chunk = chunks(chunk)%chunk_neighbours(chunk_top)
!            receiver=chunks(top_neighbour_chunk)%task
!            chunks(top_neighbour_chunk)[receiver+1]%bottom_rcv_buffer(1:size) = top_snd_buffer(1:size)
!        ENDIF
!
!
!
!        size = depth*depth
!
!        IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
!            left_top_neighbour_chunk = chunks(chunk)%chunk_neighbours(chunk_left_top)
!            receiver = chunks(left_top_neighbour_chunk)%task + 1
!            chunks(left_top_neighbour_chunk)[receiver]%right_bottom_rcv_buffer(1:size) = left_top_snd_buffer(1:size)
!        ENDIF
!
!        IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
!            right_top_neighbour_chunk = chunks(chunk)%chunk_neighbours(chunk_right_top)
!            receiver = chunks(right_top_neighbour_chunk)%task + 1
!            chunks(right_top_neighbour_chunk)[receiver]%left_bottom_rcv_buffer(1:size) = right_top_snd_buffer(1:size)
!        ENDIF
!
!
!        IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
!            right_bottom_neighbour_chunk = chunks(chunk)%chunk_neighbours(chunk_right_bottom)
!            receiver = chunks(right_bottom_neighbour_chunk)%task + 1
!            chunks(right_bottom_neighbour_chunk)[receiver]%left_top_rcv_buffer(1:size) = right_bottom_snd_buffer(1:size)
!        ENDIF
!
!        IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
!            left_bottom_neighbour_chunk = chunks(chunk)%chunk_neighbours(chunk_left_bottom)
!            receiver = chunks(left_bottom_neighbour_chunk)%task + 1
!            chunks(left_bottom_neighbour_chunk)[receiver]%right_top_rcv_buffer(1:size) = left_bottom_snd_buffer(1:size)
!        ENDIF
!
!
!  ! Wait for the messages
!#ifdef LOCAL_SYNC
!        sync images( chunks(chunk)%imageNeighbours )
!#else
!        sync all
!#endif
!
!
!        ! Unpack buffers in halo cells
!        IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
!            CALL unpack_left_buffer_seq(chunk, depth, x_inc, y_inc, left_rcv_buffer, field)
!        ENDIF
!        IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
!            CALL unpack_right_buffer_seq(chunk, depth, x_inc, y_inc, right_rcv_buffer, field)
!        ENDIF
!        IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
!            CALL unpack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, bottom_rcv_buffer, field)
!        ENDIF
!        IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
!            CALL unpack_top_buffer_seq(chunk, depth, x_inc, y_inc, top_rcv_buffer, field)
!        ENDIF
!        IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
!            CALL unpack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, left_top_rcv_buffer, field)
!        ENDIF
!        IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) ) THEN
!            CALL unpack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, right_top_rcv_buffer, field)
!        ENDIF
!        IF ( (chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
!            CALL unpack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, right_bottom_rcv_buffer, field)
!        ENDIF
!        IF ( (chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) .AND. (chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) ) THEN
!            CALL unpack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, left_bottom_rcv_buffer, field)
!        ENDIF
!
!
!        ! Wait for the messages
!#ifdef LOCAL_SYNC
!        sync images( chunks(chunk)%imageNeighbours )
!#else
!        sync all
!#endif
!        ! Unpack buffers in halo cells
!        !IF(parallel%task.EQ.chunks(chunk)%task) THEN
!        !    CALL unpack_left_right_buffers(chunks(chunk)%field%x_min,chunks(chunk)%field%x_max, &
!        !                                   chunks(chunk)%field%y_min,chunks(chunk)%field%y_max, &
!        !                                   chunks(chunk)%chunk_neighbours(chunk_left),          &
!        !                                   chunks(chunk)%chunk_neighbours(chunk_right),         &
!        !                                   external_face,                                       &
!        !                                   x_inc,y_inc,depth,size,                              &
!        !                                   field,left_rcv_buffer,right_rcv_buffer)
!
!        !  CALL unpack_top_bottom_buffers(chunks(chunk)%field%x_min,chunks(chunk)%field%x_max, &
!        !                                 chunks(chunk)%field%y_min,chunks(chunk)%field%y_max, &
!        !                                 chunks(chunk)%chunk_neighbours(chunk_bottom),        &
!        !                                 chunks(chunk)%chunk_neighbours(chunk_top),           &
!        !                                 external_face,                                       &
!        !                                 x_inc,y_inc,depth,size,                              &
!        !                                 field,bottom_rcv_buffer,top_rcv_buffer)
!
!    ENDIF
!
!END SUBROUTINE clover_exchange_message

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

SUBROUTINE clover_check_error(error)

  IMPLICIT NONE

  INTEGER :: error

  INTEGER :: maximum

  maximum=error

  CALL CO_MAX(error, maximum)

  error=maximum

END SUBROUTINE clover_check_error


END MODULE clover_module

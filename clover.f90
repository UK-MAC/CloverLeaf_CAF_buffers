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

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_neighbour_chunk)[receiver]%density0_right_rcv_buffer(1:size) = chunks(chunk)%density0_left_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        !CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
        !                                        chunks(chunk)%density1_left_snd_buffer, chunks(chunk)%density1_right_rcv_buffer, chunks(chunk)%field%density1)

        CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%density1_left_snd_buffer, chunks(chunk)%field%density1)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_neighbour_chunk)[receiver]%density1_right_rcv_buffer(1:size) = chunks(chunk)%density1_left_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        !CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
        !                                        chunks(chunk)%energy0_left_snd_buffer, chunks(chunk)%energy0_right_rcv_buffer, chunks(chunk)%field%energy0)
        CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%energy0_left_snd_buffer, chunks(chunk)%field%energy0)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_neighbour_chunk)[receiver]%energy0_right_rcv_buffer(1:size) = chunks(chunk)%energy0_left_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        !CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
        !                                        chunks(chunk)%energy1_left_snd_buffer, chunks(chunk)%energy1_right_rcv_buffer, chunks(chunk)%field%energy1)
        CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%energy1_left_snd_buffer, chunks(chunk)%field%energy1)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_neighbour_chunk)[receiver]%energy1_right_rcv_buffer(1:size) = chunks(chunk)%energy1_left_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        !CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
        !                                        chunks(chunk)%pressure_left_snd_buffer, chunks(chunk)%pressure_right_rcv_buffer, chunks(chunk)%field%pressure)
        CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%pressure_left_snd_buffer, chunks(chunk)%field%pressure)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_neighbour_chunk)[receiver]%pressure_right_rcv_buffer(1:size) = chunks(chunk)%pressure_left_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        !CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
        !                                        chunks(chunk)%viscosity_left_snd_buffer, chunks(chunk)%viscosity_right_rcv_buffer, chunks(chunk)%field%viscosity)
        CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%viscosity_left_snd_buffer, chunks(chunk)%field%viscosity)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_neighbour_chunk)[receiver]%viscosity_right_rcv_buffer(1:size) = chunks(chunk)%viscosity_left_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        !CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
        !                                        chunks(chunk)%soundspeed_left_snd_buffer, chunks(chunk)%soundspeed_right_rcv_buffer, chunks(chunk)%field%soundspeed)
        CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%soundspeed_left_snd_buffer, chunks(chunk)%field%soundspeed)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_neighbour_chunk)[receiver]%soundspeed_right_rcv_buffer(1:size) = chunks(chunk)%soundspeed_left_snd_buffer(1:size)
    ENDIF

    x_inc=1
    y_inc=1
    size=(1+(chunks(chunk)%field%y_max+y_inc+topedge)-(chunks(chunk)%field%y_min-bottomedge))*depth

    IF(fields(FIELD_XVEL0).EQ.1) THEN
        !CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, VERTEX_DATA, & 
        !                                        chunks(chunk)%xvel0_left_snd_buffer, chunks(chunk)%xvel0_right_rcv_buffer, chunks(chunk)%field%xvel0)
        CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%xvel0_left_snd_buffer, chunks(chunk)%field%xvel0)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_neighbour_chunk)[receiver]%xvel0_right_rcv_buffer(1:size) = chunks(chunk)%xvel0_left_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        !CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, VERTEX_DATA, & 
        !                                        chunks(chunk)%xvel1_left_snd_buffer, chunks(chunk)%xvel1_right_rcv_buffer, chunks(chunk)%field%xvel1)
        CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%xvel1_left_snd_buffer, chunks(chunk)%field%xvel1)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_neighbour_chunk)[receiver]%xvel1_right_rcv_buffer(1:size) = chunks(chunk)%xvel1_left_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        !CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, VERTEX_DATA, & 
        !                                        chunks(chunk)%yvel0_left_snd_buffer, chunks(chunk)%yvel0_right_rcv_buffer, chunks(chunk)%field%yvel0)
        CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%yvel0_left_snd_buffer, chunks(chunk)%field%yvel0)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_neighbour_chunk)[receiver]%yvel0_right_rcv_buffer(1:size) = chunks(chunk)%yvel0_left_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        !CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, VERTEX_DATA, & 
        !                                        chunks(chunk)%yvel1_left_snd_buffer, chunks(chunk)%yvel1_right_rcv_buffer, chunks(chunk)%field%yvel1)
        CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%yvel1_left_snd_buffer, chunks(chunk)%field%yvel1)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_neighbour_chunk)[receiver]%yvel1_right_rcv_buffer(1:size) = chunks(chunk)%yvel1_left_snd_buffer(1:size)
    ENDIF

    x_inc=1
    y_inc=0
    size=(1+(chunks(chunk)%field%y_max+y_inc+topedge)-(chunks(chunk)%field%y_min-bottomedge))*depth

    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        !CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, X_FACE_DATA, & 
        !                                        chunks(chunk)%volflux_x_left_snd_buffer, chunks(chunk)%volflux_x_right_rcv_buffer, chunks(chunk)%field%vol_flux_x)
        CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%volflux_x_left_snd_buffer, chunks(chunk)%field%vol_flux_x)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_neighbour_chunk)[receiver]%volflux_x_right_rcv_buffer(1:size) = chunks(chunk)%volflux_x_left_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        !CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, X_FACE_DATA, & 
        !                                        chunks(chunk)%massflux_x_left_snd_buffer, chunks(chunk)%massflux_x_right_rcv_buffer, chunks(chunk)%field%mass_flux_x)
        CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%massflux_x_left_snd_buffer, chunks(chunk)%field%mass_flux_x)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_neighbour_chunk)[receiver]%massflux_x_right_rcv_buffer(1:size) = chunks(chunk)%massflux_x_left_snd_buffer(1:size)
    ENDIF

    x_inc=0
    y_inc=1
    size=(1+(chunks(chunk)%field%y_max+y_inc+topedge)-(chunks(chunk)%field%y_min-bottomedge))*depth

    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        !CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, Y_FACE_DATA, & 
        !                                        chunks(chunk)%volflux_y_left_snd_buffer, chunks(chunk)%volflux_y_right_rcv_buffer, chunks(chunk)%field%vol_flux_y)
        CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%volflux_y_left_snd_buffer, chunks(chunk)%field%vol_flux_y)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_neighbour_chunk)[receiver]%volflux_y_right_rcv_buffer(1:size) = chunks(chunk)%volflux_y_left_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        !CALL clover_exchange_write_message_left(chunk, depth, receiver, topedge, bottomedge, Y_FACE_DATA, & 
        !                                        chunks(chunk)%massflux_y_left_snd_buffer, chunks(chunk)%massflux_y_right_rcv_buffer, chunks(chunk)%field%mass_flux_y)
        CALL pack_left_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%massflux_y_left_snd_buffer, chunks(chunk)%field%mass_flux_y)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
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

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_neighbour_chunk)[receiver]%density0_left_rcv_buffer(1:size) = chunks(chunk)%density0_right_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        !CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
        !                                         chunks(chunk)%density1_right_snd_buffer, chunks(chunk)%density1_left_rcv_buffer, chunks(chunk)%field%density1)
        CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%density1_right_snd_buffer, chunks(chunk)%field%density1)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_neighbour_chunk)[receiver]%density1_left_rcv_buffer(1:size) = chunks(chunk)%density1_right_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        !CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
        !                                         chunks(chunk)%energy0_right_snd_buffer, chunks(chunk)%energy0_left_rcv_buffer, chunks(chunk)%field%energy0)
        CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%energy0_right_snd_buffer, chunks(chunk)%field%energy0)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_neighbour_chunk)[receiver]%energy0_left_rcv_buffer(1:size) = chunks(chunk)%energy0_right_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        !CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
        !                                         chunks(chunk)%energy1_right_snd_buffer, chunks(chunk)%energy1_left_rcv_buffer, chunks(chunk)%field%energy1)
        CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%energy1_right_snd_buffer, chunks(chunk)%field%energy1)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_neighbour_chunk)[receiver]%energy1_left_rcv_buffer(1:size) = chunks(chunk)%energy1_right_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        !CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
        !                                         chunks(chunk)%pressure_right_snd_buffer, chunks(chunk)%pressure_left_rcv_buffer, chunks(chunk)%field%pressure)
        CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%pressure_right_snd_buffer, chunks(chunk)%field%pressure)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_neighbour_chunk)[receiver]%pressure_left_rcv_buffer(1:size) = chunks(chunk)%pressure_right_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        !CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
        !                                         chunks(chunk)%viscosity_right_snd_buffer, chunks(chunk)%viscosity_left_rcv_buffer, chunks(chunk)%field%viscosity)
        CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%viscosity_right_snd_buffer, chunks(chunk)%field%viscosity)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_neighbour_chunk)[receiver]%viscosity_left_rcv_buffer(1:size) = chunks(chunk)%viscosity_right_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        !CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, CELL_DATA, & 
        !                                         chunks(chunk)%soundspeed_right_snd_buffer, chunks(chunk)%soundspeed_left_rcv_buffer, chunks(chunk)%field%soundspeed)
        CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%soundspeed_right_snd_buffer, chunks(chunk)%field%soundspeed)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_neighbour_chunk)[receiver]%soundspeed_left_rcv_buffer(1:size) = chunks(chunk)%soundspeed_right_snd_buffer(1:size)
    ENDIF

    x_inc=1
    y_inc=1
    size=(1+(chunks(chunk)%field%y_max+y_inc+topedge)-(chunks(chunk)%field%y_min-bottomedge))*depth

    IF(fields(FIELD_XVEL0).EQ.1) THEN
        !CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, VERTEX_DATA, & 
        !                                         chunks(chunk)%xvel0_right_snd_buffer, chunks(chunk)%xvel0_left_rcv_buffer, chunks(chunk)%field%xvel0)
        CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%xvel0_right_snd_buffer, chunks(chunk)%field%xvel0)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_neighbour_chunk)[receiver]%xvel0_left_rcv_buffer(1:size) = chunks(chunk)%xvel0_right_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        !CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, VERTEX_DATA, & 
        !                                         chunks(chunk)%xvel1_right_snd_buffer, chunks(chunk)%xvel1_left_rcv_buffer, chunks(chunk)%field%xvel1)
        CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%xvel1_right_snd_buffer, chunks(chunk)%field%xvel1)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_neighbour_chunk)[receiver]%xvel1_left_rcv_buffer(1:size) = chunks(chunk)%xvel1_right_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        !CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, VERTEX_DATA, & 
        !                                         chunks(chunk)%yvel0_right_snd_buffer, chunks(chunk)%yvel0_left_rcv_buffer, chunks(chunk)%field%yvel0)
        CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%yvel0_right_snd_buffer, chunks(chunk)%field%yvel0)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_neighbour_chunk)[receiver]%yvel0_left_rcv_buffer(1:size) = chunks(chunk)%yvel0_right_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        !CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, VERTEX_DATA, & 
        !                                         chunks(chunk)%yvel1_right_snd_buffer, chunks(chunk)%yvel1_left_rcv_buffer, chunks(chunk)%field%yvel1)
        CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%yvel1_right_snd_buffer, chunks(chunk)%field%yvel1)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_neighbour_chunk)[receiver]%yvel1_left_rcv_buffer(1:size) = chunks(chunk)%yvel1_right_snd_buffer(1:size)
    ENDIF

    x_inc=1
    y_inc=0
    size=(1+(chunks(chunk)%field%y_max+y_inc+topedge)-(chunks(chunk)%field%y_min-bottomedge))*depth

    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        !CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, X_FACE_DATA, & 
        !                                         chunks(chunk)%volflux_x_right_snd_buffer, chunks(chunk)%volflux_x_left_rcv_buffer, chunks(chunk)%field%vol_flux_x)
        CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%volflux_x_right_snd_buffer, chunks(chunk)%field%vol_flux_x)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_neighbour_chunk)[receiver]%volflux_x_left_rcv_buffer(1:size) = chunks(chunk)%volflux_x_right_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        !CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, X_FACE_DATA, & 
        !                                         chunks(chunk)%massflux_x_right_snd_buffer, chunks(chunk)%massflux_x_left_rcv_buffer, chunks(chunk)%field%mass_flux_x)
        CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%massflux_x_right_snd_buffer, chunks(chunk)%field%mass_flux_x)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_neighbour_chunk)[receiver]%massflux_x_left_rcv_buffer(1:size) = chunks(chunk)%massflux_x_right_snd_buffer(1:size)
    ENDIF

    x_inc=0
    y_inc=1
    size=(1+(chunks(chunk)%field%y_max+y_inc+topedge)-(chunks(chunk)%field%y_min-bottomedge))*depth
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        !CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, Y_FACE_DATA, & 
        !                                         chunks(chunk)%volflux_y_right_snd_buffer, chunks(chunk)%volflux_y_left_rcv_buffer, chunks(chunk)%field%vol_flux_y)
        CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%volflux_y_right_snd_buffer, chunks(chunk)%field%vol_flux_y)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_neighbour_chunk)[receiver]%volflux_y_left_rcv_buffer(1:size) = chunks(chunk)%volflux_y_right_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        !CALL clover_exchange_write_message_right(chunk, depth, receiver, topedge, bottomedge, Y_FACE_DATA, & 
        !                                         chunks(chunk)%massflux_y_right_snd_buffer, chunks(chunk)%massflux_y_left_rcv_buffer, chunks(chunk)%field%mass_flux_y)
        CALL pack_right_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%massflux_y_right_snd_buffer, chunks(chunk)%field%mass_flux_y)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
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

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(bottom_neighbour_chunk)[receiver]%density0_top_rcv_buffer(1:size) = chunks(chunk)%density0_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        !CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
        !                                          chunks(chunk)%density1_bottom_snd_buffer, chunks(chunk)%density1_top_rcv_buffer, chunks(chunk)%field%density1)
        CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%density1_bottom_snd_buffer, chunks(chunk)%field%density1)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(bottom_neighbour_chunk)[receiver]%density1_top_rcv_buffer(1:size) = chunks(chunk)%density1_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        !CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
        !                                          chunks(chunk)%energy0_bottom_snd_buffer, chunks(chunk)%energy0_top_rcv_buffer, chunks(chunk)%field%energy0)
        CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%energy0_bottom_snd_buffer, chunks(chunk)%field%energy0)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(bottom_neighbour_chunk)[receiver]%energy0_top_rcv_buffer(1:size) = chunks(chunk)%energy0_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        !CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
        !                                          chunks(chunk)%energy1_bottom_snd_buffer, chunks(chunk)%energy1_top_rcv_buffer, chunks(chunk)%field%energy1)
        CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%energy1_bottom_snd_buffer, chunks(chunk)%field%energy1)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(bottom_neighbour_chunk)[receiver]%energy1_top_rcv_buffer(1:size) = chunks(chunk)%energy1_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        !CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
        !                                          chunks(chunk)%pressure_bottom_snd_buffer, chunks(chunk)%pressure_top_rcv_buffer, chunks(chunk)%field%pressure)
        CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%pressure_bottom_snd_buffer, chunks(chunk)%field%pressure)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(bottom_neighbour_chunk)[receiver]%pressure_top_rcv_buffer(1:size) = chunks(chunk)%pressure_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        !CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
        !                                          chunks(chunk)%viscosity_bottom_snd_buffer, chunks(chunk)%viscosity_top_rcv_buffer, chunks(chunk)%field%viscosity)
        CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%viscosity_bottom_snd_buffer, chunks(chunk)%field%viscosity)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(bottom_neighbour_chunk)[receiver]%viscosity_top_rcv_buffer(1:size) = chunks(chunk)%viscosity_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        !CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
        !                                          chunks(chunk)%soundspeed_bottom_snd_buffer, chunks(chunk)%soundspeed_top_rcv_buffer, chunks(chunk)%field%soundspeed)
        CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%soundspeed_bottom_snd_buffer, chunks(chunk)%field%soundspeed)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(bottom_neighbour_chunk)[receiver]%soundspeed_top_rcv_buffer(1:size) = chunks(chunk)%soundspeed_bottom_snd_buffer(1:size)
    ENDIF

    x_inc=1
    y_inc=1
    size=(1+(chunks(chunk)%field%x_max+x_inc+rightedge)-(chunks(chunk)%field%x_min-leftedge))*depth

    IF(fields(FIELD_XVEL0).EQ.1) THEN
        !CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, VERTEX_DATA, & 
        !                                          chunks(chunk)%xvel0_bottom_snd_buffer, chunks(chunk)%xvel0_top_rcv_buffer, chunks(chunk)%field%xvel0)
        CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%xvel0_bottom_snd_buffer, chunks(chunk)%field%xvel0)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(bottom_neighbour_chunk)[receiver]%xvel0_top_rcv_buffer(1:size) = chunks(chunk)%xvel0_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        !CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, VERTEX_DATA, & 
        !                                          chunks(chunk)%xvel1_bottom_snd_buffer, chunks(chunk)%xvel1_top_rcv_buffer, chunks(chunk)%field%xvel1)
        CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%xvel1_bottom_snd_buffer, chunks(chunk)%field%xvel1)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(bottom_neighbour_chunk)[receiver]%xvel1_top_rcv_buffer(1:size) = chunks(chunk)%xvel1_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        !CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, VERTEX_DATA, & 
        !                                          chunks(chunk)%yvel0_bottom_snd_buffer, chunks(chunk)%yvel0_top_rcv_buffer, chunks(chunk)%field%yvel0)
        CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%yvel0_bottom_snd_buffer, chunks(chunk)%field%yvel0)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(bottom_neighbour_chunk)[receiver]%yvel0_top_rcv_buffer(1:size) = chunks(chunk)%yvel0_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        !CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, VERTEX_DATA, & 
        !                                          chunks(chunk)%yvel1_bottom_snd_buffer, chunks(chunk)%yvel1_top_rcv_buffer, chunks(chunk)%field%yvel1)
        CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%yvel1_bottom_snd_buffer, chunks(chunk)%field%yvel1)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(bottom_neighbour_chunk)[receiver]%yvel1_top_rcv_buffer(1:size) = chunks(chunk)%yvel1_bottom_snd_buffer(1:size)
    ENDIF

    x_inc=1
    y_inc=0
    size=(1+(chunks(chunk)%field%x_max+x_inc+rightedge)-(chunks(chunk)%field%x_min-leftedge))*depth

    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        !CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, X_FACE_DATA, & 
        !                                          chunks(chunk)%volflux_x_bottom_snd_buffer, chunks(chunk)%volflux_x_top_rcv_buffer, chunks(chunk)%field%vol_flux_x)
        CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%volflux_x_bottom_snd_buffer, chunks(chunk)%field%vol_flux_x)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(bottom_neighbour_chunk)[receiver]%volflux_x_top_rcv_buffer(1:size) = chunks(chunk)%volflux_x_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        !CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, X_FACE_DATA, & 
        !                                          chunks(chunk)%massflux_x_bottom_snd_buffer, chunks(chunk)%massflux_x_top_rcv_buffer, chunks(chunk)%field%mass_flux_x)
        CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%massflux_x_bottom_snd_buffer, chunks(chunk)%field%mass_flux_x)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(bottom_neighbour_chunk)[receiver]%massflux_x_top_rcv_buffer(1:size) = chunks(chunk)%massflux_x_bottom_snd_buffer(1:size)
    ENDIF

    x_inc=0
    y_inc=1
    size=(1+(chunks(chunk)%field%x_max+x_inc+rightedge)-(chunks(chunk)%field%x_min-leftedge))*depth

    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        !CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, Y_FACE_DATA, & 
        !                                          chunks(chunk)%volflux_y_bottom_snd_buffer, chunks(chunk)%volflux_y_top_rcv_buffer, chunks(chunk)%field%vol_flux_y)
        CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%volflux_y_bottom_snd_buffer, chunks(chunk)%field%vol_flux_y)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(bottom_neighbour_chunk)[receiver]%volflux_y_top_rcv_buffer(1:size) = chunks(chunk)%volflux_y_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        !CALL clover_exchange_write_message_bottom(chunk, depth, receiver, leftedge, rightedge, Y_FACE_DATA, & 
        !                                          chunks(chunk)%massflux_y_bottom_snd_buffer, chunks(chunk)%massflux_y_top_rcv_buffer, chunks(chunk)%field%mass_flux_y)
        CALL pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%massflux_y_bottom_snd_buffer, chunks(chunk)%field%mass_flux_y)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
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

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(top_neighbour_chunk)[receiver]%density0_bottom_rcv_buffer(1:size) = chunks(chunk)%density0_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        !CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
        !                                       chunks(chunk)%density1_top_snd_buffer, chunks(chunk)%density1_bottom_rcv_buffer, chunks(chunk)%field%density1)
        CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%density1_top_snd_buffer, chunks(chunk)%field%density1)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(top_neighbour_chunk)[receiver]%density1_bottom_rcv_buffer(1:size) = chunks(chunk)%density1_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        !CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
        !                                       chunks(chunk)%energy0_top_snd_buffer, chunks(chunk)%energy0_bottom_rcv_buffer, chunks(chunk)%field%energy0)
        CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%energy0_top_snd_buffer, chunks(chunk)%field%energy0)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(top_neighbour_chunk)[receiver]%energy0_bottom_rcv_buffer(1:size) = chunks(chunk)%energy0_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        !CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
        !                                       chunks(chunk)%energy1_top_snd_buffer, chunks(chunk)%energy1_bottom_rcv_buffer, chunks(chunk)%field%energy1)
        CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%energy1_top_snd_buffer, chunks(chunk)%field%energy1)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(top_neighbour_chunk)[receiver]%energy1_bottom_rcv_buffer(1:size) = chunks(chunk)%energy1_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        !CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
        !                                       chunks(chunk)%pressure_top_snd_buffer, chunks(chunk)%pressure_bottom_rcv_buffer, chunks(chunk)%field%pressure)
        CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%pressure_top_snd_buffer, chunks(chunk)%field%pressure)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(top_neighbour_chunk)[receiver]%pressure_bottom_rcv_buffer(1:size) = chunks(chunk)%pressure_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        !CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
        !                                       chunks(chunk)%viscosity_top_snd_buffer, chunks(chunk)%viscosity_bottom_rcv_buffer, chunks(chunk)%field%viscosity)
        CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%viscosity_top_snd_buffer, chunks(chunk)%field%viscosity)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(top_neighbour_chunk)[receiver]%viscosity_bottom_rcv_buffer(1:size) = chunks(chunk)%viscosity_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        !CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, CELL_DATA, & 
        !                                       chunks(chunk)%soundspeed_top_snd_buffer, chunks(chunk)%soundspeed_bottom_rcv_buffer, chunks(chunk)%field%soundspeed)
        CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%soundspeed_top_snd_buffer, chunks(chunk)%field%soundspeed)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(top_neighbour_chunk)[receiver]%soundspeed_bottom_rcv_buffer(1:size) = chunks(chunk)%soundspeed_top_snd_buffer(1:size)
    ENDIF

    x_inc=1
    y_inc=1
    size=(1+(chunks(chunk)%field%x_max+x_inc+rightedge)-(chunks(chunk)%field%x_min-leftedge))*depth

    IF(fields(FIELD_XVEL0).EQ.1) THEN
        !CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, VERTEX_DATA, & 
        !                                       chunks(chunk)%xvel0_top_snd_buffer, chunks(chunk)%xvel0_bottom_rcv_buffer, chunks(chunk)%field%xvel0)
        CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%xvel0_top_snd_buffer, chunks(chunk)%field%xvel0)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(top_neighbour_chunk)[receiver]%xvel0_bottom_rcv_buffer(1:size) = chunks(chunk)%xvel0_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        !CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, VERTEX_DATA, & 
        !                                       chunks(chunk)%xvel1_top_snd_buffer, chunks(chunk)%xvel1_bottom_rcv_buffer, chunks(chunk)%field%xvel1)
        CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%xvel1_top_snd_buffer, chunks(chunk)%field%xvel1)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(top_neighbour_chunk)[receiver]%xvel1_bottom_rcv_buffer(1:size) = chunks(chunk)%xvel1_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        !CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, VERTEX_DATA, & 
        !                                       chunks(chunk)%yvel0_top_snd_buffer, chunks(chunk)%yvel0_bottom_rcv_buffer, chunks(chunk)%field%yvel0)
        CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%yvel0_top_snd_buffer, chunks(chunk)%field%yvel0)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(top_neighbour_chunk)[receiver]%yvel0_bottom_rcv_buffer(1:size) = chunks(chunk)%yvel0_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        !CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, VERTEX_DATA, & 
        !                                       chunks(chunk)%yvel1_top_snd_buffer, chunks(chunk)%yvel1_bottom_rcv_buffer, chunks(chunk)%field%yvel1)
        CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%yvel1_top_snd_buffer, chunks(chunk)%field%yvel1)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(top_neighbour_chunk)[receiver]%yvel1_bottom_rcv_buffer(1:size) = chunks(chunk)%yvel1_top_snd_buffer(1:size)
    ENDIF

    x_inc=1
    y_inc=0
    size=(1+(chunks(chunk)%field%x_max+x_inc+rightedge)-(chunks(chunk)%field%x_min-leftedge))*depth
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        !CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, X_FACE_DATA, & 
        !                                       chunks(chunk)%volflux_x_top_snd_buffer, chunks(chunk)%volflux_x_bottom_rcv_buffer, chunks(chunk)%field%vol_flux_x)
        CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%volflux_x_top_snd_buffer, chunks(chunk)%field%vol_flux_x)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(top_neighbour_chunk)[receiver]%volflux_x_bottom_rcv_buffer(1:size) = chunks(chunk)%volflux_x_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        !CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, X_FACE_DATA, & 
        !                                       chunks(chunk)%massflux_x_top_snd_buffer, chunks(chunk)%massflux_x_bottom_rcv_buffer, chunks(chunk)%field%mass_flux_x)
        CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%massflux_x_top_snd_buffer, chunks(chunk)%field%mass_flux_x)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(top_neighbour_chunk)[receiver]%massflux_x_bottom_rcv_buffer(1:size) = chunks(chunk)%massflux_x_top_snd_buffer(1:size)
    ENDIF

    x_inc=0
    y_inc=1
    size=(1+(chunks(chunk)%field%x_max+x_inc+rightedge)-(chunks(chunk)%field%x_min-leftedge))*depth
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        !CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, Y_FACE_DATA, & 
        !                                       chunks(chunk)%volflux_y_top_snd_buffer, chunks(chunk)%volflux_y_bottom_rcv_buffer, chunks(chunk)%field%vol_flux_y)
        CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%volflux_y_top_snd_buffer, chunks(chunk)%field%vol_flux_y)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(top_neighbour_chunk)[receiver]%volflux_y_bottom_rcv_buffer(1:size) = chunks(chunk)%volflux_y_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        !CALL clover_exchange_write_message_top(chunk, depth, receiver, leftedge, rightedge, Y_FACE_DATA, & 
        !                                       chunks(chunk)%massflux_y_top_snd_buffer, chunks(chunk)%massflux_y_bottom_rcv_buffer, chunks(chunk)%field%mass_flux_y)
        CALL pack_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%massflux_y_top_snd_buffer, chunks(chunk)%field%mass_flux_y)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
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

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_top_neighbour_chunk)[receiver]%density0_right_bottom_rcv_buffer(1:size) = chunks(chunk)%density0_left_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        !CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, CELL_DATA, & 
        !                                            chunks(chunk)%density1_left_top_snd_buffer, chunks(chunk)%density1_right_bottom_rcv_buffer, chunks(chunk)%field%density1)
        CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%density1_left_top_snd_buffer, chunks(chunk)%field%density1)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_top_neighbour_chunk)[receiver]%density1_right_bottom_rcv_buffer(1:size) = chunks(chunk)%density1_left_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        !CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, CELL_DATA, & 
        !                                            chunks(chunk)%energy0_left_top_snd_buffer, chunks(chunk)%energy0_right_bottom_rcv_buffer, chunks(chunk)%field%energy0)
        CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%energy0_left_top_snd_buffer, chunks(chunk)%field%energy0)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_top_neighbour_chunk)[receiver]%energy0_right_bottom_rcv_buffer(1:size) = chunks(chunk)%energy0_left_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        !CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, CELL_DATA, & 
        !                                            chunks(chunk)%energy1_left_top_snd_buffer, chunks(chunk)%energy1_right_bottom_rcv_buffer, chunks(chunk)%field%energy1)
        CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%energy1_left_top_snd_buffer, chunks(chunk)%field%energy1)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_top_neighbour_chunk)[receiver]%energy1_right_bottom_rcv_buffer(1:size) = chunks(chunk)%energy1_left_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        !CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, CELL_DATA, & 
        !                                            chunks(chunk)%pressure_left_top_snd_buffer, chunks(chunk)%pressure_right_bottom_rcv_buffer, chunks(chunk)%field%pressure)
        CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%pressure_left_top_snd_buffer, chunks(chunk)%field%pressure)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_top_neighbour_chunk)[receiver]%pressure_right_bottom_rcv_buffer(1:size) = chunks(chunk)%pressure_left_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        !CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, CELL_DATA, & 
        !                                            chunks(chunk)%viscosity_left_top_snd_buffer, chunks(chunk)%viscosity_right_bottom_rcv_buffer, chunks(chunk)%field%viscosity)
        CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%viscosity_left_top_snd_buffer, chunks(chunk)%field%viscosity)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_top_neighbour_chunk)[receiver]%viscosity_right_bottom_rcv_buffer(1:size) = chunks(chunk)%viscosity_left_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        !CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, CELL_DATA, & 
        !                                            chunks(chunk)%soundspeed_left_top_snd_buffer, chunks(chunk)%soundspeed_right_bottom_rcv_buffer, chunks(chunk)%field%soundspeed)
        CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%soundspeed_left_top_snd_buffer, chunks(chunk)%field%soundspeed)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_top_neighbour_chunk)[receiver]%soundspeed_right_bottom_rcv_buffer(1:size) = chunks(chunk)%soundspeed_left_top_snd_buffer(1:size)
    ENDIF

    x_inc=1
    y_inc=1
    IF(fields(FIELD_XVEL0).EQ.1) THEN
        !CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, VERTEX_DATA, & 
        !                                            chunks(chunk)%xvel0_left_top_snd_buffer, chunks(chunk)%xvel0_right_bottom_rcv_buffer, chunks(chunk)%field%xvel0)
        CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%xvel0_left_top_snd_buffer, chunks(chunk)%field%xvel0)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_top_neighbour_chunk)[receiver]%xvel0_right_bottom_rcv_buffer(1:size) = chunks(chunk)%xvel0_left_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        !CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, VERTEX_DATA, & 
        !                                            chunks(chunk)%xvel1_left_top_snd_buffer, chunks(chunk)%xvel1_right_bottom_rcv_buffer, chunks(chunk)%field%xvel1)
        CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%xvel1_left_top_snd_buffer, chunks(chunk)%field%xvel1)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_top_neighbour_chunk)[receiver]%xvel1_right_bottom_rcv_buffer(1:size) = chunks(chunk)%xvel1_left_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        !CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, VERTEX_DATA, & 
        !                                            chunks(chunk)%yvel0_left_top_snd_buffer, chunks(chunk)%yvel0_right_bottom_rcv_buffer, chunks(chunk)%field%yvel0)
        CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%yvel0_left_top_snd_buffer, chunks(chunk)%field%yvel0)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_top_neighbour_chunk)[receiver]%yvel0_right_bottom_rcv_buffer(1:size) = chunks(chunk)%yvel0_left_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        !CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, VERTEX_DATA, & 
        !                                            chunks(chunk)%yvel1_left_top_snd_buffer, chunks(chunk)%yvel1_right_bottom_rcv_buffer, chunks(chunk)%field%yvel1)
        CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%yvel1_left_top_snd_buffer, chunks(chunk)%field%yvel1)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_top_neighbour_chunk)[receiver]%yvel1_right_bottom_rcv_buffer(1:size) = chunks(chunk)%yvel1_left_top_snd_buffer(1:size)
    ENDIF

    x_inc=1
    y_inc=0
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        !CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, X_FACE_DATA, & 
        !                                            chunks(chunk)%volflux_x_left_top_snd_buffer, chunks(chunk)%volflux_x_right_bottom_rcv_buffer, chunks(chunk)%field%vol_flux_x)
        CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%volflux_x_left_top_snd_buffer, chunks(chunk)%field%vol_flux_x)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_top_neighbour_chunk)[receiver]%volflux_x_right_bottom_rcv_buffer(1:size) = chunks(chunk)%volflux_x_left_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        !CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, X_FACE_DATA, & 
        !                                            chunks(chunk)%massflux_x_left_top_snd_buffer, chunks(chunk)%massflux_x_right_bottom_rcv_buffer, chunks(chunk)%field%mass_flux_x)
        CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%massflux_x_left_top_snd_buffer, chunks(chunk)%field%mass_flux_x)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_top_neighbour_chunk)[receiver]%massflux_x_right_bottom_rcv_buffer(1:size) = chunks(chunk)%massflux_x_left_top_snd_buffer(1:size)
    ENDIF

    x_inc=0
    y_inc=1
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        !CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, Y_FACE_DATA, & 
        !                                            chunks(chunk)%volflux_y_left_top_snd_buffer, chunks(chunk)%volflux_y_right_bottom_rcv_buffer, chunks(chunk)%field%vol_flux_y)
        CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%volflux_y_left_top_snd_buffer, chunks(chunk)%field%vol_flux_y)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_top_neighbour_chunk)[receiver]%volflux_y_right_bottom_rcv_buffer(1:size) = chunks(chunk)%volflux_y_left_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        !CALL clover_exchange_write_message_left_top(chunk, depth, receiver, size, Y_FACE_DATA, & 
        !                                            chunks(chunk)%massflux_y_left_top_snd_buffer, chunks(chunk)%massflux_y_right_bottom_rcv_buffer, chunks(chunk)%field%mass_flux_y)
        CALL pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%massflux_y_left_top_snd_buffer, chunks(chunk)%field%mass_flux_y)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
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

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_top_neighbour_chunk)[receiver]%density0_left_bottom_rcv_buffer(1:size) = chunks(chunk)%density0_right_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        !CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, CELL_DATA, & 
        !                                             chunks(chunk)%density1_right_top_snd_buffer, chunks(chunk)%density1_left_bottom_rcv_buffer, chunks(chunk)%field%density1)
        CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%density1_right_top_snd_buffer, chunks(chunk)%field%density1)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_top_neighbour_chunk)[receiver]%density1_left_bottom_rcv_buffer(1:size) = chunks(chunk)%density1_right_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        !CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, CELL_DATA, & 
        !                                             chunks(chunk)%energy0_right_top_snd_buffer, chunks(chunk)%energy0_left_bottom_rcv_buffer, chunks(chunk)%field%energy0)
        CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%energy0_right_top_snd_buffer, chunks(chunk)%field%energy0)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_top_neighbour_chunk)[receiver]%energy0_left_bottom_rcv_buffer(1:size) = chunks(chunk)%energy0_right_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        !CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, CELL_DATA, & 
        !                                             chunks(chunk)%energy1_right_top_snd_buffer, chunks(chunk)%energy1_left_bottom_rcv_buffer, chunks(chunk)%field%energy1)
        CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%energy1_right_top_snd_buffer, chunks(chunk)%field%energy1)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_top_neighbour_chunk)[receiver]%energy1_left_bottom_rcv_buffer(1:size) = chunks(chunk)%energy1_right_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        !CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, CELL_DATA, & 
        !                                             chunks(chunk)%pressure_right_top_snd_buffer, chunks(chunk)%pressure_left_bottom_rcv_buffer, chunks(chunk)%field%pressure)
        CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%pressure_right_top_snd_buffer, chunks(chunk)%field%pressure)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_top_neighbour_chunk)[receiver]%pressure_left_bottom_rcv_buffer(1:size) = chunks(chunk)%pressure_right_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        !CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, CELL_DATA, & 
        !                                             chunks(chunk)%viscosity_right_top_snd_buffer, chunks(chunk)%viscosity_left_bottom_rcv_buffer, chunks(chunk)%field%viscosity)
        CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%viscosity_right_top_snd_buffer, chunks(chunk)%field%viscosity)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_top_neighbour_chunk)[receiver]%viscosity_left_bottom_rcv_buffer(1:size) = chunks(chunk)%viscosity_right_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        !CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, CELL_DATA, & 
        !                                             chunks(chunk)%soundspeed_right_top_snd_buffer, chunks(chunk)%soundspeed_left_bottom_rcv_buffer, chunks(chunk)%field%soundspeed)
        CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%soundspeed_right_top_snd_buffer, chunks(chunk)%field%soundspeed)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_top_neighbour_chunk)[receiver]%soundspeed_left_bottom_rcv_buffer(1:size) = chunks(chunk)%soundspeed_right_top_snd_buffer(1:size)
    ENDIF

    x_inc=1
    y_inc=1
    IF(fields(FIELD_XVEL0).EQ.1) THEN
        !CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, VERTEX_DATA, & 
        !                                             chunks(chunk)%xvel0_right_top_snd_buffer, chunks(chunk)%xvel0_left_bottom_rcv_buffer, chunks(chunk)%field%xvel0)
        CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%xvel0_right_top_snd_buffer, chunks(chunk)%field%xvel0)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_top_neighbour_chunk)[receiver]%xvel0_left_bottom_rcv_buffer(1:size) = chunks(chunk)%xvel0_right_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        !CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, VERTEX_DATA, & 
        !                                             chunks(chunk)%xvel1_right_top_snd_buffer, chunks(chunk)%xvel1_left_bottom_rcv_buffer, chunks(chunk)%field%xvel1)
        CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%xvel1_right_top_snd_buffer, chunks(chunk)%field%xvel1)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_top_neighbour_chunk)[receiver]%xvel1_left_bottom_rcv_buffer(1:size) = chunks(chunk)%xvel1_right_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        !CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, VERTEX_DATA, & 
        !                                             chunks(chunk)%yvel0_right_top_snd_buffer, chunks(chunk)%yvel0_left_bottom_rcv_buffer, chunks(chunk)%field%yvel0)
        CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%yvel0_right_top_snd_buffer, chunks(chunk)%field%yvel0)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_top_neighbour_chunk)[receiver]%yvel0_left_bottom_rcv_buffer(1:size) = chunks(chunk)%yvel0_right_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        !CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, VERTEX_DATA, & 
        !                                             chunks(chunk)%yvel1_right_top_snd_buffer, chunks(chunk)%yvel1_left_bottom_rcv_buffer, chunks(chunk)%field%yvel1)
        CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%yvel1_right_top_snd_buffer, chunks(chunk)%field%yvel1)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_top_neighbour_chunk)[receiver]%yvel1_left_bottom_rcv_buffer(1:size) = chunks(chunk)%yvel1_right_top_snd_buffer(1:size)
    ENDIF

    x_inc=1
    y_inc=0
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        !CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, X_FACE_DATA, & 
        !                                             chunks(chunk)%volflux_x_right_top_snd_buffer, chunks(chunk)%volflux_x_left_bottom_rcv_buffer, chunks(chunk)%field%vol_flux_x)
        CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%volflux_x_right_top_snd_buffer, chunks(chunk)%field%vol_flux_x)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_top_neighbour_chunk)[receiver]%volflux_x_left_bottom_rcv_buffer(1:size) = chunks(chunk)%volflux_x_right_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        !CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, X_FACE_DATA, & 
        !                                             chunks(chunk)%massflux_x_right_top_snd_buffer, chunks(chunk)%massflux_x_left_bottom_rcv_buffer, chunks(chunk)%field%mass_flux_x)
        CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%massflux_x_right_top_snd_buffer, chunks(chunk)%field%mass_flux_x)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_top_neighbour_chunk)[receiver]%massflux_x_left_bottom_rcv_buffer(1:size) = chunks(chunk)%massflux_x_right_top_snd_buffer(1:size)
    ENDIF

    x_inc=0
    y_inc=1
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        !CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, Y_FACE_DATA, & 
        !                                             chunks(chunk)%volflux_y_right_top_snd_buffer, chunks(chunk)%volflux_y_left_bottom_rcv_buffer, chunks(chunk)%field%vol_flux_y)
        CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%volflux_y_right_top_snd_buffer, chunks(chunk)%field%vol_flux_y)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_top_neighbour_chunk)[receiver]%volflux_y_left_bottom_rcv_buffer(1:size) = chunks(chunk)%volflux_y_right_top_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        !CALL clover_exchange_write_message_right_top(chunk, depth, receiver, size, Y_FACE_DATA, & 
        !                                             chunks(chunk)%massflux_y_right_top_snd_buffer, chunks(chunk)%massflux_y_left_bottom_rcv_buffer, chunks(chunk)%field%mass_flux_y)
        CALL pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%massflux_y_right_top_snd_buffer, chunks(chunk)%field%mass_flux_y)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
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

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_bottom_neighbour_chunk)[receiver]%density0_left_top_rcv_buffer(1:size) = chunks(chunk)%density0_right_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        !CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, CELL_DATA, & 
        !                                                chunks(chunk)%density1_right_bottom_snd_buffer, chunks(chunk)%density1_left_top_rcv_buffer, chunks(chunk)%field%density1)
        CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%density1_right_bottom_snd_buffer, chunks(chunk)%field%density1)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_bottom_neighbour_chunk)[receiver]%density1_left_top_rcv_buffer(1:size) = chunks(chunk)%density1_right_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        !CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, CELL_DATA, & 
        !                                                chunks(chunk)%energy0_right_bottom_snd_buffer, chunks(chunk)%energy0_left_top_rcv_buffer, chunks(chunk)%field%energy0)
        CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%energy0_right_bottom_snd_buffer, chunks(chunk)%field%energy0)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_bottom_neighbour_chunk)[receiver]%energy0_left_top_rcv_buffer(1:size) = chunks(chunk)%energy0_right_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        !CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, CELL_DATA, & 
        !                                                chunks(chunk)%energy1_right_bottom_snd_buffer, chunks(chunk)%energy1_left_top_rcv_buffer, chunks(chunk)%field%energy1)
        CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%energy1_right_bottom_snd_buffer, chunks(chunk)%field%energy1)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_bottom_neighbour_chunk)[receiver]%energy1_left_top_rcv_buffer(1:size) = chunks(chunk)%energy1_right_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        !CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, CELL_DATA, & 
        !                                                chunks(chunk)%pressure_right_bottom_snd_buffer, chunks(chunk)%pressure_left_top_rcv_buffer, chunks(chunk)%field%pressure)
        CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%pressure_right_bottom_snd_buffer, chunks(chunk)%field%pressure)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_bottom_neighbour_chunk)[receiver]%pressure_left_top_rcv_buffer(1:size) = chunks(chunk)%pressure_right_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        !CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, CELL_DATA, & 
        !                                                chunks(chunk)%viscosity_right_bottom_snd_buffer, chunks(chunk)%viscosity_left_top_rcv_buffer, chunks(chunk)%field%viscosity)
        CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%viscosity_right_bottom_snd_buffer, chunks(chunk)%field%viscosity)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_bottom_neighbour_chunk)[receiver]%viscosity_left_top_rcv_buffer(1:size) = chunks(chunk)%viscosity_right_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        !CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, CELL_DATA, & 
        !                                                chunks(chunk)%soundspeed_right_bottom_snd_buffer, chunks(chunk)%soundspeed_left_top_rcv_buffer, chunks(chunk)%field%soundspeed)
        CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%soundspeed_right_bottom_snd_buffer, chunks(chunk)%field%soundspeed)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_bottom_neighbour_chunk)[receiver]%soundspeed_left_top_rcv_buffer(1:size) = chunks(chunk)%soundspeed_right_bottom_snd_buffer(1:size)
    ENDIF

    x_inc=1
    y_inc=1
    IF(fields(FIELD_XVEL0).EQ.1) THEN
        !CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, VERTEX_DATA, & 
        !                                                chunks(chunk)%xvel0_right_bottom_snd_buffer, chunks(chunk)%xvel0_left_top_rcv_buffer, chunks(chunk)%field%xvel0)
        CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%xvel0_right_bottom_snd_buffer, chunks(chunk)%field%xvel0)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_bottom_neighbour_chunk)[receiver]%xvel0_left_top_rcv_buffer(1:size) = chunks(chunk)%xvel0_right_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        !CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, VERTEX_DATA, & 
        !                                                chunks(chunk)%xvel1_right_bottom_snd_buffer, chunks(chunk)%xvel1_left_top_rcv_buffer, chunks(chunk)%field%xvel1)
        CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%xvel1_right_bottom_snd_buffer, chunks(chunk)%field%xvel1)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_bottom_neighbour_chunk)[receiver]%xvel1_left_top_rcv_buffer(1:size) = chunks(chunk)%xvel1_right_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        !CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, VERTEX_DATA, & 
        !                                                chunks(chunk)%yvel0_right_bottom_snd_buffer, chunks(chunk)%yvel0_left_top_rcv_buffer, chunks(chunk)%field%yvel0)
        CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%yvel0_right_bottom_snd_buffer, chunks(chunk)%field%yvel0)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_bottom_neighbour_chunk)[receiver]%yvel0_left_top_rcv_buffer(1:size) = chunks(chunk)%yvel0_right_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        !CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, VERTEX_DATA, & 
        !                                                chunks(chunk)%yvel1_right_bottom_snd_buffer, chunks(chunk)%yvel1_left_top_rcv_buffer, chunks(chunk)%field%yvel1)
        CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%yvel1_right_bottom_snd_buffer, chunks(chunk)%field%yvel1)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_bottom_neighbour_chunk)[receiver]%yvel1_left_top_rcv_buffer(1:size) = chunks(chunk)%yvel1_right_bottom_snd_buffer(1:size)
    ENDIF

    x_inc=1
    y_inc=0
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        !CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, X_FACE_DATA, & 
        !                                                chunks(chunk)%volflux_x_right_bottom_snd_buffer, chunks(chunk)%volflux_x_left_top_rcv_buffer, chunks(chunk)%field%vol_flux_x)
        CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%volflux_x_right_bottom_snd_buffer, chunks(chunk)%field%vol_flux_x)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_bottom_neighbour_chunk)[receiver]%volflux_x_left_top_rcv_buffer(1:size) = chunks(chunk)%volflux_x_right_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        !CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, X_FACE_DATA, & 
        !                                                chunks(chunk)%massflux_x_right_bottom_snd_buffer, chunks(chunk)%massflux_x_left_top_rcv_buffer, chunks(chunk)%field%mass_flux_x)
        CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%massflux_x_right_bottom_snd_buffer, chunks(chunk)%field%mass_flux_x)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_bottom_neighbour_chunk)[receiver]%massflux_x_left_top_rcv_buffer(1:size) = chunks(chunk)%massflux_x_right_bottom_snd_buffer(1:size)
    ENDIF

    x_inc=0
    y_inc=1
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        !CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, Y_FACE_DATA, & 
        !                                                chunks(chunk)%volflux_y_right_bottom_snd_buffer, chunks(chunk)%volflux_y_left_top_rcv_buffer, chunks(chunk)%field%vol_flux_y)
        CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%volflux_y_right_bottom_snd_buffer, chunks(chunk)%field%vol_flux_y)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(right_bottom_neighbour_chunk)[receiver]%volflux_y_left_top_rcv_buffer(1:size) = chunks(chunk)%volflux_y_right_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        !CALL clover_exchange_write_message_right_bottom(chunk, depth, receiver, size, Y_FACE_DATA, & 
        !                                                chunks(chunk)%massflux_y_right_bottom_snd_buffer, chunks(chunk)%massflux_y_left_top_rcv_buffer, chunks(chunk)%field%mass_flux_y)
        CALL pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%massflux_y_right_bottom_snd_buffer, chunks(chunk)%field%mass_flux_y)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
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

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_bottom_neighbour_chunk)[receiver]%density0_right_top_rcv_buffer(1:size) = chunks(chunk)%density0_left_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_DENSITY1).EQ.1) THEN
        !CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, CELL_DATA, & 
        !                                               chunks(chunk)%density1_left_bottom_snd_buffer, chunks(chunk)%density1_right_top_rcv_buffer, chunks(chunk)%field%density1)
        CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%density1_left_bottom_snd_buffer, chunks(chunk)%field%density1)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_bottom_neighbour_chunk)[receiver]%density1_right_top_rcv_buffer(1:size) = chunks(chunk)%density1_left_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_ENERGY0).EQ.1) THEN
        !CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, CELL_DATA, & 
        !                                               chunks(chunk)%energy0_left_bottom_snd_buffer, chunks(chunk)%energy0_right_top_rcv_buffer, chunks(chunk)%field%energy0)
        CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%energy0_left_bottom_snd_buffer, chunks(chunk)%field%energy0)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_bottom_neighbour_chunk)[receiver]%energy0_right_top_rcv_buffer(1:size) = chunks(chunk)%energy0_left_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_ENERGY1).EQ.1) THEN
        !CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, CELL_DATA, & 
        !                                               chunks(chunk)%energy1_left_bottom_snd_buffer, chunks(chunk)%energy1_right_top_rcv_buffer, chunks(chunk)%field%energy1)
        CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%energy1_left_bottom_snd_buffer, chunks(chunk)%field%energy1)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_bottom_neighbour_chunk)[receiver]%energy1_right_top_rcv_buffer(1:size) = chunks(chunk)%energy1_left_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_PRESSURE).EQ.1) THEN
        !CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, CELL_DATA, & 
        !                                               chunks(chunk)%pressure_left_bottom_snd_buffer, chunks(chunk)%pressure_right_top_rcv_buffer, chunks(chunk)%field%pressure)
        CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%pressure_left_bottom_snd_buffer, chunks(chunk)%field%pressure)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_bottom_neighbour_chunk)[receiver]%pressure_right_top_rcv_buffer(1:size) = chunks(chunk)%pressure_left_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_VISCOSITY).EQ.1) THEN
        !CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, CELL_DATA, & 
        !                                               chunks(chunk)%viscosity_left_bottom_snd_buffer, chunks(chunk)%viscosity_right_top_rcv_buffer, chunks(chunk)%field%viscosity)
        CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%viscosity_left_bottom_snd_buffer, chunks(chunk)%field%viscosity)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_bottom_neighbour_chunk)[receiver]%viscosity_right_top_rcv_buffer(1:size) = chunks(chunk)%viscosity_left_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_SOUNDSPEED).EQ.1) THEN
        !CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, CELL_DATA, & 
        !                                               chunks(chunk)%soundspeed_left_bottom_snd_buffer, chunks(chunk)%soundspeed_right_top_rcv_buffer, chunks(chunk)%field%soundspeed)
        CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%soundspeed_left_bottom_snd_buffer, chunks(chunk)%field%soundspeed)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_bottom_neighbour_chunk)[receiver]%soundspeed_right_top_rcv_buffer(1:size) = chunks(chunk)%soundspeed_left_bottom_snd_buffer(1:size)
    ENDIF

    x_inc=1
    y_inc=1
    IF(fields(FIELD_XVEL0).EQ.1) THEN
        !CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, VERTEX_DATA, & 
        !                                               chunks(chunk)%xvel0_left_bottom_snd_buffer, chunks(chunk)%xvel0_right_top_rcv_buffer, chunks(chunk)%field%xvel0)
        CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%xvel0_left_bottom_snd_buffer, chunks(chunk)%field%xvel0)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_bottom_neighbour_chunk)[receiver]%xvel0_right_top_rcv_buffer(1:size) = chunks(chunk)%xvel0_left_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_XVEL1).EQ.1) THEN
        !CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, VERTEX_DATA, & 
        !                                               chunks(chunk)%xvel1_left_bottom_snd_buffer, chunks(chunk)%xvel1_right_top_rcv_buffer, chunks(chunk)%field%xvel1)
        CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%xvel1_left_bottom_snd_buffer, chunks(chunk)%field%xvel1)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_bottom_neighbour_chunk)[receiver]%xvel1_right_top_rcv_buffer(1:size) = chunks(chunk)%xvel1_left_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_YVEL0).EQ.1) THEN
        !CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, VERTEX_DATA, & 
        !                                               chunks(chunk)%yvel0_left_bottom_snd_buffer, chunks(chunk)%yvel0_right_top_rcv_buffer, chunks(chunk)%field%yvel0)
        CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%yvel0_left_bottom_snd_buffer, chunks(chunk)%field%yvel0)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_bottom_neighbour_chunk)[receiver]%yvel0_right_top_rcv_buffer(1:size) = chunks(chunk)%yvel0_left_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_YVEL1).EQ.1) THEN
        !CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, VERTEX_DATA, & 
        !                                               chunks(chunk)%yvel1_left_bottom_snd_buffer, chunks(chunk)%yvel1_right_top_rcv_buffer, chunks(chunk)%field%yvel1)
        CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%yvel1_left_bottom_snd_buffer, chunks(chunk)%field%yvel1)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_bottom_neighbour_chunk)[receiver]%yvel1_right_top_rcv_buffer(1:size) = chunks(chunk)%yvel1_left_bottom_snd_buffer(1:size)
    ENDIF

    x_inc=1
    y_inc=0
    IF(fields(FIELD_VOL_FLUX_X).EQ.1) THEN
        !CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, X_FACE_DATA, & 
        !                                               chunks(chunk)%volflux_x_left_bottom_snd_buffer, chunks(chunk)%volflux_x_right_top_rcv_buffer, chunks(chunk)%field%vol_flux_x)
        CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%volflux_x_left_bottom_snd_buffer, chunks(chunk)%field%vol_flux_x)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_bottom_neighbour_chunk)[receiver]%volflux_x_right_top_rcv_buffer(1:size) = chunks(chunk)%volflux_x_left_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_X).EQ.1) THEN
        !CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, X_FACE_DATA, & 
        !                                               chunks(chunk)%massflux_x_left_bottom_snd_buffer, chunks(chunk)%massflux_x_right_top_rcv_buffer, chunks(chunk)%field%mass_flux_x)
        CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%massflux_x_left_bottom_snd_buffer, chunks(chunk)%field%mass_flux_x)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_bottom_neighbour_chunk)[receiver]%massflux_x_right_top_rcv_buffer(1:size) = chunks(chunk)%massflux_x_left_bottom_snd_buffer(1:size)
    ENDIF

    x_inc=0
    y_inc=1
    IF(fields(FIELD_VOL_FLUX_Y).EQ.1) THEN
        !CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, Y_FACE_DATA, & 
        !                                               chunks(chunk)%volflux_y_left_bottom_snd_buffer, chunks(chunk)%volflux_y_right_top_rcv_buffer, chunks(chunk)%field%vol_flux_y)
        CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%volflux_y_left_bottom_snd_buffer, chunks(chunk)%field%vol_flux_y)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_bottom_neighbour_chunk)[receiver]%volflux_y_right_top_rcv_buffer(1:size) = chunks(chunk)%volflux_y_left_bottom_snd_buffer(1:size)
    ENDIF
    IF(fields(FIELD_MASS_FLUX_Y).EQ.1) THEN
        !CALL clover_exchange_write_message_left_bottom(chunk, depth, receiver, size, Y_FACE_DATA, & 
        !                                               chunks(chunk)%massflux_y_left_bottom_snd_buffer, chunks(chunk)%massflux_y_right_top_rcv_buffer, chunks(chunk)%field%mass_flux_y)
        CALL pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, chunks(chunk)%massflux_y_left_bottom_snd_buffer, chunks(chunk)%field%mass_flux_y)

#ifdef DEF_SYNC
        !dir$ pgas defer_sync        
#endif
        chunks(left_bottom_neighbour_chunk)[receiver]%massflux_y_right_top_rcv_buffer(1:size) = chunks(chunk)%massflux_y_left_bottom_snd_buffer(1:size)
    ENDIF

END SUBROUTINE clover_exchange_write_all_buffers_left_bottom












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
        CALL unpack_left_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%yvel0_left_rcv_buffer, chunks(chunk)%field%yvel0)
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
        CALL unpack_right_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%yvel0_right_rcv_buffer, chunks(chunk)%field%yvel0)
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
        CALL unpack_bottom_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%yvel0_bottom_rcv_buffer, chunks(chunk)%field%yvel0)
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
        CALL unpack_top_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%yvel0_top_rcv_buffer, chunks(chunk)%field%yvel0)
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
        CALL unpack_left_top_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%yvel0_left_top_rcv_buffer, chunks(chunk)%field%yvel0)
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
        CALL unpack_right_top_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%yvel0_right_top_rcv_buffer, chunks(chunk)%field%yvel0)
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
        CALL unpack_right_bottom_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%yvel0_right_bottom_rcv_buffer, chunks(chunk)%field%yvel0)
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
        CALL unpack_left_bottom_buffer_seq(chunk, depth, VERTEX_DATA, chunks(chunk)%yvel0_left_bottom_rcv_buffer, chunks(chunk)%field%yvel0)
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

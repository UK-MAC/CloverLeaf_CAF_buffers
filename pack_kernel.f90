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

!>  @brief Fortran mpi buffer packing kernel
!>  @author Wayne Gaudin
!>  @details Packs/unpacks mpi send and receive buffers

MODULE pack_kernel_module

USE definitions_module 

CONTAINS

SUBROUTINE pack_left_buffer_seq(chunk, depth, x_inc, y_inc, left_snd_buffer, field)

    IMPLICIT NONE

    INTEGER :: topedge, bottomedge, index, chunk, x_inc, y_inc, depth, j, k
    REAL(KIND=8) :: left_snd_buffer(:)
    REAL(KIND=8) :: field(-1:,-1:) 

    topedge = 0
    bottomedge = 0
    IF (chunks(chunk)%chunk_neighbours(chunk_top).EQ.external_face) THEN
      topedge = depth
    ENDIF
    IF (chunks(chunk)%chunk_neighbours(chunk_bottom).EQ.external_face) THEN
      bottomedge = depth
    ENDIF

    index=1
    DO k=chunks(chunk)%field%y_min-bottomedge,chunks(chunk)%field%y_max+y_inc+topedge
      DO j=1,depth
        left_snd_buffer(index)=field(chunks(chunk)%field%x_min+x_inc-1+j,k)
        index = index + 1
      ENDDO
    ENDDO

END SUBROUTINE pack_left_buffer_seq

SUBROUTINE pack_right_buffer_seq(chunk, depth, x_inc, y_inc, right_snd_buffer, field)

    IMPLICIT NONE

    INTEGER :: topedge, bottomedge, index, chunk, x_inc, y_inc, depth, j, k
    REAL(KIND=8) :: right_snd_buffer(:)
    REAL(KIND=8) :: field(-1:,-1:) 

    topedge = 0
    bottomedge = 0
    IF (chunks(chunk)%chunk_neighbours(chunk_top).EQ.external_face) THEN
      topedge = depth
    ENDIF
    IF (chunks(chunk)%chunk_neighbours(chunk_bottom).EQ.external_face) THEN
      bottomedge = depth
    ENDIF

    index=1
    DO k=chunks(chunk)%field%y_min-bottomedge,chunks(chunk)%field%y_max+y_inc+topedge
      DO j=1,depth
        right_snd_buffer(index)=field(chunks(chunk)%field%x_max+1-j,k)
        index = index + 1
      ENDDO
    ENDDO

END SUBROUTINE pack_right_buffer_seq

SUBROUTINE pack_bottom_buffer_seq(chunk, depth, x_inc, y_inc, bottom_snd_buffer, field)

    IMPLICIT NONE

    INTEGER :: leftedge, rightedge, index, chunk, x_inc, y_inc, depth, j, k
    REAL(KIND=8) :: bottom_snd_buffer(:)
    REAL(KIND=8) :: field(-1:,-1:)

    leftedge= 0
    rightedge= 0
    IF (chunks(chunk)%chunk_neighbours(chunk_left).EQ.external_face) THEN
      leftedge = depth
    ENDIF
    IF (chunks(chunk)%chunk_neighbours(chunk_right).EQ.external_face) THEN
      rightedge = depth
    ENDIF

    index = 1
    DO k=1,depth
      DO j=chunks(chunk)%field%x_min-leftedge,chunks(chunk)%field%x_max+x_inc+rightedge
        bottom_snd_buffer(index)=field(j,chunks(chunk)%field%y_min+y_inc-1+k)
        index = index + 1
      ENDDO
    ENDDO

END SUBROUTINE pack_bottom_buffer_seq

SUBROUTINE pack_top_buffer_seq(chunk, depth, x_inc, y_inc, top_snd_buffer, field)

    IMPLICIT NONE

    INTEGER :: leftedge, rightedge, index, chunk, x_inc, y_inc, depth, j, k
    REAL(KIND=8) :: top_snd_buffer(:)
    REAL(KIND=8) :: field(-1:,-1:)

    leftedge= 0
    rightedge= 0
    IF (chunks(chunk)%chunk_neighbours(chunk_left).EQ.external_face) THEN
      leftedge = depth
    ENDIF
    IF (chunks(chunk)%chunk_neighbours(chunk_right).EQ.external_face) THEN
      rightedge = depth
    ENDIF

    index = 1
    DO k=1,depth
      DO j=chunks(chunk)%field%x_min-leftedge,chunks(chunk)%field%x_max+x_inc+rightedge
        top_snd_buffer(index)=field(j,chunks(chunk)%field%y_max+1-k)
        index = index + 1
      ENDDO
    ENDDO

END SUBROUTINE pack_top_buffer_seq

SUBROUTINE pack_left_top_buffer_seq(chunk, depth, x_inc, y_inc, left_top_snd_buffer, field)

    IMPLICIT NONE

    INTEGER :: chunk, depth, x_inc, y_inc, index, j, k
    REAL(KIND=8) :: left_top_snd_buffer(:)
    REAL(KIND=8) :: field(-1:,-1:)

    index = 1
    DO k=1,depth
        DO j=1,depth
            left_top_snd_buffer(index) = field(chunks(chunk)%field%x_min+x_inc-1+j,chunks(chunk)%field%y_max+1-k)
            index = index + 1
        ENDDO
    ENDDO

END SUBROUTINE pack_left_top_buffer_seq

SUBROUTINE pack_right_top_buffer_seq(chunk, depth, x_inc, y_inc, right_top_snd_buffer, field)

    IMPLICIT NONE

    INTEGER :: chunk, depth, x_inc, y_inc, index, j, k
    REAL(KIND=8) :: right_top_snd_buffer(:)
    REAL(KIND=8) :: field(-1:,-1:)

    index = 1
    DO k=1,depth
        DO j=1,depth
            right_top_snd_buffer(index) = field(chunks(chunk)%field%x_max+1-j,chunks(chunk)%field%y_max+1-k)
            index = index + 1
        ENDDO
    ENDDO

END SUBROUTINE pack_right_top_buffer_seq

SUBROUTINE pack_right_bottom_buffer_seq(chunk, depth, x_inc, y_inc, right_bottom_snd_buffer, field)

    IMPLICIT NONE

    INTEGER :: chunk, depth, x_inc, y_inc, index, j, k
    REAL(KIND=8) :: right_bottom_snd_buffer(:)
    REAL(KIND=8) :: field(-1:,-1:)

    index = 1
    DO k=1,depth
        DO j=1,depth
            right_bottom_snd_buffer(index) = field(chunks(chunk)%field%x_max+1-j,chunks(chunk)%field%y_min+y_inc-1+k)
            index = index + 1
        ENDDO
    ENDDO

END SUBROUTINE pack_right_bottom_buffer_seq

SUBROUTINE pack_left_bottom_buffer_seq(chunk, depth, x_inc, y_inc, left_bottom_snd_buffer, field)

    IMPLICIT NONE

    INTEGER :: chunk, depth, x_inc, y_inc, index, j, k
    REAL(KIND=8) :: left_bottom_snd_buffer(:)
    REAL(KIND=8) :: field(-1:,-1:)

    index = 1
    DO k=1,depth
        DO j=1,depth
            left_bottom_snd_buffer(index) = field(chunks(chunk)%field%x_min+x_inc-1+j,chunks(chunk)%field%y_min+y_inc-1+k)
            index = index + 1
        ENDDO
    ENDDO


END SUBROUTINE pack_left_bottom_buffer_seq


!new unpack routines
SUBROUTINE unpack_left_buffer_seq(chunk, depth, field_type, left_rcv_buffer, field)

    IMPLICIT NONE

    INTEGER :: chunk, depth, x_inc, y_inc, index, topedge, bottomedge, j, k, field_type
    REAL(KIND=8) :: left_rcv_buffer(:)
    REAL(KIND=8) :: field(-1:,-1:)

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

    index = 1
    topedge = 0
    bottomedge = 0
    IF (chunks(chunk)%chunk_neighbours(chunk_top).EQ.external_face) THEN
      topedge = depth
    ENDIF
    IF (chunks(chunk)%chunk_neighbours(chunk_bottom).EQ.external_face) THEN
      bottomedge = depth
    ENDIF

    DO k=chunks(chunk)%field%y_min-bottomedge,chunks(chunk)%field%y_max+y_inc+topedge
      DO j=1,depth
        field(chunks(chunk)%field%x_min-j,k)=left_rcv_buffer(index)
        index = index + 1
      ENDDO
    ENDDO

END SUBROUTINE unpack_left_buffer_seq


SUBROUTINE unpack_right_buffer_seq(chunk, depth, field_type, right_rcv_buffer, field)

    IMPLICIT NONE

    INTEGER :: chunk, depth, x_inc, y_inc, index, topedge, bottomedge, j, k, field_type
    REAL(KIND=8) :: right_rcv_buffer(:)
    REAL(KIND=8) :: field(-1:,-1:)

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

    index = 1
    topedge = 0
    bottomedge = 0
    IF (chunks(chunk)%chunk_neighbours(chunk_top).EQ.external_face) THEN
      topedge = depth
    ENDIF
    IF (chunks(chunk)%chunk_neighbours(chunk_bottom).EQ.external_face) THEN
      bottomedge = depth
    ENDIF

    DO k=chunks(chunk)%field%y_min-bottomedge,chunks(chunk)%field%y_max+y_inc+topedge
      DO j=1,depth
        field(chunks(chunk)%field%x_max+x_inc+j,k)=right_rcv_buffer(index)
        index = index + 1
      ENDDO
    ENDDO

END SUBROUTINE unpack_right_buffer_seq


SUBROUTINE unpack_bottom_buffer_seq(chunk, depth, field_type, bottom_rcv_buffer, field)

    IMPLICIT NONE

    INTEGER :: chunk, depth, x_inc, y_inc, index, leftedge, rightedge, j, k, field_type
    REAL(KIND=8) :: bottom_rcv_buffer(:)
    REAL(KIND=8) :: field(-1:,-1:)

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

    index = 1
    leftedge = 0
    rightedge = 0
    IF (chunks(chunk)%chunk_neighbours(chunk_left).EQ.external_face) THEN
      leftedge = depth
    ENDIF
    IF (chunks(chunk)%chunk_neighbours(chunk_right).EQ.external_face) THEN
      rightedge = depth
    ENDIF

    DO k=1,depth
      DO j=chunks(chunk)%field%x_min-leftedge,chunks(chunk)%field%x_max+x_inc+rightedge
        field(j,chunks(chunk)%field%y_min-k)=bottom_rcv_buffer(index)
        index = index + 1
      ENDDO
    ENDDO

END SUBROUTINE unpack_bottom_buffer_seq


SUBROUTINE unpack_top_buffer_seq(chunk, depth, field_type, top_rcv_buffer, field)

    IMPLICIT NONE

    INTEGER :: chunk, depth, x_inc, y_inc, index, leftedge, rightedge, j, k, field_type
    REAL(KIND=8) :: top_rcv_buffer(:)
    REAL(KIND=8) :: field(-1:,-1:)

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

    index = 1
    leftedge = 0
    rightedge = 0
    IF (chunks(chunk)%chunk_neighbours(chunk_left).EQ.external_face) THEN
      leftedge = depth
    ENDIF
    IF (chunks(chunk)%chunk_neighbours(chunk_right).EQ.external_face) THEN
      rightedge = depth
    ENDIF

    DO k=1,depth
      DO j=chunks(chunk)%field%x_min-leftedge,chunks(chunk)%field%x_max+x_inc+rightedge
        field(j,chunks(chunk)%field%y_max+y_inc+k)=top_rcv_buffer(index)
        index = index + 1
      ENDDO
    ENDDO

END SUBROUTINE unpack_top_buffer_seq


SUBROUTINE unpack_left_top_buffer_seq(chunk, depth, field_type, left_top_rcv_buffer, field)

    IMPLICIT NONE

    INTEGER :: chunk, depth, x_inc, y_inc, index, j, k, field_type
    REAL(KIND=8) :: left_top_rcv_buffer(:)
    REAL(KIND=8) :: field(-1:,-1:)

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

    index = 1
    DO k=1,depth
        DO j=1,depth
            field(chunks(chunk)%field%x_min-j,chunks(chunk)%field%y_max+y_inc+k)=left_top_rcv_buffer(index)
            index = index + 1
        ENDDO
    ENDDO

END SUBROUTINE unpack_left_top_buffer_seq

SUBROUTINE unpack_right_top_buffer_seq(chunk, depth, field_type, right_top_rcv_buffer, field)

    IMPLICIT NONE

    INTEGER :: chunk, depth, x_inc, y_inc, index, j, k, field_type
    REAL(KIND=8) :: right_top_rcv_buffer(:)
    REAL(KIND=8) :: field(-1:,-1:)

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

    index = 1
    DO k=1,depth
        DO j=1,depth
            field(chunks(chunk)%field%x_max+x_inc+j,chunks(chunk)%field%y_max+y_inc+k)=right_top_rcv_buffer(index)
            index = index + 1
        ENDDO
    ENDDO

END SUBROUTINE unpack_right_top_buffer_seq

SUBROUTINE unpack_right_bottom_buffer_seq(chunk, depth, field_type, right_bottom_rcv_buffer, field)

    IMPLICIT NONE

    INTEGER :: chunk, depth, x_inc, y_inc, index, j, k, field_type
    REAL(KIND=8) :: right_bottom_rcv_buffer(:)
    REAL(KIND=8) :: field(-1:,-1:)

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

    index = 1
    DO k=1,depth
        DO j=1,depth
            field(chunks(chunk)%field%x_max+x_inc+j,chunks(chunk)%field%y_min-k)=right_bottom_rcv_buffer(index)
            index = index + 1
        ENDDO
    ENDDO

END SUBROUTINE unpack_right_bottom_buffer_seq

SUBROUTINE unpack_left_bottom_buffer_seq(chunk, depth, field_type, left_bottom_rcv_buffer, field)

    IMPLICIT NONE

    INTEGER :: chunk, depth, x_inc, y_inc, index, j, k, field_type
    REAL(KIND=8) :: left_bottom_rcv_buffer(:)
    REAL(KIND=8) :: field(-1:,-1:)

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

    index = 1
    DO k=1,depth
        DO j=1,depth
            field(chunks(chunk)%field%x_min-j,chunks(chunk)%field%y_min-k)=left_bottom_rcv_buffer(index)
            index = index + 1
        ENDDO
    ENDDO

END SUBROUTINE unpack_left_bottom_buffer_seq



!old pack and unpack routines
!SUBROUTINE pack_left_right_buffers(x_min,x_max,y_min,y_max,              &
!                                   chunk_left,chunk_right,external_face, &
!                                   x_inc,y_inc,depth,size,               &
!                                   field,left_snd_buffer,right_snd_buffer)
!
!  IMPLICIT NONE
!
!  INTEGER      :: x_min,x_max,y_min,y_max
!  INTEGER      :: chunk_left,chunk_right,external_face
!  INTEGER      :: x_inc,y_inc,depth,size
!
!  REAL(KIND=8) :: field(-1:,-1:) ! This seems to work for any type of mesh data
!  REAL(KIND=8) :: left_snd_buffer(:),right_snd_buffer(:)
!
!  INTEGER      :: j,k,index
!
!  IF(chunk_left.NE.external_face) THEN
!    DO k=y_min-depth,y_max+y_inc+depth
!      DO j=1,depth
!        index=j+(k+depth-1)*depth
!        left_snd_buffer(index)=field(x_min+x_inc-1+j,k)
!      ENDDO
!    ENDDO
!  ENDIF
!  IF(chunk_right.NE.external_face) THEN
!    DO k=y_min-depth,y_max+y_inc+depth
!      DO j=1,depth
!        index=j+(k+depth-1)*depth
!        right_snd_buffer(index)=field(x_max+1-j,k)
!      ENDDO
!    ENDDO
!  ENDIF
!
!END SUBROUTINE pack_left_right_buffers
!
!SUBROUTINE unpack_left_right_buffers(x_min,x_max,y_min,y_max,              &
!                                     chunk_left,chunk_right,external_face, &
!                                     x_inc,y_inc,depth,size,               &
!                                     field,left_rcv_buffer,right_rcv_buffer)
!
!  IMPLICIT NONE
!
!  INTEGER      :: x_min,x_max,y_min,y_max
!  INTEGER      :: chunk_left,chunk_right,external_face
!  INTEGER      :: x_inc,y_inc,depth,size
!
!  REAL(KIND=8) :: field(-1:,-1:) ! This seems to work for any type of mesh data
!  REAL(KIND=8) :: left_rcv_buffer(:),right_rcv_buffer(:)
!
!  INTEGER      :: j,k,index
!
!  IF(chunk_left.NE.external_face) THEN
!    DO k=y_min-depth,y_max+y_inc+depth
!      DO j=1,depth
!        index=j+(k+depth-1)*depth
!        field(x_min-j,k)=left_rcv_buffer(index)
!      ENDDO
!    ENDDO
!  ENDIF
!  IF(chunk_right.NE.external_face) THEN
!    DO k=y_min-depth,y_max+y_inc+depth
!      DO j=1,depth
!        index=j+(k+depth-1)*depth
!        field(x_max+x_inc+j,k)=right_rcv_buffer(index)
!      ENDDO
!    ENDDO
!  ENDIF
!
!END SUBROUTINE unpack_left_right_buffers
!
!SUBROUTINE pack_top_bottom_buffers(x_min,x_max,y_min,y_max,              &
!                                   chunk_bottom,chunk_top,external_face, &
!                                   x_inc,y_inc,depth,size,               &
!                                   field,bottom_snd_buffer,top_snd_buffer)
!
!  IMPLICIT NONE
!
!  INTEGER      :: x_min,x_max,y_min,y_max
!  INTEGER      :: chunk_bottom,chunk_top,external_face
!  INTEGER      :: x_inc,y_inc,depth,size
!
!  REAL(KIND=8) :: field(-1:,-1:) ! This seems to work for any type of mesh data
!  REAL(KIND=8) :: bottom_snd_buffer(:),top_snd_buffer(:)
!
!  INTEGER      :: j,k,index
!
!  IF(chunk_bottom.NE.external_face) THEN
!    DO k=1,depth
!      DO j=x_min-depth,x_max+x_inc+depth
!        index=j+depth+(k-1)*(x_max+x_inc+(2*depth))
!        bottom_snd_buffer(index)=field(j,y_min+y_inc-1+k)
!      ENDDO
!    ENDDO
!  ENDIF
!  IF(chunk_top.NE.external_face) THEN
!    DO k=1,depth
!      DO j=x_min-depth,x_max+x_inc+depth
!        index=j+depth+(k-1)*(x_max+x_inc+(2*depth))
!        top_snd_buffer(index)=field(j,y_max+1-k)
!      ENDDO
!    ENDDO
!  ENDIF
!
!END SUBROUTINE pack_top_bottom_buffers
!
!SUBROUTINE unpack_top_bottom_buffers(x_min,x_max,y_min,y_max,             &
!                                    chunk_bottom,chunk_top,external_face, &
!                                    x_inc,y_inc,depth,size,               &
!                                    field,bottom_rcv_buffer,top_rcv_buffer)
!
!  IMPLICIT NONE
!
!  INTEGER      :: x_min,x_max,y_min,y_max
!  INTEGER      :: chunk_bottom,chunk_top,external_face
!  INTEGER      :: x_inc,y_inc,depth,size
!
!  REAL(KIND=8) :: field(-1:,-1:) ! This seems to work for any type of mesh data
!  REAL(KIND=8) :: bottom_rcv_buffer(:),top_rcv_buffer(:)
!
!  INTEGER      :: j,k,index
!
!  IF(chunk_bottom.NE.external_face) THEN
!    DO k=1,depth
!      DO j=x_min-depth,x_max+x_inc+depth
!        index=j+depth+(k-1)*(x_max+x_inc+(2*depth))
!        field(j,y_min-k)=bottom_rcv_buffer(index)
!      ENDDO
!    ENDDO
!  ENDIF
!  IF(chunk_top.NE.external_face) THEN
!    DO k=1,depth
!      DO j=x_min-depth,x_max+x_inc+depth
!        index=j+depth+(k-1)*(x_max+x_inc+(2*depth))
!        field(j,y_max+y_inc+k)=top_rcv_buffer(index)
!      ENDDO
!    ENDDO
!  ENDIF
!
!END SUBROUTINE unpack_top_bottom_buffers

END MODULE pack_kernel_module

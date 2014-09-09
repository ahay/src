c GETBUF - returns pointers to beginning and end of buffer segmentc
c
c WWS, 23.9.92
c
c Arguments:
c
c buffer      --->  buffer name
c pointer    <---   pointer assigned to beginning of segment
c length      --->  length of segment
c next       <--->  pointer to beginning of current segment on
c                   call, next segment on return
c last        --->  length of buffer
c ier        <--->  usual error flag
c
c---------------------------------------------------------------
c
      subroutine getbuf(buffer,pointer,length,next,last,ipdmp,ier)

      character*(*) buffer
      integer pointer,length,next,last,ier,ipdmp

      integer idbg       ! debug flag

      if (ier.ne.0) return

      idbg=0


      if (idbg.ne.0) then
         write(ipdmp,*)' GETBUF:'
         write(ipdmp,*)' buffer  = ',buffer
         write(ipdmp,*)' pointer = ',pointer
         write(ipdmp,*)' length  = ',length
         write(ipdmp,*)' next    = ',next
         write(ipdmp,*)' last    = ',last
      end if

      pointer = next
      next = next + length

      if (next.gt.last) then
         write(ipdmp,*)' Error: GETBUF'
         write(ipdmp,*)' ran off end of buffer ',buffer
         write(ipdmp,*)' pointer = ',pointer
         write(ipdmp,*)' length  = ',length
         write(ipdmp,*)' next    = ',next
         write(ipdmp,*)' last    = ',last
         ier=99
         return
      end if

      return
      end

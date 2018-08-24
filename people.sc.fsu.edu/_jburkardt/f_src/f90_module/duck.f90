module duck
!
!  You can't declare the type of function GRINK here.
!
! integer :: grink

  contains

  subroutine flink ( toop )
    implicit none
!
!  You can't declare the type of function GRINK here!
!
    integer ink
    integer toop
    ink = 7
    call plink ( ink )
    toop = 1 + ink + grink ( ink )
    return
  end subroutine flink

  subroutine plink ( ink )
    integer ink
    ink = ink + 3
    return
  end subroutine plink

  function grink ( ink )
!
!  This, apparently, is the only place where the idiotic
!  compiler allows you to declare the type of GRINK.
!
    integer grink
    integer ink
    grink = ink + 3
    return
  end function grink

end module duck

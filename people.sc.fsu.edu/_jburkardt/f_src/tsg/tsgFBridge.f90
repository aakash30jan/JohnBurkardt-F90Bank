subroutine tsgFBridgeGetPoints ( nump, numd, bridge_points )

  use tsgFGrid, only : tsgFAliasGetPoints

  integer numd
  integer nump

  double precision bridge_points(nump*numd)
        
! print *,'Bridge nump = ', nump
        
  call tsgFAliasGetPoints ( nump, numd, bridge_points )

  return
end subroutine tsgFBridgeGetPoints

subroutine tsgFBridgeGetNeededPoints ( nump, numd, bridge_points )

  use tsgFGrid, only : tsgFAliasGetNeededPoints

  integer numd
  integer nump

  double precision bridge_points(nump*numd)
        
! print *,'Bridge nump = ', nump
        
  call tsgFAliasGetNeededPoints ( nump, numd, bridge_points )

  return
end subroutine tsgFBridgeGetNeededPoints

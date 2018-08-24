module tsgFGrid

  public tsgFMakeLocalPolynomialGrid, &
    tsgDeleteGrid, &
    tsgFGetNumPoints, &
    tsgFGetNumDimensions, &
    tsgFGetNumOutputs, &
    tsgFGetPoints, &
    tsgFGetNeededPoints, &
    tsgFAliasGetPoints, &
    tsgFAliasGetNeededPoints, &
    tsgFLoadNeededValues

  private
    double precision, pointer :: tsgf_module_poitns(:,:)
    double precision, pointer :: tsgf_module_newpoitns(:,:)
      
  contains
      
    subroutine tsgFMakeLocalPolynomialGrid ( dimensions, outputs, &
      depth, order, boundary )
      integer :: dimensions, outputs, depth, order, boundary
        
      call tsgFCMakeLocalPolynomialGrid ( dimensions, outputs, depth, &
        order, boundary )

    end subroutine tsgFMakeLocalPolynomialGrid
      
    function tsgFGetNumPoints ( ) result ( n )
      integer :: n
      call tsgFCgetNumPoints ( n )
!     print *,' tsgFGetNumPoints() n= ',n
    end function tsgFGetNumPoints
      
    function tsgFGetNumDimensions ( ) result ( d )
      integer :: d
      call tsgFCgetNumDimensions ( d )
    end function tsgFGetNumDimensions
      
    function tsgFGetNumOutputs ( ) result ( o )
      integer :: o
      call tsgFCgetNumOutputs ( o )
    end function tsgFGetNumOutputs
        
    subroutine tsgFGetPoints ( p )
      double precision, pointer :: p(:,:)
!     integer :: nump, numd
        
      call tsgFCgetPoints ( )
        
      p => tsgf_module_poitns
        
    end subroutine tsgFGetPoints
      
    subroutine tsgFGetNeededPoints ( p )
      double precision, pointer :: p(:,:)
!     integer :: nump, numd
        
      call tsgFCgetNeededPoints ( )
        
      p => tsgf_module_newpoitns
        
    end subroutine tsgFGetNeededPoints
      
    subroutine tsgFAliasGetPoints ( nump, numd, bridge_points )
      integer :: nump, numd
      double precision, target :: bridge_points(nump,numd)
        
      tsgf_module_poitns => bridge_points(1:nump,1:numd)
        
    end subroutine tsgFAliasGetPoints
      
    subroutine tsgFAliasGetNeededPoints ( nump, numd, bridge_points )
      integer :: nump, numd
      double precision, target :: bridge_points(nump,numd)
        
      tsgf_module_newpoitns => bridge_points(1:nump,1:numd)
        
    end subroutine tsgFAliasGetNeededPoints
      
    subroutine tsgFLoadNeededValues ( needed_points )
     double precision :: needed_points(:)
        
     call tsgFCLoadNeededPoints ( needed_points )
        
    end subroutine tsgFLoadNeededValues

!   subroutine tsgDeleteGrid ( )
!     call c_tsgDeleteGrid ( )
!   end subroutine tsgDeleteGrid

end module tsgFGrid

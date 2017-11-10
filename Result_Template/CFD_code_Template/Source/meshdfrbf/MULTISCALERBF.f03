! F2003
! --------- Multiscale RBF mesh deformation module ---------
!   All subroutines and functions for implementing
!    multiscale rbf mesh deformation (see test program)
!  
!   Developed using gfortran v4.7 & v4.9.2
!
!  LKEDWARD 2016



module MULTISCALERBF
  use MATRIX
  use UTILS
  implicit none
  

  contains

    !
    ! FUNCTION: Simple norm fcn.
    ! 
    function norm(x, y, z) result(n)
      real(dp), intent(in) :: x, y, z
      real(dp) :: n
      
      n = 0
      n = sqrt(x**2 + y**2 + z**2)
    
      return
    end function norm
    ! -------------------------------





    !
    ! FUNCTION: Anisotropic norm fcn.
    ! 
    function anorm(x, y, z, normx, normy, normz, zparam) result(n)
      real(dp), intent(in) :: x, y, z, normx, normy, normz, zparam
      real(dp) :: n
      
      real(dp) :: dot, r, a!, b 

      r = sqrt(x**2 + y**2 + z**2)
      r = max(r,mat_tol)
      dot = (x*normx + y*normy + z*normz)/r

      a = zparam*count([dot<0])
      n = (1+a*dot**2)*r
      
      return

    end function anorm
    ! -------------------------------



    !
    ! SUBROUTINE:  Full RBF: Control-point influence matrix
    !
    !  Generates influence matrices for RBF interpolation
    !
    subroutine fullrbf_Mmat(X, Y, Z, radius, Mmat, Pmat)
        real(dp), intent(in) :: X(:), Y(:), Z(:)
        real(dp), intent(in) :: radius
        real(dp), intent(out), allocatable :: Mmat(:,:)
        real(dp), intent(out), allocatable, optional :: Pmat(:,:)
        
        integer :: i, j, nST
        real(dp) :: normx
        real(dp) :: e
        
        nST = size(X,1)
        
        allocate(Mmat(nST,nST))
        if (present(Pmat)) then
          allocate(Pmat(4,nST))
          Pmat(1,:) = 1
          Pmat(2,:) = X(:)
          Pmat(3,:) = Y(:)
          Pmat(4,:) = Z(:)
        end if

        ! --- Populate upper triangle influences
        do j=1, nST, 1
            do i=j, nST, 1
                
                normx = norm(X(i) - X(j), Y(i) - Y(j), Z(i) - Z(j))
                e = normx/radius
                
                if (e >= 1) then
                    Mmat(i,j) = 0;
                else
                    Mmat(i,j) = ((1-e)**4)*(4*e+1)  ! Wendland C2
                end if
                
            end do
        end do

        ! --- Fill lower triangle by symmetry
        do j=1, nST, 1
          do i=1,j-1, 1
            Mmat(i,j) = Mmat(j,i)
          end do
        end do

        return
        
    end subroutine fullrbf_Mmat
    ! ------------------------------------------------------



    !
    ! SUBROUTINE: Full RBF Solve
    !  
    !  Solves a 3D RBF system using 
    !   - approximate LU decomposition for (optional) polynomial basis
    !   - Cholesky decomposition for RBF basis
    !
    subroutine fullrbf_solve(Mmat, Pmat, dX, dY, dZ, ax, ay, az, cond)
      real(dp), intent(in) :: Mmat(:,:)
      real(dp), intent(in), optional :: Pmat(:,:)
      real(dp), intent(in) :: dX(:), dY(:), dZ(:)
      real(dp), intent(out) :: ax(:), ay(:), az(:)
      real(dp), intent(out), optional :: cond

      real(dp) :: M_LU(size(Mmat,1),size(Mmat,1))
      real(dp), dimension(size(Mmat,1)) :: dXres, dYres, dZres
      integer :: N
      integer :: offset

      N = size(Mmat,1)
      
      dXres = dX
      dYres = dY
      dZres = dZ


      if (present(Pmat)) then
        call LU_decomp_solve(transpose(Pmat),dXres,ax(1:4))
        call LU_decomp_solve(transpose(Pmat),dYres,ay(1:4))
        call LU_decomp_solve(transpose(Pmat),dZres,az(1:4))
        dXres = dXres - matmul(transpose(Pmat),ax(1:4))
        dYres = dYres - matmul(transpose(Pmat),ay(1:4))
        dZres = dZres - matmul(transpose(Pmat),az(1:4))
        offset = 5
      else
        offset = 1
      end if

      if (present(cond)) then
        call CH_decomp(Mmat,M_LU,cond)
      else
        call CH_decomp(Mmat,M_LU)
      end if

      call CH_solve(M_LU, dXres, ax(offset:))
      call CH_solve(M_LU, dYres, ay(offset:))
      call CH_solve(M_LU, dZres, az(offset:))

      return

    end subroutine fullrbf_solve
    ! -----------------------------------------------------



    !
    ! SUBROUTINE: Full RBF Transfer
    !
    !  Operation-intensive transfer of RBF field
    !  Detects use of polynomial terms based on dimension of ax/ay/az coefficients
    !
    subroutine fullrbf_transfer(Xc, Yc, Zc, X, Y, Z, radius, ax, ay, az, dX, dY, dZ)
      real(dp), intent(in) :: Xc(:), Yc(:), Zc(:)
      real(dp), intent(in) :: X(:), Y(:), Z(:)
      real(dp), intent(in) :: radius
      real(dp), intent(in) :: ax(:), ay(:), az(:)
      real(dp), intent(out) :: dX(:), dY(:), dZ(:)

      integer :: i,j, Ncp, N, offset
      real(dp) :: coef, e, normx

      Ncp = size(Xc,1)
      N = size (X,1)

      dX(:) = 0
      dY(:) = 0
      dZ(:) = 0
      do i=1, N, 1
        
        ! If polynomial is included or not
        if (size(ax,1) > Ncp) then
          dX(i) = ax(1) + ax(2)*X(i) + ax(3)*Y(i) + ax(4)*Z(i)
          dY(i) = ay(1) + ay(2)*X(i) + ay(3)*Y(i) + ay(4)*Z(i)
          dZ(i) = az(1) + az(2)*X(i) + az(3)*Y(i) + az(4)*Z(i)
          offset = 4
        else
          offset = 0
        end if

        do j=1, Ncp, 1

          normx = norm(X(i) - Xc(j), Y(i) - Yc(j), Z(i) - Zc(j))
          e = normx/radius
          
          if (e < 1) then
            coef = ((1-e)**4)*(4*e+1)  ! Wendland C2
            dX(i) = dX(i) + coef*ax(j+offset)
            dY(i) = dY(i) + coef*ay(j+offset)
            dZ(i) = dZ(i) + coef*az(j+offset)
          end if
            
        end do

      end do
        
      return


    end subroutine fullrbf_transfer
    ! -------------------------------------------------------



    
    !
    ! SUBROUTINE: multiscale
    !  Processes control point mesh (surface points) to produce
    !  interpolation matrices phi_b, psi_r & L of:
    !
    !   |u_b| = |phi_b  0 | * |alpha|
    !   |u_r|   |psi_r  L |   |beta |
    ! 
    !  Also outputs arrays of new indices and support radii
    !
    subroutine multiscale(Xv, Yv, Zv, nBase, r0, phi_b, psi_r, LCRS, inew, radii, normsX, normsY, normsZ, scalings)
      real(dp), intent(in) :: Xv(:), Yv(:), Zv(:)  ! Control point coordinates
      integer, intent(in) :: nBase          ! No. of base level points
      real(dp), intent(in) :: r0
      real(dp), intent(out), allocatable :: phi_b(:,:), psi_r(:,:)
      type(CRSMAT), intent(out) :: LCRS     ! Refinement point interp. matrix in CRS format
      integer, intent(out), allocatable :: inew(:)
      real(dp), intent(out), allocatable :: radii(:)
      real(dp), intent(in) :: normsX(:), normsY(:), normsZ(:)
      real(dp), intent(out), allocatable :: scalings(:)

      ! ------ Scalars
      integer :: i, j, k, l, p
      integer :: Ncp
      integer :: nActive, nInactive
      real(dp) :: normx, minRad, coef, e, minSepDist, dot
      integer :: iMax, iLoc
      logical :: space_filling
      integer(dp) :: Noper

      ! ------ Arrays
      real(dp), allocatable :: sepDist(:), sepDistTemp(:), sepDistLoc(:), siblings(:), rho(:), phi_temp(:,:)
      integer, allocatable ::  activeList(:), inactive(:), parent(:)
      integer, allocatable :: col_ind_temp(:)
      real(dp), allocatable :: vals_temp(:)
      
      minSepDist = 0.3

      ! ---- CMD ARG: base selection mode
      space_filling = .True.
      
      Ncp = size(Xv, 1) ! Total no. of control points

      ! ------ anorm scalings
      allocate(sepDist(Ncp),scalings(Ncp))
      scalings = 1.d0
      do i=1,nCp
        sepDist(i) = 1000000
        do j=1,nCp

          dot = normsX(i)*normsX(j) + normsY(i)*normsY(j) + normsZ(i)*normsZ(j)
          !write(*,*) dot
          if (dot < 0) then
            normx = norm(Xv(i) - Xv(j), Yv(i) - Yv(j), Zv(i) - Zv(j))
            if (normx < sepDist(i)) then
              sepDist(i) = normx
              e = dot
            end if
          end if

        end do
        scalings(i) = 1.d0/sepDist(i)
        !write(*,*) sepDist(i), scalings(i), scalings(i)*aparam
        !read(*,*)
      end do


      ! ------ Initiaise
      allocate(phi_temp(Ncp, nBase))
      allocate(siblings(Ncp))
      allocate(parent(Ncp), sepDistTemp(Ncp), sepDistLoc(Ncp), rho(nbase))
      phi_temp(:,:) = 0
      sepDist(:) = 1000000
      sepDistLoc(:) = 1000000

      allocate(activeList(Ncp))
      allocate(inew(Ncp))
      activeList(1) = 1
      nActive = 1
      parent(:) = 1

      allocate(inactive(Ncp))
      i=1
      inactive = [(i, i=2, Ncp)]
      nInactive = (Ncp-1)
      sepDist(1) = -1
      sepDistLoc(1) = -1

      allocate(radii(Ncp))
      radii(1) = r0
      

      allocate(col_ind_temp(int(Ncp,dp)*int(Ncp,dp)/int(nBase,dp)))
      allocate(vals_temp(int(Ncp,dp)*int(Ncp,dp)/int(nBase,dp)))
      allocate(LCRS%row_ptr(Ncp+1-nBase))
      
      
      ! Loop until all points processed
      Noper = 0
      iMax = 1
      p = 1
      LCRS%row_ptr(1) = 1
      do while (nActive < Ncp)

        ! --- Loop through remaining inactive points
        k = activeList(nActive)   ! K = index of last pointed added to active list

        do j=1, nInactive
          l = inactive(j)           ! L = current point from inactive list
          
          ! Update separation distance of point L, based on last added point K
          normx = anorm(Xv(k) - Xv(l), Yv(k) - Yv(l), Zv(k) - Zv(l), normsX(l), normsY(l), normsZ(l), scalings(l))
          if (normx < sepDist(l)) then
            sepDist(l) = normx
            parent(l) = k                    ! ####
          end if


          ! Calculate pair influence coefficient
          normx = norm(Xv(l) - Xv(k), Yv(l) - Yv(k), Zv(l) - Zv(k))
          if (normx < sepDistLoc(l)) then
            sepDistLoc(l) = normx
          end if
          if (nActive > nBase) then
            normx = anorm(Xv(l) - Xv(k), Yv(l) - Yv(k), Zv(l) - Zv(k), normsX(k), normsY(k), normsZ(k), scalings(k))
          end if
          e = normx/radii(nActive)
          if (e >= 1) then
            coef = 0
          else
            coef = ((1-e)**4)*(4*e+1) ! Wendland C2
          end if

          ! phi_b & psi_r matrices
          ! (rows using old indices)
          if (nActive <= nBase) then
            phi_temp(l, nActive) = coef
          end if

          ! CRS(M): Compressed row storage of control-point influence matrix
          !   Save indices where influence exists
          !   (Only done for refinement points)
          if ((e < 1).AND.(nActive > nBase)) then
            col_ind_temp(p) = l
            vals_temp(p) = coef
            p = p + 1
          end if

          Noper = Noper + 1
        end do


        if (nActive <= nBase) then
          ! Self- influences in phi_b matrix (rows using old indices)
          phi_temp(k, nActive) = 1.0
        else
          ! Self-influence in LCRS (col_ind using old indices)
          col_ind_temp(p) = k
          vals_temp(p) = 1.0
          p = p + 1
          !CRS(L): Save row pointer indices
          LCRS%row_ptr(nActive-nBase+1) = p
        end if

        if (nActive < nBase)  then
          if (space_filling .eqv. .false.) then
            sepDistTemp = parent
            where (sepDistLoc < minSepDist) sepDistTemp = -1
            if (maxval(sepDistTemp) > 0) then
              iMax = maxloc( [(count(sepDistTemp==activeList(i)), i=1,nActive)],1 )
              sepDistTemp(:) = sepDistLoc(:)
              where (parent/=activeList(iMax)) sepDistTemp = -1
              iMax = maxloc(sepDistTemp,1 )
              parent(iMax) = iMax
            else
              iMax = maxloc(sepDist,1)
            end if
            
          else
            iMax = maxloc(sepDist,1) !#### sepDistLoc
          end if
          !write(*,*) sepDist(iMax)
        else
          ! Find point with max separation distance
          iMax = maxloc(sepDist,1) !#### sepDistLoc
        end if

        Noper = Noper + int(nInactive,dp)

        ! Save to active list
        activeList(nActive+1) = iMax

        ! Save sep dist as local RBF radius
        if (nActive < nBase) then
          radii(nActive+1) = r0
        else
          minRad = sepDist(iMax)  !#### sepDistLoc
          if (minRad <= mat_tol) then
            minRad = radii(nActive)
          end if
          radii(nActive+1) = minRad
        end if
        
        
        ! Increment global counter
        nActive = nActive + 1
        
        ! Remove from inactive list & set sep dist criterion to zero
        iLoc = minloc(abs(inactive - iMax), 1) ! Find index in inactive
        inactive(iLoc) = inactive(nInactive) ! (Put last element of inactive into the removed slot)
        nInactive = nInactive - 1

        sepDist(iMax) = -1
        sepDistLoc(iMax) = -1


      end do ! All points processed

      ! Last point (only influences self)
      col_ind_temp(p) = iMax
      vals_temp(p) = 1.0
      p = p + 1
      LCRS%row_ptr(Ncp+1-nBase) = p
      
      ! --- Save L CRS format
      inew = activeList
      do i=1,p-1
        col_ind_temp(i) = minloc(abs(iNew - col_ind_temp(i)),1) ! Convert old indices to new ones
      end do
      allocate(LCRS%col_ind(p-1))
      allocate(LCRS%vals(p-1))
      LCRS%col_ind(1:p-1) = col_ind_temp(1:p-1)
      LCRS%vals(1:p-1) = vals_temp(1:p-1)
      deallocate(col_ind_temp, vals_temp)

      ! Reorder temporary influence matrix (since rows stored using old indices)
      phi_temp(1:Ncp,:) = phi_temp(activeList,:)
      scalings = scalings(activeList)

      ! Split into individual matrices (base & refinement)
      allocate(psi_r(Ncp-nBase,nBase))
      allocate(phi_b(nBase,nBase))
      psi_r(1:Ncp-nBase,:) = phi_temp(nBase+1:,:)
      phi_b(1:nBase,:) = phi_temp(1:nBase,:)

      ! Fill in phi_b by symmetry
      do i=1,nBase
        do j=i,nBase
          !if (i==j) then
          !  phi_b(i,j) = 1.0
          !else
            coef = max(phi_b(i,j),phi_b(j,i))
            phi_b(i,j) = coef
            phi_b(j,i) = coef
          !end if
        end do
      end do

      return
      
    end subroutine multiscale
    ! -----------------------------------------------------------------



    !
    ! SUBROUTINE: Multiscale RBF Solve
    !  Solves for coefficent vector alpha_x/y/z
    !
    subroutine multiscale_solve(nBase, DXs, DYs, DZs, phi_b, psi_r, LCRS, &
                                ax, ay, az)
      integer, intent(in) :: nBase
      real(dp), intent(in) :: DXs(:), DYs(:), DZs(:)
      real(dp), intent(in) :: phi_b(:,:)
      real(dp), intent(in) :: psi_r(:,:)
      type(CRSMAT), intent(in) :: LCRS
      real(dp), intent(out) , allocatable:: ax(:), ay(:), az(:)
      
      
      ! ------- Scalars
      integer :: i, j, p, Ncp
      real(dp) :: coef, cond
      integer(dp) :: Noper

      ! ------- Arrays
      real(dp), allocatable :: ax_temp(:), ay_temp(:), az_temp(:)
      real(dp), allocatable :: dXres(:), dYres(:), dZres(:)
      real(dp), allocatable :: phi_b_LU(:,:)
      integer, allocatable :: ipiv(:)

      Ncp = size(DXs, 1) ! Total no. of control points
      
      ! ------ Initialise arrays
      allocate(ax(Ncp))
      allocate(ay(Ncp))
      allocate(az(Ncp))
      ax(:) = 0
      ay(:) = 0
      az(:) = 0
      allocate(dXres(Ncp))
      allocate(dYres(Ncp))
      allocate(dZres(Ncp))
      dXres(:) = DXs(:)
      dYres(:) = DYs(:)
      dZres(:) = DZs(:)
      
      
      ! ------- Solve base level
      !write(*,*) 'Solving base level...'
      allocate(ax_temp(nBase))
      allocate(ay_temp(nBase))
      allocate(az_temp(nBase))
      allocate(phi_b_LU(nBase,nBase))
      allocate(ipiv(nBase))

      Noper = int(nbase,dp)*int(nbase,dp)*int(nbase,dp)
      !call tic('Solving base level', t)
      call LU_decomp(phi_b,phi_b_LU,ipiv, cond)
      call LU_solve(phi_b_LU, ipiv, dXres(1:nBase), ax_temp)
      call LU_solve(phi_b_LU, ipiv, dYres(1:nBase), ay_temp)
      call LU_solve(phi_b_LU, ipiv, dZres(1:nBase), az_temp)
      !call toc(t)
      !call print_metric('BASE_COND_LOG ', dbl=log10(cond))

      ! Update residual for remaining points
      dXres(nBase+1:) = dXres(nBase+1:) - matmul(psi_r, ax_temp)
      dYres(nBase+1:) = dYres(nBase+1:) - matmul(psi_r, ay_temp)
      dZres(nBase+1:) = dZres(nBase+1:) - matmul(psi_r, az_temp)
      
      ax(1:nBase) = ax_temp
      ay(1:nBase) = ay_temp
      az(1:nBase) = az_temp
      deallocate(ax_temp, ay_temp, az_temp)
      
      ! ------ Loop through remaining (localised) points
      !call tic('Solving refinement points', t)
      do i=(nBase+1),Ncp
        
        ! Update localised coefficients
        ax(i) = dXres(i)
        ay(i) = dYres(i)
        az(i) = dZres(i)

        ! Update residual for remaining points
        do j = LCRS%row_ptr(i-nBase),LCRS%row_ptr(i-nBase+1)-1
          p = LCRS%col_ind(j)
          coef = LCRS%vals(j)
          dXres(p) = dXres(p) - coef*ax(i)
          dYres(p) = dYres(p) - coef*ay(i)
          dZres(p) = dZres(p) - coef*az(i)
          Noper = Noper + 1
        end do
        
      end do
      !call toc(t)

      !call print_metric('Nop_solve',dint = Noper)

      !ax(nBase+1:) = 0

      ! ------ Calculate total residuals
      !write(*,*) 'Multiscale residuals...'
      dXres(:) = DXs(:) - multiscale_eval(phi_b, psi_r, LCRS, ax) 
      dYres(:) = DYs(:) - multiscale_eval(phi_b, psi_r, LCRS, ay) 
      dZres(:) = DZs(:) - multiscale_eval(phi_b, psi_r, LCRS, az) 

      !write(*,*) 'RMS error (x) = ', sqrt(dot_product(dXres(:), dXres(:))/Ncp)
      !write(*,*) 'RMS error (y) = ', sqrt(dot_product(dYres(:), dYres(:))/Ncp)
      !write(*,*) 'RMS error (z) = ', sqrt(dot_product(dZres(:), dZres(:))/Ncp)

      !write(*,*) 'Error(x) (RMS/MAX/MAXLOC)= '
      !write(*,*) sqrt(dot_product(dXres(:), dXres(:))), maxval(abs(dXres(:))), maxloc(abs(dXres(:)))
      !write(*,*) 'Error(y) (RMS/MAX/MAXLOC)= '
      !write(*,*) sqrt(dot_product(dYres(:), dYres(:))), maxval(abs(dYres(:))), maxloc(abs(dYres(:)))
      !write(*,*) 'Error(z) (RMS/MAX/MAXLOC)= '
      !write(*,*) sqrt(dot_product(dZres(:), dZres(:))), maxval(abs(dZres(:))), maxloc(abs(dZres(:)))

      return
      
    end subroutine multiscale_solve
    ! --------------------------------------------------------



    !
    ! FUNCTION: evaluates a set of alpha values on CP mesh points
    !  (For calculating residuals)
    !
    function multiscale_eval(phi_b, psi_r, LCRS, alpha) result(X)
      real(dp), intent(in) :: phi_b(:,:)
      real(dp), intent(in) :: psi_r(:,:)
      type(CRSMAT), intent(in) :: LCRS
      real(dp), intent(in) :: alpha(:)
      real(dp) :: X(size(alpha,1))

      integer :: i, j, p, nBase, nRef, nCP
      real(dp) :: coef

      nBase = size(phi_b,1)
      nRef = size(psi_r,1)
      nCP = nBase + nRef


      X(:) = 0
      do i=1,nCp
        if (i <= nBase) then
          do j=1,nCP
            if (j <=nBase) then
              X(j) = X(j) + phi_b(j,i)*alpha(i)
            else
              X(j) = X(j) + psi_r(j-nBase,i)*alpha(i)
            end if
          end do
        else
          do j = LCRS%row_ptr(i-nBase),LCRS%row_ptr(i-nBase+1)-1
            p = LCRS%col_ind(j)
            coef = LCRS%vals(j)
            X(p) = X(p) + coef*alpha(i)
          end do
        end if
      end do

    end function multiscale_eval
    ! --------------------------------------------------------------



    !
    ! SUBROUTINE: Preprocess target mesh
    !  Generates influence matrix for volume mesh in CRS format
    ! 
    ! Uses wall distances to accelerate pre-processing
    !
    subroutine preproc_vol(Xs, Ys, Zs, Xv, Yv, Zv, wallDist, nBase, & 
                       r0, radii, psi_v, Nvolcp, normsX, normsY, normsZ, scalings)
      real(dp), intent(in) :: Xs(:), Ys(:), Zs(:)
      real(dp), intent(in) :: Xv(:), Yv(:), Zv(:)
      real(dp), intent(in) :: wallDist(:)
      integer, intent(in) :: nBase
      real(dp), intent(in) :: r0
      real(dp), intent(in) :: radii(:)
      type(CRSMAT), intent(out) :: psi_v
      integer, intent(out) :: Nvolcp(size(Xv,1))
      real(dp), intent(in) :: normsX(:), normsY(:), normsZ(:)
      real(dp), intent(in) :: scalings(:)

      ! 1GB Limit on 4-byte integer array length
      integer :: MAX_ARRAY_LENGTH = 250000000
      integer :: debug_n1, debug_n2, debug_n3

      integer :: i, j, k
      integer :: Ncp, Nt
      integer, allocatable :: col_ind_temp(:)
      !real(dp), allocatable :: vals_temp(:)
      real(dp) :: normx, e
      integer(dp) :: Noper

      
      !write(*,*) 'Preprocessing target mesh...'
      Ncp = size(Xs,1)
      Nt = size(Xv,1)
      
      allocate(col_ind_temp(MAX_ARRAY_LENGTH))
      !allocate(vals_temp(MAX_ARRAY_LENGTH))
      allocate(psi_v%row_ptr(Nt+1))
      
      psi_v%row_ptr(1) = 1
      Noper = 0
      j = 1
      debug_n1 = 0
      debug_n2 = 0
      debug_n3 = 0

      do i=1,Nt
        
        ! Find index in radii
        !CPiMax(i) = maxloc(merge(1,0,wallDist(i)<radii)*[(j,j=1,size(radii))],1)
        if (wallDist(i) > r0) then
          col_ind_temp(j) = 0
          j = j + 1
          !debug_n1 = debug_n1 + 1
          Noper = Noper + 1
        elseif (wallDist(i) > radii(nBase+1)) then
          col_ind_temp(j) = 1
          j = j+1
          Noper = Noper + 1
          !debug_n2 = debug_n2 + 1
        else
          col_ind_temp(j) = 1
          j = j+1
          
          do k=nBase+1, Ncp
            !if (blockTags(k)/=block) cycle   ! ###

            !dx = [Xv(i) - Xs(k), Yv(i) - Ys(k), Zv(i) - Zs(k)]
            !dot = (dx(1)*surfNorms(k,1) + dx(2)*surfNorms(k,2) + dx(3)*surfNorms(k,3))/sqrt(dx(1)**2+dx(2)**2+dx(3)**2)
            !dx = (1+exp(-10*dot))*dx
            !normx = norm(dx(1), dx(2), dx(3))
            !normx =  norm(Xv(i) - Xs(k), Yv(i) - Ys(k), Zv(i) - Zs(k))
            normx = anorm(Xv(i) - Xs(k), Yv(i) - Ys(k), Zv(i) - Zs(k), normsX(k), normsY(k), normsZ(k), scalings(k)) ! ####
            e = normx/radii(k)
            !if (e >= 1) then
            !  coef = 0
            !else
            !  coef = ((1-e)**4)*(4*e+1) ! Wendland C2
            !end if
            if (e < 1) then
              col_ind_temp(j) = k+1
              !vals_temp(j) = coef
              j = j+1
            end if
            Noper = Noper + 1
          end do
          !debug_n3 = debug_n3 + 1
        end if
        
        ! --- Save no. of active CPs for plotting
        Nvolcp(i) = (nBase+j-1-psi_v%row_ptr(i))
        psi_v%row_ptr(i+1) = j
      end do
      !!$omp end parallel

      allocate(psi_v%col_ind(j-1))
      !allocate(psi_v%vals(j-1))
      psi_v%col_ind(1:j-1) = col_ind_temp(1:j-1)
      !psi_v%vals(1:j-1) = vals_temp(1:j-1)
      deallocate(col_ind_temp)!, vals_temp)
      
      !write(*,*) 'N1,N2,N3 : ',debug_n1, debug_n2, debug_n3

      !memCost = (j-1+Nt+1)*4
      
      !call print_metric('Nop_ppvol',dint=Noper)
      !call print_metric('Mem_ppvol',dbl=(j-1+Nt+1)*4.d0)

      return
      !write(*,*) 'Done (NOper = ', Noper,')'
      !write(*,*) 'nnz = ', j-1, row_ptr(Nt+1)-1
      !write(*,*) ' Memory Cost :~ ', (j-1+Nt+1)*4/1000, 'kB'
      
    end subroutine preproc_vol
    ! -------------------------------------------------------------------------



    !
    ! SUBROUTINE: Multiscale RBF transfer
    ! (Transfer control point displacements to target mesh)
    !
    subroutine multiscale_transfer(Xs, Ys, Zs, Xv, Yv, Zv, radii, psi_v, nBase, ax, ay, az, DXv, DYv, DZv, &
                                                      normsX, normsY, normsZ, scalings)
      real(dp), intent(in) :: Xs(:), Ys(:), Zs(:)
      real(dp), intent(in) :: Xv(:), Yv(:), Zv(:)
      type(CRSMAT) :: psi_v
      real(dp), intent(in) :: radii(:)
      integer, intent(in) :: nBase
      real(dp), intent(in) :: ax(:), ay(:), az(:)
      real(dp), intent(out), allocatable :: DXv(:), DYv(:), DZv(:)
      real(dp), intent(in) :: normsX(:), normsY(:), normsZ(:)
      real(dp), intent(in) :: scalings(:)

      ! --- VARS
      integer :: i, j, k, Ns, Nv
      real(dp) :: normx, e, coef
      integer(dp) :: Noper

      Ns = size(Xs,1)
      Nv = size(Xv,1)
      
      allocate(DXv(Nv), DYv(Nv), DZv(Nv))

      DXv(:) = 0
      DYv(:) = 0
      DZv(:) = 0

      ! --- Transfer
      !write(*,*) 'Transfering displacements...'
      Noper = 1
      do i = 1, Nv, 1
        
        !write(*,*) i
        do k=psi_v%row_ptr(i), (psi_v%row_ptr(i+1)-1), 1
          
          if (psi_v%col_ind(k)==0) then
            ! do nothing
          elseif (psi_v%col_ind(k)==1) then
            do j=1, nBase, 1
              
              normx = norm(Xv(i) - Xs(j), Yv(i) - Ys(j), Zv(i) - Zs(j))
              !normx = anorm(Xs(j)-Xv(i), Ys(j)-Yv(i), Zs(j)-Zv(i), [surfNorms(j,1), surfNorms(j,2), surfNorms(j,3)]) ! ####
              e = normx/radii(j)

              if (e <=1) then
                coef = ((1-e)**4)*(4*e+1) ! Wendland C2
                DXv(i) = DXv(i) + coef*ax(j)
                DYv(i) = DYv(i) + coef*ay(j)
                DZv(i) = DZv(i) + coef*az(j)
              end if
              Noper = Noper + 1

            end do
          else
            j = psi_v%col_ind(k)-1
            
            !dx = [Xv(i) - Xs(j), Yv(i) - Ys(j), Zv(i) - Zs(j)]
            !dot = (dx(1)*surfNorms(j,1) + dx(2)*surfNorms(j,2) + dx(3)*surfNorms(j,3))/sqrt(dx(1)**2+dx(2)**2+dx(3)**2)
            !dx = (1+exp(-10*dot))*dx         ! dx = sqrt(1-dot+dot**2)*dx
            !normx = norm(dx(1), dx(2), dx(3))
            !normx = norm(Xv(i) - Xs(j), Yv(i) - Ys(j), Zv(i) - Zs(j))
            normx = anorm(Xv(i) - Xs(j), Yv(i) - Ys(j), Zv(i) - Zs(j), normsX(j), normsY(j), normsZ(j),scalings(j)) ! ####
            e = normx/radii(j)

            if (e <=1) then
              coef = ((1-e)**4)*(4*e+1) ! Wendland C2
              DXv(i) = DXv(i) + coef*ax(j)
              DYv(i) = DYv(i) + coef*ay(j)
              DZv(i) = DZv(i) + coef*az(j)
            end if
            Noper = Noper + 1
          end if
        end do
        
      end do
      
      !call print_metric('Nop_update',dint=Noper)
      !write(*,*) 'Noper = ', Noper, '( ', Ns*Nv/(1.0*Noper), ' reduction)'

      return

    end subroutine multiscale_transfer
    ! -------------------------------------------------------------------


    !
    ! SUBROUTINE: Greedy full point selection algorithm
    !
    subroutine greedy_full(X, Y, Z, DX, DY, DZ, nBase, r0, iBase, ax, ay, az, rCor)
      real(dp), intent(in) :: X(:), Y(:), Z(:)  ! Surface coordinates
      real(dp), intent(in) :: DX(:), DY(:), DZ(:)
      integer, intent(in) :: nBase          ! No. of base level points
      real(dp), intent(in) :: r0
      integer, intent(out), allocatable :: iBase(:)
      real(dp), intent(out), allocatable, optional :: ax(:), ay(:), az(:)
      real(dp), intent(out), optional :: rCor

      integer :: i, j, iMax
      integer(dp) :: Noper
      integer :: Ninitial = 1
      integer :: nActive, step
      real(dp), allocatable :: ax1(:), ay1(:), az1(:)
      real(dp), allocatable :: PMat(:,:), Mmat(:,:), ri(:)
      real(dp) :: cond, normx, rmax, ravg
      real(dp), allocatable :: dx_eval(:), dy_eval(:), dz_eval(:)
      real(dp) :: Ex(size(X,1)),  Ey(size(X,1)), Ez(size(X,1)), Ei(size(X,1))
      real(dp), dimension(nBase) :: Xac, Yac, Zac

      !write(*,*) 'Starting greedy selection...'

      ! --- Initial point selection ---
      Noper = 0
      i=1
      allocate(iBase(nBase))
      allocate(ax1(nBase), ay1(nBase), az1(nBase))
      allocate(dx_eval(size(X,1)), dy_eval(size(X,1)),  dz_eval(size(X,1)))
      step = size(X,1)/Ninitial
      nActive = size([(i, i=1,size(X,1),step)],1)
      iBase(1:nActive) = [(i, i=1,size(X,1),step)]
      Xac(:nActive) = X(iBase(:nActive))
      Yac(:nActive) = Y(iBase(:nActive))
      Zac(:nActive) = Z(iBase(:nActive))
      call fullrbf_Mmat(X, Y, Z, r0, Pmat, Mmat)

      do while (nActive < nBase)
        
        !write(*,*) nActive

        call fullrbf_solve(Mmat(iBase(:nActive),iBase(:nActive)), &
          dX=DX(iBase(:nActive)), dY=DY(iBase(:nActive)), dZ=DZ(iBase(:nActive)), &
          ax=ax1(:nActive), ay=ay1(:nActive), az=az1(:nActive))

        call fullrbf_transfer(Xac(:nActive), Yac(:nActive), Zac(:nActive), &
                      X, Y, Z, r0, ax1(:nActive), ay1(:nActive), az1(:nActive), dx_eval, dy_eval, dz_eval)

        Ex = DX - dx_eval
        Ey = DY - dy_eval
        Ez = DZ - dz_eval
        Ei = sqrt(Ex*Ex + Ey*Ey + Ez*Ez)

        iMax = maxloc(Ei,1)
        iBase(nActive+1) = iMax
        Xac(nActive+1) = X(iMax)
        Yac(nActive+1) = Y(iMax)
        Zac(nActive+1) = Z(iMax)
        nActive = nActive + 1

        !deallocate(ax1, ay1, az1)
        Noper = Noper + int(nActive,dp)*int(nActive,dp)*int(nActive,dp) &
          + int(nActive,dp)*(size(X,1,dp)-int(nActive,dp))
      end do

      if (present(rCor)) then

        rmax = -1000
        ravg = 0
        allocate(ri(nBase))

        do i=1,nBase
          ri(i) = 1000
          do j=1,nBase
            if (i==j) then
              cycle
            end if
            normx = norm(Xac(j)-Xac(i),Yac(j)-Yac(i),Zac(j)-Zac(i))
            ri(i) = min(ri(i),normx)
          end do

          rmax = max(rmax,ri(i))
          ravg = ravg + ri(i)
        end do
        ravg = ravg/nBase
        rCor = ravg
      end if

      if (present(ax).AND.present(ay).AND.present(az)) then
        allocate(ax(nBase))
        allocate(ay(nBase))
        allocate(az(nBase))
        call fullrbf_Mmat(X(iBase), Y(iBase), Z(iBase), r0, Pmat, Mmat)
        call fullrbf_solve(Mmat, dX=DX(iBase), dY=DY(iBase), dZ=DZ(iBase), &
                  ax=ax, ay=ay, az=az, cond=cond)

        write(*,*) 'RBF System conditioning: ', cond
      end if

      !call print_metric('Nop_ppsrf', dint=Noper)

      return


    end subroutine greedy_full
    ! ----------------------------------------------------------------


    !
    ! SUBROUTINE: Volume correction preprocessing
    !
    subroutine volcor_preproc(Xs, Ys, Zs, Xv, Yv, Zv, Rcor, vi)
      real(dp), intent(in) :: Xs(:), Ys(:), Zs(:)
      real(dp), intent(in) :: Xv(:), Yv(:), Zv(:)
      real(dp), intent(in) :: Rcor
      integer, intent(out) :: vi(:)

      integer :: i, j
      real(dp) :: normx, e
      real(dp), allocatable :: dist(:)

      allocate(dist(size(Xv,1)))

      dist(:) = 1000000
      vi(:) = -1
      do i=1,size(Xv,1)
        do j=1,size(Xs,1)

          normx = norm(Xv(i) - Xs(j), Yv(i) - Ys(j), Zv(i) - Zs(j))
          e = normx/Rcor

          if (normx<dist(i)) then
            dist(i) = normx
            if (e<=1) then
              vi(i) = j
            end if
          end if

        end do
      end do

      !call print_metric('Nop_ppvol',dint=size(Xv,1,dp)*size(Xs,1,dp))

      return

    end subroutine volcor_preproc
    ! ----------------------------------------------------------------


    !
    ! SUBROUTINE: Volume correction
    !
    subroutine volcorrect(Xv, Yv, Zv, vi, Xv2, Yv2, Zv2, Xs, Ys, Zs, DDXs, DDYs, DDZs, Rcor)
      real(dp), intent(in) :: Xv(:), Yv(:), Zv(:)
      real(dp), intent(inout) :: Xv2(:), Yv2(:), Zv2(:)
      integer, intent(in) :: vi(:)
      real(dp), intent(in) :: Xs(:), Ys(:), Zs(:)
      real(dp), intent(in) :: DDXs(:), DDYs(:), DDZs(:)
      real(dp), intent(in) :: Rcor

      integer :: i,j
      real(dp) :: normx, e, coef

      do i=1,size(Xv,1)
        j = vi(i)
        if (j>0) then
          normx = norm(Xv(i) - Xs(j), Yv(i) - Ys(j), Zv(i) - Zs(j))
          e = normx/Rcor
          coef = ((1-e)**4)*(4*e+1) ! Wendland C2
          Xv2(i) = Xv2(i) + coef*DDXs(j)
          Yv2(i) = Yv2(i) + coef*DDYs(j)
          Zv2(i) = Zv2(i) + coef*DDZs(j)

        end if

      end do

      return

    end subroutine volcorrect
    ! -----------------------------------------------------------------





end module MULTISCALERBF

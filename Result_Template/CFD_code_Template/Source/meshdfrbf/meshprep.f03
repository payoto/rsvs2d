! F2003
! --------- Multiscale rbf mesh deformation for SU2 (PREPROC) ---------
!   Multiscale mesh deformation pre-processor
!
!   Developed using gfortran v5.4.0
!
!  LKEDWARD 2017

program meshprep
  use MULTISCALERBF
  use MATRIX, only: CRSMAT, dp
  use UTILS
  use FILEIO
  implicit none


  ! --------- Config vars ---------
  integer :: rbfmode
  character(100) :: meshtype, volfile, surffile
  real(dp) :: r0
  integer :: nbase
  character(100) :: outputFile


  ! --------- General vars ---------
  logical :: fExist
  integer :: nvol, nsurf, fh
  integer, allocatable :: inew(:), Nvolcp(:)
  real(dp), allocatable, dimension(:) :: Xv, Yv, Zv, Xs, Ys, Zs, zeronorms, zerowallDs(:)
  real(dp), allocatable :: radii(:), phi_b(:,:), psi_r(:,:), scalings(:)
  real(dp), allocatable :: Mmat(:,:)
  type(CRSMAT) :: LCRS, psi_v
  type(SU2MESH) :: meshdata
  

  ! --------- Read Config. Parameters ---------
  call confFileOpen('meshdef.conf',fh)
  call confScanString(fh,'mesh',volfile)
  call confScanString(fh,'voltype',meshtype)   ! Supported meshtypes: xyz, su2
  call confScanString(fh,'surfpts',surffile)
  !if (trim(upperstr(meshtype)) == 'SU2') call confScanString(fh,'surftag',surfTag)
  call confScanInt(fh,'rbfmode',rbfmode,1)
  call confScanReal(fh,'r0',r0)
  if (rbfmode == 1) call confScanInt(fh,'nbase',nbase)
  close(fh)


  ! --------- Check input files exist ---------
  inquire(file=trim(volfile), exist=fExist)
  if (.NOT.fExist) then
    write(*,*) "(!) Volume points file: '", trim(volfile), "' not found."
    stop
  end if
  inquire(file=trim(surffile), exist=fExist)
  if (.NOT.fExist) then
    write(*,*) "(!) Surface points file: '", trim(surffile), "' not found."
    stop
  end if


  ! --------- Load volume mesh data ---------
  if (trim(upperstr(meshtype)) == 'XYZ' ) then
    ! --- Mesh type: XYZ
    call txtload_xyz(volfile,nVol, Xv, Yv, Zv)

  elseif (trim(upperstr(meshtype)) == 'SU2') then
    ! --- Mesh type: SU2
    call txtload_su2(volfile, meshdata)
    nVol = size(meshdata%XYZ,1)
    allocate(Xv(nVol))
    allocate(Yv(nVol))
    allocate(Zv(nVol))
    Xv = meshdata%XYZ(:,1)
    Yv = meshdata%XYZ(:,2)
    if (meshdata%NDIME > 2) then
      Zv = meshdata%XYZ(:,3)
    else
      Zv(:) = 0
    end if
  else
    write(*,*) "(!) Mesh type '",trim(meshtype),"' is not supported"
    stop
  end if


  ! --------- Load surface mesh data ---------
  call txtload_xyz(surffile,Nsurf,Xs,Ys,Zs)
  allocate(zeronorms(nsurf))
  zeronorms(:) = 0
  allocate(Nvolcp(nvol))
  allocate(zerowallDs(nvol))
  zerowallDs(:) = 0

  write(*,*) 'Nvol=', nvol
  write(*,*) 'Nsurf=', nsurf


  ! --------- Perform preprocessing ---------
  if (rbfmode == 1) then
    ! --- Multiscale preprocessing ---
    write(*,*) "Multiscale preprocessing..."
    call multiscale(Xs, Ys, Zs, nbase, r0, phi_b, psi_r, LCRS, inew, radii, zeronorms, zeronorms, zeronorms, scalings)
    Xs(:) = Xs(inew)
    Ys(:) = Ys(inew)
    Zs(:) = Zs(inew)

    call preproc_vol(Xs, Ys, Zs, Xv, Yv, Zv, zerowallDs, nBase, & 
                           r0, radii, psi_v, Nvolcp, zeronorms, zeronorms, zeronorms, scalings)

  elseif (rbfmode == 0) then
    ! --- Full RBF method ---
    write(*,*) "Generating full RBF interpolation matrix..."
    call fullrbf_Mmat(Xs, Ys, Zs, r0, Mmat)

  else
    write(*,*) "(!) Unrecognised rbfmode '",rbfmode,"' "
    write(*,*) "    0 = full rbf"
    write(*,*) "    1 = multiscale"
    stop

  end if


  ! --------- Save preprocessor file ---------
  outputFile = trim(volfile)//'.meshdef'
  call bin_open(fh, outputfile)

  ! --- (1) Save mesh data ---
  call binWrite_realvec(fh,Xs)
  call binWrite_realVec(fh,Ys)
  call binWrite_realVec(fh,Zs)
  if (trim(upperstr(meshtype)) == 'XYZ' ) then
    call binWrite_int(fh,0)        ! Mesh save type xyz
    call binWrite_realVec(fh,Xv)
    call binWrite_realVec(fh,Yv)
    call binWrite_realVec(fh,Zv)
  elseif (trim(upperstr(meshtype)) == 'SU2' ) then
    call binWrite_int(fh,1)        ! Mesh save type su2
    call binWrite_SU2MESH(fh,meshdata)
  end if

  ! --- (2) Save preprocessor data ---
  call binWrite_int(fh,rbfmode)
  call binWrite_real(fh,r0)
  if (rbfmode == 1) then
    write(*,*) 'Saving preprocessed (multiscale) mesh data to: ', trim(outputFile), ' ...'
    call binWrite_int(fh,nbase)
    call binWrite_realMat(fh,phi_b)
    call binWrite_realMat(fh,psi_r)
    call binWrite_realCRS(fh,LCRS)
    call binWrite_intVec(fh,inew)
    call binWrite_realVec(fh,radii)
    call binWrite_realVec(fh,scalings)
    call binWrite_realCRS(fh,psi_v)
  elseif (rbfmode == 0) then
    write(*,*) 'Saving preprocessed (full rbf) mesh data to: ', trim(outputFile), ' ...'
    call binWrite_realMat(fh,Mmat)
  end if

  close(fh)
  write(*,*) ' Done.'


end program meshprep

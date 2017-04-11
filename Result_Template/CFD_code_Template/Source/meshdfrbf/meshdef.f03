! F2003
! --------- Multiscale rbf mesh deformation for SU2 (DEFORMER) ---------
!   Multiscale mesh deformation
!
!   Developed using gfortran v5.4.0
!
!  LKEDWARD 2017

program meshdef
  use omp_lib
  use MULTISCALERBF
  use MATRIX, only: dp, CRSMAT
  use UTILS
  use FILEIO
  implicit none


  ! --------- Config vars ---------
  integer :: rbfmode
  character(100) :: preprocfile, deffile, outfile
  real(dp) :: r0
  integer :: nbase, stat

  ! --------- General Vars ---------
  logical :: fExist, zeroMode
  integer :: i, nvol, nsurf, fh, meshType
  integer, allocatable :: inew(:)
  real(dp), allocatable, dimension(:) :: Xv, Yv, Zv, Xs, Ys, Zs, zeronorms, scalings(:)
  real(dp),  allocatable, dimension(:) :: DXs, DYs, DZs, ax, ay, az
  real(dp), allocatable, dimension(:) :: DXv, DYv, DZv
  real(dp), allocatable :: radii(:), phi_b(:,:), psi_r(:,:), Mmat(:,:)
  type(CRSMAT) :: LCRS, psi_v
  type(SU2MESH) ::meshdata


  ! --------- Load Configuration ---------
  !  .meshdef preprocessor file is first command-line arg
  !  surface displacements file is second
  !  output file is third
  call GET_COMMAND_ARGUMENT(1, preprocfile,status=stat)
  if (stat /= 0) then
    write(*,*) "(!) Expecting .meshdef preprocessor filename as first argument"
    stop
  end if
  inquire(file=trim(preprocfile), exist=fExist)
  if (.NOT.fExist) then
    write(*,*) "(!) Can't find .meshdef preprocessor file (",trim(preprocfile),")"
    stop
  end if
  
  call GET_COMMAND_ARGUMENT(2, deffile,status=stat)
  if (stat /= 0) then
    write(*,*) "(!) No surface displacements filename supplied (second argument)"
    write(*,*) "    Zero-displacement mode: outputting original mesh"
    zeroMode = .True.
  else
    inquire(file=trim(preprocfile), exist=fExist)
    if (.NOT.fExist) then
      write(*,*) "(!) Can't find .meshdef preprocessor file (",trim(preprocfile),")"
      stop
    end if
    zeroMode = .False.
  end if
  call GET_COMMAND_ARGUMENT(3, outfile,status=stat)


  ! ------------ Load preprocessed data ---------
  call bin_open(fh, preprocfile)

  ! --- (1) Load mesh data ---
  call binread_realVec(fh,Xs)
  call binread_realVec(fh,Ys)
  call binread_realVec(fh,Zs)
  nsurf = size(Xs,1)

  call binread_int(fh,meshType)
  if (meshType == 0) then
    call binRead_realVec(fh,Xv)
    call binRead_realVec(fh,Yv)
    call binRead_realVec(fh,Zv)
    nvol = size(Xv,1)
  elseif (meshType == 1) then
    call binRead_SU2MESH(fh,meshdata)
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
  end if
  
  ! --- (2) Load preprocessor data ---
  call binRead_int(fh,rbfmode)
  call binRead_real(fh,r0)

  if (rbfmode == 1) then

    ! Multiscale data
    allocate(zeronorms(nsurf))
    zeronorms(:) = 0

    call binread_int(fh,nbase)

    call binread_realMat(fh,phi_b)
    call binread_realMat(fh,psi_r)
    call binread_realCRS(fh,LCRS)
    call binread_intVec(fh,inew)
    call binread_realVec(fh,radii)
    call binread_realVec(fh,scalings)
    call binread_realCRS(fh,psi_v)
    
  elseif (rbfmode == 0) then
    ! Full rbf data
    call binread_realMat(fh,Mmat)
  end if
  close(fh)


  ! --------- Perform volume mesh deformation ---------
  if (.Not.zeroMode) then

    ! --- Read deflection file ---
    ! ASCII
    ! Columns of DX DY DZ xNSURF
    allocate(DXs(nsurf))
    allocate(DYs(nsurf))
    allocate(DZs(nsurf))
    call txtload_xyz(deffile, i, DXs, DYs, DZs)
    
    ! --- Solve & transfer system ---
    allocate(ax(nsurf))
    allocate(ay(nsurf))
    allocate(az(nsurf))
    if (rbfmode == 1) then
      ! Multiscale rbfs
      DXs(:) = DXs(inew)
      DYs(:) = DYs(inew)
      DZs(:) = DZs(inew)
      call multiscale_solve(nBase, DXs, DYs, DZs, phi_b, psi_r, LCRS, &
                                    ax, ay, az)

      call multiscale_transfer(Xs, Ys, Zs, Xv, Yv, Zv, radii, &
                psi_v, nBase, ax, ay, az, DXv, DYv, DZv, zeronorms, zeronorms, zeronorms, scalings)

    elseif (rbfmode == 0) then
      ! Conventional rbfs
      call fullrbf_solve(Mmat, dx=DXs, dY=DYs, dZ=DZs, ax=ax, ay=ay, az=az)
      allocate(DXv(nvol))
      allocate(DYv(nvol))
      allocate(DZv(nvol))
      call fullrbf_transfer(Xs,Ys,Zs,Xv,Yv,Zv,r0,ax,ay,az,DXv,DYv,DZv)
    end if

  else

    ! --- Zero-displacement mode: output original undeformed mesh ---
    allocate(DXv(nvol))
    allocate(DYv(nvol))
    allocate(DZv(nvol))
    DXv(:) = 0
    DYv(:) = 0
    DZv(:) = 0

  end if

   
  ! --------- Deform & write to file ---------
  Xv = Xv + DXv
  Yv = Yv + DYv
  Zv = Zv + DZv
  if (stat/=0) then
    write(*,*) nvol
    do i=1,nvol
      write(*,*) Xv(i), Yv(i), Zv(i)
    end do
  else
    if (meshType==0) then
      open(unit=newunit(fh),file=outfile)
      write(fh,*) nvol
      do i=1,nvol
        write(fh,*) Xv(i), Yv(i), Zv(i)
      end do
      close(fh)
    else
      meshdata%XYZ(:,1) = Xv
      meshdata%XYZ(:,2) = Yv
      meshdata%XYZ(:,3) = Zv
      call txtwrite_su2(outfile,meshdata)
    end if
  end if



end program meshdef

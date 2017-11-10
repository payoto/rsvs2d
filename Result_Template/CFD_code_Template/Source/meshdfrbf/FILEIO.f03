! F2003
! --------- FILE I/O Module ---------
!
!   Assorted subroutines for reading/writing files
!
! Contains:
!  - newunit         - generates new integer for file handle
!  - txtload_xyz     - read 3D point cloud data from text file
!  - txtload_SU2     - read unstructured volume mesh from SU2 file
!  - txtwrite_SU2    - write unstructured volume mesh to SU2 file
!  
!  - binary data stream subroutines (see below for more info.)
!
!  L.KEDWARD 2017



module FILEIO
  use MATRIX, only : dp, CRSMAT
  implicit none
  public
  
  type SU2MESH
    integer :: NDIME
    integer :: NELEM
    integer, allocatable :: CONNECTIVITY(:,:)
    integer :: NPOIN
    real(dp), allocatable :: XYZ(:,:)
    integer, allocatable :: PINDICES(:)   ! Global point indices
    integer :: NMARK
    character(100), allocatable :: MARKER_TAGS(:)
    integer, allocatable :: MARKER_NELEMS(:)
    integer, allocatable :: MARKER_CONNECTIVITY(:,:,:)
  end type SU2MESH


  contains
  
  
  ! 
  ! FUNCTION: New file unit helper function
  ! 
  !  (future note: not necessary in F08)
  !
  integer function newunit(unit)
    integer, intent(out), optional :: unit
    ! local
    integer, parameter :: LUN_MIN=10, LUN_MAX=1000
    logical :: opened
    integer :: lun
    ! begin
    newunit=-1
    do lun=LUN_MIN,LUN_MAX
        inquire(unit=lun,opened=opened)
        if (.not. opened) then
            newunit=lun
            exit
        end if
    end do
    if (present(unit)) unit=newunit
  end function newunit
  ! -------------------------------------------------
    
  

  ! 
  ! SUBROUTINE: Point cloud (xyz) Loader
  !
  ! Read point cloud files of form: Npts followed by x,y,z columns
  !
  subroutine txtload_xyz(filename, Npts, X, Y, Z)
    character(*), intent(in) :: filename
    integer, intent(out) :: Npts
    real(dp), intent(out), allocatable :: X(:), Y(:), Z(:)

    real(dp), allocatable :: A(:,:)
    integer :: fh, i

    open(unit=newunit(fh),file=filename)
    read(fh,*) Npts
    allocate(A(Npts,3))
    read(fh,*) (A(i,:),i=1,Npts)
    close(fh)

    allocate(X(Npts))
    allocate(Y(Npts))
    allocate(Z(Npts))
    X(:) = A(:,1)
    Y(:) = A(:,2)
    Z(:) = A(:,3)
    deallocate(A)

    return

  end subroutine txtload_xyz



  !
  ! SUBROUTINE: SU2 Mesh Loader
  !  Loads contents of SU2 mesh text file into SU2MESH structure (see above)
  !
  ! 
  subroutine txtload_su2(filename, meshdata)
    character(*), intent(in) :: filename
    type(SU2MESH), intent(out) :: meshdata

    integer :: fh, eType, i1, i, j, mark_nelem_max
    character(100) :: linestring

    open(unit=newunit(fh), file=filename, status="old")

    ! --- Read no. of dimensions ---
    read(fh,'(A)') lineString
    i1 = index(lineString, '=')
    read(lineString(i1+1:),*) meshdata%NDIME

    ! --- Read no. of cells ---
    read(fh,'(A)') lineString
    i1 = index(lineString, '=')
    read(lineString(i1+1:),*) meshdata%NELEM

    ! --- Read cell connectivity ---
    if (meshdata%NDIME < 3) then
      allocate(meshdata%CONNECTIVITY(meshdata%NELEM,6))
    else
      allocate(meshdata%CONNECTIVITY(meshdata%NELEM,10))
    end if
    do i=1,meshdata%NELEM
      read(fh,'(A)') linestring
      read(lineString,*) eType
      if (eType == 5) then
        read(linestring,*) meshdata%CONNECTIVITY(i,1:5) ! Triangle type
      elseif (eType == 9 .OR. eType == 10) then
        read(linestring,*) meshdata%CONNECTIVITY(i,1:6) ! Quad & tet type
      elseif (eType == 12) then
        read(linestring,*) meshdata%CONNECTIVITY(i,1:10) ! Hexhedral type
      elseif (eType == 13) then
        read(linestring,*) meshdata%CONNECTIVITY(i,1:8) ! Wedge type
      elseif (eType == 14) then
        read(linestring,*) meshdata%CONNECTIVITY(i,1:7) ! Pyramid type
      else
        write(*,*) "(!) txtload_su2: Unsupported element type: ", eType
        write(*,*) " Element index ", i
        stop
      end if
    end do
    

    ! --- Read no. of points ---
    read(fh,'(A)') lineString
    i1 = index(lineString, '=')
    read(lineString(i1+1:),*) meshdata%NPOIN

    ! --- Read point coords ---
    allocate(meshdata%XYZ(meshdata%NPOIN,meshdata%NDIME))
    allocate(meshdata%PINDICES(meshdata%NPOIN))
    do i=1,meshdata%NPOIN
      read(fh,*) meshdata%XYZ(i,:), meshdata%PINDICES(i)
    end do

    ! --- Read no. of markers ---
    read(fh,'(A)') lineString
    i1 = index(lineString, '=')
    read(lineString(i1+1:),*) meshdata%NMARK

    ! --- Boundary marker data: first pass ---
    allocate(meshdata%MARKER_TAGS(meshdata%NMARK))
    allocate(meshdata%MARKER_NELEMS(meshdata%NMARK))
    mark_nelem_max = 0
    do i=1,meshdata%NMARK

      ! -- Get marker tag
      read(fh,'(A)') linestring
      i1 = index(lineString, '=')
      read(lineString(i1+1:),*) meshdata%MARKER_TAGS(i)

      ! -- Get no. of marker elements
      read(fh,'(A)') linestring
      i1 = index(lineString, '=')
      read(lineString(i1+1:),*) meshdata%MARKER_NELEMS(i)

      ! --- First pass: skip through marker connectivity ---
      do j = 1, meshdata%MARKER_NELEMS(i)
        read(fh,*)
      end do

    end do
    mark_nelem_max = maxval(meshdata%MARKER_NELEMS)

    ! --- Boundary marker data: second pass ---
    if (meshdata%NDIME < 3) then
      allocate(meshdata%MARKER_CONNECTIVITY(meshdata%NMARK,mark_nelem_max,3))
    else
      allocate(meshdata%MARKER_CONNECTIVITY(meshdata%NMARK,mark_nelem_max,5))
    end if
    rewind(fh)
    do i=1,2+meshdata%NELEM+1+meshdata%NPOIN+1
      read(fh,*)
    end do
    do i=1,meshdata%NMARK
      read(fh,*)
      read(fh,*)
      do j = 1, meshdata%MARKER_NELEMS(i)
        read(fh,'(A)') linestring
        read(linestring,*) eType
        if (eType == 3) then
          read(linestring,*) meshdata%MARKER_CONNECTIVITY(i,j,1:3) ! Line type
        elseif (eType == 5) then
          read(linestring,*) meshdata%MARKER_CONNECTIVITY(i,j,1:4) ! Triangle type
        elseif (eType == 9) then
          read(linestring,*) meshdata%MARKER_CONNECTIVITY(i,j,1:5) ! Quad type
        else
          write(*,*) "(!) txtload_su2: Unsupported boundary element type: ", eType
          stop
        end if
      end do
    end do

    close(fh)

    return

  end subroutine txtload_su2
  ! --------------------------------------------------------------------



  !
  ! SUBROUTINE: SU2 Writer
  !  Takes an SU2MESH structure (see above) and writes to .su2 file
  !
  !
  subroutine txtwrite_su2(filename, meshdata)
    character(*), intent(in) :: filename
    type(SU2MESH), intent(in) :: meshdata

    integer :: fh, i, j, eType

    open(unit=newunit(fh), file=filename, status='unknown')
    
    write(fh,'(A,I7)') 'NDIME=', meshdata%NDIME
    write(fh,'(A,I7)') 'NELEM=', meshdata%NELEM

    do i=1,meshdata%NELEM
      eType = meshdata%CONNECTIVITY(i,1)
      if (eType == 5) then
        write(fh,'(I7,I7,I7,I7,I7)') meshdata%CONNECTIVITY(i,1:5) ! Triangle type
      elseif (eType == 9 .OR. eType == 10) then
        write(fh,'(I7,I7,I7,I7,I7,I7)') meshdata%CONNECTIVITY(i,1:6) ! Quad & tet type
      elseif (eType == 12) then
        write(fh,'(I7,I7,I7,I7,I7,I7,I7,I7,I7,I7)') meshdata%CONNECTIVITY(i,1:10) ! Hexhedral type
      elseif (eType == 13) then
        write(fh,'(I7,I7,I7,I7,I7,I7,I7,I7)') meshdata%CONNECTIVITY(i,1:8) ! Wedge type
      elseif (eType == 14) then
        write(fh,'(I7,I7,I7,I7,I7,I7,I7)') meshdata%CONNECTIVITY(i,1:7) ! Pyramid type
      else
        write(*,*) "(!) txtwrite_su2: Unsupported element type: ", eType
        stop
      end if
    end do

    write(fh,'(A,I7)') 'NPOIN=', meshdata%NPOIN

    do i=1,meshdata%NPOIN
      write(fh,*) meshdata%XYZ(i,:), meshdata%PINDICES(i)
    end do

    write(fh,'(A,I7)') 'NMARK=', meshdata%NMARK
    do i=1,meshdata%NMARK

      write(fh,'(A,A)') 'MARKER_TAG= ',meshdata%MARKER_TAGS(i)

      write(fh,'(A,I7)') 'MARKER_TAG= ',meshdata%MARKER_NELEMS(i)

      do j = 1, meshdata%MARKER_NELEMS(i)
        eType = meshdata%MARKER_CONNECTIVITY(i,j,1)
        if (eType == 3) then
          write(fh,'(I7,I7,I7)') meshdata%MARKER_CONNECTIVITY(i,j,1:3) ! Line type
        elseif (eType == 5) then
          write(fh,'(I7,I7,I7,I7)') meshdata%MARKER_CONNECTIVITY(i,j,1:4) ! Triangle type
        elseif (eType == 9) then
          write(fh,'(I7,I7,I7,I7,I7)') meshdata%MARKER_CONNECTIVITY(i,j,1:5) ! Quad type
        else
          write(*,*) "(!) txtwrite_su2: Unsupported boundary element type: ", eType
          stop
        end if
      end do

    end do

    close(fh)
    return


  end subroutine txtwrite_su2
  ! --------------------------------------------------------------------


 
  
  
  
  
  !
  ! -------------- Generic binary stream file format ---------------
  ! 
  ! Format for storing arbitrary data consisting of variables of
  !  various dimension, length and type.
  ! Sequential data - variables are written and accessed in same order
  !
  ! 'Write(..)' uses machine definition of endianness and
  !  local compiler definition of double precision reals and integers,
  !  so transfer of data between machines is not guaranteed to work.
  !  IE. Preprocess on the same machine as main run
  ! (Future note: F08 can fix the type definitions)
  ! 
  ! SUBROUTINES:
  !   bin_open(fh,filename) - Open file for binary stream access (r|w)
  !   binWrite_real         - Read/write real (dp) scalar data
  !   binRead_real
  !   binWrite_int          - Read/write integer scalar data
  !   binRead_int
  !   binWrite_realVec      - Read/write real (dp) vector data
  !   binRead_realVec
  !   binWrite_intVec       - Read/write integer vector data
  !   binRead_intVec
  !   binWrite_realMat      - Read/write real (dp) matrix data
  !   binRead_realMat
  !   binWrite_realCRS      - Read/write compressed-row storage matrix
  !   binRead_realCRS
  !
  !  L.Kedward 2017
  !
  subroutine bin_open(fh,filename)
    integer, intent(inout) :: fh
    character(*), intent(in) :: filename

    open(unit=newunit(fh), file=filename, access='stream')

  end subroutine bin_open
  
  ! --- Real Scalar ---
  subroutine binWrite_real(fh, realScalar)
    integer, intent(in) :: fh
    real(dp), intent(in) :: realScalar

    integer :: int32

    int32 = 1; write(fh) int32     ! vartype = double (1)
    int32 = 0; write(fh) int32     ! ndim = 0
    int32 = 1; write(fh) int32     ! vector size (1=scalar)
    write(fh) realScalar

  end subroutine binWrite_real
  
  subroutine binRead_real(fh, realScalar)
    integer, intent(in) :: fh
    real(dp), intent(out) :: realScalar

    integer :: int32

    read(fh) int32                 ! vartype = double (1)
    if (int32 /= 1) then
      write(*,*) "(!) Binary read error: wrong type, expecting real scalar, got type=",int32
      stop
    end if
    
    read(fh) int32                  ! ndim = 0
    if (int32 /= 0) then
      write(*,*) "(!) Binary read error: wrong no. dimensions, expecting real scalar, got ndim=",int32
      stop
    end if

    read(fh) int32                ! vector size (1=scalar)
    if (int32 /= 1) then
      write(*,*) "(!) Binary read error: wrong dimension size, expecting real scalar, got dim=",int32
      stop
    end if
    
    read(fh) realScalar

  end subroutine binRead_real
  ! --------------------------------------------------------------------


  ! --- Int scalar ---
  subroutine binWrite_int(fh, intScalar)
    integer, intent(in) :: fh
    integer, intent(in) :: intScalar

    integer :: int32

    int32 = 0; write(fh) int32     ! vartype = double (1)
    int32 = 0; write(fh) int32     ! ndim = 1
    int32 = 1; write(fh) int32     ! vector size (1=scalar)
    write(fh) intScalar

  end subroutine binWrite_int

  subroutine binRead_int(fh, intScalar)
    integer, intent(in) :: fh
    integer, intent(out) :: intScalar

    integer :: int32

    read(fh) int32                 ! vartype = int (0)
    if (int32 /= 0) then
      write(*,*) "(!) Binary read error: wrong type, expecting int scalar, got type=",int32
      stop
    end if
    
    read(fh) int32                  ! ndim = 0
    if (int32 /= 0) then
      write(*,*) "(!) Binary read error: wrong no. dimensions, expecting int scalar, got ndim=",int32
      stop
    end if

    read(fh) int32                ! vector size (1=scalar)
    if (int32 /= 1) then
      write(*,*) "(!) Binary read error: wrong dimension size, expecting int scalar, got dim=",int32
      stop
    end if
    
    read(fh) intScalar

  end subroutine binRead_int
  ! --------------------------------------------------------------------


  ! --- Real vector ---
  subroutine binWrite_realvec(fh, realVector)
    integer, intent(in) :: fh
    real(dp), intent(in) :: realVector(:)

    integer :: int32

    int32 = 1; write(fh) int32                   ! vartype = double (1)
    int32 = 1; write(fh) int32                   ! ndim = 1
    int32 = size(realVector,1); write(fh) int32  ! vector size
    write(fh) realVector

  end subroutine binWrite_realvec

  subroutine binRead_realvec(fh, realVector)
    integer, intent(in) :: fh
    real(dp), intent(out), allocatable :: realVector(:)

    integer :: int32

    read(fh) int32                 ! vartype = double (1)
    if (int32 /= 1) then
      write(*,*) "(!) Binary read error: wrong type, expecting real vector, got type=",int32
      stop
    end if
    
    read(fh) int32                  ! ndim = 1
    if (int32 /= 1) then
      write(*,*) "(!) Binary read error: wrong no. dimensions, expecting real vector, got ndim=",int32
      stop
    end if

    read(fh) int32                ! vector size
    allocate(realVector(int32))
    read(fh) realVector(:)

  end subroutine binRead_realvec
  ! --------------------------------------------------------------------


  ! --- Int vector ---
  subroutine binWrite_intVec(fh, intVector)
    integer, intent(in) :: fh
    integer, intent(in) :: intVector(:)

    integer :: int32

    int32 = 0; write(fh) int32                   ! vartype = int (0)
    int32 = 1; write(fh) int32                   ! ndim = 1
    int32 = size(intVector,1); write(fh) int32  ! vector size
    write(fh) intVector

  end subroutine binWrite_intVec

  subroutine binRead_intVec(fh, intVector)
    integer, intent(in) :: fh
    integer, intent(out), allocatable :: intVector(:)

    integer :: int32

    read(fh) int32                 ! vartype = int (0)
    if (int32 /= 0) then
      write(*,*) "(!) Binary read error: wrong type, expecting int vector, got type=",int32
      stop
    end if
    
    read(fh) int32                  ! ndim = 1
    if (int32 /= 1) then
      write(*,*) "(!) Binary read error: wrong no. dimensions, expecting int vector, got ndim=",int32
      stop
    end if

    read(fh) int32                ! vector size
    allocate(intVector(int32))
    read(fh) intVector(:)

  end subroutine binRead_intVec
  ! --------------------------------------------------------------------


  ! --- Real matrix ---
  subroutine binWrite_realMat(fh, realMat)
    integer, intent(in) :: fh
    real(dp), intent(in) :: realMat(:,:)

    integer :: int32,j

    int32 = 1; write(fh) int32                  ! vartype = double (1)
    int32 = 2; write(fh) int32                  ! ndim = 2 (matrix)
    int32 = size(realMat,1); write(fh) int32    ! matrix n (rows)
    int32 = size(realMat,2); write(fh)  int32   ! matrix m (cols)

    do j = 1,size(realMat,2)
      write(fh) realMat(:,j)
    end do

  end subroutine binWrite_realMat

  subroutine binRead_realMat(fh, realMat)
    integer, intent(in) :: fh
    real(dp), intent(out), allocatable :: realMat(:,:)

    integer :: int32,j,n,m

    read(fh) int32                 ! vartype = double (1)
    if (int32 /= 1) then
      write(*,*) "(!) Binary read error: wrong type, expecting real matrix, got type=",int32
      stop
    end if
    
    read(fh) int32                  ! ndim = 2
    if (int32 /= 2) then
      write(*,*) "(!) Binary read error: wrong no. dimensions, expecting real matrix, got ndim=",int32
      stop
    end if

    read(fh) int32                ! matrix size (n rows)
    n = int32
    read(fh) int32                ! matrix size (m cols)
    m = int32
    allocate(realMat(n,m))
    do j=1,m
      read(fh) realMat(:,j)
    end do

  end subroutine binRead_realMat
  ! --------------------------------------------------------------------


  ! --- Real CRS Matrix ---
  subroutine binWrite_realCRS(fh, crsmatrix)
    integer, intent(in) :: fh
    type(CRSMAT), intent(in) :: crsmatrix

    integer :: int32, n1, n2, n3

    int32 = 1; write(fh) int32                          ! vartype = double (1)
    int32 = -1; write(fh) int32                         ! ndim = crs matrix (-1)

    n1 = 0;
    n2 = 0;
    n3 = 0;
    if (allocated(crsmatrix%vals)) n1 = size(crsmatrix%vals,1)
    if (allocated(crsmatrix%col_ind)) n2 = size(crsmatrix%col_ind,1)
    if (allocated(crsmatrix%row_ptr)) n3 = size(crsmatrix%row_ptr,1)

    write(fh) n1     ! size(vals)
    write(fh) n2     ! size(col_ind)
    write(fh) n3     ! size(row_ptr)

    if (n1 > 0) write(fh) crsmatrix%vals(:)
    if (n2 > 0) write(fh) crsmatrix%col_ind(:)
    if (n3 > 0) write(fh) crsmatrix%row_ptr(:)

  end subroutine binWrite_realCRS

  subroutine binRead_realCRS(fh, crsmatrix)
    integer, intent(in) :: fh
    type(CRSMAT), intent(out) :: crsmatrix

    integer :: int32, n1, n2, n3

    read(fh) int32                 ! vartype = double (1)
    if (int32 /= 1) then
      write(*,*) "(!) Binary read error: wrong type, expecting real crs matrix, got type=",int32
      stop
    end if
    
    read(fh) int32                  ! ndim = -1
    if (int32 /= -1) then
      write(*,*) "(!) Binary read error: wrong no. dimensions, expecting real crs matrix type (ndim=-1), got ndim=",int32
      stop
    end if

    read(fh) n1                ! n1 - size(vals)
    read(fh) n2                ! n2 - size(col_ind)
    read(fh) n3                ! n3 - size(row_ptr) = (n+1) dimension size

    if (n1 > 0) allocate(crsmatrix%vals(n1))
    if (n2 > 0) allocate(crsmatrix%col_ind(n2))
    if (n3 > 0) allocate(crsmatrix%row_ptr(n3))

    if (n1 > 0) read(fh) crsmatrix%vals(:)
    if (n2 > 0) read(fh) crsmatrix%col_ind(:)
    if (n3 > 0) read(fh) crsmatrix%row_ptr(:)
    crsmatrix%n = n3-1



  end subroutine binRead_realCRS
  ! --------------------------------------------------------------------
  
  
  ! --- SU2 Mesh Data ---
  subroutine binWrite_SU2MESH(fh, meshdata)
    integer, intent(in) :: fh
    type(SU2MESH), intent(in) :: meshdata

    integer :: j, jj, int32

    int32 = 1; write(fh) int32                          ! vartype = double (1) [redundant for SU2MESH type]
    int32 = -2; write(fh) int32                         ! ndim = SU2 data (-2)
    int32 = meshdata%NDIME; write(fh) int32     ! NDIME
    int32 = meshdata%NELEM; write(fh) int32     ! NELEM
    
    int32 = size(meshdata%CONNECTIVITY,2); write(fh) int32     ! NCols of connectivity
    do j=1,size(meshdata%CONNECTIVITY,2)
      write(fh) meshdata%CONNECTIVITY(:,j)
    end do

    int32 = meshdata%NPOIN; write(fh) int32     ! NPOIN
    write(fh) meshdata%XYZ(:,1)                 ! Xs
    write(fh) meshdata%XYZ(:,2)                 ! Ys
    write(fh) meshdata%XYZ(:,3)                 ! Zs

    int32 = meshdata%NMARK; write(fh) int32     ! NMARK

    do j=1,meshdata%NMARK                       ! Marker tags
      int32 = len(trim(meshdata%MARKER_TAGS(j))); write(fh) int32
      write(fh) trim(meshdata%MARKER_TAGS(j))
    end do

    write(fh) meshdata%MARKER_NELEMS(:)         ! Marker NELEMS
    int32 = size(meshdata%MARKER_CONNECTIVITY,2); write(fh) int32     ! Nrows of marker connectivity
    int32 = size(meshdata%MARKER_CONNECTIVITY,3); write(fh) int32     ! NCols of marker connectivity
    do j=1,meshdata%NMARK                       ! Marker tags
      do jj=1,size(meshdata%MARKER_CONNECTIVITY,3)
        write(fh) meshdata%MARKER_CONNECTIVITY(j,1:meshdata%MARKER_NELEMS(j),jj)
      end do
    end do

  end subroutine binWrite_SU2MESH

  subroutine binRead_SU2MESH(fh, meshdata)
    integer, intent(in) :: fh
    type(SU2MESH), intent(out) :: meshdata

    integer :: j, jj, int32, ncol, nrow
    character(100) :: string

    read(fh) int32                 ! vartype = double (1) [redundant for SU2MESH type]
    
    read(fh) int32                  ! ndim = -2
    if (int32 /= -2) then
      write(*,*) "(!) Binary read error: wrong no. dimensions, expecting SU2MESH type (ndim=-2), got ndim=",int32
      stop
    end if

    read(fh) meshdata%NDIME                  ! NDIME
    read(fh) meshdata%NELEM     ! NELEM
    
    read(fh) ncol     ! NCols of connectivity
    allocate(meshdata%CONNECTIVITY(meshdata%NELEM,ncol))
    do j=1,ncol
      read(fh) meshdata%CONNECTIVITY(:,j)
    end do

    read(fh) meshdata%NPOIN     ! NPOIN
    allocate(meshdata%XYZ(meshdata%NPOIN,3))
    read(fh) meshdata%XYZ(:,1)                 ! Xs
    read(fh) meshdata%XYZ(:,2)                 ! Ys
    read(fh) meshdata%XYZ(:,3)                 ! Zs

    read(fh) meshdata%NMARK     ! NMARK
    allocate(meshdata%MARKER_TAGS(meshdata%NMARK))
    do j=1,meshdata%NMARK                       ! Marker tags
      read(fh) int32
      read(fh) string(1:int32)
      meshdata%MARKER_TAGS(j) = string(1:int32)
    end do

    allocate(meshdata%MARKER_NELEMS(meshdata%NMARK))
    read(fh) meshdata%MARKER_NELEMS(:)         ! Marker NELEMS

    read(fh) nrow     ! Nrow of marker connectivity (marker_nelems max)
    read(fh) ncol     ! NCols of marker connectivity
    allocate(meshdata%MARKER_CONNECTIVITY(meshdata%NMARK,nrow,ncol))
    do j=1,meshdata%NMARK                       ! Marker tags
      do jj=1,ncol
        read(fh) meshdata%MARKER_CONNECTIVITY(j,1:meshdata%MARKER_NELEMS(j),jj)
      end do
    end do

  end subroutine binRead_SU2MESH
  ! --------------------------------------------------------------------
  
end module FILEIO

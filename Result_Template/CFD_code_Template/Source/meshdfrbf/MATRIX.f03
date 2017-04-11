! F2003
! --------- Matrix Utilities Module --------- 
!
!  Various LAPACK wrappers and custom routines
! 
! Contains:
!  - CRSMAT (type)    : Compressed-row storage for sparse matrices
!  - CH_decomp        : Cholesky decomposition for SPD matrices
!  - CH_solve         : Solve linear system using Cholesky decomposition
!  - LU_decomp        : LU decomposition of general matrix
!  - LU_solve         : Solve linear system using LU decomposition
!  - LU_decomp_solve  : LU-decompose and solve linear system
!  - CG_solve         : Solve linear system using conjugate gradient
!  - CG_solve_crs     : CG solve of linear system with A in CRS format
!  
!  - MAT2CRS          : Convert matrix into CRS format
!  - CRS2MAT          : Convert CRS matrix into matrix
!  - CRS_MATVEC       : Multiply CRS matrix with vector
!  - CRS_destroy      : Deallocate CRS matrix object memory
!
! LKEDWARD 2016


module MATRIX
  implicit none
  public

  integer, parameter:: dp=kind(0.d0) ! Machine-compiler double precision
  
  ! Matrix tolerance
  !  Used for conjugate gradient and sparse matrix tolerances
  real(dp) :: mat_tol = 0.000000001d0
  
  
  
  !
  ! TYPE: Compressed row storage matrix
  !
  ! Ref: http://netlib.org/linalg/html_templates/node91.html
  !
  type CRSMAT
    real(dp), allocatable :: vals(:)
    integer, allocatable :: col_ind(:)
    integer, allocatable :: row_ptr(:)
    integer :: n
  end type CRSMAT
  
  
  ! ------------------------- Subroutines ------------------------------
  contains
  
  
  !
  ! SUBROUTINE: Cholesky Decomposition
  !
  !  Wrapper for LAPACK routines to perform a cholesky decomposition
  !   of a symmmetric positive definite matrix
  !
  ! Refs:
  !  http://www.netlib.org/lapack/double/dpotrf.f
  !  http://www.netlib.org/lapack/double/dpocon.f
  !
  subroutine CH_decomp(A,AU, cond)
    real(dp), intent(in) :: A(:,:)
    real(dp), intent(out) :: AU(:,:)
    real(dp), intent(out), optional :: cond

    real(dp) :: work(4*size(A,2))
    integer :: Iwork(size(A,2))
    real(dp) ::ANorm, RCond
    integer :: i, n, LDA, info

    ! External procedures defined in LAPACK
    external DPOTRF
    external DPOCON
    
    ! Store A & B to prevent it from being overwritten by LAPACK
    LDA = size(A,1)
    n = size(A,2)
    AU = A
    
    ! DPOTRF computes the Cholesky factorization of a real symmetric
    ! positive definite matrix A.
    ! 
    ! The factorization has the form
    !    A = U**T * U,  if UPLO = 'U', or
    !    A = L  * L**T,  if UPLO = 'L',
    ! where U is an upper triangular matrix and L is lower triangular.'
    call DPOTRF( 'U', n, AU, LDA, info )

    if (info > 0) then
        write(*,*) 'Matrix is not positive definite (error=',info,')'
        do i=1,5
          write(*,*) A(i,1:5)
        end do
        stop
    end if

    ! DGECON Estimate reciprocal of condition number
    ! 
    if (present(cond)) then
      do i=1,n
        work(i) = sum(A(:,i))
      end do
      ANorm = maxval(work(1:n))
      
      call DPOCON( 'U', n, AU, LDA, ANorm, RCond, work, Iwork, info )
      cond = 1.0_dp/RCond
    end if

  end subroutine CH_decomp
  ! --------------------------------------------------------------------



  !
  ! SUBROUTINE: Cholesky Decomposition Solve
  !  Using existing CH decomposition
  !
  ! http://www.netlib.org/lapack/double/dpotrs.f
  !
  subroutine CH_solve(AU,B,Xout)
    real(dp), intent(in) :: AU(:,:)
    real(dp), intent(in) :: B(:)
    real(dp), intent(out) :: Xout(:)

    real(dp) :: Btemp(size(AU,2))  ! Local copy of B
    integer :: n, LDA, info
    
    ! External procedures defined in LAPACK
    external DPOTRS

    ! Store A & B to prevent it from being overwritten by LAPACK
    LDA = size(AU,1)
    n = size(AU,2)
    Btemp = B
    
    ! DGETRS solve AX=B
    ! (http://www.netlib.no/netlib/lapack/double/dgetrs.f)
    call dpotrs('U', n, 1, AU, LDA, Btemp, n, info)
    
    if (info /= 0) then
      write(*,*) '(!) CH Solve failed (err=', info, ')'
    else
      Xout = Btemp
    end if
    
    return
    
  end subroutine CH_solve
  ! --------------------------------------------------------------------



  
  !
  ! SUBROUTINE: LU Decomposition
  ! (!) Square A only currently
  !
  subroutine LU_decomp(A,ALU,ipivout, cond)
    real(dp), intent(in) :: A(:,:)
    real(dp), intent(out) :: ALU(:,:)
    integer, intent(out), optional :: ipivout(size(A,2))
    real(dp), intent(out), optional :: cond

    real(dp) :: work(4*size(A,2))
    integer :: Iwork(size(A,2))
    integer :: ipiv(size(A,2))
    real(dp) ::ANorm, RCond
    integer :: i, n, LDA, info

    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRS

    ! Store A & B to prevent it from being overwritten by LAPACK
    LDA = size(A,1)
    n = size(A,2)
    ALU = A
    
    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    !write(*,*) 'LU...'
    call DGETRF(LDA, n, ALU, LDA, ipiv, info)
    
    if (present(ipivout)) then
      ipivout = ipiv
    end if
    
    if (info /= 0) then
        stop 'Matrix is numerically singular!'
    end if

    ! DGECON Estimate reciprocal of condition number
    ! 
    if (present(cond)) then
      do i=1,n
        work(i) = sum(A(:,i))
      end do
      ANorm = maxval(work(1:n))
      
      call DGECON( '1', n, ALU, LDA, ANorm, RCond, work, Iwork, info )
      cond = 1.0_dp/RCond
      !write(*,*) '   (RCondition number = ', 1.0_dp/RCond,')'
    end if

  end subroutine LU_decomp
  ! --------------------------------------------------------------------
  

  !
  ! SUBROUTINE: Solve using LU Decomposition
  !  Using existing LU decomposition
  !
  subroutine LU_solve(ALU,ipiv,B,Xout)
    real(dp), intent(in) :: ALU(:,:)
    integer, intent(in) :: ipiv(size(ALU,2))   ! pivot indices
    real(dp), intent(in) :: B(:)
    real(dp), intent(out) :: Xout(:)

    real(dp) :: Btemp(size(ALU,2))  ! Local copy of B
    integer :: n, LDA, info
    
    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRS
    external DGECON

    ! Store A & B to prevent it from being overwritten by LAPACK
    LDA = size(ALU,1)
    n = size(ALU,2)
    Btemp = B
    
    ! DGETRS solve AX=B
    ! (http://www.netlib.no/netlib/lapack/double/dgetrs.f)
    call DGETRS('N', n, 1, ALU, LDA, ipiv, Btemp, n, info)
    
    if (info /= 0) then
      write(*,*) '(!) LU Solve failed (err=', info, ')'
    else
      Xout = Btemp
    end if
    
    return
    
  end subroutine LU_solve
  ! --------------------------------------------------------------------
  

  !
  ! SUBROUTINE: LU Decomposition & Solve
  ! 
  ! 18/11/16 - Modification to allow non-square (needs checking)
  !
  subroutine LU_decomp_solve(A,B,Xout,cond)
    real(dp), intent(in) :: A(:,:)
    real(dp), intent(in) :: B(:)
    real(dp), intent(out) :: Xout(:)
    real(dp), intent(out), optional :: cond
    
    real(dp) :: ALU( size(A,1) , size(A,2) )  ! LU decomp of A
    real(dp) :: Btemp(size(A,1))  ! Local copy of B
    integer, dimension(size(A,2)) :: ipiv   ! pivot indices
    integer :: n, LDA, info
    real(dp) :: ANorm                  ! Norm of A
    real(dp) :: RCond                  ! Reciprocal condition no. of A
    real(dp) :: work(4*size(A,2))
    integer :: Iwork(size(A,2))
    integer :: i
    
    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRS
    external DGECON

    ! Store A & B to prevent it from being overwritten by LAPACK
    LDA = size(A,1)
    n = size(A,2)
    ALU = A
    Btemp = 0
    Btemp(:size(B)) = B
    
    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call DGETRF(LDA, n, ALU, LDA, ipiv, info)
    
    if (info /= 0) then
        write(*,*) 'Matrix is numerically singular! (info=',info,')'
        stop
    end if
    
    
    ! DGETRS solve AX=B
    ! (http://www.netlib.no/netlib/lapack/double/dgetrs.f)
    call DGETRS('N', n, 1, ALU, LDA, ipiv, Btemp, n, info)
    
    if (info /= 0) then
      write(*,*) '(!) LU Solve failed (err=', info, ')'
    else
      Xout = Btemp(:size(Xout))
    end if
    
    
    ! DGECON Estimate reciprocal of condition number
    ! 
    if (present(cond)) then
      do i=1,n
        work(i) = sum(A(:,i))
      end do
      ANorm = maxval(work(1:n))
      
      call DGECON( '1', n, ALU, LDA, ANorm, RCond, work, Iwork, info )
      cond = 1.0_dp/RCond
    end if
    
    
  end subroutine LU_decomp_solve
  ! --------------------------------------------------------------------
  
  
  !
  ! SUBROUTINE: Conjugate Gradient Solve
  ! 
  !
  subroutine CG_solve(A,B,Xout, maxit, Minv)
    real(dp), intent(in) :: A(:,:)
    real(dp), intent(in) :: B(:)
    real(dp), intent(inout) :: Xout(size(A,2))
    integer, intent(in), optional :: maxit
    type(CRSMAT), intent(in), optional :: Minv
    
    ! ---------- Working Vars. ----------
    real(dp), dimension(size(A,2)) :: X
    real(dp), dimension(size(B,1)) :: r, p, Ap, s
    real(dp) :: alpha, r2, r2old, pAp, tol
    integer :: i, maxiter
    
    ! --------- Initialise ----------
    X(:) = Xout
    !X(:) = 0
    
    if (present(maxit)) then
      maxiter = maxit
    else
      maxiter = size(A,1)
    end if
    
    r = B - matmul(A,X)
    if (present(Minv)) then
      p = CRS_MATVEC(Minv, r)
    else
      p = r
    end if
    r2 = dot_product(r, p)
    write(*,*) r2
    tol = maxval([mat_tol, (mat_tol**2)*r2])
    !write(*,*) 'tol:', tol
    
    ! --------- Start  iteration ----------
    do i=1,maxiter
      
      Ap = matmul(A,p)
      pAp = dot_product(p, Ap)
      
      if (pAp < mat_tol) then
        !alpha = 0
        !write(*,*) '(!) zero denom'
        exit
      else
        alpha = r2/pAp
      end if
      
      !if (mod(i,1000) == 0) then
      !    r = B - matmul(A,X)
      !else
          r = r - alpha*Ap
      !end if
        
      X = X + alpha*p
      
      r2old = r2
      if (present(Minv)) then
        s = CRS_MATVEC(Minv, r)
      else
        s = r
      end if
      r2 = dot_product(r, s)
      
      if (r2<tol) then
        exit
      end if
      
      if (r2old > 0.0_dp) then
        p = s + (r2/r2old)*p
      else
        p = s
      end if
      
     !write(*,*) sqrt(r2)
    end do
    
    if (i.eq.(maxiter+1)) then
      write(*,*) 'Hit MaxIt. (r2 = ', r2, ')'
    else
      write(*,*) 'n = ', i, 'r2=', r2
    end if
    Xout = X
    
    
  end subroutine CG_solve
  ! --------------------------------------------------------------------
  
  
  
  
  
  
  !
  ! SUBROUTINE: Conjugate Gradient Solve
  ! (Take in A as CRS format matrix, Ac)
  !
  subroutine CG_solve_crs(Ac,B,Xout, maxit, Minv)
    type(CRSMAT), intent(in) :: Ac
    real(dp), intent(in) :: B(:)
    real(dp), intent(out) :: Xout(Ac%n)
    integer, intent(in) :: maxit
    type(CRSMAT), intent(in), optional :: Minv
    
    ! ---------- Working Vars. ----------
    real(dp), dimension(Ac%n) :: X
    real(dp), dimension(size(B,1)) :: r, p, Ap, s
    real(dp) :: alpha, r2, r2old, pAp, tol
    integer :: i
    
    ! --------- Initialise ----------
    !X(:) = Xout
    X(:) = 0
    
    r = B - CRS_MATVEC(Ac,X)
    if (present(Minv)) then
      p = CRS_MATVEC(Minv, r)
    else
      p = r
    end if
    r2 = dot_product(r, p)
    
    tol = (mat_tol**2)*r2
    !write(*,*) 'tol:', tol
    
    ! --------- Start  iteration ----------
    do i=1,maxit
      
      Ap = CRS_MATVEC(Ac,p)
      pAp = dot_product(p, Ap)
      
      if (pAp < mat_tol) then
        alpha = 0
        write(*,*) '(!) zero denom'
        exit
      else
        alpha = r2/pAp
      end if
      
      if (mod(i,50) == 0) then
          r = B - CRS_MATVEC(Ac,X)
      else
          r = r - alpha*Ap
      end if
        
      X = X + alpha*p
      
      r2old = r2
      if (present(Minv)) then
        s = CRS_MATVEC(Minv, r)
      else
        s = r
      end if
      r2 = dot_product(r, s)
      
      if (r2<tol) then
        exit
      end if
      
      if (r2old > 0.0_dp) then
        p = s + (r2/r2old)*p
      else
        p = s
      end if
      
     !write(*,*) sqrt(r2)
    end do
    
    if (i.eq.(maxit+1)) then
      write(*,*) 'Hit MaxIt. (r2 = ', r2, ')'
    else
      write(*,*) 'n = ', i
    end if
    Xout = X
    
    
  end subroutine CG_solve_crs
  ! --------------------------------------------------------------------

  

  
  
  
  !
  ! SUBROUTINE: Put A into CRS format
  !
  function MAT2CRS(A) result(ACRS)
    real(dp), intent(in) :: A(:,:)
    type(CRSMAT) :: ACRS
    
    integer :: nnz, n, i, j, counter
    
    ! Count number non-zero elements
    n = size(A,1)
    nnz = count(dabs(A) > mat_tol)
    
    ! Initialise CRS format
    allocate(ACRS%vals(nnz))
    allocate(ACRS%col_ind(nnz))
    allocate(ACRS%row_ptr(n+1))
    ACRS%n = n
    
    counter=0
    ACRS%row_ptr(1)=1
    do i=1,n
      do j=1,size(A,2)
        if(dabs(A(i,j)) > mat_tol)then
          counter = counter+1
          ACRS%vals(counter) = A(i,j)
          ACRS%col_ind(counter)=j
        end if
      end do
      
      ! Store row pointer
      ACRS%row_ptr(i+1)=counter+1
        
    end do
    
    ! By convention store nnz+1 here
    ACRS%row_ptr(n+1)=nnz+1
    
  end function MAT2CRS
  ! --------------------------------------------------------------------
  
  
  
  !
  ! SUBROUTINE: Deallocate CRSMAT object
  !
  subroutine CRS_destroy(ACRS)
    type(CRSMAT), intent(inout) :: ACRS
  
    deallocate(ACRS%vals)
    deallocate(ACRS%col_ind)
    deallocate(ACRS%row_ptr)
    
  end subroutine CRS_destroy
  ! --------------------------------------------------------------------
  
  
  
  
  !
  ! SUBROUTINE: Take A out of CRS format
  !
  function CRS2MAT(ACRS) result(A)
    type(CRSMAT), intent(in) :: ACRS
    real(dp), allocatable :: A(:,:)
    
    integer :: counter, counter2, n, i, j
    
    counter=0
    n = ACRS%n
    allocate(A(n,n))
    
    
    do counter2=1,n
      do counter=ACRS%row_ptr(counter2),ACRS%row_ptr(counter2+1)-1
        j = ACRS%col_ind(counter)
        i = counter2
        A(i,j) = ACRS%vals(counter)
      enddo
    enddo
    
  end function CRS2MAT
  ! --------------------------------------------------------------------
  
  
  !
  ! SUBROUTINE: CRS Mat Vec Multiplication
  ! (http://www.netlib.org/linalg/html_templates/node98.html)
  !
  function CRS_MATVEC(ACRS, x) result(Ax)
    type(CRSMAT), intent(in) :: ACRS
    real(dp), intent(in) :: x(:)
    real(dp) :: Ax(size(x,1))
    
    integer :: i, j, n
    
    n = ACRS%n
    
    do i=1,n,1
      Ax(i) = 0
      do j=ACRS%row_ptr(i), ACRS%row_ptr(i+1)-1
        
        Ax(i) = Ax(i) + ACRS%vals(j) * x(ACRS%col_ind(j))
        
      end do
    end do
    
  end function CRS_MATVEC
  ! --------------------------------------------------------------------
  
  
    
end module MATRIX

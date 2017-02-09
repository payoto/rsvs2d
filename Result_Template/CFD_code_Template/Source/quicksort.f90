
PROGRAM TestQuicksort

   IMPLICIT NONE

   ! ------------------------------------------------------------------------------
   ! Interface block(s)
   !
   ! For some reason, the interface block for the quicksort routine is problematic.
   ! If the interface block is contained in a module, -- it hangs at run-time.
   ! ------------------------------------------------------------------------------

   INTERFACE
   RECURSIVE SUBROUTINE Quicksort(Item, First, Last, Indices, kk)
   DOUBLE PRECISION,    DIMENSION(:), INTENT(INOUT) :: Item	! array of values
   INTEGER,               INTENT(IN)    :: First,Last,kk
   INTEGER, DIMENSION(:), INTENT(INOUT) :: Indices
END SUBROUTINE Quicksort
END INTERFACE

! ------------------------------------------------------------------------------
! Local variables
! ------------------------------------------------------------------------------

INTEGER :: GNX
INTEGER,allocatable, DIMENSION(:)	:: indarr
!DOUBLE PRECISION,    DIMENSION(GNX)	:: randvals = 0.0, randcopy = 0.0
double precision,allocatable,DIMENSION(:)  ::randvals,randcopy
INTEGER	:: i,j

print *, "Enter size of arrays:"
read(*,*) GNX

allocate(randvals(GNX))
allocate(randcopy(GNX))
allocate(indarr(GNX))

indarr = (/ (i,i=1,GNX) /)		! indexical array

call RANDOM_NUMBER(randvals)
!randvals=1.0
randcopy = randvals

print *,'   org_arr index'
if (GNX<=100) then
   do i = 1,GNX
      print '(f10.4,1x,i3)',randvals(i),indarr(i)
   enddo
end if

!call Quicksort(randvals,1,GNX,indarr,0)
call KB07AI(randvals,GNX,indarr)
if (GNX<=100) then
   print *,'       sorted  O_index'
   do i = 1,GNX
      print '(i3,1x,f10.4,1x,i3,1x,f10.4)',i,randvals(i),indarr(i),randcopy(indarr(i)) 
   enddo
end if
print *,'sorted - indexed original = ',SUM(abs(randvals - randcopy(indarr)))

END PROGRAM TestQuicksort

   !----------------------------------------------------------------------------
   !
   ! This file is based on the the routine in "Fortran 90 for Engineers & 
   ! Scientists" by Nyhoff and Leestma
   !
   ! Note: In the following subroutines, Item is an assumed-shape array
   !       so a program unit that calls these subroutines must:
   !	   1. contain this subroutine as an internal subprogram,
   !	   2. import this subroutine from a module, or
   !	   3. contain an interface block for this subroutine.
   !
   !----------------------------------------------------------------------------

   !-Quicksort------------------------------------------------------------------
   !
   ! Subroutine to sort a list using the quicksort method. Call it with 
   ! First = the lower bound on the subscripts of the array and 
   ! Last  = the upper bound. 
   !
   ! Accepts : Array "Item", array "Indices"
   ! Returns : Array "Item"    (modified) with elements in ascending order
   !           array "Indices" (modified) with elements 
   !----------------------------------------------------------------------------
   
   RECURSIVE SUBROUTINE Quicksort(Item, First, Last, Indices,kk)
   !----------------------------------------------------------------------------
   ! This routine is based on a similar routine in "Fortran 90 for Engineers & 
   ! Scientists" by Nyhoff and Leestma.  I modified it to return an integer 
   ! array sorted based on the relationship of the DOUBLE PRECISION data in "Item".
   !
   ! Example:
   ! DOUBLE PRECISION,    dimension(100) :: randvals,randcopy
   ! integer, dimension(100) :: indarr = (/ (i, i=1,100) /)
   !  ...
   ! call random_number(randvals)		! F90 intrinsic subroutine
   ! randcopy = randvals			! save for comparison
   ! call Quicksort(randvals,1,size(randvals),indarr)
   ! print *,'sorted - indexed original is ',SUM(randvals - randcopy(indarr))
   !
   ! TJH 21 Oct 1998
   !----------------------------------------------------------------------------

      DOUBLE PRECISION,    DIMENSION(:), INTENT(INOUT) :: Item	! array of values
      INTEGER,               INTENT(IN)    :: First,Last,kk
      INTEGER, DIMENSION(:), INTENT(INOUT) :: Indices
   
      !--------------------------------------------------------------------
      ! Interface block(s) & Local Variables
      !--------------------------------------------------------------------

      INTERFACE
       SUBROUTINE Split(Item, Low, High, Mid, Indices)
          DOUBLE PRECISION,    DIMENSION(:), INTENT(INOUT) :: Item
          INTEGER,               INTENT(IN)    :: Low, High
          INTEGER,               INTENT(OUT)   :: Mid
          INTEGER, DIMENSION(:), INTENT(INOUT) :: Indices
       END SUBROUTINE Split
      END INTERFACE

      INTEGER	:: Mid,kk2
      
      !--------------------------------------------------------------------
      kk2=kk+1;
      print *,kk2
      IF (First < Last) THEN				! IF list size >= 2
         CALL Split(Item, First, Last, Mid, Indices)	! Split it
         CALL Quicksort(Item, First, Mid-1, Indices,kk2)	! Sort left  half
         CALL Quicksort(Item, Mid+1, Last,  Indices,kk2)	! Sort right half
      END IF

   END SUBROUTINE Quicksort

   !-Split----------------------------------------------------------------------
   !
   ! Subroutine to split a list into two sublists, using the first element 
   ! as a pivot, and return the position of the element about which the 
   ! list was divided. Local variables used are:
   ! Left	: position of the first element
   ! Right	: position of the last element
   ! Pivot	: pivot element
   ! Swap	: used to swap elements
   !
   ! Accepts:	Array Item and positions Low and High of the first and 
   !            last elements
   ! Returns:	Array Item (modified) with elements in ascending order
   !
   ! Note:	Item is an assumed-shape array so a program unit that calls
   !		this subroutine must:
   !		1. contain this subroutine as an internal subprogram,
   !		2. import this subroutine from a module
   !		3. contain an interface block for this subroutine.
   !----------------------------------------------------------------------------

   SUBROUTINE Split(Item, Low, High, Mid, Indices)

      DOUBLE PRECISION,    DIMENSION(:), INTENT(INOUT) :: Item
      INTEGER,               INTENT(IN)    :: Low, High
      INTEGER,               INTENT(OUT)   :: Mid
      INTEGER, DIMENSION(:), INTENT(INOUT) :: Indices


      INTEGER ::   Left, Right
      DOUBLE PRECISION    ::  Pivot,  Swap
      INTEGER :: iPivot, iSwap

      Left   = Low
      Right  = High
      Pivot  = Item(Low)
      iPivot = Indices(Low)
   
      ! Repeat the following while Left and Right haven't met
   
      DO
         IF ( Left >= Right ) Exit
   
         ! Scan right to left to find element < Pivot
   
         DO
   	    IF ( Left >= Right .OR. Item(Right) < Pivot ) EXIT
   	    Right = Right - 1
         END DO

         ! Scan left to right to find element > Pivot

         DO
   	    IF ( Left>= High .or. Item(Left) > Pivot) EXIT
   	    Left = Left + 1
         END DO
   
         ! If Left and Right haven't met, exchange the items
   
         IF (Left < Right) THEN
   	    Swap        = Item(Left)		! EXCHANGE THE ARRAY ITEMS
   	    Item(Left)  = Item(Right)
   	    Item(Right) = Swap

   	    iSwap          = Indices(Left)	! EXCHANGE THE INDICES ITEMS
   	    Indices(Left)  = Indices(Right)
   	    Indices(Right) = iSwap
         END IF
   
      END DO
   
      ! Switch element in split position with pivot
   
      Item(Low)   = Item(Right)			! SWITCH ARRAY ELEMS
      Item(Right) = Pivot
      Mid         = Right

      Indices(Low)   = Indices(Right)		! SWITCH ARRAY ELEMS
      Indices(Right) = iPivot

   END SUBROUTINE Split

   SUBROUTINE KB07AI(COUNT,N,INDEX)
!
!           KB07AI      HANDLES INTEGER VARIABLES
!THE WORK-SPACE 'MARK' OF LENGTH 100 PERMITS UP TO 2**50 NUMBERS
!TO BE SORTED.

!   .. Scalar Arguments ..
      INTEGER N
!   ..
!   .. Array Arguments ..
      DOUBLE PRECISION COUNT(*)
      INTEGER INDEX(*)
!   ..
!   .. Local Scalars ..
      DOUBLE PRECISION AV,X
      INTEGER I,IF,IFK,IFKA,INT,INTEST,IP,IS,IS1,IY,J,K,K1,LA,LNGTH,M, &
     &        MLOOP
!   ..
!   .. Local Arrays ..
      INTEGER MARK(100)
!   ..
!   .. Executable Statements ..
!SET INDEX ARRAY TO ORIGINAL ORDER .
      DO 10 I = 1,N
        INDEX(I) = I
   10 CONTINUE
!CHECK THAT A TRIVIAL CASE HAS NOT BEEN ENTERED .
      IF (N.EQ.1) GO TO 200
      IF (N.GE.1) GO TO 30
      WRITE (6,FMT=20)

   20 FORMAT (/,/,/,20X,' ***KB07AI** NO NUMBERS TO BE SORTED ** ', &
     & 'RETURN TO CALLING PROGRAM')

      GO TO 200
!'M' IS THE LENGTH OF SEGMENT WHICH IS SHORT ENOUGH TO ENTER
!THE FINAL SORTING ROUTINE. IT MAY BE EASILY CHANGED.
   30 M = 12
!SET UP INITIAL VALUES.
      LA = 2
      IS = 1
      IF = N
      DO 190 MLOOP = 1,N
!IF SEGMENT IS SHORT ENOUGH SORT WITH FINAL SORTING ROUTINE .
        IFKA = IF - IS
        IF ((IFKA+1).GT.M) GO TO 70
!********* FINAL SORTING ***
!( A SIMPLE BUBBLE SORT )
        IS1 = IS + 1
        DO 60 J = IS1,IF
          I = J
   40     IF (COUNT(I-1).LT.COUNT(I)) GO TO 60
          IF (COUNT(I-1).GT.COUNT(I)) GO TO 50
          IF (INDEX(I-1).LT.INDEX(I)) GO TO 60
   50     AV = COUNT(I-1)
          COUNT(I-1) = COUNT(I)
          COUNT(I) = AV
          INT = INDEX(I-1)
          INDEX(I-1) = INDEX(I)
          INDEX(I) = INT
          I = I - 1
          IF (I.GT.IS) GO TO 40
   60   CONTINUE
        LA = LA - 2
        GO TO 170
!           *******  QUICKSORT  ********
!SELECT THE NUMBER IN THE CENTRAL POSITION IN THE SEGMENT AS
!THE TEST NUMBER.REPLACE IT WITH THE NUMBER FROM THE SEGMENT'S
!HIGHEST ADDRESS.
   70   IY = (IS+IF)/2
        X = COUNT(IY)
        INTEST = INDEX(IY)
        COUNT(IY) = COUNT(IF)
        INDEX(IY) = INDEX(IF)
!THE MARKERS 'I' AND 'IFK' ARE USED FOR THE BEGINNING AND END
!OF THE SECTION NOT SO FAR TESTED AGAINST THE PRESENT VALUE
!OF X .
        K = 1
        IFK = IF
!WE ALTERNATE BETWEEN THE OUTER LOOP THAT INCREASES I AND THE
!INNER LOOP THAT REDUCES IFK, MOVING NUMBERS AND INDICES AS
!NECESSARY, UNTIL THEY MEET .
        DO 110 I = IS,IF
          IF (X.GT.COUNT(I)) GO TO 110
          IF (X.LT.COUNT(I)) GO TO 80
          IF (INTEST.GT.INDEX(I)) GO TO 110
   80     IF (I.GE.IFK) GO TO 120
          COUNT(IFK) = COUNT(I)
          INDEX(IFK) = INDEX(I)
          K1 = K
          DO 100 K = K1,IFKA
            IFK = IF - K
            IF (COUNT(IFK).GT.X) GO TO 100
            IF (COUNT(IFK).LT.X) GO TO 90
            IF (INTEST.LE.INDEX(IFK)) GO TO 100
   90       IF (I.GE.IFK) GO TO 130
            COUNT(I) = COUNT(IFK)
            INDEX(I) = INDEX(IFK)
            GO TO 110

  100     CONTINUE
          GO TO 120

  110   CONTINUE
!RETURN THE TEST NUMBER TO THE POSITION MARKED BY THE MARKER
!WHICH DID NOT MOVE LAST. IT DIVIDES THE INITIAL SEGMENT INTO
!2 PARTS. ANY ELEMENT IN THE FIRST PART IS LESS THAN OR EQUAL
!TO ANY ELEMENT IN THE SECOND PART, AND THEY MAY NOW BE SORTED
!INDEPENDENTLY .
  120   COUNT(IFK) = X
        INDEX(IFK) = INTEST
        IP = IFK
        GO TO 140

  130   COUNT(I) = X
        INDEX(I) = INTEST
        IP = I
!STORE THE LONGER SUBDIVISION IN WORKSPACE.
  140   IF ((IP-IS).GT. (IF-IP)) GO TO 150
        MARK(LA) = IF
        MARK(LA-1) = IP + 1
        IF = IP - 1
        GO TO 160

  150   MARK(LA) = IP - 1
        MARK(LA-1) = IS
        IS = IP + 1
!FIND THE LENGTH OF THE SHORTER SUBDIVISION.
  160   LNGTH = IF - IS
        IF (LNGTH.LE.0) GO TO 180
!IF IT CONTAINS MORE THAN ONE ELEMENT SUPPLY IT WITH WORKSPACE .
        LA = LA + 2
        GO TO 190

  170   IF (LA.LE.0) GO TO 200
!OBTAIN THE ADDRESS OF THE SHORTEST SEGMENT AWAITING QUICKSORT
  180   IF = MARK(LA)
        IS = MARK(LA-1)
  190 CONTINUE
  200 RETURN

      END
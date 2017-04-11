! F2003
! --------- UTILS Module ---------
!
!   General utilities and helper routines
!
! Contains:
!  - argReadKV         : Read key-value parameter from command line args
!  - confFileOpen      : Find and open configuration file, returns handle fh
!  - confScan          : Config parser core routine (see below for info)
!  - confScanString    : Config parser wrapper for string parameter
!  - confScanReal      :   "      "    wrapper for float parameter
!  - confScanInt       :   "      "    wrapper for integer parameter
!  - upperstr          : Convert string to upper-case
!
!  L.KEDWARD 2016

! #TODO: Theres a BUG in the config parser!
!        You can't currently have a parameter name which starts with another parameter name
!        E.G. mesh=foo   &    meshtype=su2
!

module UTILS
  use MATRIX, only: dp
  use FILEIO, only: newunit
  implicit none
  public
  
  contains
  
  
  !
  ! SUBROUTINE: argReadKV
  !  Read Key-Value parameter from command-line arguments
  !
  !  Example: command: 'program.exe --NCores 8'
  !           code usage: call argReadKV('--NCores',N)
  !
  !  IN:  [string] argKey     : key string used to identify parameter
  !  OUT: [string] argValue   : parameter value as string
  !       [bool]   found      : (optional) indicates whether argument
  !                              was found
  !
  subroutine argReadKV(argKey, argValue, found)
    character(*), intent(in) :: argKey
    character(*), intent(out), optional :: argValue
    logical, intent(out), optional :: found

    integer :: count, i
    character(100) :: arg

    if (present(found)) found = .False.

    count = command_argument_count()
    if (present(argValue)) count = count - 1

    do i=1,count
      call get_command_argument(i,arg)
      if (upperstr(trim(arg))==upperstr(argKey)) then
        if (present(found)) found = .True.
        if (present(argValue)) call get_command_argument(i+1,argValue)
        exit
      end if
    end do

    return

  end subroutine argReadKV
  ! --------------------------------------------------------------------
  
  
  
  !
  ! ----------------------- Config Parser ------------------------------
  !
  !  Subroutines for getting configuration parameter inputs.
  !  Parameters can be defined in a conf file and/or overridden using
  !   command line parameters.
  !  Default values can also be specified in code. In the absence of a
  !   parameter defined by config file, by command line or by default,
  !   then prompts user. (Therefore if you need to avoid prompts, always 
  !   provide default values)
  !
  !  Example usage: 
  !   - configFile.txt:
  !      NCORES = 8
  !      FILE=mesh.xyz
  !      reynolds =1E6
  !   - Code:
  !      open(unit=newunit(fh), 'configFile.txt')
  !      call confScanString(fh,'file',meshfile)
  !      call confScanInt(fh,'ncores',N,1)         ! default to 1 here
  !      call confScanReal(fh,'reynolds',RE)
  !   - Notes:
  !      - parameter names are case insensitive
  !      - spaces before/after names, values & equals-sign are ignored
  !      - parameters do not need to read in the same order as file
  !      - can override config file values from command line by:
  !          >> program.exe --NCORES 4 --FILE mesh2.xyz
  !      - to only use command line args without config file, set fh = -1
  !      - can reference other parameter values using $(varName):
  !          GEOM = naca0012
  !          FILE = mesh_$(geom).xyz
  !
  !  Subroutines:
  !    confFileOpen     - Find and open configuration file, returns handle fh
  !    confScan         - Core code for extracting string data
  !    confScanString | - Wrappers for specific types also providing
  !    confScanReal   |    default value and prompting functionality
  !    confScanInt    |   
  !  
  

  !
  ! SUBROUTINE: confFileOpen
  !  Code to check for config file: First checks if default filename
  !  exists, then for any filename specified in first command-line argument.
  !  If no config file is found, then assume all conf parameters are
  !  passed in command-line args and return fh=-1
  !
  ! IN:  [string]  defaultName : Name of default config file
  ! OUT: [integer] fh          : If config file is found, then valid file handle
  !                               else fh = -1
  !
  subroutine confFileOpen(defaultName, fh)
    character(*), intent(in) :: defaultName
    integer, intent(out) :: fh

    character(100) :: confFile
    logical :: fExist

    confFile = defaultName
    inquire(file=confFile, exist=fExist)
    if (fExist) then
      ! --- Found default conf file
      open(unit=newunit(fh), file=confFile, status="old")
    else
      ! --- Default conf file not found: check first command argument
      if (COMMAND_ARGUMENT_COUNT() == 0) then
        write(*,*) '(!) Configuration filename is missing'
        stop
      else
        call GET_COMMAND_ARGUMENT(1, confFile)
        inquire(file=confFile, exist=fExist)
        if (fExist) then
          ! --- Configuration filename found in first argument
          open(unit=newunit(fh), file=confFile, status="old")
        else
          ! --- No configuration file found
          !  Proceed assuming parameters are passed as command line args
          fh = -1
        end if
      end if
    end if


  end subroutine confFileOpen
  ! --------------------------------------------------------------------


  !
  ! SUBROUTINE: confScan
  !  Core code for reading and parsing parameter inputs
  !
  ! IN:  [integer] fh          : File unit handle for config file
  !                               (Set to -1 to not use config file)
  !      [string]  keystring   : Parameter name, e.g. 'NCORES'
  ! OUT: [string]  datastring  : Parameter value as string
  !      [bool]    foundStat   : (optional) indicates whether parameter
  !                               was found in config file or cmd args
  !
  recursive subroutine confScan(fh, keystring, datastring, foundStat)
    integer, intent(in) :: fh
    character(*), intent(in) :: keystring
    character(*), intent(out) :: datastring
    logical, intent(out), optional :: foundStat
    
    character(100) :: linestring, varName, varValue
    integer :: iloc, jloc, kloc, iostatus
    logical :: found, hasConfFile, varFound
    
    integer :: a1loc, a2loc, i1
    
    found = .False.

    ! --- FIRST: Look in cmd args ---
    ! (Command line argument will over-rule config file)
    call argReadKV('--'//keystring,datastring,found)
    
    
    ! --- SECOND: Look in conf file if provided ---
    inquire(fh, opened=hasConfFile)
    if (.Not.found .and. hasConfFile) then
      read(fh, '(A)', iostat=IOstatus) linestring
      do while (IOstatus == 0)
      
        iloc = index(adjustl(upperstr(linestring)), & 
                     trim(adjustl(upperstr(keystring))) )
        jloc = index(linestring, '=')
        kloc = index(adjustl(linestring), '#')
        
        if ( (iloc == 1).AND.(jloc > 0).AND.(kloc /= 1) ) then
          found = .True.
          datastring = trim(adjustl(linestring(jloc+1:)))
          exit
        end if
          
        read(fh, '(A)', iostat=IOstatus) linestring
      end do
      rewind(fh)   ! Set back to beginning for subsequent calls
    end if
    
    
    ! ------ Check for and parse variables -----
    i1 = 1
    do while ( i1 < len(datastring))
      a1loc = index(datastring(i1:), '$(')
      a2loc = index(datastring(i1:), ')')
      if ((a1loc > 0).AND.(a2loc > 0).AND.(a2loc>a1loc)) then
        varName = datastring(i1+a1loc+1:i1+a2loc-2)
        call confScan(fh, varName, varValue, varFound)
        if (.Not.varFound) varValue = varName
        datastring = datastring(:i1+a1loc-2) &
                        //trim(adjustl(varValue))//datastring(i1+a2loc:)
        i1 = i1+a1loc+1!-2 + len(varValue)
      else
        i1 = len(datastring)
      end if 
    end do
    
    if (present(foundStat)) foundStat = found
    
  end subroutine confScan
  ! --------------------------------------------------------------------
  
  
  !
  ! SUBROUTINE: confScanString
  !  Get string parameter value
  !
  ! IN:  [integer] fh            : File unit handle for config file
  !                                 (Set to -1 to not use config file)
  !      [string]  key           : String parameter name, e.g. 'label'
  !      [string]  dataValue     : String parameter value
  !      [string]  defaultValue  : (Optional) Default string value
  !      [bool]    interact      : Behaviour when input is missing and
  !                                 no default value is specified:
  !                                 TRUE  - Prompt user for value
  !                                 FALSE - Do not prompt user and halt
  !      [bool]    foundStat     : (optional) TRUE if value found in 
  !                                 conf file or cmd args. FALSE if
  !                                 value from default or user prompt
  !
  subroutine confScanString(fh, key, dataValue, defaultValue, interact, foundStat)
    integer, intent(in) :: fh
    character(*), intent(in) :: key
    character(*), intent(out) :: dataValue
    character(*), intent(in), optional :: defaultValue
    logical, intent(in), optional :: interact
    logical, intent(out), optional :: foundStat
    
    logical :: found
    
    call confScan(fh, key, dataValue, found)
    
    ! --- Check if found in args/conf
    if (.Not.found) then
      if (present(defaultValue)) then
        ! Use provided defaul value
        dataValue = defaultValue
      else
        ! All else has failed :-(
        if (present(interact).AND.interact) then
          ! Prompt user if interactive = True
          write(*,*) "Enter a value for parameter, '", key, "' :"
          read(*,*) dataValue
        else
          ! Otherwise fail
          write(*,*) "(!) Required parameter, '", key, "' is missing"
          write(*,*) " Stopping."
          stop
        end if
      end if
    end if
    
    if (present(foundStat)) foundStat = found
    
  end subroutine confScanString
  ! --------------------------------------------------------------------
  
  
  !
  ! SUBROUTINE: confScanReal
  !  Get float parameter value
  !
  ! IN:  [integer] fh            : File unit handle for config file
  !                                 (Set to -1 to not use config file)
  !      [string]  key           : Float parameter name, e.g. 'label'
  !      [double]  dataValue     : Float parameter value
  !      [double]  defaultValue  : Default float parameter value
  !      [bool]    foundStat     : (optional) TRUE if value found in 
  !                                 conf file or cmd args. FALSE if
  !                                 value from default or user prompt
  !
  subroutine confScanReal(fh, key, dataValue, defaultValue, interact, foundStat)
    integer, intent(in) :: fh
    character(*), intent(in) :: key
    real(dp), intent(out) :: dataValue
    real(dp), intent(in), optional :: defaultValue
    logical, intent(in), optional :: interact
    logical, intent(out), optional :: foundStat
    
    character(100) :: datastring
    logical :: found
    
    call confScan(fh, key, datastring, found)
    
    if (found) then
      ! Parse datastring for real data
      read(datastring,*) dataValue
    else
      if (present(defaultValue)) then
        ! Use provided defaul value
        dataValue = defaultValue
      else
        ! All else has failed :-(
        if (present(interact).AND.interact) then
          ! Prompt user if interactive = True
          write(*,*) "Enter a value for parameter, '", key, "' :"
          read(*,*) dataValue
        else
          ! Otherwise fail
          write(*,*) "(!) Required parameter, '", key, "' is missing"
          write(*,*) " Stopping."
          stop
      end if
      end if
    end if
    
    if (present(foundStat)) foundStat = found
    
  end subroutine confScanReal
  ! --------------------------------------------------------------------
  
  
  !
  ! SUBROUTINE: confScanInt
  !  Get integer parameter value
  !
  ! IN:  [integer] fh            : File unit handle for config file
  !                                 (Set to -1 to not use config file)
  !      [string]  key           : Integer parameter name, e.g. 'label'
  !      [integer]  dataValue    : Integer parameter value
  !      [integer]  defaultValue : Default integer parameter value
  !      [bool]    foundStat     : (optional) TRUE if value found in 
  !                                 conf file or cmd args. FALSE if
  !                                 value from default or user prompt
  !
  subroutine confScanInt(fh, key, dataValue, defaultValue, interact, foundStat)
    integer, intent(in) :: fh
    character(*), intent(in) :: key
    integer, intent(out) :: dataValue
    integer, intent(in), optional :: defaultValue
    logical, intent(in), optional :: interact
    logical, intent(out), optional :: foundStat
    
    character(100) :: datastring
    logical :: found
    
    call confScan(fh, key, datastring, found)
    
    if (found) then
      ! Parse datastring for real data
      read(datastring,*) dataValue
    else
      if (present(defaultValue)) then
        ! Use provided defaul value
        dataValue = defaultValue
      else
        ! All else has failed :-(
        if (present(interact).AND.interact) then
          ! Prompt user if interactive = True
          write(*,*) "Enter a value for parameter, '", key, "' :"
          read(*,*) dataValue
        else
          ! Otherwise fail
          write(*,*) "(!) Required parameter, '", key, "' is missing"
          write(*,*) " Stopping."
          stop
      end if
      end if
    end if
    
    if (present(foundStat)) foundStat = found
    
  end subroutine confScanInt
  ! --------------------------------------------------------------------
  
  
  
  !
  ! FUNCTION: upperstr
  !  Return copy of string converted to uppercase
  !  Used for case-insensitive string comparison
  !
  ! 1996, John S. Urban http://fortranwiki.org/fortran/show/ufpp
  character(len=len(linei)) function upperstr(linei)
    character(len=*),intent(in) :: linei               ! input string to convert to uppercase

    intrinsic ichar, char, len
    integer :: inlen ! number of characters in trimmed input string
    integer :: i10 ! counter to increment through input and output string
    integer :: ilet ! current character being converted represented using ASCII Decimal Equivalent
    
    inlen=len_trim(linei) ! number of characters to convert to uppercase
    upperstr=' '  ! initialize output string to all blanks

    if(inlen.gt.len(upperstr))then ! make sure there is room to store the output characters
      write(*,'(a)')'*ufpp* FATAL - OUTPUT TOO LONG TO CONVERT TO UPPERCASE:'
    endif
    
    ! loop through each character in input string
    do i10=1,inlen,1                                   
      ilet=ichar(linei(i10:i10))                ! current character in input to convert to output converted to ADE
      if( (ilet.ge.97) .and. (ilet.le.122))then ! lowercase a-z in ASCII is 97 to 122; uppercase A-Z in ASCII is 65 to 90
         upperstr(i10:i10)=char(ilet-32)        ! convert lowercase a-z to uppercase A-Z
      else
         upperstr(i10:i10)=linei(i10:i10)       ! character is not a lowercase a-z, just put it in output
      endif
    enddo
  end function upperstr
  ! --------------------------------------------------------------------
  
  
  
end module UTILS

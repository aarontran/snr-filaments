c Short program to start learning FORTRAN77 syntax
c Following http://web.stanford.edu/class/me200c/tutorial_77/

c Aaron Tran
c Summer 2014 (57 years later...)

c ===================
c Some starting notes
c ===================

c Another good reference: http://www.star.le.ac.uk/~cgp/prof77.html

c Column layout in FORTRAN (fixed format for FORTRAN 77)
c Col. 1    : blank
c Col. 1-5  : statement label (optional)
c Col. 6    : continuation character (anything, prefer [+&0-9])
c Col. 7-72 : statements
c Col. 73-80: sequence number (optional, rarely used)

c Conveniently, vim auto-indents statements 6 spaces,
c and cuts statements off at 72 columns

c FORTRAN character set: 26 letters, 10 digits, 13 special characters
c Special characters are: + (plus), - (minus), * (asterisk), / (slash)
c   (blank), = (equals), ( (left paren), ) (right paren), . (decimal pt)
c , (comma), ' (apostrophe), : (colon), $ (currency symbol)

c Inline comments with exclamation marks (!) are not FORTRAN 77 standard

c FORTRAN 90 standard
c 1. Free form layout (up to 132 char)
c 2. Variable names up to 31 char (with underscores)
c 3. implicit none, ! comments, more...

c ============
c Main program
c ============

c Only one program <name> in every fortran file

      program main
          ! Declarations must precede executable statements
          implicit none  ! This is not ANSI FORTRAN?
          real r, area, pi, avgdro
          logical truth
          ! Parameter declarations can't be changed later
          parameter (pi = 3.1415926535897932384626)
          parameter (avgdro = 6.022E22)
          parameter (truth = .TRUE.)

          write(unit=*, fmt=*) 'Hello world!'

          write(*,*) 'Radius:'
          read(*,*) r
          area = pi*r**2
          write(*,*) area

          write(*,*) 'Avogadro''s number is:', avgdro

          if (5 .LT. 10) write(*,*) 'Five is less than ten'

          if (.NOT. 5 .GE. 10) then
              write(*,*) 'Five is not greater than ten'
              write(*,*) 'If statements are great, if they are great'
          else
              write(*,*) 'I''ll never be printed!'
          endif

          call submain
          call pssref

          stop
      end

c Possible type declarations are: integer, real, double precision,
c complex, logical, character.  integer, real are 4 bytes (32 bit)

c =================================
c Subroutines, functions, and loops
c =================================

      subroutine submain
          implicit none
          integer i, n, fact, exp2
          real circum

          write(*,*) 'Compute factorial for integer:'
          read(*,*) n

          fact = 1
          ! 10 is a statement label, typically consec. multiples of 10
          do 10 i = 1, n
              fact = fact * i
              write (*,*) 'i =', i
              write (*,*) 'fact =', fact
   10     continue  ! Many F77 compilers accept enddo, but it's not
          ! part of ANSI FORTRAN

          ! ANSI FORTRAN doesn't have WHILE, FOR loops
          ! must implement using DO, IF, GOTO

          exp2 = 2
   20     if (exp2 .LE. 100) then
              exp2 = exp2 * 2
              goto 20
          endif

          write(*,*) 'The smallest power of 2 greater than 100 is', exp2

          write(*,*) 'The circumference of a circle with radius 1 is'
          write(*,*) circum(1.0)
          
          return
      end

      real function circum(r)
          implicit none
          real pi, r  ! Don't need to declare circum
          parameter (pi = 3.1415926535897932384626)

          circum = 2*pi*r

          return
      end

c FORTRAN 77 doesn't explicitly allow recursion, but compiler might
c FORTRAN is said to be pass-by-reference, here's an example

      subroutine iswap(a,b)
          implicit none
          integer a, b, tmp

          tmp = a
          a = b
          b = tmp

          return
      end


      subroutine pssref
          implicit none
          integer m, n

          m = 1
          n = 9001
          
          write(*,*) 'Demonstration of pass-by-reference'
          write(*,*) 'Before:', m, n
          call iswap(m,n)
          write(*,*) 'After: ', m, n
          ! The variables were passed in by reference
          ! reference (memory address) value is updated by iswap
          ! so variables in caller's scope are modified

          return
      end

c ==============
c FORTRAN arrays
c ==============

      subroutine arrfnc(m,n)
          implicit none
          integer i, j, m, n
          real mat(m, n)

          do 30 j = 1, n      ! Inner-most loop should iterate over row
              do 40 i = 1, m  ! index, more efficient memory access
                  mat(i,j) = i*j
   40         continue
   30    continue
         write(*,*) 'Array element at 2,2 is:', mat(2,2)

         return
      end

c FORTRAN 77 does not allocate memory dynamically, so preallocate arrays
c larger than needed.  Arrays are stored in column-major order (i.e.,
c elements ordered (1,1), (2,1), (3,1), ... (m,1), (1,2), ...
c So allocate more columns than needed, but not more rows

c In subroutines, all but the last matrix dimension MUST be specified
c e.g., declare `real mat(m, n, *)` and not `mat(*,*,*)`.
c However, typically arrays should be declared in the main program and
c passed to subroutines to modify.


c Common blocks and more

c Avoid common blocks (basically, global), but the syntax is pretty
c simple.  data statement for initializing variables/arrays isn't that
c useful?

c ========
c File I/O
c ========

      subroutine fileio(fname)
          implicit none
          character(len=*) fname

          open(unit=1, file=fname)  ! unit specifies the file...
          ! open(1, file=fname)     ! works as well

          read(1,*) n
          write(1,*) 'lalalala'

          ! Formatting statements

          write(1, 100) 100, 2.71828
          write(*, 200) 'Wrote', 100, 2.71828
  100     format(I4, F8.3)
  200     format(A, I4, F8.3)

          close(1)

          return
      end


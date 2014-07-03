c Short program to start learning FORTRAN77 syntax
c Following http://web.stanford.edu/class/me200c/tutorial_77/

c Aaron Tran
c Summer 2014 (57 years later...)

c ===================
c Some starting notes
c ===================

c Another good reference: http://www.star.le.ac.uk/~cgp/prof77.html

c Column layout in FORTRAN
c Col. 1    : blank
c Col. 1-5  : statement label (optional) (what is this?)
c Col. 6    : continuation character (anything, prefer [+&0-9])
c Col. 7-72 : statements
c Col. 73-80: sequence number (optional, rarely used)

c Conveniently, vim auto-indents statements 6 spaces,
c and cuts statements off at 72 columns

c FORTRAN character set: 26 letters, 10 digits, 13 special characters
c Special characters are: + (plus), - (minus), * (asterisk), / (slash)
c   (blank), = (equals), ( (left paren), ) (right paren), . (decimal pt)
c , (comma), ' (apostrophe), : (colon), $ (currency symbol)

c Inline comments with exclamation marks (!) are not FORTRAN standard

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

c Coming soon!

DEFINT H-N
' Models for the simulations of the color pattern on the shells of mollusks
' see also: Meinhardt,H. and Klingler,M. (1987) J. theor. Biol 126, 63-69
' see also: H.Meinhardt: "Algorithmic beauty of sea shells"
' (Springer Verlag) (c) H.Meinhardt, T¸bingen

'This is a short version of a program for the simulations of the color
'patterns on tropical sea shells, here 'Oliva porphyria'.
'An autocatalytic activator a(i) leads to a burst-like activation
'that is regulated back by the action of an inhibitor b(i). The life
'time of the inhibitor is regulated via a hormone c, that is
'homogeneously distributed along the growing edge. Whenever the number
'of activated cells cells become too small, active cells remain activated
'until backwards waves are triggered

'The program runs with the interpreter QBASIC, but this is very slow.
'Better are the compiler Power Basic, Microsoft QB 4.5, Professional
'Basic 7.1, and Visual Basic fÅr DOS. A compiled version is included

' i = 1...kx < imax = cells at the growing edge

imax = 640: DIM a(imax), b(imax), z(imax)
RANDOMIZE TIMER ' By different fluctuations,
                 'simulation will be slightly different
KT = 455  'Number of displays
'          it*ip = number of iterations in total
KP = 12   'number of iterations between the displays ( = lines on the screen)
dx = 1    'width of a cell in pixel;   with kp=6 ; kx=315 and dx=2 =>
kx = 630  'Number of cells               simulation in a smaller field
kx = 630  'number of cells
DA = .015 'Diffusion of the activator
RA = .1   'Decay rate of the inhibitor
BA = .1   'Basic production of the activator
SA = .25  'Saturation of the autocatalysis
DB = 0    'Diffusion of the inhibitor
RB = .014 'Decay rate of the inhibitor
SB = .1   'Michaelis-Menten constant of the inhibition
RC = .1   'Decay rate of the hormone

start:
REM ----------- Initial condition  --------------------------
FOR i = 1 TO kx
   a(i) = 0    'Activator, general initial concentration
   b(i) = .1   'Inhibitor, general initial concentration
   z(i) = RA * (.96 + .08 * RND)'Fluctuation of the autocatalysis
NEXT i
C = .5 'Hormone-concentration, homogeneous in all cells
i = 10: FOR j = 1 TO 20  'initially active cells
a(i) = 1: i = i + 50 * RND + 10: IF i > imax THEN EXIT FOR
NEXT
DAC = 1! - RA - 2! * DA ' These constant factors are used again and again
DBC = 1! - RB - 2! * DB ' therefore, they are calculated only once
DBCC = DBC'             ' at the begin of the calculation
SCREEN 12
WINDOW (1, 1)-(640, 480)
continuo:
v1 = TIMER
LINE (1, 1)-(640, 480), 1, BF 'background blue
x0 = 5: y1 = 477       'Initial position of the first line
itot = 0               'Number of iteration calculated so far
FOR itot = 0 TO KT
FOR iprint% = 1 TO KP  ' Begin of the iteration
REM -----  --- Boundary impermeable
A1 = a(1) '    a1 is the concentration  of the actual cell. since this
B1 = b(1) '    concentration is identical, no diffusion through the border.
a(kx + 1) = a(kx) '          Concentration in a virtual right cell
b(kx + 1) = b(kx)
BSA = 0!  '    This will carry the sum of all activations of all cells
REM ---------- Reactions  ------
FOR i = 1 TO kx' i = actual cell, kx = right cell
   AF = a(i) 'local activator concentration
   BF = b(i) 'local inhibitor concentration
   AQ = z(i) * AF * AF / (1! + SA * AF * AF)  'Saturating autocatalysis
' Calculation of the new activator and inhibitor concentration in cell i:
   a(i) = AF * DAC + DA * (A1 + a(i + 1)) + AQ / (SB + BF)
   ' 1/BF => Action of the inhibitor; SB = Michaelis-Menten constant
   b(i) = BF * DBCC + DB * (B1 + b(i + 1)) + AQ 'new inhibitor conc.
   BSA = BSA + RC * AF 'Hormone production -> Sum of activations
   A1 = AF '    actual concentration will be concentration in left cell
   B1 = BF '    in the concentration change of the next cell
   NEXT i
   C = C * (1! - RC) + BSA / kx ' New hormone concentration , 1/kx=normalization
   RBB = RB / C                 ' on total number of cells
   'RBB => Effective life time of the inhibitor
   DBCC = 1! - 2! * DB - RBB      ' Change in a cell by diffusion
   NEXT iprint%                   ' and decay. Must be recalculated since
                                  ' lifetime of the inhibitor changes
   REM ----------------Plot -------------
   y1 = y1 - 1 'Next plot, one line below
   LINE (x0, y1)-(x0 + dx * kx, y1), 15 'Background white
   FOR ix% = 1 TO kx   'Pigment is drawn when a is higher than a threshold
   IF a(ix%) > .5 THEN LINE (x0 + dx * (ix% - 1), y1)-(x0 + dx * ix%, y1), 6
   NEXT ix%
   IF INKEY$ > "" THEN EXIT FOR
 NEXT itot
     v2 = v1 - TIMER
LOCATE 30, 1:  PRINT "c = continue; s = a new start, all other keys = End (c) H.Meinhardt      ";
PRINT USING "##.##"; TIMER - v1;
DO: a$ = INKEY$: IF a$ > "" THEN EXIT DO
LOOP
IF a$ = "c" GOTO continuo
IF a$ = "s" GOTO start
END


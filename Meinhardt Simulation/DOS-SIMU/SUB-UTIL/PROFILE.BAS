DECLARE SUB init (inewscreen%, igraph%, ivorcol%, ibackcol%, ivgashift%)
DECLARE SUB cmessage (itype%, iline%, irow%, ilettercol%, iletterback%, ibackcol%, text$, resp$)
DECLARE SUB linep (ihp%, x1!, y1!, x2!, y2!, ipat%, iwidth%, ico%, ityp%)
DECLARE SUB lptopen (lpt$)
DECLARE SUB switchscreen (inz%, iactive%, ibackcol%)
DECLARE SUB zeingabe (igt%, iquest%, inz%, i%, f!, LABEL$, text$)
DECLARE SUB prof (a!(), il%, dw$, ipattern%, inz%, iendm%, DX!, ff!, ja%, js%, iy%, xo!, yo!)
DECLARE SUB showprof (inz%, iendm%, dw$, a!(), ja%, js%, icc%)
'XT-Subroutines for the simulation of pigmentation patterns on sea shells
' (c) by Hans Meinhardt, 1994

DEFINT H-N
' $INCLUDE: 'sub-util\scommon.bi' nur fr Microsoft-Basic gebraucht
' $INCLUDE: 'sub-util\errlist.bas' nur fr Microsoft-Basic gebraucht

SUB manipulationxt (pcontrol2$, ja, js, kinsert, imxl, imxl2, icc) STATIC
'for saving or manipulation of distributions, growth etc.
SHARED ivgashift, displaywas$, iparam, itot&, K1, igt, iactive
SHARED a(), kax, KX, kay, ky, KD, KN, imx
SHARED inv(), ivorcol, ibackcol, igraph, LABEL$()
'pcontrol2$: what to do; ja, js leftmost, rightmost cell
'kinsert: where to insert a new cell, imxl: max number of cells

IF ja = 0 THEN
   mess$ = "left cell is zero, run program first"
  CALL cmessage(6, -1, 1, 15, 4, ibackcol, mess$, "")
  itot& = 99999
  EXIT SUB
END IF
mamessage$ = ""
manipulationxtlabel:
SELECT CASE pcontrol2$
CASE "Q", CHR$(27): EXIT SUB
CASE "S"'to store present concentrations in memory arrays
  GOSUB storesubroutine
  mamessage$ = " Distributions are stored"
  GOTO mandisplay
CASE "R"'to restore former state of cells from memory
  mamessage$ = " Distributions are restored"
  GOSUB restoresubroutine:
  GOTO mandisplay

CASE "C"'to change margins of the field
mess$ = "Store present or Restore previous concentrations: S/R/RETURN=nothing"
CALL zeingabe(igt, 0, 3, kanew, dummy, resp$, mess$)
  SELECT CASE resp$
  CASE CHR$(27): EXIT SUB
  CASE "s", "S"
    GOSUB storesubroutine
CASE "r", "R"
    GOSUB restoresubroutine
  END SELECT
ksnewinput:
  kanew = ja: ksnew = js
  mess$ = "Input new (smaller) field size: Present field size is:" + STR$(ja) + " - " + STR$(js)
  CALL cmessage(5, -2, 1, 4, ibackcol, ibackcol, mess$, "")
CALL zeingabe(igt, 1, 1, kanew, dummy, dummy$, "New left-most cell")
CALL zeingabe(igt, 1, 1, ksnew, dummy, dummy$, "New right-most cell")
  IF ksnew = 0 THEN ksnew = js
  IF ksnew > js THEN
  ksnew = js
   mess$ = "New field size must be smaller than the old one, left unchanged"
  resp$ = "ok"
  CALL cmessage(4, -1, 1, 15, 4, ibackcol, mess$, resp$)
  END IF
  IF kanew < ja THEN kanew = ja
  js = ksnew
  ja = kanew
  icc = icc + 3
  CALL cmessage(3, -2, 1, 4, ibackcol, ibackcol, SPACE$(80), "")
  GOTO mandisplay
CASE "I"'insertion of new cells
CALL zeingabe(igt, 1, 1, kinsert, dummy, dummy$, " Position of the new cell, RETURN=right margin")
  IF kinsert = 0 THEN kinsert = js
  GOSUB growth

CASE "G"
  GOSUB growthpreparation
CASE "4", "5"
  inz = VAL(pcontrol2$)
  CALL showprof(inz, 1, LEFT$("sabcdefgh", KN + 1) + "h", a(), ja, js, icc)
  iparam = 1
  EXIT SUB
CASE "A"
  itot& = 0
CALL zeingabe(igt, 0, 1, il, dummy, dummy$, "Array to be changed, 0=source, 1=a,2=b....")
  IF il > KN THEN
  mess$ = "This substance does not exist, must be <= KN = " + STR$(KN)
  resp$ = "ok": CALL cmessage(5, -1, 1, 15, 4, ibackcol, mess$, resp$)
  EXIT SUB
  END IF
 mess$ = "Present field size is " + STR$(ja) + " - " + STR$(js)
    CALL cmessage(5, -2, 1, 4, ibackcol, ibackcol, mess$, "")
  i1 = ja: i2 = js
CALL zeingabe(igt, 1, 1, i1, dummy, dummy$, "First cell to be changed")
CALL zeingabe(igt, 1, 1, i2, dummy, dummy$, "Last cell to be changed")
CALL zeingabe(igt, 0, 2, i2, cfactor, dummy$, " factor [1 for no change] ")
    CALL cmessage(3, -2, 1, ibackcol, ibackcol, ibackcol, SPACE$(80), "")
  IF i1 < ja THEN i1 = ja
  IF i2 = 0 OR i2 > js THEN i2 = js
  FOR i = i1 TO i2: a(il, i) = a(il, i) * cfactor: NEXT
  GOTO mandisplay
CASE "T"
  itot& = 0
CALL zeingabe(igt, 0, 3, kanew, dummy, resp$, "Store present or Restore previous concentrations S/R/RETURN=nothing")
  SELECT CASE resp$
  CASE "s", "S":   GOSUB storesubroutine
  CASE "r", "R":   GOSUB restoresubroutine
  END SELECT
  ja = 1: i = 0
  mess$ = "Input new (smaller) field size: Present field size is:" + STR$(ja) + " - " + STR$(js)
  CALL cmessage(3, -2, 1, 4, ibackcol, ibackcol, mess$, "")
fromto:
ifrom = kamem: ito = ksmem
CALL zeingabe(igt, 1, 1, ifrom, dummy, dummy$, "Left-most cell of the piece to be added, 0=finish")
IF ifrom = 0 THEN
mamessage$ = "Transplantation complete"
icc = 1': LINE (0, 0)-(640, 40), 15, BF
GOTO mandisplay
END IF
CALL zeingabe(igt, 1, 1, ito, dummy, dummy$, "Right-most cell of the piece to be added")
IF ito = 0 THEN ito = ksmem
  IF ito > ksmem OR ifrom < kamem THEN
CALL zeingabe(igt, 0, 3, i, dummy, resp$, "Sorry, this piece is larger than the field, type RETURN")
    GOTO fromto
END IF
  idirection = 1
  IF ifrom > ito THEN idirection = -1
  FOR im = ifrom TO ito STEP idirection
    i = i + 1
IF i >= imxl THEN
CALL cmessage(6, -1, 1, 15, 4, ibackcol, " Sorry, XT-array too large!", "OK")
    END IF
    FOR il = 0 TO KN: a(il, i) = a(il, im + imxl2): NEXT
  IF genedw$ > "" THEN
IF i >= imxl THEN
CALL cmessage(6, -1, 1, 15, 4, ibackcol, " Sorry, XY-array too large!", "OK")
    END IF
    FOR il = 0 TO KN: FOR iy = 1 TO ky(il)
    ax(il, i, iy) = ax(il, im, iy + kyadd)
    NEXT iy:  NEXT il
  END IF
  NEXT im
  js = i: FOR il = 0 TO KN: KX(il) = i: NEXT
  GOTO fromto
CASE ELSE
  CLS
  LOCATE 10, 1
  PRINT "S = Save current results in a save buffer"
  PRINT "R = Restore previous results from the save buffer"
  PRINT "C = Change field size (for simulation of regeneration)"
  PRINT "T = Transplantation fo several pieces"
  PRINT "I = Insert a new cell, also growth at marinal positions"
  PRINT "A = Change activator etc.  [a()] in some cells"
  PRINT "Q = quit without action"
CALL zeingabe(igt, 0, 3, i, dummy, pcontrol2$, "Type of manipulation: ")
IF PCONTROL$ = "" THEN EXIT SUB
pcontrol2$ = UCASE$(pcontrol2$)
  GOTO manipulationxtlabel
END SELECT'  'select of manipulation
EXIT SUB

mandisplay:
IF igraph = 9 THEN CALL switchscreen(1, iactive, ibackcol)
CALL showprof(KD, 1, displaywas$, a(), ja, js, icc)
CALL cmessage(7, -1, 1, 1, 15, ibackcol, mamessage$, mamessage$)
iparam = 2
EXIT SUB

storesubroutine:
ksmem = js
kamem = ja
IF js + imxl2 > imxl THEN
  mess$ = "sorry, array to large to be stored, max = " + STR$(imxl2)
resp$ = "ok"
    CALL cmessage(4, -1, 1, 15, 4, ibackcol, mess$, resp$)
  EXIT SUB
END IF
FOR il = 0 TO KN
  FOR i = ja TO js
    a(il, i + imxl2) = a(il, i)
  NEXT i: NEXT il
RETURN

restoresubroutine:
IF ksmem = 0 THEN
  mess$ = "sorry, nothing was stored in the memory"
    CALL cmessage(4, -1, 1, 15, 4, ibackcol, mess$, "OK")
  EXIT SUB
END IF
itot& = 0
js = ksmem

ja = kamem
FOR il = 0 TO KN
  FOR i = ja TO js
    a(il, i) = a(il, i + imxl2)
  NEXT i: NEXT il
RETURN

growthpreparation:
SELECT CASE K1'Where to insert additional cell
      CASE 0'terminal growth
	IF ja > 1 THEN
	  ja = ja - 1
	  FOR il = 1 TO KN: a(il, ja) = a(il, ja + 1): NEXT il
	END IF
	kinsert = js
      CASE 1'two insertion, one at each half
	km = ja + (js - ja) / 2
	kinsert = ja + RND * (km - ja)'insertion in the first half
	IF ja > 1 THEN
	  ja = ja - 1
	  FOR il = 1 TO KN
	    FOR ix = ja TO kinsert: a(il, ix) = a(il, ix + 1): NEXT ix
	  NEXT il
	END IF
	kinsert = km + RND * (js - km)'insertion in the second half
      CASE 2: kinsert = ja + (js - ja) / 2'central
      CASE 3: kinsert = ja + (js - ja) / 3 'one-third
      CASE 4: kinsert = ja + (js - ja) * RND'random position
	IF kinsert > js THEN kinsert = js
      CASE ELSE
	kinsert = js
END SELECT

growth: 'subroutine
IF js = imyl THEN
    mess$ = " No further growth is possible, too many cells"
    CALL cmessage(1, -1, 1, 1, ibackcol, ibackcol, mess$, "")
ELSE
  FOR il = 1 TO KN
    FOR i = js TO kinsert STEP -1'still the old js is being used!
      a(il, i + 1) = a(il, i)
    NEXT i: NEXT il
  js = js + 1' js is adjusted to the new field size
END IF
RETURN
END SUB

SUB prof (a(), il, dw$, ipattern, inz, iendm, DX, ff, ja, js, iy, xo, yo) STATIC
'this routine makes the drawings----------------------------------------
'a(ix,il): array to be plotted    il: substance, 1=a, 2=b, ...
'dw$: display what (abc..)        ipattern: controlling some patterns
'inz: type of the plot            iendm: = 99 if plot is to be finished
'ff: normalization factor         ja: left cell js: right cell
'iy: counting index, =1 -> initialization, xo, yo: origin
'------------------------------------------------------------------------
SHARED igraph, x1, y1, dy, ivorcol, ibackcol, krot, per
SHARED icol(), iwi(), ityp()
SHARED ivgashift, msiluet(), msiluetmax, ihp
DIM sp$(11), ipat(1 TO 4, 1 TO 5), icolf(10), col(3, 10)
icolor = icol(il)
f10 = 2.5 * ff
f1 = 40 * ff
'colours for false colour plot:
icolf(0) = 15: icolf(1) = 14: icolf(2) = 11: icolf(3) = 10: icolf(4) = 7: icolf(5) = 13
icolf(6) = 3: icolf(7) = 2: icolf(8) = 6: icolf(9) = 5: icolf(10) = 0
twas$ = "sabcdefghij"
tv$(il) = MID$(twas$, il + 1, 1)
iwhatarray = 1: iwhatarray = INSTR(iwhatarray, dw$, tv$(il)) - 1
IF iwhatarray = -1 THEN EXIT SUB
iplottotal = LEN(dw$)

profilestart:
IF iendm = 99 GOTO profileframe
IF iy = 1 AND ahidden = 0 THEN '----------- Initialization ---------------------
PALETTE
  SELECT CASE inz
  CASE 1, 2  '3D-plots
    msiluetm5 = msiluetmax / iplottotal - 1
    ahidden = -1: IF inz = 2 THEN ahidden = 1000
    FOR i = 1 TO msiluetmax: msiluet(i) = ahidden: NEXT i
    IF krot = 0 THEN krot = 2
    angle = ATN(1! / (krot))
    pdx = DX * COS(angle)
    px = pdx / krot
    pyy = DX * per * COS(angle)
    pdyx = DX * SIN(angle) * per
    pyx = pdyx / krot
    f1 = 40 * ff: fmodus = f1 / krot
  CASE 3, 18, 19'------new curves are drawn
  CASE 4, 5'----- schematic plot with letters --------
    spp$ = "N .,:;°±²ÛÛ"
    FOR i% = 0 TO 9: sp$(i%) = MID$(spp$, i% + 1, 1): NEXT i%
  CASE 6 TO 14, 22'pixelpattern for pixel plots
    ipat(1, 1) = &H4444
    ipat(1, 2) = &H1111
    ipat(1, 3) = &H8888
    ipat(1, 4) = &H2222
    ipat(2, 1) = &HAAAA
    ipat(2, 2) = &H5555
    ipat(2, 3) = &H3333
    ipat(2, 4) = &HCCCC
    ipat(3, 1) = &HEEEE
    ipat(3, 2) = &HBBBB
    ipat(3, 3) = &HDDDD
    ipat(3, 4) = &H7777
    ipat(4, 1) = &HFFFF
    ipat(4, 2) = &HFFFF
    ipat(4, 3) = &HFFFF
    ipat(4, 4) = &HFFFF
    'for false colour,inz=6-8
    icolf(0) = 15: icolf(1) = 14: icolf(2) = 11: icolf(3) = 10: icolf(4) = 7: icolf(5) = 13
    icolf(6) = 3: icolf(7) = 2: icolf(8) = 6: icolf(9) = 5: icolf(10) = 0
CASE 15: alpha = 0: iradius = 50: wi = -.01 * DX'for circular plots
CASE 16'(used for plots around a cell surface)
f1 = 40 * ff: wi = 6.28 / (js - ja + 1): r = 60
CASE 17'for plots in different brown colours, works only with MS BASIC
isc = 5
col(1, 1) = .3: col(1, isc) = 1
col(2, 1) = 0: col(2, isc) = .8
col(3, 1) = 0: col(3, isc) = .4
FOR ic = 1 TO 3
FOR isci = 1 TO isc
col(ic, isci) = col(ic, 1) + (isci - 1) * (col(ic, isc) - col(ic, 1)) / (isc - 1)
ipat(ic, isci) = 63 * col(ic, isci)
NEXT isci: NEXT ic
FOR i = 1 TO isc
PALETTE isc + 1 - i, 65536 * ipat(3, i) + 256 * ipat(2, i) + ipat(1, i)
NEXT i
CASE 20, 21, 22, 23 ' no special conditions required
CASE ELSE
END SELECT
END IF

' ------ Drawing the curves--------
mess$ = "max. array size for hidden lines reached, use a smaller KROT"
SELECT CASE inz
CASE 1 '3D plot: time = background -> foreground, with hidden lines
  f1 = 40 * ff: fmodus = f1 / krot
  x2 = px * (ja * krot + iy) + xo
  ya = -pyx * ja * krot + pyy * iy + yo
  y2 = ya + a(il, ja) * f1
  j = msiluetm5 * iwhatarray + iy + (ja - 1) * krot
  IF j > msiluetm5 * (iwhatarray + 1) THEN
  CALL cmessage(5, -1, 1, 4, ibackcol, ibackcol, mess$, "")
  EXIT SUB
  END IF

  s2 = msiluet(j)'y-Coordinate of the still visible line
  IF y2 > s2 THEN 'left point is visible
    yb = ya: IF s2 > yb THEN yb = s2
    'PRESET (x2, yb)'Initialtion at the left side
    LINE (x2, yb)-(x2, y2), icolor
    msiluet(j) = y2'new horizon for hidden lines
    ivisible = 1
  END IF
   xanf = x2
   yanf = y2
  FOR ix = ja + 1 TO js
    dy = (a(il, ix) - a(il, ix - 1)) * fmodus - pyx
    xr = 0
    yr = 0
    FOR i = 1 TO krot'------Interpolation-------
      j = j + 1'new x-coordinate of the horizon
      IF j > msiluetm5 * (iwhatarray + 1) THEN
	CALL cmessage(5, -1, 1, 4, ibackcol, ibackcol, mess$, "")
	EXIT SUB
      END IF
      s1 = s2'new y-coordinate of the horizon at point 1
      y1 = y2'old left point becomes the new right point
      x1 = x2
      x2 = x2 + px
      y2 = y1 + dy
      s2 = msiluet(j)'new horizon
      IF y1 > s1 THEN 'IF ivisible > 0 THEN     ' line is not hidden
	IF y2 >= s2 THEN ' point 1 and 2 have been visible
	  xr = xr + px
	  yr = yr + dy
	  msiluet(j) = y2
	ELSE
	  xs = px * (s1 - y1) / (y2 + s1 - y1 - s2 + .000001)' line becomes invisible
	  ys = xs * (y2 - y1) / px'point where kt becomes hidden
	  xr = xr + xs
	  yr = yr + ys
	  'LINE -STEP(xr, yr), icolor
     xend = xanf + xr
     YEND = yanf + yr
     LINE (xanf, yanf)-(xend, YEND), icolor'  line after krot-loop
     xanf = xend
     yanf = YEND
	  xr = 0!
	  yr = 0!
	END IF
      ELSE
	IF y2 > s2 THEN ' point 2 becomes visible
	  xs = px * (s1 - y1) / (y2 + s1 - y1 - s2 + .001)
	  ys = y1 + xs * (y2 - y1) / px
	  xs = x1 + xs
	  xanf = xs
	  yanf = ys
	  'PRESET (xs, ys)
	  xr = x2 - xs
	  yr = y2 - ys
	  msiluet(j) = y2
	END IF
      END IF
    NEXT i'end of krot loop
     IF xr > 0 THEN
     xend = xanf + xr
     YEND = yanf + yr
     LINE (xanf, yanf)-(xend, YEND), icolor'  line after krot-loop
     xanf = xend
     yanf = YEND
     END IF
  NEXT ix' end IX- loop
  ya = -pdyx * (js + ipattern) + pyy * iy + yo
  IF s2 > ya AND s2 < YEND THEN ya = s2
  LINE (xend, YEND)-(xend, ya), icolor

CASE 2'3D plot: time = background -> foreground
  f1 = 40 * ff: fmodus = f1 / krot
  IF iy = ibaseline + 1 THEN 'remove last baseline during continuation
    x2 = px * (ja * krot + iy - 1) + xo
    y2 = pdyx * ja - pyy * (iy - 1) + yo
    LINE (x2, y2 - 1)-(x2 + dxtot, y2 + dytot - 1), ibackcol
  END IF
  x2 = px * (ja * krot + iy) + xo
  ya = pdyx * ja - pyy * iy + yo
  y2 = ya + a(il, ja) * f1
  yma = ya: ym2 = y2
  xm2 = x2
  j = msiluetm5 * iwhatarray + iy + (ja - 1) * krot
  'j = (il - 1) *  * iwhatarray + iy + (ja - 1) * krot
  IF j > msiluetm5 * (iwhatarray + 1) THEN
	CALL cmessage(5, -1, 1, 4, ibackcol, ibackcol, mess$, "")
    EXIT SUB
  END IF
  IF ya < msiluet(j) THEN msiluet(j) = ya
  s2 = msiluet(j)
  
  FOR ix = ja + 1 TO js
    dy = (a(il, ix) - a(il, ix - 1)) * fmodus + pyx
    'Change by difference per step
    ddy = dy / px'+ shift in y-direction
    xmmodus = x2
    ymmodus = y2
    xr = 0
    yr = 0
    FOR i = 1 TO krot'------Interpolation-------
      j = j + 1'new x-coordinate of the horizon
      IF j > msiluetm5 * (iwhatarray + 1) THEN
	CALL cmessage(5, -1, 1, 4, ibackcol, ibackcol, mess$, "")
	EXIT SUB
      END IF
      s1 = s2'new y-coordinate of the horizon at point 1
      y1 = y2'old left point becomes the new right point
      x1 = x2
      x2 = x1 + px
      y2 = y1 + dy
      s2 = msiluet(j)'new horizon
      IF y1 >= s1 THEN ' point 1 has to be removed
	IF y2 > s2 THEN ' point 2 has to be removed too
	  xl = x1
	  yl = y1
	  DO WHILE xl <= x2
	    LINE (xl, yl)-(xl, yma), ibackcol
	    xl = xl + 1
	    yl = yl + ddy
	  LOOP

	ELSE 'remove from point 1 to intersection
	  xs = px * (s1 - y1) / (y2 + s1 - y1 - s2 + .000001)' intersection
	  ys = xs * (y2 - y1) / px'point where kt becomes visible
	  xr = x1 + xs
	  xl = x1
	  yl = y1
	  DO WHILE xl <= xr
	    LINE (xl, yl)-(xl, yma), ibackcol
	    xl = xl + 1
	    yl = yl + ddy
	  LOOP
	END IF
      ELSE ' point 1 was visible
	IF y2 > s2 THEN '     line at point 2 must be removed
	  xs = px * (s1 - y1) / (y2 + s1 - y1 - s2 + .001)'point 2 becomes visible in
	  yr = y1 + xs * (y2 - y1) / px
	  xr = x1 + xs
	  xl = xr
	  yl = yr
	  DO WHILE xl <= x2
	    LINE (xl, yl)-(xl, yma), ibackcol
	    xl = xl + 1
	    yl = yl + ddy
	  LOOP
	END IF
      END IF
      msiluet(j) = y2
    NEXT i' end krot-loop
    LINE (xmmodus, ymmodus)-(x2, y2), icolor
  NEXT ix'' end ix-loop
  msiluet(j) = yma
  xmb = x2
  ymb = yma + pdyx * (js - ja)
  LINE -(x2, ymb), icolor'left border
  LINE (xm2, ym2)-(xm2, yma), icolor'left border

CASE 3 'plot of single curves
type3:
  f1 = 40! * ff: y150 = yo - 150
  ahidden = 0
  ya = a(il, ja) * f1 + yo
  xoo = xo + DX / 2
IF icolor > 19 GOTO stepcurves
  xa = xoo
  IF il = 0 THEN   'special for source density = dots
    FOR i = ja TO js
      ya = a(il, i) * f1 + yo
      LINE (xa, ya)-(xa + iwi(il), ya + iwi(il)), icolor, BF
      xa = xa + DX
    NEXT i
  ELSE
    FOR i = ja + 1 TO js
      xb = xa + DX
      yb = a(il, i) * f1 + yo
      IF yb < y150 THEN yb = y150
      CALL linep(ihp, xa, ya, xb, yb, 0, iwi(il), icolor, ityp(il))
      xa = xb: ya = yb
    NEXT i
  END IF
EXIT SUB
stepcurves:
icol20 = icolor - 20
iwi2 = iwi(il) / 2
  ya = a(il, ja) * f1 + yo
  xa = xoo - .5 * DX
  xb = xa + DX
      LINE (xa + iwi2, ya)-(xb - iwi2, ya + iwi(il)), icol20, BF
    FOR i = ja + 1 TO js
      xa = xb - iwi
      yb = a(il, i) * f1 + yo
      IF yb > ya THEN
      LINE (xa - iwi2, ya)-(xa + iwi2, yb + iwi(il)), icol20, BF
      ELSE
      LINE (xa - iwi2, ya + iwi(il))-(xa + iwi2, yb), icol20, BF
      END IF
      xb = xa + DX
      LINE (xa + iwi2, yb)-(xb - iwi2, yb + iwi(il)), icol20, BF
      ya = yb
    NEXT i


CASE 4, 5'plot with symbols
  amx = .00000001#
  FOR ix = ja TO js
    IF a(il, ix) > amx THEN amx = a(il, ix)
  NEXT
  PRINT #8, tv$(il); : PRINT #8, USING "#####.###"; amx;
  am10 = 8 / amx: dummy$ = " º"
  FOR ix = ja TO js
    nx = a(il, ix) * am10 + 1
    IF nx < 0 THEN nx = 0
    dummy$ = dummy$ + sp$(nx)
  NEXT ix
    dummy$ = dummy$ + "º"
  PRINT #8, MID$(dummy$, 1, 69)


'----------vertical pixel pattern (especially for shell patterns)----
CASE 6, 8, 9, 11, 12, 13, 14, 17, 22
  x1 = xo: idelta = iy: IF inz = 22 THEN idelta = -iy
  y1 = yo - idelta + 1
  IF y1 < 2 THEN EXIT SUB
  LINE (x1 - 2, y1)-(x1 - 1, y1), ivorcol
  SELECT CASE inz
  CASE 6, 8'false colour
    FOR ix = ja TO js
      jlocal& = a(il, ix) * f10
      IF jlocal& > 10 THEN jlocal& = 10
      IF jlocal& >= krot THEN LINE (x1, y1)-(x1 + DX, y1), icolf(jlocal&)
      x1 = x1 + DX
    NEXT ix
  CASE 9, 11, 12, 13, 22
    krnd = RND * 3 + 1
    FOR ix = ja TO js
      jlocal& = a(il, ix) * ff
      IF jlocal& >= 1 THEN
	'IF jlocal& >= krot THEN
	IF jlocal& > 4 THEN jlocal& = 4
	LINE (x1, y1)-(x1 + DX, y1), icolor, , ipat(jlocal&, krnd)
	krnd = krnd + 1: IF krnd > 4 THEN krnd = 1
      END IF
      x1 = x1 + DX
    NEXT ix
  CASE 14 'full line
    FOR ix = ja TO js
      jlocal& = a(il, ix) * ff
      IF jlocal& > 10 THEN jlocal& = 10
      IF jlocal& >= 1 THEN LINE (x1, y1)-(x1 + DX, y1), icolor
      x1 = x1 + DX
    NEXT ix
  REM
CASE 17 'brownish colour
    FOR ix = ja TO js
      jlocal& = a(il, ix) * ff
      IF jlocal& > 5 THEN jlocal& = 5
      IF jlocal& >= 1 THEN LINE (x1, y1)-(x1 + DX, y1), jlocal&
      x1 = x1 + DX
    NEXT ix
  REM
  END SELECT
  LINE (x1 + 1, y1)-(x1 + 2, y1), ivorcol

CASE 7, 10 'landscape orientation
  y1 = yo
  x1 = xo + iy - 1
  IF x1 > 636 THEN EXIT SUB
  LINE (x1, y1 - 2)-(x1, y1 - 1), ivorcol
  SELECT CASE inz
  CASE 7 'false colour
    FOR ix = ja TO js
      jlocal& = a(il, ix) * f10
      IF jlocal& > 10 THEN jlocal& = 10
      IF jlocal& >= krot THEN LINE (x1, y1)-(x1, y1 + DX), icolf(jlocal&)
      y1 = y1 + DX
    '  x1 = x1 + dx
    NEXT ix
  CASE ELSE'pixel
    krnd = RND * 3 + 1
    FOR ix = ja TO js
      jlocal& = a(il, ix) * ff
      IF jlocal& >= 1 THEN
	IF jlocal& > 4 THEN jlocal& = 4
	LINE (x1, y1)-(x1, y1 + DX), icolor, , ipat(jlocal&, krnd)
	krnd = krnd + 1: IF krnd > 4 THEN krnd = 1
      END IF
    y1 = y1 + DX
  NEXT ix
  END SELECT
  LINE (x1, y1)-(x1, y1 + 1), ivorcol

CASE 15'-drawing on circles, Nautilus
  'wi controls the angular displacement
  krnd = RND * 3 + 1
  FOR ix = ja TO js
    jlocal& = a(il, ix) * ff
    IF jlocal& >= 1 THEN
      IF jlocal& > 4 THEN jlocal& = 4
      FOR ddx = 0 TO DX STEP .3
      x1 = 320 + (iradius + (DX * ix + ddx)) * COS(alpha)
      y1 = 240 + (iradius + (DX * ix + ddx)) * SIN(alpha)
      x2 = 320 + (iradius + (DX * ix + ddx)) * COS(alpha + wi)
      y2 = 240 + (iradius + (DX * ix + ddx)) * SIN(alpha + wi)
      LINE (x1, y1)-(x2, y2), icolor
      NEXT ddx
GOTO auweiter
'ELSE
      x1 = 320 + DX * (iradius + ix + ddx) * COS(alpha)
      y1 = 240 + DX * (iradius + ix + ddx) * SIN(alpha)
      x2 = 320 + DX * (iradius + ix + DX) * COS(alpha)
      y2 = 240 + DX * (iradius + ix + ddx) * SIN(alpha)
      LINE (x1, y1)-(x2, y2), icolor
'      END IF
auweiter:
    END IF
  NEXT ix
  alpha = alpha + .5 * wi
CASE 16
xo = 320: yo = 240: r = 60
wi1 = (ja - 1) * wi
ixa = xo + (r + a(il, ja) * f1) * COS(wi1)
iya = yo + per * (r + a(il, ja) * f1) * SIN(wi1)
PSET (ixa, iya), 6 + il
FOR ix = ja + 1 TO js
wi1 = (ix - 1) * wi
ixb = xo + (r + a(il, ix) * f1) * COS(wi1)
iyb = yo + per * (r + a(il, ix) * f1) * SIN(wi1)
LINE -(ixb, iyb), 6 + il
NEXT
LINE -(ixa, iya), 6 + il
PAINT (xo, yo), icolor, 6 + il

CASE 18, 19
      xoo = xo + DX / 2
      IF icolor > 39 THEN
      icolor = icolor - 40
      LINE (xoo, yo)-(xoo + (js - ja) * DX, yo), icolor
      GOTO type3
      END IF
      IF icolor > 19 GOTO rectangles
      yoo = yo + 3: IF f1 < 0 THEN yoo = yo - 3
      FOR ibordercol = 7 TO 13 STEP 6'first gray, then pink
      xa = xoo
      ya = yoo + a(il, ja) * f1
      LINE (xa, yo)-(xa, ya), ibordercol
      FOR ix = ja + 1 TO js
      yb = yoo + a(il, ix) * f1
      xb = xa + DX
      LINE (xa, ya)-(xb, yb), ibordercol
      xa = xb
      ya = yb
      NEXT ix
      'ybaseline=
      LINE (xa, ya)-(xa, yo), ibordercol
      IF ibordercol = 7 THEN
      LINE (xa, yo)-(xoo, yo), ibordercol
      PAINT (xoo + 2, yo + (yoo - yo) / 2), icolor, 7
      'PAINT (xoo + 2, yoo - yo), icolor, 7
      END IF
      NEXT
      EXIT SUB
rectangles:
      icolor = icolor - 20
  ya = a(il, ja) * f1 + yo
  xa = xoo - .5 * DX
    FOR i = ja TO js
      xb = xa + DX
      ya = a(il, i) * f1 + yo
      LINE (xa, yo)-(xb, ya), icolor, BF
      xa = xb
    NEXT i

CASE 20
  f1 = 40! * ff
  ya = a(il, ja) * f1 + yo
  IF iy = 1 THEN
    PRESET (xo, ya), icolor
    ym(il) = ya
  ELSE
    LINE (xo - 1, ym(il) + 1)-(xo, ya), icolor, BF
  END IF
  ym(il) = ya
  LINE (xo, yo - 1)-(xo, yo - 2), ivorcol'Basis line

CASE 21'separate rectangles
  f1 = 40! * ff
  xa = xo + .5 * DX
    FOR i = ja TO js
      xb = xa + DX * .8
      ya = a(il, i) * f1 + yo
      CALL linep(ihp, xa, yo, xb, ya, 2, 1, i, 1)
      xa = xa + DX
    NEXT i
END SELECT

EXIT SUB

profileframe: '-------------finishing frame after simulation
ahidden = 0
SELECT CASE inz
CASE 1
  x1 = px * ((ja + ipattern) * krot + 1) + xo
  y1 = -pdyx * (ja + ipattern) + pyy + yo
  x2 = x1 + pdx * (js - ja)
  y2 = y1 - pdyx * (js - ja)
  LINE (x1, y1)-(x2, y2), icolor
  x3 = x2 + (iy - 1) * px
  y3 = y2 + pyy * (iy - 1)
  LINE -(x3, y3), icolor
  x4 = x3 - pdx * (js - ja)
  y4 = y3 + pdyx * (js - ja)
  LINE -(x4, y4), icolor, , &H1010
  LINE -(x1, y1), icolor, , &H1010
  yver = y4 + a(il, ja) * f1
  LINE (x4, y4)-(x4, yver), icolor, , &H1010

CASE 2
  ja = ja + ipattern
  x1 = px * (ja * krot + 1) + xo
  y1 = pdyx * ja - pyy + yo
  x2 = px * (ja * krot + iy) + xo
  y2 = pdyx * ja - pyy * iy + yo
  x3 = x2 + (js - ja) * px * krot
  y3 = y2 + pdyx * (js - ja)
  LINE (x3, y3)-(x3 - x2 + x1, y3 - y2 + y1), icolor, , &H1010
  LINE -(x1, y1), icolor, , &H1010
  LINE -(x2, y2), icolor
  LINE (x2, y2 - 1)-(x3, y3 - 1), icolor
  dxtot = x3 - x2
  dytot = y3 - y2
  ibaseline = iy


CASE 6, 8, 9, 11, 12, 13, 14, 17
  LINE (xo - 2, y1 - 1)-(xo + DX * (js - ja + 1) + 2, y1 - 1), ivorcol
  LINE (xo - 2, y1 - 2)-(xo + DX * (js - ja + 1) + 2, y1 - 2), ivorcol
CASE 7, 10
  LINE (x1 + 1, yo - 2)-(x1 + 1, yo + (js - ja + 1) * DX + 1), ivorcol
  LINE (x1, yo - 2)-(x1, yo + (js - ja + 1) * DX + 1), ivorcol
CASE 22
  LINE (xo - 2, y1 + 1)-(xo + DX * (js - ja + 1) + 2, y1 + 1), ivorcol
  LINE (xo - 2, y1 + 2)-(xo + DX * (js - ja + 1) + 2, y1 + 2), ivorcol
END SELECT
EXIT SUB
END SUB

SUB showprof (inz, iendm, dw$, a(), ja, js, icc) STATIC
'calls the prof(ile)- subroutine for each substance
'inz=type of the plot, iendm=99 to finish the plot
'dw$: substance to be displayed
'a(il,ix) two-dimensional array, substance and position
'ja, js left,right cell, icc counting index, icc=1:initiation
'
SHARED KX, ky, KD, KN, KG, icol, K1, K2, K3, K4
SHARED DX, dy, ip$, krot, per, iparam, iactive
SHARED xa(), fa(), am()
SHARED ya(), f1(), icol(), ityp(), iwi(), msiluet(), msiluetmax
SHARED itot&, LABEL$(), displaywas$, icheckerror
SHARED ihp, dxp, dyp, fkasten, ivgashift, ianz, ireaddata
SHARED ivorcol, ibackcol, igraph, istelz, iletztk, ihintr
ON ERROR GOTO 0

IF ja = 0 THEN
   mess$ = "left cell is zero, run program first"
	resp$ = "ok"
	CALL cmessage(4, -1, 1, 15, 4, ibackcol, mess$, resp$)
   itot& = 99999
   EXIT SUB
END IF

SELECT CASE inz

CASE 1, 2
  CALL prof(a(), 2, dw$, 0, inz, iendm, DX, fa(2), ja, js, icc, xa(2), ya(2))
  CALL prof(a(), 3, dw$, 0, inz, iendm, DX, fa(3), ja, js, icc, xa(3), ya(3))
  CALL prof(a(), 1, dw$, 0, inz, iendm, DX, fa(1), ja, js, icc, xa(1), ya(1))

CASE 3
  IF iendm = 99 THEN EXIT SUB
  IF icc = 1 OR KG > 0 THEN
  y2 = 440: yo = ya(1)
  jam = ja:
  xa = 320 - ((js - ja + 1) / 2) * DX
  yo = 120
  END IF
  iwi2 = iwi(1) / 2
  y3 = yo - iwi(1): y4 = yo - iwi(1) - iwi2
  xo = xa + DX * (ja - jam) 'to have correct position after change of field
LINE (1, 1)-(650, 480), ibackcol, BF '= CLS
  jja = ja: jjs = js: GOSUB type3frame
GOSUB drawnow

CASE 4, 5'----- schematic plot with letters --------
  iparm = 1
  IF iendm = 99 THEN EXIT SUB
  IF icc = 1 THEN CLS
  LOCATE 1, 1
  OPEN "SCRN:" FOR OUTPUT AS #8
  GOSUB showprint45
  IF inz = 5 THEN
  CALL lptopen("LPT1:"): IF icheckerror = 1 THEN EXIT SUB
  GOSUB showprint45
  END IF
GOTO esubb

showprint45:
  PRINT #8, "Gr”áe"; : PRINT #8, USING "##"; ja; : PRINT #8, "-"; : PRINT #8, USING "###"; js;
  dummy$ = "É" + STRING$(js - ja + 1, "Í") + "»"
  PRINT #8, MID$(dummy$, 1, 68)
  FOR il = 0 TO 8'kn
   CALL prof(a(), il, dw$, 0, inz, iendm, DX, gB, ja, js, 1, xo, yo)
  NEXT il
  PRINT #8, "IT"; : PRINT #8, USING "########"; itot&;
  dummy$ = " È" + STRING$(js - ja + 1, "Í") + "¼"
  IF (js - ja) > 65 THEN dummy$ = " È" + STRING$(34, "Í") + " Das ist nur der linke Teil...   "
  PRINT #8, MID$(dummy$, 1, 68)
  PRINT #8, : PRINT #8,
  CLOSE (8)
  RETURN


CASE 6 TO 14, 17, 22'----- pixel-plot----------12,13 print in addition
  '22 t-axis upwards, for phyllotaxis
  IF icc = 1 THEN
    SELECT CASE igraph'units=screen-pixel
    CASE 2, 3
      xmax = 720: ymax = 348
    CASE 9
      xmax = 640: ymax = 350
    CASE 12
      xmax = 640: ymax = 480
    END SELECT
    WINDOW (1, 1)-(xmax, ymax)
    END IF
  FOR ilv = LEN(displaywas$) TO 1 STEP -1
    il = ASC(MID$(displaywas$, ilv, 1)) - 96'a=asc(97)
    IF il = 19 THEN il = 0 's=source
    IF icc = 1 THEN
    y1 = ymax - 4: IF inz = 22 THEN y1 = 25
    x1 = 1 + DX * ja + xa(il)
    x2 = 1 + DX * (js + 1) + xa(il)
    SELECT CASE KD
    CASE 7, 10
      x1 = 4
      y1 = 35 + DX * ja
      y2 = 35 + DX * (js + 1)
      LINE (x1 - 2, y1 - 2)-(x1 - 2, y2 + 1), ivorcol
      LINE (x1 - 1, y1 - 2)-(x1 - 1, y2 + 1), ivorcol

    CASE 8, 11 'Plot of the spatial stable pattern
      y1 = y1 - 50
      LINE (x1 - 1, y1)-(x2 + 1, y1 + 50), ivorcol, B
      LINE (x1 - 2, y1 - 1)-(x2 + 2, y1 + 51), ivorcol, B
      amz = 0
      amz = 45: xm = x1 + .5 * DX: ym = y1 + a(8, ja) * amz
      FOR i = ja + 1 TO js
	xz = xm + DX: yz = y1 + a(8, i) * amz
	'LINE (xz - 1, yz - 1)-(xm + 1, ym + 1), ivorcol, BF
	LINE (xz, yz)-(xm, ym), ivorcol
	LINE (xz + 1, yz)-(xm + 1, ym), ivorcol
	LINE (xz + 1, yz + .9)-(xm + 1, ym + .9), ivorcol
	LINE (xz + 1, yz + 1.5)-(xm + 1, ym + 1.5), ivorcol
	xm = xz: ym = yz
      NEXT i
      y1 = y1 - 10
      LINE (x1 - 2, y1 + 2)-(x2 + 2, y1 + 2), ivorcol
      LINE (x1 - 2, y1 + 1)-(x2 + 2, y1 + 1), ivorcol
    CASE ELSE
      LINE (x1 - 2, y1 + 2)-(x2 + 2, y1 + 2), ivorcol
      LINE (x1 - 2, y1 + 1)-(x2 + 2, y1 + 1), ivorcol
    END SELECT
  END IF

  
  IF KD = 7 OR KD = 10 THEN
    y1 = 35 + DX * ja
  ELSE
    x1 = 1 + DX * ja + xa(il)
  END IF
    CALL prof(a(), il, dw$, kv, inz, iendm, DX, fa(il), ja, js, icc, x1, y1)
  NEXT ilv

  IF inz < 12 OR inz > 13 THEN EXIT SUB
  IF iopen12 = 0 THEN
    OPEN "pixel.tmp" FOR OUTPUT AS #12
    iopen12 = 1
  END IF

  IF iendm < 99 THEN 'Storage of the concentrations for later printout
    resp$ = ""
    SELECT CASE inz
    CASE 12
      FOR ix = ja TO js
	temp$ = "0"
	jlocal& = a(1, ix) * fa(1)
	IF jlocal& > 5 THEN jlocal& = 5
	IF jlocal& > 0 THEN temp$ = MID$(STR$(jlocal&), 2, 1)
	resp$ = resp$ + temp$
      NEXT ix
    CASE 13'array a and c are plotted with thresholds
      FOR ix = ja TO js
	temp$ = "0"
	IF a(3, ix) * fa(3) > 2 THEN temp$ = "1"
	IF a(1, ix) * fa(1) > 2 THEN temp$ = "4"
	resp$ = resp$ + temp$
      NEXT ix
    END SELECT
    PRINT #12, resp$
    ELSE 'final plott, iendm=99
    CLOSE #12: iopen12 = 0
    mess$ = "You like to print this picture on a HP printer, y or RETURN"
    CALL zeingabe(igt, 0, 4, i, dummy, resp$, mess$)
    IF resp$ = "Y" OR resp$ = "y" OR resp$ = "j" THEN
     ''the NEC needle printer is no longer supported
	OPEN "pixel.tmp" FOR INPUT AS #7
	CALL lptopen("LPT1:"): IF icheckerror = 1 THEN CLOSE : EXIT SUB
	WIDTH #8, 255
	CALL cmessage(3, -1, 1, 15, 2, 0, "... es wird gedruckt ...", "")
	PRINT #8, CHR$(27); "&a5L" 'left border
	PRINT #8, CHR$(27); "&k4S"; ' Elite
	PRINT #8, CHR$(27); "*t150R";
	ianzahl = (js - ja + 1) / 2
	anzahl$ = "*b" + MID$(STR$(ianzahl + 2), 2, 255) + "W"
	PRINT #8, CHR$(27); "*rA";
	PRINT #8, CHR$(27); anzahl$; CHR$(3);
	FOR ix = 1 TO ianzahl: PRINT #8, CHR$(255); : NEXT ix
	PRINT #8, CHR$(192);

	DO UNTIL EOF(7)
	  INPUT #7, resp$
	  PRINT #8, CHR$(27); anzahl$; CHR$(3);
	  FOR ixl = 0 TO ianzahl - 1
	    ix = ja + 2 * ixl
	    jlocal& = VAL(MID$(resp$, ix, 1))
	    IF jlocal& > 4 THEN jlocal& = 4
	    jlocal& = jlocal& + id
	    SELECT CASE jlocal&'pixel for HP-CL
	    CASE 0, 5, 10, 15, 6, 16: iwas = 0
	    CASE 1: iwas = 16
	    CASE 11: iwas = 64
	    CASE 2, 12: iwas = 128
	    CASE 7, 17: iwas = 32
	    CASE 3, 13: iwas = 160
	    CASE 8, 18: iwas = 80
	    CASE 4, 9, 14, 19: iwas = 240
	    END SELECT
	    jlocal& = VAL(MID$(resp$, ix, 1))
	    IF jlocal& > 4 THEN jlocal& = 4
	    jlocal& = jlocal& + id
	    SELECT CASE jlocal&
	    CASE 0, 1, 5, 10, 11, 15:
	    CASE 6: iwas = iwas + 1
	    CASE 16: iwas = iwas + 4
	    CASE 2, 12: iwas = iwas + 8
	    CASE 7, 17: iwas = iwas + 2
	    CASE 3, 13: iwas = iwas + 10
	    CASE 8, 18: iwas = iwas + 5
	    CASE 4, 9, 14, 19: iwas = iwas + 15
	    CASE ELSE
	      PRINT jlocal&
	      INPUT resp$
	    END SELECT
	    PRINT #8, CHR$(iwas);
	  NEXT ixl
	  PRINT #8, CHR$(192);
	  id = id + 5
	  IF id = 20 THEN id = 0
	IF INKEY$ = CHR$(27) THEN EXIT DO
	LOOP
	PRINT #8, CHR$(27); anzahl$; CHR$(3);  'final line
	FOR ix = 1 TO ianzahl: PRINT #8, CHR$(255); : NEXT ix
	PRINT #8, CHR$(192); CHR$(27); "*rB": PRINT #8, : PRINT #8,
	CLOSE #8
	CLOSE #7
	CALL cmessage(6, -1, 1, 15, 2, 0, "eject page with command FF, if you like", "")
  END IF
END IF
CASE 15
    il = 1
    CALL prof(a(), il, dw$, 0, inz, iendm, DX, fa(il), ja, js, icc, xo, yo)

CASE 16
IF iendm = 99 THEN EXIT SUB
WINDOW (1, 1)-(640, 480)
LINE (1, 1)-(650, 480), ibackcol, BF 'erase frame
GOSUB drawnow
r = 60 '10! * dx
wi = 6.28 / (js - ja + 1)
FOR ix = ja TO js
wi1 = wi * ix
xo = 320: yo = 240
xa = xo + r * COS(wi1)
ya = yo + per * r * SIN(wi1)
xm = xo + .8 * r * COS(wi1)
ym = yo + per * .8 * r * SIN(wi1)
CALL linep(ihp, xm, ym, xa, ya, 0, iwi(1), 1, 1)
xb = xo + r * COS(wi1 + wi)
yb = yo + per * r * SIN(wi1 + wi)
CALL linep(ihp, xa, ya, xb, yb, 0, iwi(1), 1, 1)
NEXT ix
PAINT (xo, yo), 15, 1
IF K2 > 0 AND dy > 0 THEN 'arrow
wff = 6.28 * (K2 - 1) / (js - ja + 1)
aw = 3! * r
x = xo + aw * COS(wff)
y = yo + per * aw * SIN(wff)
xa = xo + aw * COS(wff) * .6
ya = yo + per * aw * SIN(wff) * .6
CALL linep(ihp, x, y, xa, ya, 0, iwi(1), ivorcol, 1)
x = xo + aw * COS(wff + .03) * .69
y = yo + per * aw * SIN(wff + .03) * .69
CALL linep(ihp, x, y, xa, ya, 0, iwi(1), ivorcol, 1)
x = xo + aw * COS(wff - .03) * .69
y = yo + per * aw * SIN(wff - .03) * .69
CALL linep(ihp, x, y, xa, ya, 0, iwi(1), ivorcol, 1)
END IF

CASE 18, 19 'Individual plots, area filled, 19 with new screen
  IF iendm = 99 THEN EXIT SUB
  IF inz = 19 OR xa(1) = 0 THEN
  IF icc = 1 OR KG > 0 THEN
  xa = 320 - ((js - ja + 1) / 2) * DX
  jam = ja
  END IF
  xo = xa + DX * (ja - jam) 'to have correct position after change of field
  END IF
  LINE (1, 1)-(650, 480), ibackcol, BF 'wipe screen
GOSUB drawnow

CASE 20'for display of oscillation within a single cell
  IF icc = 1 THEN
    yo = ya(1): y2 = yo + dy: xo = 10
    LINE (xo, y2)-(xo, yo - 2), ivorcol
    xo = xo + 1: LINE (xo, y2)-(xo, yo - 2), ivorcol
  ELSE
    xo = xo + 1
    IF xo > 630 THEN
      yo = yo - dy - 30: y2 = yo + dy: xo = 10: icc = 1
      LINE (xo, y2)-(xo, yo - 2), ivorcol
      xo = xo + 1: LINE (xo, y2)-(xo, yo - 2), ivorcol
    END IF
  END IF
  GOSUB drawnow

CASE 21' cells plotted as separate squares, for gene activation
  IF iendm = 99 THEN EXIT SUB
  LINE (1, 1)-(650, 480), ibackcol, BF 'wipe screen
  xo3 = 320 - (js + 1) / 2 * DX
  xo = xo3 + (ja - 1) * DX
  x2 = xo + (js - ja + 2) * DX
    yo = ya(3): y2 = yo + a(3, 1) * fa(3) * 40
    LINE (xo, yo)-(x2, y2), icol(3), BF
    CALL prof(a(), 1, dw$, 0, inz, iendm, DX, fa(1), ja, js, icc, xo, ya(1))
    CALL prof(a(), 2, dw$, 0, inz, iendm, DX, fa(2), ja, js, icc, xo, ya(2))



CASE 23  'two plots side by side (for regeneration)
  IF iendm = 99 THEN EXIT SUB
  IF icc = 1 OR KG > 0 THEN
  yo = 120: y2 = 440: iwi2 = iwi(1) / 2
  y3 = yo - iwi(1): y4 = yo - iwi(1) - iwi2
    END IF
  LINE (1, 1)-(640, 480), ibackcol, BF 'erease frame
  xa = 320 - ((js - ja + 1) / 2) * DX
  jja = ja: jjs = K4 - 1: xo = xa + (jja - ja) * DX: GOSUB type3frame
  FOR il = KN TO 0 STEP -1
    CALL prof(a(), il, dw$, 0, 3, iendm, DX, fa(il), ja, K4 - 1, icc, xo, yo)
  NEXT il
jja = K4 + 1: jjs = js: xo = xa + (jja - ja) * DX: GOSUB type3frame
  FOR il = KN TO 0 STEP -1
    CALL prof(a(), il, dw$, 0, 3, iendm, DX, fa(il), K4 + 1, js, icc, xo, yo)
  NEXT il


CASE ELSE
mess$ = "for this kd no graphic display is yet installed"
resp$ = "ok"
CALL cmessage(5, -1, 1, 15, 4, 0, mess$, resp$): EXIT SUB
END SELECT

esubb:
IF igraph = 9 THEN
IF LABEL$(0) > "" THEN
CALL cmessage(5, -1, 1, 1, 15, ibackcol, LABEL$(0), "")
END IF

SELECT CASE KD
CASE 3, 4, 5, 16, 18, 19, 21, 22, 23
CALL switchscreen(2, iactive, ibackcol)

END SELECT
END IF
EXIT SUB

type3frame:
  x2 = xo + (jjs - jja + 1) * DX
    CALL linep(ihp, xo, y4, xo, y2 + iwi2, 0, iwi(1), ivorcol, 1)
    CALL linep(ihp, xo, y3 - 1, x2, y3 - 1, 0, iwi(1), ivorcol, 1)
    CALL linep(ihp, x2, y4, x2, y2 + 2, 0, iwi(1), ivorcol, 1)
    x3 = xo
    FOR i = jja TO jjs + 1
      CALL linep(ihp, x3, y4, x3, y4 - 5, 0, iwi(1), ivorcol, 1)
      x3 = x3 + DX
    NEXT i
RETURN
drawnow:
  FOR ilv = LEN(displaywas$) TO 1 STEP -1
    il = ASC(MID$(displaywas$, ilv, 1)) - 96'a=asc(97)
    IF il = 19 THEN il = 0
    IF inz = 18 OR inz = 19 OR inz = 22 THEN
    CALL prof(a(), il, dw$, 0, inz, iendm, DX, fa(il), ja, js, icc, xo, ya(il))
    ELSE
    CALL prof(a(), il, dw$, 0, inz, iendm, DX, fa(il), ja, js, icc, xo, yo)
    END IF
  NEXT ilv
RETURN
END SUB


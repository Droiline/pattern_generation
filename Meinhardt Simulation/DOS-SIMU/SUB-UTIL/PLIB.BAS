DECLARE SUB cmessage (itype%, iline%, irow%, ilettercol%, iletterback%, ibackcol%, text$, resp$)
DECLARE SUB eingabe (igt%, inz%, iaf%(), f!(), LABEL$(), text$, ianf%, ianz%)
DECLARE SUB gmanage (wasnun$)
DECLARE SUB gtourcontrol (PCONTROL$, igt%, LABEL$())
DECLARE SUB gtreturn (inz%, igt%)
DECLARE SUB init (inewscreen%, igraph%, ivorcol%, ibackcol%, ivgashift%)
DECLARE SUB lptopen (lpt$)
DECLARE SUB prtscreen (scrfile$)
DECLARE SUB rwparameterdisk (rw$, ifilenumber%)
DECLARE SUB SCREENCHECK (igraph%)
DECLARE SUB smenu (irfilename$, iline%, headline$)
DECLARE SUB sorts (s$(), icount%)
DECLARE SUB switchscreen (inz%, iactive%, ibackcol%)
DECLARE SUB titles (inz%, title$, mess$)
DECLARE SUB whatprmfiles (programname$, progrtmp$, exten$, drive$, irfilename$)
DECLARE SUB zeingabe (igt%, iquest%, inz%, i%, f!, LABEL$, text$)
DECLARE SUB ColorPrint (text$, row%, col%, fc%, bc%)
'Subroutines for the simulation of pigmentation patterns on sea shells
' (c) by Hans Meinhardt, 1994

DEFINT H-N
' $INCLUDE: 'sub-util\scommon.bi' nur fr Microsoft-Basic gebraucht
' $INCLUDE: 'sub-util\errlist.bas' nur fr Microsoft-Basic gebraucht

SUB cmessage (iwhat, iline, irow, ilettercol, iletterback, ibackcol, text$, resp$)
SHARED ivgashift
' Prints a message on the screen, >3 centered
' iwhat = 2 or 4: a RETURN is required
' iline: The line on which the message appears, VGA-Corrected
' if iline <0 'counted from lower border
ilinea = iline: IF iline < 0 THEN ilinea = 26 + ivgashift + iline
IF ilinea < ivgashift + 25 THEN linea = ivgashift + 25
IF iwhat > 2 THEN
iblank = 81: IF ilinea = 25 + ivgashift THEN iblank = 80
CALL ColorPrint(SPACE$(iblank - irow), ilinea, 1, ilettercol, iletterback)
END IF
IF irow > 79 THEN irow = 79
tex$ = LEFT$(text$, 80 - irow)
IF iwhat > 3 THEN 'centered
   isp = LEN(tex$)
   IF isp + irow < 72 THEN
   isp = (79 - isp) / 2:
   tex$ = LEFT$(SPACE$(isp) + tex$, 80)
   END IF
   END IF
CALL ColorPrint(tex$, ilinea, irow, ilettercol, iletterback)
  SELECT CASE iwhat
   CASE 2, 4
   IF resp$ > "" THEN
   CALL ColorPrint("[" + resp$ + "]", ilinea, 78 - LEN(resp$), 15, 13)
   END IF
   resp$ = ""
   DO UNTIL resp$ > "": resp$ = INKEY$: LOOP
   CALL ColorPrint(SPACE$(79), ilinea, 1, ibackcol, ibackcol)
   CASE 5
   IF resp$ > "" THEN
   CALL ColorPrint("[" + resp$ + "]", ilinea, 78 - LEN(resp$), 15, 13)
   END IF
   CASE 6, 7
   jdelay = 2: IF iwhat = 7 THEN jdelay = .6
   resp$ = ""
   t1 = TIMER: DO UNTIL TIMER - t1 > jdelay: IF INKEY$ > "" THEN EXIT DO:
   LOOP
      CALL ColorPrint(SPACE$(79), ilinea, 1, ibackcol, ibackcol)
   END SELECT
END SUB

SUB eingabe (igt, inz, iaf(), f(), LABEL$(), text$, ianf, ianz)
SHARED ivgashift, displaywas$, igtdelay, ibackcol
'Input routine for mupltiple changes of the plot parameters
CALL cmessage(1, -1, 1, 15, 15, ibackcol, SPACE$(79), resp$)
twas$ = "sabcdefgh"
FOR i = ianf TO ianz
  ipu = 1: ipu = INSTR(ipu, displaywas$, MID$(twas$, i + 1, 1))
  IF ipu > 0 THEN
    textprint$ = "Substance  : " + text$ + "  OLD = "
  IF inz = 1 THEN textprint$ = textprint$ + STR$(iaf(i))
  IF inz = 2 THEN textprint$ = textprint$ + STR$(f(i))
  IF inz = 3 THEN textprint$ = textprint$ + LABEL$(i)
  textprint$ = textprint$ + "  NEW = "
CALL cmessage(1, -1, 1, 2, 15, ibackcol, textprint$, resp$)
CALL cmessage(1, -1, 11, 4, 15, ibackcol, MID$(twas$, i + 1, 1), resp$)
ilength = LEN(textprint$) + 1
fac$ = ""
CALL cmessage(1, -1, ilength, 12, 15, ibackcol, fac$ + "_ ", resp$)

facinput:
SELECT CASE igt
CASE 0, -1 'No GUIDED TOUR active
IF inz = 3 THEN
LINE INPUT ; " "; fac$
ELSE
fac$ = ""
DO
dummy$ = INKEY$
REM
SELECT CASE dummy$
CASE ""
CASE "g", "G", "-", "+", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "."
fac$ = fac$ + dummy$
CALL cmessage(1, -1, ilength, 12, 15, ibackcol, fac$ + "_  ", resp$)
CASE CHR$(8)
lentext = LEN(fac$)
CALL cmessage(1, -1, ilength, 12, 15, ibackcol, fac$ + "_  ", resp$)
IF LEN(fac$) > 0 THEN fac$ = LEFT$(fac$, LEN(fac$) - 1)
CALL cmessage(1, -1, ilength, 12, 15, ibackcol, fac$ + "_  ", resp$)
CASE "x", CHR$(27):
CALL cmessage(1, -1, 1, 12, 15, ibackcol, SPACE$(79), "")
GOTO nexti
CASE "q": GOTO einende
CASE CHR$(13):
CALL cmessage(1, -1, 1, 12, 15, ibackcol, SPACE$(79), "")
EXIT DO
CASE ELSE: BEEP
END SELECT
LOOP
END IF

CASE 1 TO 1000: LINE INPUT #11, fac$
IF igtdelay > 1 THEN CALL gtreturn(5, igt)
CALL cmessage(1, -1, ilength, 12, 15, ibackcol, fac$ + "_    ", resp$)
IF igtdelay > 1 THEN CALL gtreturn(5, igt)
IF fac$ = "q" GOTO einende
IF fac$ = "x" GOTO nexti
END SELECT
    IF LEFT$(fac$, 1) = "g" OR LEFT$(fac$, 1) = "G" THEN
      fac$ = MID$(fac$, 2, 20)
      FOR il = i TO ianz
	IF inz = 1 THEN iaf(il) = VAL(fac$)'      Global change
	IF inz = 2 THEN f(il) = VAL(fac$)
	IF inz = 3 THEN LABEL$(il) = fac$
      NEXT il
      GOTO einende
    END IF

IF inz = 1 THEN iaf(i) = VAL(fac$)
IF inz = 2 THEN f(i) = VAL(fac$)
IF inz = 3 THEN LABEL$(i) = fac$
nexti:
  END IF
NEXT i
einende:
CALL cmessage(1, -1, 1, 15, 15, ibackcol, SPACE$(79), resp$)
END SUB

SUB gmanage (wasnun$) STATIC
'for changes of the graphic parameters
SHARED flv(), inv(), igt
SHARED xa(), fa(), am(), ya(), f1(), icol()
SHARED ityp(), iwi(), LABEL$()
SHARED KN, K4
SHARED ihp, dxp, dyp, fkasten, DX, per, krot, ivgashift
SHARED ivorcol, ibackcol, igraph
ianf = VAL(RIGHT$(wasnun$, 1))
IF ianf > 0 THEN
  lenwas = LEN(wasnun$)
  wasnun$ = MID$(wasnun$, 1, lenwas - 1)
END IF

SELECT CASE wasnun$
CASE "PER": CALL zeingabe(igt, 1, 2, i, per, resp$, "Perspective")
CASE "KROT": CALL zeingabe(igt, 1, 1, krot, D, resp$, "kROTATION")
IF krot > 7 THEN
    mess$ = "KROT was too large, it is changed to 7"
    CALL cmessage(6, -1, 1, 15, 4, ibackcol, mess$, "")
krot = 7
END IF
CASE "F": CALL eingabe(igt, 2, ityp(), fa(), LABEL$(), "Factor", ianf, KN)
CASE "X": CALL eingabe(igt, 2, ityp(), xa(), LABEL$(), "Xa", ianf, KN)
CASE "Y": CALL eingabe(igt, 2, ityp(), ya(), LABEL$(), "Ya", ianf, KN)

CASE "ITYP", "ITYPE":
  CALL eingabe(igt, 1, ityp(), fa(), LABEL$(), "ityp", ianf, KN)
CASE "ICOL": CALL eingabe(igt, 1, icol(), fa(), LABEL$(), "icolor", ianf, KN)
CASE "IWI": CALL eingabe(igt, 1, iwi(), fa(), LABEL$(), "iwidth", ianf, KN)
CASE "BCOL"
  mess$ = "Background colour, 15 = white, 0 = black, etc"
  CALL zeingabe(igt, 1, 1, ibackcol, dummy, resp$, mess$)
  mess$ = "Foreground colour, 15 = white, 0 = black, etc"
  CALL zeingabe(igt, 1, 1, ivorcol, dummy, resp$, mess$)
  IF ibackcol = 0 AND ivorcol = 0 THEN ivorcol = 15
  IF ibackcol = 15 AND ivorcol = 15 THEN ivorcol = 0
CASE ELSE: EXIT SUB
END SELECT
wasnun$ = ""
EXIT SUB
END SUB

SUB gtourcontrol (PCONTROL$, igt, LABEL$()) STATIC
'control of GUIDED TOURS
SHARED igtdelay, inewscreen, ibackcol, title$, fdelay
DO WHILE igt > 0
IF istop = 1 GOTO istopagain' used if calculation is continued
readnextgt:
ON EOF(11) GOTO ENDGT
LINE INPUT #11, PCONTROL$
SELECT CASE LEFT$(PCONTROL$, 1)
  CASE " "
    LABEL$(0) = PCONTROL$
    CALL titles(1, title$, PCONTROL$)
    GOTO readnextgt
  CASE "%":  GOTO readnextgt
  CASE "<"
    idelay = VAL(MID$(PCONTROL$, 2, 7))
    CALL gtreturn(idelay * 10, igt)
    IF igt > 0 GOTO readnextgt
  CASE "+":
istopagain: ' used if calculation has been continued
    istop = 0
	IF igtdelay < 0 THEN
	CALL gtreturn(-igtdelay, igt)' delay if igtdelay is negative
	ELSE
	CALL gtreturn(100, igt)
	IF PCONTROL$ = "C" OR PCONTROL$ = "S" OR PCONTROL$ = "N" THEN
	istop = 1: EXIT SUB
	END IF
	END IF
    IF igt > 0 GOTO readnextgt
END SELECT
SELECT CASE PCONTROL$
    CASE "%"
    CASE "-", "_": LABEL$(0) = "": title$ = ""
    CASE "#":
	  LOCATE 25 + ivgashift - 5, 1
	  FOR i = 1 TO 4: LINE INPUT #11, LABEL$(i)
	  ON EOF(11) GOTO ENDGT
	  IF LABEL$(i) = "" THEN EXIT FOR
	  IF LABEL$(i) = "." THEN LABEL$(i) = ""
	  NEXT i
	  IF igtdelay > 0 THEN
	  FOR ii = 1 TO i - 1
	  CALL cmessage(3, ii - i - 1, 1, 4, 15, ibackcol, " " + LABEL$(ii), "")
	  'CALL cmessage(3, ii - i, 1, 4, 15, ibackcol, " " + label$(ii), "")
	  NEXT
	  CALL cmessage(3, ii - i - 1, 1, 4, 15, ibackcol, SPACE$(78), "")
	  CALL gtreturn(100, igt)
	  END IF
	  LABEL$(0) = ""
    CASE "&":
	  INPUT #11, ierease, izeile, idelay, inewscreen
	  ON EOF(11) GOTO ENDGT
	  IF ierease = 1 THEN
	  CALL SCREENCHECK(igraph)
	  CALL init(1, igraph, 0, 15, ivgashift)
	  LINE (1, 1)-(640, 480), 15, BF
	  END IF
	  FOR i = 1 TO 4: LINE INPUT #11, LABEL$(i)
	  ON EOF(11) GOTO ENDGT
	  IF LABEL$(i) = "" THEN EXIT FOR
	  IF LABEL$(i) = "." THEN LABEL$(i) = ""
	  NEXT i
	  FOR ii = 1 TO i - 1
	  CALL cmessage(5, ii - i - 1, 1, 4, 15, ibackcol, " " + LABEL$(ii), "")
	  NEXT ii
	  LABEL$(0) = ""
	  CALL gtreturn(idelay * 10, igt)
    CASE "EOF-GT"
ENDGT:
	 igt = 0: CLOSE (11): LABEL$(0) = ""
	  CLOSE (27)
	  fdelay = .5
	  mess$ = "-------   End of the GUIDED TOUR  ------"
	  CALL cmessage(6, -1, 1, 15, 4, ibackcol, mess$, "")
    CASE ELSE:
	IF igtdelay > 5 AND igtdelay < 80 THEN
	CALL cmessage(1, -1, 1, 15, 1, ibackcol, PCONTROL$ + "      ", "")
	CALL gtreturn(5, igt) 'short delay
	END IF
    EXIT SUB
END SELECT
LOOP
END SUB

SUB gtreturn (inz, igt)
'makes delays or asks for <RETURN> in GUIDED TOURS
SHARED ivgashift, PCONTROL$, ibackcol, readdatafile$, fdelay
SELECT CASE inz
CASE 0
CASE 1 TO 50  'delay in .1 seconds
  t1 = TIMER
  DO UNTIL TIMER - t1 > .1 * inz
  resp$ = INKEY$
  IF resp$ = "q" THEN GOSUB terminateTour
  IF resp$ > "" THEN EXIT DO
  LOOP
CASE 51 TO 1000 'requires RETURN
  PCONTROL$ = MID$(PCONTROL$, 3, 80)
  IF PCONTROL$ > "" THEN
  PCONTROL$ = MID$(PCONTROL$ + SPACE$(77), 1, 80)
  CALL ColorPrint(PCONTROL$, 25 + ivgashift, 1, ibackcol, ibackcol)
  END IF
   CALL ColorPrint("[ok]", 25 + ivgashift, 76, 15, 13)
   resp$ = "":   DO UNTIL resp$ > "": resp$ = INKEY$: LOOP
   CALL ColorPrint(SPACE$(79), 25 + ivgashift, 1, 15, 15)

IF igt > 0 THEN
SELECT CASE resp$
CASE "c": PCONTROL$ = "C"
CASE "s": PCONTROL$ = "S"
CASE "n": PCONTROL$ = "N"
CASE "q"
GOSUB terminateTour
CASE "-", "_"
  mess$ = " GUIDED TOUR suspended, back with command GT"
  CALL cmessage(6, -1, 1, 15, 4, ibackcol, mess$, "")
  igt = -1: PCONTROL$ = ""
CASE ELSE: PCONTROL$ = ""
END SELECT
END IF
END SELECT
EXIT SUB

terminateTour:
  mess$ = "GUIDED TOUR terminated"
  CALL cmessage(6, -1, 1, 15, 4, ibackcol, mess$, "")
  CLOSE (11)
  igt = 0: fdelay = 0.5
RETURN

END SUB

SUB init (inewscreen, igraph, ivorcol, ibackcol, ivgashift) STATIC
SHARED iwidthm, iactive, icol()
ivgashift = 0

SELECT CASE igraph
CASE 9
CALL switchscreen(1, iactive, ibackcol)
  WINDOW (1, 1)-(640, 480)
  IF inewscreen = 1 THEN LINE (1, 1)-(640, 480), ibackcol, BF
  ivgashift = 0
CASE 12
  SCREEN 12
  ivgashift = 5
  WINDOW (1, 1)-(640, 480)
  IF inewscreen = 1 THEN LINE (1, 1)-(640, 480), ibackcol, BF
CASE 2, 3      'Hercules Monochrome
  SCREEN igraph
  FOR i = 1 TO 7: icol(i) = 1: NEXT
  ivorcol = 1
  ibackcol = 0
  WINDOW (1, 1)-(720, 480)
  IF inewscreen = 1 THEN LINE (1, 1)-(720, 480), ibackcol, BF
END SELECT'igraph
inewscreen = 1
EXIT SUB
END SUB

SUB linep (ihp, x1, y1, x2, y2, ipat, iwidth, ico, ityp) STATIC
' allows also to use line width iwidth > 1
' ico= color (0-15)
' ipat: 0 = line, 1 = rectangle, 2 =filled
' ityp = line type (- - -; -- -- -- etc)
' (ihp: switch for HPGL-output, not in use here)
DIM ilinetype(4)
linetype(1) = &HFFFF
linetype(2) = &H4040
linetype(3) = &HF0F0
linetype(4) = &H8888
SHARED iwidthm
SELECT CASE ipat
CASE 0' line
 REM
 SELECT CASE iwidth
 CASE 0, 1
  SELECT CASE ityp'linetype
  CASE 0: LINE (x1, y1)-(x2, y2), ico
  CASE 1: LINE (x1, y1)-(x2, y2), ico, , linetype(ityp)
  CASE 2:
  LINE (x1, y1)-(x2, y2), ico, , &H4040
  CASE 3: LINE (x1, y1)-(x2, y2), ico, , &HF0F0
  CASE 4: LINE (x1, y1)-(x2, y2), ico, , &H8888
  CASE ELSE
    LINE (x1, y1)-(x2, y2), ico
  END SELECT

 CASE ELSE 'iwidth > 1 'Thick lines
 wi = ATN((x2 - x1) / (y2 - y1 + .00000012#))
 pix = .5 * COS(wi): piy = -.5 * SIN(wi)
 pimax = iwidth * pix
 pimay = iwidth * piy
 xa1 = x1 - pimax: xa2 = x2 - pimax: ya1 = y1 - pimay: ya2 = y2 - pimay
 IF iwidth < 4 THEN
   FOR iwiac = 1 TO iwidth
   LINE (xa1, ya1)-(xa2, ya2), ico, , linetype(ityp)
   xa1 = xa1 + pix: xa2 = xa2 + pix: ya1 = ya1 + piy: ya2 = ya2 + piy
   NEXT
 ELSE 'broader lines, filled by PAINT
    xb1 = x1 + pimax: xb2 = x2 + pimax: yb1 = y1 + pimay: yb2 = y2 + pimay
    LINE (xa1, ya1)-(xa2, ya2), ico
    LINE -(xb2, yb2), ico
    LINE -(xb1, yb1), ico
    LINE -(xa1, ya1), ico
    PAINT (x1 + (x2 - x1) / 2, y1 + (y2 - y1) / 2), ico
 END IF
 END SELECT

CASE 1
  SELECT CASE ityp' rectangle
  CASE 1: LINE (x1, y1)-(x2, y2), ico, B
  CASE 2: LINE (x1, y1)-(x2, y2), ico, B, &H4040
  CASE 3: LINE (x1, y1)-(x2, y2), ico, B, &HF0F0
  CASE 4: LINE (x1, y1)-(x2, y2), ico, B, &H8888
  END SELECT
  
CASE 2 'filled rectangle
  LINE (x1, y1)-(x2, y2), ico, BF

CASE 3
  ddy = y2 - y1
  idx = x2 - x1
  idy = ddy
  idx1 = idy - 1
  idx2 = idy - 2
  idx3 = idy - 3
  x21 = x2 - 1
  SELECT CASE ityp
      CASE 10 TO 32000
	LINE (x1, y1)-(x21, y2), ico, BF
      CASE 7 TO 9
	FOR iyd = y1 TO y1 + idx1 STEP 1
	  LINE (x1, iyd)-(x21, iyd), ico, , &H5555
	  iyd = iyd + 1
	  LINE (x1, iyd)-(x21, iyd), ico, , &HAAAA
	NEXT iyd
      CASE 4 TO 6
	FOR iyd = y1 TO y1 + idx2 STEP 2
	  LINE (x1, iyd)-(x21, iyd), ico, , &H5555
	  iyd = iyd + 2
	  LINE (x1, iyd)-(x21, iyd), ico, , &HAAAA
	NEXT iyd
      CASE 1 TO 3
	FOR iyd = y1 TO y1 + idx2 STEP 2
	  LINE (x1, iyd)-(x21, iyd), ico, , &H1111
	  iyd = iyd + 2
	  LINE (x1, iyd)-(x21, iyd), ico, , &H4444
	NEXT iyd
      END SELECT 'shading

END SELECT'ipat
END SUB

SUB lptopen (lpt$)
' checks and if OK opens printer as #8
'            if not, icheckerror=1
SHARED icheckerror
icheckerror = 0
ON ERROR GOTO printerproblem
DEF SEG = 0
POKE &H478, 1'reducing time for printer time out
DEF SEG
OPEN lpt$ FOR OUTPUT AS #8
PRINT #8, CHR$(8)
IF icheckerror = 1 THEN CLOSE (8)
ON ERROR GOTO 0
DEF SEG = 0
POKE &H478, 255'elongating time for printer time out
DEF SEG
END SUB

SUB paramenu STATIC'=======================================================-
'Allows to change floating point parameters (for calculation)
'               or of integer parameters (for program flow)
'  or inputting a command. Control variable:   PCONTROL$
'  e.g. name of the parameter to be changed or command
'--------------------------------------------------------------------------
'   inv(): temporary storage of integer values KT, KP .. KG, K1..
'   flv():        ..    floating point values DA, RA ..
'   parnam$(): contains the names of the variables for identification



SHARED parnam$(), flv(), inv(), igt, iparam
SHARED KP, KT, KX, ky, KD, KI, KE, KR, KN, KG
SHARED K1, K2, K3, K4, DX, DY, DZ, diffmax, imxl, imyl
SHARED DA, ra, ba, sA, ca, aA, gA, DB, rB, bb, sb, cB, aB, gB
SHARED DC, rc, bc, sc, cC, aC, gC, DD, rd, bd, sd, cd, aD, gD
SHARED DE, re, be, se, cE, aE, ge, DF, rf, bf, sf, cF, aF, gf
SHARED DG, rg, bg, sG, cG, aG, gg, DH, rH, bH, sH, cH, aH, gH
SHARED xa(), ya(), fa(), ityp(), icol(), iwi(), inewscreen
SHARED dxp, dyp, per, fkasten, itot&, igtdelay, fdelay, ihp, iactive
SHARED krot, ivorcol, ibackcol, istelz, ihintr, openfilename$, iwohin, ifilenr
SHARED PCONTROL$, pcontrol2$, comment$, icheckerror, ivgashift, programname$
SHARED igraph, ipcontrol, displaywas$, displaywastmp$, filemessage$, basfilename$
SHARED LABEL$(), s$(), ilineparam, exposure, pdir$, title$


IF igraph = 0 THEN 'initialization after program start
  t1 = TIMER
  CALL SCREENCHECK(igraph)

  a$ = "KT,KP,KX,KY,KD,KI,KE,KR,KN,KG,"
  a$ = a$ + "K1,K2,K3,K4,DX,DY,DZ,"
  a$ = a$ + "DA,RA,BA,SA,CA,AA,GA,DB,RB,BB,SB,CB,AB,GB,"
  a$ = a$ + "DC,RC,BC,SC,CC,AC,GC,DD,RD,BD,SD,CD,AD,GD,"
  a$ = a$ + "DE,RE,BE,SE,CE,AE,GE,DF,RF,BF,SF,CF,AF,GF,"
  a$ = a$ + "DG,RG,BG,SG,CG,AG,GG,DH,RH,BH,SH,CH,AH,GH,"
  FOR i = 1 TO 73
    parnam$(i) = MID$(a$, (i - 1) * 3 + 1, 2)
  NEXT i'param$() contains now the names of the variables
  KEY 1, "MENU" + CHR$(13)
  KEY 2, "R" + CHR$(13)
  KEY 3, "RN" + CHR$(13)
  KEY 4, "m4" + CHR$(13)
  KEY 5, "m5" + CHR$(13)
  KEY 6, "PE" + CHR$(13)
  KEY 7, "PP" + CHR$(13)
  KEY 8, "W0" + CHR$(13)
  KEY 9, "R0" + CHR$(13)
  KEY 10, "SN" + CHR$(13)
  KEY 30, "GT" + CHR$(13)
  KEY 31, "RR" + CHR$(13)

  irfilenr$ = programname$ + "1"
  progrtmp$ = programname$
    IF COMMAND$ > "" THEN
     PCONTROL$ = COMMAND$
     GOTO usepcontrol
  END IF
GOTO readfilenow
END IF

parametermenu:
IF iparam > 0 GOTO showparameterlist
CALL prtscreen("sub-util\" + programname$ + ".scr")
showparameterlist: '--------- formatted print of parameters on screen
IF iparam = 2 GOTO parameinput'only bottom line is shown
IF inv(9) < 1 THEN inv(9) = 2
IF inv(9) > 8 THEN
inv(9) = 8: mess$ = "please not more then 8 substances (KN<=8)"
CALL cmessage(4, -1, 1, 15, 4, ibackcol, mess$, "OK")
END IF
GOSUB screenprintparameter
parameinput:
pdirectinput = -1
IF igt > 0 THEN
CALL gtourcontrol(PCONTROL$, igt, LABEL$())  'if GUIDED TOUR is on
IF igt = 0 THEN GOTO parameinput
ELSE
PCONTROL$ = ""
 CALL zeingabe(igt, 0, 3, i, dummy, PCONTROL$, "->")
 pdirectinput = -1  'to find and input of the type "KT = 1"
 t1 = TIMER

END IF
usepcontrol:
PCONTROL$ = UCASE$(PCONTROL$)'------------------ all in uppercase letters
INS = 0: ipcontrol = 0: pcontrol2$ = ""
ipu = INSTR(PCONTROL$, "=")
IF ipu > 0 THEN  'e.g., KT=12 is also a legal input, used for menu
pdirectinput = VAL(MID$(PCONTROL$, ipu + 1, 12))
PCONTROL$ = MID$(PCONTROL$, 1, ipu - 1)
END IF
FOR i = 1 TO 73'----- Find the parameter name in parnam$(INS)
  IF PCONTROL$ = parnam$(i) THEN INS = i
NEXT i
SELECT CASE INS
CASE 1 TO 14'----- input new integer constant
  IF pdirectinput = -1 THEN
  CALL zeingabe(igt, 1, 1, inv(INS), dummy, dummy$, parnam$(INS))
  ELSE
  inv(INS) = pdirectinput
  END IF
  modified$ = ", modified             ": comment$ = ""
  IF INS = 9 THEN iparam = 0: 'number of substances has been changed
GOSUB setparameter
  GOTO showparameterlist
CASE 15 TO 73'-input new floating point constant:
  IF pdirectinput = -1 THEN
  CALL zeingabe(igt, 1, 2, idummy, flv(INS - 14), dummy$, parnam$(INS))
  ELSE
  flv(INS - 14) = pdirectinput
  END IF
  IF LEFT$(PCONTROL$, 1) = "D" AND INS > 17 AND flv(INS - 14) > diffmax THEN
    mess$ = "to avoid numerical instability, diffusion must be < " + STR$(diffmax) + ", corrected"
    CALL cmessage(6, -1, 1, 15, 4, ibackcol, mess$, "")
    flv(INS - 14) = diffmax
  END IF
  modified$ = ", modifiziert            ": comment$ = ""
  GOSUB setparameter
  GOTO showparameterlist
END SELECT

SELECT CASE PCONTROL$
CASE ""'--- input "RETURN" for full menu
  iparam = iparam - 1: GOTO parametermenu

CASE "S", "C", "N", "I", "II", "M", "CS", "D"
  GOTO leave' fixes parameters, jump into the main program
CASE "%", " " 'ignore
CASE "QQ": END
CASE "Q"'That is the END
resp$ = "Y / N": mess$ = " You like to terminate the program?"
    CALL cmessage(4, -1, 1, 15, 4, ibackcol, mess$, resp$)
IF resp$ = "y" OR resp$ = "Y" OR resp$ = "j" THEN CLS : END
CASE "DW"
  CALL zeingabe(igt, 1, 3, i, D, displaywas$, "display what (sab...) ")
  GOTO showparameterlist

CASE "IGTDELAY"
CALL zeingabe(igt, 1, 1, igtdelay, D, dummy$, "Delay in the GUIDED TOUR; in .1 seconds")
CASE "FDELAY"
	 mess$ = "Delay in the Movie, only at SCREEN 9, 10 = 1sec"
	 CALL zeingabe(igt, 1, 2, i, fdelay, resp$, mess$)
CASE "M4", "M5"
SELECT CASE programname$
CASE "SP", "XT"
    IF PCONTROL$ = "M4" THEN pcontrol2$ = "4"
    IF PCONTROL$ = "M5" THEN pcontrol2$ = "5"
  PCONTROL$ = "M"
  GOTO leave

CASE ELSE
    displaywastmp$ = ""
    IF PCONTROL$ = "M4" THEN ipcontrol = 12
    IF PCONTROL$ = "M5" THEN ipcontrol = 13
    PCONTROL$ = "D"
    GOTO leave
END SELECT


CASE "EGA": igraph = 9: CALL SCREENCHECK(igraph)
CASE "VGA": igraph = 12: CALL SCREENCHECK(igraph)
CASE "PR" 'Druckt Bemerkung (LPT1:)
CALL lptopen("LPT1:"): IF icheckerror = 1 GOTO parametermenu
GOSUB remark

CASE "PE"
  such$ = "SELECT CASE KE'EQUATION"
  lpt$ = "SCRN:": OPEN lpt$ FOR OUTPUT AS #8
  mess$ = "press any key, or `p` for printing"
  GOSUB FindEquation
  GOSUB screenprintparameter
  iparam = 1
  resp$ = ""
  CALL cmessage(4, -1, 1, 15, 2, ibackcol, mess$, resp$)
  IF ipu = 0 OR resp$ <> "p" GOTO parametermenu
  lpt$ = "LPT1:": CALL lptopen(lpt$):
  IF icheckerror = 1 THEN CLOSE : GOTO parametermenu
  GOSUB FindEquation
  CALL cmessage(6, -1, 1, 15, 2, 0, "eject page with command FF, if you like", "")
  GOTO parametermenu
CASE "NEC": CALL cmessage(6, -1, 1, 15, 2, 0, "The NEC printer is no longer supported", "")
CASE "GT"
  IF igt < 0 THEN
  igt = 1
  ELSE
gtlist:
  CALL whatprmfiles(programname$, programname$, "GT", pdir$, irfilenr$)
  CALL SCREENCHECK(igraph)
  iparam = 0
  IF irfilenr$ = "q" OR irfilenr$ = "p" OR irfilenr$ = "" GOTO parametermenu
  gtfilename$ = pdir$ + irfilenr$ + ".GT"
readgtnow:
  icheckerror = 0
  igraph = 12: CALL SCREENCHECK(igraph)
  LINE (1, 1)-(640, 480), 15, BF

  ON ERROR GOTO checkerror
  OPEN gtfilename$ FOR INPUT AS #11
  IF icheckerror = 1 THEN
    'ivgashift = 0
    mess$ = " the file " + gtfilename$ + " does not exist. You will get a list"
    resp$ = "ok": CALL cmessage(2, -1, 1, 15, 4, ibackcol, mess$, resp$)
    GOTO gtlist
  END IF
  LINE INPUT #11, resp$
  igt = 1
END IF
iparam = 2

CASE "PP" 'prints parameter
CALL lptopen("LPT1:"): IF icheckerror = 1 GOTO parametermenu
PCONTROL$ = ""
GOSUB printparameter
CALL cmessage(6, -1, 1, 15, 2, 0, "eject page with command FF, if you like", "")

CASE "FF":
CALL lptopen("LPT1:"): IF icheckerror = 1 GOTO parametermenu
PRINT #8, CHR$(12): CLOSE (8)

CASE "PDIR"
  CALL zeingabe(igt, 1, 3, i, D, pdir$, "Directory of parameter files ")
  IF RIGHT$(pdir$, 1) <> "\" THEN pdir$ = pdir$ + "\"

CASE "CLS":  iparam = 2:
igraph = 12: CALL SCREENCHECK(igraph)
  LINE (1, 1)-(640, 480), ibackcol, BF


CASE "MENU"
  PCONTROL$ = ""
  CALL whatprmfiles(programname$, "HELP", "HLM", "sub-util\", irfilenr$)
  IF irfilenr$ <> "q" THEN
  OPEN "sub-util\" + irfilenr$ + ".HLM" FOR INPUT AS #2
  LINE INPUT #2, mess$: i = 1
  DO UNTIL EOF(2)
    LINE INPUT #2, s$(i)
    IF s$(i) = "" THEN EXIT DO
    i = i + 1
  LOOP
  CLOSE #2
  s$(i) = ""
  CALL smenu(PCONTROL$, iline, mess$)
  IF PCONTROL$ = "q" THEN PCONTROL$ = ""
  IF RIGHT$(PCONTROL$, 1) = "." THEN
    COLOR 15, 4
    LOCATE 19, 1:
    PRINT SPACE$(10) + "ÚÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄ¿"
    PRINT SPACE$(10) + "³  substance to which the change should apply, a,b,...     ³"
    PRINT SPACE$(10) + "ÀÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÄÙ"
    resp$ = "": DO UNTIL resp$ > "": resp$ = INKEY$: LOOP
    PCONTROL$ = LEFT$(PCONTROL$, 1) + resp$
  END IF
  END IF
  CALL SCREENCHECK(igraph)
  iparam = 0
  LOCATE 20 - inv(9) + ivgashift, 1
  OPEN "SCRN:" FOR OUTPUT AS #8
  GOSUB printparameter
  GOTO usepcontrol


CASE ELSE
pcleft$ = LEFT$(PCONTROL$, 1)'number behind display can be used
pcontrol2$ = MID$(PCONTROL$, 2, 4)
ipcontrol = VAL(pcontrol2$)
SELECT CASE pcleft$
CASE "M"
  PCONTROL$ = pcleft$
  GOTO leave
CASE "W", "R", "S"
  IF MID$(PCONTROL$, 1, 2) = "RR" GOTO readfilenow1
  IF MID$(PCONTROL$, 1, 2) = "WW" GOTO writefilenow
  IF MID$(PCONTROL$, 2, 1) = "N" THEN
  ilineparam = ilineparam + 1
  ipu = INSTR(s$(ilineparam), " ")
  IF ipu > 0 THEN
  filenametmp$ = pdir$ + MID$(s$(ilineparam), 1, ipu) + ".prm"
  GOTO readfilenow1
  ELSE
  mess$ = "either last simulation of this chapter, or use F2 first"
  CALL cmessage(6, -1, 1, 15, 4, ibackcol, mess$, "")
  GOTO parameinput
  END IF
  END IF
  irfilenr$ = MID$(PCONTROL$, 2, 6)
  ipu = INSTR(irfilenr$, ".")
  IF ipu > 0 THEN irfilenr$ = MID$(irfilenr$, 1, ipu - 1) + MID$(irfilenr$, ipu + 1, 8)
  IF irfilenr$ > "" THEN
  irfilenr$ = progrtmp$ + irfilenr$: GOTO readfilenow
  ELSE
  IF pcleft$ = "W" THEN
	mess$ = "please use the command Wxxx, xxx=filename or number"
	CALL cmessage(6, -1, 1, 15, 4, ibackcol, mess$, "")
  GOTO parametermenu
  END IF
  END IF
  mess$ = "parameter for which chapter (2-9 or 12 <RETURN> = all), or ESC"
  CALL zeingabe(igt, 0, 4, i, dummy, resp$, mess$)
  IF resp$ = "q" OR resp$ = CHR$(27) THEN GOTO parametermenu
  CLS
  progrtmp$ = programname$ + resp$
whatprmfileslabel:
CALL whatprmfiles(programname$, progrtmp$, "prm", pdir$, irfilenr$)
  CALL SCREENCHECK(igraph)
  IF irfilenr$ = "q" THEN iparam = 0: progrtmp$ = programname$: GOTO parametermenu
readfilenow:
  filenametmp$ = pdir$ + irfilenr$ + ".prm"
  progrtmp$ = programname$
  IF pcleft$ = "W" GOTO writefilenow
readfilenow1:
  icheckerror = 0
  ON ERROR GOTO checkerror
  OPEN filenametmp$ FOR INPUT AS #9
  IF icheckerror = 1 THEN
    mess$ = " the file " + filenametmp$ + " does not exist. Use F2 for a list"
    CALL cmessage(6, -1, 1, 15, 4, ibackcol, mess$, "")
GOTO parametermenu
  END IF
  filename$ = filenametmp$
  CALL rwparameterdisk("R", 9)
  CLOSE (9)
  comment$ = "read from file "
  modified$ = ""
  LABEL$(0) = filemessage$
  title$ = filemessage$
  ON ERROR GOTO 0
  IF igt = 0 THEN iparam = 0
  IF pcleft$ = "S" THEN
  PCONTROL$ = LEFT$(PCONTROL$, 1)
  GOSUB leave'direct start after reading a parameter file
  END IF
  GOSUB setparameter
  GOTO parametermenu

writefilenow:
  mess$ = "New title of this set of parameter:"
IF LEN(filemessage$) > 40 THEN
  mess$ = mess$ + " [" + MID$(filemessage$, 1, 40) + "...." + "]"
  ELSE
  mess$ = mess$ + " [" + filemessage$ + "]"
  END IF
  CALL cmessage(3, -2, 1, 15, 2, ibackcol, mess$, "")
  CALL zeingabe(igt, 0, 4, i, dummy, response$, "")
  CALL cmessage(3, -2, 1, 0, 0, ibackcol, SPACE$(79), "")
  IF response$ = CHR$(27) OR response$ = "q" GOTO parametermenu
  IF response$ > "" THEN filemessage$ = response$
  filename$ = filenametmp$
  OPEN filename$ FOR OUTPUT AS #9
  CALL rwparameterdisk("W", 9)
  ON ERROR GOTO 0
  CLOSE (9)
  comment$ = "written to file "
  modified$ = ""
  iparam = 0
  GOTO parametermenu
CASE "G"
gtfilename$ = MID$(PCONTROL$, 3, 5)
IF gtfilename$ > "" THEN
  gtfilename$ = pdir$ + programname$ + gtfilename$ + ".GT"
  GOTO readgtnow
  END IF
END SELECT
checkpcontrol:
CALL gmanage(PCONTROL$)
IF PCONTROL$ = "" GOTO parametermenu
mess$ = PCONTROL$ + " is a nonsense input, try again or use F1"
CALL cmessage(6, -1, 1, 15, 4, ibackcol, mess$, "")
END SELECT
GOTO parametermenu

setparameter:
' transferring the array contents into single parameters used by the equations
KT = inv(1): KP = inv(2): KX = inv(3): ky = inv(4): KD = inv(5)
KI = inv(6): KE = inv(7): KR = inv(8): KN = inv(9): KG = inv(10)
K1 = inv(11): K2 = inv(12): K3 = inv(13): K4 = inv(14)
DX = flv(1): DY = flv(2): DZ = flv(3)
DA = flv(4): ra = flv(5): ba = flv(6): sA = flv(7): ca = flv(8): aA = flv(9): gA = flv(10)
DB = flv(11): rB = flv(12): bb = flv(13): sb = flv(14): cB = flv(15): aB = flv(16): gB = flv(17)
DC = flv(18): rc = flv(19): bc = flv(20): sc = flv(21): cC = flv(22): aC = flv(23): gC = flv(24)
DD = flv(25): rd = flv(26): bd = flv(27): sd = flv(28): cd = flv(29): aD = flv(30): gD = flv(31)
DE = flv(32): re = flv(33): be = flv(34): se = flv(35): cE = flv(36): aE = flv(37): ge = flv(38)
DF = flv(39): rf = flv(40): bf = flv(41): sf = flv(42): cF = flv(43): aF = flv(44): gf = flv(45)
DG = flv(46): rg = flv(47): bg = flv(48): sG = flv(49): cG = flv(50): aG = flv(51): gg = flv(52)
DH = flv(53): rH = flv(54): bH = flv(55): sH = flv(56): cH = flv(57): aH = flv(58): gH = flv(59)
  IF inv(4) > imyl - 1 THEN
    inv(4) = imyl - 1: ky = inv(4)
    mess$ = "KY (Number of cells) too large (<=" + STR$(imyl - 1) + "), will be corrected"
    resp$ = "ok": CALL cmessage(4, -1, 1, 15, 4, ibackcol, mess$, resp$)
  END IF
  IF inv(3) > imxl - 1 THEN
    inv(3) = imxl - 1: KX = inv(3)
    mess$ = "Kx (Number of cells) too large (<=" + STR$(imxl - 1) + "), will be corrected"
    resp$ = "ok": CALL cmessage(4, -1, 1, 15, 4, ibackcol, mess$, resp$)
  END IF
RETURN

leave: '-------------to be executed before  going back to the main program----
GOSUB setparameter
IF KP = 0 OR KT = 0 OR KX = 0 OR ky = 0 OR KE = 0 OR DX = 0 THEN
 mess$ = "list of parameters incomplete! check KT, KP, KX, KY, KE or DX"
 resp$ = "ok": CALL cmessage(4, -1, 1, 15, 4, ibackcol, mess$, resp$)
 GOTO parameinput
END IF
EXIT SUB

screenprintparameter:
FOR i = 19 - inv(9) + ivgashift TO 25 + ivgashift
LOCATE i, 1: PRINT SPACE$(80); : NEXT i
LOCATE 20 - inv(9) + ivgashift, 1
OPEN "SCRN:" FOR OUTPUT AS #8

printparameter:
PRINT #8, "Progr.: "; basfilename$; "; "; comment$; filename$; modified$; ""
ishift = (74 - (LEN(filemessage$))) / 2: IF ishift < 1 THEN ishift = 1
PRINT #8, SPACE$(ishift) + filemessage$
FOR i = 1 TO 2: PRINT #8, USING "####&& "; inv(i); "-"; parnam$(i); : NEXT i
FOR i = 3 TO 10: PRINT #8, USING "###&& "; inv(i); "-"; parnam$(i); : NEXT i
PRINT #8,
FOR i = 11 TO 14
  PRINT #8, USING "###&& "; inv(i); "-"; parnam$(i);
NEXT i: PRINT #8, "     ";
FOR k = 1 TO 3
  PRINT #8, USING "###.####&&"; flv(k); "-"; parnam$(k + 14);
NEXT k
PRINT #8, "  " + LEFT$(displaywas$ + "-DW        ", 10)
FOR k = 1 TO (7 * inv(9))
  PRINT #8, USING "###.####&&"; flv(k + 3); "-"; parnam$(k + 17);
  IF k MOD 7 = 0 THEN PRINT #8,
NEXT k
CLOSE (8)
RETURN

remark:
mess$ = "Remark to printer, stop with <RETURN> or <ESC>"
CALL cmessage(3, -2, 1, 15, 2, ibackcol, mess$, "")
  CALL zeingabe(igt, 0, 4, i, dummy, resp$, "Text: ")
IF resp$ = "" OR resp$ = CHR$(27) THEN
CALL cmessage(6, -1, 1, 15, 2, 0, "eject page with command FF, if you like", "")
'CALL cmessage(3, -2, 1, ibackcol, ibackcol, ibackcol, SPACE$(80), "")
CLOSE #8: GOTO parametermenu
ELSE
  PRINT #8, resp$: GOTO remark
END IF
RETURN


FindEquation:
  icheckerror = 0
  ON ERROR GOTO checkerror
OPEN basfilename$ FOR INPUT AS #2
IF icheckerror = 1 GOTO parametermenu
DO UNTIL EOF(2)
  LINE INPUT #2, s$
  s$ = UCASE$(s$)
  ipu = 1: ipu = INSTR(ipu, s$, such$)
  IF ipu > 0 THEN
  EXIT DO
  END IF
  LOOP
DO UNTIL EOF(2)
  LINE INPUT #2, s$
  s$ = UCASE$(s$)
  ipu = 1: ipu = INSTR(ipu, s$, "CASE" + STR$(inv(7)) + " ")
  IF ipu > 0 THEN EXIT DO
LOOP
IF ipu = 0 THEN
mess$ = "sorry, no equation of this type found, press any key [ok]"
ELSE
IF lpt$ = "SCRN:" THEN
  CLS : LOCATE 5, 1
  END IF
	PRINT #8, s$
	DO UNTIL EOF(2)
	  LINE INPUT #2, s$
	  ipu = 1: ipu = INSTR(ipu, s$, "CASE")
	  IF ipu = 0 THEN PRINT #8, s$
	  IF ipu > 0 THEN EXIT DO
	LOOP
END IF
CLOSE (2)
CLOSE (8)
RETURN

END SUB'===========end of sub paramenu()  ========================

SUB prtscreen (scrfile$)
SHARED ivgashift, ibackcol
OPEN scrfile$ FOR INPUT AS #2
izeile = 1
DO
LINE INPUT #2, mess$
IF EOF(2) THEN EXIT DO
CALL cmessage(3, izeile, 1, 15, 1, ibackcol, mess$, resp$)
IF izeile = 23 THEN EXIT DO
izeile = izeile + 1
LOOP
CLOSE (2)
END SUB

SUB rwparameterdisk (rw$, ifilenumber)
'reads (rw$="R") or writes parameters from or to disk
SHARED flv(), inv(), KN
SHARED xa(), ya(), fa(), ityp(), icol(), iwi()
SHARED dxp, dyp, per, fkasten
SHARED krot, ivorcol, ibackcol, istelz, ihintr
SHARED PCONTROL$, pcontrol2$, comment$, icheckerror, ivgashift, programname$
SHARED igraph, ipcontrol, displaywas$, filemessage$

SELECT CASE rw$
CASE "R"
  IF icheckerror = 2 THEN EXIT SUB
  LINE INPUT #ifilenumber, filemessage$
  FOR i = 1 TO 14
    'ON ERROR GOTO 0
    INPUT #ifilenumber, inv(i)
    IF i MOD 10 = 0 THEN INPUT #ifilenumber, dummy
  NEXT i
  FOR k = 1 TO 3: INPUT #ifilenumber, flv(k): NEXT k
  INPUT #ifilenumber, dummy
  KN = inv(9):
  FOR k = 1 TO 7 * KN
    INPUT #ifilenumber, flv(k + 3)
    IF k MOD 7 = 0 THEN INPUT #ifilenumber, dummy
  NEXT k

  FOR i = 0 TO inv(9)
  INPUT #ifilenumber, xa(i), ya(i), fa(i), ityp(i), icol(i), iwi(i): NEXT
  INPUT #ifilenumber, displaywas$, DX, dxp, dyp, per, fkasten
  INPUT #ifilenumber, krot, ivorcol, ibackcol, istelz, ihintr
  IF EOF(ifilenumber) THEN
	igraph = 12
	ELSE
	INPUT #ifilenumber, igraph
  END IF
  CALL SCREENCHECK(igraph)
CASE "W"
  PRINT #ifilenumber, filemessage$
  FOR i = 1 TO 10: PRINT #ifilenumber, USING "####"; inv(i); : PRINT #ifilenumber, ","; : NEXT i
  PRINT #ifilenumber,
  FOR i = 11 TO 14: PRINT #ifilenumber, USING "####"; inv(i); : PRINT #ifilenumber, ","; : NEXT i
  FOR k = 1 TO 3: PRINT #ifilenumber, USING "###.####"; flv(k); : PRINT #ifilenumber, ","; : NEXT k
  PRINT #ifilenumber,
  FOR k = 1 TO 7 * inv(9)
    PRINT #ifilenumber, USING "###.####"; flv(k + 3); : PRINT #ifilenumber, ",";
    IF k MOD 7 = 0 THEN PRINT #ifilenumber,
  NEXT k
  FOR i = 0 TO inv(9)
  WRITE #ifilenumber, xa(i), ya(i), fa(i), ityp(i), icol(i), iwi(i): NEXT
  WRITE #ifilenumber, displaywas$, DX, dxp, dyp, per, fkasten
  WRITE #ifilenumber, krot, ivorcol, ibackcol, istelz, ihintr
  WRITE #ifilenumber, igraph
END SELECT
modified$ = ""
END SUB

SUB SCREENCHECK (igraph)
' checks which graphic card is installed; sets windows
SHARED tryanother, ivgashift
ON ERROR GOTO wrongscreen
IF igraph = 0 THEN igraph = 12
tryanother = 0
FOR i = igraph TO 1 STEP -1
  SCREEN i
  IF tryanother = 0 THEN EXIT FOR
  tryanother = 0
NEXT i
IF i = 1 THEN
  PRINT " no graphic cardfound. You can use"
  PRINT " the display type KD=4 only. Somethink seems to be wrong"
ELSE
  igraph = i
END IF
ivgashift = 0
  WINDOW (1, 1)-(640, 480)
  IF igraph = 12 THEN ivgashift = 5
ON ERROR GOTO 0
END SUB

SUB smenu (irfilename$, iline, headline$)
'Produces a menu for selection, s$() contains the items to choose
'iline is returned
SHARED s$(), ivgashift, icheckerror
SCREEN 0
ivgashift = 0
CLS
'length of menu lines
izlenmax = LEN(headline$) + 4
FOR ial = 1 TO 150
izlen = LEN(s$(ial))
IF izlen = 0 THEN EXIT FOR
IF izlen > 75 THEN izlen = 75: s$(ial) = LEFT$(s$(ial), 72) + "..."
IF izlen > izlenmax THEN izlenmax = izlen
NEXT ial
iarmax = ial - 1: iline = 1: ivon = 1:
smenunew:
COLOR 15, 1
LOCATE 1, 1
PRINT "ÚÄ  " + headline$ + " " + STRING$(izlenmax - LEN(headline$) - 2, 196) + "¿"
ibis = ivon + 21: IF ibis > iarmax THEN ibis = iarmax
FOR ial = ivon TO ibis
PRINT "³ " + s$(ial) + SPACE$(izlenmax - LEN(s$(ial))) + " ³"
NEXT ial
IF ibis = iarmax THEN
PRINT "ÀÄ" + STRING$(izlenmax, 196) + "ÄÙ";
ELSE
PRINT "ÀÄ" + STRING$(izlenmax - 9, 196) + " more ...ÄÙ";
END IF
COLOR 15, 2
DO
LOCATE 25, 3
LOCATE 25, 3: PRINT "Use arrows, Pg-UP and PG-DOWN to select, RETURN to take, ESC to quit, p=print";
'selectline
COLOR 1, 15'2, 0
LOCATE 2 + iline - ivon, 3
ioldline = iline
PRINT s$(iline) + SPACE$(izlenmax - LEN(s$(iline)))
i$ = ""
WHILE i$ = ""
i$ = INKEY$
WEND
i = ASC(RIGHT$(i$, 1))
SELECT CASE i
CASE 27, 113: irfilename$ = "q": EXIT SUB
CASE 80, 32: IF iline < iarmax THEN iline = iline + 1
IF iline > ibis THEN
ivon = ivon + 22
IF ivon > iarmax - 21 THEN ivon = iarmax - 21
IF ivon < 1 THEN ivon = 1
iline = ivon
END IF
GOTO smenunew
CASE 72: IF iline > 1 THEN iline = iline - 1
IF iline < ivon THEN ivon = ivon - 22
IF ivon < 1 THEN ivon = 1
GOTO smenunew
CASE 81'pg-down
ivon = ivon + 22
iline = ivon
IF ivon > iarmax - 21 THEN
ivon = iarmax - 21
IF ivon < 1 THEN ivon = 1
iline = iarmax
END IF
GOTO smenunew
CASE 73'pg-up
ivon = ivon - 22
IF ivon < 1 THEN ivon = 1
iline = ivon
GOTO smenunew
CASE 13
irfilename$ = s$(iline)
ipu = INSTR(s$(iline), " ")
IF ipu > 0 THEN irfilename$ = LEFT$(s$(iline), ipu - 1)
EXIT SUB
CASE 112
   CALL lptopen("LPT1:"):
   IF icheckerror = 1 THEN
   irfilename$ = "q":
   EXIT SUB
   END IF
PRINT #8, "ÚÄ  " + headline$ + " " + STRING$(izlenmax - LEN(headline$) - 2, 196) + "¿"
FOR ial = 1 TO iarmax
PRINT #8, "³ " + s$(ial) + SPACE$(izlenmax - LEN(s$(ial))) + " ³"
NEXT ial
PRINT #8, "ÀÄ" + STRING$(izlenmax, 196) + "ÄÙ"
   CLOSE #8
GOTO smenunew
CASE ELSE
BEEP
END SELECT
LOCATE 2 + ioldline - ivon, 3
PRINT s$(ioldline) + SPACE$(izlenmax - LEN(s$(ioldline)))
LOOP
END SUB

SUB sorts (s$(), icount)
' to sort e.g. files for the menu, icount = number of elements
ICA = 0
FOR icc = 0 TO icount
VAM$ = "zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz"
FOR ICD = ICA TO icount
IF s$(ICD) > VAM$ GOTO SNEXTICD
IVAM = ICD
VAM$ = s$(ICD)
SNEXTICD:
NEXT ICD
s$(IVAM) = s$(ICA)
s$(ICA) = VAM$
ICA = ICA + 1
NEXT icc
EXIT SUB
 END SUB

SUB switchscreen (inz, iactive, ibackcol) STATIC
SHARED KD, programname$, igraph, fdelay
'switch active and displayed screens 9
'inz=1 :initiation; inz=2: toggle
'inz=3 termination, active screen is displayed
IF igraph <> 9 THEN EXIT SUB
IF programname$ = "XT" OR programname$ = "SP" THEN
SELECT CASE KD
CASE 3, 4, 5, 16, 18, 19, 21, 22, 23'only applicable for that modes
CASE ELSE: EXIT SUB
END SELECT
END IF

SELECT CASE inz
CASE 1'initiation
	SELECT CASE iactive
	CASE 0: SCREEN 9, , 1, 0: iactive = 1
	CASE ELSE: SCREEN 9, , 0, 1: iactive = 0
	END SELECT
CASE 2'switch
	SELECT CASE iactive
	CASE 0: SCREEN 9, , 1, 0: iactive = 1
	CASE 1: SCREEN 9, , 0, 1: iactive = 0
	END SELECT
	IF fdelay > 0 THEN DO UNTIL TIMER - t1 > fdelay * .1: LOOP
	t1 = TIMER
CASE 3 'end
    SELECT CASE iactive
	CASE 0: SCREEN 9, , 1, 1: iactive = 1
	CASE 1: SCREEN 9, , 0, 0: iactive = 0
      END SELECT
CASE 4'initiation, visible
	SELECT CASE iactive
	CASE 0: SCREEN 9, , 0, 0: iactive = 0
	CASE ELSE: SCREEN 9, , 1, 1: iactive = 1
	END SELECT
END SELECT
END SUB

SUB titles (inz, title$, mess$)
  ishift = (74 - (LEN(mess$))) / 2: IF ishift < 1 THEN ishift = 1
  title$ = SPACE$(ishift) + mess$
END SUB

SUB whatprmfiles (programname$, progrtmp$, exten$, drive$, irfilename$)
' scans the disk for a particular type of files
' assures the correct order in the menu
' s$() contains afterwards the file names and the first line of that files
' i.e. the head lines
SHARED s$(), ilineparam, ivgashift, icheckerror
IF progrtmp$ = "" THEN
  mess$ = "Program name to be scanned (ESC to skip)  [" + programname$ + "]"
  CALL zeingabe(igt, 0, 4, i, dummy, progrtmp$, mess$)
  IF progrtmp$ = "q" OR progrtmp$ = CHR$(27) THEN irfilename$ = "q": EXIT SUB
  IF progrtmp$ = "" THEN progrtmp$ = programname$
END IF
exten$ = UCASE$(exten$)
ilr = LEN(progrtmp$)
dymmy$ = drive$ + progrtmp$ + "*." + exten$
s1$ = DIR$(drive$ + progrtmp$ + "*." + exten$)
IF s1$ = "" THEN
resp$ = "ok"
CALL cmessage(4, -1, 1, 15, 4, 0, "sorry, no such files found", resp$)
irfilename$ = "q": EXIT SUB
END IF
CALL cmessage(3, -1, 1, 15, 2, 0, "... reading headlines....", "")
icount = 0: ilinecount = 0
DO
IF s1$ > "" THEN
OPEN drive$ + s1$ FOR INPUT AS #2
LINE INPUT #2, si$
CLOSE #2
icount = icount + 1
ipu = 1: ipu = INSTR(ipu, s1$, "." + exten$)
s1$ = MID$(s1$, LEN(progrtmp$) + 1, ipu - LEN(progrtmp$) - 1)
s1$ = s1$ + "        "
s$(icount) = MID$(s1$, 1, 7) + si$
IF VAL(s$(icount)) < 10 THEN
s$(icount) = ".." + s$(icount)
 ELSEIF VAL(s$(icount)) < 100 THEN
s$(icount) = LEFT$(s$(icount), 1) + "." + MID$(s$(icount), 2, 250)

 END IF
END IF
s1$ = DIR$
LOOP UNTIL s1$ = ""
s$(icount + 1) = ""
CALL sorts(s$(), icount)
mess$ = "FILES: " + drive$ + progrtmp$ + "*." + exten$ + "     TITLES"
IF icount > 0 THEN
FOR i = 1 TO icount
IF LEFT$(s$(i), 2) = ".." THEN s$(i) = MID$(s$(i), 3, 250)
IF MID$(s$(i), 2, 1) = "." THEN s$(i) = LEFT$(s$(i), 1) + MID$(s$(i), 3, 250)
s$(i) = progrtmp$ + s$(i)
NEXT i
smenuagain:
CALL smenu(irfilename$, ilineparam, mess$)
EXIT SUB
  ELSE 'icount=0
   PRINT "sorry, no such files"
   irfilename$ = "q"
 END IF
END SUB

SUB zeingabe (igt, iquest, inz, i, f, LABEL$, text$)
'allows inputting single constants, if igt>O from file (guided tour)
'iquest =1: original value is prompted
'inz=1 integer i: inz=2 floating point f, inz=3: string label$
'text$ is the message displayed with the request
'---------------------------------------------------------
SHARED ivgashift, igtdelay
textprint$ = text$ + " "
IF iquest = 1 THEN
  textprint$ = textprint$ + " OLD = "
  IF inz = 1 THEN textprint$ = textprint$ + STR$(i)
  IF inz = 2 THEN textprint$ = textprint$ + STR$(f)
  IF inz = 3 OR inz = 4 THEN textprint$ = textprint$ + LABEL$
textprint$ = textprint$ + "  NEW = "
END IF
CALL cmessage(1, -1, 1, 2, 15, ibackcol, textprint$, resp$)
ilength = LEN(textprint$) + 1
fac$ = ""
CALL cmessage(1, -1, ilength, 12, 15, ibackcol, fac$ + "_ ", resp$)
SELECT CASE igt
CASE -1000 TO 0
fac$ = ""
DO
DO: dummy$ = INKEY$: LOOP UNTIL dummy$ > ""
IF dummy$ = CHR$(13) THEN EXIT DO
IF dummy$ = CHR$(8) THEN
  lentext = LEN(fac$)
  IF LEN(fac$) > 0 THEN fac$ = LEFT$(fac$, LEN(fac$) - 1)
CALL cmessage(1, -1, ilength, 12, 15, ibackcol, fac$ + "_  ", resp$)
ELSE
SELECT CASE inz
CASE 1, 2
SELECT CASE dummy$
CASE "-", "+", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "."
fac$ = fac$ + dummy$
CALL cmessage(1, -1, ilength, 12, 15, ibackcol, fac$ + "_  ", resp$)
CASE CHR$(27): GOTO zende
CASE ELSE
BEEP
END SELECT
CASE 3, 4
IF dummy$ = CHR$(27) AND inz = 3 GOTO zende
IF dummy$ = CHR$(27) AND inz = 4 THEN LABEL$ = CHR$(27): GOTO zende
fac$ = fac$ + dummy$
CALL cmessage(1, -1, ilength, 12, 15, ibackcol, fac$ + "_  ", resp$)
END SELECT
END IF
LOOP

CASE 1 TO 1000: LINE INPUT #11, fac$ 'only for GUIDED TOURS
IF igtdelay > 10 THEN CALL gtreturn(8, igt)
CALL cmessage(1, -1, ilength, 12, 15, ibackcol, fac$ + "_  ", resp$)
IF igtdelay > 1 THEN CALL gtreturn(8, igt)
IF fac$ = "x" GOTO zende
END SELECT
IF inz = 1 THEN i = VAL(fac$)
IF inz = 2 THEN f = VAL(fac$)
IF inz = 3 OR inz = 4 THEN LABEL$ = fac$
zende:
CALL cmessage(1, -1, 1, 15, 15, ibackcol, SPACE$(79), resp$)
END SUB


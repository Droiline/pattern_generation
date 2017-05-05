'Program for the simulation of pattern formation in two-dimensional fields
' (c) by Hans Meinhardt, 1995-2002

'-----------------------------------------------------------------------
' A word of precaution...                                                                 ³
'  The software was written with care but the customer knows that
'  software cannot be produced without errors. Springer-Verlag and the
'  author disclaim any liability for software errors or their
'  consequences.
'
'  All user rights to the software are held exclusively by the author.
'  Springer-Verlag does not author data and programs but makes them
'  available. Any commercial use or distribution independently of the
'  book is not permitted
'-----------------------------------------------------------------------

' This version is appropriate for the compilation with POWER-BASIC;
' use XYMS.BAS for MS - BASIC

$ERROR ALL ON '    de-activate for MS-BASIC!
DEFINT H-N
DECLARE SUB cmessage (itype%, iline%, irow%, ilettercol%, iletterback%, ibackcol%, text$, resp$)
DECLARE SUB message (itype%, iline%, text$)
DECLARE SUB init (ihp%, igraph%, ivorcol%, ibackcol%, ivgashift%)
DECLARE SUB switchscreen (inz%, iactive%, ibackcol%)
DECLARE SUB gmanage (wasnun$)
DECLARE SUB paramenu ()
DECLARE SUB manipulation (mwasnun$, iinsert)
DECLARE SUB manipulate ()
DECLARE SUB showdisplay (ax!(), iwas%, kdisplay%, displaywas$)
DECLARE SUB inip (igraph%, ivorcol%, ibackcol%, ivgashift%)
'$INCLUDE: 'scommon.bi'
imxl = 36: ilm = 9: imx = 36: imy = 36
'-----------------------------for Power-Basic----------------
pbilm = 9: pbimx = 82: pbimy = 82
$INCLUDE "sub-util\plibxy.inc" '   de-activate this for MS-BASIC!
DIM ax(ilm, imx, imy) 'for very large arrays, this statement must be de-activated...
''$INCLUDE "dimhuge.inc" ' ...and this must activated, becomes slower
imxl = imx: imyl = imy
programname$ = "XY"
pdir$ = "PXY\"
basfilename$ = "XY.BAS"
DIM parnam$(75), flv(60), inv(15)' to store parameters
DIM label$(10)
DIM asu(ilm), u(ilm, imx), al(imy), difa(ilm), afc(ilm)
DIM kax(9), KX(9), kay(9), ky(9), xa(9), fa(9), am(9)
DIM ya(9), f1(9), icol(9), ityp(9), iwi(9)
DIM tv$(9), afa(9), itot&(9), s$(150)
DIM a(1, 1), B(1, 1), C(1, 1), D(1, 1), sa(1, 1)
DIM mempos(imx), icolm(9), ibm(16), a2(9), ansp(9)
msiluetmax = 500
DIM msiluet(msiluetmax)
ihp = 1: fdelay =0.5
diffmax = .2 'Maximum numerical constant for diffusion
CLS : RANDOMIZE TIMER

'--- set parameter and check screen
parameter:
CALL paramenu

Mainstart:
mess$ = ""
SELECT CASE PCONTROL$
CASE "D" 'To display a single figure
iinit = 1:
   CALL init(ihp, igraph, ivorcol, ibackcol, ivgashift)
IF ipcontrol = 0 THEN
  CALL showdisplay(ax(), iinit, KD, displaywas$)
  ELSE
  CALL showdisplay(ax(), iinit, ipcontrol, displaywastmp$)
  END IF
IF igraph = 9 THEN CALL switchscreen(2, iactive, ibackcol)
GOTO EndSimulation

CASE "M"
       iinsert = 0
	CALL manipulation(pcontrol2$, iinsert)
	GOTO makedraw

CASE "S", "N", "C", "I", "II"
ON ERROR GOTO erroverflow
      icount = 0
      igrowth = 1
      IF PCONTROL$ = "N" GOTO screeninit
      IF PCONTROL$ = "C" GOTO continue
      itot& = 0      ' will start right after end select
   CASE ELSE
CALL cmessage(6, -1, 1, 15, 4, ibackcol, "Kein legaler Befehl, ggf. F1 benutzen", "ok")
END SELECT

'--- This is the start where you enter with pcontrol$="S"

js = KX: jy = ky: kax = 1: kay = 1
FOR il = 0 TO ilm: kax(il) = 1: KX(il) = js: kay(il) = 1: ky(il) = jy: NEXT
fran = KR: fran = fran / 100: fran5 = .5 * fran
sra = ra: IF ca > 0 THEN sra = ca
IF PCONTROL$ <> "II" THEN
	FOR iy = 1 TO imy
	FOR ix = 1 TO imx      'source density =random fluctuation
	ax(0, ix, iy) = sra * (1! + RND * fran - fran5)
	ax(9, ix, iy) = 1! + RND * fran - fran5
	NEXT ix: NEXT iy
END IF

FOR iy = 1 TO jy
FOR ix = 1 TO js     'set all cell to initial condtions
	ax(1, ix, iy) = gA
	ax(2, ix, iy) = gB
	IF ilm >= 3 THEN ax(3, ix, iy) = gC
	IF ilm >= 4 THEN ax(4, ix, iy) = gD
	IF ilm >= 5 THEN ax(5, ix, iy) = ge
	IF ilm >= 6 THEN ax(6, ix, iy) = gf
	IF ilm >= 7 THEN ax(7, ix, iy) = gg
	IF ilm >= 8 THEN ax(8, ix, iy) = gH
	NEXT ix: NEXT iy
jxx = (js + 1) / 2
jyy = (jy + 1) / 2
SELECT CASE KI
   CASE 1
      ax(1, jxx, jyy) = aA
      ax(3, jxx, jyy) = aC
   CASE 2
      ax(0, jxx, jyy) = aA * ax(0, jxx, jyy)
   CASE 3
      ax(1, jxx - 2, jyy - 2) = aA
      ax(4, jxx - 2, jyy - 2) = aD
   CASE 12  'for Hydra hypostome and tentakles
      FOR ix = 1 TO js
      FOR iy = 1 TO jy
      ax(3, ix, iy) = gC * EXP(-aC * (jxx - ix) ^ 2) * EXP(-aC * (jyy - iy) ^ 2)
      NEXT iy: NEXT ix
   CASE 13 ' fr Notochord
     FOR iy = k2 TO ky: ax(1, jxx, iy) = aa: NEXT
     ax(3, jxx, ky) = ac
     ax(3, jxx, 1) = ad

CASE 14:  'for insertion of new maxima during growth
   ax(1, 3, 3) = aA
   ax(1, js - 2, 3) = aA
   ax(1, js - 2, ky(1) - 2) = aA
   ax(1, 3, ky(1) - 2) = aA

   CASE 21
      ax(4, 3, 3) = aD
END SELECT
screeninit:
CALL init(ihp, igraph, ivorcol, ibackcol, ivgashift)
refresh:
CALL showdisplay(ax(), 1, KD, displaywas$)
IF label$(0) > "" THEN CALL cmessage(5, -1, 1, 15, 4, ibackcol, label$(0), "")

IF igraph = 9 THEN CALL switchscreen(2, iactive, ibackcol)
'--- Here you will enter with pcontrol$= "C"
continue:
IF icount >= KT OR PCONTROL$ = "D" OR PCONTROL$ = "I" OR PCONTROL$ = "II" GOTO EndSimulation
dra = 1 - ra - 4 * DA: drb = 1 - rB - 4 * DB
drc = 1 - rc - 4 * DC: drd = 1 - rd - 4 * DD
dre = 1 - re - 4 * DE: drf = 1 - rf - 4 * DF
drg = 1 - rg - 4 * DG: drh = 1 - rH - 4 * DH
IF KG > 0 THEN 'growth by insertion of additional cells
 igrowth = igrowth + 1
 IF igrowth > KG THEN
   igrowth = 1
     SELECT CASE K1
     CASE 0:
     kinsert = js
     CASE 1: kinsert = kax ' foot in hydra
     CASE 2: kinsert = (js + 1) / 2
     CASE 3
       IF js < 3 THEN
       kinsert = 1
     ELSE
       kinsert = js - 2   'sub-hypostome
     END IF
     CASE 4: kinsert = js: iyinsert = jy
     CASE 5:
	kinsert = (js + 1) / 2
	iyinsert = (kay(1) + ky(1)) / 2
     CASE ELSE:
CALL cmessage(6, -1, 1, 15, 4, ibackcol, "K1=Mode of growth, 0-5", "OK")
     END SELECT
       CALL manipulation("GX", kinsert)
       IF K1 > 3 THEN CALL manipulation("GY", iyinsert)
  END IF
END IF
  FOR iprint = 1 TO KP '-- SECOND LOOP: kp iterations until next plot ------
  ' Boundary conditions: impermeable; this is achieved by
  ' giving the cells kax - 1 and js+1 the same values as cells kax resp.js
  ' virtual left border cells:
   nx = js + 1
   ny = jy + 1
ixleft = 1
ixright = js

'boundary conditions
SELECT CASE K3
CASE 0, 1 'normal impermeable
   FOR il = 1 TO KN
	   FOR ix = kax TO js
	     ax(il, ix, ny) = ax(il, ix, jy)
	     u(il, ix) = ax(il, ix, 1): NEXT ix: NEXT il
CASE 2, 3, 4'Cylinder, wrapped around the y-axis
   FOR il = 1 TO KN
	   FOR ix = kax TO js
	     ax(il, ix, ny) = ax(il, ix, 1)
	     u(il, ix) = ax(il, ix, jy): NEXT ix: NEXT il
END SELECT
   FOR il = 1 TO KN  'boundary along y-axis
	   FOR iy = 1 TO jy
	   ax(il, nx, iy) = ax(il, ixright, iy)
	   NEXT iy: NEXT il
  '------ THIRD LOOP: cells kax...js get new values --------------

SELECT CASE KE'EQUATIONS 'selection of the equation according to ke

CASE 21 ' Activator-inhibitor system, sa=saturation
FOR iy = kay TO jy: GOSUB StoreLeftCell
FOR ix = 1 TO js: GOSUB olddecay
	aq = s * (a * a + ba)
	ax(1, ix, iy) = olddecaydiffA + aq / B
	ax(2, ix, iy) = olddecaydiffB + aq
NEXT: NEXT

CASE 23 ' Activator-inhibitor system, sa=saturation
'bb is Michaelis-Menten constant
FOR iy = kay TO jy: GOSUB StoreLeftCell
FOR ix = 1 TO js: GOSUB olddecay
	ax(1, ix, iy) = olddecaydiffA + s * (a * a / ((sb + B) * (1 + sa * a * a)) + ba)
	ax(2, ix, iy) = olddecaydiffB + s * a * a + bb
	NEXT: NEXT

CASE 24 ' activator-depleted substance
FOR iy = kay TO jy: GOSUB StoreLeftCell
FOR ix = 1 TO js: GOSUB olddecay
	aq = s * B * (a * a + ba)
	ax(1, ix, iy) = olddecaydiffA + aq / (1 + sa * a * a)
	ax(2, ix, iy) = olddecaydiffB + bb - aq
	IF ax(2, ix, iy) < 0 THEN ax(2, ix, iy) = 0
	NEXT: NEXT

CASE 28 ' Inhibition if an inhibition
FOR iy = kay TO jy: GOSUB StoreLeftCell
FOR ix = 1 TO js: GOSUB olddecay
	ax(1, ix, iy) = olddecaydiffA + s / (sa + C * C) + ba
	ax(2, ix, iy) = olddecaydiffB + rB * a
	ax(3, ix, iy) = olddecaydiffC + rc / (sc + a * a / (sb + B * B))
	NEXT: NEXT

CASE 122 'LI Net-like structure    (a,b) Activator-inhibitor
FOR iy = kay TO jy: GOSUB StoreLeftCell   '(c) substrate (d) differentiation
FOR ix = 1 TO js: GOSUB olddecay
	aq = s * C * a * a
	dq = D * D
	ax(1, ix, iy) = olddecaydiffA + aq / B + ba * D
	ax(2, ix, iy) = olddecaydiffB + aq + bb * D
	ax(3, ix, iy) = olddecaydiffC + bc - cC * C * D
	ax(4, ix, iy) = olddecaydiffD + rd * dq / (1 + sd * dq) + bd * a
	NEXT: NEXT

CASE 126 'Hydra head (a, b); source density (c), tentacles (d, e)
          ' Foot (f,g)
FOR iy = kay TO jy: GOSUB StoreLeftCell
FOR ix = 1 TO js: GOSUB olddecay
	aq = C * s * (a * a + ba)
	ax(1, ix, iy) = olddecaydiffA + aq / B
	ax(2, ix, iy) = olddecaydiffB + aq
	'Quelldichte
	ax(3, ix, iy) = olddecaydiffC + rc * a + bc - sc * C * f
	'Tentakeln
	dfq = D * D + bd
	dq = rd * C * ax(9, ix, iy) * dfq / (1 + sd * dfq) / (1 + cE * a)
	ax(4, ix, iy) = olddecaydiffD + dq / e
	ax(5, ix, iy) = olddecaydiffE + dq + be
	NEXT: NEXT

CASE 132 ' elongation of a stripe by a moving spot
FOR iy = kay TO jy: GOSUB StoreLeftCell
FOR ix = 1 TO js: GOSUB olddecay
	aq = s * (a * a + ba * c)
	ax(1, ix, iy) = olddecaydiffA + aq / B / (1 + sa * a * a)
	ax(2, ix, iy) = olddecaydiffB + s * a * a + bb
	cq = e * ax(9, ix, iy) * (c * c + bc * a)
	ax(3, ix, iy) = olddecaydiffC + rc * e * c * c / d
	ax(4, ix, iy) = olddecaydiffD + rd * cq + bd
	ax(5, ix, iy) = olddecaydiffE + be - se * e * a
	NEXT: NEXT

CASE 202 'closed loops by patch (ab) and stripe system (c,d)
  FOR iy = kay TO jy: GOSUB StoreLeftCell
  FOR ix = 1 TO js: GOSUB olddecay
	aq = s * B * (a * a + ba)
	ax(1, ix, iy) = olddecaydiffA + aq / (1 + sa * a * a)
	ax(2, ix, iy) = olddecaydiffB + bb - aq
	cq = rc * B ^ 2 * ax(9, ix, iy)  * (c * c + bc) / (d * (1 + sc * c * c))
	ax(3, ix, iy) = olddecaydiffC + cq
	ax(4, ix, iy) = olddecaydiffD + rc * c * c
	NEXT: NEXT


CASE ELSE
IF igraph = 9 THEN CALL switchscreen(3, iactive, ibackcol)
CALL cmessage(6, -1, 1, 15, 4, ibackcol, "No equation of this type, illegal KE ", "OK")
GOTO parameter
END SELECT 'ke
'-------------- End of equations ---------------------------------------------
IF INKEY$ = CHR$(27) GOTO makedraw
itot& = itot& + 1
NEXT iprint
icount = icount + 1
GOTO refresh

erroverflow:
CALL switchscreen(3, iactive, ibackcol)
mess$ = "Error, presumable an overflow, num. instability, zero-division"
CALL cmessage(6, -1, 1, 15, 4, ibackcol, mess$, "OK")
RESUME parameter

makedraw:
    CALL showdisplay(ax(), 1, KD, displaywas$)
    CALL switchscreen(2, iactive, ibackcol)

EndSimulation:
iparam = 2: ipcontrol = 0
CALL switchscreen(3, iactive, ibackcol)
GOTO parameter

END  '==================  End of Main-program  ==========

StoreLeftCell:  ' concentration in the left cell is stored in al...
    al = ax(1, ixleft, iy): bl = ax(2, ixleft, iy)
    cl = ax(3, ixleft, iy): dl = ax(4, ixleft, iy)
    IF KN > 4 THEN el = ax(5, ixleft, iy): fl = ax(6, ixleft, iy):
    gl = ax(7, ixleft, iy):   zhl = ax(8, ixleft, iy)
    RETURN

olddecay:
   a = ax(1, ix, iy): B = ax(2, ix, iy): s = ax(0, ix, iy)
  olddecaydiffA = dra * a + DA * (u(1, ix) + ax(1, ix, iy + 1) + al + ax(1, ix + 1, iy))
  olddecaydiffB = drb * B + DB * (u(2, ix) + ax(2, ix, iy + 1) + bl + ax(2, ix + 1, iy))
  al = a: bl = B: u(1, ix) = a: u(2, ix) = B
   IF KN > 2 THEN
     C = ax(3, ix, iy): D = ax(4, ix, iy)
     olddecaydiffC = drc * C + DC * (u(3, ix) + ax(3, ix, iy + 1) + cl + ax(3, ix + 1, iy))
     olddecaydiffD = drd * D + DD * (u(4, ix) + ax(4, ix, iy + 1) + dl + ax(4, ix + 1, iy))
     cl = C: dl = D: u(3, ix) = C: u(4, ix) = D
   IF KN > 4 THEN
     e = ax(5, ix, iy): f = ax(6, ix, iy): g = ax(7, ix, iy): zh = ax(8, ix, iy)
     olddecaydiffE = dre * e + DE * (u(5, ix) + ax(5, ix, iy + 1) + el + ax(5, ix + 1, iy))
     olddecaydiffF = drf * f + DF * (u(6, ix) + ax(6, ix, iy + 1) + fl + ax(6, ix + 1, iy))
     olddecaydiffG = drg * g + DG * (u(7, ix) + ax(7, ix, iy + 1) + gl + ax(7, ix + 1, iy))
     olddecaydiffH = drh * zh + DH * (u(8, ix) + ax(8, ix, iy + 1) + zhl + ax(8, ix + 1, iy))
     el = e: fl = f: gl = g: zhl = h
     u(5, ix) = e: u(6, ix) = f: u(7, ix) = g: u(8, ix) = g
   END IF
   END IF
 RETURN


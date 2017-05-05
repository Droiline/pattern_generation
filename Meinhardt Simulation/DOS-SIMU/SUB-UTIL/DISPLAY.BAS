DECLARE SUB lptopen (lpt$)
DECLARE SUB cmessage (itype%, iline%, irow%, ilettercol%, iletterback%, ibackcol%, text$, resp$)
DECLARE SUB init (ihp%, igraph%, ivorcol%, ibackcol%, ivgashift%)
DECLARE SUB linep (ihp%, x1!, y1!, x2!, y2!, ipat%, iwidth%, ico%, ityp%)
DECLARE SUB rwparameterdisk (rw$, ifilenumber%)
DECLARE SUB switchscreen (INZ%, iactive%, ibackcol%)
DECLARE SUB whatprmfiles (programname$, progrtmp$, exten$, drive$, irfilename$)
DECLARE SUB zeingabe (igt%, iquest%, INZ%, i%, f!, LABEL$, text$)
DECLARE SUB backfront (INZ%, x0!, x1!, y4!, y1!, s0r!, s1r!, s0l!, s1l!, iwi%, icol%, ihp%)
DECLARE SUB disp (il%, ax!(), krot%, INZ%, icol%, iwi%, ff!, xo!, yo!, kax%, KX%, kay%, ky%, tv$)
DECLARE SUB leftright (x1!, x2!, x3!, y1!, y2!, y3!, s1r!, s2r!, s3r!, s1l!, s2l!, s3l!, iwi%, icol%, ihp%)
DECLARE SUB perpendicular (x0!, y0!, ybase0!, s0l!, s0r!, iwi%, icol%, ihp%)
DECLARE SUB rl (ihp%, INZ%, x2!, y2!, s2!, x1!, yto!, s1!, y3!, msiluet%(), isi%, iwi%, icol%, iy%, jyto%)
DEFINT I-N
' $INCLUDE: 'sub-util\scommon.bi'
' $INCLUDE: 'sub-util\errlist.bas'

SUB backfront (INZ, x0, x1, y4, y1, s0r, s1r, s0l, s1l, iwi, icol, ihp)
SHARED p4, px, pyx, pdx, pdyx, pyy
IF y4 > s0l THEN 'line can start
  IF y1 > s1r THEN 'full
    CALL linep(ihp, x0, y4, x1, y1, 0, iwi, icol, 1)
  ELSE 'becomes invisible
    IF INZ < 5 OR INZ > 6 THEN
      p4 = px / (y4 - y1 + s1r - s0l)
      xs = (s1r - y1) * p4
      ys = y1 + xs * (y4 - y1) / px
    ELSE
      xs = (s1r - y1) * p4
      ys = y1 + xs * pyy / px
      CALL linep(ihp, x0, y4, x1 + xs, ys, 0, iwi, icol, 1)
    END IF
    CALL linep(ihp, x0, y4, x1 + xs, ys, 0, iwi, icol, 1)
  END IF
END IF
EXIT SUB
END SUB

SUB disp (il, ax(), krot, INZ, icol, iwi, ff, xo, yo, kax, KX, kay, ky, tv$) STATIC
SHARED icolm(), ibm(), a2(), ansp(), lpt$, icc, KG
SHARED ihp, itot&, DX, DY, per, LABEL$()
SHARED readdatafile$, writedatafile$, ibackcol, ivorcol, ivgashift
SHARED istelz, iletztk, ihintr
SHARED p4, px, py, pdx, pdyx, pyy, DZ
SHARED msiluet(), msiluetmax, ilm, icheckerror
DIM ty$(14)
IF KX = 0 OR ky = 0 THEN
mess$ = "Presumably you ask for a display before a first simulation, use S first"
CALL cmessage(6, -1, 1, 15, 4, ibackcol, mess$, "OK")
  EXIT SUB
END IF
IF il > ilm THEN
mess$ = "Too many substances, il or kn  <=" + STR$(ilm)
CALL cmessage(4, -1, 1, 15, 4, ibackcol, mess$, "OK")
  EXIT SUB
END IF
SELECT CASE INZ
CASE 0, 1, 2, 3, 4, 5, 6
  FOR i = 1 TO msiluetmax: msiluet(i) = 0: NEXT
  IF krot = 0 OR INZ = 5 OR INZ = 6 THEN krot = 2
  IF krot * KX + ky > msiluetmax THEN
    krot = (msiluetmax - ky) / KX
    mess$ = "KROT too large, korrected to " + STR$(krot)
    CALL cmessage(6, -1, 1, 15, 6, ibackcol, mess$, resp$)
  END IF
  winkel = ATN(1 / krot)
  pdx = DX * COS(winkel)
  px = pdx / krot
  pdyx = per * SIN(winkel) * DX
  pyx = pdyx / krot
  pyy = per * SQR(DX * DX - px * px)
  f1 = 40! * ff
  fmodus = f1 / krot
END SELECT

SELECT CASE INZ
CASE 0 TO 4' ================ Hidden line plot, 2=reverse
  IF krot = 0 THEN krot = 2
  jyfrom = kay: jyto = ky: iystep = 1
  IF INZ = 2 THEN
    pdyx = -pdyx: pyx = -pyx: pyy = -pyy: jyfrom = ky: jyto = kay: iystep = -1
  END IF
  kxfirstmem = krot * (KX - 1) + kay
  kxlastmem = KX
  FOR iy = jyfrom TO jyto STEP iystep
    kxlast = KX
    y2 = ax(il, KX, iy) * f1
    IF INZ <> 3 THEN 'skip zero-concentration in all other modes
      DO WHILE y2 = 0 AND kxlast > kax
	kxlast = kxlast - 1
	y2 = ax(il, kxlast, iy) * f1
      LOOP
    END IF
    IF kxlastmem <= kxlast AND kxlastmem > kax AND INZ = 1 THEN '-------basic lines after zero-skip,
    'isi = krot * (kxlast - 1) + iy
    's2 = msiluet(isi)
    'hier mu0 noch die Bremse hinein wenn es keliner geworden war und wieder gr”áer wurde
      ybaseto = yo - pdyx * (kxlast) + (iy - 1) * pyy
      x2 = xo + px * iy + pdx * (kxlast - 1)
      y1 = ybaseto + (kxlast - kxlastmem) * pdyx
      x1 = x2 - (kxlast - kxlastmem) * pdx
      IF kxlastmem < kxlast THEN
	CALL linep(ihp, x2, ybaseto, x1, y1, 0, iwi, icol, 1)'inner corner
	CALL linep(ihp, x1, y1, x1, y1 + ax(il, kxlastmem, iy) * f1, 0, iwi, icol, 1)
	CALL backfront(INZ, x1, x1 - px, y1, y1 - pyy, 0, 0, 0, 0, iwi, icol, ihp)
      ELSE '= no change of size
       IF iy > jyfrom THEN CALL backfront(INZ, x1, x1 - px, y1, y1 - pyy, 0, 0, 0, 0, iwi, icol, ihp)
      END IF
    END IF
    kxlastmem = kxlast
    'IF kxlast > kxlastmem THEN kxlastmem = kxlast
    ybaseto = yo - pdyx * (kxlast) + (iy - 1) * pyy
    y2 = y2 + ybaseto
    isi = krot * (kxlast - 1) + iy
    x2 = xo + px * iy + pdx * (kxlast - 1)
    s2 = msiluet(isi)
    CALL perpendicular(x2, y2, ybaseto, s2, s2, iwi, icol, ihp)
    ivisible = 1'may change if inz=3
    IF y2 > msiluet(isi) THEN msiluet(isi) = y2
    FOR ix = kxlast - 1 TO kax STEP -1
   IF ihintr < 1 THEN
   IF ax(1, ix, iy) > .5 THEN icol = 2
   IF ax(2, ix, iy) > .5 THEN icol = 1
   IF ax(3, ix, iy) > .5 OR ax(3, ix + 1, iy) > .5 THEN icol = 12
   END IF
      IF INZ > 2 THEN
	IF iy = jyfrom THEN ' front
	  CALL perpendicular(x2, y2, ybaseto, s2, s2, iwi, icol, ihp)
	ELSE
	  yto = ax(il, ix + 1, iy - 1) * f1
	  IF INZ = 4 AND yto = 0 THEN
	    CALL perpendicular(x2, y2, ybaseto, s2, s2, iwi, icol, ihp)
	  ELSE
	    y3 = ybaseto - pyy + yto
	    s1 = msiluet(isi - 1)
	    CALL backfront(INZ, x2, x2 - px, y2, y3, s2, s1, s2, s1, iwi, icol, ihp)
	  END IF
	END IF
      END IF
      'ybaseto = ybaseto + pdyx
      yto = ax(il, ix, iy) * f1
      IF yto = 0 THEN
       IF INZ <> 3 THEN
	IF INZ < 3 THEN CALL perpendicular(x2, y2, ybaseto, s2, s2, iwi, icol, ihp)
	kxfirst = isi
	GOTO nextiycase2
      END IF
      END IF
      ybaseto = ybaseto + pdyx
      yto = yto + ybaseto
      CALL rl(ihp, INZ, x2, y2, s2, x1, yto, s1, y3, msiluet(), isi, iwi, icol, iy, jyto)
    NEXT ix
    IF INZ > 2 AND iy > jyfrom THEN
      yto = ax(il, kax, iy - 1) * f1
      IF INZ = 4 AND yto = 0 THEN
	CALL perpendicular(x2, y2, ybaseto, s2, s2, iwi, icol, ihp)
      ELSE
	y3 = ybaseto - pyy + yto
	s1 = msiluet(isi - 1)
	CALL backfront(INZ, x2, x2 - px, y2, y3, s2, s1, s2, s1, iwi, icol, ihp)
      END IF
    ELSE
      CALL perpendicular(x2, y2, ybaseto, s2, s2, iwi, icol, ihp)
    END IF
    kxfirst = isi
nextiycase2:
    IF INZ = 4 AND kxfirst < kxfirstmem AND iy > kay THEN
      idiff = (kxfirstmem - kxfirst)
      x1 = x2 + px * idiff
      y1 = ybaseto - pyx * idiff
      CALL linep(ihp, x2, ybaseto, x1, y1, 0, iwi, icol, 1)
    END IF
    IF kxfirst < kxfirstmem THEN kxfirstmem = kxfirst
  NEXT iy
  irx = 0
  iry = 0
  GOTO frame


CASE 5, 6 'column-graphic
  ddx = .01 * DX
  w9 = ATN(.5)
  px = DX * COS(w9) * .5
  pyy = per * SQR(DX * DX - px * px)
  py = per * SIN(w9) * DX * .5
  pdx = 2! * px
  pdy = 2! * py
  sch = 2!
  p4 = px / (pyy + py)
  FOR iy = kay TO ky
    isi = 2 * KX + iy
    x1 = xo - 1.5 * px + px * isi
    yos = yo - (pyy + pdyx) / 2
    ybase1 = yos - pdy * (KX + 1) + (iy - 1) * pyy
    FOR ix = KX TO kax STEP -1
      x0 = x1 + px
      x2 = x1 - px
      x3 = x1 - pdx
      ybase3 = ybase1 + pdy
      ybase0 = ybase1 + pyy

      yhight = ax(il, ix, iy) * f1
      y1 = ybase1 + yhight
      y2 = y1 + py
      y3 = y1 + pdy
      y4 = y1 + pyy
      y5 = y4 + py
      y6 = y4 + pdy
      s0r = msiluet(isi + 1)
      s1r = msiluet(isi)
      s2r = msiluet(isi - 1)
      s3r = msiluet(isi - 2)
      isil = isi + 200
      s0l = msiluet(isil + 1)
      s1l = msiluet(isil)
      s2l = msiluet(isil - 1)
      s3l = msiluet(isil - 2)
      '   IF inz = 6 AND yhight > sch THEN
      '   ybase2 = ybase1 + dy
      '   CALL xleftright(x1, x2, x3, ybase1, ybase2, ybase3, s1r, s2r, s3r, s1l, s2l, s3l,iwi,icol,ihp)
      '   CALL xbackfront(x0, x1, ybase0, ybase1, s0r, s1r, s0l, s1l,iwi,icol,ihp)
      '   END IF
      
      IF INZ = 5 OR yhight > sch THEN
	CALL leftright(x1, x2, x3, y1, y2, y3, s1r, s2r, s3r, s1l, s2l, s3l, iwi, icol, ihp)
	CALL leftright(x0, x1, x2, y4, y5, y6, s0r, s1r, s2r, s0l, s1l, s2l, iwi, icol, ihp)
	CALL backfront(INZ, x0, x1, y4, y1, s0r, s1r, s0l, s1l, iwi, icol, ihp)
	CALL backfront(INZ, x2, x3, y6, y3, s2r, s3r, s2l, s3l, iwi, icol, ihp)
	CALL perpendicular(x0, y4, ybase0, s0l, s0r, iwi, icol, ihp)
	CALL perpendicular(x1, y1, ybase1, s1l, s1r, iwi, icol, ihp)
	CALL perpendicular(x3, y3, ybase3, s3l, s3r, iwi, icol, ihp)
      END IF
      'new silhouette
      IF y4 > msiluet(isil + 1) THEN msiluet(isil + 1) = y4
      IF y5 > msiluet(isi) THEN msiluet(isi) = y5
      IF y5 > msiluet(isil) THEN msiluet(isil) = y5
      IF y6 > msiluet(isi - 1) THEN msiluet(isi - 1) = y6
      IF y6 > msiluet(isil - 1) THEN msiluet(isil - 1) = y6
      IF y3 > msiluet(isi - 2) THEN msiluet(isi - 2) = y3
      x1 = x3
      ybase1 = ybase3
      isi = isi - 2
    NEXT ix
  NEXT iy
  irx = 1
  iry = 1
  GOTO frame

CASE 7
zylinder:  '-----------------------------------------------------
  iwil = iwi
  ra = DX * (ky - kay + 1) / 6.28
  fpy = 4 * DX * ff
  sq = DX * (1 - per * per)
  xpl = xo
  ypl = yo - sq * KX / 2
  IF krot = 5 THEN ypl = yo
  IF krot = 6 THEN ypl = yo - sq * KX
  aaa = 0: k1t = kax: kxt = KX
  IF krot >= 3 THEN
    aaa = ax(il, KX, 1)'Spitze
    kxt = KX - 1
  END IF
  aab = 0
  IF krot >= 4 THEN
    aab = ax(il, kax, 1)'Spitze
    k1t = kax + 1
  END IF
  wq = 6.28 / (ky - kay + 1)
  ' drawing of the circles
  FOR ix = k1t TO kxt
    '       nx = kxt + kax - ix
    iwile = iwil' thickness of lines, higher in the foreground
    ym = ypl + sq * (ix)
    angle = DY * wq
    pr = ax(il, ix, kay) * fpy + ra
    y1 = ym + per * pr * COS(angle)
    x1 = xpl + pr * SIN(angle)
    FOR iy = kay + 1 TO ky
      angle = ((iy - kay) + DY) * wq
      pr = ax(il, ix, iy) * fpy + ra
      y2 = ym + per * pr * COS(angle)
      x2 = xpl + pr * SIN(angle)
      IF COS(angle) < 0 OR ix = k1t THEN iwile = iwil + 1
      CALL linep(ihp, x1, y1, x2, y2, 0, iwile, icol, 1)
      'IF iwile > iwil THEN
      'END IF
      x1 = x2
      y1 = y2
    NEXT iy

    '       nx = kxt + k1t - ix
    '       nx = kxt + kax - ix
    angle = DY * wq
    pr = ax(il, ix, kay) * fpy + ra
    y2 = ym + per * pr * COS(angle)
    x2 = xpl + pr * SIN(angle)
    CALL linep(ihp, x1, y1, x2, y2, 0, iwile, icol, 1)
  NEXT ix

  ' drawing the lines
  'nx = kxt + kax - k1t
  yma = ypl + sq * (kxt) + aaa * fpy * (1 - per ^ 2)
  ymb = ypl + sq * k1t - aab * fpy * (1 - per ^ 2)
  FOR iy = kay TO ky
    angle = (iy - kay + DY) * wq
    iwile = iwil
    IF COS(angle) < 0 THEN iwile = iwil + 1
    x1 = xpl'lines start at the centre
    y1 = yma

    FOR ix = kxt TO k1t STEP -1
      '   nx = kxt + kax - ix
      ym = ypl + sq * (ix)
      angle = (iy - kay + DY) * wq
      pr = ax(il, ix, iy) * fpy + ra
      y2 = ym + per * pr * COS(angle)
      x2 = xpl + pr * SIN(angle)
      iwiler = iwile
      IF krot < 2 AND ix = KX THEN
	x1 = x2: y1 = y2
      ELSE
	CALL linep(ihp, x1, y1, x2, y2, 0, iwile, icol, 1)
      END IF
      x1 = x2
      y1 = y2
    NEXT ix
    IF krot >= 4 THEN CALL linep(ihp, x1, y1, xpl, ymb, 0, iwile, icol, 1)'second tip
  NEXT iy'

frame:
  '---front and left
 xos = xo: IF INZ = 5 OR INZ = 6 THEN xos = xo - 1.5 * px
 yos = yo: IF INZ = 5 OR INZ = 6 THEN yos = yo - (pyy + pdyx) / 2
  x1 = xos + px * (krot * (kax - 1) + kay)
  y1 = yos - pdyx * kax + pyy * (kay - 1)
  x2 = x1 + pdx * (KX - kax + irx)
  y2 = y1 - pdyx * (KX - kax + irx)
  'behind and right--------------
  x3 = x2 + px * (ky - kay + irx)
  IF INZ = 7 THEN x3 = x2 + px * (ky - kay)
  y3 = y2 + pyy * (ky - kay + iry)
  x4 = x1 + px * (ky - kay + irx)
  y4 = y1 + pyy * (ky - kay + iry)
  SELECT CASE INZ
  CASE 0, 1, 3, 4, 5, 6, 10
  CALL linep(ihp, x1, y1, x2, y2, 0, iwi, icol, 1)
    CALL linep(ihp, x2, y2, x3, y3, 0, iwi, icol, 1)
   IF ihintr < 1 THEN EXIT SUB
    IF ihintr > 0 THEN CALL linep(ihp, x3, y3, x4, y4, 0, 1, icol, 4)
    ver = ax(il, kax, ky) * f1
    IF ihintr > 1 THEN CALL linep(ihp, x4, y4, x4, y4 + ver, 0, 1, icol, 4)
    IF ihintr > 0 THEN CALL linep(ihp, x4, y4, x1, y1, 0, 1, icol, 4)

    IF INZ = 6 THEN
      CALL linep(ihp, x3, y3, x4, y4, 0, iwi, icol, 2)'
      CALL linep(ihp, x4, y4, x1, y1, 0, iwi, icol, 2)'
    END IF
    IF INZ = 10 THEN
      CALL linep(ihp, x3, y3, x4, y4, 0, iwi, icol, 1)'
      CALL linep(ihp, x4, y4, x1, y1, 0, iwi, icol, 1)'
    END IF
  CASE 2
'iwi>1 geht nicht
    IF ihintr = 1 THEN CALL linep(ihp, x2, y2, x1, y1, 0, 1, icol, 5)'back
    IF ihintr = 2 THEN CALL linep(ihp, x2, y2, x1, y1, 0, 1, icol, 1)
    CALL linep(ihp, x4, y4, x3, y3, 0, 1, icol, 1)'front
    CALL linep(ihp, x4, y4, x1, y1, 0, 1, icol, 1) 'left
    CALL linep(ihp, x2, y2, x3, y3, 0, 1, icol, 5) 'right
  END SELECT
  IF istelz > 0 AND iletztk > 10 THEN
    stelzm = iletztk - y1
    CALL linep(ihp, x1, y1, x1, y1 + stelzm, 0, 1, icol, 2)'Stelzen
    CALL linep(ihp, x2, y2, x2, y2 + stelzm, 0, 1, icol, 2)'Stelzen
    CALL linep(ihp, x3, y3, x3, y3 + stelzm, 0, 1, icol, 2)'Stelzen
    CALL linep(ihp, x4, y4, x4, y4 + stelzm, 0, 1, icol, 2)'Stelzen
  END IF
  iletztk = y1
  EXIT SUB


CASE 8, 9, 19'-------------------------------------pixel-------------------------------------
  dxt = (KX - kax + 1) * DX
  ddy = DX
  IF DY > 0 THEN ddy = DY
  dyt = (ky - kay + 1) * ddy
  IF xo = 0 THEN
  IF icc = 1 OR KG > 0 THEN
  xoa = 320 - (KX / 2 + kax - 1) * DX
  jam = kax
  END IF
  xa = xoa + DX * (kax - jam)
  ELSE
  xa = xo + DX * (kax - 1)
  END IF
  CALL linep(ihp, xa - 1, yo - 1, xa + dxt, yo + dyt + 1, 1, iwi, ivorcol, 1)
  CALL linep(ihp, xa - 2, yo - 2, xa + dxt + 1, yo + dyt + 2, 1, iwi, ivorcol, 1)
  FOR ix = kax TO KX
    FOR iy = kay TO ky
      icoly = icol
      IF INZ = 9 OR INZ = 19 THEN
      icoly = iy + icol - 1: IF icoly > 14 THEN icoly = icoly - 14
      icoly = ibm(icoly)
      END IF
      x1 = xa + (ix - kax) * DX
      y1 = yo + (iy - kay) * ddy
      x2 = x1 + DX
      y2 = y1 + ddy
      n = 8 * ax(il, ix, iy) * ff
      IF n < krot THEN n = 0
  IF n > 0 THEN CALL linep(ihp, x1, y1, x2, y2, 3, iwi, icoly, n)
    NEXT iy
  NEXT ix

CASE 12, 13
    am = -1
    FOR iy = kay TO ky
      FOR ix = kax TO KX
	IF ax(il, ix, iy) > am THEN am = ax(il, ix, iy)
      NEXT ix
    NEXT iy
    OPEN "SCRN:" FOR OUTPUT AS #8
    GOSUB ploxprint
    IF INZ = 13 THEN
      CALL lptopen("LPT1:"):
      IF icheckerror = 1 THEN EXIT SUB
      GOSUB ploxprint
    END IF

CASE ELSE

END SELECT' -----------------------------------------------------------------------

EXIT SUB



ploxprint: 'subroutine which can be reached with gosub
'Symbols for plot-labelling
ul$ = "É": um$ = "Í": ur$ = "»": ml$ = "º": mr$ = "º": ll$ = "È"
lm$ = "Í": lr$ = "¼"
ty$(0) = "-": ty$(1) = " ": ty$(2) = ".": ty$(3) = ":": ty$(4) = "°"
ty$(5) = "±": ty$(6) = "²": ty$(7) = "Û": ty$(8) = "Û": ty$(9) = "X"
ty$(10) = "X": ty$(11) = "X"

amx = 7 / (am + .00001): ram100 = 100 / (am + .0000001)
PRINT #8,
PRINT #8, tv$; : PRINT #8, "-Maximum = ";
IF am >= 1000 THEN
  PRINT #8, USING "#######.##"; am;
ELSE
  PRINT #8, USING "###.####"; am;
END IF
PRINT #8, USING "#######"; itot&; : PRINT #8, " iterations   ";
PRINT #8, " Field size: ";
PRINT #8, USING "###"; kax; KX; kay; ky
IF am < .0001 GOTO clploxprint
'------upper frame-----
PRINT #8, ul$;
FOR ix = kax TO KX: PRINT #8, um$;
NEXT ix
PRINT #8, ur$
IF KX - kax < 16 THEN zwraum$ = SPACE$(22 - KX + kax)
FOR iy = ky TO kay STEP -1
  PRINT #8, mr$;
  FOR ix = kax TO KX
    n = ax(il, ix, iy) * amx + 1.0001
    IF n <= 0 THEN n = 0
    PRINT #8, ty$(n);
    msiluet(ix) = ax(il, ix, iy) * ram100 + .4
  NEXT ix
  PRINT #8, mr$;
  IF KX - kax > 16 THEN
    PRINT #8,
  ELSE
    PRINT #8, zwraum$;
    FOR ix = kax TO KX: PRINT #8, USING "###"; msiluet(ix);
    NEXT ix
    PRINT #8,
  END IF
NEXT iy
PRINT #8, ll$;
FOR ix = kax TO KX: PRINT #8, lm$;
NEXT ix
PRINT #8, lr$
clploxprint:
CLOSE #8
RETURN



END SUB

SUB leftright (x1, x2, x3, y1, y2, y3, s1r, s2r, s3r, s1l, s2l, s3l, iwi, icol, ihp)
SHARED p4, px, py, pdx, pdy, pyy
IF y1 > s1l THEN 'line starts, point 1 (or 4) is visible
  IF y3 > s3r THEN 'point 3 (or 6) is visible
    CALL linep(ihp, x1, y1, x3, y3, 0, iwi, icol, 1)
  ELSE '3 not visible
    IF y2 > s2r THEN 'point 2 (or 5) is visible
      CALL linep(ihp, x1, y1, x2, y2, 0, iwi, icol, 1)
    END IF
  END IF
  
ELSE 'point 1 is not visible,
  IF y2 > s2r THEN 'point 2 (or 5) is visible
    IF y2 > s2l THEN
      xs = (s1l - y1) * p4
      ys = y1 + xs * py / px
    ELSE
      xs = px: ys = y2
    END IF
    IF y3 > s3r THEN 'point 3 (or 6) is visible
      CALL linep(ihp, x1 - xs, ys, x3, y3, 0, iwi, icol, 1)
    ELSE
      CALL linep(ihp, x1 - xs, ys, x2, y2, 0, iwi, icol, 1)
    END IF
  ELSE 'point 2 was not visible
    
    'CALL linep(ihp, x2 + 1, y2, x2 + 1, s2r, 2, iwi, 13, 1)
    IF y3 > s3r THEN 'point 3 (or 6) is visible
      IF y2 > s2l THEN
	CALL linep(ihp, x2, y2, x3, y3, 0, iwi, icol, 1)
      ELSE
	xs = (s2l - y2) * p4
	ys = y2 + xs * py / px
	CALL linep(ihp, x2 - xs, ys, x3, y3, 0, iwi, icol, 1)
	' PSET (x2, s2l), 11
	'PSET (x2, y2), 13
      END IF
      'PRINT y2, s2l, xs, px, y3
      'INPUT resp$
      
    END IF
  END IF
END IF
EXIT SUB
END SUB

SUB manipulation (mwasnun$, iinsert) STATIC
SHARED ax(), kax, KX, kay, jy, js, KN, imx, imy, kax(), KX(), kay(), ky()
SHARED mempos(), igt, ivgashift, icc
IF mwasnun$ = "" THEN
CLS
manipulationlabel:
iinsert = 0
LOCATE 1, 1
PRINT " What shall be done?"
PRINT
PRINT "  To ...                                            type"
PRINT
PRINT "     store the current state                         s"
PRINT "     restore the stored state                        r"
PRINT "     change field size                               c"
PRINT "     set local alterations of array s, a, b...       a"
PRINT "     double all cells in X-dir.                      d"
PRINT "     growth at the right margin                      gm"
PRINT "     wipe a and d in a part of the field             0"
PRINT "     insert a column in the field (ix)               gx"
PRINT "     insert a row in the field (iy)                  gy"
PRINT "     remove a piece from the field                   de"
PRINT "     transplant or rearrange pieces in X-dir.       t"
PRINT "     restore a graft combination                     sg"
PRINT
PRINT "     go back to the parameter menu                   RETURN"
PRINT
INPUT resp$: resp$ = UCASE$(resp$)

ELSE
resp$ = mwasnun$
END IF
SELECT CASE resp$

CASE "GM"     'growth at the right margin
      iinsert = js
      GOSUB growthxy

CASE "CA", "A"
il = -1: ix = 0: iy = 0
CALL zeingabe(igt, 0, 1, il, dummy, resp$, "Substance to be changed, 1= a, 2=b, ect.")
IF il = -1 THEN EXIT SUB
CALL zeingabe(igt, 0, 1, ix, dummy, resp$, "ix of the cell to be changed")
IF ix = 0 THEN EXIT SUB
CALL zeingabe(igt, 0, 1, iy, dummy, resp$, "iy of the cell to be changed")
IF iy = 0 THEN EXIT SUB
CALL zeingabe(igt, 1, 2, ix, ax(il, ix, iy), resp$, "new value")

CASE "S"        'to store present state of cells in memory arrays
      GOSUB storearray

CASE "R"        'to restore former state of cells from memory
      GOSUB restorearray
CASE "C"'to change margins of the field
CALL zeingabe(igt, 0, 3, kanew, dummy, resp$, "Store present or Restore previous concentrations S/R/RETURN=nothing")
  SELECT CASE resp$
  CASE "S", "s"
      GOSUB storearray
  CASE "r", "R"
      GOSUB restorearray
  END SELECT
  kanew = kax: ksnew = js
mess$ = " Input new (smaller) field size: Present field size is: " + STR$(kax) + " - " + STR$(js)
CALL cmessage(3, -2, 1, 15, 4, ibackcol, mess$, "OK")
CALL zeingabe(igt, 1, 1, kanew, dummy, dummy$, "New left-most cell")
CALL zeingabe(igt, 1, 1, ksnew, dummy, dummy$, "New right-most cell")
CALL cmessage(3, -2, 1, 15, 4, ibackcol, "", "OK")
  IF ksnew = 0 THEN ksnew = js
  IF kanew < kax THEN kanew = kax
  js = ksnew
  kax = kanew
      FOR il = 1 TO KN: KX(il) = js: kax(il) = kax: NEXT


CASE "DE"      'to remove a piece (delete)
  irem1 = kax: rem2 = js
  mess$ = " present field size is now " + STR$(kax) + "; " + STR$(js)
  CALL cmessage(3, -1, 1, 15, 4, ibackcol, mess$, "OK")
  mess$ = " Position of leftmost cell of piece to be removed:  "
  CALL zeingabe(igt, 0, 1, irem1, dummy, dummy$, mess$)
  mess$ = " Position rightmost cell of piece to be removed:  "
  CALL zeingabe(igt, 0, 1, irem2, dummy, dummy$, mess$)
      iremnumber = irem2 - irem1 + 1 'number of cells to be removed
      FOR il = 1 TO KN
      KX(il) = js - iremnumber
      FOR ix = (irem2 + 1) TO js
	 FOR iy = 1 TO ky(il)
	    ax(il, ix - iremnumber, iy) = ax(il, ix, iy)
      NEXT iy: NEXT ix: NEXT il
      js = js - iremnumber

CASE "0", "0C" 'to set a and d (activator and differentiation to zero)
mess$ = " X-region to be wiped out; Present field size is: " + STR$(kax) + " - " + STR$(js)
kanew = 1: ksnew = js: kaynew = 1: kynew = jy
CALL cmessage(6, -1, 1, 15, 4, ibackcol, mess$, "OK")
CALL zeingabe(igt, 1, 1, kanew, dummy, dummy$, "X-left-most cell")
CALL zeingabe(igt, 1, 1, ksnew, dummy, dummy$, "X-right-most cell")
  IF ksnew = 0 OR ksnew > js THEN ksnew = js
  IF kanew < 1 THEN kanew = 1

mess$ = "Y-region to be wiped out; Present field size is: " + STR$(kay) + " - " + STR$(jy)
CALL cmessage(6, -1, 1, 15, 4, ibackcol, mess$, "OK")
CALL zeingabe(igt, 1, 1, kaynew, dummy, dummy$, "Y-left-most cell")
CALL zeingabe(igt, 1, 1, kynew, dummy, dummy$, "Y-right-most cell")
  IF kynew < 1 OR kynew > jy THEN kynew = jy
  IF kaynew < 1 THEN kaynew = kay

FOR ix = kanew TO ksnew
FOR iy = kaynew TO kynew
ax(1, ix, iy) = 0
IF resp$ = "0" THEN
 ax(4, ix, iy) = 0:
ELSE
 ax(3, ix, iy) = 0
END IF
NEXT iy: NEXT ix

CASE "D"
jinsert = iinsert
FOR jinsert = js TO 1 STEP -1
GOSUB growthxy
NEXT

CASE "GX"
IF iinsert = 0 THEN
 mess$ = " A row of cell will be inserted; position?, 0=centre"
CALL zeingabe(igt, 1, 1, iinsert, dummy, dummy$, mess$)
IF iinsert = 0 THEN iinsert = (kax(1) + KX(1)) / 2
END IF
GOSUB growthxy


CASE "GY"
IF iinsert = 0 THEN
mess$ = " A column of cells will be inserted; position?, 0=centre"
CALL zeingabe(igt, 1, 1, iinsert, dummy, dummy$, mess$)
IF iinsert = 0 THEN iinsert = (kay(1) + ky(1)) / 2
END IF
IF jy < imy - 1 THEN
FOR il = 1 TO KN
FOR iy = jy TO iinsert STEP -1          'still the old js is being used!
FOR ix = 1 TO KX(il)
      ax(il, ix, iy + 1) = ax(il, ix, iy)
   NEXT ix:  NEXT iy
   ky(il) = ky(il) + 1
   NEXT il

jy = jy + 1                             'finally js is adjusted to the new state
ELSE
mess$ = "array becomes too large,  js>imx [" + STR$(js) + ", " + STR$(imx) + " ]"
CALL cmessage(3, -1, 1, 15, 4, ibackcol, mess$, "OK")
END IF

CASE "GR", "T"'graft pieces together
  CALL zeingabe(igt, 0, 3, kanew, dummy, resp$, "Store present or Restore previous concentrations S/R/RETURN=nothing")
  SELECT CASE resp$
  CASE "s", "S"
      GOSUB storearray
  CASE "r", "R"
      GOSUB restorearray
  END SELECT
      IF ksmem = 0 THEN GOSUB storearray
      FOR ix = 1 TO imx
	 mempos(ix) = 0 'first mop up anything left from previous manipulations
      NEXT ix
      jxright = 0
      IF igt = 0 THEN
      FOR i = 18 + ivgashift TO 24 + ivgashift:
      LOCATE i, 1: PRINT SPACE$(80); : NEXT
      LOCATE 18 + ivgashift, 1
      PRINT " You can combine cells from the stored state to create a new field."
      PRINT " Each cell can be used in the new field several times."
      PRINT " The stored field has"; ksmem; "cells."
      PRINT " If you want to take just one cell, return two identical positions. ";
      PRINT "Stored field size is"; kamem; " - "; ksmem; ""
      END IF
graftarray:
  mess$ = " starting cell of the piece      (or RETURN = end of grafting process)"
  CALL zeingabe(igt, 0, 1, ifrom, dummy, dummy$, mess$)
  IF ifrom = 0 THEN GOTO graftassembly
  mess$ = "  end cell of the piece"
  CALL zeingabe(igt, 0, 1, ito, dummy, dummy$, mess$)
      itofromdiff = ABS(ito - ifrom)
	 IF (ileftx + 1 + itofromdiff) > imx THEN
CALL cmessage(3, -1, 1, 15, 4, ibackcol, " Piece too long, has been cuttet", "OK")
	    itofromdiff = imx - ileftx - 1: ito = ifrom + ifromdiff
	 END IF
	 istep = 1: IF ifrom > ito THEN istep = -1
	 FOR ix = ifrom TO ito STEP istep
	 jxright = jxright + 1
	    mempos(jxright) = ix
	 NEXT ix
      IF igt = 0 THEN
      FOR i = 18 + ivgashift TO 24 + ivgashift:
      LOCATE i, 1: PRINT SPACE$(80); : NEXT
      LOCATE 18 + ivgashift, 1
      PRINT "Stored field size is"; kamem; " - "; ksmem; ""
	 PRINT "present field size is 1 - "; jxright;
	 PRINT " The field contains for far the following cells:"; jxright; "in total"
	 FOR jxx = 1 TO jxright
	    PRINT mempos(jxx);
	 NEXT jxx
	 END IF
	 js = jxright
      GOTO graftarray
CASE "SG"
      GOTO graftassembly

   CASE ELSE
      FOR il = 1 TO KN: KX(il) = js: NEXT
      EXIT SUB
END SELECT
EXIT SUB


growthxy:
IF js < imx THEN
FOR il = 1 TO KN:
FOR ix = KX(il) TO iinsert STEP -1          'still the old js is being used!
     FOR iy = 1 TO ky(il)
      ax(il, ix + 1, iy) = ax(il, ix, iy)
   NEXT iy:  NEXT ix
   KX(il) = KX(il) + 1
   NEXT il
js = js + 1
ELSE
mess$ = "array becomes too large,  js>imx [" + STR$(js) + ", " + STR$(imx) + " ]"
CALL cmessage(3, -1, 1, 15, 4, ibackcol, mess$, "OK")
END IF
RETURN


storearray:
IF 2 * jy > imy THEN
mess$ = "Array too large to be stored, (jy too large or imy too small)"
CALL cmessage(3, -1, 1, 15, 12, ibackcol, mess$, "OK")
ELSE
ksmem = js: kymem = jy: kamem = kax
kyadd = imy / 2
      FOR il = 1 TO KN: FOR ix = kax TO js: FOR iy = kay TO jy
      ax(il, ix, iy + kyadd) = ax(il, ix, iy)
      NEXT iy: NEXT ix: NEXT il
CALL cmessage(3, -1, 1, 15, 10, ibackcol, " State has been stored", "OK")
END IF
RETURN


restorearray:
IF ksmem = 0 THEN
mess$ = "sorry, nothing was stored"
CALL cmessage(6, -1, 1, 15, 4, ibackcol, mess$, "OK")
ELSE
js = ksmem: jy = kymem: kax = kamem: icc = 1
      FOR il = 1 TO KN:
      kax(il) = kamem: KX(il) = js: ky(il) = jy
      FOR ix = kax TO js: FOR iy = 1 TO jy
      ax(il, ix, iy) = ax(il, ix, iy + kyadd)
	 NEXT iy: NEXT ix: NEXT il
mess$ = " State has been restored"
CALL cmessage(3, -1, 1, 15, 2, ibackcol, mess$, "OK")
END IF
RETURN

graftassembly:
       'Now all pieces are being assembled
      js = jxright
      KX(1) = js
      jy = kymem
      FOR ix = 1 TO js
      FOR iy = 1 TO jy
      FOR il = 1 TO KN
	    ax(il, ix, iy) = ax(il, mempos(ix), iy + kyadd)
	 NEXT il
      NEXT iy
      NEXT ix
      FOR il = 1 TO KN: KX(il) = js: NEXT
      EXIT SUB
       'Now you can see what you have done

END SUB

SUB perpendicular (x0, y0, ybase0, s0l, s0r, iwi, icol, ihp)
SHARED p4, px, py, pdx, pdy, pyy


IF y0 > s0r THEN 'line starts
  yba = ybase0
  IF s0l > yba THEN yba = s0l
  IF s0r > yba THEN yba = s0r
  IF yba < y0 THEN CALL linep(ihp, x0, y0, x0, yba, 0, iwi, icol, 1)
END IF
END SUB

SUB rl (ihp, INZ, x2, y2, s2, x1, yto, s1, y3, msiluet(), isi, iwi, icol, iy, jyto)
SHARED px, pyx, pdx, pdyx, pyy, krot, ihintr
'line from left to right
dym = (yto - y2) / krot
IF y2 >= s2 THEN xstart = x2: ystart = y2
FOR i = 1 TO krot'------Interpolation-------
  isi = isi - 1'new x-coordinate of the horizon
  s1 = msiluet(isi)'new y-coordinate of the horizon at point 1
  x1 = x2 - px
  y1 = y2 + dym
  IF y2 >= s2 THEN ' line is not hidden
    msiluet(isi + 1) = y2
    IF i = 1 AND y3 > y1 AND INZ > 2 THEN
      ivisible = 0
      'xstart = x1: ystart = y1
      GOTO nextimodus
    END IF
    ivisible = 1
    IF y1 < s1 THEN
      xs = px * (s1 - y1) / (y2 + s1 - y1 - s2 + .000001)' line becomes invisible
      ys = y1 + xs * (y2 - y1) / px'point where kp becomes hidden
      xs = x1 + xs:
      ivisible = 0
      CALL linep(ihp, xstart, ystart, xs, ys, 0, iwi, icol, 1)'    , icol
    ELSE
      msiluet(isi) = y1
    END IF
  ELSE 'point 2 was invisible
    IF y1 > s1 THEN 'point 1 becomes visible
      xs = px * (s1 - y1) / (y2 + s1 - y1 - s2 + .000001)' line becomes invisible
      ystart = y1 + xs * (y2 - y1) / px'point where kp becomes hidden
      xstart = x1 + xs:
      msiluet(isi) = y1
      ivisible = 1

    ELSE 'point 1 and 2was invisible,last kurve marked by ihintr>2
      IF iy = jyto AND ihintr > 2 THEN
	CALL linep(ihp, x2, y2, x2, y2, 0, iwi, icol, 1)'    , icol'  line after krot-loop
      END IF
    END IF
  END IF
nextimodus:
  y2 = y1'old left point becomes the new right point
  x2 = x1
  s2 = s1
NEXT i
IF ivisible = 1 THEN CALL linep(ihp, xstart, ystart, x1, y1, 0, iwi, icol, 1)
EXIT SUB
END SUB

SUB showdisplay (ax(), iwas, kdisplay, displaywas$) STATIC
SHARED kax(), KX(), kay(), ky(), xa(), fa(), am(), ya(), f1(), icol()
SHARED icolm(), ibm(), a2(), ansp(), lpt$, js, jy
SHARED ityp(), iwi(), tv$(), afa(), ito(), LABEL$(), title$, msiluet()
SHARED ihp, dxp, dyp, DX, DY, per, krot, ivgashift, igraph, ianz, ireaddata
SHARED istelz, iletztk, ihintr, iparam, iactive, iwidthm
SHARED ivorcol, ibackcol, dataname$, readdatafile$, writedatafile$
SHARED kax, KX, kay, ky, itot&, KN, igt
ibm(1) = 1: ibm(2) = 2: ibm(3) = 12: ibm(4) = 4: ibm(5) = 5: ibm(6) = 12: ibm(7) = 13
ibm(8) = 8: ibm(9) = 9: ibm(10) = 10
ibm(11) = 11: ibm(12) = 12: ibm(13) = 13: ibm(14) = 14: ibm(15) = 15
WINDOW (0, 0)-(639, 479)
D = DX: iletztk = 0
ipw = INSTR(displaywas$, "w")
IF displaywas$ = "" AND kdisplay <> 27 THEN
IF igraph = 9 THEN CALL switchscreen(3, iactive, ibackcol)
CALL zeingabe(igt, 0, 3, i, D, displaywas$, "Substances to be displayed, e.g, abd..")
displaywas$ = LCASE$(displaywas$)
END IF


twas$ = "sabcdefgh"
FOR il = 0 TO KN
  ipu = 1: ipu = INSTR(ipu, displaywas$, MID$(twas$, il + 1, 1))
  IF ipu > 0 THEN
    tv$(il) = MID$(twas$, il + 1, 1)
    ity(il) = ityp(il): ist = il
    IF kdisplay > 0 THEN ity(il) = kdisplay
  ELSE
    ity(il) = 0
  END IF
NEXT il

SELECT CASE kdisplay
CASE 0 TO 9, 12, 13
    IF kdisplay < 12 THEN
      LINE (0, 0)-(640, 480), ibackcol, BF'initialization
      ELSE
    LOCATE 1, 1
    END IF


twas$ = "sabcdefghij"
FOR ilv = LEN(displaywas$) TO 1 STEP -1
    il = ASC(MID$(displaywas$, ilv, 1)) - 96'a=asc(97)
    IF il <> 23 THEN
    IF il = 19 THEN il = 0
    tv$(il) = MID$(twas$, il + 1, 1)
    ity(il) = ityp(il): ist = il
    IF kdisplay > 0 THEN ity(il) = kdisplay
    CALL disp(il, ax(), krot, ity(il), icol(il), iwi(il), fa(il), xa(il), ya(il), kax(il), KX(il), kay(il), ky(il), tv$(il))
    END IF
  NEXT ilv

CASE 10 ' for netlike Strukture
  LINE (0, 0)-(640, 480), ibackcol, BF
  dyp = 400: d2 = 400! / ky
  d21 = d2 - 1!: d22 = d2 / 2!
  y2 = 50!: x2 = 300! - KX / 2 * d2
  FOR ix = 1 TO KX
    FOR iy = 1 TO ky
      dxt = x2 + d2 * (ix - 1)
      dyt = y2 + d2 * (iy - 1)
      IF ax(1, ix, iy) > 1.8 THEN LINE (dxt + 1, dyt + 1)-(dxt + d21, dyt + d21), icol(1), BF
      IF ax(4, ix, iy) > .8 THEN LINE (dxt, dyt)-(dxt + d2, dyt + d2), icol(4), B
      IF ax(3, ix, iy) > .8 THEN
	LINE (dxt, dyt)-(dxt + d22, dyt + d2), icol(3)
	LINE (dxt + d22, dyt + d2)-(dxt + d2, dyt), icol(3)
      END IF
    NEXT: NEXT
ipu = 1: ipu = INSTR(ipu, displaywas$, "h")

IF ipu > 0 THEN
  FOR ix = 1 TO KX
    FOR iy = 1 TO ky
      dxt = x2 + d2 * (ix - 1)
      dyt = y2 + d2 * (iy - 1)
      IF ax(5, ix, iy) > 1.8 THEN LINE (dxt + 1, dyt + 1)-(dxt + d21, dyt + d21), icol(5), BF
      IF ax(8, ix, iy) > .8 THEN LINE (dxt, dyt)-(dxt + d2, dyt + d2), icol(8), B
      IF ax(7, ix, iy) > .8 THEN
	LINE (dxt, dyt + d2)-(dxt + d22, dyt), icol(7)
	LINE (dxt + d22, dyt)-(dxt + d2, dyt + d2), icol(7)
      END IF
    NEXT: NEXT
END IF
  LINE (x2 - 2, y2 - 2)-(x2 + d2 * KX + 2, y2 + d2 * ky + 2), ivorcol, B
  LINE (x2 - 3, y2 - 3)-(x2 + d2 * KX + 3, y2 + d2 * ky + 3), ivorcol, B




CASE 11 'colored random dots, density proportional to level
  xo = xa(ist): yo = ya(ist)
  LINE (0, 0)-(640, 480), ibackcol, BF
  inkd = DX: kax = 1: kay = 1
  dxx = (KX - kax + 1) * inkd
  dyy = (ky - kay + 1) * inkd
  CALL linep(ihp, xo - 1, yo - 1, xo + dxx + 2, yo + dyy + 2, 2, 1, ibackcol, 1)
  CALL linep(ihp, xo - 1, yo - 1, xo + dxx + 2, yo + dyy + 2, 1, 1, ivorcol, 1)
  CALL linep(ihp, xo, yo, xo + dxx + 1, yo + dyy + 1, 1, 1, ivorcol, 1)
  yb = yo + (iy - kay) * inkd
  inc = 0
  FOR ix = kax TO KX
    FOR iy = kay TO ky
      xb = xo + (ix - kax) * DX + 1
      yb = yo + (iy - kay) * DX + 1
      su = 0
      inc = 0
      FOR il = 0 TO KN
	f1 = fa(il) * .3
	IF ity(il) = 11 THEN
	  inc = inc + 1
	  icolm(inc) = icol(il)
	  a2(inc) = ax(il, ix, iy) * f1
	  su = su + a2(inc)
	END IF
      NEXT il

      IF su > .01 THEN
	anc = 0
	FOR ibc = 1 TO inc'normieren
	  IF su < 1! THEN
	    anc = anc + a2(ibc)
	  ELSE
	    anc = anc + a2(ibc) / su
	  END IF
	  ansp(ibc) = anc'anteil summativ
	NEXT ibc
	FOR ikx = 1 TO inkd
	  ybb = yb
	  FOR iky = 1 TO inkd
	    rra = RND
	    IF rra < ansp(inc) THEN
	      i = 1
	      DO UNTIL rra < ansp(i)
		i = i + 1
	      LOOP
	      PSET (xb, ybb), icolm(i)
	    END IF
	    ybb = ybb + 1!
	  NEXT iky
	  xb = xb + 1!
	NEXT ikx
      END IF
    NEXT iy: NEXT ix

END SELECT
EXIT SUB
END SUB


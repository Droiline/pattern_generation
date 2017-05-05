SUB ColorPrint (text$, row%, col%, fc%, bc%) STATIC
SHARED igraph
IF igraph < 12 THEN
IF bc% < 8 THEN
COLOR bc% + 8
ELSE
COLOR bc%
END IF
LOCATE row%, col%: PRINT text$;
COLOR 15
EXIT SUB
END IF
'This creates the background
c$ = STRING$(LEN(text$), 219)'Defines a Big Block BackGround(really foreground)
reg 1, &H1300
reg 2, bc%
reg 3, LEN(text$)
reg 4, 256 * (row% - 1) + (col% - 1)
reg 9, strseg(c$)     'es
reg 7, strptr(c$)     'bp
call interrupt &h10
reg 1, &H1300
reg 2, (fc% XOR bc%) + &H80 'xor the previous color block with new characters
reg 3, LEN(text$)
reg 4, 256 * (row% - 1) + (col% - 1)
reg 9, strseg(text$)
reg 7, strptr(text$)
call interrupt &h10
END SUB


$INCLUDE "sub-util\PROFILE.bas"
$segment
$INCLUDE "sub-util\plib.bas"
$segment

'======================================================

GOTO jumpover
'=========the next line is necessary and allowed only for POWER-BASIC======
$include "sub-util\errlist.bas"
jumpover:




SUB ColorPrint (text$, row%, col%, fc%, bc%) STATIC
'SHARED igraph
IF fc% < 8 THEN
COLOR fc% + 8
ELSE
COLOR fc%
END IF
IF col% > 78 THEN col% = 78
LOCATE row%, col%: PRINT text$;
COLOR 15
END SUB


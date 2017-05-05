DECLARE SUB message (itype%, iline%, text$)
printerproblem:
mess$="  =======   Printer needs attention,  ============"
CALL cmessage (4,-1,1,15,4, 0,mess$, "Type Return")
icheckerror = 1
RESUME NEXT

checkerror:
icheckerror = 1
RESUME NEXT


wrongscreen:
tryanother = 1
RESUME NEXT


nosuchfile:
RESUME NEXT

checkfile:
icheckerror = 1
SELECT CASE ERR
CASE 25
mess$= readdatafile$ + " Hardware-error "
CASE 52
mess$= readdatafile$ + " Wrong file name or number ": icheckerror = 2
CASE 53
  mess$= readdatafile$ + " File does not exist "
CASE 62
    mess$= readdatafile$  + " End of file is reached, file closed "
END SELECT
mess$=mess$ + "  Err.-code= " + STR$(ERR)
CALL cmessage (4,-1,1,15,4, ibackcol,mess$, "Type Return")
CLOSE #27
readdatafile$ = ""
RESUME NEXT


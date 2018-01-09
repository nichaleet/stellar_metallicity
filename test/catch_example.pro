PRO CATCH_EXAMPLE
 
   ; Define variable A:
   A = FLTARR(10)
 
   ; Establish error handler. When errors occur, the index of the
   ; error is returned in the variable Error_status:
   CATCH, Error_status
 
   ;This statement begins the error handler:
   IF Error_status NE 0 THEN BEGIN
      PRINT, 'Error index: ', Error_status
      PRINT, 'Error message: ', !ERROR_STATE.MSG
      ; Handle the error by extending A:
      A=FLTARR(12)
      CATCH, /CANCEL
   ENDIF
   stop
   ; Cause an error:
   A[11]=12
 
   ; Even though an error occurs in the line above, program
   ; execution continues to this point because the event handler
   ; extended the definition of A so that the statement can be 
   ; re-executed.
   HELP, A
END

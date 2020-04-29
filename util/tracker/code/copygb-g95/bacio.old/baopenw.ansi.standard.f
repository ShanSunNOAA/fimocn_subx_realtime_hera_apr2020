C-----------------------------------------------------------------------
      SUBROUTINE BAOPENWT(LU,CFN,IRET)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM: BAOPENWT
C
C$$$
      CHARACTER CFN*(*)
      OPEN(UNIT=LU,FILE=CFN,form='unformatted')
      rewind(LU)
      return
      END

      SUBROUTINE BAOPENWA(LU,CFN,IRET)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C
C SUBPROGRAM: BAOPENWA
C
C$$$
      CHARACTER CFN*(*)
      OPEN(UNIT=LU,FILE=CFN,form='unformatted',position="append")
      rewind(LU)
      return
      end

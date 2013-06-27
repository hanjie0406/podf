C*****precision > double
      SUBROUTINE DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
*     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
*     ..
*
*  Purpose
*  =======
*
*  DGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X',
*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*
*              TRANSA = 'T' or 't',  op( A ) = A'.
*
*              TRANSA = 'C' or 'c',  op( A ) = A'.
*
*           Unchanged on exit.
*
*  TRANSB - CHARACTER*1.
*           On entry, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*
*              TRANSB = 'T' or 't',  op( B ) = B'.
*
*              TRANSB = 'C' or 'c',  op( B ) = B'.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      DOUBLE PRECISION   TEMP
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Executable Statements ..
*
*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
*     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
*     and  columns of  A  and the  number of  rows  of  B  respectively.
*
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
*
*     Test the input parameters.
*
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.
     $         ( .NOT.LSAME( TRANSB, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     And if  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
*
*     Start the operations.
*
      IF( NOTB )THEN
         IF( NOTA )THEN
*
*           Form  C := alpha*A*B + beta*C.
*
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B + beta*C
*
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
      ELSE
         IF( NOTA )THEN
*
*           Form  C := alpha*A*B' + beta*C
*
            DO 170, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 130, I = 1, M
                     C( I, J ) = ZERO
  130             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 140, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  140             CONTINUE
               END IF
               DO 160, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 150, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  150                CONTINUE
                  END IF
  160          CONTINUE
  170       CONTINUE
         ELSE
*
*           Form  C := alpha*A'*B' + beta*C
*
            DO 200, J = 1, N
               DO 190, I = 1, M
                  TEMP = ZERO
                  DO 180, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  180             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  190          CONTINUE
  200       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DGEMM .
*
      END
      SUBROUTINE DGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
     $                   BETA, Y, INCY )
*     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA, BETA
      INTEGER            INCX, INCY, LDA, M, N
      CHARACTER*1        TRANS
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  DGEMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n matrix.
*
*  Parameters
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*
*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
*
*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry with BETA non-zero, the incremented array Y
*           must contain the vector y. On exit, Y is overwritten by the
*           updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 1
      ELSE IF( M.LT.0 )THEN
         INFO = 2
      ELSE IF( N.LT.0 )THEN
         INFO = 3
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
      IF( LSAME( TRANS, 'N' ) )THEN
         LENX = N
         LENY = M
      ELSE
         LENX = M
         LENY = N
      END IF
      IF( INCX.GT.0 )THEN
         KX = 1
      ELSE
         KX = 1 - ( LENX - 1 )*INCX
      END IF
      IF( INCY.GT.0 )THEN
         KY = 1
      ELSE
         KY = 1 - ( LENY - 1 )*INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y.
*
      IF( BETA.NE.ONE )THEN
         IF( INCY.EQ.1 )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 10, I = 1, LENY
                  Y( I ) = ZERO
   10          CONTINUE
            ELSE
               DO 20, I = 1, LENY
                  Y( I ) = BETA*Y( I )
   20          CONTINUE
            END IF
         ELSE
            IY = KY
            IF( BETA.EQ.ZERO )THEN
               DO 30, I = 1, LENY
                  Y( IY ) = ZERO
                  IY      = IY   + INCY
   30          CONTINUE
            ELSE
               DO 40, I = 1, LENY
                  Y( IY ) = BETA*Y( IY )
                  IY      = IY           + INCY
   40          CONTINUE
            END IF
         END IF
      END IF
      IF( ALPHA.EQ.ZERO )
     $   RETURN
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  y := alpha*A*x + y.
*
         JX = KX
         IF( INCY.EQ.1 )THEN
            DO 60, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  DO 50, I = 1, M
                     Y( I ) = Y( I ) + TEMP*A( I, J )
   50             CONTINUE
               END IF
               JX = JX + INCX
   60       CONTINUE
         ELSE
            DO 80, J = 1, N
               IF( X( JX ).NE.ZERO )THEN
                  TEMP = ALPHA*X( JX )
                  IY   = KY
                  DO 70, I = 1, M
                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
                     IY      = IY      + INCY
   70             CONTINUE
               END IF
               JX = JX + INCX
   80       CONTINUE
         END IF
      ELSE
*
*        Form  y := alpha*A'*x + y.
*
         JY = KY
         IF( INCX.EQ.1 )THEN
            DO 100, J = 1, N
               TEMP = ZERO
               DO 90, I = 1, M
                  TEMP = TEMP + A( I, J )*X( I )
   90          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  100       CONTINUE
         ELSE
            DO 120, J = 1, N
               TEMP = ZERO
               IX   = KX
               DO 110, I = 1, M
                  TEMP = TEMP + A( I, J )*X( IX )
                  IX   = IX   + INCX
  110          CONTINUE
               Y( JY ) = Y( JY ) + ALPHA*TEMP
               JY      = JY      + INCY
  120       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DGEMV .
*
      END
      SUBROUTINE DGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
*     .. Scalar Arguments ..
      DOUBLE PRECISION   ALPHA
      INTEGER            INCX, INCY, LDA, M, N
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  DGER   performs the rank 1 operation
*
*     A := alpha*x*y' + A,
*
*  where alpha is a scalar, x is an m element vector, y is an n element
*  vector and A is an m by n matrix.
*
*  Parameters
*  ==========
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( m - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the m
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y.
*           Unchanged on exit.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients. On exit, A is
*           overwritten by the updated matrix.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JY, KX
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( M.LT.0 )THEN
         INFO = 1
      ELSE IF( N.LT.0 )THEN
         INFO = 2
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 5
      ELSE IF( INCY.EQ.0 )THEN
         INFO = 7
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGER  ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
     $   RETURN
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF( INCY.GT.0 )THEN
         JY = 1
      ELSE
         JY = 1 - ( N - 1 )*INCY
      END IF
      IF( INCX.EQ.1 )THEN
         DO 20, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               DO 10, I = 1, M
                  A( I, J ) = A( I, J ) + X( I )*TEMP
   10          CONTINUE
            END IF
            JY = JY + INCY
   20    CONTINUE
      ELSE
         IF( INCX.GT.0 )THEN
            KX = 1
         ELSE
            KX = 1 - ( M - 1 )*INCX
         END IF
         DO 40, J = 1, N
            IF( Y( JY ).NE.ZERO )THEN
               TEMP = ALPHA*Y( JY )
               IX   = KX
               DO 30, I = 1, M
                  A( I, J ) = A( I, J ) + X( IX )*TEMP
                  IX        = IX        + INCX
   30          CONTINUE
            END IF
            JY = JY + INCY
   40    CONTINUE
      END IF
*
      RETURN
*
*     End of DGER  .
*
      END
      SUBROUTINE DGETF2( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETF2 computes an LU factorization of a general m-by-n matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 2 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the m by n matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
*               has been completed, but the factor U is exactly
*               singular, and division by zero will occur if it is used
*               to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J, JP
*     ..
*     .. External Functions ..
      INTEGER            IDAMAX
      EXTERNAL           IDAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGER, DSCAL, DSWAP, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
      DO 10 J = 1, MIN( M, N )
*
*        Find pivot and test for singularity.
*
         JP = J - 1 + IDAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
*
*           Apply the interchange to columns 1:N.
*
            IF( JP.NE.J )
     $         CALL DSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
*
*           Compute elements J+1:M of J-th column.
*
            IF( J.LT.M )
     $         CALL DSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
*
         ELSE IF( INFO.EQ.0 ) THEN
*
            INFO = J
         END IF
*
         IF( J.LT.MIN( M, N ) ) THEN
*
*           Update trailing submatrix.
*
            CALL DGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA,
     $                 A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
*
*     End of DGETF2
*
      END
      SUBROUTINE DGETRF( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETRF computes an LU factorization of a general M-by-N matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 3 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*                has been completed, but the factor U is exactly
*                singular, and division by zero will occur if it is used
*                to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DGETF2, DLASWP, DTRSM, XERBLA
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'DGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
*
*        Use unblocked code.
*
         CALL DGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
*
*        Use blocked code.
*
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
*
*           Factor diagonal and subdiagonal blocks and test for exact
*           singularity.
*
            CALL DGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
*
*           Adjust INFO and the pivot indices.
*
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
*
*           Apply interchanges to columns 1:J-1.
*
            CALL DLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
*
            IF( J+JB.LE.N ) THEN
*
*              Apply interchanges to columns J+JB:N.
*
               CALL DLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,
     $                      IPIV, 1 )
*
*              Compute block row of U.
*
               CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,
     $                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),
     $                     LDA )
               IF( J+JB.LE.M ) THEN
*
*                 Update trailing submatrix.
*
                  CALL DGEMM( 'No transpose', 'No transpose', M-J-JB+1,
     $                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
     $                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
     $                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
*
*     End of DGETRF
*
      END
      SUBROUTINE DGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), WORK( LWORK )
*     ..
*
*  Purpose
*  =======
*
*  DGETRI computes the inverse of a matrix using the LU factorization
*  computed by DGETRF.
*
*  This method inverts U and then computes inv(A) by solving the system
*  inv(A)*L = inv(U) for inv(A).
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the factors L and U from the factorization
*          A = P*L*U as computed by DGETRF.
*          On exit, if INFO = 0, the inverse of the original matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices from DGETRF; for 1<=i<=N, row i of the
*          matrix was interchanged with row IPIV(i).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK)
*          On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,N).
*          For optimal performance LWORK >= N*NB, where NB is
*          the optimal blocksize returned by ILAENV.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
*                singular and its inverse could not be computed.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IWS, J, JB, JJ, JP, LDWORK, NB, NBMIN, NN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEMM, DGEMV, DSWAP, DTRSM, DTRTRI, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      WORK( 1 ) = MAX( N, 1 )
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -3
      ELSE IF( LWORK.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRI', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Form inv(U).  If INFO > 0 from DTRTRI, then U is singular,
*     and the inverse is not computed.
*
      CALL DTRTRI( 'Upper', 'Non-unit', N, A, LDA, INFO )
      IF( INFO.GT.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'DGETRI', ' ', N, -1, -1, -1 )
      NBMIN = 2
      LDWORK = N
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
         IWS = MAX( LDWORK*NB, 1 )
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'DGETRI', ' ', N, -1, -1, -1 ) )
         END IF
      ELSE
         IWS = N
      END IF
*
*     Solve the equation inv(A)*L = inv(U) for inv(A).
*
      IF( NB.LT.NBMIN .OR. NB.GE.N ) THEN
*
*        Use unblocked code.
*
         DO 20 J = N, 1, -1
*
*           Copy current column of L to WORK and replace with zeros.
*
            DO 10 I = J + 1, N
               WORK( I ) = A( I, J )
               A( I, J ) = ZERO
   10       CONTINUE
*
*           Compute current column of inv(A).
*
            IF( J.LT.N )
     $         CALL DGEMV( 'No transpose', N, N-J, -ONE, A( 1, J+1 ),
     $                     LDA, WORK( J+1 ), 1, ONE, A( 1, J ), 1 )
   20    CONTINUE
      ELSE
*
*        Use blocked code.
*
         NN = ( ( N-1 ) / NB )*NB + 1
         DO 50 J = NN, 1, -NB
            JB = MIN( NB, N-J+1 )
*
*           Copy current block column of L to WORK and replace with
*           zeros.
*
            DO 40 JJ = J, J + JB - 1
               DO 30 I = JJ + 1, N
                  WORK( I+( JJ-J )*LDWORK ) = A( I, JJ )
                  A( I, JJ ) = ZERO
   30          CONTINUE
   40       CONTINUE
*
*           Compute current block column of inv(A).
*
            IF( J+JB.LE.N )
     $         CALL DGEMM( 'No transpose', 'No transpose', N, JB,
     $                     N-J-JB+1, -ONE, A( 1, J+JB ), LDA,
     $                     WORK( J+JB ), LDWORK, ONE, A( 1, J ), LDA )
            CALL DTRSM( 'Right', 'Lower', 'No transpose', 'Unit', N, JB,
     $                  ONE, WORK( J ), LDWORK, A( 1, J ), LDA )
   50    CONTINUE
      END IF
*
*     Apply column interchanges.
*
      DO 60 J = N - 1, 1, -1
         JP = IPIV( J )
         IF( JP.NE.J )
     $      CALL DSWAP( N, A( 1, J ), 1, A( 1, JP ), 1 )
   60 CONTINUE
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of DGETRI
*
      END
      SUBROUTINE DGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          TRANS
      INTEGER            INFO, LDA, LDB, N, NRHS
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DGETRS solves a system of linear equations
*     A * X = B  or  A' * X = B
*  with a general N-by-N matrix A using the LU factorization computed
*  by DGETRF.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations:
*          = 'N':  A * X = B  (No transpose)
*          = 'T':  A'* X = B  (Transpose)
*          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The factors L and U from the factorization A = P*L*U
*          as computed by DGETRF.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices from DGETRF; for 1<=i<=N, row i of the
*          matrix was interchanged with row IPIV(i).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOTRAN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASWP, DTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NOTRAN = LSAME( TRANS, 'N' )
      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
     $    LSAME( TRANS, 'C' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NRHS.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DGETRS', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. NRHS.EQ.0 )
     $   RETURN
*
      IF( NOTRAN ) THEN
*
*        Solve A * X = B.
*
*        Apply row interchanges to the right hand sides.
*
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
*
*        Solve L*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
*
*        Solve U*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
     $               NRHS, ONE, A, LDA, B, LDB )
      ELSE
*
*        Solve A' * X = B.
*
*        Solve U'*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS,
     $               ONE, A, LDA, B, LDB )
*
*        Solve L'*X = B, overwriting B with X.
*
         CALL DTRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE,
     $               A, LDA, B, LDB )
*
*        Apply row interchanges to the solution vectors.
*
         CALL DLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
      END IF
*
      RETURN
*
*     End of DGETRS
*
      END
      SUBROUTINE DLASWP( N, A, LDA, K1, K2, IPIV, INCX )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DLASWP performs a series of row interchanges on the matrix A.
*  One row interchange is initiated for each of rows K1 through K2 of A.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the matrix of column dimension N to which the row
*          interchanges will be applied.
*          On exit, the permuted matrix.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*
*  K1      (input) INTEGER
*          The first element of IPIV for which a row interchange will
*          be done.
*
*  K2      (input) INTEGER
*          The last element of IPIV for which a row interchange will
*          be done.
*
*  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
*          The vector of pivot indices.  Only the elements in positions
*          K1 through K2 of IPIV are accessed.
*          IPIV(K) = L implies rows K and L are to be interchanged.
*
*  INCX    (input) INTEGER
*          The increment between successive values of IPIV.  If IPIV
*          is negative, the pivots are applied in reverse order.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IP, IX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSWAP
*     ..
*     .. Executable Statements ..
*
*     Interchange row I with row IPIV(I) for each of rows K1 through K2.
*
      IF( INCX.EQ.0 )
     $   RETURN
      IF( INCX.GT.0 ) THEN
         IX = K1
      ELSE
         IX = 1 + ( 1-K2 )*INCX
      END IF
      IF( INCX.EQ.1 ) THEN
         DO 10 I = K1, K2
            IP = IPIV( I )
            IF( IP.NE.I )
     $         CALL DSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
   10    CONTINUE
      ELSE IF( INCX.GT.1 ) THEN
         DO 20 I = K1, K2
            IP = IPIV( IX )
            IF( IP.NE.I )
     $         CALL DSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   20    CONTINUE
      ELSE IF( INCX.LT.0 ) THEN
         DO 30 I = K2, K1, -1
            IP = IPIV( IX )
            IF( IP.NE.I )
     $         CALL DSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
            IX = IX + INCX
   30    CONTINUE
      END IF
*
      RETURN
*
*     End of DLASWP
*
      END
      SUBROUTINE DTRMM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   B, LDB )
*     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      DOUBLE PRECISION   ALPHA
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DTRMM  performs one of the matrix-matrix operations
*
*     B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
*
*  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*
*     op( A ) = A   or   op( A ) = A'.
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry,  SIDE specifies whether  op( A ) multiplies B from
*           the left or right as follows:
*
*              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
*
*              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix A is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n'   op( A ) = A.
*
*              TRANSA = 'T' or 't'   op( A ) = A'.
*
*              TRANSA = 'C' or 'c'   op( A ) = A'.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit triangular
*           as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of B. M must be at
*           least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of B.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*           zero then  A is not referenced and  B need not be set before
*           entry.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*           upper triangular part of the array  A must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           A is not referenced.
*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*           lower triangular part of the array  A must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*           A  are not referenced either,  but are assumed to be  unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*           then LDA must be at least max( 1, n ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
*           Before entry,  the leading  m by n part of the array  B must
*           contain the matrix  B,  and  on exit  is overwritten  by the
*           transformed matrix.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            LSIDE, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      DOUBLE PRECISION   TEMP
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
*
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.
     $         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.
     $         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
     $         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRMM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
*     And when  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
*
*     Start the operations.
*
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*A*B.
*
            IF( UPPER )THEN
               DO 50, J = 1, N
                  DO 40, K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*B( K, J )
                        DO 30, I = 1, K - 1
                           B( I, J ) = B( I, J ) + TEMP*A( I, K )
   30                   CONTINUE
                        IF( NOUNIT )
     $                     TEMP = TEMP*A( K, K )
                        B( K, J ) = TEMP
                     END IF
   40             CONTINUE
   50          CONTINUE
            ELSE
               DO 80, J = 1, N
                  DO 70 K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        TEMP      = ALPHA*B( K, J )
                        B( K, J ) = TEMP
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )*A( K, K )
                        DO 60, I = K + 1, M
                           B( I, J ) = B( I, J ) + TEMP*A( I, K )
   60                   CONTINUE
                     END IF
   70             CONTINUE
   80          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*B*A'.
*
            IF( UPPER )THEN
               DO 110, J = 1, N
                  DO 100, I = M, 1, -1
                     TEMP = B( I, J )
                     IF( NOUNIT )
     $                  TEMP = TEMP*A( I, I )
                     DO 90, K = 1, I - 1
                        TEMP = TEMP + A( K, I )*B( K, J )
   90                CONTINUE
                     B( I, J ) = ALPHA*TEMP
  100             CONTINUE
  110          CONTINUE
            ELSE
               DO 140, J = 1, N
                  DO 130, I = 1, M
                     TEMP = B( I, J )
                     IF( NOUNIT )
     $                  TEMP = TEMP*A( I, I )
                     DO 120, K = I + 1, M
                        TEMP = TEMP + A( K, I )*B( K, J )
  120                CONTINUE
                     B( I, J ) = ALPHA*TEMP
  130             CONTINUE
  140          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*B*A.
*
            IF( UPPER )THEN
               DO 180, J = N, 1, -1
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 150, I = 1, M
                     B( I, J ) = TEMP*B( I, J )
  150             CONTINUE
                  DO 170, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*A( K, J )
                        DO 160, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  160                   CONTINUE
                     END IF
  170             CONTINUE
  180          CONTINUE
            ELSE
               DO 220, J = 1, N
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 190, I = 1, M
                     B( I, J ) = TEMP*B( I, J )
  190             CONTINUE
                  DO 210, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*A( K, J )
                        DO 200, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  200                   CONTINUE
                     END IF
  210             CONTINUE
  220          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*B*A'.
*
            IF( UPPER )THEN
               DO 260, K = 1, N
                  DO 240, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = ALPHA*A( J, K )
                        DO 230, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  230                   CONTINUE
                     END IF
  240             CONTINUE
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( K, K )
                  IF( TEMP.NE.ONE )THEN
                     DO 250, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  250                CONTINUE
                  END IF
  260          CONTINUE
            ELSE
               DO 300, K = N, 1, -1
                  DO 280, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = ALPHA*A( J, K )
                        DO 270, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  270                   CONTINUE
                     END IF
  280             CONTINUE
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( K, K )
                  IF( TEMP.NE.ONE )THEN
                     DO 290, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  290                CONTINUE
                  END IF
  300          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DTRMM .
*
      END
      SUBROUTINE DTRMV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
*     .. Scalar Arguments ..
      INTEGER            INCX, LDA, N
      CHARACTER*1        DIAG, TRANS, UPLO
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), X( * )
*     ..
*
*  Purpose
*  =======
*
*  DTRMV  performs one of the matrix-vector operations
*
*     x := A*x,   or   x := A'*x,
*
*  where x is an n element vector and  A is an n by n unit, or non-unit,
*  upper or lower triangular matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   x := A*x.
*
*              TRANS = 'T' or 't'   x := A'*x.
*
*              TRANS = 'C' or 'c'   x := A'*x.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U' or 'u', the leading n by n
*           upper triangular part of the array A must contain the upper
*           triangular matrix and the strictly lower triangular part of
*           A is not referenced.
*           Before entry with UPLO = 'L' or 'l', the leading n by n
*           lower triangular part of the array A must contain the lower
*           triangular matrix and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u', the diagonal elements of
*           A are not referenced either, but are assumed to be unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, n ).
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x. On exit, X is overwritten with the
*           tranformed vector x.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*
*  Level 2 Blas routine.
*
*  -- Written on 22-October-1986.
*     Jack Dongarra, Argonne National Lab.
*     Jeremy Du Croz, Nag Central Office.
*     Sven Hammarling, Nag Central Office.
*     Richard Hanson, Sandia National Labs.
*
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER        ( ZERO = 0.0D+0 )
*     .. Local Scalars ..
      DOUBLE PRECISION   TEMP
      INTEGER            I, INFO, IX, J, JX, KX
      LOGICAL            NOUNIT
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF     ( .NOT.LSAME( UPLO , 'U' ).AND.
     $         .NOT.LSAME( UPLO , 'L' )      )THEN
         INFO = 1
      ELSE IF( .NOT.LSAME( TRANS, 'N' ).AND.
     $         .NOT.LSAME( TRANS, 'T' ).AND.
     $         .NOT.LSAME( TRANS, 'C' )      )THEN
         INFO = 2
      ELSE IF( .NOT.LSAME( DIAG , 'U' ).AND.
     $         .NOT.LSAME( DIAG , 'N' )      )THEN
         INFO = 3
      ELSE IF( N.LT.0 )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
         INFO = 6
      ELSE IF( INCX.EQ.0 )THEN
         INFO = 8
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRMV ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
      NOUNIT = LSAME( DIAG, 'N' )
*
*     Set up the start point in X if the increment is not unity. This
*     will be  ( N - 1 )*INCX  too small for descending loops.
*
      IF( INCX.LE.0 )THEN
         KX = 1 - ( N - 1 )*INCX
      ELSE IF( INCX.NE.1 )THEN
         KX = 1
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
      IF( LSAME( TRANS, 'N' ) )THEN
*
*        Form  x := A*x.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 20, J = 1, N
                  IF( X( J ).NE.ZERO )THEN
                     TEMP = X( J )
                     DO 10, I = 1, J - 1
                        X( I ) = X( I ) + TEMP*A( I, J )
   10                CONTINUE
                     IF( NOUNIT )
     $                  X( J ) = X( J )*A( J, J )
                  END IF
   20          CONTINUE
            ELSE
               JX = KX
               DO 40, J = 1, N
                  IF( X( JX ).NE.ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX
                     DO 30, I = 1, J - 1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      + INCX
   30                CONTINUE
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX + INCX
   40          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 60, J = N, 1, -1
                  IF( X( J ).NE.ZERO )THEN
                     TEMP = X( J )
                     DO 50, I = N, J + 1, -1
                        X( I ) = X( I ) + TEMP*A( I, J )
   50                CONTINUE
                     IF( NOUNIT )
     $                  X( J ) = X( J )*A( J, J )
                  END IF
   60          CONTINUE
            ELSE
               KX = KX + ( N - 1 )*INCX
               JX = KX
               DO 80, J = N, 1, -1
                  IF( X( JX ).NE.ZERO )THEN
                     TEMP = X( JX )
                     IX   = KX
                     DO 70, I = N, J + 1, -1
                        X( IX ) = X( IX ) + TEMP*A( I, J )
                        IX      = IX      - INCX
   70                CONTINUE
                     IF( NOUNIT )
     $                  X( JX ) = X( JX )*A( J, J )
                  END IF
                  JX = JX - INCX
   80          CONTINUE
            END IF
         END IF
      ELSE
*
*        Form  x := A'*x.
*
         IF( LSAME( UPLO, 'U' ) )THEN
            IF( INCX.EQ.1 )THEN
               DO 100, J = N, 1, -1
                  TEMP = X( J )
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 90, I = J - 1, 1, -1
                     TEMP = TEMP + A( I, J )*X( I )
   90             CONTINUE
                  X( J ) = TEMP
  100          CONTINUE
            ELSE
               JX = KX + ( N - 1 )*INCX
               DO 120, J = N, 1, -1
                  TEMP = X( JX )
                  IX   = JX
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 110, I = J - 1, 1, -1
                     IX   = IX   - INCX
                     TEMP = TEMP + A( I, J )*X( IX )
  110             CONTINUE
                  X( JX ) = TEMP
                  JX      = JX   - INCX
  120          CONTINUE
            END IF
         ELSE
            IF( INCX.EQ.1 )THEN
               DO 140, J = 1, N
                  TEMP = X( J )
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 130, I = J + 1, N
                     TEMP = TEMP + A( I, J )*X( I )
  130             CONTINUE
                  X( J ) = TEMP
  140          CONTINUE
            ELSE
               JX = KX
               DO 160, J = 1, N
                  TEMP = X( JX )
                  IX   = JX
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 150, I = J + 1, N
                     IX   = IX   + INCX
                     TEMP = TEMP + A( I, J )*X( IX )
  150             CONTINUE
                  X( JX ) = TEMP
                  JX      = JX   + INCX
  160          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DTRMV .
*
      END
      SUBROUTINE DTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   B, LDB )
*     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      DOUBLE PRECISION   ALPHA
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DTRSM  solves one of the matrix equations
*
*     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
*
*  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*
*     op( A ) = A   or   op( A ) = A'.
*
*  The matrix X is overwritten on B.
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry, SIDE specifies whether op( A ) appears on the left
*           or right of X as follows:
*
*              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
*
*              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix A is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n'   op( A ) = A.
*
*              TRANSA = 'T' or 't'   op( A ) = A'.
*
*              TRANSA = 'C' or 'c'   op( A ) = A'.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit triangular
*           as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of B. M must be at
*           least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of B.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*           zero then  A is not referenced and  B need not be set before
*           entry.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*           upper triangular part of the array  A must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           A is not referenced.
*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*           lower triangular part of the array  A must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*           A  are not referenced either,  but are assumed to be  unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*           then LDA must be at least max( 1, n ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
*           Before entry,  the leading  m by n part of the array  B must
*           contain  the  right-hand  side  matrix  B,  and  on exit  is
*           overwritten by the solution matrix  X.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*
*  Level 3 Blas routine.
*
*
*  -- Written on 8-February-1989.
*     Jack Dongarra, Argonne National Laboratory.
*     Iain Duff, AERE Harwell.
*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
*     Sven Hammarling, Numerical Algorithms Group Ltd.
*
*
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     .. Local Scalars ..
      LOGICAL            LSIDE, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      DOUBLE PRECISION   TEMP
*     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
*
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.
     $         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.
     $         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
     $         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRSM ', INFO )
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
*     And when  alpha.eq.zero.
*
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
*
*     Start the operations.
*
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*inv( A )*B.
*
            IF( UPPER )THEN
               DO 60, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 30, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   30                CONTINUE
                  END IF
                  DO 50, K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 40, I = 1, K - 1
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   40                   CONTINUE
                     END IF
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 100, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 70, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   70                CONTINUE
                  END IF
                  DO 90 K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 80, I = K + 1, M
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   80                   CONTINUE
                     END IF
   90             CONTINUE
  100          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*inv( A' )*B.
*
            IF( UPPER )THEN
               DO 130, J = 1, N
                  DO 120, I = 1, M
                     TEMP = ALPHA*B( I, J )
                     DO 110, K = 1, I - 1
                        TEMP = TEMP - A( K, I )*B( K, J )
  110                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  120             CONTINUE
  130          CONTINUE
            ELSE
               DO 160, J = 1, N
                  DO 150, I = M, 1, -1
                     TEMP = ALPHA*B( I, J )
                     DO 140, K = I + 1, M
                        TEMP = TEMP - A( K, I )*B( K, J )
  140                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  150             CONTINUE
  160          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
*
*           Form  B := alpha*B*inv( A ).
*
            IF( UPPER )THEN
               DO 210, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 170, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  170                CONTINUE
                  END IF
                  DO 190, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 180, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  180                   CONTINUE
                     END IF
  190             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 200, I = 1, M
                        B( I, J ) = TEMP*B( I, J )
  200                CONTINUE
                  END IF
  210          CONTINUE
            ELSE
               DO 260, J = N, 1, -1
                  IF( ALPHA.NE.ONE )THEN
                     DO 220, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  220                CONTINUE
                  END IF
                  DO 240, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 230, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  230                   CONTINUE
                     END IF
  240             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 250, I = 1, M
                       B( I, J ) = TEMP*B( I, J )
  250                CONTINUE
                  END IF
  260          CONTINUE
            END IF
         ELSE
*
*           Form  B := alpha*B*inv( A' ).
*
            IF( UPPER )THEN
               DO 310, K = N, 1, -1
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 270, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  270                CONTINUE
                  END IF
                  DO 290, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 280, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  280                   CONTINUE
                     END IF
  290             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 300, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  300                CONTINUE
                  END IF
  310          CONTINUE
            ELSE
               DO 360, K = 1, N
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 320, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  320                CONTINUE
                  END IF
                  DO 340, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 330, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  330                   CONTINUE
                     END IF
  340             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 350, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  350                CONTINUE
                  END IF
  360          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of DTRSM .
*
      END
      SUBROUTINE DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DTRTI2 computes the inverse of a real upper or lower triangular
*  matrix.
*
*  This is the Level 2 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the matrix A is upper or lower triangular.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  DIAG    (input) CHARACTER*1
*          Specifies whether or not the matrix A is unit triangular.
*          = 'N':  Non-unit triangular
*          = 'U':  Unit triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the triangular matrix A.  If UPLO = 'U', the
*          leading n by n upper triangular part of the array A contains
*          the upper triangular matrix, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading n by n lower triangular part of the array A contains
*          the lower triangular matrix, and the strictly upper
*          triangular part of A is not referenced.  If DIAG = 'U', the
*          diagonal elements of A are also not referenced and are
*          assumed to be 1.
*
*          On exit, the (triangular) inverse of the original matrix, in
*          the same storage format.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT, UPPER
      INTEGER            J
      DOUBLE PRECISION   AJJ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL, DTRMV, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTRTI2', -INFO )
         RETURN
      END IF
*
      IF( UPPER ) THEN
*
*        Compute inverse of upper triangular matrix.
*
         DO 10 J = 1, N
            IF( NOUNIT ) THEN
               A( J, J ) = ONE / A( J, J )
               AJJ = -A( J, J )
            ELSE
               AJJ = -ONE
            END IF
*
*           Compute elements 1:j-1 of j-th column.
*
            CALL DTRMV( 'Upper', 'No transpose', DIAG, J-1, A, LDA,
     $                  A( 1, J ), 1 )
            CALL DSCAL( J-1, AJJ, A( 1, J ), 1 )
   10    CONTINUE
      ELSE
*
*        Compute inverse of lower triangular matrix.
*
         DO 20 J = N, 1, -1
            IF( NOUNIT ) THEN
               A( J, J ) = ONE / A( J, J )
               AJJ = -A( J, J )
            ELSE
               AJJ = -ONE
            END IF
            IF( J.LT.N ) THEN
*
*              Compute elements j+1:n of j-th column.
*
               CALL DTRMV( 'Lower', 'No transpose', DIAG, N-J,
     $                     A( J+1, J+1 ), LDA, A( J+1, J ), 1 )
               CALL DSCAL( N-J, AJJ, A( J+1, J ), 1 )
            END IF
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of DTRTI2
*
      END
      SUBROUTINE DTRTRI( UPLO, DIAG, N, A, LDA, INFO )
*
*  -- LAPACK routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     March 31, 1993
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  DTRTRI computes the inverse of a real upper or lower triangular
*  matrix A.
*
*  This is the Level 3 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  A is upper triangular;
*          = 'L':  A is lower triangular.
*
*  DIAG    (input) CHARACTER*1
*          = 'N':  A is non-unit triangular;
*          = 'U':  A is unit triangular.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the triangular matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of the array A contains
*          the upper triangular matrix, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of the array A contains
*          the lower triangular matrix, and the strictly upper
*          triangular part of A is not referenced.  If DIAG = 'U', the
*          diagonal elements of A are also not referenced and are
*          assumed to be 1.
*          On exit, the (triangular) inverse of the original matrix, in
*          the same storage format.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
*               matrix is singular and its inverse can not be computed.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT, UPPER
      INTEGER            J, JB, NB, NN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DTRMM, DTRSM, DTRTI2, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTRTRI', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Check for singularity if non-unit.
*
      IF( NOUNIT ) THEN
         DO 10 INFO = 1, N
            IF( A( INFO, INFO ).EQ.ZERO )
     $         RETURN
   10    CONTINUE
         INFO = 0
      END IF
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'DTRTRI', UPLO // DIAG, N, -1, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
*
*        Use unblocked code
*
         CALL DTRTI2( UPLO, DIAG, N, A, LDA, INFO )
      ELSE
*
*        Use blocked code
*
         IF( UPPER ) THEN
*
*           Compute inverse of upper triangular matrix
*
            DO 20 J = 1, N, NB
               JB = MIN( NB, N-J+1 )
*
*              Compute rows 1:j-1 of current block column
*
               CALL DTRMM( 'Left', 'Upper', 'No transpose', DIAG, J-1,
     $                     JB, ONE, A, LDA, A( 1, J ), LDA )
               CALL DTRSM( 'Right', 'Upper', 'No transpose', DIAG, J-1,
     $                     JB, -ONE, A( J, J ), LDA, A( 1, J ), LDA )
*
*              Compute inverse of current diagonal block
*
               CALL DTRTI2( 'Upper', DIAG, JB, A( J, J ), LDA, INFO )
   20       CONTINUE
         ELSE
*
*           Compute inverse of lower triangular matrix
*
            NN = ( ( N-1 ) / NB )*NB + 1
            DO 30 J = NN, 1, -NB
               JB = MIN( NB, N-J+1 )
               IF( J+JB.LE.N ) THEN
*
*                 Compute rows j+jb:n of current block column
*
                  CALL DTRMM( 'Left', 'Lower', 'No transpose', DIAG,
     $                        N-J-JB+1, JB, ONE, A( J+JB, J+JB ), LDA,
     $                        A( J+JB, J ), LDA )
                  CALL DTRSM( 'Right', 'Lower', 'No transpose', DIAG,
     $                        N-J-JB+1, JB, -ONE, A( J, J ), LDA,
     $                        A( J+JB, J ), LDA )
               END IF
*
*              Compute inverse of current diagonal block
*
               CALL DTRTI2( 'Lower', DIAG, JB, A( J, J ), LDA, INFO )
   30       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of DTRTRI
*
      END
      double precision function dasum(n,dx,incx)
c
c     takes the sum of the absolute values.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      double precision dx(1),dtemp
      integer i,incx,m,mp1,n,nincx
c
      dasum = 0.0d0
      dtemp = 0.0d0
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dtemp = dtemp + dabs(dx(i))
   10 continue
      dasum = dtemp
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dabs(dx(i))
   30 continue
      if( n .lt. 6 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,6
        dtemp = dtemp + dabs(dx(i)) + dabs(dx(i + 1)) + dabs(dx(i + 2))
     *  + dabs(dx(i + 3)) + dabs(dx(i + 4)) + dabs(dx(i + 5))
   50 continue
   60 dasum = dtemp
      return
      end
      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end
      subroutine  dcopy(n,dx,incx,dy,incy)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
        dy(i) = dx(i)
        dy(i + 1) = dx(i + 1)
        dy(i + 2) = dx(i + 2)
        dy(i + 3) = dx(i + 3)
        dy(i + 4) = dx(i + 4)
        dy(i + 5) = dx(i + 5)
        dy(i + 6) = dx(i + 6)
   50 continue
      return
      end
      double precision function ddot(n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end
      subroutine dgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
      integer lda,n,ml,mu,ipvt(1)
      double precision abd(lda,1),z(1)
      double precision rcond
c
c     dgbco factors a double precision band matrix by gaussian
c     elimination and estimates the condition of the matrix.
c
c     if  rcond  is not needed, dgbfa is slightly faster.
c     to solve  a*x = b , follow dgbco by dgbsl.
c     to compute  inverse(a)*c , follow dgbco by dgbsl.
c     to compute  determinant(a) , follow dgbco by dgbdi.
c
c     on entry
c
c        abd     double precision(lda, n)
c                contains the matrix in band storage.  the columns
c                of the matrix are stored in the columns of  abd  and
c                the diagonals of the matrix are stored in rows
c                ml+1 through 2*ml+mu+1 of  abd .
c                see the comments below for details.
c
c        lda     integer
c                the leading dimension of the array  abd .
c                lda must be .ge. 2*ml + mu + 1 .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c                0 .le. ml .lt. n .
c
c        mu      integer
c                number of diagonals above the main diagonal.
c                0 .le. mu .lt. n .
c                more efficient if  ml .le. mu .
c
c     on return
c
c        abd     an upper triangular matrix in band storage and
c                the multipliers which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        rcond   double precision
c                an estimate of the reciprocal condition of  a .
c                for the system  a*x = b , relative perturbations
c                in  a  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  a  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.
c
c        z       double precision(n)
c                a work vector whose contents are usually unimportant.
c                if  a  is close to a singular matrix, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c
c     band storage
c
c           if  a  is a band matrix, the following program segment
c           will set up the input.
c
c                   ml = (band width below the diagonal)
c                   mu = (band width above the diagonal)
c                   m = ml + mu + 1
c                   do 20 j = 1, n
c                      i1 = max0(1, j-mu)
c                      i2 = min0(n, j+ml)
c                      do 10 i = i1, i2
c                         k = i - j + m
c                         abd(k,j) = a(i,j)
c                10    continue
c                20 continue
c
c           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
c           in addition, the first  ml  rows in  abd  are used for
c           elements generated during the triangularization.
c           the total number of rows needed in  abd  is  2*ml+mu+1 .
c           the  ml+mu by ml+mu  upper left triangle and the
c           ml by ml  lower right triangle are not referenced.
c
c     example..  if the original matrix is
c
c           11 12 13  0  0  0
c           21 22 23 24  0  0
c            0 32 33 34 35  0
c            0  0 43 44 45 46
c            0  0  0 54 55 56
c            0  0  0  0 65 66
c
c      then  n = 6, ml = 1, mu = 2, lda .ge. 5  and abd should contain
c
c            *  *  *  +  +  +  , * = not used
c            *  * 13 24 35 46  , + = used for pivoting
c            * 12 23 34 45 56
c           11 22 33 44 55 66
c           21 32 43 54 65  *
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     linpack dgbfa
c     blas daxpy,ddot,dscal,dasum
c     fortran dabs,dmax1,max0,min0,dsign
c
c     internal variables
c
      double precision ddot,ek,t,wk,wkm
      double precision anorm,s,dasum,sm,ynorm
      integer is,info,j,ju,k,kb,kp1,l,la,lm,lz,m,mm
c
c
c     compute 1-norm of a
c
      anorm = 0.0d0
      l = ml + 1
      is = l + mu
      do 10 j = 1, n
         anorm = dmax1(anorm,dasum(l,abd(is,j),1))
         if (is .gt. ml + 1) is = is - 1
         if (j .le. mu) l = l + 1
         if (j .ge. n - ml) l = l - 1
   10 continue
c
c     factor
c
      call dgbfa(abd,lda,n,ml,mu,ipvt,info)
c
c     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c     estimate = norm(z)/norm(y) where  a*z = y  and  trans(a)*y = e .
c     trans(a)  is the transpose of a .  the components of  e  are
c     chosen to cause maximum local growth in the elements of w  where
c     trans(u)*w = e .  the vectors are frequently rescaled to avoid
c     overflow.
c
c     solve trans(u)*w = e
c
      ek = 1.0d0
      do 20 j = 1, n
         z(j) = 0.0d0
   20 continue
      m = ml + mu + 1
      ju = 0
      do 100 k = 1, n
         if (z(k) .ne. 0.0d0) ek = dsign(ek,-z(k))
         if (dabs(ek-z(k)) .le. dabs(abd(m,k))) go to 30
            s = dabs(abd(m,k))/dabs(ek-z(k))
            call dscal(n,s,z,1)
            ek = s*ek
   30    continue
         wk = ek - z(k)
         wkm = -ek - z(k)
         s = dabs(wk)
         sm = dabs(wkm)
         if (abd(m,k) .eq. 0.0d0) go to 40
            wk = wk/abd(m,k)
            wkm = wkm/abd(m,k)
         go to 50
   40    continue
            wk = 1.0d0
            wkm = 1.0d0
   50    continue
         kp1 = k + 1
         ju = min0(max0(ju,mu+ipvt(k)),n)
         mm = m
         if (kp1 .gt. ju) go to 90
            do 60 j = kp1, ju
               mm = mm - 1
               sm = sm + dabs(z(j)+wkm*abd(mm,j))
               z(j) = z(j) + wk*abd(mm,j)
               s = s + dabs(z(j))
   60       continue
            if (s .ge. sm) go to 80
               t = wkm - wk
               wk = wkm
               mm = m
               do 70 j = kp1, ju
                  mm = mm - 1
                  z(j) = z(j) + t*abd(mm,j)
   70          continue
   80       continue
   90    continue
         z(k) = wk
  100 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
c
c     solve trans(l)*y = w
c
      do 120 kb = 1, n
         k = n + 1 - kb
         lm = min0(ml,n-k)
         if (k .lt. n) z(k) = z(k) + ddot(lm,abd(m+1,k),1,z(k+1),1)
         if (dabs(z(k)) .le. 1.0d0) go to 110
            s = 1.0d0/dabs(z(k))
            call dscal(n,s,z,1)
  110    continue
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
  120 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
c
      ynorm = 1.0d0
c
c     solve l*v = y
c
      do 140 k = 1, n
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
         lm = min0(ml,n-k)
         if (k .lt. n) call daxpy(lm,t,abd(m+1,k),1,z(k+1),1)
         if (dabs(z(k)) .le. 1.0d0) go to 130
            s = 1.0d0/dabs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
  130    continue
  140 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
c
c     solve  u*z = w
c
      do 160 kb = 1, n
         k = n + 1 - kb
         if (dabs(z(k)) .le. dabs(abd(m,k))) go to 150
            s = dabs(abd(m,k))/dabs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
  150    continue
         if (abd(m,k) .ne. 0.0d0) z(k) = z(k)/abd(m,k)
         if (abd(m,k) .eq. 0.0d0) z(k) = 1.0d0
         lm = min0(k,m) - 1
         la = m - lm
         lz = k - lm
         t = -z(k)
         call daxpy(lm,t,abd(la,k),1,z(lz),1)
  160 continue
c     make znorm = 1.0
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
c
      if (anorm .ne. 0.0d0) rcond = ynorm/anorm
      if (anorm .eq. 0.0d0) rcond = 0.0d0
      return
      end
      subroutine dgbfa(abd,lda,n,ml,mu,ipvt,info)
      integer lda,n,ml,mu,ipvt(1),info
      double precision abd(lda,1)
c
c     dgbfa factors a double precision band matrix by elimination.
c
c     dgbfa is usually called by dgbco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c
c     on entry
c
c        abd     double precision(lda, n)
c                contains the matrix in band storage.  the columns
c                of the matrix are stored in the columns of  abd  and
c                the diagonals of the matrix are stored in rows
c                ml+1 through 2*ml+mu+1 of  abd .
c                see the comments below for details.
c
c        lda     integer
c                the leading dimension of the array  abd .
c                lda must be .ge. 2*ml + mu + 1 .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c                0 .le. ml .lt. n .
c
c        mu      integer
c                number of diagonals above the main diagonal.
c                0 .le. mu .lt. n .
c                more efficient if  ml .le. mu .
c     on return
c
c        abd     an upper triangular matrix in band storage and
c                the multipliers which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgbsl will divide by zero if
c                     called.  use  rcond  in dgbco for a reliable
c                     indication of singularity.
c
c     band storage
c
c           if  a  is a band matrix, the following program segment
c           will set up the input.
c
c                   ml = (band width below the diagonal)
c                   mu = (band width above the diagonal)
c                   m = ml + mu + 1
c                   do 20 j = 1, n
c                      i1 = max0(1, j-mu)
c                      i2 = min0(n, j+ml)
c                      do 10 i = i1, i2
c                         k = i - j + m
c                         abd(k,j) = a(i,j)
c                10    continue
c                20 continue
c
c           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
c           in addition, the first  ml  rows in  abd  are used for
c           elements generated during the triangularization.
c           the total number of rows needed in  abd  is  2*ml+mu+1 .
c           the  ml+mu by ml+mu  upper left triangle and the
c           ml by ml  lower right triangle are not referenced.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c     fortran max0,min0
c
c     internal variables
c
      double precision t
      integer i,idamax,i0,j,ju,jz,j0,j1,k,kp1,l,lm,m,mm,nm1
c
c
      m = ml + mu + 1
      info = 0
c
c     zero initial fill-in columns
c
      j0 = mu + 2
      j1 = min0(n,m) - 1
      if (j1 .lt. j0) go to 30
      do 20 jz = j0, j1
         i0 = m + 1 - jz
         do 10 i = i0, ml
            abd(i,jz) = 0.0d0
   10    continue
   20 continue
   30 continue
      jz = j1
      ju = 0
c
c     gaussian elimination with partial pivoting
c
      nm1 = n - 1
      if (nm1 .lt. 1) go to 130
      do 120 k = 1, nm1
         kp1 = k + 1
c
c        zero next fill-in column
c
         jz = jz + 1
         if (jz .gt. n) go to 50
         if (ml .lt. 1) go to 50
            do 40 i = 1, ml
               abd(i,jz) = 0.0d0
   40       continue
   50    continue
c
c        find l = pivot index
c
         lm = min0(ml,n-k)
         l = idamax(lm+1,abd(m,k),1) + m - 1
         ipvt(k) = l + k - m
c
c        zero pivot implies this column already triangularized
c
         if (abd(l,k) .eq. 0.0d0) go to 100
c
c           interchange if necessary
c
            if (l .eq. m) go to 60
               t = abd(l,k)
               abd(l,k) = abd(m,k)
               abd(m,k) = t
   60       continue
c
c           compute multipliers
c
            t = -1.0d0/abd(m,k)
            call dscal(lm,t,abd(m+1,k),1)
c
c           row elimination with column indexing
c
            ju = min0(max0(ju,mu+ipvt(k)),n)
            mm = m
            if (ju .lt. kp1) go to 90
            do 80 j = kp1, ju
               l = l - 1
               mm = mm - 1
               t = abd(l,j)
               if (l .eq. mm) go to 70
                  abd(l,j) = abd(mm,j)
                  abd(mm,j) = t
   70          continue
               call daxpy(lm,t,abd(m+1,k),1,abd(mm+1,j),1)
   80       continue
   90       continue
         go to 110
  100    continue
            info = k
  110    continue
  120 continue
  130 continue
      ipvt(n) = n
      if (abd(m,n) .eq. 0.0d0) info = n
      return
      end
      subroutine dgbsl(abd,lda,n,ml,mu,ipvt,b,job)
      integer lda,n,ml,mu,ipvt(1),job
      double precision abd(lda,1),b(1)
c
c     dgbsl solves the double precision band system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgbco or dgbfa.
c
c     on entry
c
c        abd     double precision(lda, n)
c                the output from dgbco or dgbfa.
c
c        lda     integer
c                the leading dimension of the array  abd .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c
c        mu      integer
c                number of diagonals above the main diagonal.
c
c        ipvt    integer(n)
c                the pivot vector from dgbco or dgbfa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b , where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgbco has set rcond .gt. 0.0
c        or dgbfa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c     fortran min0
c
c     internal variables
c
      double precision ddot,t
      integer k,kb,l,la,lb,lm,m,nm1
c
      m = mu + ml + 1
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve l*y = b
c
         if (ml .eq. 0) go to 30
         if (nm1 .lt. 1) go to 30
            do 20 k = 1, nm1
               lm = min0(ml,n-k)
               l = ipvt(k)
               t = b(l)
               if (l .eq. k) go to 10
                  b(l) = b(k)
                  b(k) = t
   10          continue
               call daxpy(lm,t,abd(m+1,k),1,b(k+1),1)
   20       continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/abd(m,k)
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = -b(k)
            call daxpy(lm,t,abd(la,k),1,b(lb),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = ddot(lm,abd(la,k),1,b(lb),1)
            b(k) = (b(k) - t)/abd(m,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (ml .eq. 0) go to 90
         if (nm1 .lt. 1) go to 90
            do 80 kb = 1, nm1
               k = n - kb
               lm = min0(ml,n-k)
               b(k) = b(k) + ddot(lm,abd(m+1,k),1,b(k+1),1)
               l = ipvt(k)
               if (l .eq. k) go to 70
                  t = b(l)
                  b(l) = b(k)
                  b(k) = t
   70          continue
   80       continue
   90    continue
  100 continue
      return
      end
      subroutine dgeco(a,lda,n,ipvt,rcond,z)
      integer lda,n,ipvt(1)
      double precision a(lda,1),z(1)
      double precision rcond
c
c     dgeco factors a double precision matrix by gaussian elimination
c     and estimates the condition of the matrix.
c
c     if  rcond  is not needed, dgefa is slightly faster.
c     to solve  a*x = b , follow dgeco by dgesl.
c     to compute  inverse(a)*c , follow dgeco by dgesl.
c     to compute  determinant(a) , follow dgeco by dgedi.
c     to compute  inverse(a) , follow dgeco by dgedi.
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        rcond   double precision
c                an estimate of the reciprocal condition of  a .
c                for the system  a*x = b , relative perturbations
c                in  a  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  a  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.
c
c        z       double precision(n)
c                a work vector whose contents are usually unimportant.
c                if  a  is close to a singular matrix, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     linpack dgefa
c     blas daxpy,ddot,dscal,dasum
c     fortran dabs,dmax1,dsign
c
c     internal variables
c
      double precision ddot,ek,t,wk,wkm
      double precision anorm,s,dasum,sm,ynorm
      integer info,j,k,kb,kp1,l
c
c
c     compute 1-norm of a
c
      anorm = 0.0d0
      do 10 j = 1, n
         anorm = dmax1(anorm,dasum(n,a(1,j),1))
   10 continue
c
c     factor
c
      call dgefa(a,lda,n,ipvt,info)
c
c     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c     estimate = norm(z)/norm(y) where  a*z = y  and  trans(a)*y = e .
c     trans(a)  is the transpose of a .  the components of  e  are
c     chosen to cause maximum local growth in the elements of w  where
c     trans(u)*w = e .  the vectors are frequently rescaled to avoid
c     overflow.
c
c     solve trans(u)*w = e
c
      ek = 1.0d0
      do 20 j = 1, n
         z(j) = 0.0d0
   20 continue
      do 100 k = 1, n
         if (z(k) .ne. 0.0d0) ek = dsign(ek,-z(k))
         if (dabs(ek-z(k)) .le. dabs(a(k,k))) go to 30
            s = dabs(a(k,k))/dabs(ek-z(k))
            call dscal(n,s,z,1)
            ek = s*ek
   30    continue
         wk = ek - z(k)
         wkm = -ek - z(k)
         s = dabs(wk)
         sm = dabs(wkm)
         if (a(k,k) .eq. 0.0d0) go to 40
            wk = wk/a(k,k)
            wkm = wkm/a(k,k)
         go to 50
   40    continue
            wk = 1.0d0
            wkm = 1.0d0
   50    continue
         kp1 = k + 1
         if (kp1 .gt. n) go to 90
            do 60 j = kp1, n
               sm = sm + dabs(z(j)+wkm*a(k,j))
               z(j) = z(j) + wk*a(k,j)
               s = s + dabs(z(j))
   60       continue
            if (s .ge. sm) go to 80
               t = wkm - wk
               wk = wkm
               do 70 j = kp1, n
                  z(j) = z(j) + t*a(k,j)
   70          continue
   80       continue
   90    continue
         z(k) = wk
  100 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
c
c     solve trans(l)*y = w
c
      do 120 kb = 1, n
         k = n + 1 - kb
         if (k .lt. n) z(k) = z(k) + ddot(n-k,a(k+1,k),1,z(k+1),1)
         if (dabs(z(k)) .le. 1.0d0) go to 110
            s = 1.0d0/dabs(z(k))
            call dscal(n,s,z,1)
  110    continue
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
  120 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
c
      ynorm = 1.0d0
c
c     solve l*v = y
c
      do 140 k = 1, n
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
         if (k .lt. n) call daxpy(n-k,t,a(k+1,k),1,z(k+1),1)
         if (dabs(z(k)) .le. 1.0d0) go to 130
            s = 1.0d0/dabs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
  130    continue
  140 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
c
c     solve  u*z = v
c
      do 160 kb = 1, n
         k = n + 1 - kb
         if (dabs(z(k)) .le. dabs(a(k,k))) go to 150
            s = dabs(a(k,k))/dabs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
  150    continue
         if (a(k,k) .ne. 0.0d0) z(k) = z(k)/a(k,k)
         if (a(k,k) .eq. 0.0d0) z(k) = 1.0d0
         t = -z(k)
         call daxpy(k-1,t,a(1,k),1,z(1),1)
  160 continue
c     make znorm = 1.0
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
c
      if (anorm .ne. 0.0d0) rcond = ynorm/anorm
      if (anorm .eq. 0.0d0) rcond = 0.0d0
      return
      end
      subroutine dgedi(a,lda,n,ipvt,det,work,job)
      integer lda,n,ipvt(1),job
      double precision a(lda,1),det(2),work(1)
c
c     dgedi computes the determinant and inverse of a matrix
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        work    double precision(n)
c                work vector.  contents destroyed.
c
c        job     integer
c                = 11   both determinant and inverse.
c                = 01   inverse only.
c                = 10   determinant only.
c
c     on return
c
c        a       inverse of original matrix if requested.
c                otherwise unchanged.
c
c        det     double precision(2)
c                determinant of original matrix if requested.
c                otherwise not referenced.
c                determinant = det(1) * 10.0**det(2)
c                with  1.0 .le. dabs(det(1)) .lt. 10.0
c                or  det(1) .eq. 0.0 .
c
c     error condition
c
c        a division by zero will occur if the input factor contains
c        a zero on the diagonal and the inverse is requested.
c        it will not occur if the subroutines are called correctly
c        and if dgeco has set rcond .gt. 0.0 or dgefa has set
c        info .eq. 0 .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,dswap
c     fortran dabs,mod
c
c     internal variables
c
      double precision t
      double precision ten
      integer i,j,k,kb,kp1,l,nm1
c
c
c     compute determinant
c
      if (job/10 .eq. 0) go to 70
         det(1) = 1.0d0
         det(2) = 0.0d0
         ten = 10.0d0
         do 50 i = 1, n
            if (ipvt(i) .ne. i) det(1) = -det(1)
            det(1) = a(i,i)*det(1)
c        ...exit
            if (det(1) .eq. 0.0d0) go to 60
   10       if (dabs(det(1)) .ge. 1.0d0) go to 20
               det(1) = ten*det(1)
               det(2) = det(2) - 1.0d0
            go to 10
   20       continue
   30       if (dabs(det(1)) .lt. ten) go to 40
               det(1) = det(1)/ten
               det(2) = det(2) + 1.0d0
            go to 30
   40       continue
   50    continue
   60    continue
   70 continue
c
c     compute inverse(u)
c
      if (mod(job,10) .eq. 0) go to 150
         do 100 k = 1, n
            a(k,k) = 1.0d0/a(k,k)
            t = -a(k,k)
            call dscal(k-1,t,a(1,k),1)
            kp1 = k + 1
            if (n .lt. kp1) go to 90
            do 80 j = kp1, n
               t = a(k,j)
               a(k,j) = 0.0d0
               call daxpy(k,t,a(1,k),1,a(1,j),1)
   80       continue
   90       continue
  100    continue
c
c        form inverse(u)*inverse(l)
c
         nm1 = n - 1
         if (nm1 .lt. 1) go to 140
         do 130 kb = 1, nm1
            k = n - kb
            kp1 = k + 1
            do 110 i = kp1, n
               work(i) = a(i,k)
               a(i,k) = 0.0d0
  110       continue
            do 120 j = kp1, n
               t = work(j)
               call daxpy(n,t,a(1,j),1,a(1,k),1)
  120       continue
            l = ipvt(k)
            if (l .ne. k) call dswap(n,a(1,k),1,a(1,l),1)
  130    continue
  140    continue
  150 continue
      return
      end
      subroutine dgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(1),info
      double precision a(lda,1)
c
c     dgefa factors a double precision matrix by gaussian elimination.
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c
c     internal variables
c
      double precision t
      integer idamax,j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end
      subroutine dgesl(a,lda,n,ipvt,b,job)
      integer lda,n,ipvt(1),job
      double precision a(lda,1),b(1)
c
c     dgesl solves the double precision system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgeco has set rcond .gt. 0.0
c        or dgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c
c     internal variables
c
      double precision ddot,t
      integer k,kb,l,nm1
c
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call daxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            t = ddot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end
      subroutine  dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      double precision da,dx(1)
      integer i,incx,m,mp1,n,nincx
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end
      subroutine  dswap (n,dx,incx,dy,incy)
c
c     interchanges two vectors.
c     uses unrolled loops for increments equal one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dx(ix)
        dx(ix) = dy(iy)
        dy(iy) = dtemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
c
c       clean-up loop
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        dtemp = dx(i)
        dx(i) = dy(i)
        dy(i) = dtemp
        dtemp = dx(i + 1)
        dx(i + 1) = dy(i + 1)
        dy(i + 1) = dtemp
        dtemp = dx(i + 2)
        dx(i + 2) = dy(i + 2)
        dy(i + 2) = dtemp
   50 continue
      return
      end
      integer function idamax(n,dx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      double precision dx(1),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end
C*****END precision > double
C
C*****precision > single
C      SUBROUTINE SGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
C     $                   BETA, C, LDC )
C*     .. Scalar Arguments ..
C      CHARACTER*1        TRANSA, TRANSB
C      INTEGER            M, N, K, LDA, LDB, LDC
C      REAL               ALPHA, BETA
C*     .. Array Arguments ..
C      REAL               A( LDA, * ), B( LDB, * ), C( LDC, * )
C*     ..
C*
C*  Purpose
C*  =======
C*
C*  SGEMM  performs one of the matrix-matrix operations
C*
C*     C := alpha*op( A )*op( B ) + beta*C,
C*
C*  where  op( X ) is one of
C*
C*     op( X ) = X   or   op( X ) = X',
C*
C*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
C*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
C*
C*  Parameters
C*  ==========
C*
C*  TRANSA - CHARACTER*1.
C*           On entry, TRANSA specifies the form of op( A ) to be used in
C*           the matrix multiplication as follows:
C*
C*              TRANSA = 'N' or 'n',  op( A ) = A.
C*
C*              TRANSA = 'T' or 't',  op( A ) = A'.
C*
C*              TRANSA = 'C' or 'c',  op( A ) = A'.
C*
C*           Unchanged on exit.
C*
C*  TRANSB - CHARACTER*1.
C*           On entry, TRANSB specifies the form of op( B ) to be used in
C*           the matrix multiplication as follows:
C*
C*              TRANSB = 'N' or 'n',  op( B ) = B.
C*
C*              TRANSB = 'T' or 't',  op( B ) = B'.
C*
C*              TRANSB = 'C' or 'c',  op( B ) = B'.
C*
C*           Unchanged on exit.
C*
C*  M      - INTEGER.
C*           On entry,  M  specifies  the number  of rows  of the  matrix
C*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
C*           Unchanged on exit.
C*
C*  N      - INTEGER.
C*           On entry,  N  specifies the number  of columns of the matrix
C*           op( B ) and the number of columns of the matrix C. N must be
C*           at least zero.
C*           Unchanged on exit.
C*
C*  K      - INTEGER.
C*           On entry,  K  specifies  the number of columns of the matrix
C*           op( A ) and the number of rows of the matrix op( B ). K must
C*           be at least  zero.
C*           Unchanged on exit.
C*
C*  ALPHA  - REAL            .
C*           On entry, ALPHA specifies the scalar alpha.
C*           Unchanged on exit.
C*
C*  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
C*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
C*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
C*           part of the array  A  must contain the matrix  A,  otherwise
C*           the leading  k by m  part of the array  A  must contain  the
C*           matrix A.
C*           Unchanged on exit.
C*
C*  LDA    - INTEGER.
C*           On entry, LDA specifies the first dimension of A as declared
C*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
C*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
C*           least  max( 1, k ).
C*           Unchanged on exit.
C*
C*  B      - REAL             array of DIMENSION ( LDB, kb ), where kb is
C*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
C*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
C*           part of the array  B  must contain the matrix  B,  otherwise
C*           the leading  n by k  part of the array  B  must contain  the
C*           matrix B.
C*           Unchanged on exit.
C*
C*  LDB    - INTEGER.
C*           On entry, LDB specifies the first dimension of B as declared
C*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
C*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
C*           least  max( 1, n ).
C*           Unchanged on exit.
C*
C*  BETA   - REAL            .
C*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
C*           supplied as zero then C need not be set on input.
C*           Unchanged on exit.
C*
C*  C      - REAL             array of DIMENSION ( LDC, n ).
C*           Before entry, the leading  m by n  part of the array  C must
C*           contain the matrix  C,  except when  beta  is zero, in which
C*           case C need not be set on entry.
C*           On exit, the array  C  is overwritten by the  m by n  matrix
C*           ( alpha*op( A )*op( B ) + beta*C ).
C*
C*  LDC    - INTEGER.
C*           On entry, LDC specifies the first dimension of C as declared
C*           in  the  calling  (sub)  program.   LDC  must  be  at  least
C*           max( 1, m ).
C*           Unchanged on exit.
C*
C*
C*  Level 3 Blas routine.
C*
C*  -- Written on 8-February-1989.
C*     Jack Dongarra, Argonne National Laboratory.
C*     Iain Duff, AERE Harwell.
C*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
C*     Sven Hammarling, Numerical Algorithms Group Ltd.
C*
C*
C*     .. External Functions ..
C      LOGICAL            LSAME
C      EXTERNAL           LSAME
C*     .. External Subroutines ..
C      EXTERNAL           XERBLA
C*     .. Intrinsic Functions ..
C      INTRINSIC          MAX
C*     .. Local Scalars ..
C      LOGICAL            NOTA, NOTB
C      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
C      REAL               TEMP
C*     .. Parameters ..
C      REAL               ONE         , ZERO
C      PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
C*     ..
C*     .. Executable Statements ..
C*
C*     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
C*     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
C*     and  columns of  A  and the  number of  rows  of  B  respectively.
C*
C      NOTA  = LSAME( TRANSA, 'N' )
C      NOTB  = LSAME( TRANSB, 'N' )
C      IF( NOTA )THEN
C         NROWA = M
C         NCOLA = K
C      ELSE
C         NROWA = K
C         NCOLA = M
C      END IF
C      IF( NOTB )THEN
C         NROWB = K
C      ELSE
C         NROWB = N
C      END IF
C*
C*     Test the input parameters.
C*
C      INFO = 0
C      IF(      ( .NOT.NOTA                 ).AND.
C     $         ( .NOT.LSAME( TRANSA, 'C' ) ).AND.
C     $         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
C         INFO = 1
C      ELSE IF( ( .NOT.NOTB                 ).AND.
C     $         ( .NOT.LSAME( TRANSB, 'C' ) ).AND.
C     $         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
C         INFO = 2
C      ELSE IF( M  .LT.0               )THEN
C         INFO = 3
C      ELSE IF( N  .LT.0               )THEN
C         INFO = 4
C      ELSE IF( K  .LT.0               )THEN
C         INFO = 5
C      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
C         INFO = 8
C      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
C         INFO = 10
C      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
C         INFO = 13
C      END IF
C      IF( INFO.NE.0 )THEN
C         CALL XERBLA( 'SGEMM ', INFO )
C         RETURN
C      END IF
C*
C*     Quick return if possible.
C*
C      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
C     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
C     $   RETURN
C*
C*     And if  alpha.eq.zero.
C*
C      IF( ALPHA.EQ.ZERO )THEN
C         IF( BETA.EQ.ZERO )THEN
C            DO 20, J = 1, N
C               DO 10, I = 1, M
C                  C( I, J ) = ZERO
C   10          CONTINUE
C   20       CONTINUE
C         ELSE
C            DO 40, J = 1, N
C               DO 30, I = 1, M
C                  C( I, J ) = BETA*C( I, J )
C   30          CONTINUE
C   40       CONTINUE
C         END IF
C         RETURN
C      END IF
C*
C*     Start the operations.
C*
C      IF( NOTB )THEN
C         IF( NOTA )THEN
C*
C*           Form  C := alpha*A*B + beta*C.
C*
C            DO 90, J = 1, N
C               IF( BETA.EQ.ZERO )THEN
C                  DO 50, I = 1, M
C                     C( I, J ) = ZERO
C   50             CONTINUE
C               ELSE IF( BETA.NE.ONE )THEN
C                  DO 60, I = 1, M
C                     C( I, J ) = BETA*C( I, J )
C   60             CONTINUE
C               END IF
C               DO 80, L = 1, K
C                  IF( B( L, J ).NE.ZERO )THEN
C                     TEMP = ALPHA*B( L, J )
C                     DO 70, I = 1, M
C                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
C   70                CONTINUE
C                  END IF
C   80          CONTINUE
C   90       CONTINUE
C         ELSE
C*
C*           Form  C := alpha*A'*B + beta*C
C*
C            DO 120, J = 1, N
C               DO 110, I = 1, M
C                  TEMP = ZERO
C                  DO 100, L = 1, K
C                     TEMP = TEMP + A( L, I )*B( L, J )
C  100             CONTINUE
C                  IF( BETA.EQ.ZERO )THEN
C                     C( I, J ) = ALPHA*TEMP
C                  ELSE
C                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
C                  END IF
C  110          CONTINUE
C  120       CONTINUE
C         END IF
C      ELSE
C         IF( NOTA )THEN
C*
C*           Form  C := alpha*A*B' + beta*C
C*
C            DO 170, J = 1, N
C               IF( BETA.EQ.ZERO )THEN
C                  DO 130, I = 1, M
C                     C( I, J ) = ZERO
C  130             CONTINUE
C               ELSE IF( BETA.NE.ONE )THEN
C                  DO 140, I = 1, M
C                     C( I, J ) = BETA*C( I, J )
C  140             CONTINUE
C               END IF
C               DO 160, L = 1, K
C                  IF( B( J, L ).NE.ZERO )THEN
C                     TEMP = ALPHA*B( J, L )
C                     DO 150, I = 1, M
C                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
C  150                CONTINUE
C                  END IF
C  160          CONTINUE
C  170       CONTINUE
C         ELSE
C*
C*           Form  C := alpha*A'*B' + beta*C
C*
C            DO 200, J = 1, N
C               DO 190, I = 1, M
C                  TEMP = ZERO
C                  DO 180, L = 1, K
C                     TEMP = TEMP + A( L, I )*B( J, L )
C  180             CONTINUE
C                  IF( BETA.EQ.ZERO )THEN
C                     C( I, J ) = ALPHA*TEMP
C                  ELSE
C                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
C                  END IF
C  190          CONTINUE
C  200       CONTINUE
C         END IF
C      END IF
C*
C      RETURN
C*
C*     End of SGEMM .
C*
C      END
C      SUBROUTINE SGEMV ( TRANS, M, N, ALPHA, A, LDA, X, INCX,
C     $                   BETA, Y, INCY )
C*     .. Scalar Arguments ..
C      REAL               ALPHA, BETA
C      INTEGER            INCX, INCY, LDA, M, N
C      CHARACTER*1        TRANS
C*     .. Array Arguments ..
C      REAL               A( LDA, * ), X( * ), Y( * )
C*     ..
C*
C*  Purpose
C*  =======
C*
C*  SGEMV  performs one of the matrix-vector operations
C*
C*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
C*
C*  where alpha and beta are scalars, x and y are vectors and A is an
C*  m by n matrix.
C*
C*  Parameters
C*  ==========
C*
C*  TRANS  - CHARACTER*1.
C*           On entry, TRANS specifies the operation to be performed as
C*           follows:
C*
C*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
C*
C*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
C*
C*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
C*
C*           Unchanged on exit.
C*
C*  M      - INTEGER.
C*           On entry, M specifies the number of rows of the matrix A.
C*           M must be at least zero.
C*           Unchanged on exit.
C*
C*  N      - INTEGER.
C*           On entry, N specifies the number of columns of the matrix A.
C*           N must be at least zero.
C*           Unchanged on exit.
C*
C*  ALPHA  - REAL            .
C*           On entry, ALPHA specifies the scalar alpha.
C*           Unchanged on exit.
C*
C*  A      - REAL             array of DIMENSION ( LDA, n ).
C*           Before entry, the leading m by n part of the array A must
C*           contain the matrix of coefficients.
C*           Unchanged on exit.
C*
C*  LDA    - INTEGER.
C*           On entry, LDA specifies the first dimension of A as declared
C*           in the calling (sub) program. LDA must be at least
C*           max( 1, m ).
C*           Unchanged on exit.
C*
C*  X      - REAL             array of DIMENSION at least
C*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
C*           and at least
C*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
C*           Before entry, the incremented array X must contain the
C*           vector x.
C*           Unchanged on exit.
C*
C*  INCX   - INTEGER.
C*           On entry, INCX specifies the increment for the elements of
C*           X. INCX must not be zero.
C*           Unchanged on exit.
C*
C*  BETA   - REAL            .
C*           On entry, BETA specifies the scalar beta. When BETA is
C*           supplied as zero then Y need not be set on input.
C*           Unchanged on exit.
C*
C*  Y      - REAL             array of DIMENSION at least
C*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
C*           and at least
C*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
C*           Before entry with BETA non-zero, the incremented array Y
C*           must contain the vector y. On exit, Y is overwritten by the
C*           updated vector y.
C*
C*  INCY   - INTEGER.
C*           On entry, INCY specifies the increment for the elements of
C*           Y. INCY must not be zero.
C*           Unchanged on exit.
C*
C*
C*  Level 2 Blas routine.
C*
C*  -- Written on 22-October-1986.
C*     Jack Dongarra, Argonne National Lab.
C*     Jeremy Du Croz, Nag Central Office.
C*     Sven Hammarling, Nag Central Office.
C*     Richard Hanson, Sandia National Labs.
C*
C*
C*     .. Parameters ..
C      REAL               ONE         , ZERO
C      PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
C*     .. Local Scalars ..
C      REAL               TEMP
C      INTEGER            I, INFO, IX, IY, J, JX, JY, KX, KY, LENX, LENY
C*     .. External Functions ..
C      LOGICAL            LSAME
C      EXTERNAL           LSAME
C*     .. External Subroutines ..
C      EXTERNAL           XERBLA
C*     .. Intrinsic Functions ..
C      INTRINSIC          MAX
C*     ..
C*     .. Executable Statements ..
C*
C*     Test the input parameters.
C*
C      INFO = 0
C      IF     ( .NOT.LSAME( TRANS, 'N' ).AND.
C     $         .NOT.LSAME( TRANS, 'T' ).AND.
C     $         .NOT.LSAME( TRANS, 'C' )      )THEN
C         INFO = 1
C      ELSE IF( M.LT.0 )THEN
C         INFO = 2
C      ELSE IF( N.LT.0 )THEN
C         INFO = 3
C      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
C         INFO = 6
C      ELSE IF( INCX.EQ.0 )THEN
C         INFO = 8
C      ELSE IF( INCY.EQ.0 )THEN
C         INFO = 11
C      END IF
C      IF( INFO.NE.0 )THEN
C         CALL XERBLA( 'SGEMV ', INFO )
C         RETURN
C      END IF
C*
C*     Quick return if possible.
C*
C      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
C     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
C     $   RETURN
C*
C*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
C*     up the start points in  X  and  Y.
C*
C      IF( LSAME( TRANS, 'N' ) )THEN
C         LENX = N
C         LENY = M
C      ELSE
C         LENX = M
C         LENY = N
C      END IF
C      IF( INCX.GT.0 )THEN
C         KX = 1
C      ELSE
C         KX = 1 - ( LENX - 1 )*INCX
C      END IF
C      IF( INCY.GT.0 )THEN
C         KY = 1
C      ELSE
C         KY = 1 - ( LENY - 1 )*INCY
C      END IF
C*
C*     Start the operations. In this version the elements of A are
C*     accessed sequentially with one pass through A.
C*
C*     First form  y := beta*y.
C*
C      IF( BETA.NE.ONE )THEN
C         IF( INCY.EQ.1 )THEN
C            IF( BETA.EQ.ZERO )THEN
C               DO 10, I = 1, LENY
C                  Y( I ) = ZERO
C   10          CONTINUE
C            ELSE
C               DO 20, I = 1, LENY
C                  Y( I ) = BETA*Y( I )
C   20          CONTINUE
C            END IF
C         ELSE
C            IY = KY
C            IF( BETA.EQ.ZERO )THEN
C               DO 30, I = 1, LENY
C                  Y( IY ) = ZERO
C                  IY      = IY   + INCY
C   30          CONTINUE
C            ELSE
C               DO 40, I = 1, LENY
C                  Y( IY ) = BETA*Y( IY )
C                  IY      = IY           + INCY
C   40          CONTINUE
C            END IF
C         END IF
C      END IF
C      IF( ALPHA.EQ.ZERO )
C     $   RETURN
C      IF( LSAME( TRANS, 'N' ) )THEN
C*
C*        Form  y := alpha*A*x + y.
C*
C         JX = KX
C         IF( INCY.EQ.1 )THEN
C            DO 60, J = 1, N
C               IF( X( JX ).NE.ZERO )THEN
C                  TEMP = ALPHA*X( JX )
C                  DO 50, I = 1, M
C                     Y( I ) = Y( I ) + TEMP*A( I, J )
C   50             CONTINUE
C               END IF
C               JX = JX + INCX
C   60       CONTINUE
C         ELSE
C            DO 80, J = 1, N
C               IF( X( JX ).NE.ZERO )THEN
C                  TEMP = ALPHA*X( JX )
C                  IY   = KY
C                  DO 70, I = 1, M
C                     Y( IY ) = Y( IY ) + TEMP*A( I, J )
C                     IY      = IY      + INCY
C   70             CONTINUE
C               END IF
C               JX = JX + INCX
C   80       CONTINUE
C         END IF
C      ELSE
C*
C*        Form  y := alpha*A'*x + y.
C*
C         JY = KY
C         IF( INCX.EQ.1 )THEN
C            DO 100, J = 1, N
C               TEMP = ZERO
C               DO 90, I = 1, M
C                  TEMP = TEMP + A( I, J )*X( I )
C   90          CONTINUE
C               Y( JY ) = Y( JY ) + ALPHA*TEMP
C               JY      = JY      + INCY
C  100       CONTINUE
C         ELSE
C            DO 120, J = 1, N
C               TEMP = ZERO
C               IX   = KX
C               DO 110, I = 1, M
C                  TEMP = TEMP + A( I, J )*X( IX )
C                  IX   = IX   + INCX
C  110          CONTINUE
C               Y( JY ) = Y( JY ) + ALPHA*TEMP
C               JY      = JY      + INCY
C  120       CONTINUE
C         END IF
C      END IF
C*
C      RETURN
C*
C*     End of SGEMV .
C*
C      END
C      SUBROUTINE SGER  ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA )
C*     .. Scalar Arguments ..
C      REAL               ALPHA
C      INTEGER            INCX, INCY, LDA, M, N
C*     .. Array Arguments ..
C      REAL               A( LDA, * ), X( * ), Y( * )
C*     ..
C*
C*  Purpose
C*  =======
C*
C*  SGER   performs the rank 1 operation
C*
C*     A := alpha*x*y' + A,
C*
C*  where alpha is a scalar, x is an m element vector, y is an n element
C*  vector and A is an m by n matrix.
C*
C*  Parameters
C*  ==========
C*
C*  M      - INTEGER.
C*           On entry, M specifies the number of rows of the matrix A.
C*           M must be at least zero.
C*           Unchanged on exit.
C*
C*  N      - INTEGER.
C*           On entry, N specifies the number of columns of the matrix A.
C*           N must be at least zero.
C*           Unchanged on exit.
C*
C*  ALPHA  - REAL            .
C*           On entry, ALPHA specifies the scalar alpha.
C*           Unchanged on exit.
C*
C*  X      - REAL             array of dimension at least
C*           ( 1 + ( m - 1 )*abs( INCX ) ).
C*           Before entry, the incremented array X must contain the m
C*           element vector x.
C*           Unchanged on exit.
C*
C*  INCX   - INTEGER.
C*           On entry, INCX specifies the increment for the elements of
C*           X. INCX must not be zero.
C*           Unchanged on exit.
C*
C*  Y      - REAL             array of dimension at least
C*           ( 1 + ( n - 1 )*abs( INCY ) ).
C*           Before entry, the incremented array Y must contain the n
C*           element vector y.
C*           Unchanged on exit.
C*
C*  INCY   - INTEGER.
C*           On entry, INCY specifies the increment for the elements of
C*           Y. INCY must not be zero.
C*           Unchanged on exit.
C*
C*  A      - REAL             array of DIMENSION ( LDA, n ).
C*           Before entry, the leading m by n part of the array A must
C*           contain the matrix of coefficients. On exit, A is
C*           overwritten by the updated matrix.
C*
C*  LDA    - INTEGER.
C*           On entry, LDA specifies the first dimension of A as declared
C*           in the calling (sub) program. LDA must be at least
C*           max( 1, m ).
C*           Unchanged on exit.
C*
C*
C*  Level 2 Blas routine.
C*
C*  -- Written on 22-October-1986.
C*     Jack Dongarra, Argonne National Lab.
C*     Jeremy Du Croz, Nag Central Office.
C*     Sven Hammarling, Nag Central Office.
C*     Richard Hanson, Sandia National Labs.
C*
C*
C*     .. Parameters ..
C      REAL               ZERO
C      PARAMETER        ( ZERO = 0.0E+0 )
C*     .. Local Scalars ..
C      REAL               TEMP
C      INTEGER            I, INFO, IX, J, JY, KX
C*     .. External Subroutines ..
C      EXTERNAL           XERBLA
C*     .. Intrinsic Functions ..
C      INTRINSIC          MAX
C*     ..
C*     .. Executable Statements ..
C*
C*     Test the input parameters.
C*
C      INFO = 0
C      IF     ( M.LT.0 )THEN
C         INFO = 1
C      ELSE IF( N.LT.0 )THEN
C         INFO = 2
C      ELSE IF( INCX.EQ.0 )THEN
C         INFO = 5
C      ELSE IF( INCY.EQ.0 )THEN
C         INFO = 7
C      ELSE IF( LDA.LT.MAX( 1, M ) )THEN
C         INFO = 9
C      END IF
C      IF( INFO.NE.0 )THEN
C         CALL XERBLA( 'SGER  ', INFO )
C         RETURN
C      END IF
C*
C*     Quick return if possible.
C*
C      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )
C     $   RETURN
C*
C*     Start the operations. In this version the elements of A are
C*     accessed sequentially with one pass through A.
C*
C      IF( INCY.GT.0 )THEN
C         JY = 1
C      ELSE
C         JY = 1 - ( N - 1 )*INCY
C      END IF
C      IF( INCX.EQ.1 )THEN
C         DO 20, J = 1, N
C            IF( Y( JY ).NE.ZERO )THEN
C               TEMP = ALPHA*Y( JY )
C               DO 10, I = 1, M
C                  A( I, J ) = A( I, J ) + X( I )*TEMP
C   10          CONTINUE
C            END IF
C            JY = JY + INCY
C   20    CONTINUE
C      ELSE
C         IF( INCX.GT.0 )THEN
C            KX = 1
C         ELSE
C            KX = 1 - ( M - 1 )*INCX
C         END IF
C         DO 40, J = 1, N
C            IF( Y( JY ).NE.ZERO )THEN
C               TEMP = ALPHA*Y( JY )
C               IX   = KX
C               DO 30, I = 1, M
C                  A( I, J ) = A( I, J ) + X( IX )*TEMP
C                  IX        = IX        + INCX
C   30          CONTINUE
C            END IF
C            JY = JY + INCY
C   40    CONTINUE
C      END IF
C*
C      RETURN
C*
C*     End of SGER  .
C*
C      END
C      SUBROUTINE SGETF2( M, N, A, LDA, IPIV, INFO )
C*
C*  -- LAPACK routine (version 1.1) --
C*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C*     Courant Institute, Argonne National Lab, and Rice University
C*     June 30, 1992
C*
C*     .. Scalar Arguments ..
C      INTEGER            INFO, LDA, M, N
C*     ..
C*     .. Array Arguments ..
C      INTEGER            IPIV( * )
C      REAL               A( LDA, * )
C*     ..
C*
C*  Purpose
C*  =======
C*
C*  SGETF2 computes an LU factorization of a general m-by-n matrix A
C*  using partial pivoting with row interchanges.
C*
C*  The factorization has the form
C*     A = P * L * U
C*  where P is a permutation matrix, L is lower triangular with unit
C*  diagonal elements (lower trapezoidal if m > n), and U is upper
C*  triangular (upper trapezoidal if m < n).
C*
C*  This is the right-looking Level 2 BLAS version of the algorithm.
C*
C*  Arguments
C*  =========
C*
C*  M       (input) INTEGER
C*          The number of rows of the matrix A.  M >= 0.
C*
C*  N       (input) INTEGER
C*          The number of columns of the matrix A.  N >= 0.
C*
C*  A       (input/output) REAL array, dimension (LDA,N)
C*          On entry, the m by n matrix to be factored.
C*          On exit, the factors L and U from the factorization
C*          A = P*L*U; the unit diagonal elements of L are not stored.
C*
C*  LDA     (input) INTEGER
C*          The leading dimension of the array A.  LDA >= max(1,M).
C*
C*  IPIV    (output) INTEGER array, dimension (min(M,N))
C*          The pivot indices; for 1 <= i <= min(M,N), row i of the
C*          matrix was interchanged with row IPIV(i).
C*
C*  INFO    (output) INTEGER
C*          = 0: successful exit
C*          < 0: if INFO = -k, the k-th argument had an illegal value
C*          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
C*               has been completed, but the factor U is exactly
C*               singular, and division by zero will occur if it is used
C*               to solve a system of equations.
C*
C*  =====================================================================
C*
C*     .. Parameters ..
C      REAL               ONE, ZERO
C      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
C*     ..
C*     .. Local Scalars ..
C      INTEGER            J, JP
C*     ..
C*     .. External Functions ..
C      INTEGER            ISAMAX
C      EXTERNAL           ISAMAX
C*     ..
C*     .. External Subroutines ..
C      EXTERNAL           SGER, SSCAL, SSWAP, XERBLA
C*     ..
C*     .. Intrinsic Functions ..
C      INTRINSIC          MAX, MIN
C*     ..
C*     .. Executable Statements ..
C*
C*     Test the input parameters.
C*
C      INFO = 0
C      IF( M.LT.0 ) THEN
C         INFO = -1
C      ELSE IF( N.LT.0 ) THEN
C         INFO = -2
C      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
C         INFO = -4
C      END IF
C      IF( INFO.NE.0 ) THEN
C         CALL XERBLA( 'SGETF2', -INFO )
C         RETURN
C      END IF
C*
C*     Quick return if possible
C*
C      IF( M.EQ.0 .OR. N.EQ.0 )
C     $   RETURN
C*
C      DO 10 J = 1, MIN( M, N )
C*
C*        Find pivot and test for singularity.
C*
C         JP = J - 1 + ISAMAX( M-J+1, A( J, J ), 1 )
C         IPIV( J ) = JP
C         IF( A( JP, J ).NE.ZERO ) THEN
C*
C*           Apply the interchange to columns 1:N.
C*
C            IF( JP.NE.J )
C     $         CALL SSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
C*
C*           Compute elements J+1:M of J-th column.
C*
C            IF( J.LT.M )
C     $         CALL SSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
C*
C         ELSE IF( INFO.EQ.0 ) THEN
C*
C            INFO = J
C         END IF
C*
C         IF( J.LT.MIN( M, N ) ) THEN
C*
C*           Update trailing submatrix.
C*
C            CALL SGER( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ), LDA,
C     $                 A( J+1, J+1 ), LDA )
C         END IF
C   10 CONTINUE
C      RETURN
C*
C*     End of SGETF2
C*
C      END
C      SUBROUTINE SGETRF( M, N, A, LDA, IPIV, INFO )
C*
C*  -- LAPACK routine (version 1.1) --
C*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C*     Courant Institute, Argonne National Lab, and Rice University
C*     March 31, 1993
C*
C*     .. Scalar Arguments ..
C      INTEGER            INFO, LDA, M, N
C*     ..
C*     .. Array Arguments ..
C      INTEGER            IPIV( * )
C      REAL               A( LDA, * )
C*     ..
C*
C*  Purpose
C*  =======
C*
C*  SGETRF computes an LU factorization of a general M-by-N matrix A
C*  using partial pivoting with row interchanges.
C*
C*  The factorization has the form
C*     A = P * L * U
C*  where P is a permutation matrix, L is lower triangular with unit
C*  diagonal elements (lower trapezoidal if m > n), and U is upper
C*  triangular (upper trapezoidal if m < n).
C*
C*  This is the right-looking Level 3 BLAS version of the algorithm.
C*
C*  Arguments
C*  =========
C*
C*  M       (input) INTEGER
C*          The number of rows of the matrix A.  M >= 0.
C*
C*  N       (input) INTEGER
C*          The number of columns of the matrix A.  N >= 0.
C*
C*  A       (input/output) REAL array, dimension (LDA,N)
C*          On entry, the M-by-N matrix to be factored.
C*          On exit, the factors L and U from the factorization
C*          A = P*L*U; the unit diagonal elements of L are not stored.
C*
C*  LDA     (input) INTEGER
C*          The leading dimension of the array A.  LDA >= max(1,M).
C*
C*  IPIV    (output) INTEGER array, dimension (min(M,N))
C*          The pivot indices; for 1 <= i <= min(M,N), row i of the
C*          matrix was interchanged with row IPIV(i).
C*
C*  INFO    (output) INTEGER
C*          = 0:  successful exit
C*          < 0:  if INFO = -i, the i-th argument had an illegal value
C*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
C*                has been completed, but the factor U is exactly
C*                singular, and division by zero will occur if it is used
C*                to solve a system of equations.
C*
C*  =====================================================================
C*
C*     .. Parameters ..
C      REAL               ONE
C      PARAMETER          ( ONE = 1.0E+0 )
C*     ..
C*     .. Local Scalars ..
C      INTEGER            I, IINFO, J, JB, NB
C*     ..
C*     .. External Subroutines ..
C      EXTERNAL           SGEMM, SGETF2, SLASWP, STRSM, XERBLA
C*     ..
C*     .. External Functions ..
C      INTEGER            ILAENV
C      EXTERNAL           ILAENV
C*     ..
C*     .. Intrinsic Functions ..
C      INTRINSIC          MAX, MIN
C*     ..
C*     .. Executable Statements ..
C*
C*     Test the input parameters.
C*
C      INFO = 0
C      IF( M.LT.0 ) THEN
C         INFO = -1
C      ELSE IF( N.LT.0 ) THEN
C         INFO = -2
C      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
C         INFO = -4
C      END IF
C      IF( INFO.NE.0 ) THEN
C         CALL XERBLA( 'SGETRF', -INFO )
C         RETURN
C      END IF
C*
C*     Quick return if possible
C*
C      IF( M.EQ.0 .OR. N.EQ.0 )
C     $   RETURN
C*
C*     Determine the block size for this environment.
C*
C      NB = ILAENV( 1, 'SGETRF', ' ', M, N, -1, -1 )
C      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
C*
C*        Use unblocked code.
C*
C         CALL SGETF2( M, N, A, LDA, IPIV, INFO )
C      ELSE
C*
C*        Use blocked code.
C*
C         DO 20 J = 1, MIN( M, N ), NB
C            JB = MIN( MIN( M, N )-J+1, NB )
C*
C*           Factor diagonal and subdiagonal blocks and test for exact
C*           singularity.
C*
C            CALL SGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
C*
C*           Adjust INFO and the pivot indices.
C*
C            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
C     $         INFO = IINFO + J - 1
C            DO 10 I = J, MIN( M, J+JB-1 )
C               IPIV( I ) = J - 1 + IPIV( I )
C   10       CONTINUE
C*
C*           Apply interchanges to columns 1:J-1.
C*
C            CALL SLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
C*
C            IF( J+JB.LE.N ) THEN
C*
C*              Apply interchanges to columns J+JB:N.
C*
C               CALL SLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,
C     $                      IPIV, 1 )
C*
C*              Compute block row of U.
C*
C               CALL STRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,
C     $                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),
C     $                     LDA )
C               IF( J+JB.LE.M ) THEN
C*
C*                 Update trailing submatrix.
C*
C                  CALL SGEMM( 'No transpose', 'No transpose', M-J-JB+1,
C     $                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
C     $                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
C     $                        LDA )
C               END IF
C            END IF
C   20    CONTINUE
C      END IF
C      RETURN
C*
C*     End of SGETRF
C*
C      END
C      SUBROUTINE SGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
C*
C*  -- LAPACK routine (version 1.1) --
C*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C*     Courant Institute, Argonne National Lab, and Rice University
C*     March 31, 1993
C*
C*     .. Scalar Arguments ..
C      INTEGER            INFO, LDA, LWORK, N
C*     ..
C*     .. Array Arguments ..
C      INTEGER            IPIV( * )
C      REAL               A( LDA, * ), WORK( LWORK )
C*     ..
C*
C*  Purpose
C*  =======
C*
C*  SGETRI computes the inverse of a matrix using the LU factorization
C*  computed by SGETRF.
C*
C*  This method inverts U and then computes inv(A) by solving the system
C*  inv(A)*L = inv(U) for inv(A).
C*
C*  Arguments
C*  =========
C*
C*  N       (input) INTEGER
C*          The order of the matrix A.  N >= 0.
C*
C*  A       (input/output) REAL array, dimension (LDA,N)
C*          On entry, the factors L and U from the factorization
C*          A = P*L*U as computed by SGETRF.
C*          On exit, if INFO = 0, the inverse of the original matrix A.
C*
C*  LDA     (input) INTEGER
C*          The leading dimension of the array A.  LDA >= max(1,N).
C*
C*  IPIV    (input) INTEGER array, dimension (N)
C*          The pivot indices from SGETRF; for 1<=i<=N, row i of the
C*          matrix was interchanged with row IPIV(i).
C*
C*  WORK    (workspace) REAL array, dimension (LWORK)
C*          On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
C*
C*  LWORK   (input) INTEGER
C*          The dimension of the array WORK.  LWORK >= max(1,N).
C*          For optimal performance LWORK >= N*NB, where NB is
C*          the optimal blocksize returned by ILAENV.
C*
C*  INFO    (output) INTEGER
C*          = 0:  successful exit
C*          < 0:  if INFO = -i, the i-th argument had an illegal value
C*          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
C*                singular and its inverse could not be computed.
C*
C*  =====================================================================
C*
C*     .. Parameters ..
C      REAL               ZERO, ONE
C      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
C*     ..
C*     .. Local Scalars ..
C      INTEGER            I, IWS, J, JB, JJ, JP, LDWORK, NB, NBMIN, NN
C*     ..
C*     .. External Functions ..
C      INTEGER            ILAENV
C      EXTERNAL           ILAENV
C*     ..
C*     .. External Subroutines ..
C      EXTERNAL           SGEMM, SGEMV, SSWAP, STRSM, STRTRI, XERBLA
C*     ..
C*     .. Intrinsic Functions ..
C      INTRINSIC          MAX, MIN
C*     ..
C*     .. Executable Statements ..
C*
C*     Test the input parameters.
C*
C      INFO = 0
C      WORK( 1 ) = MAX( N, 1 )
C      IF( N.LT.0 ) THEN
C         INFO = -1
C      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
C         INFO = -3
C      ELSE IF( LWORK.LT.MAX( 1, N ) ) THEN
C         INFO = -6
C      END IF
C      IF( INFO.NE.0 ) THEN
C         CALL XERBLA( 'SGETRI', -INFO )
C         RETURN
C      END IF
C*
C*     Quick return if possible
C*
C      IF( N.EQ.0 )
C     $   RETURN
C*
C*     Form inv(U).  If INFO > 0 from STRTRI, then U is singular,
C*     and the inverse is not computed.
C*
C      CALL STRTRI( 'Upper', 'Non-unit', N, A, LDA, INFO )
C      IF( INFO.GT.0 )
C     $   RETURN
C*
C*     Determine the block size for this environment.
C*
C      NB = ILAENV( 1, 'SGETRI', ' ', N, -1, -1, -1 )
C      NBMIN = 2
C      LDWORK = N
C      IF( NB.GT.1 .AND. NB.LT.N ) THEN
C         IWS = MAX( LDWORK*NB, 1 )
C         IF( LWORK.LT.IWS ) THEN
C            NB = LWORK / LDWORK
C            NBMIN = MAX( 2, ILAENV( 2, 'SGETRI', ' ', N, -1, -1, -1 ) )
C         END IF
C      ELSE
C         IWS = N
C      END IF
C*
C*     Solve the equation inv(A)*L = inv(U) for inv(A).
C*
C      IF( NB.LT.NBMIN .OR. NB.GE.N ) THEN
C*
C*        Use unblocked code.
C*
C         DO 20 J = N, 1, -1
C*
C*           Copy current column of L to WORK and replace with zeros.
C*
C            DO 10 I = J + 1, N
C               WORK( I ) = A( I, J )
C               A( I, J ) = ZERO
C   10       CONTINUE
C*
C*           Compute current column of inv(A).
C*
C            IF( J.LT.N )
C     $         CALL SGEMV( 'No transpose', N, N-J, -ONE, A( 1, J+1 ),
C     $                     LDA, WORK( J+1 ), 1, ONE, A( 1, J ), 1 )
C   20    CONTINUE
C      ELSE
C*
C*        Use blocked code.
C*
C         NN = ( ( N-1 ) / NB )*NB + 1
C         DO 50 J = NN, 1, -NB
C            JB = MIN( NB, N-J+1 )
C*
C*           Copy current block column of L to WORK and replace with
C*           zeros.
C*
C            DO 40 JJ = J, J + JB - 1
C               DO 30 I = JJ + 1, N
C                  WORK( I+( JJ-J )*LDWORK ) = A( I, JJ )
C                  A( I, JJ ) = ZERO
C   30          CONTINUE
C   40       CONTINUE
C*
C*           Compute current block column of inv(A).
C*
C            IF( J+JB.LE.N )
C     $         CALL SGEMM( 'No transpose', 'No transpose', N, JB,
C     $                     N-J-JB+1, -ONE, A( 1, J+JB ), LDA,
C     $                     WORK( J+JB ), LDWORK, ONE, A( 1, J ), LDA )
C            CALL STRSM( 'Right', 'Lower', 'No transpose', 'Unit', N, JB,
C     $                  ONE, WORK( J ), LDWORK, A( 1, J ), LDA )
C   50    CONTINUE
C      END IF
C*
C*     Apply column interchanges.
C*
C      DO 60 J = N - 1, 1, -1
C         JP = IPIV( J )
C         IF( JP.NE.J )
C     $      CALL SSWAP( N, A( 1, J ), 1, A( 1, JP ), 1 )
C   60 CONTINUE
C*
C      WORK( 1 ) = IWS
C      RETURN
C*
C*     End of SGETRI
C*
C      END
C      SUBROUTINE SGETRS( TRANS, N, NRHS, A, LDA, IPIV, B, LDB, INFO )
C*
C*  -- LAPACK routine (version 1.1) --
C*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C*     Courant Institute, Argonne National Lab, and Rice University
C*     March 31, 1993
C*
C*     .. Scalar Arguments ..
C      CHARACTER          TRANS
C      INTEGER            INFO, LDA, LDB, N, NRHS
C*     ..
C*     .. Array Arguments ..
C      INTEGER            IPIV( * )
C      REAL               A( LDA, * ), B( LDB, * )
C*     ..
C*
C*  Purpose
C*  =======
C*
C*  SGETRS solves a system of linear equations
C*     A * X = B  or  A' * X = B
C*  with a general N-by-N matrix A using the LU factorization computed
C*  by SGETRF.
C*
C*  Arguments
C*  =========
C*
C*  TRANS   (input) CHARACTER*1
C*          Specifies the form of the system of equations:
C*          = 'N':  A * X = B  (No transpose)
C*          = 'T':  A'* X = B  (Transpose)
C*          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
C*
C*  N       (input) INTEGER
C*          The order of the matrix A.  N >= 0.
C*
C*  NRHS    (input) INTEGER
C*          The number of right hand sides, i.e., the number of columns
C*          of the matrix B.  NRHS >= 0.
C*
C*  A       (input) REAL array, dimension (LDA,N)
C*          The factors L and U from the factorization A = P*L*U
C*          as computed by SGETRF.
C*
C*  LDA     (input) INTEGER
C*          The leading dimension of the array A.  LDA >= max(1,N).
C*
C*  IPIV    (input) INTEGER array, dimension (N)
C*          The pivot indices from SGETRF; for 1<=i<=N, row i of the
C*          matrix was interchanged with row IPIV(i).
C*
C*  B       (input/output) REAL array, dimension (LDB,NRHS)
C*          On entry, the right hand side matrix B.
C*          On exit, the solution matrix X.
C*
C*  LDB     (input) INTEGER
C*          The leading dimension of the array B.  LDB >= max(1,N).
C*
C*  INFO    (output) INTEGER
C*          = 0:  successful exit
C*          < 0:  if INFO = -i, the i-th argument had an illegal value
C*
C*  =====================================================================
C*
C*     .. Parameters ..
C      REAL               ONE
C      PARAMETER          ( ONE = 1.0E+0 )
C*     ..
C*     .. Local Scalars ..
C      LOGICAL            NOTRAN
C*     ..
C*     .. External Functions ..
C      LOGICAL            LSAME
C      EXTERNAL           LSAME
C*     ..
C*     .. External Subroutines ..
C      EXTERNAL           SLASWP, STRSM, XERBLA
C*     ..
C*     .. Intrinsic Functions ..
C      INTRINSIC          MAX
C*     ..
C*     .. Executable Statements ..
C*
C*     Test the input parameters.
C*
C      INFO = 0
C      NOTRAN = LSAME( TRANS, 'N' )
C      IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'T' ) .AND. .NOT.
C     $    LSAME( TRANS, 'C' ) ) THEN
C         INFO = -1
C      ELSE IF( N.LT.0 ) THEN
C         INFO = -2
C      ELSE IF( NRHS.LT.0 ) THEN
C         INFO = -3
C      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
C         INFO = -5
C      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
C         INFO = -8
C      END IF
C      IF( INFO.NE.0 ) THEN
C         CALL XERBLA( 'SGETRS', -INFO )
C         RETURN
C      END IF
C*
C*     Quick return if possible
C*
C      IF( N.EQ.0 .OR. NRHS.EQ.0 )
C     $   RETURN
C*
C      IF( NOTRAN ) THEN
C*
C*        Solve A * X = B.
C*
C*        Apply row interchanges to the right hand sides.
C*
C         CALL SLASWP( NRHS, B, LDB, 1, N, IPIV, 1 )
C*
C*        Solve L*X = B, overwriting B with X.
C*
C         CALL STRSM( 'Left', 'Lower', 'No transpose', 'Unit', N, NRHS,
C     $               ONE, A, LDA, B, LDB )
C*
C*        Solve U*X = B, overwriting B with X.
C*
C         CALL STRSM( 'Left', 'Upper', 'No transpose', 'Non-unit', N,
C     $               NRHS, ONE, A, LDA, B, LDB )
C      ELSE
C*
C*        Solve A' * X = B.
C*
C*        Solve U'*X = B, overwriting B with X.
C*
C         CALL STRSM( 'Left', 'Upper', 'Transpose', 'Non-unit', N, NRHS,
C     $               ONE, A, LDA, B, LDB )
C*
C*        Solve L'*X = B, overwriting B with X.
C*
C         CALL STRSM( 'Left', 'Lower', 'Transpose', 'Unit', N, NRHS, ONE,
C     $               A, LDA, B, LDB )
C*
C*        Apply row interchanges to the solution vectors.
C*
C         CALL SLASWP( NRHS, B, LDB, 1, N, IPIV, -1 )
C      END IF
C*
C      RETURN
C*
C*     End of SGETRS
C*
C      END
C      SUBROUTINE SLASWP( N, A, LDA, K1, K2, IPIV, INCX )
C*
C*  -- LAPACK auxiliary routine (version 1.1) --
C*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C*     Courant Institute, Argonne National Lab, and Rice University
C*     October 31, 1992
C*
C*     .. Scalar Arguments ..
C      INTEGER            INCX, K1, K2, LDA, N
C*     ..
C*     .. Array Arguments ..
C      INTEGER            IPIV( * )
C      REAL               A( LDA, * )
C*     ..
C*
C*  Purpose
C*  =======
C*
C*  SLASWP performs a series of row interchanges on the matrix A.
C*  One row interchange is initiated for each of rows K1 through K2 of A.
C*
C*  Arguments
C*  =========
C*
C*  N       (input) INTEGER
C*          The number of columns of the matrix A.
C*
C*  A       (input/output) REAL array, dimension (LDA,N)
C*          On entry, the matrix of column dimension N to which the row
C*          interchanges will be applied.
C*          On exit, the permuted matrix.
C*
C*  LDA     (input) INTEGER
C*          The leading dimension of the array A.
C*
C*  K1      (input) INTEGER
C*          The first element of IPIV for which a row interchange will
C*          be done.
C*
C*  K2      (input) INTEGER
C*          The last element of IPIV for which a row interchange will
C*          be done.
C*
C*  IPIV    (input) INTEGER array, dimension (M*abs(INCX))
C*          The vector of pivot indices.  Only the elements in positions
C*          K1 through K2 of IPIV are accessed.
C*          IPIV(K) = L implies rows K and L are to be interchanged.
C*
C*  INCX    (input) INTEGER
C*          The increment between successive values of IPIV.  If IPIV
C*          is negative, the pivots are applied in reverse order.
C*
C* =====================================================================
C*
C*     .. Local Scalars ..
C      INTEGER            I, IP, IX
C*     ..
C*     .. External Subroutines ..
C      EXTERNAL           SSWAP
C*     ..
C*     .. Executable Statements ..
C*
C*     Interchange row I with row IPIV(I) for each of rows K1 through K2.
C*
C      IF( INCX.EQ.0 )
C     $   RETURN
C      IF( INCX.GT.0 ) THEN
C         IX = K1
C      ELSE
C         IX = 1 + ( 1-K2 )*INCX
C      END IF
C      IF( INCX.EQ.1 ) THEN
C         DO 10 I = K1, K2
C            IP = IPIV( I )
C            IF( IP.NE.I )
C     $         CALL SSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
C   10    CONTINUE
C      ELSE IF( INCX.GT.1 ) THEN
C         DO 20 I = K1, K2
C            IP = IPIV( IX )
C            IF( IP.NE.I )
C     $         CALL SSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
C            IX = IX + INCX
C   20    CONTINUE
C      ELSE IF( INCX.LT.0 ) THEN
C         DO 30 I = K2, K1, -1
C            IP = IPIV( IX )
C            IF( IP.NE.I )
C     $         CALL SSWAP( N, A( I, 1 ), LDA, A( IP, 1 ), LDA )
C            IX = IX + INCX
C   30    CONTINUE
C      END IF
C*
C      RETURN
C*
C*     End of SLASWP
C*
C      END
C      SUBROUTINE STRMM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
C     $                   B, LDB )
C*     .. Scalar Arguments ..
C      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
C      INTEGER            M, N, LDA, LDB
C      REAL               ALPHA
C*     .. Array Arguments ..
C      REAL               A( LDA, * ), B( LDB, * )
C*     ..
C*
C*  Purpose
C*  =======
C*
C*  STRMM  performs one of the matrix-matrix operations
C*
C*     B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
C*
C*  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
C*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
C*
C*     op( A ) = A   or   op( A ) = A'.
C*
C*  Parameters
C*  ==========
C*
C*  SIDE   - CHARACTER*1.
C*           On entry,  SIDE specifies whether  op( A ) multiplies B from
C*           the left or right as follows:
C*
C*              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
C*
C*              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
C*
C*           Unchanged on exit.
C*
C*  UPLO   - CHARACTER*1.
C*           On entry, UPLO specifies whether the matrix A is an upper or
C*           lower triangular matrix as follows:
C*
C*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
C*
C*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
C*
C*           Unchanged on exit.
C*
C*  TRANSA - CHARACTER*1.
C*           On entry, TRANSA specifies the form of op( A ) to be used in
C*           the matrix multiplication as follows:
C*
C*              TRANSA = 'N' or 'n'   op( A ) = A.
C*
C*              TRANSA = 'T' or 't'   op( A ) = A'.
C*
C*              TRANSA = 'C' or 'c'   op( A ) = A'.
C*
C*           Unchanged on exit.
C*
C*  DIAG   - CHARACTER*1.
C*           On entry, DIAG specifies whether or not A is unit triangular
C*           as follows:
C*
C*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
C*
C*              DIAG = 'N' or 'n'   A is not assumed to be unit
C*                                  triangular.
C*
C*           Unchanged on exit.
C*
C*  M      - INTEGER.
C*           On entry, M specifies the number of rows of B. M must be at
C*           least zero.
C*           Unchanged on exit.
C*
C*  N      - INTEGER.
C*           On entry, N specifies the number of columns of B.  N must be
C*           at least zero.
C*           Unchanged on exit.
C*
C*  ALPHA  - REAL            .
C*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
C*           zero then  A is not referenced and  B need not be set before
C*           entry.
C*           Unchanged on exit.
C*
C*  A      - REAL             array of DIMENSION ( LDA, k ), where k is m
C*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
C*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
C*           upper triangular part of the array  A must contain the upper
C*           triangular matrix  and the strictly lower triangular part of
C*           A is not referenced.
C*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
C*           lower triangular part of the array  A must contain the lower
C*           triangular matrix  and the strictly upper triangular part of
C*           A is not referenced.
C*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
C*           A  are not referenced either,  but are assumed to be  unity.
C*           Unchanged on exit.
C*
C*  LDA    - INTEGER.
C*           On entry, LDA specifies the first dimension of A as declared
C*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
C*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
C*           then LDA must be at least max( 1, n ).
C*           Unchanged on exit.
C*
C*  B      - REAL             array of DIMENSION ( LDB, n ).
C*           Before entry,  the leading  m by n part of the array  B must
C*           contain the matrix  B,  and  on exit  is overwritten  by the
C*           transformed matrix.
C*
C*  LDB    - INTEGER.
C*           On entry, LDB specifies the first dimension of B as declared
C*           in  the  calling  (sub)  program.   LDB  must  be  at  least
C*           max( 1, m ).
C*           Unchanged on exit.
C*
C*
C*  Level 3 Blas routine.
C*
C*  -- Written on 8-February-1989.
C*     Jack Dongarra, Argonne National Laboratory.
C*     Iain Duff, AERE Harwell.
C*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
C*     Sven Hammarling, Numerical Algorithms Group Ltd.
C*
C*
C*     .. External Functions ..
C      LOGICAL            LSAME
C      EXTERNAL           LSAME
C*     .. External Subroutines ..
C      EXTERNAL           XERBLA
C*     .. Intrinsic Functions ..
C      INTRINSIC          MAX
C*     .. Local Scalars ..
C      LOGICAL            LSIDE, NOUNIT, UPPER
C      INTEGER            I, INFO, J, K, NROWA
C      REAL               TEMP
C*     .. Parameters ..
C      REAL               ONE         , ZERO
C      PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
C*     ..
C*     .. Executable Statements ..
C*
C*     Test the input parameters.
C*
C      LSIDE  = LSAME( SIDE  , 'L' )
C      IF( LSIDE )THEN
C         NROWA = M
C      ELSE
C         NROWA = N
C      END IF
C      NOUNIT = LSAME( DIAG  , 'N' )
C      UPPER  = LSAME( UPLO  , 'U' )
C*
C      INFO   = 0
C      IF(      ( .NOT.LSIDE                ).AND.
C     $         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
C         INFO = 1
C      ELSE IF( ( .NOT.UPPER                ).AND.
C     $         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
C         INFO = 2
C      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
C     $         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
C     $         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
C         INFO = 3
C      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
C     $         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
C         INFO = 4
C      ELSE IF( M  .LT.0               )THEN
C         INFO = 5
C      ELSE IF( N  .LT.0               )THEN
C         INFO = 6
C      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
C         INFO = 9
C      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
C         INFO = 11
C      END IF
C      IF( INFO.NE.0 )THEN
C         CALL XERBLA( 'STRMM ', INFO )
C         RETURN
C      END IF
C*
C*     Quick return if possible.
C*
C      IF( N.EQ.0 )
C     $   RETURN
C*
C*     And when  alpha.eq.zero.
C*
C      IF( ALPHA.EQ.ZERO )THEN
C         DO 20, J = 1, N
C            DO 10, I = 1, M
C               B( I, J ) = ZERO
C   10       CONTINUE
C   20    CONTINUE
C         RETURN
C      END IF
C*
C*     Start the operations.
C*
C      IF( LSIDE )THEN
C         IF( LSAME( TRANSA, 'N' ) )THEN
C*
C*           Form  B := alpha*A*B.
C*
C            IF( UPPER )THEN
C               DO 50, J = 1, N
C                  DO 40, K = 1, M
C                     IF( B( K, J ).NE.ZERO )THEN
C                        TEMP = ALPHA*B( K, J )
C                        DO 30, I = 1, K - 1
C                           B( I, J ) = B( I, J ) + TEMP*A( I, K )
C   30                   CONTINUE
C                        IF( NOUNIT )
C     $                     TEMP = TEMP*A( K, K )
C                        B( K, J ) = TEMP
C                     END IF
C   40             CONTINUE
C   50          CONTINUE
C            ELSE
C               DO 80, J = 1, N
C                  DO 70 K = M, 1, -1
C                     IF( B( K, J ).NE.ZERO )THEN
C                        TEMP      = ALPHA*B( K, J )
C                        B( K, J ) = TEMP
C                        IF( NOUNIT )
C     $                     B( K, J ) = B( K, J )*A( K, K )
C                        DO 60, I = K + 1, M
C                           B( I, J ) = B( I, J ) + TEMP*A( I, K )
C   60                   CONTINUE
C                     END IF
C   70             CONTINUE
C   80          CONTINUE
C            END IF
C         ELSE
C*
C*           Form  B := alpha*B*A'.
C*
C            IF( UPPER )THEN
C               DO 110, J = 1, N
C                  DO 100, I = M, 1, -1
C                     TEMP = B( I, J )
C                     IF( NOUNIT )
C     $                  TEMP = TEMP*A( I, I )
C                     DO 90, K = 1, I - 1
C                        TEMP = TEMP + A( K, I )*B( K, J )
C   90                CONTINUE
C                     B( I, J ) = ALPHA*TEMP
C  100             CONTINUE
C  110          CONTINUE
C            ELSE
C               DO 140, J = 1, N
C                  DO 130, I = 1, M
C                     TEMP = B( I, J )
C                     IF( NOUNIT )
C     $                  TEMP = TEMP*A( I, I )
C                     DO 120, K = I + 1, M
C                        TEMP = TEMP + A( K, I )*B( K, J )
C  120                CONTINUE
C                     B( I, J ) = ALPHA*TEMP
C  130             CONTINUE
C  140          CONTINUE
C            END IF
C         END IF
C      ELSE
C         IF( LSAME( TRANSA, 'N' ) )THEN
C*
C*           Form  B := alpha*B*A.
C*
C            IF( UPPER )THEN
C               DO 180, J = N, 1, -1
C                  TEMP = ALPHA
C                  IF( NOUNIT )
C     $               TEMP = TEMP*A( J, J )
C                  DO 150, I = 1, M
C                     B( I, J ) = TEMP*B( I, J )
C  150             CONTINUE
C                  DO 170, K = 1, J - 1
C                     IF( A( K, J ).NE.ZERO )THEN
C                        TEMP = ALPHA*A( K, J )
C                        DO 160, I = 1, M
C                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
C  160                   CONTINUE
C                     END IF
C  170             CONTINUE
C  180          CONTINUE
C            ELSE
C               DO 220, J = 1, N
C                  TEMP = ALPHA
C                  IF( NOUNIT )
C     $               TEMP = TEMP*A( J, J )
C                  DO 190, I = 1, M
C                     B( I, J ) = TEMP*B( I, J )
C  190             CONTINUE
C                  DO 210, K = J + 1, N
C                     IF( A( K, J ).NE.ZERO )THEN
C                        TEMP = ALPHA*A( K, J )
C                        DO 200, I = 1, M
C                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
C  200                   CONTINUE
C                     END IF
C  210             CONTINUE
C  220          CONTINUE
C            END IF
C         ELSE
C*
C*           Form  B := alpha*B*A'.
C*
C            IF( UPPER )THEN
C               DO 260, K = 1, N
C                  DO 240, J = 1, K - 1
C                     IF( A( J, K ).NE.ZERO )THEN
C                        TEMP = ALPHA*A( J, K )
C                        DO 230, I = 1, M
C                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
C  230                   CONTINUE
C                     END IF
C  240             CONTINUE
C                  TEMP = ALPHA
C                  IF( NOUNIT )
C     $               TEMP = TEMP*A( K, K )
C                  IF( TEMP.NE.ONE )THEN
C                     DO 250, I = 1, M
C                        B( I, K ) = TEMP*B( I, K )
C  250                CONTINUE
C                  END IF
C  260          CONTINUE
C            ELSE
C               DO 300, K = N, 1, -1
C                  DO 280, J = K + 1, N
C                     IF( A( J, K ).NE.ZERO )THEN
C                        TEMP = ALPHA*A( J, K )
C                        DO 270, I = 1, M
C                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
C  270                   CONTINUE
C                     END IF
C  280             CONTINUE
C                  TEMP = ALPHA
C                  IF( NOUNIT )
C     $               TEMP = TEMP*A( K, K )
C                  IF( TEMP.NE.ONE )THEN
C                     DO 290, I = 1, M
C                        B( I, K ) = TEMP*B( I, K )
C  290                CONTINUE
C                  END IF
C  300          CONTINUE
C            END IF
C         END IF
C      END IF
C*
C      RETURN
C*
C*     End of STRMM .
C*
C      END
C      SUBROUTINE STRMV ( UPLO, TRANS, DIAG, N, A, LDA, X, INCX )
C*     .. Scalar Arguments ..
C      INTEGER            INCX, LDA, N
C      CHARACTER*1        DIAG, TRANS, UPLO
C*     .. Array Arguments ..
C      REAL               A( LDA, * ), X( * )
C*     ..
C*
C*  Purpose
C*  =======
C*
C*  STRMV  performs one of the matrix-vector operations
C*
C*     x := A*x,   or   x := A'*x,
C*
C*  where x is an n element vector and  A is an n by n unit, or non-unit,
C*  upper or lower triangular matrix.
C*
C*  Parameters
C*  ==========
C*
C*  UPLO   - CHARACTER*1.
C*           On entry, UPLO specifies whether the matrix is an upper or
C*           lower triangular matrix as follows:
C*
C*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
C*
C*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
C*
C*           Unchanged on exit.
C*
C*  TRANS  - CHARACTER*1.
C*           On entry, TRANS specifies the operation to be performed as
C*           follows:
C*
C*              TRANS = 'N' or 'n'   x := A*x.
C*
C*              TRANS = 'T' or 't'   x := A'*x.
C*
C*              TRANS = 'C' or 'c'   x := A'*x.
C*
C*           Unchanged on exit.
C*
C*  DIAG   - CHARACTER*1.
C*           On entry, DIAG specifies whether or not A is unit
C*           triangular as follows:
C*
C*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
C*
C*              DIAG = 'N' or 'n'   A is not assumed to be unit
C*                                  triangular.
C*
C*           Unchanged on exit.
C*
C*  N      - INTEGER.
C*           On entry, N specifies the order of the matrix A.
C*           N must be at least zero.
C*           Unchanged on exit.
C*
C*  A      - REAL             array of DIMENSION ( LDA, n ).
C*           Before entry with  UPLO = 'U' or 'u', the leading n by n
C*           upper triangular part of the array A must contain the upper
C*           triangular matrix and the strictly lower triangular part of
C*           A is not referenced.
C*           Before entry with UPLO = 'L' or 'l', the leading n by n
C*           lower triangular part of the array A must contain the lower
C*           triangular matrix and the strictly upper triangular part of
C*           A is not referenced.
C*           Note that when  DIAG = 'U' or 'u', the diagonal elements of
C*           A are not referenced either, but are assumed to be unity.
C*           Unchanged on exit.
C*
C*  LDA    - INTEGER.
C*           On entry, LDA specifies the first dimension of A as declared
C*           in the calling (sub) program. LDA must be at least
C*           max( 1, n ).
C*           Unchanged on exit.
C*
C*  X      - REAL             array of dimension at least
C*           ( 1 + ( n - 1 )*abs( INCX ) ).
C*           Before entry, the incremented array X must contain the n
C*           element vector x. On exit, X is overwritten with the
C*           tranformed vector x.
C*
C*  INCX   - INTEGER.
C*           On entry, INCX specifies the increment for the elements of
C*           X. INCX must not be zero.
C*           Unchanged on exit.
C*
C*
C*  Level 2 Blas routine.
C*
C*  -- Written on 22-October-1986.
C*     Jack Dongarra, Argonne National Lab.
C*     Jeremy Du Croz, Nag Central Office.
C*     Sven Hammarling, Nag Central Office.
C*     Richard Hanson, Sandia National Labs.
C*
C*
C*     .. Parameters ..
C      REAL               ZERO
C      PARAMETER        ( ZERO = 0.0E+0 )
C*     .. Local Scalars ..
C      REAL               TEMP
C      INTEGER            I, INFO, IX, J, JX, KX
C      LOGICAL            NOUNIT
C*     .. External Functions ..
C      LOGICAL            LSAME
C      EXTERNAL           LSAME
C*     .. External Subroutines ..
C      EXTERNAL           XERBLA
C*     .. Intrinsic Functions ..
C      INTRINSIC          MAX
C*     ..
C*     .. Executable Statements ..
C*
C*     Test the input parameters.
C*
C      INFO = 0
C      IF     ( .NOT.LSAME( UPLO , 'U' ).AND.
C     $         .NOT.LSAME( UPLO , 'L' )      )THEN
C         INFO = 1
C      ELSE IF( .NOT.LSAME( TRANS, 'N' ).AND.
C     $         .NOT.LSAME( TRANS, 'T' ).AND.
C     $         .NOT.LSAME( TRANS, 'C' )      )THEN
C         INFO = 2
C      ELSE IF( .NOT.LSAME( DIAG , 'U' ).AND.
C     $         .NOT.LSAME( DIAG , 'N' )      )THEN
C         INFO = 3
C      ELSE IF( N.LT.0 )THEN
C         INFO = 4
C      ELSE IF( LDA.LT.MAX( 1, N ) )THEN
C         INFO = 6
C      ELSE IF( INCX.EQ.0 )THEN
C         INFO = 8
C      END IF
C      IF( INFO.NE.0 )THEN
C         CALL XERBLA( 'STRMV ', INFO )
C         RETURN
C      END IF
C*
C*     Quick return if possible.
C*
C      IF( N.EQ.0 )
C     $   RETURN
C*
C      NOUNIT = LSAME( DIAG, 'N' )
C*
C*     Set up the start point in X if the increment is not unity. This
C*     will be  ( N - 1 )*INCX  too small for descending loops.
C*
C      IF( INCX.LE.0 )THEN
C         KX = 1 - ( N - 1 )*INCX
C      ELSE IF( INCX.NE.1 )THEN
C         KX = 1
C      END IF
C*
C*     Start the operations. In this version the elements of A are
C*     accessed sequentially with one pass through A.
C*
C      IF( LSAME( TRANS, 'N' ) )THEN
C*
C*        Form  x := A*x.
C*
C         IF( LSAME( UPLO, 'U' ) )THEN
C            IF( INCX.EQ.1 )THEN
C               DO 20, J = 1, N
C                  IF( X( J ).NE.ZERO )THEN
C                     TEMP = X( J )
C                     DO 10, I = 1, J - 1
C                        X( I ) = X( I ) + TEMP*A( I, J )
C   10                CONTINUE
C                     IF( NOUNIT )
C     $                  X( J ) = X( J )*A( J, J )
C                  END IF
C   20          CONTINUE
C            ELSE
C               JX = KX
C               DO 40, J = 1, N
C                  IF( X( JX ).NE.ZERO )THEN
C                     TEMP = X( JX )
C                     IX   = KX
C                     DO 30, I = 1, J - 1
C                        X( IX ) = X( IX ) + TEMP*A( I, J )
C                        IX      = IX      + INCX
C   30                CONTINUE
C                     IF( NOUNIT )
C     $                  X( JX ) = X( JX )*A( J, J )
C                  END IF
C                  JX = JX + INCX
C   40          CONTINUE
C            END IF
C         ELSE
C            IF( INCX.EQ.1 )THEN
C               DO 60, J = N, 1, -1
C                  IF( X( J ).NE.ZERO )THEN
C                     TEMP = X( J )
C                     DO 50, I = N, J + 1, -1
C                        X( I ) = X( I ) + TEMP*A( I, J )
C   50                CONTINUE
C                     IF( NOUNIT )
C     $                  X( J ) = X( J )*A( J, J )
C                  END IF
C   60          CONTINUE
C            ELSE
C               KX = KX + ( N - 1 )*INCX
C               JX = KX
C               DO 80, J = N, 1, -1
C                  IF( X( JX ).NE.ZERO )THEN
C                     TEMP = X( JX )
C                     IX   = KX
C                     DO 70, I = N, J + 1, -1
C                        X( IX ) = X( IX ) + TEMP*A( I, J )
C                        IX      = IX      - INCX
C   70                CONTINUE
C                     IF( NOUNIT )
C     $                  X( JX ) = X( JX )*A( J, J )
C                  END IF
C                  JX = JX - INCX
C   80          CONTINUE
C            END IF
C         END IF
C      ELSE
C*
C*        Form  x := A'*x.
C*
C         IF( LSAME( UPLO, 'U' ) )THEN
C            IF( INCX.EQ.1 )THEN
C               DO 100, J = N, 1, -1
C                  TEMP = X( J )
C                  IF( NOUNIT )
C     $               TEMP = TEMP*A( J, J )
C                  DO 90, I = J - 1, 1, -1
C                     TEMP = TEMP + A( I, J )*X( I )
C   90             CONTINUE
C                  X( J ) = TEMP
C  100          CONTINUE
C            ELSE
C               JX = KX + ( N - 1 )*INCX
C               DO 120, J = N, 1, -1
C                  TEMP = X( JX )
C                  IX   = JX
C                  IF( NOUNIT )
C     $               TEMP = TEMP*A( J, J )
C                  DO 110, I = J - 1, 1, -1
C                     IX   = IX   - INCX
C                     TEMP = TEMP + A( I, J )*X( IX )
C  110             CONTINUE
C                  X( JX ) = TEMP
C                  JX      = JX   - INCX
C  120          CONTINUE
C            END IF
C         ELSE
C            IF( INCX.EQ.1 )THEN
C               DO 140, J = 1, N
C                  TEMP = X( J )
C                  IF( NOUNIT )
C     $               TEMP = TEMP*A( J, J )
C                  DO 130, I = J + 1, N
C                     TEMP = TEMP + A( I, J )*X( I )
C  130             CONTINUE
C                  X( J ) = TEMP
C  140          CONTINUE
C            ELSE
C               JX = KX
C               DO 160, J = 1, N
C                  TEMP = X( JX )
C                  IX   = JX
C                  IF( NOUNIT )
C     $               TEMP = TEMP*A( J, J )
C                  DO 150, I = J + 1, N
C                     IX   = IX   + INCX
C                     TEMP = TEMP + A( I, J )*X( IX )
C  150             CONTINUE
C                  X( JX ) = TEMP
C                  JX      = JX   + INCX
C  160          CONTINUE
C            END IF
C         END IF
C      END IF
C*
C      RETURN
C*
C*     End of STRMV .
C*
C      END
C      SUBROUTINE STRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
C     $                   B, LDB )
C*     .. Scalar Arguments ..
C      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
C      INTEGER            M, N, LDA, LDB
C      REAL               ALPHA
C*     .. Array Arguments ..
C      REAL               A( LDA, * ), B( LDB, * )
C*     ..
C*
C*  Purpose
C*  =======
C*
C*  STRSM  solves one of the matrix equations
C*
C*     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
C*
C*  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
C*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
C*
C*     op( A ) = A   or   op( A ) = A'.
C*
C*  The matrix X is overwritten on B.
C*
C*  Parameters
C*  ==========
C*
C*  SIDE   - CHARACTER*1.
C*           On entry, SIDE specifies whether op( A ) appears on the left
C*           or right of X as follows:
C*
C*              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
C*
C*              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
C*
C*           Unchanged on exit.
C*
C*  UPLO   - CHARACTER*1.
C*           On entry, UPLO specifies whether the matrix A is an upper or
C*           lower triangular matrix as follows:
C*
C*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
C*
C*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
C*
C*           Unchanged on exit.
C*
C*  TRANSA - CHARACTER*1.
C*           On entry, TRANSA specifies the form of op( A ) to be used in
C*           the matrix multiplication as follows:
C*
C*              TRANSA = 'N' or 'n'   op( A ) = A.
C*
C*              TRANSA = 'T' or 't'   op( A ) = A'.
C*
C*              TRANSA = 'C' or 'c'   op( A ) = A'.
C*
C*           Unchanged on exit.
C*
C*  DIAG   - CHARACTER*1.
C*           On entry, DIAG specifies whether or not A is unit triangular
C*           as follows:
C*
C*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
C*
C*              DIAG = 'N' or 'n'   A is not assumed to be unit
C*                                  triangular.
C*
C*           Unchanged on exit.
C*
C*  M      - INTEGER.
C*           On entry, M specifies the number of rows of B. M must be at
C*           least zero.
C*           Unchanged on exit.
C*
C*  N      - INTEGER.
C*           On entry, N specifies the number of columns of B.  N must be
C*           at least zero.
C*           Unchanged on exit.
C*
C*  ALPHA  - REAL            .
C*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
C*           zero then  A is not referenced and  B need not be set before
C*           entry.
C*           Unchanged on exit.
C*
C*  A      - REAL             array of DIMENSION ( LDA, k ), where k is m
C*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
C*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
C*           upper triangular part of the array  A must contain the upper
C*           triangular matrix  and the strictly lower triangular part of
C*           A is not referenced.
C*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
C*           lower triangular part of the array  A must contain the lower
C*           triangular matrix  and the strictly upper triangular part of
C*           A is not referenced.
C*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
C*           A  are not referenced either,  but are assumed to be  unity.
C*           Unchanged on exit.
C*
C*  LDA    - INTEGER.
C*           On entry, LDA specifies the first dimension of A as declared
C*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
C*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
C*           then LDA must be at least max( 1, n ).
C*           Unchanged on exit.
C*
C*  B      - REAL             array of DIMENSION ( LDB, n ).
C*           Before entry,  the leading  m by n part of the array  B must
C*           contain  the  right-hand  side  matrix  B,  and  on exit  is
C*           overwritten by the solution matrix  X.
C*
C*  LDB    - INTEGER.
C*           On entry, LDB specifies the first dimension of B as declared
C*           in  the  calling  (sub)  program.   LDB  must  be  at  least
C*           max( 1, m ).
C*           Unchanged on exit.
C*
C*
C*  Level 3 Blas routine.
C*
C*
C*  -- Written on 8-February-1989.
C*     Jack Dongarra, Argonne National Laboratory.
C*     Iain Duff, AERE Harwell.
C*     Jeremy Du Croz, Numerical Algorithms Group Ltd.
C*     Sven Hammarling, Numerical Algorithms Group Ltd.
C*
C*
C*     .. External Functions ..
C      LOGICAL            LSAME
C      EXTERNAL           LSAME
C*     .. External Subroutines ..
C      EXTERNAL           XERBLA
C*     .. Intrinsic Functions ..
C      INTRINSIC          MAX
C*     .. Local Scalars ..
C      LOGICAL            LSIDE, NOUNIT, UPPER
C      INTEGER            I, INFO, J, K, NROWA
C      REAL               TEMP
C*     .. Parameters ..
C      REAL               ONE         , ZERO
C      PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
C*     ..
C*     .. Executable Statements ..
C*
C*     Test the input parameters.
C*
C      LSIDE  = LSAME( SIDE  , 'L' )
C      IF( LSIDE )THEN
C         NROWA = M
C      ELSE
C         NROWA = N
C      END IF
C      NOUNIT = LSAME( DIAG  , 'N' )
C      UPPER  = LSAME( UPLO  , 'U' )
C*
C      INFO   = 0
C      IF(      ( .NOT.LSIDE                ).AND.
C     $         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
C         INFO = 1
C      ELSE IF( ( .NOT.UPPER                ).AND.
C     $         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
C         INFO = 2
C      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
C     $         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
C     $         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
C         INFO = 3
C      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
C     $         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
C         INFO = 4
C      ELSE IF( M  .LT.0               )THEN
C         INFO = 5
C      ELSE IF( N  .LT.0               )THEN
C         INFO = 6
C      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
C         INFO = 9
C      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
C         INFO = 11
C      END IF
C      IF( INFO.NE.0 )THEN
C         CALL XERBLA( 'STRSM ', INFO )
C         RETURN
C      END IF
C*
C*     Quick return if possible.
C*
C      IF( N.EQ.0 )
C     $   RETURN
C*
C*     And when  alpha.eq.zero.
C*
C      IF( ALPHA.EQ.ZERO )THEN
C         DO 20, J = 1, N
C            DO 10, I = 1, M
C               B( I, J ) = ZERO
C   10       CONTINUE
C   20    CONTINUE
C         RETURN
C      END IF
C*
C*     Start the operations.
C*
C      IF( LSIDE )THEN
C         IF( LSAME( TRANSA, 'N' ) )THEN
C*
C*           Form  B := alpha*inv( A )*B.
C*
C            IF( UPPER )THEN
C               DO 60, J = 1, N
C                  IF( ALPHA.NE.ONE )THEN
C                     DO 30, I = 1, M
C                        B( I, J ) = ALPHA*B( I, J )
C   30                CONTINUE
C                  END IF
C                  DO 50, K = M, 1, -1
C                     IF( B( K, J ).NE.ZERO )THEN
C                        IF( NOUNIT )
C     $                     B( K, J ) = B( K, J )/A( K, K )
C                        DO 40, I = 1, K - 1
C                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
C   40                   CONTINUE
C                     END IF
C   50             CONTINUE
C   60          CONTINUE
C            ELSE
C               DO 100, J = 1, N
C                  IF( ALPHA.NE.ONE )THEN
C                     DO 70, I = 1, M
C                        B( I, J ) = ALPHA*B( I, J )
C   70                CONTINUE
C                  END IF
C                  DO 90 K = 1, M
C                     IF( B( K, J ).NE.ZERO )THEN
C                        IF( NOUNIT )
C     $                     B( K, J ) = B( K, J )/A( K, K )
C                        DO 80, I = K + 1, M
C                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
C   80                   CONTINUE
C                     END IF
C   90             CONTINUE
C  100          CONTINUE
C            END IF
C         ELSE
C*
C*           Form  B := alpha*inv( A' )*B.
C*
C            IF( UPPER )THEN
C               DO 130, J = 1, N
C                  DO 120, I = 1, M
C                     TEMP = ALPHA*B( I, J )
C                     DO 110, K = 1, I - 1
C                        TEMP = TEMP - A( K, I )*B( K, J )
C  110                CONTINUE
C                     IF( NOUNIT )
C     $                  TEMP = TEMP/A( I, I )
C                     B( I, J ) = TEMP
C  120             CONTINUE
C  130          CONTINUE
C            ELSE
C               DO 160, J = 1, N
C                  DO 150, I = M, 1, -1
C                     TEMP = ALPHA*B( I, J )
C                     DO 140, K = I + 1, M
C                        TEMP = TEMP - A( K, I )*B( K, J )
C  140                CONTINUE
C                     IF( NOUNIT )
C     $                  TEMP = TEMP/A( I, I )
C                     B( I, J ) = TEMP
C  150             CONTINUE
C  160          CONTINUE
C            END IF
C         END IF
C      ELSE
C         IF( LSAME( TRANSA, 'N' ) )THEN
C*
C*           Form  B := alpha*B*inv( A ).
C*
C            IF( UPPER )THEN
C               DO 210, J = 1, N
C                  IF( ALPHA.NE.ONE )THEN
C                     DO 170, I = 1, M
C                        B( I, J ) = ALPHA*B( I, J )
C  170                CONTINUE
C                  END IF
C                  DO 190, K = 1, J - 1
C                     IF( A( K, J ).NE.ZERO )THEN
C                        DO 180, I = 1, M
C                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
C  180                   CONTINUE
C                     END IF
C  190             CONTINUE
C                  IF( NOUNIT )THEN
C                     TEMP = ONE/A( J, J )
C                     DO 200, I = 1, M
C                        B( I, J ) = TEMP*B( I, J )
C  200                CONTINUE
C                  END IF
C  210          CONTINUE
C            ELSE
C               DO 260, J = N, 1, -1
C                  IF( ALPHA.NE.ONE )THEN
C                     DO 220, I = 1, M
C                        B( I, J ) = ALPHA*B( I, J )
C  220                CONTINUE
C                  END IF
C                  DO 240, K = J + 1, N
C                     IF( A( K, J ).NE.ZERO )THEN
C                        DO 230, I = 1, M
C                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
C  230                   CONTINUE
C                     END IF
C  240             CONTINUE
C                  IF( NOUNIT )THEN
C                     TEMP = ONE/A( J, J )
C                     DO 250, I = 1, M
C                       B( I, J ) = TEMP*B( I, J )
C  250                CONTINUE
C                  END IF
C  260          CONTINUE
C            END IF
C         ELSE
C*
C*           Form  B := alpha*B*inv( A' ).
C*
C            IF( UPPER )THEN
C               DO 310, K = N, 1, -1
C                  IF( NOUNIT )THEN
C                     TEMP = ONE/A( K, K )
C                     DO 270, I = 1, M
C                        B( I, K ) = TEMP*B( I, K )
C  270                CONTINUE
C                  END IF
C                  DO 290, J = 1, K - 1
C                     IF( A( J, K ).NE.ZERO )THEN
C                        TEMP = A( J, K )
C                        DO 280, I = 1, M
C                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
C  280                   CONTINUE
C                     END IF
C  290             CONTINUE
C                  IF( ALPHA.NE.ONE )THEN
C                     DO 300, I = 1, M
C                        B( I, K ) = ALPHA*B( I, K )
C  300                CONTINUE
C                  END IF
C  310          CONTINUE
C            ELSE
C               DO 360, K = 1, N
C                  IF( NOUNIT )THEN
C                     TEMP = ONE/A( K, K )
C                     DO 320, I = 1, M
C                        B( I, K ) = TEMP*B( I, K )
C  320                CONTINUE
C                  END IF
C                  DO 340, J = K + 1, N
C                     IF( A( J, K ).NE.ZERO )THEN
C                        TEMP = A( J, K )
C                        DO 330, I = 1, M
C                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
C  330                   CONTINUE
C                     END IF
C  340             CONTINUE
C                  IF( ALPHA.NE.ONE )THEN
C                     DO 350, I = 1, M
C                        B( I, K ) = ALPHA*B( I, K )
C  350                CONTINUE
C                  END IF
C  360          CONTINUE
C            END IF
C         END IF
C      END IF
C*
C      RETURN
C*
C*     End of STRSM .
C*
C      END
C      SUBROUTINE STRTI2( UPLO, DIAG, N, A, LDA, INFO )
C*
C*  -- LAPACK routine (version 1.1) --
C*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C*     Courant Institute, Argonne National Lab, and Rice University
C*     February 29, 1992
C*
C*     .. Scalar Arguments ..
C      CHARACTER          DIAG, UPLO
C      INTEGER            INFO, LDA, N
C*     ..
C*     .. Array Arguments ..
C      REAL               A( LDA, * )
C*     ..
C*
C*  Purpose
C*  =======
C*
C*  STRTI2 computes the inverse of a real upper or lower triangular
C*  matrix.
C*
C*  This is the Level 2 BLAS version of the algorithm.
C*
C*  Arguments
C*  =========
C*
C*  UPLO    (input) CHARACTER*1
C*          Specifies whether the matrix A is upper or lower triangular.
C*          = 'U':  Upper triangular
C*          = 'L':  Lower triangular
C*
C*  DIAG    (input) CHARACTER*1
C*          Specifies whether or not the matrix A is unit triangular.
C*          = 'N':  Non-unit triangular
C*          = 'U':  Unit triangular
C*
C*  N       (input) INTEGER
C*          The order of the matrix A.  N >= 0.
C*
C*  A       (input/output) REAL array, dimension (LDA,N)
C*          On entry, the triangular matrix A.  If UPLO = 'U', the
C*          leading n by n upper triangular part of the array A contains
C*          the upper triangular matrix, and the strictly lower
C*          triangular part of A is not referenced.  If UPLO = 'L', the
C*          leading n by n lower triangular part of the array A contains
C*          the lower triangular matrix, and the strictly upper
C*          triangular part of A is not referenced.  If DIAG = 'U', the
C*          diagonal elements of A are also not referenced and are
C*          assumed to be 1.
C*
C*          On exit, the (triangular) inverse of the original matrix, in
C*          the same storage format.
C*
C*  LDA     (input) INTEGER
C*          The leading dimension of the array A.  LDA >= max(1,N).
C*
C*  INFO    (output) INTEGER
C*          = 0: successful exit
C*          < 0: if INFO = -k, the k-th argument had an illegal value
C*
C*  =====================================================================
C*
C*     .. Parameters ..
C      REAL               ONE
C      PARAMETER          ( ONE = 1.0E+0 )
C*     ..
C*     .. Local Scalars ..
C      LOGICAL            NOUNIT, UPPER
C      INTEGER            J
C      REAL               AJJ
C*     ..
C*     .. External Functions ..
C      LOGICAL            LSAME
C      EXTERNAL           LSAME
C*     ..
C*     .. External Subroutines ..
C      EXTERNAL           SSCAL, STRMV, XERBLA
C*     ..
C*     .. Intrinsic Functions ..
C      INTRINSIC          MAX
C*     ..
C*     .. Executable Statements ..
C*
C*     Test the input parameters.
C*
C      INFO = 0
C      UPPER = LSAME( UPLO, 'U' )
C      NOUNIT = LSAME( DIAG, 'N' )
C      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
C         INFO = -1
C      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
C         INFO = -2
C      ELSE IF( N.LT.0 ) THEN
C         INFO = -3
C      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
C         INFO = -5
C      END IF
C      IF( INFO.NE.0 ) THEN
C         CALL XERBLA( 'STRTI2', -INFO )
C         RETURN
C      END IF
C*
C      IF( UPPER ) THEN
C*
C*        Compute inverse of upper triangular matrix.
C*
C         DO 10 J = 1, N
C            IF( NOUNIT ) THEN
C               A( J, J ) = ONE / A( J, J )
C               AJJ = -A( J, J )
C            ELSE
C               AJJ = -ONE
C            END IF
C*
C*           Compute elements 1:j-1 of j-th column.
C*
C            CALL STRMV( 'Upper', 'No transpose', DIAG, J-1, A, LDA,
C     $                  A( 1, J ), 1 )
C            CALL SSCAL( J-1, AJJ, A( 1, J ), 1 )
C   10    CONTINUE
C      ELSE
C*
C*        Compute inverse of lower triangular matrix.
C*
C         DO 20 J = N, 1, -1
C            IF( NOUNIT ) THEN
C               A( J, J ) = ONE / A( J, J )
C               AJJ = -A( J, J )
C            ELSE
C               AJJ = -ONE
C            END IF
C            IF( J.LT.N ) THEN
C*
C*              Compute elements j+1:n of j-th column.
C*
C               CALL STRMV( 'Lower', 'No transpose', DIAG, N-J,
C     $                     A( J+1, J+1 ), LDA, A( J+1, J ), 1 )
C               CALL SSCAL( N-J, AJJ, A( J+1, J ), 1 )
C            END IF
C   20    CONTINUE
C      END IF
C*
C      RETURN
C*
C*     End of STRTI2
C*
C      END
C      SUBROUTINE STRTRI( UPLO, DIAG, N, A, LDA, INFO )
C*
C*  -- LAPACK routine (version 1.1) --
C*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
C*     Courant Institute, Argonne National Lab, and Rice University
C*     March 31, 1993
C*
C*     .. Scalar Arguments ..
C      CHARACTER          DIAG, UPLO
C      INTEGER            INFO, LDA, N
C*     ..
C*     .. Array Arguments ..
C      REAL               A( LDA, * )
C*     ..
C*
C*  Purpose
C*  =======
C*
C*  STRTRI computes the inverse of a real upper or lower triangular
C*  matrix A.
C*
C*  This is the Level 3 BLAS version of the algorithm.
C*
C*  Arguments
C*  =========
C*
C*  UPLO    (input) CHARACTER*1
C*          = 'U':  A is upper triangular;
C*          = 'L':  A is lower triangular.
C*
C*  DIAG    (input) CHARACTER*1
C*          = 'N':  A is non-unit triangular;
C*          = 'U':  A is unit triangular.
C*
C*  N       (input) INTEGER
C*          The order of the matrix A.  N >= 0.
C*
C*  A       (input/output) REAL array, dimension (LDA,N)
C*          On entry, the triangular matrix A.  If UPLO = 'U', the
C*          leading N-by-N upper triangular part of the array A contains
C*          the upper triangular matrix, and the strictly lower
C*          triangular part of A is not referenced.  If UPLO = 'L', the
C*          leading N-by-N lower triangular part of the array A contains
C*          the lower triangular matrix, and the strictly upper
C*          triangular part of A is not referenced.  If DIAG = 'U', the
C*          diagonal elements of A are also not referenced and are
C*          assumed to be 1.
C*          On exit, the (triangular) inverse of the original matrix, in
C*          the same storage format.
C*
C*  LDA     (input) INTEGER
C*          The leading dimension of the array A.  LDA >= max(1,N).
C*
C*  INFO    (output) INTEGER
C*          = 0: successful exit
C*          < 0: if INFO = -i, the i-th argument had an illegal value
C*          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
C*               matrix is singular and its inverse can not be computed.
C*
C*  =====================================================================
C*
C*     .. Parameters ..
C      REAL               ONE, ZERO
C      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
C*     ..
C*     .. Local Scalars ..
C      LOGICAL            NOUNIT, UPPER
C      INTEGER            J, JB, NB, NN
C*     ..
C*     .. External Functions ..
C      LOGICAL            LSAME
C      INTEGER            ILAENV
C      EXTERNAL           LSAME, ILAENV
C*     ..
C*     .. External Subroutines ..
C      EXTERNAL           STRMM, STRSM, STRTI2, XERBLA
C*     ..
C*     .. Intrinsic Functions ..
C      INTRINSIC          MAX, MIN
C*     ..
C*     .. Executable Statements ..
C*
C*     Test the input parameters.
C*
C      INFO = 0
C      UPPER = LSAME( UPLO, 'U' )
C      NOUNIT = LSAME( DIAG, 'N' )
C      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
C         INFO = -1
C      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
C         INFO = -2
C      ELSE IF( N.LT.0 ) THEN
C         INFO = -3
C      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
C         INFO = -5
C      END IF
C      IF( INFO.NE.0 ) THEN
C         CALL XERBLA( 'STRTRI', -INFO )
C         RETURN
C      END IF
C*
C*     Quick return if possible
C*
C      IF( N.EQ.0 )
C     $   RETURN
C*
C*     Check for singularity if non-unit.
C*
C      IF( NOUNIT ) THEN
C         DO 10 INFO = 1, N
C            IF( A( INFO, INFO ).EQ.ZERO )
C     $         RETURN
C   10    CONTINUE
C         INFO = 0
C      END IF
C*
C*     Determine the block size for this environment.
C*
C      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
C      IF( NB.LE.1 .OR. NB.GE.N ) THEN
C*
C*        Use unblocked code
C*
C         CALL STRTI2( UPLO, DIAG, N, A, LDA, INFO )
C      ELSE
C*
C*        Use blocked code
C*
C         IF( UPPER ) THEN
C*
C*           Compute inverse of upper triangular matrix
C*
C            DO 20 J = 1, N, NB
C               JB = MIN( NB, N-J+1 )
C*
C*              Compute rows 1:j-1 of current block column
C*
C               CALL STRMM( 'Left', 'Upper', 'No transpose', DIAG, J-1,
C     $                     JB, ONE, A, LDA, A( 1, J ), LDA )
C               CALL STRSM( 'Right', 'Upper', 'No transpose', DIAG, J-1,
C     $                     JB, -ONE, A( J, J ), LDA, A( 1, J ), LDA )
C*
C*              Compute inverse of current diagonal block
C*
C               CALL STRTI2( 'Upper', DIAG, JB, A( J, J ), LDA, INFO )
C   20       CONTINUE
C         ELSE
C*
C*           Compute inverse of lower triangular matrix
C*
C            NN = ( ( N-1 ) / NB )*NB + 1
C            DO 30 J = NN, 1, -NB
C               JB = MIN( NB, N-J+1 )
C               IF( J+JB.LE.N ) THEN
C*
C*                 Compute rows j+jb:n of current block column
C*
C                  CALL STRMM( 'Left', 'Lower', 'No transpose', DIAG,
C     $                        N-J-JB+1, JB, ONE, A( J+JB, J+JB ), LDA,
C     $                        A( J+JB, J ), LDA )
C                  CALL STRSM( 'Right', 'Lower', 'No transpose', DIAG,
C     $                        N-J-JB+1, JB, -ONE, A( J, J ), LDA,
C     $                        A( J+JB, J ), LDA )
C               END IF
C*
C*              Compute inverse of current diagonal block
C*
C               CALL STRTI2( 'Lower', DIAG, JB, A( J, J ), LDA, INFO )
C   30       CONTINUE
C         END IF
C      END IF
C*
C      RETURN
C*
C*     End of STRTRI
C*
C      END
C      real function sasum(n,sx,incx)
Cc
Cc     takes the sum of the absolute values.
Cc     uses unrolled loops for increment equal to one.
Cc     jack dongarra, linpack, 3/11/78.
Cc     modified 3/93 to return if incx .le. 0.
Cc
C      real sx(1),stemp
C      integer i,incx,m,mp1,n,nincx
Cc
C      sasum = 0.0e0
C      stemp = 0.0e0
C      if( n.le.0 .or. incx.le.0 )return
C      if(incx.eq.1)go to 20
Cc
Cc        code for increment not equal to 1
Cc
C      nincx = n*incx
C      do 10 i = 1,nincx,incx
C        stemp = stemp + abs(sx(i))
C   10 continue
C      sasum = stemp
C      return
Cc
Cc        code for increment equal to 1
Cc
Cc
Cc        clean-up loop
Cc
C   20 m = mod(n,6)
C      if( m .eq. 0 ) go to 40
C      do 30 i = 1,m
C        stemp = stemp + abs(sx(i))
C   30 continue
C      if( n .lt. 6 ) go to 60
C   40 mp1 = m + 1
C      do 50 i = mp1,n,6
C        stemp = stemp + abs(sx(i)) + abs(sx(i + 1)) + abs(sx(i + 2))
C     *  + abs(sx(i + 3)) + abs(sx(i + 4)) + abs(sx(i + 5))
C   50 continue
C   60 sasum = stemp
C      return
C      end
C      subroutine saxpy(n,sa,sx,incx,sy,incy)
Cc
Cc     constant times a vector plus a vector.
Cc     uses unrolled loop for increments equal to one.
Cc     jack dongarra, linpack, 3/11/78.
Cc
C      real sx(1),sy(1),sa
C      integer i,incx,incy,ix,iy,m,mp1,n
Cc
C      if(n.le.0)return
C      if (sa .eq. 0.0) return
C      if(incx.eq.1.and.incy.eq.1)go to 20
Cc
Cc        code for unequal increments or equal increments
Cc          not equal to 1
Cc
C      ix = 1
C      iy = 1
C      if(incx.lt.0)ix = (-n+1)*incx + 1
C      if(incy.lt.0)iy = (-n+1)*incy + 1
C      do 10 i = 1,n
C        sy(iy) = sy(iy) + sa*sx(ix)
C        ix = ix + incx
C        iy = iy + incy
C   10 continue
C      return
Cc
Cc        code for both increments equal to 1
Cc
Cc
Cc        clean-up loop
Cc
C   20 m = mod(n,4)
C      if( m .eq. 0 ) go to 40
C      do 30 i = 1,m
C        sy(i) = sy(i) + sa*sx(i)
C   30 continue
C      if( n .lt. 4 ) return
C   40 mp1 = m + 1
C      do 50 i = mp1,n,4
C        sy(i) = sy(i) + sa*sx(i)
C        sy(i + 1) = sy(i + 1) + sa*sx(i + 1)
C        sy(i + 2) = sy(i + 2) + sa*sx(i + 2)
C        sy(i + 3) = sy(i + 3) + sa*sx(i + 3)
C   50 continue
C      return
C      end
C      subroutine scopy(n,sx,incx,sy,incy)
Cc
Cc     copies a vector, x, to a vector, y.
Cc     uses unrolled loops for increments equal to 1.
Cc     jack dongarra, linpack, 3/11/78.
Cc
C      real sx(1),sy(1)
C      integer i,incx,incy,ix,iy,m,mp1,n
Cc
C      if(n.le.0)return
C      if(incx.eq.1.and.incy.eq.1)go to 20
Cc
Cc        code for unequal increments or equal increments
Cc          not equal to 1
Cc
C      ix = 1
C      iy = 1
C      if(incx.lt.0)ix = (-n+1)*incx + 1
C      if(incy.lt.0)iy = (-n+1)*incy + 1
C      do 10 i = 1,n
C        sy(iy) = sx(ix)
C        ix = ix + incx
C        iy = iy + incy
C   10 continue
C      return
Cc
Cc        code for both increments equal to 1
Cc
Cc
Cc        clean-up loop
Cc
C   20 m = mod(n,7)
C      if( m .eq. 0 ) go to 40
C      do 30 i = 1,m
C        sy(i) = sx(i)
C   30 continue
C      if( n .lt. 7 ) return
C   40 mp1 = m + 1
C      do 50 i = mp1,n,7
C        sy(i) = sx(i)
C        sy(i + 1) = sx(i + 1)
C        sy(i + 2) = sx(i + 2)
C        sy(i + 3) = sx(i + 3)
C        sy(i + 4) = sx(i + 4)
C        sy(i + 5) = sx(i + 5)
C        sy(i + 6) = sx(i + 6)
C   50 continue
C      return
C      end
C      real function sdot(n,sx,incx,sy,incy)
Cc
Cc     forms the dot product of two vectors.
Cc     uses unrolled loops for increments equal to one.
Cc     jack dongarra, linpack, 3/11/78.
Cc
C      real sx(1),sy(1),stemp
C      integer i,incx,incy,ix,iy,m,mp1,n
Cc
C      stemp = 0.0e0
C      sdot = 0.0e0
C      if(n.le.0)return
C      if(incx.eq.1.and.incy.eq.1)go to 20
Cc
Cc        code for unequal increments or equal increments
Cc          not equal to 1
Cc
C      ix = 1
C      iy = 1
C      if(incx.lt.0)ix = (-n+1)*incx + 1
C      if(incy.lt.0)iy = (-n+1)*incy + 1
C      do 10 i = 1,n
C        stemp = stemp + sx(ix)*sy(iy)
C        ix = ix + incx
C        iy = iy + incy
C   10 continue
C      sdot = stemp
C      return
Cc
Cc        code for both increments equal to 1
Cc
Cc
Cc        clean-up loop
Cc
C   20 m = mod(n,5)
C      if( m .eq. 0 ) go to 40
C      do 30 i = 1,m
C        stemp = stemp + sx(i)*sy(i)
C   30 continue
C      if( n .lt. 5 ) go to 60
C   40 mp1 = m + 1
C      do 50 i = mp1,n,5
C        stemp = stemp + sx(i)*sy(i) + sx(i + 1)*sy(i + 1) +
C     *   sx(i + 2)*sy(i + 2) + sx(i + 3)*sy(i + 3) + sx(i + 4)*sy(i + 4)
C   50 continue
C   60 sdot = stemp
C      return
C      end
C      subroutine sgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
C      integer lda,n,ml,mu,ipvt(1)
C      real abd(lda,1),z(1)
C      real rcond
Cc
Cc     sgbco factors a real band matrix by gaussian
Cc     elimination and estimates the condition of the matrix.
Cc
Cc     if  rcond  is not needed, sgbfa is slightly faster.
Cc     to solve  a*x = b , follow sgbco by sgbsl.
Cc     to compute  inverse(a)*c , follow sgbco by sgbsl.
Cc     to compute  determinant(a) , follow sgbco by sgbdi.
Cc
Cc     on entry
Cc
Cc        abd     real(lda, n)
Cc                contains the matrix in band storage.  the columns
Cc                of the matrix are stored in the columns of  abd  and
Cc                the diagonals of the matrix are stored in rows
Cc                ml+1 through 2*ml+mu+1 of  abd .
Cc                see the comments below for details.
Cc
Cc        lda     integer
Cc                the leading dimension of the array  abd .
Cc                lda must be .ge. 2*ml + mu + 1 .
Cc
Cc        n       integer
Cc                the order of the original matrix.
Cc
Cc        ml      integer
Cc                number of diagonals below the main diagonal.
Cc                0 .le. ml .lt. n .
Cc
Cc        mu      integer
Cc                number of diagonals above the main diagonal.
Cc                0 .le. mu .lt. n .
Cc                more efficient if  ml .le. mu .
Cc
Cc     on return
Cc
Cc        abd     an upper triangular matrix in band storage and
Cc                the multipliers which were used to obtain it.
Cc                the factorization can be written  a = l*u  where
Cc                l  is a product of permutation and unit lower
Cc                triangular matrices and  u  is upper triangular.
Cc
Cc        ipvt    integer(n)
Cc                an integer vector of pivot indices.
Cc
Cc        rcond   real
Cc                an estimate of the reciprocal condition of  a .
Cc                for the system  a*x = b , relative perturbations
Cc                in  a  and  b  of size  epsilon  may cause
Cc                relative perturbations in  x  of size  epsilon/rcond .
Cc                if  rcond  is so small that the logical expression
Cc                           1.0 + rcond .eq. 1.0
Cc                is true, then  a  may be singular to working
Cc                precision.  in particular,  rcond  is zero  if
Cc                exact singularity is detected or the estimate
Cc                underflows.
Cc
Cc        z       real(n)
Cc                a work vector whose contents are usually unimportant.
Cc                if  a  is close to a singular matrix, then  z  is
Cc                an approximate null vector in the sense that
Cc                norm(a*z) = rcond*norm(a)*norm(z) .
Cc
Cc     band storage
Cc
Cc           if  a  is a band matrix, the following program segment
Cc           will set up the input.
Cc
Cc                   ml = (band width below the diagonal)
Cc                   mu = (band width above the diagonal)
Cc                   m = ml + mu + 1
Cc                   do 20 j = 1, n
Cc                      i1 = max0(1, j-mu)
Cc                      i2 = min0(n, j+ml)
Cc                      do 10 i = i1, i2
Cc                         k = i - j + m
Cc                         abd(k,j) = a(i,j)
Cc                10    continue
Cc                20 continue
Cc
Cc           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
Cc           in addition, the first  ml  rows in  abd  are used for
Cc           elements generated during the triangularization.
Cc           the total number of rows needed in  abd  is  2*ml+mu+1 .
Cc           the  ml+mu by ml+mu  upper left triangle and the
Cc           ml by ml  lower right triangle are not referenced.
Cc
Cc     example..  if the original matrix is
Cc
Cc           11 12 13  0  0  0
Cc           21 22 23 24  0  0
Cc            0 32 33 34 35  0
Cc            0  0 43 44 45 46
Cc            0  0  0 54 55 56
Cc            0  0  0  0 65 66
Cc
Cc      then  n = 6, ml = 1, mu = 2, lda .ge. 5  and abd should contain
Cc
Cc            *  *  *  +  +  +  , * = not used
Cc            *  * 13 24 35 46  , + = used for pivoting
Cc            * 12 23 34 45 56
Cc           11 22 33 44 55 66
Cc           21 32 43 54 65  *
Cc
Cc     linpack. this version dated 08/14/78 .
Cc     cleve moler, university of new mexico, argonne national lab.
Cc
Cc     subroutines and functions
Cc
Cc     linpack sgbfa
Cc     blas saxpy,sdot,sscal,sasum
Cc     fortran abs,amax1,max0,min0,sign
Cc
Cc     internal variables
Cc
C      real sdot,ek,t,wk,wkm
C      real anorm,s,sasum,sm,ynorm
C      integer is,info,j,ju,k,kb,kp1,l,la,lm,lz,m,mm
Cc
Cc
Cc     compute 1-norm of a
Cc
C      anorm = 0.0e0
C      l = ml + 1
C      is = l + mu
C      do 10 j = 1, n
C         anorm = amax1(anorm,sasum(l,abd(is,j),1))
C         if (is .gt. ml + 1) is = is - 1
C         if (j .le. mu) l = l + 1
C         if (j .ge. n - ml) l = l - 1
C   10 continue
Cc
Cc     factor
Cc
C      call sgbfa(abd,lda,n,ml,mu,ipvt,info)
Cc
Cc     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
Cc     estimate = norm(z)/norm(y) where  a*z = y  and  trans(a)*y = e .
Cc     trans(a)  is the transpose of a .  the components of  e  are
Cc     chosen to cause maximum local growth in the elements of w  where
Cc     trans(u)*w = e .  the vectors are frequently rescaled to avoid
Cc     overflow.
Cc
Cc     solve trans(u)*w = e
Cc
C      ek = 1.0e0
C      do 20 j = 1, n
C         z(j) = 0.0e0
C   20 continue
C      m = ml + mu + 1
C      ju = 0
C      do 100 k = 1, n
C         if (z(k) .ne. 0.0e0) ek = sign(ek,-z(k))
C         if (abs(ek-z(k)) .le. abs(abd(m,k))) go to 30
C            s = abs(abd(m,k))/abs(ek-z(k))
C            call sscal(n,s,z,1)
C            ek = s*ek
C   30    continue
C         wk = ek - z(k)
C         wkm = -ek - z(k)
C         s = abs(wk)
C         sm = abs(wkm)
C         if (abd(m,k) .eq. 0.0e0) go to 40
C            wk = wk/abd(m,k)
C            wkm = wkm/abd(m,k)
C         go to 50
C   40    continue
C            wk = 1.0e0
C            wkm = 1.0e0
C   50    continue
C         kp1 = k + 1
C         ju = min0(max0(ju,mu+ipvt(k)),n)
C         mm = m
C         if (kp1 .gt. ju) go to 90
C            do 60 j = kp1, ju
C               mm = mm - 1
C               sm = sm + abs(z(j)+wkm*abd(mm,j))
C               z(j) = z(j) + wk*abd(mm,j)
C               s = s + abs(z(j))
C   60       continue
C            if (s .ge. sm) go to 80
C               t = wkm - wk
C               wk = wkm
C               mm = m
C               do 70 j = kp1, ju
C                  mm = mm - 1
C                  z(j) = z(j) + t*abd(mm,j)
C   70          continue
C   80       continue
C   90    continue
C         z(k) = wk
C  100 continue
C      s = 1.0e0/sasum(n,z,1)
C      call sscal(n,s,z,1)
Cc
Cc     solve trans(l)*y = w
Cc
C      do 120 kb = 1, n
C         k = n + 1 - kb
C         lm = min0(ml,n-k)
C         if (k .lt. n) z(k) = z(k) + sdot(lm,abd(m+1,k),1,z(k+1),1)
C         if (abs(z(k)) .le. 1.0e0) go to 110
C            s = 1.0e0/abs(z(k))
C            call sscal(n,s,z,1)
C  110    continue
C         l = ipvt(k)
C         t = z(l)
C         z(l) = z(k)
C         z(k) = t
C  120 continue
C      s = 1.0e0/sasum(n,z,1)
C      call sscal(n,s,z,1)
Cc
C      ynorm = 1.0e0
Cc
Cc     solve l*v = y
Cc
C      do 140 k = 1, n
C         l = ipvt(k)
C         t = z(l)
C         z(l) = z(k)
C         z(k) = t
C         lm = min0(ml,n-k)
C         if (k .lt. n) call saxpy(lm,t,abd(m+1,k),1,z(k+1),1)
C         if (abs(z(k)) .le. 1.0e0) go to 130
C            s = 1.0e0/abs(z(k))
C            call sscal(n,s,z,1)
C            ynorm = s*ynorm
C  130    continue
C  140 continue
C      s = 1.0e0/sasum(n,z,1)
C      call sscal(n,s,z,1)
C      ynorm = s*ynorm
Cc
Cc     solve  u*z = w
Cc
C      do 160 kb = 1, n
C         k = n + 1 - kb
C         if (abs(z(k)) .le. abs(abd(m,k))) go to 150
C            s = abs(abd(m,k))/abs(z(k))
C            call sscal(n,s,z,1)
C            ynorm = s*ynorm
C  150    continue
C         if (abd(m,k) .ne. 0.0e0) z(k) = z(k)/abd(m,k)
C         if (abd(m,k) .eq. 0.0e0) z(k) = 1.0e0
C         lm = min0(k,m) - 1
C         la = m - lm
C         lz = k - lm
C         t = -z(k)
C         call saxpy(lm,t,abd(la,k),1,z(lz),1)
C  160 continue
Cc     make znorm = 1.0
C      s = 1.0e0/sasum(n,z,1)
C      call sscal(n,s,z,1)
C      ynorm = s*ynorm
Cc
C      if (anorm .ne. 0.0e0) rcond = ynorm/anorm
C      if (anorm .eq. 0.0e0) rcond = 0.0e0
C      return
C      end
C      subroutine sgbfa(abd,lda,n,ml,mu,ipvt,info)
C      integer lda,n,ml,mu,ipvt(1),info
C      real abd(lda,1)
Cc
Cc     sgbfa factors a real band matrix by elimination.
Cc
Cc     sgbfa is usually called by sgbco, but it can be called
Cc     directly with a saving in time if  rcond  is not needed.
Cc
Cc     on entry
Cc
Cc        abd     real(lda, n)
Cc                contains the matrix in band storage.  the columns
Cc                of the matrix are stored in the columns of  abd  and
Cc                the diagonals of the matrix are stored in rows
Cc                ml+1 through 2*ml+mu+1 of  abd .
Cc                see the comments below for details.
Cc
Cc        lda     integer
Cc                the leading dimension of the array  abd .
Cc                lda must be .ge. 2*ml + mu + 1 .
Cc
Cc        n       integer
Cc                the order of the original matrix.
Cc
Cc        ml      integer
Cc                number of diagonals below the main diagonal.
Cc                0 .le. ml .lt. n .
Cc
Cc        mu      integer
Cc                number of diagonals above the main diagonal.
Cc                0 .le. mu .lt. n .
Cc                more efficient if  ml .le. mu .
Cc     on return
Cc
Cc        abd     an upper triangular matrix in band storage and
Cc                the multipliers which were used to obtain it.
Cc                the factorization can be written  a = l*u  where
Cc                l  is a product of permutation and unit lower
Cc                triangular matrices and  u  is upper triangular.
Cc
Cc        ipvt    integer(n)
Cc                an integer vector of pivot indices.
Cc
Cc        info    integer
Cc                = 0  normal value.
Cc                = k  if  u(k,k) .eq. 0.0 .  this is not an error
Cc                     condition for this subroutine, but it does
Cc                     indicate that sgbsl will divide by zero if
Cc                     called.  use  rcond  in sgbco for a reliable
Cc                     indication of singularity.
Cc
Cc     band storage
Cc
Cc           if  a  is a band matrix, the following program segment
Cc           will set up the input.
Cc
Cc                   ml = (band width below the diagonal)
Cc                   mu = (band width above the diagonal)
Cc                   m = ml + mu + 1
Cc                   do 20 j = 1, n
Cc                      i1 = max0(1, j-mu)
Cc                      i2 = min0(n, j+ml)
Cc                      do 10 i = i1, i2
Cc                         k = i - j + m
Cc                         abd(k,j) = a(i,j)
Cc                10    continue
Cc                20 continue
Cc
Cc           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
Cc           in addition, the first  ml  rows in  abd  are used for
Cc           elements generated during the triangularization.
Cc           the total number of rows needed in  abd  is  2*ml+mu+1 .
Cc           the  ml+mu by ml+mu  upper left triangle and the
Cc           ml by ml  lower right triangle are not referenced.
Cc
Cc     linpack. this version dated 08/14/78 .
Cc     cleve moler, university of new mexico, argonne national lab.
Cc
Cc     subroutines and functions
Cc
Cc     blas saxpy,sscal,isamax
Cc     fortran max0,min0
Cc
Cc     internal variables
Cc
C      real t
C      integer i,isamax,i0,j,ju,jz,j0,j1,k,kp1,l,lm,m,mm,nm1
Cc
Cc
C      m = ml + mu + 1
C      info = 0
Cc
Cc     zero initial fill-in columns
Cc
C      j0 = mu + 2
C      j1 = min0(n,m) - 1
C      if (j1 .lt. j0) go to 30
C      do 20 jz = j0, j1
C         i0 = m + 1 - jz
C         do 10 i = i0, ml
C            abd(i,jz) = 0.0e0
C   10    continue
C   20 continue
C   30 continue
C      jz = j1
C      ju = 0
Cc
Cc     gaussian elimination with partial pivoting
Cc
C      nm1 = n - 1
C      if (nm1 .lt. 1) go to 130
C      do 120 k = 1, nm1
C         kp1 = k + 1
Cc
Cc        zero next fill-in column
Cc
C         jz = jz + 1
C         if (jz .gt. n) go to 50
C         if (ml .lt. 1) go to 50
C            do 40 i = 1, ml
C               abd(i,jz) = 0.0e0
C   40       continue
C   50    continue
Cc
Cc        find l = pivot index
Cc
C         lm = min0(ml,n-k)
C         l = isamax(lm+1,abd(m,k),1) + m - 1
C         ipvt(k) = l + k - m
Cc
Cc        zero pivot implies this column already triangularized
Cc
C         if (abd(l,k) .eq. 0.0e0) go to 100
Cc
Cc           interchange if necessary
Cc
C            if (l .eq. m) go to 60
C               t = abd(l,k)
C               abd(l,k) = abd(m,k)
C               abd(m,k) = t
C   60       continue
Cc
Cc           compute multipliers
Cc
C            t = -1.0e0/abd(m,k)
C            call sscal(lm,t,abd(m+1,k),1)
Cc
Cc           row elimination with column indexing
Cc
C            ju = min0(max0(ju,mu+ipvt(k)),n)
C            mm = m
C            if (ju .lt. kp1) go to 90
C            do 80 j = kp1, ju
C               l = l - 1
C               mm = mm - 1
C               t = abd(l,j)
C               if (l .eq. mm) go to 70
C                  abd(l,j) = abd(mm,j)
C                  abd(mm,j) = t
C   70          continue
C               call saxpy(lm,t,abd(m+1,k),1,abd(mm+1,j),1)
C   80       continue
C   90       continue
C         go to 110
C  100    continue
C            info = k
C  110    continue
C  120 continue
C  130 continue
C      ipvt(n) = n
C      if (abd(m,n) .eq. 0.0e0) info = n
C      return
C      end
C      subroutine sgbsl(abd,lda,n,ml,mu,ipvt,b,job)
C      integer lda,n,ml,mu,ipvt(1),job
C      real abd(lda,1),b(1)
Cc
Cc     sgbsl solves the real band system
Cc     a * x = b  or  trans(a) * x = b
Cc     using the factors computed by sgbco or sgbfa.
Cc
Cc     on entry
Cc
Cc        abd     real(lda, n)
Cc                the output from sgbco or sgbfa.
Cc
Cc        lda     integer
Cc                the leading dimension of the array  abd .
Cc
Cc        n       integer
Cc                the order of the original matrix.
Cc
Cc        ml      integer
Cc                number of diagonals below the main diagonal.
Cc
Cc        mu      integer
Cc                number of diagonals above the main diagonal.
Cc
Cc        ipvt    integer(n)
Cc                the pivot vector from sgbco or sgbfa.
Cc
Cc        b       real(n)
Cc                the right hand side vector.
Cc
Cc        job     integer
Cc                = 0         to solve  a*x = b ,
Cc                = nonzero   to solve  trans(a)*x = b , where
Cc                            trans(a)  is the transpose.
Cc
Cc     on return
Cc
Cc        b       the solution vector  x .
Cc
Cc     error condition
Cc
Cc        a division by zero will occur if the input factor contains a
Cc        zero on the diagonal.  technically this indicates singularity
Cc        but it is often caused by improper arguments or improper
Cc        setting of lda .  it will not occur if the subroutines are
Cc        called correctly and if sgbco has set rcond .gt. 0.0
Cc        or sgbfa has set info .eq. 0 .
Cc
Cc     to compute  inverse(a) * c  where  c  is a matrix
Cc     with  p  columns
Cc           call sgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
Cc           if (rcond is too small) go to ...
Cc           do 10 j = 1, p
Cc              call sgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
Cc        10 continue
Cc
Cc     linpack. this version dated 08/14/78 .
Cc     cleve moler, university of new mexico, argonne national lab.
Cc
Cc     subroutines and functions
Cc
Cc     blas saxpy,sdot
Cc     fortran min0
Cc
Cc     internal variables
Cc
C      real sdot,t
C      integer k,kb,l,la,lb,lm,m,nm1
Cc
C      m = mu + ml + 1
C      nm1 = n - 1
C      if (job .ne. 0) go to 50
Cc
Cc        job = 0 , solve  a * x = b
Cc        first solve l*y = b
Cc
C         if (ml .eq. 0) go to 30
C         if (nm1 .lt. 1) go to 30
C            do 20 k = 1, nm1
C               lm = min0(ml,n-k)
C               l = ipvt(k)
C               t = b(l)
C               if (l .eq. k) go to 10
C                  b(l) = b(k)
C                  b(k) = t
C   10          continue
C               call saxpy(lm,t,abd(m+1,k),1,b(k+1),1)
C   20       continue
C   30    continue
Cc
Cc        now solve  u*x = y
Cc
C         do 40 kb = 1, n
C            k = n + 1 - kb
C            b(k) = b(k)/abd(m,k)
C            lm = min0(k,m) - 1
C            la = m - lm
C            lb = k - lm
C            t = -b(k)
C            call saxpy(lm,t,abd(la,k),1,b(lb),1)
C   40    continue
C      go to 100
C   50 continue
Cc
Cc        job = nonzero, solve  trans(a) * x = b
Cc        first solve  trans(u)*y = b
Cc
C         do 60 k = 1, n
C            lm = min0(k,m) - 1
C            la = m - lm
C            lb = k - lm
C            t = sdot(lm,abd(la,k),1,b(lb),1)
C            b(k) = (b(k) - t)/abd(m,k)
C   60    continue
Cc
Cc        now solve trans(l)*x = y
Cc
C         if (ml .eq. 0) go to 90
C         if (nm1 .lt. 1) go to 90
C            do 80 kb = 1, nm1
C               k = n - kb
C               lm = min0(ml,n-k)
C               b(k) = b(k) + sdot(lm,abd(m+1,k),1,b(k+1),1)
C               l = ipvt(k)
C               if (l .eq. k) go to 70
C                  t = b(l)
C                  b(l) = b(k)
C                  b(k) = t
C   70          continue
C   80       continue
C   90    continue
C  100 continue
C      return
C      end
C      subroutine sgeco(a,lda,n,ipvt,rcond,z)
C      integer lda,n,ipvt(1)
C      real a(lda,1),z(1)
C      real rcond
Cc
Cc     sgeco factors a real matrix by gaussian elimination
Cc     and estimates the condition of the matrix.
Cc
Cc     if  rcond  is not needed, sgefa is slightly faster.
Cc     to solve  a*x = b , follow sgeco by sgesl.
Cc     to compute  inverse(a)*c , follow sgeco by sgesl.
Cc     to compute  determinant(a) , follow sgeco by sgedi.
Cc     to compute  inverse(a) , follow sgeco by sgedi.
Cc
Cc     on entry
Cc
Cc        a       real(lda, n)
Cc                the matrix to be factored.
Cc
Cc        lda     integer
Cc                the leading dimension of the array  a .
Cc
Cc        n       integer
Cc                the order of the matrix  a .
Cc
Cc     on return
Cc
Cc        a       an upper triangular matrix and the multipliers
Cc                which were used to obtain it.
Cc                the factorization can be written  a = l*u  where
Cc                l  is a product of permutation and unit lower
Cc                triangular matrices and  u  is upper triangular.
Cc
Cc        ipvt    integer(n)
Cc                an integer vector of pivot indices.
Cc
Cc        rcond   real
Cc                an estimate of the reciprocal condition of  a .
Cc                for the system  a*x = b , relative perturbations
Cc                in  a  and  b  of size  epsilon  may cause
Cc                relative perturbations in  x  of size  epsilon/rcond .
Cc                if  rcond  is so small that the logical expression
Cc                           1.0 + rcond .eq. 1.0
Cc                is true, then  a  may be singular to working
Cc                precision.  in particular,  rcond  is zero  if
Cc                exact singularity is detected or the estimate
Cc                underflows.
Cc
Cc        z       real(n)
Cc                a work vector whose contents are usually unimportant.
Cc                if  a  is close to a singular matrix, then  z  is
Cc                an approximate null vector in the sense that
Cc                norm(a*z) = rcond*norm(a)*norm(z) .
Cc
Cc     linpack. this version dated 08/14/78 .
Cc     cleve moler, university of new mexico, argonne national lab.
Cc
Cc     subroutines and functions
Cc
Cc     linpack sgefa
Cc     blas saxpy,sdot,sscal,sasum
Cc     fortran abs,amax1,sign
Cc
Cc     internal variables
Cc
C      real sdot,ek,t,wk,wkm
C      real anorm,s,sasum,sm,ynorm
C      integer info,j,k,kb,kp1,l
Cc
Cc
Cc     compute 1-norm of a
Cc
C      anorm = 0.0e0
C      do 10 j = 1, n
C         anorm = amax1(anorm,sasum(n,a(1,j),1))
C   10 continue
Cc
Cc     factor
Cc
C      call sgefa(a,lda,n,ipvt,info)
Cc
Cc     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
Cc     estimate = norm(z)/norm(y) where  a*z = y  and  trans(a)*y = e .
Cc     trans(a)  is the transpose of a .  the components of  e  are
Cc     chosen to cause maximum local growth in the elements of w  where
Cc     trans(u)*w = e .  the vectors are frequently rescaled to avoid
Cc     overflow.
Cc
Cc     solve trans(u)*w = e
Cc
C      ek = 1.0e0
C      do 20 j = 1, n
C         z(j) = 0.0e0
C   20 continue
C      do 100 k = 1, n
C         if (z(k) .ne. 0.0e0) ek = sign(ek,-z(k))
C         if (abs(ek-z(k)) .le. abs(a(k,k))) go to 30
C            s = abs(a(k,k))/abs(ek-z(k))
C            call sscal(n,s,z,1)
C            ek = s*ek
C   30    continue
C         wk = ek - z(k)
C         wkm = -ek - z(k)
C         s = abs(wk)
C         sm = abs(wkm)
C         if (a(k,k) .eq. 0.0e0) go to 40
C            wk = wk/a(k,k)
C            wkm = wkm/a(k,k)
C         go to 50
C   40    continue
C            wk = 1.0e0
C            wkm = 1.0e0
C   50    continue
C         kp1 = k + 1
C         if (kp1 .gt. n) go to 90
C            do 60 j = kp1, n
C               sm = sm + abs(z(j)+wkm*a(k,j))
C               z(j) = z(j) + wk*a(k,j)
C               s = s + abs(z(j))
C   60       continue
C            if (s .ge. sm) go to 80
C               t = wkm - wk
C               wk = wkm
C               do 70 j = kp1, n
C                  z(j) = z(j) + t*a(k,j)
C   70          continue
C   80       continue
C   90    continue
C         z(k) = wk
C  100 continue
C      s = 1.0e0/sasum(n,z,1)
C      call sscal(n,s,z,1)
Cc
Cc     solve trans(l)*y = w
Cc
C      do 120 kb = 1, n
C         k = n + 1 - kb
C         if (k .lt. n) z(k) = z(k) + sdot(n-k,a(k+1,k),1,z(k+1),1)
C         if (abs(z(k)) .le. 1.0e0) go to 110
C            s = 1.0e0/abs(z(k))
C            call sscal(n,s,z,1)
C  110    continue
C         l = ipvt(k)
C         t = z(l)
C         z(l) = z(k)
C         z(k) = t
C  120 continue
C      s = 1.0e0/sasum(n,z,1)
C      call sscal(n,s,z,1)
Cc
C      ynorm = 1.0e0
Cc
Cc     solve l*v = y
Cc
C      do 140 k = 1, n
C         l = ipvt(k)
C         t = z(l)
C         z(l) = z(k)
C         z(k) = t
C         if (k .lt. n) call saxpy(n-k,t,a(k+1,k),1,z(k+1),1)
C         if (abs(z(k)) .le. 1.0e0) go to 130
C            s = 1.0e0/abs(z(k))
C            call sscal(n,s,z,1)
C            ynorm = s*ynorm
C  130    continue
C  140 continue
C      s = 1.0e0/sasum(n,z,1)
C      call sscal(n,s,z,1)
C      ynorm = s*ynorm
Cc
Cc     solve  u*z = v
Cc
C      do 160 kb = 1, n
C         k = n + 1 - kb
C         if (abs(z(k)) .le. abs(a(k,k))) go to 150
C            s = abs(a(k,k))/abs(z(k))
C            call sscal(n,s,z,1)
C            ynorm = s*ynorm
C  150    continue
C         if (a(k,k) .ne. 0.0e0) z(k) = z(k)/a(k,k)
C         if (a(k,k) .eq. 0.0e0) z(k) = 1.0e0
C         t = -z(k)
C         call saxpy(k-1,t,a(1,k),1,z(1),1)
C  160 continue
Cc     make znorm = 1.0
C      s = 1.0e0/sasum(n,z,1)
C      call sscal(n,s,z,1)
C      ynorm = s*ynorm
Cc
C      if (anorm .ne. 0.0e0) rcond = ynorm/anorm
C      if (anorm .eq. 0.0e0) rcond = 0.0e0
C      return
C      end
C      subroutine sgedi(a,lda,n,ipvt,det,work,job)
C      integer lda,n,ipvt(1),job
C      real a(lda,1),det(2),work(1)
Cc
Cc     sgedi computes the determinant and inverse of a matrix
Cc     using the factors computed by sgeco or sgefa.
Cc
Cc     on entry
Cc
Cc        a       real(lda, n)
Cc                the output from sgeco or sgefa.
Cc
Cc        lda     integer
Cc                the leading dimension of the array  a .
Cc
Cc        n       integer
Cc                the order of the matrix  a .
Cc
Cc        ipvt    integer(n)
Cc                the pivot vector from sgeco or sgefa.
Cc
Cc        work    real(n)
Cc                work vector.  contents destroyed.
Cc
Cc        job     integer
Cc                = 11   both determinant and inverse.
Cc                = 01   inverse only.
Cc                = 10   determinant only.
Cc
Cc     on return
Cc
Cc        a       inverse of original matrix if requested.
Cc                otherwise unchanged.
Cc
Cc        det     real(2)
Cc                determinant of original matrix if requested.
Cc                otherwise not referenced.
Cc                determinant = det(1) * 10.0**det(2)
Cc                with  1.0 .le. abs(det(1)) .lt. 10.0
Cc                or  det(1) .eq. 0.0 .
Cc
Cc     error condition
Cc
Cc        a division by zero will occur if the input factor contains
Cc        a zero on the diagonal and the inverse is requested.
Cc        it will not occur if the subroutines are called correctly
Cc        and if sgeco has set rcond .gt. 0.0 or sgefa has set
Cc        info .eq. 0 .
Cc
Cc     linpack. this version dated 08/14/78 .
Cc     cleve moler, university of new mexico, argonne national lab.
Cc
Cc     subroutines and functions
Cc
Cc     blas saxpy,sscal,sswap
Cc     fortran abs,mod
Cc
Cc     internal variables
Cc
C      real t
C      real ten
C      integer i,j,k,kb,kp1,l,nm1
Cc
Cc
Cc     compute determinant
Cc
C      if (job/10 .eq. 0) go to 70
C         det(1) = 1.0e0
C         det(2) = 0.0e0
C         ten = 10.0e0
C         do 50 i = 1, n
C            if (ipvt(i) .ne. i) det(1) = -det(1)
C            det(1) = a(i,i)*det(1)
Cc        ...exit
C            if (det(1) .eq. 0.0e0) go to 60
C   10       if (abs(det(1)) .ge. 1.0e0) go to 20
C               det(1) = ten*det(1)
C               det(2) = det(2) - 1.0e0
C            go to 10
C   20       continue
C   30       if (abs(det(1)) .lt. ten) go to 40
C               det(1) = det(1)/ten
C               det(2) = det(2) + 1.0e0
C            go to 30
C   40       continue
C   50    continue
C   60    continue
C   70 continue
Cc
Cc     compute inverse(u)
Cc
C      if (mod(job,10) .eq. 0) go to 150
C         do 100 k = 1, n
C            a(k,k) = 1.0e0/a(k,k)
C            t = -a(k,k)
C            call sscal(k-1,t,a(1,k),1)
C            kp1 = k + 1
C            if (n .lt. kp1) go to 90
C            do 80 j = kp1, n
C               t = a(k,j)
C               a(k,j) = 0.0e0
C               call saxpy(k,t,a(1,k),1,a(1,j),1)
C   80       continue
C   90       continue
C  100    continue
Cc
Cc        form inverse(u)*inverse(l)
Cc
C         nm1 = n - 1
C         if (nm1 .lt. 1) go to 140
C         do 130 kb = 1, nm1
C            k = n - kb
C            kp1 = k + 1
C            do 110 i = kp1, n
C               work(i) = a(i,k)
C               a(i,k) = 0.0e0
C  110       continue
C            do 120 j = kp1, n
C               t = work(j)
C               call saxpy(n,t,a(1,j),1,a(1,k),1)
C  120       continue
C            l = ipvt(k)
C            if (l .ne. k) call sswap(n,a(1,k),1,a(1,l),1)
C  130    continue
C  140    continue
C  150 continue
C      return
C      end
C      subroutine sgefa(a,lda,n,ipvt,info)
C      integer lda,n,ipvt(1),info
C      real a(lda,1)
Cc
Cc     sgefa factors a real matrix by gaussian elimination.
Cc
Cc     sgefa is usually called by sgeco, but it can be called
Cc     directly with a saving in time if  rcond  is not needed.
Cc     (time for sgeco) = (1 + 9/n)*(time for sgefa) .
Cc
Cc     on entry
Cc
Cc        a       real(lda, n)
Cc                the matrix to be factored.
Cc
Cc        lda     integer
Cc                the leading dimension of the array  a .
Cc
Cc        n       integer
Cc                the order of the matrix  a .
Cc
Cc     on return
Cc
Cc        a       an upper triangular matrix and the multipliers
Cc                which were used to obtain it.
Cc                the factorization can be written  a = l*u  where
Cc                l  is a product of permutation and unit lower
Cc                triangular matrices and  u  is upper triangular.
Cc
Cc        ipvt    integer(n)
Cc                an integer vector of pivot indices.
Cc
Cc        info    integer
Cc                = 0  normal value.
Cc                = k  if  u(k,k) .eq. 0.0 .  this is not an error
Cc                     condition for this subroutine, but it does
Cc                     indicate that sgesl or sgedi will divide by zero
Cc                     if called.  use  rcond  in sgeco for a reliable
Cc                     indication of singularity.
Cc
Cc     linpack. this version dated 08/14/78 .
Cc     cleve moler, university of new mexico, argonne national lab.
Cc
Cc     subroutines and functions
Cc
Cc     blas saxpy,sscal,isamax
Cc
Cc     internal variables
Cc
C      real t
C      integer isamax,j,k,kp1,l,nm1
Cc
Cc
Cc     gaussian elimination with partial pivoting
Cc
C      info = 0
C      nm1 = n - 1
C      if (nm1 .lt. 1) go to 70
C      do 60 k = 1, nm1
C         kp1 = k + 1
Cc
Cc        find l = pivot index
Cc
C         l = isamax(n-k+1,a(k,k),1) + k - 1
C         ipvt(k) = l
Cc
Cc        zero pivot implies this column already triangularized
Cc
C         if (a(l,k) .eq. 0.0e0) go to 40
Cc
Cc           interchange if necessary
Cc
C            if (l .eq. k) go to 10
C               t = a(l,k)
C               a(l,k) = a(k,k)
C               a(k,k) = t
C   10       continue
Cc
Cc           compute multipliers
Cc
C            t = -1.0e0/a(k,k)
C            call sscal(n-k,t,a(k+1,k),1)
Cc
Cc           row elimination with column indexing
Cc
C            do 30 j = kp1, n
C               t = a(l,j)
C               if (l .eq. k) go to 20
C                  a(l,j) = a(k,j)
C                  a(k,j) = t
C   20          continue
C               call saxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
C   30       continue
C         go to 50
C   40    continue
C            info = k
C   50    continue
C   60 continue
C   70 continue
C      ipvt(n) = n
C      if (a(n,n) .eq. 0.0e0) info = n
C      return
C      end
C      subroutine sgesl(a,lda,n,ipvt,b,job)
C      integer lda,n,ipvt(1),job
C      real a(lda,1),b(1)
Cc
Cc     sgesl solves the real system
Cc     a * x = b  or  trans(a) * x = b
Cc     using the factors computed by sgeco or sgefa.
Cc
Cc     on entry
Cc
Cc        a       real(lda, n)
Cc                the output from sgeco or sgefa.
Cc
Cc        lda     integer
Cc                the leading dimension of the array  a .
Cc
Cc        n       integer
Cc                the order of the matrix  a .
Cc
Cc        ipvt    integer(n)
Cc                the pivot vector from sgeco or sgefa.
Cc
Cc        b       real(n)
Cc                the right hand side vector.
Cc
Cc        job     integer
Cc                = 0         to solve  a*x = b ,
Cc                = nonzero   to solve  trans(a)*x = b  where
Cc                            trans(a)  is the transpose.
Cc
Cc     on return
Cc
Cc        b       the solution vector  x .
Cc
Cc     error condition
Cc
Cc        a division by zero will occur if the input factor contains a
Cc        zero on the diagonal.  technically this indicates singularity
Cc        but it is often caused by improper arguments or improper
Cc        setting of lda .  it will not occur if the subroutines are
Cc        called correctly and if sgeco has set rcond .gt. 0.0
Cc        or sgefa has set info .eq. 0 .
Cc
Cc     to compute  inverse(a) * c  where  c  is a matrix
Cc     with  p  columns
Cc           call sgeco(a,lda,n,ipvt,rcond,z)
Cc           if (rcond is too small) go to ...
Cc           do 10 j = 1, p
Cc              call sgesl(a,lda,n,ipvt,c(1,j),0)
Cc        10 continue
Cc
Cc     linpack. this version dated 08/14/78 .
Cc     cleve moler, university of new mexico, argonne national lab.
Cc
Cc     subroutines and functions
Cc
Cc     blas saxpy,sdot
Cc
Cc     internal variables
Cc
C      real sdot,t
C      integer k,kb,l,nm1
Cc
C      nm1 = n - 1
C      if (job .ne. 0) go to 50
Cc
Cc        job = 0 , solve  a * x = b
Cc        first solve  l*y = b
Cc
C         if (nm1 .lt. 1) go to 30
C         do 20 k = 1, nm1
C            l = ipvt(k)
C            t = b(l)
C            if (l .eq. k) go to 10
C               b(l) = b(k)
C               b(k) = t
C   10       continue
C            call saxpy(n-k,t,a(k+1,k),1,b(k+1),1)
C   20    continue
C   30    continue
Cc
Cc        now solve  u*x = y
Cc
C         do 40 kb = 1, n
C            k = n + 1 - kb
C            b(k) = b(k)/a(k,k)
C            t = -b(k)
C            call saxpy(k-1,t,a(1,k),1,b(1),1)
C   40    continue
C      go to 100
C   50 continue
Cc
Cc        job = nonzero, solve  trans(a) * x = b
Cc        first solve  trans(u)*y = b
Cc
C         do 60 k = 1, n
C            t = sdot(k-1,a(1,k),1,b(1),1)
C            b(k) = (b(k) - t)/a(k,k)
C   60    continue
Cc
Cc        now solve trans(l)*x = y
Cc
C         if (nm1 .lt. 1) go to 90
C         do 80 kb = 1, nm1
C            k = n - kb
C            b(k) = b(k) + sdot(n-k,a(k+1,k),1,b(k+1),1)
C            l = ipvt(k)
C            if (l .eq. k) go to 70
C               t = b(l)
C               b(l) = b(k)
C               b(k) = t
C   70       continue
C   80    continue
C   90    continue
C  100 continue
C      return
C      end
C      subroutine sscal(n,sa,sx,incx)
Cc
Cc     scales a vector by a constant.
Cc     uses unrolled loops for increment equal to 1.
Cc     jack dongarra, linpack, 3/11/78.
Cc     modified 3/93 to return if incx .le. 0.
Cc
C      real sa,sx(1)
C      integer i,incx,m,mp1,n,nincx
Cc
C      if( n.le.0 .or. incx.le.0 )return
C      if(incx.eq.1)go to 20
Cc
Cc        code for increment not equal to 1
Cc
C      nincx = n*incx
C      do 10 i = 1,nincx,incx
C        sx(i) = sa*sx(i)
C   10 continue
C      return
Cc
Cc        code for increment equal to 1
Cc
Cc
Cc        clean-up loop
Cc
C   20 m = mod(n,5)
C      if( m .eq. 0 ) go to 40
C      do 30 i = 1,m
C        sx(i) = sa*sx(i)
C   30 continue
C      if( n .lt. 5 ) return
C   40 mp1 = m + 1
C      do 50 i = mp1,n,5
C        sx(i) = sa*sx(i)
C        sx(i + 1) = sa*sx(i + 1)
C        sx(i + 2) = sa*sx(i + 2)
C        sx(i + 3) = sa*sx(i + 3)
C        sx(i + 4) = sa*sx(i + 4)
C   50 continue
C      return
C      end
C      subroutine sswap (n,sx,incx,sy,incy)
Cc
Cc     interchanges two vectors.
Cc     uses unrolled loops for increments equal to 1.
Cc     jack dongarra, linpack, 3/11/78.
Cc
C      real sx(1),sy(1),stemp
C      integer i,incx,incy,ix,iy,m,mp1,n
Cc
C      if(n.le.0)return
C      if(incx.eq.1.and.incy.eq.1)go to 20
Cc
Cc       code for unequal increments or equal increments not equal
Cc         to 1
Cc
C      ix = 1
C      iy = 1
C      if(incx.lt.0)ix = (-n+1)*incx + 1
C      if(incy.lt.0)iy = (-n+1)*incy + 1
C      do 10 i = 1,n
C        stemp = sx(ix)
C        sx(ix) = sy(iy)
C        sy(iy) = stemp
C        ix = ix + incx
C        iy = iy + incy
C   10 continue
C      return
Cc
Cc       code for both increments equal to 1
Cc
Cc
Cc       clean-up loop
Cc
C   20 m = mod(n,3)
C      if( m .eq. 0 ) go to 40
C      do 30 i = 1,m
C        stemp = sx(i)
C        sx(i) = sy(i)
C        sy(i) = stemp
C   30 continue
C      if( n .lt. 3 ) return
C   40 mp1 = m + 1
C      do 50 i = mp1,n,3
C        stemp = sx(i)
C        sx(i) = sy(i)
C        sy(i) = stemp
C        stemp = sx(i + 1)
C        sx(i + 1) = sy(i + 1)
C        sy(i + 1) = stemp
C        stemp = sx(i + 2)
C        sx(i + 2) = sy(i + 2)
C        sy(i + 2) = stemp
C   50 continue
C      return
C      end
C      integer function isamax(n,sx,incx)
Cc
Cc     finds the index of element having max. absolute value.
Cc     jack dongarra, linpack, 3/11/78.
Cc     modified 3/93 to return if incx .le. 0.
Cc
C      real sx(1),smax
C      integer i,incx,ix,n
Cc
C      isamax = 0
C      if( n.lt.1 .or. incx.le.0 ) return
C      isamax = 1
C      if(n.eq.1)return
C      if(incx.eq.1)go to 20
Cc
Cc        code for increment not equal to 1
Cc
C      ix = 1
C      smax = abs(sx(1))
C      ix = ix + incx
C      do 10 i = 2,n
C         if(abs(sx(ix)).le.smax) go to 5
C         isamax = i
C         smax = abs(sx(ix))
C    5    ix = ix + incx
C   10 continue
C      return
Cc
Cc        code for increment equal to 1
Cc
C   20 smax = abs(sx(1))
C      do 30 i = 2,n
C         if(abs(sx(i)).le.smax) go to 30
C         isamax = i
C         smax = abs(sx(i))
C   30 continue
C      return
C      end
C*****END precision > single
C
      INTEGER          FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3,
     $                 N4 )
*
*  -- LAPACK auxiliary routine (preliminary version) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 20, 1992
*
*     .. Scalar Arguments ..
      CHARACTER*( * )    NAME, OPTS
      INTEGER            ISPEC, N1, N2, N3, N4
*     ..
*
*  Purpose
*  =======
*
*  ILAENV is called from the LAPACK routines to choose problem-dependent
*  parameters for the local environment.  See ISPEC for a description of
*  the parameters.
*
*  This version provides a set of parameters which should give good,
*  but not optimal, performance on many of the currently available
*  computers.  Users are encouraged to modify this subroutine to set
*  the tuning parameters for their particular machine using the option
*  and problem size information in the arguments.
*
*  This routine will not function correctly if it is converted to all
*  lower case.  Converting it to all upper case is allowed.
*
*  Arguments
*  =========
*
*  ISPEC   (input) INTEGER
*          Specifies the parameter to be returned as the value of
*          ILAENV.
*          = 1: the optimal blocksize; if this value is 1, an unblocked
*               algorithm will give the best performance.
*          = 2: the minimum block size for which the block routine
*               should be used; if the usable block size is less than
*               this value, an unblocked routine should be used.
*          = 3: the crossover point (in a block routine, for N less
*               than this value, an unblocked routine should be used)
*          = 4: the number of shifts, used in the nonsymmetric
*               eigenvalue routines
*          = 5: the minimum column dimension for blocking to be used;
*               rectangular blocks must have dimension at least k by m,
*               where k is given by ILAENV(2,...) and m by ILAENV(5,...)
*          = 6: the crossover point for the SVD (when reducing an m by n
*               matrix to bidiagonal form, if max(m,n)/min(m,n) exceeds
*               this value, a QR factorization is used first to reduce
*               the matrix to a triangular form.)
*          = 7: the number of processors
*          = 8: the crossover point for the multishift QR and QZ methods
*               for nonsymmetric eigenvalue problems.
*
*  NAME    (input) CHARACTER*(*)
*          The name of the calling subroutine, in either upper case or
*          lower case.
*
*  OPTS    (input) CHARACTER*(*)
*          The character options to the subroutine NAME, concatenated
*          into a single character string.  For example, UPLO = 'U',
*          TRANS = 'T', and DIAG = 'N' for a triangular routine would
*          be specified as OPTS = 'UTN'.
*
*  N1      (input) INTEGER
*  N2      (input) INTEGER
*  N3      (input) INTEGER
*  N4      (input) INTEGER
*          Problem dimensions for the subroutine NAME; these may not all
*          be required.
*
* (ILAENV) (output) INTEGER
*          >= 0: the value of the parameter specified by ISPEC
*          < 0:  if ILAENV = -k, the k-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  The following conventions have been used when calling ILAENV from the
*  LAPACK routines:
*  1)  OPTS is a concatenation of all of the character options to
*      subroutine NAME, in the same order that they appear in the
*      argument list for NAME, even if they are not used in determining
*      the value of the parameter specified by ISPEC.
*  2)  The problem dimensions N1, N2, N3, N4 are specified in the order
*      that they appear in the argument list for NAME.  N1 is used
*      first, N2 second, and so on, and unused problem dimensions are
*      passed a value of -1.
*  3)  The parameter value returned by ILAENV is checked for validity in
*      the calling subroutine.  For example, ILAENV is used to retrieve
*      the optimal blocksize for STRTRI as follows:
*
*      NB = ILAENV( 1, 'STRTRI', UPLO // DIAG, N, -1, -1, -1 )
*      IF( NB.LE.1 ) NB = MAX( 1, N )
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            CNAME, SNAME
      CHARACTER*1        C1
      CHARACTER*2        C2, C4
      CHARACTER*3        C3
      CHARACTER*6        SUBNAM
      INTEGER            I, IC, IZ, NB, NBMIN, NX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          CHAR, ICHAR, INT, MIN, REAL
*     ..
*     .. Executable Statements ..
*
      GO TO ( 100, 100, 100, 400, 500, 600, 700, 800 ) ISPEC
*
*     Invalid value for ISPEC
*
      ILAENV = -1
      RETURN
*
  100 CONTINUE
*
*     Convert NAME to upper case if the first character is lower case.
*
      ILAENV = 1
      SUBNAM = NAME
      IC = ICHAR( SUBNAM( 1:1 ) )
      IZ = ICHAR( 'Z' )
      IF( IZ.EQ.90 .OR. IZ.EQ.122 ) THEN
*
*        ASCII character set
*
         IF( IC.GE.97 .AND. IC.LE.122 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 10 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.97 .AND. IC.LE.122 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   10       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.233 .OR. IZ.EQ.169 ) THEN
*
*        EBCDIC character set
*
         IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $       ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $       ( IC.GE.162 .AND. IC.LE.169 ) ) THEN
            SUBNAM( 1:1 ) = CHAR( IC+64 )
            DO 20 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( ( IC.GE.129 .AND. IC.LE.137 ) .OR.
     $             ( IC.GE.145 .AND. IC.LE.153 ) .OR.
     $             ( IC.GE.162 .AND. IC.LE.169 ) )
     $            SUBNAM( I:I ) = CHAR( IC+64 )
   20       CONTINUE
         END IF
*
      ELSE IF( IZ.EQ.218 .OR. IZ.EQ.250 ) THEN
*
*        Prime machines:  ASCII+128
*
         IF( IC.GE.225 .AND. IC.LE.250 ) THEN
            SUBNAM( 1:1 ) = CHAR( IC-32 )
            DO 30 I = 2, 6
               IC = ICHAR( SUBNAM( I:I ) )
               IF( IC.GE.225 .AND. IC.LE.250 )
     $            SUBNAM( I:I ) = CHAR( IC-32 )
   30       CONTINUE
         END IF
      END IF
*
      C1 = SUBNAM( 1:1 )
      SNAME = C1.EQ.'S' .OR. C1.EQ.'D'
      CNAME = C1.EQ.'C' .OR. C1.EQ.'Z'
      IF( .NOT.( CNAME .OR. SNAME ) )
     $   RETURN
      C2 = SUBNAM( 2:3 )
      C3 = SUBNAM( 4:6 )
      C4 = C3( 2:3 )
*
      GO TO ( 110, 200, 300 ) ISPEC
*
  110 CONTINUE
*
*     ISPEC = 1:  block size
*
*     In these examples, separate code is provided for setting NB for
*     real and complex.  We assume that NB will take the same value in
*     single or double precision.
*
      NB = 1
*
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $            C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NB = 32
            ELSE
               NB = 32
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'PO' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( SNAME .AND. C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            NB = 64
         ELSE IF( C3.EQ.'TRD' ) THEN
            NB = 1
         ELSE IF( C3.EQ.'GST' ) THEN
            NB = 64
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NB = 32
            END IF
         END IF
      ELSE IF( C2.EQ.'GB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N4.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'PB' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            ELSE
               IF( N2.LE.64 ) THEN
                  NB = 1
               ELSE
                  NB = 32
               END IF
            END IF
         END IF
      ELSE IF( C2.EQ.'TR' ) THEN
         IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( C2.EQ.'LA' ) THEN
         IF( C3.EQ.'UUM' ) THEN
            IF( SNAME ) THEN
               NB = 64
            ELSE
               NB = 64
            END IF
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'ST' ) THEN
         IF( C3.EQ.'EBZ' ) THEN
            NB = 1
         END IF
      END IF
      ILAENV = NB
      RETURN
*
  200 CONTINUE
*
*     ISPEC = 2:  minimum block size
*
      NBMIN = 2
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( C3.EQ.'TRI' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( C3.EQ.'TRF' ) THEN
            IF( SNAME ) THEN
               NBMIN = 2
            ELSE
               NBMIN = 2
            END IF
         ELSE IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NBMIN = 2
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         ELSE IF( C3( 1:1 ).EQ.'M' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NBMIN = 2
            END IF
         END IF
      END IF
      ILAENV = NBMIN
      RETURN
*
  300 CONTINUE
*
*     ISPEC = 3:  crossover point
*
      NX = 0
      IF( C2.EQ.'GE' ) THEN
         IF( C3.EQ.'QRF' .OR. C3.EQ.'RQF' .OR. C3.EQ.'LQF' .OR.
     $       C3.EQ.'QLF' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'HRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         ELSE IF( C3.EQ.'BRD' ) THEN
            IF( SNAME ) THEN
               NX = 128
            ELSE
               NX = 128
            END IF
         END IF
      ELSE IF( C2.EQ.'SY' ) THEN
         IF( SNAME .AND. C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'HE' ) THEN
         IF( C3.EQ.'TRD' ) THEN
            NX = 1
         END IF
      ELSE IF( SNAME .AND. C2.EQ.'OR' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      ELSE IF( CNAME .AND. C2.EQ.'UN' ) THEN
         IF( C3( 1:1 ).EQ.'G' ) THEN
            IF( C4.EQ.'QR' .OR. C4.EQ.'RQ' .OR. C4.EQ.'LQ' .OR.
     $          C4.EQ.'QL' .OR. C4.EQ.'HR' .OR. C4.EQ.'TR' .OR.
     $          C4.EQ.'BR' ) THEN
               NX = 128
            END IF
         END IF
      END IF
      ILAENV = NX
      RETURN
*
  400 CONTINUE
*
*     ISPEC = 4:  number of shifts (used by xHSEQR)
*
      ILAENV = 6
      RETURN
*
  500 CONTINUE
*
*     ISPEC = 5:  minimum column dimension (not used)
*
      ILAENV = 2
      RETURN
*
  600 CONTINUE
*
*     ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
*
      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )
      RETURN
*
  700 CONTINUE
*
*     ISPEC = 7:  number of processors (not used)
*
      ILAENV = 1
      RETURN
*
  800 CONTINUE
*
*     ISPEC = 8:  crossover point for multishift (used by xHSEQR)
*
      ILAENV = 50
      RETURN
*
*     End of ILAENV
*
      END
      LOGICAL          FUNCTION LSAME( CA, CB )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          CA, CB
*     ..
*
*  Purpose
*  =======
*
*  LSAME returns .TRUE. if CA is the same letter as CB regardless of
*  case.
*
*  Arguments
*  =========
*
*  CA      (input) CHARACTER*1
*  CB      (input) CHARACTER*1
*          CA and CB specify the single characters to be compared.
*
*     .. Intrinsic Functions ..
      INTRINSIC          ICHAR
*     ..
*     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
*     ..
*     .. Executable Statements ..
*
*     Test if the characters are equal
*
      LSAME = CA.EQ.CB
      IF( LSAME )
     $   RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      ZCODE = ICHAR( 'Z' )
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
*
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
*
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.
     $       INTA.GE.145 .AND. INTA.LE.153 .OR.
     $       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.
     $       INTB.GE.145 .AND. INTB.LE.153 .OR.
     $       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
*
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
*
*     RETURN
*
*     End of LSAME
*
      END
      SUBROUTINE XERBLA ( SRNAME, INFO )
*     ..    Scalar Arguments ..
      INTEGER            INFO
      CHARACTER*6        SRNAME
*     ..
*
*  Purpose
*  =======
*
*  XERBLA  is an error handler for the Level 2 BLAS routines.
*
*  It is called by the Level 2 BLAS routines if an input parameter is
*  invalid.
*
*  Installers should consider modifying the STOP statement in order to
*  call system-specific exception-handling facilities.
*
*  Parameters
*  ==========
*
*  SRNAME - CHARACTER*6.
*           On entry, SRNAME specifies the name of the routine which
*           called XERBLA.
*
*  INFO   - INTEGER.
*           On entry, INFO specifies the position of the invalid
*           parameter in the parameter-list of the calling routine.
*
*
*  Auxiliary routine for Level 2 Blas.
*
*  Written on 20-July-1986.
*
*     .. Executable Statements ..
*
      WRITE (*,99999) SRNAME, INFO
*
      STOP
*
99999 FORMAT ( ' ** On entry to ', A6, ' parameter number ', I2,
     $         ' had an illegal value' )
*
*     End of XERBLA.
*
      END
      SUBROUTINE ADDCOLSQ( M, N, A, LDA, Q, LDQ, K, P, TAU, WORK, INFO)
*
*     Craig Lucas, University of Manchester
*     March, 2004
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LDQ, M, N, P
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), Q( LDQ, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ADDCOLSQ generates an m by m real matrix Q with orthogonal columns,
*  which is defined as the product of Q_B, elementary  reflectors and 
*  Givens rotations
*
*     Q = Q_B * H(k) * H(k+1) *...* H(k+p-1) * G(k+p-1,k+p) *...
*         *G(k,k+1) * G(k+p,k+p+1) *...* G(k+2p-2,k+2p-1)
*
*  where the H(j) and G(i,j) are as returned by ADDCOLS, such that
*  C = Q * R and C is the matrix B = Q_B * R_B, with p columns added
*  from the kth column onwards.
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
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the elements below the diagonal in columns
*          K:K+P-1 (if M > M-P+1) must contain the vector which defines 
*          the elementary reflector H(J). The elements above these 
*          vectors and below the diagonal store the scalars such that 
*          the Givens rotations can be constructed, as returned by
*          ADDCOLS.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  Q       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the matrix Q_B.
*          On exit, the matrix Q.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.  LDQ >= M.
*
*  K       (input) INTEGER
*          The postion of first column added to B.  
*          0 < K <= N-P+1.
* 
*  P       (input) INTEGER
*          The number columns added.  P > 0.
*
*  TAU     (output) DOUBLE PRECISION array, dimension(N-K+1)
*          The scalar factors of the elementary reflectors.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (2*N)
*          Work space.
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -I, the I-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   DTEMP
      INTEGER            COL, I, INC, ISTART, J, JSTOP
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARF, DLASR, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
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
      ELSE IF( K.GT.N-P+1 .OR. K.LE.0 ) THEN
         INFO = -5
      ELSE IF( P.LE.0 ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ADDCLQ', -INFO )
         RETURN
      END IF
*
*     We did a QR factorization on rows below N-P+1
*
      IF( M.GT.N-P+1 ) THEN
*
         COL = N - P + 1
         DO 10 J = K, K + P - 1
*
            DTEMP = A( COL, J )
            A( COL, J ) = ONE
*
*           If  N+P > M-N we have only factored the first M-N columns.
*
            IF( M-COL+1.LE.0 )
     $         GO TO 10
            CALL DLARF( 'R', M, M-COL+1, A( COL, J ), 1, TAU( J-K+1 ),
     $                  Q( 1, COL ), LDQ, WORK )
*
            A( COL, J ) = DTEMP
            COL = COL + 1
*
   10    CONTINUE
      END IF
*
*     If K not equal to number of columns in B then there was
*     some elimination by Givens
*
      IF( K+P-1.LT.N .AND. K.LE.M-1 ) THEN
*
*        Allow for M < N, i.e DO P wide unless hit the bottom first
*
         JSTOP = MIN( P+K-1, M-1 )
         DO 30 J = K, JSTOP
*
            ISTART = MIN( N-P+J-K+1, M )
            INC = ISTART - J
*
*           Compute vectors of C and S for rotations
*
            DO 20 I = ISTART, J + 1, -1
*
               IF( A( I, J ).EQ.ONE ) THEN
                  WORK( INC ) = ZERO
                  WORK( N+INC ) = ONE
               ELSE IF( ABS( A( I, J ) ).LT.ONE ) THEN
                  WORK( N+INC ) = A( I, J )
                  WORK( INC ) = SQRT( ( 1-A( I, J )**2 ) )
               ELSE
                  WORK( INC ) = ONE / A( I, J )
                  WORK( N+INC ) = SQRT( ( 1-WORK( INC )**2 ) )
               END IF
               INC = INC - 1
   20       CONTINUE
*
*           Apply rotations to the Jth column from the right
*
            CALL DLASR( 'R', 'V', 'b', M, ISTART-I+1, WORK( 1 ),
     $                  WORK( N+1 ), Q( 1, I ), LDQ )
*
   30    CONTINUE
*
      END IF
      RETURN
*
*     End of ADDCOLS
*
      END

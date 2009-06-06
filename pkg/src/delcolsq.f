      SUBROUTINE DELCOLSQ( M, N, A, LDA, Q, LDQ, K, P, TAU, WORK, INFO )
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
*  DELCOLSQ generates an m by m real matrix Q with orthogonal columns,
*  which is defined as the product of Q_B and elementary reflectors 
*
*        Q = Q_B * H(k) * H(k+1) *...* H(last), last = min( m-1, n ) .
*
*  where the H(j) are as returned by DELCOLSQ, such that C = Q * R and
*  C is the matrix B = Q_B * R_B, with p columns deleted from the
*  kth column onwards.
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
*          On entry, the elements below the diagonal in columns k:n
*          must contain the vector which defines the elementary 
*          reflector H(J) as returned by DELCOLS.
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
*          The position of the first column deleted from B.  
*          0 < K <= N+P.
*
*  P       (input) INTEGER
*          The number of columns deleted from B.  P > 0.
*
*  TAU     (input) DOUBLE PRECISION array, dimension(N-K+1)
*          TAU(J) must contain the scalar factor of the elementary
*          reflector H(J), as returned by DELCOLS.
*
*  WORK    DOUBLE PRECISION array, dimension (P+1)
*          Work space.
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -I, the I-th argument had an illegal value.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   AJJ
      INTEGER            J, LAST, LENH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARF, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
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
      ELSE IF( K.GT.N+P .OR. K.LE.0 ) THEN
         INFO = -5
      ELSE IF( P.LE.0 ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DELCOLSQ', -INFO )
         RETURN
      END IF
*
      LAST = MIN( M-1, N )
*
      DO 10 J = K, LAST
*
         LENH = MIN( P+1, M-J+1 )
*
*        Apply H(J) from right
*
         AJJ = A( J, J )
         A( J, J ) = ONE
*
         CALL DLARF( 'R', M, LENH, A( J, J ), 1, TAU( J-K+1 ),
     $               Q( 1, J ), LDQ, WORK )
*
         A( J, J ) = AJJ
*
   10 CONTINUE
*
      RETURN
*
*     End of DELCOLSQ
*
      END

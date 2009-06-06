      SUBROUTINE DELCOLS( M, N, A, LDA, K, P, TAU, WORK, INFO )
*
*     Craig Lucas, University of Manchester
*     March, 2004
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N, P
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  Given a real m by (n+p) matrix, B, and the QR factorization
*  B = Q_B * R_B, DELCOLS computes the QR factorization
*  C = Q * R where C is the matrix B with p columns deleted
*  from the kth column onwards.
*
*  The input to this routine is Q_B' * C
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix C.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the matrix Q_B' * C. The elements in columns
*          1:K-1 are not referenced.
*
*          On exit, the elements on and above the diagonal contain
*          the n by n upper triangular part of the matrix R. The
*          elements below the diagonal in columns k:n, together with
*          TAU represent the orthogonal matrix Q as a product of
*          elementary reflectors (see Further Details).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  K       (input) INTEGER
*          The position of the first column deleted from B. 
*          0 < K <= N+P.
*
*  P       (input) INTEGER
*          The number of columns deleted from B.  P > 0.
*
*  TAU     (output) DOUBLE PRECISION array, dimension(N-K+1)
*          The scalar factors of the elementary reflectors
*          (see Further Details).
*
*  WORK    DOUBLE PRECISION array, dimension (P+1)
*          Work space.
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -I, the I-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of Q_B and elementary 
*  reflectors
*
*     Q = Q_B * H(k) * H(k+1) *...* H(last), last = min( m-1, n ).
*
*  Each H(j) has the form
*
*     H(j) = I - tau*v*v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:j-1) = 0, v(j) = 1, v(j+1:j+lenh-1), lenh = min( p+1, m-j+1 ),
*  stored on exit in A(j+1:j+lenh-1,j) and v(j+lenh:m) = 0, tau is
*  stored in  TAU(j).
*
*  The matrix Q can be formed with DELCOLSQ
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
      EXTERNAL           DLARF, DLARFG, XERBLA
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
         CALL XERBLA( 'DELCOLS', -INFO )
         RETURN
      END IF
*
      LAST = MIN( M-1, N )
*
      DO 10 J = K, LAST
*
*        Generate elementary reflector H(J) to annihilate the nonzero
*        entries below A(J,J)
*
         LENH = MIN( P+1, M-J+1 )
         CALL DLARFG( LENH, A( J, J ), A( J+1, J ), 1, TAU( J-K+1 ) )
*
         IF( J.LT.N ) THEN
*
*           Apply H(J) to trailing matrix from left
*
            AJJ = A( J, J )
            A( J, J ) = ONE
            CALL DLARF( 'L', LENH, N-J, A( J, J ), 1, TAU( J-K+1 ),
     $                  A( J, J+1 ), LDA, WORK )
            A( J, J ) = AJJ
*
         END IF
*
   10 CONTINUE
*
      RETURN
*
*     End of DELCOLS
*
      END

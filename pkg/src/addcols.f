      SUBROUTINE ADDCOLS( M, N, A, LDA, K, P, TAU, WORK, LWORK, INFO )
*
*     Craig Lucas, University of Manchester
*     March, 2004
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N, P
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  Given a real m by (n-p) matrix, B, and the QR factorization 
*  B = Q_B * R_B, ADDCOLS computes the QR factorization
*  C = Q * R where C is the matrix B with p columns added
*  in the kth column onwards.
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
*          elements below the diagonal in columns K:N, together with
*          TAU represent the orthogonal matrix Q as a product of
*          elementary reflectors and Givens rotations. 
*          (see Further Details).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  K       (input) INTEGER
*          The position of the first column added to B.  
*          0 < K <= N-P+1.
*
*  P       (input) INTEGER
*          The number of columns added to B.  P > 0.
*
*  TAU     (output) DOUBLE PRECISION array, dimension(P)
*          The scalar factors of the elementary reflectors
*          (see Further Details).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension ( LWORK )
*          Work space.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= P.
*          For optimal performance LWORK >= P*NB, where NB is the
*          optimal block size.
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -I, the I-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of Q_B, elementary 
*  reflectors and Givens rotations
*
*     Q = Q_B * H(k) * H(k+1) *...* H(k+p-1) * G(k+p-1,k+p) *...
*         *G(k,k+1) * G(k+p,k+p+1) *...* G(k+2p-2,k+2p-1)
*  
*  Each H(j) has the form
*
*     H(j) = I - tau*v*v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:n-p-j+1) = 0, v(j) = 1, and v(j+1:m) stored on exit in 
*  A(j+1:m,j), tau is stored in TAU(j).
*
*  Each G(i,j) has the form
*
*                  i-1  i
*              [ I          ]
*              [   c   -s   ] i-1
*     G(i,j) = [   s    c   ] i
*              [          I ] 
*
*  and zero A(i,j), where c and s are encoded in scalar and 
*  stored in A(i,j) and
*     
*     IF A(i,j) = 1, c = 0, s = 1
*     ELSE IF | A(i,j) | < 1, s = A(i,j), c = sqrt(1-s**2)
*     ELSE c = 1 / A(i,j), s = sqrt(1-c**2)
*
*  The matrix Q can be formed with ADDCOLSQ
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION   C, S
      INTEGER            I, INC, ISTART, J, JSTOP, UPLEN
*     ..
*     .. External Subroutines ..
      EXTERNAL           DGEQRF, DLASR, DROT, DROTG, XERBLA
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
      ELSE IF( K.GT.N-P+1 .OR. K.LE.0 ) THEN
         INFO = -5
      ELSE IF( P.LE.0 ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ADDCOLS', -INFO )
         RETURN
      END IF
*
*     Do a QR factorization on rows below N-P, if there is more than one
*
      IF( M.GT.N-P+1 ) THEN
*
*        Level 3 QR factorization
*
         CALL DGEQRF( M-N+P, P, A( N-P+1, K ), LDA, TAU, WORK, LWORK,
     $                INFO )
*
      END IF
*
*     If K not equal to number of columns in B and not <= M-1 then 
*     there is some elimination by Givens to do
*
      IF( K+P-1.NE.N .AND. K.LE.M-1) THEN
*
*        Zero out the rest with Givens
*        Allow for M < N
*
         JSTOP = MIN( P+K-1, M-1 )
         DO 20 J = K, JSTOP
*
*           Allow for M < N
*
            ISTART = MIN( N-P+J-K+1, M )
            UPLEN = N - K - P - ISTART + J + 1
*
            INC = ISTART - J
*
            DO 10 I = ISTART, J + 1, -1
*
*              Recall DROTG updates A( I-1, J ) and 
*              stores C and S encoded as scalar in A( I, J )
*
               CALL DROTG( A( I-1, J ), A( I, J ), C, S )
               WORK( INC ) = C
               WORK( N+INC ) = S
*
*              Update nonzero rows of R
*              Do the next two line this way round because
*              A( I-1, N-UPLEN+1 ) gets updated
*
               A( I, N-UPLEN ) = -S*A( I-1, N-UPLEN )
               A( I-1, N-UPLEN ) = C*A( I-1, N-UPLEN )
*
               CALL DROT( UPLEN, A( I-1, N-UPLEN+1 ), LDA,
     $                    A( I, N-UPLEN+1 ), LDA, C, S )
*
               UPLEN = UPLEN + 1
               INC = INC - 1
*
   10       CONTINUE
*
*           Update inserted columns in one go
*           Max number of rotations is N-1, we've allowed N
*
            IF( J.LT.P+K-1 ) THEN
*
           
               CALL DLASR( 'L', 'V', 'B', ISTART-J+1, K+P-1-J,
     $                     WORK( 1 ), WORK( N+1 ), A( J, J+1 ), LDA )
*
            END IF
*
   20    CONTINUE
*
      END IF
      RETURN
*
*     End of ADDCOLS
*
      END

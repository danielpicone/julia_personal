void C_Richardson(int P, int Q, double U[Q+1][P+1], double r[Q-1][P-1], 
                  int nu, double omega, double Dx, double Dy, double b[Q-1][P-1]) 
{
    double rdx2, rdy2, resid;
    int k, p, q;

    rdx2 = 1 / ( Dx*Dx );
    rdy2 = 1 / ( Dy*Dy );
    for (k=0; k<nu; k++) {
        for (q=1; q<Q; q++) {
            for (p=1; p<P; p++) {
                r[q-1][p-1] = b[q-1][p-1] - (
                          rdx2 * ( -U[q][p+1] + 2*U[q][p] - U[q][p-1] )
                        + rdy2 * ( -U[q+1][p] + 2*U[q][p] - U[q-1][p] ) );
                U[q][p] += omega * r[q-1][p-1];
            }
        }
    }
}

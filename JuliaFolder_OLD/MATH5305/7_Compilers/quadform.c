double quadform(int n, double A[n][n], double v[n])
{
    double ans, inner_sum;
    int i, j;

    ans = 0.0;
    for (i=0; i<n; i++) {
        inner_sum = 0.0;
        for (j=0; j<n; j++) {
            inner_sum += A[i][j] * v[j];
        }
        ans += v[i] * inner_sum;
    }
    return ans;
}

void cblas_daxpy(int n, double da, double* dx, int incx, 
		                   double* dy, int incy)
{
    int ix, iy, i;

    ix = 0;
    iy = 0;
    for(i=0; i < n; i++) {
        dy[iy] += da * dx[ix];
        ix += incx;
        iy += incy;
    }
}


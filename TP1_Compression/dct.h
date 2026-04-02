void gen_dct(double **c, int N) ;
void dct1dim(double *v, double *y, int N);
void dct1dim_inv(double *vt, double *vrec, int N);
void dct2dim(double **x, double **y, int M, int N);
void dct2dim_inv(double **xt, double **xrec, int M, int N);

void dct2dim_bloc(double **x, double **y, int M, int N, int Bx, int By, int step) ; 
void dct2dim_bloc_inv(double **x, double **y, int M, int N, int Bx, int By);

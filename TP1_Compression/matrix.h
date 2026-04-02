//operations avec matrices


typedef struct {
  int x;
  int y;  } COORD; 

void mult_mat(double **a,double **b,double **c,
unsigned int ma,unsigned int na,unsigned int mb,unsigned int nb);
void mult_mvd(double **a,double *b,double *c,
unsigned int ma,unsigned int na,unsigned int mb);
void add_mat(double **a,double **b,double **c,
	     unsigned int m,unsigned int n);
void add_mat_ct(double **a,double ct,double **c,
		unsigned int m,unsigned int n);
void sub_mat(double **a,double **b,double **c,
	     unsigned int m,unsigned int n);
void sub_mat_i(int **a,int **b,int **c,
	     unsigned int m,unsigned int n);
void ch_sig_mat(double **a,double **b,
unsigned int m,unsigned int n);
int min_mi(int **x,unsigned int M,unsigned int N);
int max_mi(int **x,unsigned int M,unsigned int N);
void zero_vd(double *h,int n);
void zero_vi(int *h,int n);
void zero_vuc(unsigned char *h,int n);
void abs_mati(int **x, int **y, unsigned int m, unsigned int n);

//transposition d'une matrice

void transmd(double **a,unsigned int m,unsigned int n);
void transmc(char **a,unsigned int m,unsigned int n);
void transmi(int **a,unsigned int m,unsigned int n);
void transml(long **a,unsigned int m,unsigned int n);
void transmf(float **a,unsigned int m,unsigned int n);

void get_col(int **x,int *v, int M, int indx_col);
void put_col(int **x,int *v, int M, int indx_col);

// egalite des matrices
void egal_dc(double **a,char **b,unsigned int m,unsigned int n);
void egal_cd(char **a,double **b,unsigned int m,unsigned int n);
void egal_ucd(unsigned char **a,double **b,unsigned int m,unsigned int n);
void egal_duc(double **a,unsigned char **b,unsigned int m,unsigned int n);
void egal_ic(int **a,char **b,unsigned int m,unsigned int n);
void egal_id(int **a,double **b,unsigned int m,unsigned int n);
void egal_ci(char **a,int **b,unsigned int m,unsigned int n);
void egal_uci(unsigned char **a,int **b,unsigned int m,unsigned int n);
void egal_iuc(int **a,unsigned char **b,unsigned int m,unsigned int n);
void eg_dd(double **a,double **b,unsigned int m,unsigned int n);
void eg_cc(char **a,char **b,unsigned int m,unsigned int n);
void eg_ucuc(unsigned char **a,unsigned char **b,unsigned int m,unsigned int n);
void eg_ff(float **a,float **b,unsigned int m,unsigned int n);
void eg_ii(int **a,int **b,unsigned int m,unsigned int n);
void eg_ll(long **a,long **b,unsigned int m,unsigned int n);

//egalite des vecteurs:
void eg_dd_v(double *v,double *u,unsigned int n);

// allocation des matrices
double ** alocamd(unsigned int dx,unsigned int dy);
char ** alocamc(unsigned int dx,unsigned int dy);
signed char ** alocamsc(unsigned int dx,unsigned int dy);
unsigned char ** alocamuc(unsigned int dx,unsigned int dy);
int ** alocami(unsigned int dx,unsigned int dy);
long ** alocaml(unsigned int dx,unsigned int dy);
float ** alocamf(unsigned int dx,unsigned int dy);

// allocations des vecteurs
double *alocavd(unsigned int n);
char *alocavc(unsigned int n);
int *alocavi(unsigned int n);
long *alocavl(unsigned int n);
float *alocavf(unsigned int n);

// des-allocation des matrices
void dalocd(double ** m,unsigned int dx);
void dalocc(char ** m,unsigned int dx);
void dalocsc(signed char ** m,unsigned int dx);
void dalocuc(unsigned char ** m,unsigned int dx);
void daloci(int ** m,unsigned int dx);
void dalocl(long ** m,unsigned int dx);
void dalocf(float ** m,unsigned int dx);
double **** alocam_md(unsigned int dx,unsigned int dy, unsigned int mdx,
		      unsigned int mdy);
double *** alocam_vd(unsigned int dx,unsigned int dy, unsigned int nr_el);


// des-allocations des vecteurs
void freevd(double *v);
void freevc(char *v);
void freevi(int *v);
void freevl(long *v);
void freevf(float *v);
void daloc_md(double **** m,unsigned int dx,unsigned int dy, unsigned int mdx);
void daloc_vd(double *** m,unsigned int dx, unsigned int dy);

void afis_matd(double **a,unsigned int m, unsigned int n);


void concat_vect(int *v1, int N1, int *v2, int N2);
void concat_vectf(float *v1, int N1, float *v2, int N2);
COORD *alocavcoord(unsigned int N);
void freevcoord(COORD *c);

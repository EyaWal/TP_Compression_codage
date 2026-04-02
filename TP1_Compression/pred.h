int quantiz(double v,int step);

double calc_entropie(int **x,unsigned int M,unsigned int N);
double entropie(double *h, int nr); 
double*  calcprob(int **x,unsigned int M, unsigned int N,int *nr); 

int codeur_adapt(unsigned char **x, int **err, int H, int W, int step); 
int decodeur_adapt(int **err, unsigned char **xrec, int H, int W); 
int codeur_DPCM(unsigned char **x, int **err, int H, int W, int step); 
int decodeur_DPCM(int **err, unsigned char **xrec, int H, int W); 

int codeur_adapt_forward(unsigned char **x, int **err, int H, int W, int step); 
int codeur_DPCM_forward(unsigned char **x, int **err, int H, int W, int step); 

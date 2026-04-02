#define NBANDES 3



int lecture_dim(const char *nom, int *N, int *M);
void lecture_pgm(char *nom, unsigned char **x);
void ecriture_pgm(char *nom, unsigned char **x, int N, int M);
void lecture_ppm(char *nom, unsigned char **r, unsigned char **g, unsigned char **b);
void ecriture_ppm(char *nom, unsigned char **r, unsigned char **g, unsigned char **b, int N, int M);

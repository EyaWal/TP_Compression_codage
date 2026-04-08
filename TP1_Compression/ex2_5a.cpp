#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fichiers.h"
#include "matrix.h"
#include "pred.h"
#include "dct.h"

// ============================================================
// Utilitaires sauvegarde image entière en PGM
// ============================================================

int SaveIntImage_pgm(char *nom, int **im, int Height, int Width)
{
 int i, j;
 unsigned char **ima_uc = alocamuc(Height, Width);
 int maxim = 0;
 for (i = 0; i < Height; i++)
 for (j = 0; j < Width; j++)
 if (abs(im[i][j]) > maxim) maxim = abs(im[i][j]);
 if (maxim == 0) maxim = 1;
 for (i = 0; i < Height; i++)
 for (j = 0; j < Width; j++)
 ima_uc[i][j] = (unsigned char)((abs(im[i][j]) * 255) / maxim);
 ecriture_pgm(nom, ima_uc, Width, Height);
 dalocuc(ima_uc, Height);
 return 1;
}

int SaveIntImage_pgm_tronc(char *nom, int **im, int Height, int Width)
{
 int i, j;
 unsigned char **ima_uc = alocamuc(Height, Width);
 for (i = 0; i < Height; i++)
 for (j = 0; j < Width; j++)
 ima_uc[i][j] = (abs(im[i][j]) >= 255) ? 255 : (unsigned char)abs(im[i][j]);
 ecriture_pgm(nom, ima_uc, Width, Height);
 dalocuc(ima_uc, Height);
 return 1;
}

// ============================================================
// Exercice 1a – DPCM avec boucle de rétroaction (closed-loop)
//
// Prédicteur : x_hat(i,j) = xrec(i,j-1)
// Le codeur reconstruit en temps réel et utilise xrec comme
// prédicteur, exactement comme le décodeur => pas de dérive.
// ============================================================

int my_codeurDPCM(unsigned char **x, int **err, int H, int W, int step)
{
 int i, j, pred, xrec_val, e;
 for (i = 0; i < H; i++)
 {
 // Premier pixel de la ligne : prédicteur = 128 (DC médian)
 pred = 128;
 e = quantiz((double)((int)x[i][0] - pred), step);
 err[i][0] = e;
 xrec_val = pred + e;
 if (xrec_val < 0) xrec_val = 0;
 if (xrec_val > 255) xrec_val = 255;

 for (j = 1; j < W; j++)
 {
 pred = xrec_val; // pixel reconstruit (boucle fermée)
 e = quantiz((double)((int)x[i][j] - pred), step);
 err[i][j] = e;
 xrec_val = pred + e;
 if (xrec_val < 0) xrec_val = 0;
 if (xrec_val > 255) xrec_val = 255;
 }
 }
 return 1;
}

int my_decodeurDPCM(int **err, unsigned char **xrec, int H, int W)
{
 int i, j, pred, xrec_val;
 for (i = 0; i < H; i++)
 {
 pred = 128;
 xrec_val = pred + err[i][0];
 if (xrec_val < 0) xrec_val = 0;
 if (xrec_val > 255) xrec_val = 255;
 xrec[i][0] = (unsigned char)xrec_val;

 for (j = 1; j < W; j++)
 {
 pred = xrec_val;
 xrec_val = pred + err[i][j];
 if (xrec_val < 0) xrec_val = 0;
 if (xrec_val > 255) xrec_val = 255;
 xrec[i][j] = (unsigned char)xrec_val;
 }
 }
 return 1;
}

// ============================================================
// Exercice 1b – DPCM sans boucle de rétroaction (open-loop)
//
// Prédicteur codeur : x_hat(i,j) = x_ORIGINAL(i,j-1)
// Le décodeur ne dispose que de ses valeurs reconstruites
// => erreurs de quantification s'accumulent différemment.
// ============================================================

int my_codeurDPCM_forward(unsigned char **x, int **err, int H, int W, int step)
{
 int i, j;
 for (i = 0; i < H; i++)
 {
 err[i][0] = quantiz((double)((int)x[i][0] - 128), step);
 for (j = 1; j < W; j++)
 {
 int pred = (int)x[i][j - 1]; // pixel ORIGINAL
 err[i][j] = quantiz((double)((int)x[i][j] - pred), step);
 }
 }
 return 1;
}

int my_decodeurDPCM_forward(int **err, unsigned char **xrec, int H, int W)
{
 int i, j, xrec_val;
 for (i = 0; i < H; i++)
 {
 xrec_val = 128 + err[i][0];
 if (xrec_val < 0) xrec_val = 0;
 if (xrec_val > 255) xrec_val = 255;
 xrec[i][0] = (unsigned char)xrec_val;

 for (j = 1; j < W; j++)
 {
 xrec_val = (int)xrec[i][j - 1] + err[i][j];
 if (xrec_val < 0) xrec_val = 0;
 if (xrec_val > 255) xrec_val = 255;
 xrec[i][j] = (unsigned char)xrec_val;
 }
 }
 return 1;
}

// ============================================================
// Exercice 2 – Prédiction adaptative (closed-loop)
//
// Voisins (raster-scan, valeurs RECONSTRUITES) :
// B | C
// -----
// A | X
//
// X_hat = A si |B-C| <= |B-A|
// C sinon
// ============================================================

int my_codeur_adapt(unsigned char **x, int **err, int H, int W, int step)
{
 int i, j;
 int **xrec_enc = alocami(H, W); // reconstruction locale au codeur

 for (i = 0; i < H; i++)
 {
 for (j = 0; j < W; j++)
 {
 int pred;
 if (i == 0 && j == 0)
 {
 pred = 128;
 }
 else if (i == 0)
 {
 pred = xrec_enc[i][j - 1];
 }
 else if (j == 0)
 {
 pred = xrec_enc[i - 1][j];
 }
 else
 {
 int A = xrec_enc[i][j - 1];
 int B = xrec_enc[i - 1][j - 1];
 int C = xrec_enc[i - 1][j];
 pred = (abs(B - C) <= abs(B - A)) ? A : C;
 }

 int e = quantiz((double)((int)x[i][j] - pred), step);
 err[i][j] = e;
 int rv = pred + e;
 if (rv < 0) rv = 0;
 if (rv > 255) rv = 255;
 xrec_enc[i][j] = rv;
 }
 }
 daloci(xrec_enc, H);
 return 1;
}

int my_decodeur_adapt(int **err, unsigned char **xrec, int H, int W)
{
 int i, j;
 for (i = 0; i < H; i++)
 {
 for (j = 0; j < W; j++)
 {
 int pred;
 if (i == 0 && j == 0)
 {
 pred = 128;
 }
 else if (i == 0)
 {
 pred = (int)xrec[i][j - 1];
 }
 else if (j == 0)
 {
 pred = (int)xrec[i - 1][j];
 }
 else
 {
 int A = (int)xrec[i][j - 1];
 int B = (int)xrec[i - 1][j - 1];
 int C = (int)xrec[i - 1][j];
 pred = (abs(B - C) <= abs(B - A)) ? A : C;
 }

 int rv = pred + err[i][j];
 if (rv < 0) rv = 0;
 if (rv > 255) rv = 255;
 xrec[i][j] = (unsigned char)rv;
 }
 }
 return 1;
}

// ============================================================
// MAIN
// ============================================================

int main(int argc, char *argv[])
{
 char nom[200], nom_out[200], nom_err[300];
 int W, H, i, j;

 if (argc != 3)
 {
 fprintf(stderr, "Utilisation: %s <nom_image_pgm> <pas_quantification>\n", argv[0]);
 exit(0);
 }

 strcpy(nom, argv[1]);
 int step = atoi(argv[2]);

 lecture_dim(nom, &W, &H);
 fprintf(stderr, "Image : %s | Width=%d Height=%d step=%d\n\n", nom, W, H, step);

 unsigned char **x = alocamuc(H, W);
 unsigned char **xrec = alocamuc(H, W);
 int **err = alocami(H, W);

 lecture_pgm(nom, x);

 // --- Entropie de l'image originale ---
 for (i = 0; i < H; i++)
 for (j = 0; j < W; j++)
 err[i][j] = (int)x[i][j];
 double h_orig = calc_entropie(err, H, W);
 fprintf(stderr, "[ORIGINAL] entropie = %.4f bits/pixel\n\n", h_orig);

 // ==========================================================
 // Exercice 1a : DPCM avec boucle de rétroaction
 // ==========================================================
 fprintf(stderr, "--- EX1a : DPCM avec boucle de retroaction ---\n");
 my_codeurDPCM(x, err, H, W, step);
 double h1a = calc_entropie(err, H, W);
 fprintf(stderr, " entropie erreurs = %.4f bits/pixel\n", h1a);

 my_decodeurDPCM(err, xrec, H, W);
 sprintf(nom_out, "%s_dpcm_closed_step%d.pgm", nom, step);
 ecriture_pgm(nom_out, xrec, W, H);
 sprintf(nom_err, "%s_dpcm_closed_err_step%d.pgm", nom, step);
 SaveIntImage_pgm(nom_err, err, H, W);
 fprintf(stderr, " image reconstruite -> %s\n\n", nom_out);

 // ==========================================================
 // Exercice 1b : DPCM sans boucle de rétroaction
 // ==========================================================
 fprintf(stderr, "--- EX1b : DPCM sans boucle de retroaction ---\n");
 my_codeurDPCM_forward(x, err, H, W, step);
 double h1b = calc_entropie(err, H, W);
 fprintf(stderr, " entropie erreurs = %.4f bits/pixel\n", h1b);

 my_decodeurDPCM_forward(err, xrec, H, W);
 sprintf(nom_out, "%s_dpcm_open_step%d.pgm", nom, step);
 ecriture_pgm(nom_out, xrec, W, H);
 sprintf(nom_err, "%s_dpcm_open_err_step%d.pgm", nom, step);
 SaveIntImage_pgm(nom_err, err, H, W);
 fprintf(stderr, " image reconstruite -> %s\n\n", nom_out);

 // ==========================================================
 // Exercice 2 : Prédiction adaptative
 // ==========================================================
 fprintf(stderr, "--- EX2 : Prediction adaptative ---\n");
 my_codeur_adapt(x, err, H, W, step);
 double h2 = calc_entropie(err, H, W);
 fprintf(stderr, " entropie erreurs = %.4f bits/pixel\n", h2);

 my_decodeur_adapt(err, xrec, H, W);
 sprintf(nom_out, "%s_adapt_step%d.pgm", nom, step);
 ecriture_pgm(nom_out, xrec, W, H);
 sprintf(nom_err, "%s_adapt_err_step%d.pgm", nom, step);
 SaveIntImage_pgm(nom_err, err, H, W);
 fprintf(stderr, " image reconstruite -> %s\n\n", nom_out);

 // ==========================================================
 // Exercice 3 : Résumé comparatif
 // ==========================================================
 fprintf(stderr, "======================================================\n");
 fprintf(stderr, " BILAN COMPARATIF (step = %d)\n", step);
 fprintf(stderr, "======================================================\n");
 fprintf(stderr, " Original : %.4f bits/pixel\n", h_orig);
 fprintf(stderr, " DPCM avec boucle : %.4f bits/pixel\n", h1a);
 fprintf(stderr, " DPCM sans boucle : %.4f bits/pixel\n", h1b);
 fprintf(stderr, " Prediction adaptative : %.4f bits/pixel\n", h2);
 fprintf(stderr, "======================================================\n\n");

 // ==========================================================
 // Exercice 4 : DCT globale (step recommandés: 5,10,20,30)
 // ==========================================================
 {
 fprintf(stderr, "--- EX4 : DCT globale (step=%d) ---\n", step);
 double **tdct = alocamd(H, W);
 double **xd = alocamd(H, W);
 double **xrecd = alocamd(H, W);

 for (i = 0; i < H; i++)
 for (j = 0; j < W; j++)
 xd[i][j] = (double)x[i][j];

 dct2dim(xd, tdct, H, W);

 // Quantification : step pour les AC (i>0 ET j>0), pas=1 pour 1ere ligne/colonne
 for (i = 1; i < H; i++)
 for (j = 1; j < W; j++)
 { err[i][j] = quantiz(tdct[i][j], step); tdct[i][j] = (double)err[i][j]; }
 for (j = 0; j < W; j++)
 { err[0][j] = quantiz(tdct[0][j], 1); tdct[0][j] = (double)err[0][j]; }
 for (i = 0; i < H; i++)
 { err[i][0] = quantiz(tdct[i][0], 1); tdct[i][0] = (double)err[i][0]; }

 double h4 = calc_entropie(err, H, W);
 fprintf(stderr, " entropie coeff. DCT globale = %.4f bits/pixel\n", h4);

 dct2dim_inv(tdct, xrecd, H, W);
 for (i = 0; i < H; i++)
 for (j = 0; j < W; j++)
 {
 if (xrecd[i][j] < 0.0) xrec[i][j] = 0;
 else if (xrecd[i][j] > 255.0) xrec[i][j] = 255;
 else xrec[i][j] = (unsigned char)xrecd[i][j];
 }
 sprintf(nom_out, "%s_dct_step%d.pgm", nom, step);
 ecriture_pgm(nom_out, xrec, W, H);
 sprintf(nom_err, "%s_dct_err_step%d.pgm", nom, step);
 SaveIntImage_pgm_tronc(nom_err, err, H, W);
 fprintf(stderr, " image reconstruite -> %s\n\n", nom_out);

 dalocd(tdct, H); dalocd(xd, H); dalocd(xrecd, H);
 }

 // ==========================================================
 // Exercice 5a : DCT par blocs (8x8, 16x16, 32x32) – quantification uniforme
 // ==========================================================
 {
 int block_sizes[] = {8, 16, 32};
 for (int bs = 0; bs < 3; bs++)
 {
 int Bx = block_sizes[bs], By = block_sizes[bs];
 fprintf(stderr, "--- EX5a : DCT par blocs %dx%d uniforme (step=%d) ---\n", Bx, By, step);

 double **tdct = alocamd(H, W);
 double **xd = alocamd(H, W);
 double **xrecd = alocamd(H, W);

 for (i = 0; i < H; i++)
 for (j = 0; j < W; j++)
 xd[i][j] = (double)x[i][j];

 dct2dim_bloc(xd, tdct, H, W, Bx, By, step);

 for (i = 0; i < H; i++)
 for (j = 0; j < W; j++)
 err[i][j] = (int)tdct[i][j];

 fprintf(stderr, " entropie coeff. DCT blocs %dx%d = %.4f bits/pixel\n",
 Bx, By, calc_entropie(err, H, W));

 dct2dim_bloc_inv(tdct, xrecd, H, W, Bx, By);

 for (i = 0; i < H; i++)
 for (j = 0; j < W; j++)
 {
 if (xrecd[i][j] < 0.0) xrec[i][j] = 0;
 else if (xrecd[i][j] > 255.0) xrec[i][j] = 255;
 else xrec[i][j] = (unsigned char)xrecd[i][j];
 }

 sprintf(nom_out, "%s_dct_bloc%dx%d_step%d.pgm", nom, Bx, By, step);
 ecriture_pgm(nom_out, xrec, W, H);
 fprintf(stderr, " image reconstruite -> %s\n\n", nom_out);

 dalocd(tdct, H); dalocd(xd, H); dalocd(xrecd, H);
 }
 }

 // ==========================================================
 // Exercice 5a : DCT 8x8 – quantification variable (Figure 2 du TP)
 // Pas de quantification = Qmat[u][v] * delta / 8, delta in {1, 2, 10}
 // ==========================================================
 {
 static const int Qmat[8][8] = {
 { 8, 17, 18, 19, 21, 23, 25, 27},
 {17, 18, 19, 21, 23, 25, 27, 28},
 {20, 21, 22, 23, 24, 26, 28, 30},
 {21, 22, 23, 24, 26, 28, 30, 32},
 {22, 23, 24, 26, 28, 30, 32, 35},
 {23, 24, 26, 28, 30, 32, 35, 38},
 {25, 26, 28, 30, 32, 35, 38, 41},
 {27, 28, 30, 32, 35, 38, 41, 45}
 };
 int deltas[] = {1, 2, 10};

 double **xd = alocamd(H, W);
 for (i = 0; i < H; i++)
 for (j = 0; j < W; j++)
 xd[i][j] = (double)x[i][j];

 for (int di = 0; di < 3; di++)
 {
 int delta = deltas[di];
 fprintf(stderr, "--- EX5a : DCT 8x8 quantif. variable (delta=%d) ---\n", delta);

 double **xrecd = alocamd(H, W);
 double **blk = alocamd(8, 8);
 double **blkd = alocamd(8, 8);
 double **blkr = alocamd(8, 8);

 int bi, bj, u, v;
 for (bi = 0; bi + 8 <= H; bi += 8)
 {
 for (bj = 0; bj + 8 <= W; bj += 8)
 {
 for (u = 0; u < 8; u++)
 for (v = 0; v < 8; v++)
 blk[u][v] = xd[bi + u][bj + v];

 dct2dim(blk, blkd, 8, 8);

 for (u = 0; u < 8; u++)
 for (v = 0; v < 8; v++)
 {
 int qstep = (int)round(Qmat[u][v] * delta / 8.0);
 if (qstep < 1) qstep = 1;
 int q = quantiz(blkd[u][v], qstep);
 err[bi + u][bj + v] = q;
 blkd[u][v] = (double)q;
 }

 dct2dim_inv(blkd, blkr, 8, 8);

 for (u = 0; u < 8; u++)
 for (v = 0; v < 8; v++)
 xrecd[bi + u][bj + v] = blkr[u][v];
 }
 }

 fprintf(stderr, " entropie coeff. DCT 8x8 varquant delta=%d = %.4f bits/pixel\n",
 delta, calc_entropie(err, H, W));

 for (i = 0; i < H; i++)
 for (j = 0; j < W; j++)
 {
 if (xrecd[i][j] < 0.0) xrec[i][j] = 0;
 else if (xrecd[i][j] > 255.0) xrec[i][j] = 255;
 else xrec[i][j] = (unsigned char)xrecd[i][j];
 }

 sprintf(nom_out, "%s_dct_8x8_varquant_delta%d.pgm", nom, delta);
 ecriture_pgm(nom_out, xrec, W, H);
 fprintf(stderr, " image reconstruite -> %s\n\n", nom_out);

 dalocd(xrecd, H);
 dalocd(blk, 8); dalocd(blkd, 8); dalocd(blkr, 8);
 }
 dalocd(xd, H);
 }

 // Libération mémoire
 dalocuc(x, H);
 dalocuc(xrec, H);
 daloci(err, H);

 return 0;
}

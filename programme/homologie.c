#include <stdio.h>
#include <conio.h>
#include <math.h>
#include <stdlib.h>

const Max = 220,                                 // Maximale Groesse der Matrizen
                                                 // (reicht fuer 3 <= h <= 4)

// *********************************************************************************

      MatrAnzahl = 2,                            // Beispiel: h=3, c=1, Komponente 2
      Dim_der_Komponente = 3,
      Anz_Zeilen [MatrAnzahl] = {5,20},
      Anz_Spalten[MatrAnzahl] = {20,16};

char huge D [MatrAnzahl][Max][Max] = {

{
{1,0,0,0,-1,0,0,0,0,1,0,0,-1,0,0,0,0,0,0,0},     // Matrix 2 -> 1
{0,-1,0,-1,0,0,-1,0,-1,0,0,0,0,0,0,0,0,1,-1,0},
{0,0,-1,0,-1,-1,0,-1,0,0,0,0,0,0,0,0,0,-1,0,1},
{0,0,0,0,0,0,0,0,0,0,-1,0,-1,0,-1,0,-1,1,0,0},
{0,0,0,0,0,0,0,0,0,-1,0,-1,0,-1,0,-1,0,0,0,1}},
{
{0,0,0,0,1,0,0,0,0,-1,0,0,0,0,0,0},              // Matrix 3 -> 2
{0,1,0,0,0,0,0,0,0,0,-1,0,0,0,0,0},
{1,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0},
{0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,0},
{0,0,0,0,0,0,0,0,0,0,0,1,0,-1,0,0},
{0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1},
{0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0},
{-1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0},
{0,-1,0,-1,0,0,0,0,0,0,0,0,0,0,0,0},
{0,0,0,0,-1,0,0,0,0,0,0,1,0,0,0,0},
{0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0},
{0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,0,-1,0,0,0,1,0,0},
{0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1},
{0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0},
{0,0,0,0,0,-1,0,-1,0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,-1,0,-1,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0},
{0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,0},
{0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1}}
};

// *********************************************************************************

const Verschiebung = Dim_der_Komponente - MatrAnzahl + 1;

int   n, m, dim, torsion = 0;

typedef char huge Matrix[Max][Max];

                  Matrix E1, F;

int sgn(int z)
{
   if (z > 0)  return(1);
   if (z == 0) return(0);
               return(-1);
}

void ZeilenVertauschen(int i, int j)             // vertausche Zeilen i und j von D[dim]
{
   int l, zwischen;

   for (l=0; l<m; l++)
   {
      zwischen = D[dim][i][l]; D[dim][i][l] = D[dim][j][l]; D[dim][j][l] = zwischen;
   }
   for (l=0; l<n; l++)                           // vertausche Spalten j und i von E1
   {
      zwischen = E1[l][i]; E1[l][i] = E1[l][j]; E1[l][j] = zwischen;
   }
   if (dim > 0)
      for (l=0; l < Anz_Zeilen[dim-1]; l++)      // vertausche Spalten j, i von D[dim-1]
      {
         zwischen = D[dim-1][l][i];
         D[dim-1][l][i] = D[dim-1][l][j];
         D[dim-1][l][j] = zwischen;
      }
}

void SpaltenVertauschen(int i, int j)            // vertausche Spalten i und j von D[dim]
{
   int l, zwischen;

   for (l=0; l<n; l++)
   {
      zwischen = D[dim][l][i]; D[dim][l][i] = D[dim][l][j]; D[dim][l][j] = zwischen;
   }
   for (l=0; l<m; l++)
   {
      zwischen = F[l][i]; F[l][i] = F[l][j]; F[l][j] = zwischen;
   }
}

void ZeilenAddieren(int i, int j, int lambda)    // Zeile j += lambda*Zeile i in D[dim]
{
   int l;

   for (l=0; l<m; l++)  D[dim][j][l] += lambda *  D[dim][i][l];
                                                 // Spalte i -= lambda*Spalte j in E1
   for (l=0; l<n; l++)  E1[l][i] -= lambda * E1[l][j];
   if (dim > 0)                                  // Spalte i -= lambda*Spalte j in D[dim-1]
      for (l=0; l < Anz_Zeilen[dim-1]; l++)
         D[dim-1][l][i] -= lambda * D[dim-1][l][j];
}

void SpaltenAddieren(int i, int j, int lambda)   // Spalte j += lambda*Spalte i in D[dim]
{
   int l;

   for (l=0; l<n; l++)  D[dim][l][j] += lambda *  D[dim][l][i];
   for (l=0; l<m; l++)  F[l][j] += lambda *  F[l][i];
}

void SpalteInvertieren(int i)
{
   SpaltenAddieren(i,i,-2);
}

void SpalteVerschwindet(int k)                   // Spalte k besteht nur noch aus dem
{                                                // Eintrag auf der Diagonale
   int i, q, r, eps;

   for (i=k+1; i<n; i++)
      if (D[dim][i][k] != 0)
      {
         q = abs(D[dim][i][k]) / abs(D[dim][k][k]);
         r = abs(D[dim][i][k]) % abs(D[dim][k][k]);
         eps = -sgn(D[dim][k][k])*sgn(D[dim][i][k]);
         ZeilenAddieren(k, i, eps*q);
         if (r != 0)  { ZeilenVertauschen(k,i); i = k; }
      }
}

void ZeileVerschwindet(int k)                    // Zeile k besteht nur noch aus dem
{                                                // Eintrag auf der Diagonale
   int j, q, r, eps;

   for (j=k+1; j<m; j++)
      if (D[dim][k][j] != 0)
      {
         q = abs(D[dim][k][j]) / abs(D[dim][k][k]);
         r = abs(D[dim][k][j]) % abs(D[dim][k][k]);
         eps = -sgn(D[dim][k][k])*sgn(D[dim][k][j]);
         SpaltenAddieren(k, j, eps*q);
         if (r != 0)
         {
            SpaltenVertauschen(k, j);
            SpalteVerschwindet(k);
            j = k;
         }
      }
}

void Normalform(void)                            // bringt die Matrix D[dim] in
{                                                // Normalform ueber Z
   int i, j, k, min, imin, jmin, q, eps;

   if (torsion > 0)
      for (i=0; i < m-torsion; i++)  SpaltenVertauschen(i, i+torsion);

   for (k=0; ((k<n) && (k<m)); k++)
   {
      min = 0;                                   // finde minimalen Eintrag
      for (i=k; i<n; i++)
         for (j=k; j<m; j++)
            if ((abs(D[dim][i][j]) > 0) && ((abs(D[dim][i][j]) < min) || (min == 0)))
            {
               min = abs(D[dim][i][j]);
               imin = i; jmin = j;
            }
      if (min == 0) return;                      // D[dim] ist diagonalisiert

      ZeilenVertauschen(k, imin);                // bringe den kleinsten Eintrag nach
      SpaltenVertauschen(k, jmin);               // links oben
      SpalteVerschwindet(k);
      ZeileVerschwindet(k);
      for (i=k+1; i<n; i++)
         for (j=k+1; j<m; j++)
            if ((abs(D[dim][i][j]) % abs(D[dim][k][k])) != 0)
            {
               q = abs(D[dim][i][j]) / abs(D[dim][k][k]);
               eps = -sgn(D[dim][k][k])*sgn(D[dim][i][j]);
               ZeilenAddieren(k, i, eps*q);
               SpaltenAddieren(k, j, 1);
               ZeilenVertauschen(k, i);
               SpaltenVertauschen(k, j);

               SpalteVerschwindet(k);
               ZeileVerschwindet(k);
               i = k; j = m;
            }
      if (D[dim][k][k] < 0)  SpalteInvertieren(k);
   }
}

void Diagonale (void)                            // druckt die Diagonale der
{                                                // diagonalisierten Matrix D[dim]
   int i;

   for (i=0; i < torsion; i++)                   printf("* ");
   for (i=0; ((i < n) && (i < m-torsion)); i++)  printf("%d ", D[dim][i][i]);
   for (i=n; i < m-torsion; i++)                 printf("+ ");
   printf("\n");
}

void Erzeuger(int nr, int vielfach)
{
   int i;

   if (vielfach != 0)
   {
      for (i=1; i <= abs(vielfach); i++)
      {
         if (vielfach>0) printf("+"); else printf("-");
      }
      printf("%d ", nr);
   }
}

void Homologie (void)                            // ermittelt aus D[dim] in Diagonalform
{                                                // die Homologie
   int i = 0, j, k;

   printf("\n");
   while ((i<m-torsion) && (i<n) && (D[dim][i][i] != 0)) i++;

   if (i < m-torsion)                            // freier Anteil in Dimension dim
   {
      printf("Z^%d in dim %d, Erzeuger:\n", m-torsion - i, dim+Verschiebung);
      for (j=i; j < m-torsion; j++)
      {
         for (k=0; k<m; k++)  Erzeuger(k+1, F[k][j]);
         printf("\n");
      }
   }
   i = 0;                                        // Torsion in Dimension dim-1 
   while ((i<m-torsion) && (i<n) && (D[dim][i][i] != 0))
   {
      if (D[dim][i][i] > 1)
      {
         printf("Z_%d in dim %d, Erzeuger:\n", D[dim][i][i], dim-1 + Verschiebung);
         printf("d(");
         for (k=0; k<m; k++)  Erzeuger(k+1, F[k][i]);
         printf(")\n= %d*(", D[dim][i][i]);
         for (k=0; k<n; k++)  Erzeuger(k+1, E1[k][i]);
         printf(")\n");
      }
      i++;
   }
   torsion = i;
                                                 // freier Anteil in Dimension dim-1,
   if ((dim == 0) && (i<n))                      // wenn D[dim] die letzte Matrix ist
   {
      printf("Z^%d in dim %d, Erzeuger:\n", n-i, dim-1 + Verschiebung);
      for (j=i; j<n; j++)
      {
         for (k=0; k<n; k++)  Erzeuger(k+1, E1[k][j]);
         printf("\n");
      }
   }
   printf("\n");
}

void initialisiere(Matrix M, int groesse)        // M := Einheitsmatrix
{
   int i, j;

   for (i=0; i<groesse; i++)
      for (j=0; j<groesse; j++) M[i][j] = 0;
   for (i=0; i<groesse; i++)    M[i][i] = 1;
}

void main (void)
{
   int i,j;

   clrscr();
   for (dim=MatrAnzahl-1; dim>=0; dim--)
   {
      n = Anz_Zeilen [dim];
      m = Anz_Spalten[dim];
      if (dim == MatrAnzahl-1) initialisiere(F, m);
      else for (i=0; i<m; i++)
              for (j=0; j<m; j++)  F[i][j] = E1[i][j];
      initialisiere(E1, n);

      Normalform();
      Diagonale();
      Homologie();
   }
   getch();
}

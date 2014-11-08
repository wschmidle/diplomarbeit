#include <stdio.h>
#include <conio.h>
#include <math.h>
#include <stdlib.h>

const Max_Zellen = 500,                          // max. Zellenanzahl pro Dimension
                                                 // (reicht fuer h<=4)

// ********************************************************************************

               h = 3,                            // g = (h-c)/2
               c = 1,

// ********************************************************************************

       Dimension = 2*h-2-(h-c)/2;                // Dimension von P^sym(h,c)

typedef int        Nummern [Dimension+1];

typedef char Schichten     [h+1],
             Permutationen      [2*h+1],
             SchichtenPerm [h+1][2*h+1],
              huge Tabelle [Dimension+1][Max_Zellen][h+1][2*h+1],
                   Angaben [Dimension+1][Max_Zellen],
              huge Matrix  [Max_Zellen][Max_Zellen];

        int      dim, n, m, komponente;          // m = Schlitzhaufen, n = Schichten
   long int      insgesamt;                      // Gesamtzahl aller Fixpunkt-Zellen

   Nummern       nr, komp_groesse;               // Anzahlen in jeder Dimension
   Schichten     h_frei;                         // Anzahl der noch uebrigen Schlitzpaare
   Permutationen taulambda, zw, gezaehlt;
   SchichtenPerm involution, pi, pi1;            // Permutation einer Schicht und Inverse
   Tabelle       liste;                          // Liste aller Fixpunkt-Zellen
   Angaben       anz_m, anz_n, zusammen;         // zusammen: Nr. der Komponente einer Zelle
   Matrix        inzidenz;                       // Inzidenzmatrizen der Randabbildungen

void Neue_Schicht(void);
void Permutation(void);

void Init(void)
{
   int i, j, k, l;

   for (i=0; i <= Dimension; i++)
      for (j=1; j < Max_Zellen; j++)
         for (k=1; k <= h; k++)
            for (l=1; l <= 2*h; l++)  liste[i][j][k][l] = 0;
   for (i=1; i <= Dimension; i++)  nr[i] = 0;
   h_frei[0] = h;
   insgesamt = 0;
}

int Kreiszahl(void)                              // ermittelt die Kreiszahl c
{
   int j, k, num = -1;

   taulambda[n] = 1; for (j=1; j<n; j++) taulambda[j] = j+1;
   for (k=1; k<=m; k++)
   {
      for (j=1; j<=n; j++)  zw[j] = taulambda[pi[k][j]];
      for (j=1; j<=n; j++)  taulambda[j] = zw[j];
   }
   for (j=1; j<=n; j++)  gezaehlt[j] = 0;
   for (j=1; j<=n; j++)
      if (gezaehlt[j] == 0)
      {
         k = j;
         do
         {
            gezaehlt[k] = 1;
            k = taulambda[k];
         }
         while (k != j);
         num++;
      }
   return(num);
}

void Ausgeben(void)
{
   int j, k;

   printf("%ld:n=%d,m=%d,dim=%d", insgesamt, n, m, dim);
   for (j=1; j<=m; j++)
   {
      printf(" ");
      for (k=1; k<=n; k++)  printf("%d", pi[j][k]);
   }
   printf("\n");
}

void Eintragen(void)
{
   int j, k;

   nr[dim]++;
   anz_m[dim][nr[dim]] = m;
   anz_n[dim][nr[dim]] = n;
   for (j=1; j<=m; j++)
      for (k=1; k<=n; k++)  liste[dim][nr[dim]][j][k] = pi[j][k];
}

int Anz_Paare(int j, int schlitzhaufen)          // ermittelt die Anzahl der Schlitzpaare
{                                                // in Schicht j; ein Zykel der Laenge q
   int i, k, num = 0;                            // in pi[j] entspricht q-1 Schlitzpaaren

   for (i=1; i<=schlitzhaufen; i++)  gezaehlt[i] = 0;
   for (i=1; i<=schlitzhaufen; i++)
      if (gezaehlt[i] == 0)
      {
         k = i;
         do
         {
            gezaehlt[k] = 1;
            k = pi[j][k];
            num++;
         }
         while (k != i);
         num--;
      }
   return(num);
}

int Freie_Schlitzhaufen(void)                    // ermittelt die Anzahl der noch freien
{                                                // Schlitzhaufen
   int j, k, frei = 0;

   for (j=1; j<=n; j++)
      {
         k=1;
         while ((k <= m) && (pi[k][j] == j))  k++;
         if (k > m)  frei++;
      }
   return(frei);
}

int Suchen(int schichten, int schlitzhaufen)     // sucht einen Fixpunkt in der Liste
{
   int i = 1, j, k, gefunden;

   do
   {
      if ((anz_m[dim-1][i] == schichten) && (anz_n[dim-1][i] == schlitzhaufen))
      {
         gefunden = 1;
         for (j=1; j <= schichten; j++)
            for (k=1; k <= schlitzhaufen; k++)
               if (pi[j][k] != liste[dim-1][i][j][k])
               {
                  gefunden = 0;
                  j = schichten; k = schlitzhaufen;
               }
         if (gefunden == 1)  return(i);
      }
   }
   while (++i <= nr[dim-1]);
   return(0);
}

void Neue_Schicht(void)                          // erzeugt rekursiv alle Fixpunkte
{                                                // mit n Schlitzhaufen
   int j, k;

   if (m == 1)  for (j=1; j<=n; j++)  involution[1][j] = n+1 - j;
   else         for (j=1; j<=n; j++)  involution[m][j] = involution[m-1][pi[m-1][j]];
   for (j=1; j<=n; j++)
   {
      pi [m][j] = -1;
      pi1[m][j] = -1;
   }
   Permutation();
}

void Permutation(void)                           // erzeugt in der aktuellen Schicht
{                                                // rekursiv alle erlaubten Permutationen
   int j = 1, l;

   while ((j<=n) && (pi[m][j] != -1))  j++;
   if (j<=n)
   {
      for (l=1; l<=n; l++)                       // vervollstaendige die Permutation
         if (pi1[m][l] == -1)
         {
            pi [m][j] =  l; pi [m][ involution[m][l] ] = involution[m][j];
            pi1[m][l] =  j; pi1[m][ involution[m][j] ] = involution[m][l];
            Permutation();
            pi [m][j] = -1; pi [m][ involution[m][l] ] = -1;
            pi1[m][l] = -1; pi1[m][ involution[m][j] ] = -1;
         }
   }
   else
   {
      h_frei[m] = h_frei[m-1] - Anz_Paare(m, n);
      if ((h_frei[m] < h_frei[m-1]) && (Freie_Schlitzhaufen() <= 2*h_frei[m]))
      {
         if (h_frei[m] == 0)
         {
            if (Kreiszahl() == c)
            {
               insgesamt++;
               dim = (n-2)/2 + m-1;
               Eintragen();
               Ausgeben();
            }
         }
         else
         {
            m++;
            Neue_Schicht();
            m--;
         }
      }
   }
}

int D_1(int i, int verschiebung)
{
   int j, l, q, s, i_max;

   for (j=1; j <= m; j++)
      for (l=1; l <= n-verschiebung; l++)  pi1[j][pi[j][l]] = l;
   for (j=m; j>=1; j--)
      if (pi[j][i] != i)  i_max = j;
   for (j=m; j >= i_max; j--)
      if (pi[j][i+1] != i+1)
      {
         s = i;
         for (q=i_max; q<=j-1; q++)
         {
            s = pi1[q][s];
            if (s == i+1)  return (0);
         }
         q = s;
         do
         {
            q = pi[j][q];
            if (q == i+1)  return (0);
         }
         while (q != s);
         pi[j][pi1[j][s]]   = pi[j][i+1];
         pi[j][pi1[j][i+1]] = s;
         pi[j][i+1]        = i+1;
      }
   for (j=1; j<=m; j++)                          // die Reste der beiden Schlitzhaufen
   {                                             // werden zusammengeschoben
      if (pi[j][i] == i)  pi[j][i] = pi[j][i+1];
      for (q=i+1; q < n-verschiebung; q++)  pi[j][q] = pi[j][q+1];
      for (q=1; q < n-verschiebung; q++)
         if (pi[j][q] > i)  pi[j][q]--;
   }
   return(1);
}

int D_1_sym(int i, int num)
{
   int j, k, nummer;

   for (j=1; j<=m; j++)
      for (k=1; k<=n; k++)  pi[j][k] = liste[dim][num][j][k];
   nummer = D_1(n-i, 0);
   if (nummer == 0)  return(0);
   if (2*i < n)  nummer = D_1(i, 1);
   if (nummer == 0)  return(0);
   if (2*i ==n)  nummer = Suchen(m, n-1);
   else          nummer = Suchen(m, n-2);
   return(nummer);
}

int D_2(int j, int num)
{
   int k, l, finger;

   for (k=1; k<=m; k++)
      for (l=1; l<=n; l++)  pi[k][l] = liste[dim][num][k][l];
   finger = Anz_Paare(j, n) + Anz_Paare(j+1, n);
   for (l=1; l<=n; l++)  zw[l] = pi[j][pi[j+1][l]];
   for (l=1; l<=n; l++)  pi[j][l] = zw[l];
   if (Anz_Paare(j, n) < finger)  return(0);
   for (k=j+1; k<m; k++)
      for (l=1; l<=n; l++)  pi[k][l] = pi[k+1][l];
   return(Suchen(m-1, n));
}

void Verbinden(int num, int inz)
{
   int i, j, min, max;

   if (zusammen[dim-1][inz] == 0)  zusammen[dim-1][inz] = zusammen[dim][num];
   else
   {
      if (zusammen[dim-1][inz] != zusammen[dim][num])
      {
         if (zusammen[dim-1][inz] < zusammen[dim][num])
         {
            min = zusammen[dim-1][inz];
            max = zusammen[dim][num];
         }
         else
         {
            min = zusammen[dim][num];
            max = zusammen[dim-1][inz];
         }
         for (i=1; i<=dim; i++)
            for (j=1; j<=nr[dim]; j++)
            {
               if (zusammen[i][j] == max)  zusammen[i][j] = min;
               if (zusammen[i][j] >  max)  zusammen[i][j]--;
            }
         komponente--;
      }
   }
}

void Komponenten(void)                           // funktioniert nur dann richtig, wenn
{                                                // es in Dimension 0 nur den degenerierten
   int i, k, num, inz;                           // Punkt W gibt (also fuer h >= 3)

   for (dim=1; dim<=Dimension; dim++)
      for (num=1; num<=nr[dim]; num++)  zusammen[dim][num] = 0;
   komponente = 0;
   for (dim=2; dim<=Dimension; dim++)
      for (num=1; num<=nr[dim]; num++)
      {
         if (zusammen[dim][num] == 0)  zusammen[dim][num] = ++komponente;
         n = anz_n[dim][num];
         m = anz_m[dim][num];
         k = n/2;
         for (i=1; i<=k; i++)
         {
            inz = D_1_sym(i,num);
            if (inz != 0)  Verbinden(num,inz);
         }
         for (i=1; i<=m-1; i++)
         {
            inz = D_2(i,num);
            if (inz != 0)  Verbinden(num,inz);
         }
      }
}

void Komp_Groessen (void)                        // gibt die Zellenanzahlen der Komponenten
{                                                // in jeder Dimension aus
   int i, j, k, groesse;

   printf("\n%d Komponenten\n", komponente);
   for (i=1; i<=komponente; i++)
   {
      groesse = 0;
      for (j=1; j <= Dimension; j++)  komp_groesse[j] = 0;
      for (j=1; j <= Dimension; j++)
         for (k=1; k <= nr[j]; k++)
            if (zusammen[j][k] == i)
            {
               groesse++;
               komp_groesse[ anz_m[j][k]-1 + (anz_n[j][k]-2)/2 ]++;
            }
      printf("Komponente Nr. %d: %d, ", i, groesse);
      for (j=1; j <= Dimension; j++)  printf(" %d: %d, ", j, komp_groesse[j]);
      printf("\n");
   }
}

void Rand(int num)                               // bestimmt den Rand einer Zelle
{
   int i, j, k, vorzeichen = 1, inz;

   n = anz_n[dim][num];
   m = anz_m[dim][num];
   k = n/2;
   for (i=1; i<=k; i++)
   {
      inz = D_1_sym(i,num);
      if (inz != 0)  inzidenz[num][inz] += vorzeichen;
      vorzeichen *= -1;
   }
   if (k > 0)  vorzeichen *= -1;
   for (j=1; j<=m-1; j++)
   {
      inz = D_2(j,num);
      if (inz != 0)  inzidenz[num][inz] += vorzeichen;
      vorzeichen *= -1;
   }
}

void Delta(int komp)                             // gibt eine Inzidenzmatrix einer
{                                                // Komponente aus
   int i, j, dim1 = 0, dim2 = 0, d1, d2;

   for (i=1; i <= nr[dim];   i++)  if (zusammen[dim][i]   == komp)  dim1++;
   if (dim1 == 0) return;
   for (j=1; j <= nr[dim-1]; j++)  if (zusammen[dim-1][j] == komp)  dim2++;
   if (dim2 == 0) return;

   for (i=1; i <= nr[dim]; i++)
      for (j=1; j <= nr[dim-1]; j++)  inzidenz[i][j] = 0;
   for (i=1; i <= nr[dim]; i++)
      if (zusammen[dim][i] == komp)  Rand(i);
   printf("\n{\n");
   d2 = 0;
   for (j=1; j <= nr[dim-1]; j++)
      if (zusammen[dim-1][j] == komp)
      {
         d2++; d1 = 0;
         printf("{");
         for (i=1; i <= nr[dim]; i++)
            if (zusammen[dim][i] == komp)
            {
               d1++;
               printf("%d", inzidenz[i][j]);
               if (d1 == dim1)
               {
                  printf("}");
                  if (d2 == dim2)  printf("}");
               }
               printf(",");
            }
         printf("\n");
      }
}

void Liste_ausgeben(void)                        // gibt alle Fixpunkt-Zellen nach
{                                                // Komponenten sortiert aus
   int i, j, k, komp, num;

   for (komp=1; komp <= komponente; komp++)
   {
      printf("\nKomponente %d\n\n", komp);
      for (dim=1; dim <= Dimension; dim++)
      {
         num = 0;
         for (i=1; i <= nr[dim]; i++)
            if (zusammen[dim][i] == komp)
            {
               num++;
               printf("d=%d, Nr.%2d (%d): ", dim, num, i);
               for (j=1; j<=anz_m[dim][i]; j++)
               {
                  printf(" ");
                  for (k=1; k<=anz_n[dim][i]; k++)  printf("%d", liste[dim][i][j][k]);
               }
               printf("\n");
            }
      }
   }
}

void main (void)
{
   int i;

   clrscr();
   printf("h = %d, c = %d\n\n", h, c);
   Init();
   for (n=2; n <= 2*h; n++)                      // erzeuge alle Fixpunkte
   {
      m = 1;
      Neue_Schicht();
   }
   Komponenten();                                // bestimme die Komponenten
   Komp_Groessen();                              // und ihre Groessen
   Liste_ausgeben();
   for (i=1; i <= komponente; i++)               // gib die Inzidenzmatrizen aus
   {
      printf("\n\nKomponente %d\n", i);
      for (dim=2; dim <= Dimension; dim++)  Delta(i);
      printf("};\n");
   }
   getch();
}

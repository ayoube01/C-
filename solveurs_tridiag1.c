/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* %%%            Validation pour Solveur_Chol_TriDiag            %%% */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* %%%  Auteur: FL, dernière modification: 02/11/2014             %%% */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */

/* ================================================================== */
/* ===               INCLUSION DES LIBRAIRIES DU C                === */
/* ================================================================== */
#include <math.h>   /* sqrt()...         */
#include <stdio.h>  /* printf()...       */
#include <stdlib.h> /* malloc()...       */
#include <stddef.h> /* NULL...           */
#include <string.h> /* atof,memcpy...    */
#include <time.h>   /* clock()...        */
#include <limits.h> /* bornes limites... */
#include <assert.h> /* assert()...       */

/* ================================================================== */
/* ===                  DEFINITION DES MACROS                     === */
/* ================================================================== */
#define min(a,b) ( (a)<(b) ? (a) : (b) )
#define max(a,b) ( (a)>(b) ? (a) : (b) )
#define CARRE(x) ( (x)*(x) )
#define  CUBE(x) ( (x)*(x)*(x) )
#define   ABS(x) ( (x)>=0 ? (x) : -(x) )
#define  SIGN(x) ( (x)>0 ? (+1) : ( (x)<0 ? (-1) : (0) ) )

#define Inf  (1.0/0.0)
#define NaN  (0.0/0.0)
#define PI    3.141592653589793e+000

/* ================================================================== */
/* ===                  DECLARATION DES STRUCTURES                === */
/* ================================================================== */

typedef struct {
    int     length; /* nombre de coeff.  du vecteur */
    double *tab;    /* les coefficients  du vecteur */
} Vecteur;

typedef struct {
   int dim;
   char type;
   Vecteur *lA;
   Vecteur *dA;
   Vecteur *uA;
}Tridiag;






/* ================================================================== */
/* ===                  DECLARATION DES ENTETES                   === */
/* ================================================================== */

double *alloc_dtab_1d(int n);
double *free_dtab_1d(double *v);
Vecteur *alloc_Vecteur(int n);
Vecteur *free_Vecteur(Vecteur *V);
void   disp_Vecteur(char *message, Vecteur *V);
double dot(double *v, double *w, int d, int f);
double NORM    (Vecteur *x);
double NORM_INF(Vecteur *x);
void Solveur_LU_TriDiag(Vecteur *dA, Vecteur *lA, Vecteur *uA, Vecteur *b, Vecteur *x);
void Solveur_TriDiag_Chol(Vecteur *dA, Vecteur *lA, Vecteur *b, Vecteur *x);


/* ================================================================== */
/* ===                   DEFINITION DES FONCTIONS                 === */
/* ================================================================== */

/* ===================== */
/* allocation dynamique: */
/* ===================== */
/* ------------------------------------------------------------------------- */
double *alloc_dtab_1d(int n)
{ assert( n>0 );
  { double *v=NULL;

    v = (double *) calloc (n, sizeof(double));
    if (v == NULL) {
      printf("plus de place pour alloc_dtab_1d()\n"); exit(1);
    }
    v = v-1;

    return v;
  }
}

double *free_dtab_1d(double *v)
{ if (v!=NULL) free(v+1);
  return NULL;
}

Vecteur *alloc_Vecteur(int n)
{ assert( n>0 );
  { Vecteur *V=NULL;

    V = (Vecteur *) malloc(sizeof(Vecteur));
    if (V == NULL) {
      printf("plus de place pour alloc_Vecteur()\n"); exit(1);
    };
    V->length = n;
    V->tab    = alloc_dtab_1d(n);
    return V;
  }
}
Vecteur *alloc_Tridiag(int n ,char type)
{    Tridiag *A;
      A->dim=n;
      A->type=type;
    if (type=='G')
    {
        A->lA=alloc_Vecteur(n-1);
        A->dA=alloc_Vecteur(n);
        A->uA=alloc_Vecteur(n-1);
    }
      if (type=='L')
    {
        A->lA=alloc_Vecteur(n-1);
        A->dA=alloc_Vecteur(n);

    }
      if (type=='U')
    {
        A->dA=alloc_Vecteur(n);
        A->uA=alloc_Vecteur(n-1);
    }
      if (type=='S')
    {
        A->lA=alloc_Vecteur(n-1);
        A->dA=alloc_Vecteur(n);


    }
      return A;
}
Vecteur *free_Vecteur(Vecteur *V)
{ if (V!=NULL) {
    free_dtab_1d(V->tab);
    free(V);
  }
  return NULL;
}


/* ===================== */
/* ====     disp:   ==== */
/* ===================== */
/* ------------------------------------------------------------------------- */
void disp_Vecteur(char *message, Vecteur *V)
{ assert(V!=NULL);
  { int i;
    printf("%s\n",message);
    for (i=1; i<=V->length; i++) printf("%24.16e\n",V->tab[i]);
    printf("\n");
  }
}


/* ============ */
/* === BLAS === */
/* ============ */
/* -------------------------- BLAS de bas niveau --------------------------- */
double dot(double *v, double *w, int d, int f)
/* v(d:f)'*w(d:f) */
{ assert(v!=NULL); assert(w!=NULL);
  { int    i;
    double s=0.0;
    for (i=d; i<=f; i++) s += v[i]*w[i];
    return s;
  }
}

/* -------------------------- BLAS de haut niveau -------------------------- */
double NORM(Vecteur *x)
/* ||x||_2 = sqrt( (x|x) ) */
{ assert(x!=NULL);
  return sqrt(dot(x->tab, x->tab, 1, x->length));
}

double NORM_INF(Vecteur *x)
/* ||x||_inf = max_{1<=i<=n} |x_i| */
{ assert(x!=NULL);
  { int i;
    double nx=0.0;
    for (i=1; i<=x->length; i++) nx = max(nx, fabs(x->tab[i]));
    return nx;
  }
}




/* ====================== */
/* ====  Solveurs:   ==== */
/* ====================== */
/* ------------------------------------------------------------------------- */
void Solveur_LU_TriDiag(Vecteur *dA, Vecteur *lA, Vecteur *uA, Vecteur *b, Vecteur *x)
{     Vecteur *lL;
    Vecteur *dU,*uU;

    int i;
 lL=alloc_Vecteur(lA->length);
 dU=alloc_Vecteur(dA->length);
 uU=alloc_Vecteur(uA->length);
 /*initialisation de notre L U */
   dU->tab[1]=dA->tab[1];
   for (i=1;i<(dA->length);i++)
   {
       lL->tab[i]=lA->tab[i]/dU->tab[i];
       uU->tab[i]=uA->tab[i];
       dU->tab[i+1]=dA->tab[i+1]-lL->tab[i]*uU->tab[i];
   }
   /*descente*/
   Vecteur *Y;
   Y=alloc_Vecteur(dA->length);
   Y->tab[1]=b->tab[1];
   for(i=1;i<(Y->length);i++)
   {
       Y->tab[i+1]=b->tab[i+1]-(lL->tab[i]*Y->tab[i]);
   }
   /*remontee*/
   x->tab[x->length]=Y->tab[Y->length]/dU->tab[dU->length];

   for(i=(x->length)-1;i>1;i--)
   {
      x->tab[i]=(Y->tab[i]-(uA->tab[i]*x->tab[i+1]))/(dU->tab[i]);

   }

}


void Solveur_Chol_TriDiag(Vecteur *dA, Vecteur *lA, Vecteur *b, Vecteur *x)
{          /*decomposition de la matrice*/
            Vecteur *dL;
            Vecteur *lL;
            dL=alloc_Vecteur(dA->length);
            lL=alloc_Vecteur(lA->length);
            int i;
            dL->tab[1]=sqrt(dA->tab[1]);
            for(i=1;i<dA->length;i++)
            {
                lL->tab[i]=lA->tab[i]/dL->tab[i];
               dL->tab[i+1]=sqrt(dA->tab[i+1]-(lL->tab[i]*lL->tab[i]));

            }
            /*résolution de L*Y = b*/
           Vecteur *Y;
            Y=alloc_Vecteur(dA->length);
            Y->tab[1]=b->tab[1]/dL->tab[1];
            for(i=2;i<dA->length+1;i++)
            {
              Y->tab[i]=(b->tab[i]-lL->tab[i-1]*Y->tab[i-1])/dL->tab[i] ;

            }
            /*resolution de L(transposé)*x=Y*/
           x->tab[dA->length]=Y->tab[dA->length]/dL->tab[dA->length];
           for(i=dA->length-1;i>0;i--)
           {
              x->tab[i]=(Y->tab[i]-lL->tab[i]*x->tab[i+1])/dL->tab[i];
           }


}
void Solveur_LU_TriDiag1(Tridiag *A, Vecteur *b, Vecteur *x)
{   Vecteur *lL;
    Vecteur *dU,*uU;
if(A->type=='G'){
    int i;
 lL=alloc_Vecteur(A->dim-1);
 dU=alloc_Vecteur(A->dim);
 uU=alloc_Vecteur(A->dim-1);
 /*initialisation de notre L U */
   dU->tab[1]=A->dA->tab[1];
   for (i=1;i<(A->dim);i++)
   {
       lL->tab[i]=A->lA->tab[i]/dU->tab[i];
       uU->tab[i]=A->uA->tab[i];
       dU->tab[i+1]=A->dA->tab[i+1]-lL->tab[i]*uU->tab[i];
   }
   /*descente*/
   Vecteur *Y;
   Y=alloc_Vecteur(A->dim);
   Y->tab[1]=b->tab[1];
   for(i=1;i<(Y->length);i++)
   {
       Y->tab[i+1]=b->tab[i+1]-(lL->tab[i]*Y->tab[i]);
   }
   /*remontee*/
   x->tab[x->length]=Y->tab[Y->length]/dU->tab[dU->length];

   for(i=(x->length)-1;i>0;i--)
   {
      x->tab[i]=(Y->tab[i]-(A->uA->tab[i]*x->tab[i+1]))/(dU->tab[i]);

   }
}
if(A->type=='S'){
    int i;
 lL=alloc_Vecteur(A->dim-1);
 dU=alloc_Vecteur(A->dim);
 uU=alloc_Vecteur(A->dim-1);
 /*initialisation de notre L U */
   dU->tab[1]=A->dA->tab[1];
   for (i=1;i<(A->dim);i++)
   {
       lL->tab[i]=A->lA->tab[i]/dU->tab[i];
       uU->tab[i]=A->lA->tab[i];
       dU->tab[i+1]=A->dA->tab[i+1]-lL->tab[i]*uU->tab[i];
   }
   /*descente*/
   Vecteur *Y;
   Y=alloc_Vecteur(A->dim);
   Y->tab[1]=b->tab[1];
   for(i=1;i<(Y->length);i++)
   {
       Y->tab[i+1]=b->tab[i+1]-(lL->tab[i]*Y->tab[i]);
   }
   /*remontee*/
   x->tab[x->length]=Y->tab[Y->length]/dU->tab[dU->length];

   for(i=(x->length)-1;i>0;i--)
   {
      x->tab[i]=(Y->tab[i]-(A->lA->tab[i]*x->tab[i+1]))/(dU->tab[i]);

   }
}
if (A->type=='L')
{   int i;
    Vecteur *Y;
   Y=alloc_Vecteur(A->dim);
     x->tab[1]=b->tab[1]/A->dA->tab[1];
    for(i=1;i<x->length+1;i++)
   {
      x->tab[i]=(b->tab[i]-(A->lA->tab[i-1]*x->tab[i-1]))/(A->dA->tab[i]);

   }

}
if (A->type=='U')
{   int i;
    Vecteur *Y;
   Y=alloc_Vecteur(A->dim);
     x->tab[x->length]=b->tab[Y->length]/A->dA->tab[A->dim];


     for(i=x->length-1;i>1;i--)
   {
      x->tab[i]=(b->tab[i]-(A->uA->tab[i]*x->tab[i+1]))/(A->dA->tab[i]);

   }

}

}
void Solveur_Chol_TriDiag1(Tridiag *A, Vecteur *b, Vecteur *x)
{
    if(A->type=='G' || A->type=='S' )/*n'est applicable qu'aux matrice symétrique*/
  {
    /*decomposition de la matrice*/
            Vecteur *dL;
            Vecteur *lL;
            dL=alloc_Vecteur(A->dA->length);
            lL=alloc_Vecteur(A->lA->length);
            int i;
            dL->tab[1]=sqrt(A->dA->tab[1]);
            for(i=1;i<A->dA->length;i++)
            {
                lL->tab[i]=A->lA->tab[i]/dL->tab[i];
               dL->tab[i+1]=sqrt(A->dA->tab[i+1]-(lL->tab[i]*lL->tab[i]));

            }
            /*résolution de L*Y = b*/
           Vecteur *Y;
            Y=alloc_Vecteur(A->dA->length);
            Y->tab[1]=b->tab[1]/dL->tab[1];
            for(i=2;i<A->dA->length+1;i++)
            {
              Y->tab[i]=(b->tab[i]-lL->tab[i-1]*Y->tab[i-1])/dL->tab[i] ;

            }
            /*resolution de L(transposé)*x=Y*/
           x->tab[A->dA->length]=Y->tab[A->dA->length]/dL->tab[A->dA->length];
           for(i=A->dA->length-1;i>0;i--)
           {
              x->tab[i]=(Y->tab[i]-lL->tab[i]*x->tab[i+1])/dL->tab[i];
           }
}



}
//fonction qui calcule le produit matricielle
void produit_tridiag_vecteur(Tridiag *A , Vecteur *X,Vecteur *Y)
{   int i;
    if(A->type=='G')
    {
        Y->tab[1]=A->dA->tab[1]*X->tab[1]+A->uA->tab[1]*X->tab[2];
        for (i=2;i<Y->length;i++)
        {
            Y->tab[i]=A->lA->tab[i-1]*X->tab[i-1]+A->dA->tab[i]*X->tab[i]+A->uA->tab[i]*X->tab[i+1];
        }
        Y->tab[Y->length]=A->lA->tab[Y->length-1]*X->tab[Y->length-1]+A->dA->tab[Y->length]*X->tab[Y->length];
    }
    if(A->type=='S')
    {
        Y->tab[1]=A->dA->tab[1]*X->tab[1]+A->lA->tab[1]*X->tab[2];
        for (i=2;i<Y->length;i++)
        {
            Y->tab[i]=A->lA->tab[i-1]*X->tab[i-1]+A->dA->tab[i]*X->tab[i]+A->lA->tab[i]*X->tab[i+1];
        }
        Y->tab[Y->length]=A->lA->tab[Y->length-1]*X->tab[Y->length-1]+A->dA->tab[Y->length]*X->tab[Y->length];
    }
 if(A->type=='U')
    {
        Y->tab[1]=A->dA->tab[1]*X->tab[1]+A->uA->tab[1]*X->tab[2];
        for (i=2;i<Y->length;i++)
        {
            Y->tab[i]=A->dA->tab[i]*X->tab[i]+A->uA->tab[i]*X->tab[i+1];
        }
        Y->tab[Y->length]=A->dA->tab[Y->length]*X->tab[Y->length];
    }
     if(A->type=='L')
    {
        Y->tab[1]=A->dA->tab[1]*X->tab[1];
        for (i=2;i<Y->length;i++)
        {
            Y->tab[i]=A->lA->tab[i-1]*X->tab[i-1]+A->dA->tab[i]*X->tab[i];
        }
        Y->tab[Y->length]=A->lA->tab[Y->length-1]*X->tab[Y->length-1]+A->dA->tab[Y->length]*X->tab[Y->length];
    }
}
double residu(Vecteur *b,Tridiag *A,Vecteur *X)
{   double resultat;
int i;
    Vecteur *Y;
    Y=alloc_Vecteur(X->length);
    produit_tridiag_vecteur(A,X,Y);
    for (i=1;i<b->length+1;i++)
    {
        b->tab[i]=b->tab[i]-Y->tab[i];
    }
    resultat=NORM_INF(b);
    return resultat;
}
/* ============================================= */
/* ====   Test pour Solveur_Chol_TriDiag:   ==== */
/* ============================================= */
/* --------------------------------------------------------------- */
int main(int narg, char **arg)
{/*saisie de notre matrice A et b et x*/


/*
Vecteur* dA;
     dA= alloc_Vecteur(5);
     dA->tab[1]=2;
     dA->tab[2]=4;
     dA->tab[3]=4;
     dA->tab[4]=4;
     dA->tab[5]=2;
Vecteur* lA;
     lA= alloc_Vecteur(4);
     lA->tab[1]=-1;
     lA->tab[2]=-1;
     lA->tab[3]=-1;
     lA->tab[4]=-1;



Vecteur* uA;
     uA= alloc_Vecteur(4);
     uA->tab[1]=-1;
     uA->tab[2]=-1;
     uA->tab[3]=-1;
     uA->tab[4]=-1;



Vecteur* x;
     x= alloc_Vecteur(5);
     x->tab[1]=9999;
     x->tab[2]=4;
     x->tab[3]=14;
     x->tab[4]=9999;
     x->tab[5]=9999;

Vecteur* b;
     b= alloc_Vecteur(5);
     b->tab[1]=0;
     b->tab[2]=4;
     b->tab[3]=6;
     b->tab[4]=8;
     b->tab[5]=6;

    Solveur_LU_TriDiag(dA, lA, uA, b, x);
    disp_Vecteur("\n",x);


    Solveur_Chol_TriDiag(dA,lA,b,x);

    disp_Vecteur("\n",x);
*/

/* exemple avec A general*/


/*
Tridiag *A;
A=alloc_Tridiag(5,'G');
(A->dA)->tab[1]=2;
(A->dA)->tab[2]=4;
(A->dA)->tab[3]=4;
(A->dA)->tab[4]=4;
(A->dA)->tab[5]=2;

(A->lA)->tab[1]=-1;
(A->lA)->tab[2]=-1;
(A->lA)->tab[3]=-1;
(A->lA)->tab[4]=-1;

(A->uA)->tab[1]=-1;
(A->uA)->tab[2]=-1;
(A->uA)->tab[3]=-1;
(A->uA)->tab[4]=-1;
Vecteur* x;
     x= alloc_Vecteur(5);
     x->tab[1]=9;
     x->tab[2]=78;
     x->tab[3]=18;
     x->tab[4]=8;
     x->tab[5]=6;
Vecteur* b;
     b= alloc_Vecteur(5);
     b->tab[1]=0;
     b->tab[2]=4;
     b->tab[3]=6;
     b->tab[4]=8;
     b->tab[5]=6;
Solveur_LU_TriDiag1(A,b,x);
disp_Vecteur("\n",x);
Solveur_Chol_TriDiag1(A,b,x);

printf("%d",A->dim);
disp_Vecteur("\n",x);
*/

/*
Tridiag *A;
A=alloc_Tridiag(3,'L');
(A->dA)->tab[1]=2;
(A->dA)->tab[2]=2;
(A->dA)->tab[3]=2;



(A->lA)->tab[1]=1;
(A->lA)->tab[2]=1;


Vecteur* x;
     x= alloc_Vecteur(3);
     x->tab[1]=7;
     x->tab[2]=7;
     x->tab[3]=8;

Vecteur* b;
     b= alloc_Vecteur(3);
     b->tab[1]=2;
     b->tab[2]=3;
     b->tab[3]=3;

Solveur_LU_TriDiag1(A,b,x);

disp_Vecteur("\n",x);
*/


/*
Tridiag *A;
A=alloc_Tridiag(3,'U');
(A->dA)->tab[1]=2;
(A->dA)->tab[2]=2;
(A->dA)->tab[3]=2;



(A->uA)->tab[1]=1;
(A->uA)->tab[2]=1;


Vecteur* x;
     x= alloc_Vecteur(3);
     x->tab[1]=5;
     x->tab[2]=2;
     x->tab[3]=5;

Vecteur* b;
     b= alloc_Vecteur(3);
     b->tab[1]=3;
     b->tab[2]=3;
     b->tab[3]=2;

Solveur_LU_TriDiag1(A,b,x);

disp_Vecteur("\n",x);
*/

/*
Tridiag *A;
A=alloc_Tridiag(5,'S');
(A->dA)->tab[1]=2;
(A->dA)->tab[2]=4;
(A->dA)->tab[3]=4;
(A->dA)->tab[4]=4;
(A->dA)->tab[5]=2;

(A->lA)->tab[1]=-1;
(A->lA)->tab[2]=-1;
(A->lA)->tab[3]=-1;
(A->lA)->tab[4]=-1;


Vecteur* x;
     x= alloc_Vecteur(5);
     x->tab[1]=9;
     x->tab[2]=78;
     x->tab[3]=18;
     x->tab[4]=8;
     x->tab[5]=6;
Vecteur* b;
     b= alloc_Vecteur(5);
     b->tab[1]=0;
     b->tab[2]=4;
     b->tab[3]=6;
     b->tab[4]=8;
     b->tab[5]=6;
 Solveur_LU_TriDiag1( A, b, x);
 disp_Vecteur("\n",x);
 Solveur_Chol_TriDiag1( A, b, x);
disp_Vecteur("\n",x);
*/


/*

Tridiag *A;
A=alloc_Tridiag(5,'G');
(A->dA)->tab[1]=2;
(A->dA)->tab[2]=4;
(A->dA)->tab[3]=4;
(A->dA)->tab[4]=4;
(A->dA)->tab[5]=2;

(A->lA)->tab[1]=-1;
(A->lA)->tab[2]=-1;
(A->lA)->tab[3]=-1;
(A->lA)->tab[4]=-1;

(A->uA)->tab[1]=-1;
(A->uA)->tab[2]=-1;
(A->uA)->tab[3]=-1;
(A->uA)->tab[4]=-1;
Vecteur* x;
     x= alloc_Vecteur(5);
     x->tab[1]=9;
     x->tab[2]=78;
     x->tab[3]=18;
     x->tab[4]=8;
     x->tab[5]=6;
Vecteur* b;
     b= alloc_Vecteur(5);
     b->tab[1]=0;
     b->tab[2]=4;
     b->tab[3]=6;
     b->tab[4]=8;
     b->tab[5]=6;
Solveur_LU_TriDiag1(A,b,x);
  disp_Vecteur("\n",x);
Vecteur* Y;
Y= alloc_Vecteur(5);
produit_tridiag_vecteur(A ,x,Y);
    disp_Vecteur("\n",Y);

printf("le residu est de : %f",residu(b,A,x));


*/

    return EXIT_SUCCESS;
}

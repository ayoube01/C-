/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* %%%  Résolution du PB Elliptique 1D par Différences Finies:    %%% */
/* %%%      -(K(x)u')' + p(x) u = f(x)  dans [0,L]                %%% */
/* %%%                     u(0) = a    (Dirichlet)                %%% */
/* %%%                     u(L) = b    (Dirichlet)                %%% */
/* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
/* %%%  F.Lefèvre, dernière modification: 09/11/2014              %%% */
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
#define PI      3.141592653589793e+000
#define CAS_PB  1

/* ================================================================== */
/* ===                  DECLARATION DES STRUCTURES                === */
/* ================================================================== */

typedef struct {
    int     length; /* nombre de coeff.  du vecteur */
    double *tab;    /* les coefficients  du vecteur */
} Vecteur;

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

double fct_kappa(double x);
double fct_f    (double x);
double fct_p    (double x);
double fct_u    (double x);


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
  { return sqrt(dot(x->tab, x->tab, 1, x->length));
  }
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
/* ====   Solveur:   ==== */
/* ====================== */
/* ------------------------------------------------------------------------- */
void Solveur_LU_TriDiag(Vecteur *dA, Vecteur *lA, Vecteur *uA, Vecteur *b, Vecteur *x)
{
                   /* ============================ */
                   /*  ...   A  COMPLETER   ....   */
                   /* ============================ */
}

/* ------------------------------------------------------------------------- */
void Solveur_TriDiag_Chol(Vecteur *dA, Vecteur *lA, Vecteur *b, Vecteur *x)
{
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


/* ============================ */
/* ====   Pb Elliptique:   ==== */
/* ============================ */
/* ------------------------------------------------------------------------- */
double fct_u(double x)
{ double z;
/////////////cas 1
  z = sin(x);
/////////////cas 2
 /* double e;
  e = fct_kappa(0.0);
  z = exp(2/sqrt(e))/(exp(2/sqrt(e))-1)*exp(-x/sqrt(e)) + 1/(1-exp(2/sqrt(e)))*exp(x/sqrt(e));
  */
///////////cas 3
 // z = CARRE(x)*(1-x);


  return z;
}

double fct_kappa(double x)
{ double z;

////////////////cas 1
  z = 1;
/////////cas 2

  //z = 0.00001;
////////////cas 3
 // z = 1;

  return z;

}

double fct_f(double x)
{ double z;

///////////////cas 1
  z = 2.0*sin(x);


/////////////cas 2
  //z = 0;
//////////cas 3
  //z = 2*x-2;


  return z;
}

double fct_p(double x)
{ double z;
/////////////cas 1
  z = 1;
  //////////////cas 2
  //z = 1;
/////////////////////cas 3
 // z = 0;


  return z;
}



/* --------------------------------------------------------------- */
int main(int narg, char **arg)
{   int i;
    int N;
    double h;
    double L=1;
    double alpha =0;
    double beta=-1;
    printf("tapez s'il vous plait le nombre de subdivision N :");
    scanf("%d",&N);
    h=L/N;

///////////vecteur x de pas
    Vecteur *x;
    x=alloc_Vecteur(N+1);
    x->tab[1]=0;
    for(i=2;i<N;i++)
    {
        x->tab[i]=x->tab[0]+i*h;
    }
    x->tab[N+1]=L;
    /////////////vecteur p
    Vecteur *p;
    p=alloc_Vecteur(N-1);
    for (i=1;i<N;i++)
    {
       p->tab[i]=fct_p(x->tab[i]);
    }

 disp_Vecteur("\n",x);

///////vecteur u
Vecteur *u;
u= alloc_Vecteur(N-1);

//////////notre vecteur F
Vecteur *b;
b=alloc_Vecteur(N-1);
b->tab[1]=fct_f(x->tab[1])+fct_kappa(x->tab[1]+h/2)*alpha/(h*h);
for (i=2;i<N-1;i++)
{
    b->tab[i]=fct_f(x->tab[i]);
}
b->tab[N-1]=fct_f(x->tab[N-1])+fct_kappa(L-h/2)*alpha/(h*h);
disp_Vecteur("\n",u);
///////////////////matrice K
Vecteur *dk;
dk= alloc_Vecteur(N-1);
for (i=1;i<N;i++)
{
    dk->tab[i]=(1/(h*h))*fct_kappa(x->tab[i]+h/2)+fct_kappa(x->tab[i+1]+h/2);
}
disp_Vecteur("\n",dk);

Vecteur *lk;
lk=alloc_Vecteur(N-2);
for (i=1;i<N-1;i++)
{
    lk->tab[i]=-1*(1/(h*h))*fct_kappa(x->tab[i+1]+h/2);
}
Vecteur *dA;
dA=alloc_Vecteur(N-1);
for (i=1;i<N;i++)
{
    dA->tab[i]=p->tab[i]+dk->tab[i];
}
disp_Vecteur("\n",lk);
Solveur_TriDiag_Chol(dA,lk, b, u);
disp_Vecteur("\n",u);
  return EXIT_SUCCESS;


}

#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#define pi=3.14
double get_last(FILE *f_in)
{
    double c;
  c=0;
  if ((f_in = fopen(f_in,"r")) == NULL)
    {
      fprintf(stderr, "\nErreur: Impossible de lire le fichier %s\n",f_in);
      return(EXIT_FAILURE);
    }
     fscanf(f_in, "%lf", &c);

    while(fgetc(f_in)!=EOF)

          fscanf(f_in, "%lf", &c);

 fclose(f_in);
return c;
  }
int compteur(FILE *f_in)
{
    int c;
   int i=0;
  if ((f_in = fopen(f_in,"r")) == NULL)
    {
      fprintf(stderr, "\nErreur: Impossible de lire le fichier %s\n",f_in);
      return(EXIT_FAILURE);
    }

  while ((c = fgetc(f_in)) != EOF)
      {
          if (c==' ')
          i=i+1;

  }

  fclose(f_in);
  return i;
  }

double f1(double x)
{
    return 4*sqrt(1-x*x);
}
double A_milieu(double(*f)(double),double a,double b,int n)
{   double x;
    double s;
    double h;
    h=(b-a)/n;
    s=0.0;
    int i ;
    for (i=0;i<n;i++)
    {
        x=a+i*h;
        s+=h*f(x+h/2.0);

    }
    return s;
}
double int_montecarlo_save(double(*f)(double),double a, double b,int n)
{    int i;
    double res=0.0;

    int c;
    FILE *fichier=NULL;
    fichier = fopen("resultat.txt", "a");
    res=get_last("resultat.txt");
    for ( i=compteur("resultat.txt");i<n;i++)

    {

        res+=f(a+rand()*1.0/RAND_MAX*(b-a));
         fprintf(fichier,"%f",res);
         fputc(' ',fichier);
    }


return res/n;
}
double int_montecarlo(double(*f)(double),double a, double b,int n)
{    int i;
    double res=0.0;

    int c;

    for ( i=0;i<n;i++)

    {

        res+=f(a+rand()*1.0/RAND_MAX*(b-a));


    }


return res/n;
}
double simpson(double(*f)(double),double a, double b, int n)
{   double res;
    res=0.0;
    int i;
    double h;
    h=(b-a)/n;
    res=f(a)+f(b)+4*f(a+h/2);
    for (i=0;i<n;i++)
    {
        res+=2*f(a+i*h)+4*f(a+(i+.5)*h);
    }
    return res*h/6;
}
double trapez(double(*f)(double),double a, double b, int n)
{   double res;
    res=0.0;
    int i;
    double h;
    h=(b-a)/n;
    res=(f(a)+f(b))/2;
    for (i=0;i<n;i++)
    {
        res+=f(a+i*h);
    }
    return res*h;
}
int main()
{ int a = 0;
double pi;
double s;
double b = 1;
int n = 1000;
//printf("par la methode du point du milieu :    %f \n",A_milieu(f1,a,b,n));

//printf("par la methode de montecarlo :     %f \n",int_montecarlo_save(f1,a,b,n));
//printf("par la methode  de simpson :    %f \n",simpson(f1,a,b,n));
//printf("par la methode de trapez :     %f \n",trapez(f1,a,b,n));
//printf("%f",get_last("resultat.txt"));
//printf("par la methode de montecarlo :     %f \n",int_montecarlo(f1,a,b,n));
    return 0;
}


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void histogram(int N, double x[N], double xmin, double xmax, int nbin){
    int i,bin;
    FILE *f;
    f=fopen("brownian.txt","w");
    double hist[nbin];
    double dx = (xmax-xmin)/nbin;
    // initialisation de l'histogramme a 0
    for(bin=0;bin<nbin;bin++){hist[bin]=0.0;}
    // on ajoute chacune des N valeur de x dans le bonn nbin de l'histogramme
    for(i=0;i<N;i++){
        if (x[i]<xmax && x[i]>xmin){
            bin = (int) ((x[i]-xmin)/dx) ;
            hist[bin] += 1.0/N/dx ;
        }
    }
    // on ecrit les valeurs de l'histogram dans un fichier
    for(bin=0;bin<nbin;bin++){
        fprintf(f,"%13.6e %13.6e \n",xmin+dx*(bin+0.5),hist[bin]);
    }
    fclose(f);
}

double uniform(){
    double num = (double)rand()/(double)RAND_MAX;
    return num;
}

/* generate a random value weighted within the normal (gaussian) distribution */
double gauss(){
    double u = uniform();
    double v = uniform();
    double x = sqrtf(-2 * log(u)) * cos(2 * M_PI * v);
    return x;
}

/* convert non-standard normal distribution to standard */
double convert(double m, double sig){
    return gauss()*sqrtf(sig)+m;
}

double* brownian1D(int N, double m, double sig, double tmax){
    FILE *f;
    f=fopen("brownian1D.txt","w");
    double epsilon = tmax/(double)(N-1);
    double dis=0;
    double tab[N];
    tab[0] = 0;
    double X[N];
    X[0]=0;
    fprintf(f,"%f %f\n",tab[0],X[0]);
    for(int i=1; i<N; i++){
        tab[i] = tab[0] + epsilon * i;
        X[i] = X[i-1] + convert(0,tab[i]-tab[i-1]);
        fprintf(f,"%f %f\n",tab[i],X[i]);
    }
    fclose(f);
}

double* brownian2D(int N, double m, double sig, double tmax){
    FILE *f;
    f=fopen("brownian2D.txt","w");
    double epsilon = tmax/(double)(N-1);
    double dis=0;
    double tab[N];
    tab[0] = 0;
    double X[N], Y[N];
    X[0]=0;
    Y[0]=0;
    fprintf(f,"%f %f\n",X[0],Y[0]);
    for(int i=1; i<N; i++){
        tab[i] = tab[0] + epsilon * i;
        X[i] = X[i-1] + convert(0,tab[i]-tab[i-1]);
        Y[i] = Y[i-1] + convert(0,tab[i]-tab[i-1]);
        fprintf(f,"%f %f\n",X[i],Y[i]);
    }
    fclose(f);
}

double* brownian3D(int N, double m, double sig, double tmax){
    FILE *f;
    f=fopen("brownian3D.txt","w");
    double epsilon = tmax/(double)(N-1);
    double dis=0;
    double tab[N];
    tab[0] = 0;
    double X[N], Y[N], Z[N];
    X[0]=0;
    Y[0]=0;
    Z[0]=0;
    fprintf(f,"%f %f %f\n",X[0],Y[0],Z[0]);
    for(int i=1; i<N; i++){
        tab[i] = tab[0] + epsilon * i;
        X[i] = X[i-1] + convert(0,tab[i]-tab[i-1]);
        Y[i] = Y[i-1] + convert(0,tab[i]-tab[i-1]);
        Z[i] = Z[i-1] + convert(0,tab[i]-tab[i-1]);
        fprintf(f,"%f %f %F\n",X[i],Y[i],Z[i]);
    }
    fclose(f);
}

int main(void) {
    int N = pow(10,3);
    srand(time(NULL));
    //tableaux tab;
    double tab[N];
    double Y[N];
    /*for(int i=0; i<N;i++){
        tab[i] = brownian1D(Y,0,10,N);
    }*/
    brownian3D(N,0,1,100);
    //histogram(N,Y,50,200,200);
}
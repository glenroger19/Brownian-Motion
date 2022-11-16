#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<time.h>

double uniform(double a, double b){
    double num = (double)rand()/(double)RAND_MAX;
    return a+(b-a)*num;
}

/* generate a random value weighted within the normal (gaussian) distribution */
double gauss(double m, double sig){
    //double u = (double)random() / RAND_MAX;
    //double v = (double)random() / RAND_MAX;
    double u = uniform(m,sig);
    double v = uniform(m,sig);
    double z = sqrt(-2 * log(u)) * cos(2 * M_PI * v);
    return z;
}

void histogram(int N, double x[N], double xmin, double xmax, int nbin){
    int i,bin;
    FILE *f;
    f=fopen("normal.txt","w");
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

int main(void) {
    int N = pow(10,6);
    srand(time(NULL));
    double tab[N];
    for(int i=0; i<N;i++){
        tab[i] = gauss(0,1);
        //printf("%f\n",tab[i]);
    }
    histogram(N,tab,-5,5,100);
}
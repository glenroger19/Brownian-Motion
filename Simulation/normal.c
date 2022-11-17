#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double uniform(){
    //double a = m-3*sig;
    //double b = m+3*sig;
    //double num = a+((b-a)*(double)rand())/(double)RAND_MAX;
    double num = (double)rand()/(double)RAND_MAX;
    printf("%f\n",num);
    return num;
}

double convert(double m, double sig){
    //printf("%f\n",uniform()*sqrtf(sig)+m);
    //printf("%f\n",fabs((uniform()-m)/sqrtf(sig)));
    return uniform()*sqrt(sig)+m;
    //return (uniform()-m)/sqrtf(sig);
}

/* generate a random value weighted within the normal (gaussian) distribution */
double gauss(double m, double sig){
    double u = convert(m,sig);
    double v = convert(m,sig);
    //printf("%f %f\n",u,v);
    double z = -2 * log(u) * cos(2 * M_PI * v);
    //printf("%f\n",-2 * log(u));
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
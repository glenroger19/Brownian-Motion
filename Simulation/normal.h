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

void brownian2D(int N, double m, double sig, double tmax){
    FILE*f1,*f2,*f3,*f4,*f5;
    f1=fopen("fich1.txt","a");
    f2=fopen("fich2.txt","a");
    f3=fopen("fich3.txt","a");
    f4=fopen("fich4.txt","a");
    f5=fopen("fich5.txt","a");
    double epsilon = tmax/(double)(N-1);
    double dis=0;
    double tab[N];
    tab[0] = 0;
    double X[N], Y[N];
    //fprintf(f1,"%f %f\n",X[0],Y[0]);
    for(int i=0; i<N; i++){
        X[0]=0;
        Y[0]=0;
        for(int j=0; j<N; j++){
            tab[j] = tab[0] + (epsilon * j);
            X[j] = X[j-1] + convert(0,tab[j]-tab[j-1]);
            Y[j] = Y[j-1] + convert(0,tab[j]-tab[j-1]);
            if(j==0){
                fprintf(f1,"%f %f\n",X[j],Y[j]);
            }
            if(j==9){
                fprintf(f2,"%f %f\n",X[j],Y[j]);
            }
            if(j==99){
                fprintf(f3,"%f %f\n",X[j],Y[j]);
            }
        }
    }
    fclose(f1);
    fclose(f2);
    fclose(f3);
}

double* brownian2D_plus(int N, double m, double sig, double tmax){
    FILE *f;
    f=fopen("brownian2D_plus.txt","w");
    double epsilon = tmax/(double)(N-1);
    double dis=0;
    double tab[N];
    tab[0] = 0;
    double X[N],Y[N],X1[N],Y1[N],X2[N],Y2[N],X3[N],Y3[N];
    X[0]=0;
    Y[0]=0;
    X1[0]=0;
    Y1[0]=0;
    X2[0]=0;
    Y2[0]=0;
    X3[0]=0;
    Y3[0]=0;
    fprintf(f,"%f %f %f %f %f %f %f %f\n",X[0],Y[0],X1[0],Y1[0],X2[0],Y2[0],X3[0],Y3[0]);
    for(int i=1; i<N; i++){
        tab[i] = tab[0] + epsilon * i;
        X[i] = X[i-1] + convert(0,tab[i]-tab[i-1]);
        Y[i] = Y[i-1] + convert(0,tab[i]-tab[i-1]);
        X1[i] = X1[i-1] + convert(0,tab[i]-tab[i-1]);
        Y1[i] = Y1[i-1] + convert(0,tab[i]-tab[i-1]);
        X2[i] = X2[i-1] + convert(0,tab[i]-tab[i-1]);
        Y2[i] = Y2[i-1] + convert(0,tab[i]-tab[i-1]);
        X3[i] = X3[i-1] + convert(0,tab[i]-tab[i-1]);
        Y3[i] = Y3[i-1] + convert(0,tab[i]-tab[i-1]);
        fprintf(f,"%f %f %f %f %f %f %f %f\n",X[i],Y[i],X1[i],Y1[i],X2[i],Y2[i],X3[i],Y3[i]);
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

/*int main(void) {
    int N = 9*pow(10,4);
    srand(time(NULL));
    //tableaux tab;
    double tab[N];
    double Y[N];
    brownian2D_plus(N,0,1,10000000000000);
    //histogram(N,Y,50,200,200);
}*/
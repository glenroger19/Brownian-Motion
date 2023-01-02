#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void histogram(int N, double x[N], double xmin, double xmax, int nbin){
    int i,bin;
    FILE *f;
    f=fopen("test.txt","w");
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
    double num = rand()/(double)RAND_MAX;
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

double* brownian1D(int N, double m, double tmax){                                   //mouvement brownien en 1D
    FILE *f;
    f=fopen("brownian1D.txt","w");
    double epsilon = tmax/(double)(N-1);
    double tab[N];
    tab[0] = 0;                                         //initialise le temps à t=0
    double X[N];
    X[0]=0;                                             //condition initiale car la particule part de la position x=0
    fprintf(f,"%f %f\n",tab[0],X[0]);
    for(int i=1; i<N; i++){
        tab[i] = tab[0] + epsilon * i;                  //incrémentation du temps 
        X[i] = X[i-1] + convert(0,tab[i]-tab[i-1]);     //la position suivante dépend de la dernière position
        fprintf(f,"%f %f\n",tab[i],X[i]);   
    }
    fclose(f);
}

void brownian2D(int N, double m, double tmax){                                      //mouvement brownien en 2D
    FILE*f1;
    f1=fopen("brownian2D.txt","w");
    double epsilon = tmax/(double)(N-1);
    double tab[N];
    tab[0] = 0;
    double X[N], Y[N], dist[N];
    X[0]=0;
    Y[0]=0;
    dist[0]=0;
    double dis_moy = 0;
    fprintf(f1,"%f %f %f\n",X[0],Y[0],dist[0]);
    for(int i=1; i<N; i++){
        tab[i] = tab[0] + epsilon * i;
        X[i] = X[i-1] + convert(0,tab[i]-tab[i-1]);
        Y[i] = Y[i-1] + convert(0,tab[i]-tab[i-1]);
        dist[i] = sqrt(pow(X[i]-X[i-1],2)+pow(Y[i]-Y[i-1],2));
        fprintf(f1,"%f %f %f\n",X[i],Y[i],dist[i]);
    }
    fclose(f1);
}

double* brownian3D(int N, double m, double tmax){                                   //mouvement brownien en 3D
    FILE *f;
    f=fopen("brownian3D.txt","w");
    double epsilon = tmax/(double)(N-1);
    double tab[N], X[N], Y[N], Z[N];
    tab[0] = 0;
    X[0]=0;
    Y[0]=0;
    Z[0]=0;
    fprintf(f,"%f %f %f\n",X[0],Y[0],Z[0]);
    for(int i=1; i<N; i++){
        tab[i] = tab[0] + epsilon * i;
        X[i] = X[i-1] + convert(0,tab[i]-tab[i-1]);
        Y[i] = Y[i-1] + convert(0,tab[i]-tab[i-1]);
        Z[i] = Z[i-1] + convert(0,tab[i]-tab[i-1]);
        fprintf(f,"%f %f %f\n",X[i],Y[i],Z[i]);
    }
    fclose(f);
}

double ink(int N, double m, double tmax, int fois){
    FILE* t0=fopen("t0.txt","w");                                                   //création d'un fichier à t=0
    FILE* t10=fopen("t10.txt","w");                                                 //création d'un fichier à t=10
    FILE* t100=fopen("t100.txt","w");                                               //création d'un fichier à t=100                                
    double epsilon = tmax/(double)(N-1);
    double tab[N];
    double X[N], Y[N], distance[N];
    double dist = 0;
    double q = 0;
    for(int j=0; j<fois; j++){
        X[0]=0;                                                                     //condition initiale sur x
        Y[0]=0;                                                                     //condition initiale sur y
        tab[0] = 0;
        distance[0]=0;                                                              //condition initiale sur la distance
        fprintf(t0,"%f %f\n",X[0],Y[0]);
        for(int i=1; i<N+1; i++){
            tab[i] = tab[0] + epsilon * i;
            X[i] = X[i-1] + convert(0,tab[i]-tab[i-1]);
            Y[i] = Y[i-1] + convert(0,tab[i]-tab[i-1]);
            distance[i] = sqrt(pow(X[i]-X[i-1],2)+pow(Y[i]-Y[i-1],2));              //calcule la distance de la particule i par rapport sa position précédente
            q += pow(distance[i]-distance[i-1],2);                                  //numérateur de la distance quadratique moyenne
        }
        fprintf(t10,"%f %f\n",X[10],Y[10]);
        fprintf(t100,"%f %f\n",X[100],Y[100]);
        dist += sqrt(pow(X[N]-X[0],2)+pow(Y[N]-Y[0],2));                            //somme de la distance de la dernière particule par rapport à sa position initiale
    }
    fclose(t0);
    fclose(t10);
    fclose(t100);
    printf("%f\n",dist/fois);                                                       //affiche la distance moyenne 
    printf("%f\n",q/fois);                                                          //affiche la distance quadratique moyenne
    return q/fois;
}

double echappement(int N, double rayon, int tmax, int fois){
    double epsilon = tmax/(double)(N-1);
    double tab[N], X[N], Y[N], Z[N];                                                //on se place en 3D avec x,y,z
    tab[0] = 0;
    X[0]=0;
    Y[0]=0;
    Z[0]=0;
    double p_ext = 0;
    double p_tot = 0;
    double distance = 0;
    for(int j=1; j<fois+1; j++){                                                    //on a ici 1000 particules
        for(int i=1; i<N+1; i++){                                                   
            tab[i] = tab[0] + epsilon * i;
            X[i] = X[i-1] + convert(0,tab[i]-tab[i-1]);
            Y[i] = Y[i-1] + convert(0,tab[i]-tab[i-1]);
            Z[i] = Z[i-1] + convert(0,tab[i]-tab[i-1]);
        }
        distance = sqrt(pow(X[N]-X[0],2)+pow(Y[N]-Y[0],2)+pow(Z[N]-Z[0],2));        //calcule de la distance après N-itérations
        if(distance>rayon){                                                         //regarde si la particule est à l'extérieure du sphère
            p_ext++;                                                                //compte le nombre de particule à l'extérieur
        }
    }
    return p_ext/(double)fois;
}

void plot_echap(int N, double rayon){                                               //calcule de la probabilité pour différents rayon
    FILE* fich = fopen("echappement.txt","w");
    double p_ext = 0;
    for(int i=1;i<rayon+1;i++){
        p_ext = echappement(N,i,100,1000);
        fprintf(fich,"%i %f\n",i,p_ext);
    }
    fclose(fich);  
}

int main(){
    int N = pow(10,3);
    srand(time(NULL));
    /*double tab[N];
    for(int i=0; i<N; i++){                                                         //affichage d'une courbe gaussienne
        tab[i] = convert(0,1);
    }
    histogram(N,tab,-5,5,100);                                                      //fin d'affichage
    brownian1D(N,0,100);
    brownian2D(N,0,100);
    brownian3D(N,0,100);*/
    ink(N,0,100,100);
    /*echappement(N,10,100,10000);
    plot_echap(N,100);*/
}
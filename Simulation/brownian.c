#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

int main(){
    srand(time(NULL));
    FILE *fich;
    fich = fopen("brownian.txt","w");
    int sigma = 2;
    int tmax = 10;
    int N = 2000;
    double epsilon = (double)tmax / N;
    double T[N];
    double X[N];
    X[0] = 0;
    for(int i=0; i<N; i++){
        T[i] = epsilon * i;
        X[i+1] = X[i] + sqrtf(sigma * epsilon) * (rand()/(double)RAND_MAX);
        fprintf(fich,"%f %f\n",T[i],X[i]);
    }
    fclose(fich);
}
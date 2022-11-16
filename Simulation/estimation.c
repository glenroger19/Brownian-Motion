#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<math.h>

double f(double x){
    double n = pow(x,3);
    return n;
}

double aleatoire(int nbMarches, int nbPas){
    FILE *fich;
    fich = fopen("estimation.txt","w");
    double y[nbMarches];
    for(int i = 0; i<nbMarches; i++){
        y[0] = 0;
        for(int i = 1; i<nbPas ; i++){
            double n = 1 + (((double)rand()/RAND_MAX) * (-2));
            printf("%f\n",n);
            y[i] = (y[i-1] + n);
            fprintf(fich,"%i %f\n",i,y[i]);
        }
    }
    fclose(fich);
}

double main(){
    srand(time(NULL));
    /*FILE *fich;
    fich = fopen("estimation.txt","w");
    int N;
    printf("Nombre de points ?\n");
    scanf("%d",&N);
    double somme=0;
    for(int i = 0;i<N;i++){
        somme += f(rand()%RAND_MAX);
        fprintf(fich,"%f\n",somme/(double)N);
    }
    fclose(fich);*/
    aleatoire(1000,100);
}
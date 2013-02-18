#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define X 4
#define Y 4
#define Z 4
#define T 4
#define nx 200
#define ny 200
#define nz 200
#define nt 800 //this is also the maximum iteration times
#define C 1 //velocity speed

//kernel code performing the iteration
void kernel_4(float*** __restrict__ Uc, float*** __restrict__ Up, float* __restrict__ c);

int main(){
    time_t start,end;
    time(&start);
    
    //building data array Up
    int i,j,k;
    float*** Up_tmp_1 = (float***)calloc(nx+3,sizeof(float**));
    if(Up_tmp_1==NULL){
        fprintf(stderr, "Could not allocate memory.\n");
        exit(1);
    }
    float** Up_tmp_2 = (float**)calloc((nx+3)*(ny+3), sizeof(float*));
    if(Up_tmp_2==NULL){
        fprintf(stderr, "Could not allocate memory.\n");
        exit(1);
    }
    float* Up_tmp_3 = (float*)calloc((nx+3)*(ny+3)*(nz+3),sizeof(float));
    if(Up_tmp_3==NULL){
        fprintf(stderr, "Could not allocate memory,\n");
        exit(1);
    }
    
    Up_tmp_1[0] = &Up_tmp_2[1];
    for(i=1;i<nx+3;i++)
        Up_tmp_1[i] = Up_tmp_1[i-1] + ny + 3;
    Up_tmp_2[0] = &Up_tmp_3[1];
    for(j=1,k=(nx+3)*(ny+3);j<k;j++)
        Up_tmp_2[j] = Up_tmp_2[j-1] + nz + 3;
    
    float*** Up = &Up_tmp_1[1];
    
    //building data array Uc
    float*** Uc_tmp_1 = (float***)calloc(nx+3,sizeof(float**));
    if(Uc_tmp_1==NULL){
        fprintf(stderr, "Could not allocate memory.\n");
        exit(1);
    }
    float** Uc_tmp_2 = (float**)calloc((nx+3)*(ny+3), sizeof(float*));
    if(Uc_tmp_2==NULL){
        fprintf(stderr, "Could not allocate memory.\n");
        exit(1);
    }
    float* Uc_tmp_3 = (float*)calloc((nx+3)*(ny+3)*(nz+3),sizeof(float));
    if(Uc_tmp_3==NULL){
        fprintf(stderr, "Could not allocate memory,\n");
        exit(1);
    }
    
    Uc_tmp_1[0] = &Uc_tmp_2[1];
    for(i=1;i<nx+3;i++)
        Uc_tmp_1[i] = Uc_tmp_1[i-1] + ny + 3;
    Uc_tmp_2[0] = &Uc_tmp_3[1];
    for(j=1,k=(nx+3)*(ny+3);j<k;j++)
        Uc_tmp_2[j] = Uc_tmp_2[j-1] + nz + 3;
    
    float*** Uc = &Uc_tmp_1[1];
    
    //initialize u0
    Uc[nx/2][ny/2][nz/2] = 1;
    
    //iteration process...
    //====================
    
    int n = 0; //iteration index
        
    float c[7];
    c[1] = C*C*T*T*nx*nx*4;
    c[1] /= X*X*nt*nt*3*2;
    c[2] = C*C*T*T*ny*ny*4;
    c[2] /= Y*Y*nt*nt*3*2;
    c[3] = C*C*T*T*nz*nz*4;
    c[3] /= Z*Z*nt*nt*3*2;
    c[4] = -c[1]/16;
    c[5] = -c[2]/16;
    c[6] = -c[3]/16;
    c[0] = -15*(c[1]+c[2]+c[3])/8+1;
    
    //calculate U(1,x,y,z)
    kernel_4(Uc, Up, c); 
    
    n++;
    
    //swipe and iterate
    for(i=0;i<7;i++)
        c[i] *= 2;
    
    float*** tmp;
    while(n+2<=nt){
        
        kernel_4(Up, Uc, c);
        
        tmp = Uc;
        Uc = Up;
        Up = tmp;
    
        n++;
    }
    
    time(&end);
    //Output results...
    //=================
    float dif = difftime(end, start);
    printf("The time duration for the initilization and the loop section is %.6f seconds.\n",dif);
    
    FILE* ofp;
    char filename[20];
    sprintf(filename, "results_extend_1.m");
    ofp = fopen(filename, "w");
    if(ofp==NULL){
        printf("Could not allocate memory.\n");
        exit(1);
    }
    
    fprintf(ofp,"U = zeros(%d,%d,%d);\n",nx+1,ny+1,nz+1);
    for(k=0;k<nz+1;k++){
        fprintf(ofp,"U(:,:,%d) = [\n", k+1);
        
        for(i=0;i<nx+1;i++){
            for(j=0;j<ny+1;j++)
                fprintf(ofp, "%f ", Up[i][j][k]);
            
            fprintf(ofp,";\n");
        }
        
        fprintf(ofp,"];\n");
    }
    fclose(ofp);
    
    
    //clean up...
    //===========
    free(Up_tmp_1);
    free(Up_tmp_2);
    free(Up_tmp_3);
    free(Uc_tmp_1);
    free(Uc_tmp_2);
    free(Uc_tmp_3);
    
}


void kernel_4(float*** __restrict__ Uc, float*** __restrict__ Up, float* __restrict__ c){
    int i,j,k;
    for(i = 1;i < nx;i++)
        for(j = 1;j < ny;j++)
            for(k = 1;k < nz; k++)
                Up[i][j][k] = -Up[i][j][k]
                + c[0]*Uc[i][j][k]
                + c[1]*(Uc[i+1][j][k]+Uc[i-1][j][k])
                + c[2]*(Uc[i][j+1][k]+Uc[i][j-1][k])
                + c[3]*(Uc[i][j][k+1]+Uc[i][j][k-1])
                + c[4]*(Uc[i+2][j][k]+Uc[i-2][j][k])
                + c[5]*(Uc[i][j+2][k]+Uc[i][j-2][k])
                + c[6]*(Uc[i][j][k+2]+Uc[i][j][k-2]);
    
    for(j = 1; j < ny; j++)
        for(k = 1; k < nz; k++){
            Up[-1][j][k] = -Up[1][j][k];
            Up[nx+1][j][k] = -Up[nx-1][j][k];
        }
    
    for(i = 1; i < nx; i++)
        for(k = 1; k < nz; k++){
            Up[i][-1][k] = -Up[i][1][k];
            Up[i][ny+1][k] = -Up[i][ny-1][k];
        }
    
    for(i = 1; i < nx; i++)
        for(j = 1; j < ny; j++){
            Up[i][j][-1] = -Up[i][j][1];
            Up[i][j][nz+1] = -Up[i][j][nz-1];
        }
}

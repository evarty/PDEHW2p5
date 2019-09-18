#include <stdio.h>
#include <math.h>
#include <stdlib.h>


void ThomasAlgo(double *r, double *unp1, double *a, double *b, double *c, int length);

int main(void){

    printf("print-1\n");

    //double P = 1;
    double deltax = 0.05;
    double deltat = 0.001;
    double xfinal = 1;
    double tfinal = 0.005;
    double ICvalue[(int)(xfinal/deltax)+1];
    //double X0BC = 0;
    double X1BC = 0;
    double mu = deltat/deltax/deltax;
    int length = (int)(xfinal/deltax);
    double r[length + 1];
    double g = 2;
    double theta = 0.5;
    double a[length + 1];
    double bdiag[length + 1];
    double c[length + 1];
    double unp1[length + 1];

    FILE *Output;
    Output = fopen("./Output.txt","w");
    fprintf(Output,"j n ujn\n");

    //Init zero array over the entire solution space
    //There are probably more efficient ways of handling the result, but this should work
    double ujn[length + 1][(int)(tfinal/deltat)+1];

    for(int j = 0; j < length + 1; j++){
        for(int n = 0; n < tfinal/deltat + 1; n++){
            ujn[j][n] = 0;
        }
    }
    printf("print0\n");

    //Set up the initial conditions
    for(int j = 0; j < length + 1; j++){
        if(j*deltax <= 0.5){
            ICvalue[j] = 2*(double)j*deltax;
        }else{
            ICvalue[j] = 2-2*(double)j*deltax;
        }
    }
    double IC[length + 1];
    for(int j = 0; j < length + 1; j++){
        IC[j] = ICvalue[j];
        ujn[j][0] = IC[j];
        //printf("%lf\n",ujn[j][0]);
    }
    printf("print1\n");

    //Init ujn with the boundary condition at x=1
    double BCx1[(int)(tfinal/deltat)];
    for(int n = 0; n < tfinal/deltat + 1; n++){
        BCx1[n] = X1BC;
        ujn[length][n] = BCx1[n];
        //printf("%lf\n",ujn[length][n]);
    }
    printf("print3\n");


    //Construct all ujn for the steps in the interior
    //Set up all the matrix values for thomas algo.
    //Solve via thomas algo.
    for(int n = 1; n < tfinal/deltat + 1; n++){
        for(int j = 0; j < length; j++){

            //Construct RHS of the equation.
            if(j==0){
                r[j] = -g*deltax*mu*theta;
                printf("r[0]=%lf\n",r[j]);
            }else if(j<length - 1){
                r[j] = ujn[j][n-1] + (1 - theta)*mu*(ujn[j-1][n-1] - 2*ujn[j][n-1] + ujn[j+1][n-1]);
                printf("r[%d]=%lf\n",j,r[j]);
            }else{
                r[j] = ujn[j][n-1] + (1 - theta)*mu*(ujn[j-1][n-1] - 2*ujn[j][n-1] + ujn[j+1][n-1]) + mu*theta*ujn[j+1][n];
                printf("r[%d]=%lf\n",j,r[j]);
            }

            //Construct a vector (subdiagonal)
            if(j==0){
                
                //All excepth the first and the last are the same.
            }else if(j<length){
                a[j] = -mu*theta;
            }

            //Construct b vector (diagonal)
            if(j==0){
                bdiag[j] = mu*theta;
            }else if(j<length){
                bdiag[j] = 1+2*mu*theta;
            }

            //Construct c vector (superdiagonal)
            if(j<length - 1){
                c[j] = -mu*theta;
            }
        }

        ThomasAlgo(r, unp1, a, bdiag, c, length);
        for(int j = 0; j < length + 1; j++){
            ujn[j][n] = unp1[j];
            //printf("%lf\n",unp1[j]);
        }

    }

    for(int n = 0; n < (int)(tfinal/deltat) + 1; n++){
        for(int j = 0; j < length + 1; j++){
            fprintf(Output,"%d %d %lf\n",j,n,ujn[j][n]);
        }
    }
    fclose(Output);

}

void ThomasAlgo(double *r, double *unp1, double *a, double *b, double *c, int length){

    double cprime[length];
    double rdoubleprime[length];
    double bprime[length];
    double rprime[length];
    double rtripleprime[length];

    cprime[0] = c[0]/b[0];
    rdoubleprime[0] = r[0]/b[0];

    for(int j = 1; j < length ; j++){

        bprime[j] = b[j] - a[j]*cprime[j-1]; 
        rprime[j] = r[j] - a[j]*rdoubleprime[j-1];
        cprime[j] = c[j]/bprime[j];
        rdoubleprime[j] = rprime[j]/bprime[j];
    }

    bprime[length] = b[length] - a[length]*cprime[length - 1];
    rprime[length] = r[length] - a[length]*rdoubleprime[length - 1];
    rtripleprime[length] = rprime[length]/bprime[length];

    for(int j = length - 1; j >= 0; j--){
        rtripleprime[j] = rdoubleprime[j] - cprime[j]*rtripleprime[j+1];
    }

    *unp1 = *rtripleprime;
    //for(int j = 0; j < length + 1;j++){
    //    printf("%lf\n",unp1[j]);
    //}



}

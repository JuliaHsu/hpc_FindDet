#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"
#include <time.h>
#include <math.h>
#include "mpfr.h"
#define Size 5000
#include "mpi.h"
#define MASTER 0 
#define From_Master 1
#define From_Worker 2
int main(int argc, char* argv[]){
    double START,END;
    double first, last;
    START=clock();
    int numTasks,taskId,numWorkers,source,dest,mtype,
        rows,averow,extra,offset;
    MPI_Status status;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskId); 
    MPI_Comm_size(MPI_COMM_WORLD, &numTasks);
    numWorkers = numTasks-1;
    double max;
    double tmp;
    double resDet;
    mpfr_t det;
    double mulTmp; 
    mpfr_t tmp1;
    mpfr_t tmp2;
    int i,j,k,maxRow,sign;
    char f_name[50];
    double logAbsDet;
    double* buffer;
    buffer = malloc(Size*Size*sizeof(double));
    mpfr_init(det);
    mpfr_init(tmp1);
    mpfr_init(tmp2);
    sprintf(f_name,"m5000x5000.bin");
    FILE *datafile=fopen(f_name,"rb");
    for(i=0;i<Size;i++){
        for(j=0;j<Size;j++){
            *(buffer+Size*i+j)=0; 
        }
    }
    if(taskId==MASTER){

        for(i=0;i<Size;i++){
            for(j=0;j<Size;j++){
                fread(buffer+(Size*i+j),sizeof(double),1,datafile);
                printf("%lf\t",*(buffer+(Size*i+j)));
            }
            printf("\n");
        }
        first=MPI_Wtime();
        averow = 5000/numWorkers;
        extra = 5000%numWorkers;
        offset = 0;
        mtype = From_Master;
        for (dest=1; dest<=numWorkers; dest++)
        {
            rows = (dest <= extra) ? averow+1 : averow;      
            printf("Sending %d rows to task %d offset=%d\n",rows,dest,offset);
            MPI_Send(&offset, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
            MPI_Send(&rows, 1, MPI_INT, dest, mtype, MPI_COMM_WORLD);
            MPI_Send(buffer+Size*offset, rows*5000, MPI_DOUBLE, dest, mtype,
                    MPI_COMM_WORLD);
            offset = offset + rows;
        }
        /* Receive results from worker tasks */
        mtype = From_Worker;
        for (i=1; i<=numWorkers; i++)
        {
            source = i;
            MPI_Recv(&offset, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
            MPI_Recv(&rows, 1, MPI_INT, source, mtype, MPI_COMM_WORLD, &status);
            MPI_Recv(&resDet,1, MPI_DOUBLE, source, mtype, 
                    MPI_COMM_WORLD, &status);
            printf("Received results from task %d\n",source);
            printf("det= %e",resDet);
        }
    }
    if(taskId > MASTER){
        mtype = From_Master;
        MPI_Recv(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
        MPI_Recv(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD, &status);
        MPI_Recv(buffer, rows*5000, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD, &status);
        for(i=0;i<Size-1;i++){
            maxRow = i;
            max = *(buffer+(Size*i+i));
            for(j=i+1;j<Size;j++){
                if(fabs(max)<fabs(*(buffer+Size*j+i))){
                    max = *(buffer+Size*j+i);
                    maxRow=j;
                }
            }
            if(fabs(max)<0.000001){
                resDet =0;
                printf("det=%e\n",resDet);
                exit(1);
            }
            if(maxRow!=i){
                sign++;
                for(j=i;j<Size;j++){
                    tmp = *(buffer + Size*i+j);
                    *(buffer+Size*i+j) = *(buffer + Size*maxRow+j);
                    *(buffer+Size*maxRow+j) = tmp;
                }
            }
            for(j=i+1;j<Size;j++){
                mulTmp = *(buffer+Size*j+i) / *(buffer+Size*i+i);
                for(k=i;k<Size;k++){

                    *(buffer+Size*j+k) = *(buffer+Size*j+k)-*(buffer + Size*i+k)*mulTmp;
                }
            }
        }
        mpfr_set_d(tmp1,1.0,GMP_RNDN);
        for(i=0;i<Size;i++){
            mpfr_set_d(tmp2,*(buffer+Size*i+i),GMP_RNDN);
            mpfr_mul(tmp1,tmp2,tmp1,GMP_RNDN);
            logAbsDet=logAbsDet+log10(fabs(*(buffer+Size*i+i)));
        }
        if(sign%2==0){
            mpfr_set(det,tmp1,GMP_RNDN);
            resDet=mpfr_get_d(det,GMP_RNDN);
        }
        else{
            mpfr_neg(tmp1,tmp1,GMP_RNDN);
            mpfr_set(det,tmp1,GMP_RNDN);

            resDet=mpfr_get_d(det,GMP_RNDN);
        }

        mtype = From_Worker;
        MPI_Send(&offset, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
        MPI_Send(&rows, 1, MPI_INT, MASTER, mtype, MPI_COMM_WORLD);
        MPI_Send(&resDet, 1, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);
    }
    MPI_Finalize();

    printf("det= %e\n",resDet);
    printf("log(abs(det)= %e\n",logAbsDet);

    END = clock();
    printf("%f\n",(END-START)/CLOCKS_PER_SEC);

    mpfr_clear(det);
    mpfr_clear(tmp1);
    mpfr_clear(tmp2);
    free(buffer);
    return 0;
}


#include <R.h>
#include <Rinternals.h>


void component_adj(double *ComponentLRT_prior, int *prior1, int *prior2, double *ComponentLRT_post, int *post2, int *rowindex, int *colindex, double *cv, int *n_adj, int *simu)
{
int i;
int j,x,nonsignificant;
for(i=0;i<*simu;++i){
   for(j=0;j<*post2;++j){
      nonsignificant=0;
      for(x=0;x<*n_adj;++x){
         if( *(ComponentLRT_prior + i * (*prior1) * (*prior2) + (*(colindex+x * (*post2)+j)-1) * (*prior1) + *(rowindex+x)-1) < *(cv+j * (*n_adj)+x) ){
             ++nonsignificant;
         }
         else{break;}
      }
      if(nonsignificant>1) *(ComponentLRT_post + i * (*post2) + j)=0;
   }                     
}
}

/*
void partition_adj(double *ComponentLRT, int *Comp1, int *Comp2, double *LRT, int *LRTrow, int *h, double *cv, int *cvlength, int *rowind, int *colind, int *simu)
{
   int x,significant;
   int i;
   for(i=0;i<*simu;++i){
      significant=*cvlength;
      for(x=0;x<*cvlength;++x){
         if( *(ComponentLRT+i* *Comp1 * *Comp2 +(*(colind+x)-1)* *Comp1 +*(rowind+x)-1)<*(cv+x) ) --significant;               
      }
      if(significant<*cvlength-1) *(LRT+i* *LRTrow +*h-1)=0;                                                      
   }                     
}


void sorting_order(double *vector, int *n, int *order){
    int i, j=0, k, taken[*n];
    double temp;
    for(i=0;i<*n;++i)taken[i]=0;
    for(i=0;i<*n;++i){
        for(k=0;k<*n;++k)if(taken[k]==0){j=k;break;}
        temp=vector[j];
        order[i]=j+1;
        for(k=j+1;k<*n;++k){
            if(taken[k]!=0)continue;
            else if(vector[k]<temp){
                order[i]=k+1;
                temp=vector[k];
            }
        }
        taken[order[i]-1]=1;
    }
}

double ConstraintMLE(double *ybar, int *trtnum, int *n){
    int i,j,order[*trtnum];
    double tempsum, tempn, meansq=0;
    ytilda=(double *)malloc((*trtnum+1)*sizeof(double));
    for(i=0;i<*trtnum+1;++i)ytilda[i]=*(ybar+i);
    sorting_order(ybar,trtnum,order);
    for(i=*trtnum-1;i>=0;++i){
        if(ybar[order[i]-1]<=ybar[*trtnum]) break;
        for(j=*trtnum-1,tempsum=ybar[*trtnum]*n[*trtnum],tempn=n[*trtnum];j>=i;--j){
            tempsum+=ybar[order[j]-1]*n[order[j]-1];
            tempn+=n[order[j]-1];
        }
        ybar[*trtnum]=tempsum/tempn;
        for(j=*trtnum-1;j>=i;--j)ybar[order[j]-1]=tempsum/tempn;
    }
    for(i=0;i<*trtnum+1;++i)meansq+=(ybar[i]-ytilda[i])*(ybar[i]-ytilda[i]);
    free(ytilda);
    return meansq;
}

void LR(double *label, int *dose, double *ybarP, double *ybarS, double *LRT, int *simu){
    int j,pj,sj,i;
    ycurrentP=(double *)malloc((*dose+1)*sizeof(double));
    ycurrentS=(double *)malloc((*dose+1)*sizeof(double));
    int *P,*S;

    for(j=0,pj=0,sj=0;j<*dose;++j){
        if(label[j]==1){
            if(P==NULL)P=(int *)malloc(1*sizeof(int));
            else P=(int *)realloc(P, (pj+1) * sizeof(int));
            P[pj]=j;++pj;
        }
        else if(label[j]==2){
            if(S==NULL)S=(int *)malloc(1*sizeof(int));
            else S=(int *)realloc(S, (sj+1) * sizeof(int));
            S[sj]=j;++sj;
        }
    }
    
    for(i=0;i<*simu;++i){
        for(j=0;j<=*dose;++j){ycurrentP[j]=*(ybarP+(*dose+1)*i+j);ycurrentS[j]=*(ybarS+(*dose+1)*i+j);}
        meansq
    }
}


void LRT(int* label, double* ybar, double* sigma, int* nu, int* part){
   int i;
   double s;
   for(i=0;i<*part;++i){
      s = *sigma * *sigma * *nu;

   }
}

int main(){
   double C_prior[12]={1.395022, 1.076493, 1.000000, 3.534913, 1.000000, 1.101371, 1.000000, 1.000000, 1.000000, 1.034383, 1.000000, 1.000000};
   long prior1=2, prior2=3;
   double C_post[6]={1.395022, 1.395022, 1.000000, 1.000000, 1.000000, 1.000000};
   long post2=3;
   long rowindex[2]={1,1};
   long colindex[6]={1,1,2,2,3,3};
   double cv[6]={6.51771,6.51771,6.51771,6.51771,6.51771,6.51771};
   long n_adj=2;
   long simu=2;
   consonant_adj(&C_prior[0],&prior1,&prior2,&C_post[0],&post2,&rowindex[0],&colindex[0],&cv[0],&n_adj,&simu);
   long i,j;
   for(i=0;i<6;++i){
         printf("%d,%f\t,%f,\t %f\n", i,*(C_prior+2*i),*(C_prior+2*i+1),*(C_post+i));    
   }
   return 0;
}
*/

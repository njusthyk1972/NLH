#include <math.h>
#include "mex.h"
#include <stdio.h>

// 声明输入变量
double *N1,*N2,*N3,*Ns;
float *Nstep;
double *im_p;
int    M,N;
// 声明输出变量
float *imr,*sumd;


void FindMin(float* OrgArray, int OrgNum, int* ResIndexArray, int ResNum)
{
    int i,j,k;
    float* ResArray = new float[ResNum];
    for ( i = 0; i < ResNum; i++ )  *(ResArray+i) = 10000000.0;

    for ( i = 0; i < OrgNum; i++ )
        {
            if ( *(OrgArray+i) >= *(ResArray+ResNum-1) )   continue;
            else if ( *(OrgArray+i) <= *ResArray )
                {
                    for ( j = ResNum-1; j > 0; j-- )
                        {
                            *(ResArray+j) = *(ResArray+j-1);
                            *(ResIndexArray+j) = *(ResIndexArray+j-1);
                        }
                   *ResArray = *(OrgArray+i);
                   *ResIndexArray = i;
                   continue;
                }
            else
                {
                    for ( j = 0; j < ResNum-1; j++ )
                        {
                            if ( *(ResArray+j) < *(OrgArray+i) && *(OrgArray+i) < *(ResArray+j+1) )
                                {
                                    for ( k = ResNum-1; k > j+1; k-- )
                                        {
                                            *(ResArray+k) = *(ResArray+k-1);
                                            *(ResIndexArray+k) = *(ResIndexArray+k-1);
                                        }
                                    *(ResArray+j+1) = *(OrgArray+i);
                                    *(ResIndexArray+j+1) = i;
                                    break;
                                }
                        }
                }
        }
    delete [] ResArray;
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //参数个数检查
    if ( nrhs != 6 ) 
		mexErrMsgTxt("6 input arguments required."); 
	if ( nlhs != 2 )
		mexErrMsgTxt("2 output arguments required.");
    // 定义输入变量
    N1       = mxGetPr(prhs[0]);
    N2       = mxGetPr(prhs[1]);
    N3       = mxGetPr(prhs[2]);
    Ns       = mxGetPr(prhs[3]);
    Nstep    = (float*)mxGetPr(prhs[4]);
    im_p     = mxGetPr(prhs[5]);
    M        = mxGetM(prhs[5]);
    N        = mxGetN(prhs[5]);
  
  // 定义输出变量
    plhs[0] = mxCreateNumericMatrix(M,N,mxGetClassID(prhs[4]),mxREAL);
    imr     = (float*)mxGetPr(plhs[0]);
    plhs[1] = mxCreateNumericMatrix(1,1,mxGetClassID(prhs[4]),mxREAL);
    sumd     = (float*)mxGetPr(plhs[1]);   
    
    int N12 = int(*N1);
    int N22 = int(*N2);
    int N32 = int(*N3);
    int Ns1 = int(*Ns);
    int Nstep1  = int(*Nstep);
    
        
    int i,j,k,l,r,u;
    float weight;
   // 中间矩阵空间定义
    // im_ps，M行M列
    float** im_ps;
    im_ps = new float*[M+2*Ns1];
    for ( i = 0; i < M+2*Ns1; i++ )   im_ps[i] = new float[N+2*Ns1]; 
     // im_r，M行M列
    float** im_r;
    im_r = new float*[M+2*Ns1];
    for ( i = 0; i < M+2*Ns1; i++ )   im_r[i] = new float[N+2*Ns1];
    
     // W，M行M列
    float** W;
    W = new float*[M+2*Ns1];
    for ( i = 0; i < M+2*Ns1; i++ )   W[i] = new float[N+2*Ns1]; 
    // Cood，N2行3列
   int** Cood;
   Cood = new int*[N22];
   for ( i = 0; i < N22; i++ )   Cood[i] = new int[3];
    // group，N1N12N22
    float*** group;
    group = new float**[N12];
    for ( i = 0; i < N12; i++ )
        {
            group[i] = new float*[N12];
            for ( j = 0; j < N12; j++ )
                group[i][j] = new float[N22];
        }
    float** groupb;
    groupb = new float*[N12*N12];
    for ( i = 0; i < N12*N12; i++ )   groupb[i] = new float[N22];
    
    float** groups;
    groups = new float*[N32];
    for ( i = 0; i < N32; i++ )       groups[i] = new float[N22];
    float** groupss;
    groupss = new float*[N32];
    for ( i = 0; i < N32; i++ )       groupss[i] = new float[N22];             
          
    float** groupss_t;
    groupss_t = new float*[N32];
    for ( i = 0; i < N32; i++ )   groupss_t[i] = new float[N22];
    float** groups_t;
    groups_t = new float*[N32];
    for ( i = 0; i < N32; i++ )   groups_t[i] = new float[N22];
    
    float** DD;
    DD = new float*[N12*N12];
    for ( i = 0; i < N12*N12; i++ )   DD[i] = new float[N22];
    
    float** WW;
    WW = new float*[N12*N12];
    for ( i = 0; i < N12*N12; i++ )   WW[i] = new float[N22];
    
    // im_n，N1行N1列
    float** im_n;
    im_n = new float*[N12];
    for ( i = 0; i < N12; i++ )   im_n[i] = new float[N12];
  
       
    // dis1，Ns*Ns-1行3列
    float** dis1;
    dis1 = new float*[(Ns1)*(Ns1)-1];
    for ( i = 0; i < (Ns1)*(Ns1)-1; i++ )   dis1[i] = new float[3];
    // dis2，N2-1行3列
    float** dis2;
    dis2 = new float*[N22-1];
    for ( i = 0; i < N22-1; i++ )   dis2[i] = new float[3];
    
    float** dis3;
    dis3 = new float*[N12*N12];
    for ( i = 0; i < N12*N12; i++ )   dis3[i] = new float[N12*N12];
    
     // dis4，N1*N1行2列
    float* dis4 = new float[N12*N12];
    
	// dis1Array
	float* dis1Array = new float[(Ns1)*(Ns1)-1];
	// IndexArray
    int* IndexArray = new int[N22-1];
    // dis1Array
	float* dis1Arrays = new float[N12*N12];
	// IndexArray
    int* IndexArrays = new int[N32];
    
     //d3
    int* dis5 = new int[N32];
    
    
    // x,y
    int x[409600],y[409600];
    // 输出变量初始化
	for ( i = 0; i < M*N; i++ )
		{
			*(imr+i) = 0.0;
		}
    for ( i = 0; i < 1; i++ )
		{
			*(sumd+i) = 0.0;
		}     
        
    
   for ( i = 0; i < M+2*Ns1; i++ )
          for ( j = 0; j < N+2*Ns1; j++ )
              {
              im_ps[i][j] = 0.0;
              im_r[i][j]  = 0.0;
              W[i][j]     = 0.0;
              }
   
       for ( i = Ns1; i < M+Ns1; i++ )
          for ( j = Ns1; j < N+Ns1; j++ )
              {
               im_ps[i][j] = *(im_p+(j-Ns1)*M+(i-Ns1));
              }

  // 周期延拓
    for ( i = 0; i < Ns1; i++ )
          for ( j = 0; j < N+2*Ns1; j++ )
              {
              im_ps[i][j] = im_ps[i+Ns1][j];
              im_ps[i+M+Ns1][j] = im_ps[i+M][j];
              }
    for ( i = 0; i < M+2*Ns1; i++ )
          for ( j = 0; j < Ns1; j++ )
              {
              im_ps[i][j] = im_ps[i][j+Ns1];
              im_ps[i][j+N+Ns1] = im_ps[i][j+N];
              }
 
    
    
                                                 
                         int h = 0;
                         
                          
                         
	// 算法主体
    for ( int a = Ns1/2+1; a <= M+Ns1; a += Nstep1 )
        for ( int b = Ns1/2+1; b <= N+Ns1; b += Nstep1 )
            {
            	for ( i = 0; i < N22; i++ )
                    for ( j = 0; j < 3; j++ )
                        Cood[i][j] = 0.0;  
                for ( i = 0; i < N12; i++ )
                    for ( j = 0; j < N12; j++ )
                        for ( k = 0; k < N22; k++ )
                            group[i][j][k] = 0.0;
               
                for ( i = a; i <= a+N12-1; i++ )
                    for ( j = b; j <= b+N12-1; j++ )
                        {
                            im_n[i-a][j-b] = im_ps[i-1][j-1];
                            group[i-a][j-b][0] = im_n[i-a][j-b];
                        }
                Cood[0][0] = 0.0;
                Cood[0][1] = a-1;
                Cood[0][2] = b-1;
                int Q = 0;
                for ( k = 1; k <= (Ns1-1)/2; k++ )   // 半径
                    {
                        for ( i = -k; i <= k-1; i++ )
                            {
                                x[Q] = a-k-1;
                                y[Q] = b+i-1;
                                Q++;
                            }
                        for ( i = -k; i <= k-1; i++ )
                            {
                                x[Q] = a+i-1;
                                y[Q] = b+k-1;
                                Q++;
                            }
                        for ( i = -k; i <= k-1; i++ )
                            {
                                x[Q] = a+k-1;
                                y[Q] = b-i-1;
                                Q++;
                            }
                        for ( i = -k; i <= k-1; i++ )
                            {
                                x[Q] = a-i-1;
                                y[Q] = b-k-1;
                                Q++;
                            }
                    }
//                 printf("%d \n",Q);  
                 int xk,yk; 
                for ( k = 0; k < Q; k++ )
                    {
                     float  sum  = 0.0,minu = 0.0; 
                          xk = x[k];
                          yk = y[k];
                       for ( i = xk; i < xk+N12; i++ )
                            for ( j = yk; j < yk+N12; j++ )
                          {
                              minu = im_ps[i][j]-im_n[i-xk][j-yk];
                                 sum += minu*minu;
                            }
                        
                        dis1[k][0] = xk;
                        dis1[k][1] = yk;
                        dis1[k][2] = sum;
                    }
                
                
                for ( i = 0; i < (Ns1)*(Ns1)-1; i++ )   dis1Array[i] = dis1[i][2];
                FindMin(dis1Array,(Ns1)*(Ns1)-1,IndexArray,N22-1);
               
                
                for ( i = 0; i < N22-1; i++ )
                    {
                        dis2[i][0] = dis1[IndexArray[i]][0];
                        dis2[i][1] = dis1[IndexArray[i]][1];
                        dis2[i][2] = dis1[IndexArray[i]][2];
                    }
                
           
                int x,y;                
                for ( k = 0; k < N22-1; k++ )
                    {    
                         x = dis2[k][0];
						 y = dis2[k][1];
                        for ( i = x; i <= x+N12-1; i++ )
                             for ( j = y; j <= y+N12-1; j++ )
							     group[i-x][j-y][k+1] = im_ps[i][j]; 
                        Cood[k+1][0] = k+1;
                        Cood[k+1][1] = x;
                        Cood[k+1][2] = y;
                    }
                
              
      
 //将每组中的所有块按列扫描给groupb         
                            
                      for ( k = 0; k < N12; k++ )
                           for ( j = 0; j < N12; j++ )
                                for ( i = 0; i < N22; i++ )
                                      groupb[k*N12+j][i] = group[k][j][i];
       
                             
        for (k = 0; k < N12*N12; k++)
               for ( i = 0; i < N12*N12; i++) 
                     dis3[i][k] = 0.0;  
                
       
               for ( k = 0; k < N12*N12; k++)
                      for ( i = 0; i < N12*N12; i++)
                            { 
                            float sum = 0.0,minu = 0.0;
                              for ( j = 0; j < N22; j++)
                            {
                                minu = groupb[k][j]-groupb[i][j];
                                sum += minu*minu;
                              }
                             dis3[k][i] = sum;
                           }
//           printf("%d,%d,%d,%d\n",dis3[0][0],dis3[0][1],dis3[5][9],dis3[10][11]); 
             
   
                
                
     for ( i = 0; i < N12*N12; i++ )
          for ( j = 0; j < N22; j++ ) 
             {     
                 DD[i][j] = 0.0;
                 WW[i][j] = 0.00000001;
             } 
 // 开始每一组内的行匹配                 
 for (k = 0; k < N12*N12; k++)    //开始循环
                {  
            
                  
    
     for ( i = 0; i < N32; i++ )
          for ( j = 0; j < N22; j++ ) 
              { 
                  groups[i][j] = 0.0;
                  groups_t[i][j] = 0.0;
                  groupss[i][j] = 0.0;
                  groupss_t[i][j] = 0.0;
              }  
           
     
           for (i = 0; i < N12*N12; i++)
                    dis4[i] = 0.0;
     
           for (i = 0; i < N12*N12; i++)
                    dis4[i] = dis3[i][k];
         

               
         for ( i = 0; i < N12*N12; i++ )   dis1Arrays[i] = dis4[i];   
             FindMin(dis1Arrays,N12*N12,IndexArrays,N32);
  
        float sum = 0.0; 
/////////////////////////////////////////////////////////////////////
        
      //     if (dis4[IndexArrays[1]] == 0.0)
      //        dis4[IndexArrays[1]] = 0.00001;
              
           sum += dis4[IndexArrays[1]]+0.00000001;
           sum = sqrt(sum/N22);
           *sumd += sum;
           h++;
/////////////////////////////////////////////////////////////////////          
           for ( i = 0; i < N32; i++ )
                   dis5[i] = IndexArrays[i];

                   
             
            for ( i = 0; i < N32; i++)
                  for ( j = 0; j < N22; j++)
                           groups[i][j] = groupb[dis5[i]][j];
  

                    
     
              /*   
               for ( i = 0; i < N32; i++ )
                    for ( j = 0; j < N22; j++ )
                    {
                         float sum = 0.0;
                                for ( r = 0; r < N22; r++ )
                                     sum += groups[i][r];
                                     if (sum != 0.0)
                                    { 
                                    DD[dis5[i]][j] += dis4[IndexArrays[1]];
                                    WW[dis5[i]][j] += 1.0;
                                     }
                                   else 
                                     { 
                                    DD[dis5[i]][j] += 0.0;
                                    WW[dis5[i]][j] += 1.0;
                                     }  
                                }
              */  
           
           for ( i = 0; i < N32; i++ )
                    for ( j = 0; j < N22; j++ )
                    {
                      DD[dis5[i]][j] += dis4[IndexArrays[1]];
                      WW[dis5[i]][j] += 1.0;
                                   
                     }
           
           
                 
       }
                
     /*   
       for ( i = 0; i < N12*N12; i++ )
                    for ( j = 0; j < N22; j++ )
                        if (WW[i][j] == 0.0)
                        WW[i][j] = 0.00000001;
     */           
         for ( i = 0; i < N22; i++ )
                 for ( k = 0; k < N12; k++ )
                           for ( j = 0; j < N12; j++ )
                                   group[k][j][i] = DD[k*N12+j][i]/ WW[k*N12+j][i];              
  

                for ( k = 0; k < N22; k++ )
                     for ( i = Cood[k][1]; i <= Cood[k][1]+N12-1; i++ )
                            for ( j = Cood[k][2]; j <= Cood[k][2]+N12-1; j++ )
                                {
                                 im_r[i][j] += group[i-Cood[k][1]][j-Cood[k][2]][k];
                                 W[i][j]   += 1.0f;
                                 } 
              
        }

  if ( h == 0.0)
     *sumd = 0.0;
    else                       
    *sumd /= h;    
 //   *sumd /= (((M/Nstep1)*(N/Nstep1))*8.0); 
   *sumd *= 255.0;
   *sumd /= 0.8;
   for ( i = Ns1; i < M+Ns1; i++ )
          for ( j = Ns1; j < N+Ns1; j++ )
              *(imr+(j-Ns1)*M+(i-Ns1)) = im_r[i][j] / W[i][j];  
      

    
    // 矩阵空间释放
    for ( i = 0; i < M+2*Ns1; i++ )   delete [] im_ps[i];
    delete [] im_ps;
    for ( i = 0; i < M+2*Ns1; i++ )   delete [] W[i];
    delete [] W;
    for ( i = 0; i < M+2*Ns1; i++ )   delete [] im_r[i];
    delete [] im_r;
    
    for ( i = 0; i < N22; i++ )       delete [] Cood[i];
    delete [] Cood;
    for ( i = 0; i < N12; i++ )
        for ( j = 0; j < N12; j++ )
            delete [] group[i][j];
    for ( i = 0; i < N12; i++ )  delete [] group[i];
    delete [] group;
     
    for ( i = 0; i < N12*N12; i++ )   delete [] groupb[i];
    delete [] groupb;
     for ( i = 0; i < N32; i++ )   delete [] groups[i];
    delete [] groups;
    for ( i = 0; i < N32; i++ )   delete [] groupss[i];
    delete [] groupss;
    for ( i = 0; i < N32; i++ )   delete [] groups_t[i];
    delete [] groups_t;
    for ( i = 0; i < N32; i++ )   delete [] groupss_t[i];
    delete [] groupss_t;
    for ( i = 0; i < N12*N12; i++ )   delete [] DD[i];
    delete [] DD;
    
    for ( i = 0; i < N12*N12; i++ )   delete [] WW[i];
    delete [] WW;
    
    for ( i = 0; i < N12; i++ )   delete [] im_n[i];
    delete [] im_n;
    for ( i = 0; i < N12; i++ )   delete [] dis1[i];
    delete [] dis1;
    for ( i = 0; i < N22-1; i++ )   delete [] dis2[i];
    delete [] dis2;
    for ( i = 0; i < N12*N12; i++ )   delete [] dis3[i];
    delete [] dis3;
    
    delete [] dis4;
    
    delete [] dis5;
    delete [] dis1Array;
    delete [] IndexArray;
    delete [] dis1Arrays;
    delete [] IndexArrays;
   
    
}

















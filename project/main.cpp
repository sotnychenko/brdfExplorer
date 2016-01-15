// my first program in C++
#include <iostream>
#include <cstdlib>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <eigen3/Eigen/Sparse>
#include "cnpy.h"

#include "redsvd.hpp"
#include "redsvdFile.hpp"
#include <dirent.h>



#define BRDF_SAMPLING_RES_THETA_H       90
#define BRDF_SAMPLING_RES_THETA_D       90
#define BRDF_SAMPLING_RES_PHI_D         360

#define EPS 0.01

using namespace std;
using namespace Eigen;

 float sgn(float val)
  {
    if(val<0.0) return -1.0;
 
  return 1.0;

}

 float* brdfData;
 bool write_brdf(const char* filename )
 {
     // read in the MERL BRDF data
     FILE *f = fopen(filename, "wb");
     if (!f)
             return false;

     int dims[3];
     dims[0]= BRDF_SAMPLING_RES_THETA_H;
     dims[1]= BRDF_SAMPLING_RES_THETA_D;
     dims[2]= BRDF_SAMPLING_RES_PHI_D / 2;
     int  numBRDFSamples=dims[0] * dims[1] * dims[2];

     if (fwrite(dims, sizeof(int), 3, f) != 3) {
         fprintf(stderr, "write error\n");
         fclose(f);
         return false;
     }

     double* brdf = new double[ numBRDFSamples * 3 ];

     for( int i = 0; i < numBRDFSamples; i++ )
     {

             brdf[i*3 + 0] = (brdfData[i*3 + 0])/(1.00/1500.0);
             brdf[i*3 + 1] = (brdfData[i*3 + 1])/(1.15/1500.0);
             brdf[i*3 + 2] = (brdfData[i*3 + 2])/(1.66/1500.0);
     }
    // cout<< brdf[0]<<endl<<brdf[1]<<endl<<brdf[2]<<endl<<brdf[3]<<endl<<brdf[4]<<endl<<brdf[5];
     double (&reshaped)[180][90][90][3] = *reinterpret_cast<double(*)[180][90][90][3]>(brdf);
     double* brdfFinal = new double[ numBRDFSamples * 3 ];
     int ind=0;

     for(int i=0; i<3;i++)
              for(int j=0; j<90;j++)
                    for(int k=0; k<90;k++)
                          for(int l=0; l<180;l++)
                               brdfFinal[ind++]=reshaped[l][j][k][i];




     // cout<< brdfFinal[0]<<endl<<brdfFinal[1]<<endl<<brdfFinal[2]<<endl<<brdfFinal[3]<<endl<<brdfFinal[4]<<endl<<brdfFinal[5];
     if (fwrite(brdfFinal, sizeof(double), 3* numBRDFSamples, f) != size_t(3* numBRDFSamples)) {
         fprintf(stderr, "write error\n");
         fclose(f);
     return false;
     }
     fclose(f);
     free(brdf);

 }

bool read_brdf( const char* filename )
{


    // read in the MERL BRDF data
    FILE *f = fopen(filename, "rb");
    if (!f)
            return false;

    int dims[3];
    if (fread(dims, sizeof(int), 3, f) != 3) {
        fprintf(stderr, "read error\n");
        fclose(f);
        return false;
    }
    int numBRDFSamples = dims[0] * dims[1] * dims[2];
    if (numBRDFSamples != BRDF_SAMPLING_RES_THETA_H *
                            BRDF_SAMPLING_RES_THETA_D *
                            BRDF_SAMPLING_RES_PHI_D / 2)
    {
        fprintf(stderr, "Dimensions don't match\n");
        fclose(f);
        return false;
    }

    // read the data
    double* brdf = (double*) malloc (sizeof(double)*3*numBRDFSamples);
    fread(brdf, sizeof(double), 3*numBRDFSamples, f);

    fclose(f);

    // now transform it to RGBA floats

    brdfData = new float[ numBRDFSamples * 3 ];
    double (&reshaped)[3][90][90][180] = *reinterpret_cast<double(*)[3][90][90][180]>(brdf);
    int ind;
 for(int i=0; i<3;i++){
     ind = 0;
          for(int j=0; j<180;j++)
                for(int k=0; k<90;k++)
                      for(int l=0; l<90;l++)
                      {

                           brdfData[ind*3 + i]=reshaped[i][k][l][j];
                           ind++;

                      }

    }
 for(int i=0; i<numBRDFSamples;i++)
 {
      brdfData[i*3 + 0]*=1.00/1500;
      brdfData[i*3 + 1]*=1.15/1500;
      brdfData[i*3 + 2]*=1.66/1500;
 }

//cout<<reshaped[0][85][0][13]<<endl;
 //cout<<reshaped[0][60][0][0]<<endl;
    // free( brdfData );
    //  brdfData = new float[ numBRDFSamples * 3 ];

   // for( int i = 0; i < 10; i++ ) cout<<reshaped[0][50][25][i]<<endl;
    // now we can dump the old data
    free( brdf );
   // free(reshaped);


    return true;
}

float*  reshape(double* CosineMap, bool* MaskMap,int mask_size,int Qsize, float* my_npy)
{
    float* result = new float[Qsize*3];
    int ind = 0;
    for(int i=0; i<mask_size;i++)
        if(MaskMap[i]) {
            result[ind*3 + 0] = my_npy[i*3+0]*CosineMap[ind];
            result[ind*3 + 1] = my_npy[i*3+1]*CosineMap[ind];
            result[ind*3 + 2] = my_npy[i*3+2]*CosineMap[ind];
            ind++;
        }
    

    return result;
}

void  MapBRDF(float* reshapedBRDF,float* median,int Qsize )
{
     for(int i=0; i<Qsize;i++)
     { reshapedBRDF[i*3+0]=  log((reshapedBRDF[i*3+0]+0.001)/(median[i]+0.001));
       reshapedBRDF[i*3+1]=  log((reshapedBRDF[i*3+1]+0.001)/(median[i]+0.001));
       reshapedBRDF[i*3+2]=  log((reshapedBRDF[i*3+2]+0.001)/(median[i]+0.001));
     }
}

void UnmapBRDF(double* CosineMap,float* mappedData,bool* MaskMap,float* median,int mask_size)

{
    int ind = 0;
     for(int i=0; i<mask_size;i++)
     if(MaskMap[i])
    {
         brdfData[i*3+0]=exp(mappedData[ind*3+0])*(median[ind]+0.001)-0.001;
         brdfData[i*3+1]=exp(mappedData[ind*3+1])*(median[ind]+0.001)-0.001;
         brdfData[i*3+2]=exp(mappedData[ind*3+2])*(median[ind]+0.001)-0.001;
        
          brdfData[i*3+0]/=CosineMap[ind];
          brdfData[i*3+1]/=CosineMap[ind];
          brdfData[i*3+2]/=CosineMap[ind];
             ind++;
    }
     else
         brdfData[i*3 + 0] = brdfData[i*3 + 1]=brdfData[i*3 + 2]=-1.0 ;

}
void rgb2labv2( float R, float G, float B, float & l_s, float &a_s, float &b_s )
{

    float var_R=R;
    float var_G=G;
    float var_B=B;
    


    //Observer. = 2°, Illuminant = D65
    float X = var_R * 0.4124 + var_G * 0.3576 + var_B * 0.1805;
    float Y = var_R * 0.2126 + var_G * 0.7152 + var_B * 0.0722;
    float Z = var_R * 0.0193 + var_G * 0.1192 + var_B * 0.9505;


    float var_X = X / 0.95047 ;         //ref_X =  95.047   Observer= 2°, Illuminant= D65
    float var_Y = Y;          //ref_Y = 100.000
    float var_Z = Z / 1.08883;          //ref_Z = 108.883

    if ( var_X > 0.008856 ) var_X = pow(var_X , ( 1./3. ) );
    else                    var_X = ( 7.787 * var_X ) + ( 16. / 116. );
    if ( var_Y > 0.008856 ) var_Y = pow(var_Y , ( 1./3. ));
    else                    var_Y = ( 7.787 * var_Y ) + ( 16. / 116. );
    if ( var_Z > 0.008856 ) var_Z = pow(var_Z , ( 1./3. ));
    else                    var_Z = ( 7.787 * var_Z ) + ( 16. / 116. );

    l_s = ( 116. * var_Y ) - 16.;
    a_s = 500. * ( var_X - var_Y );
    b_s = 200. * ( var_Y - var_Z );


}

//http://www.easyrgb.com/index.php?X=MATH&H=01#text1
void lab2rgbv2( float l_s, float a_s, float b_s, float& R, float& G, float& B )
{
    float var_Y = ( l_s + 16. ) / 116.;
    float var_X = a_s / 500. + var_Y;
    float var_Z = var_Y - b_s / 200.;

    if ( pow(var_Y,3) > 0.008856 ) var_Y = pow(var_Y,3);
    else                      var_Y = ( var_Y - 16. / 116. ) / 7.787;
    if ( pow(var_X,3) > 0.008856 ) var_X = pow(var_X,3);
    else                      var_X = ( var_X - 16. / 116. ) / 7.787;
    if ( pow(var_Z,3) > 0.008856 ) var_Z = pow(var_Z,3);
    else                      var_Z = ( var_Z - 16. / 116. ) / 7.787;

    float X = 95.047 * var_X ;    //ref_X =  95.047     Observer= 2°, Illuminant= D65
    float Y = 100.000 * var_Y  ;   //ref_Y = 100.000
    float Z = 108.883 * var_Z ;    //ref_Z = 108.883


    var_X = X / 100. ;       //X from 0 to  95.047      (Observer = 2°, Illuminant = D65)
    var_Y = Y / 100. ;       //Y from 0 to 100.000
    var_Z = Z / 100. ;      //Z from 0 to 108.883

     R = var_X *  3.2406 + var_Y * -1.5372 + var_Z * -0.4986;
     G = var_X * -0.9689 + var_Y *  1.8758 + var_Z *  0.0415;
     B = var_X *  0.0557 + var_Y * -0.2040 + var_Z *  1.0570;

 
}
// using http://www.easyrgb.com/index.php?X=MATH&H=01#text1
void rgb2lab( float var_R, float var_G, float var_B, float & l_s, float &a_s, float &b_s )
{

  /*  if ( var_R > 0.04045 ) var_R = pow( (( var_R + 0.055 ) / 1.055 ), 2.4 );
    else                   var_R = var_R / 12.92;
    if ( var_G > 0.04045 ) var_G = pow( ( ( var_G + 0.055 ) / 1.055 ), 2.4);
    else                   var_G = var_G / 12.92;
    if ( var_B > 0.04045 ) var_B = pow( ( ( var_B + 0.055 ) / 1.055 ), 2.4);
    else                   var_B = var_B / 12.92;

    var_R = var_R * 100.;
    var_G = var_G * 100.;
    var_B = var_B * 100.;
*/

                
                

    //Observer. = 2°, Illuminant = D65
    float X = var_R * 0.412453 + var_G * 0.357580 + var_B * 0.180423;
    float Y = var_R * 0.212671 + var_G * 0.715160 + var_B * 0.072169;
    float Z = var_R * 0.019334 + var_G * 0.119193 + var_B * 0.950227;


    float var_X = X / 0.95046 ;         //ref_X =  95.047   Observer= 2°, Illuminant= D65
    float var_Y = Y;          //ref_Y = 100.000
    float var_Z = Z /  1.088754;          //ref_Z = 108.883
     float var_Yf=var_Y;
    if ( var_X > 0.008856 ) var_X = sgn(var_X)*pow(abs(var_X) , ( 1./3. ) );
    else                    var_X = ( 7.787 * var_X ) + ( 16. / 116. );
    if ( var_Y > 0.008856 ) var_Y = sgn(var_Y)*pow(abs(var_Y ), ( 1./3. ));
    else                    var_Y = ( 7.787 * var_Y ) + ( 16. / 116. );
    if ( var_Z > 0.008856 ) var_Z = sgn(var_Z)*pow(abs(var_Z) , ( 1./3. ));
    else                    var_Z = ( 7.787 * var_Z ) + ( 16. / 116. );

     if ( var_Yf > 0.008856 ) l_s = ( 116.0 * var_Y ) - 16.0;
      else l_s  = (903.3*var_Yf);
    a_s = 500. * ( var_X - var_Y );
    b_s = 200. * ( var_Y - var_Z );


}

//http://www.easyrgb.com/index.php?X=MATH&H=01#text1
void lab2rgb( float l_s, float a_s, float b_s, float& R, float& G, float& B )
{
    float  T1 = 0.008856;
    float T2 = 0.206893;
    
   float fY = pow((l_s + 16.0)/116.0,1.0/3.0);
   
   float fY_temp=fY;
   float Y;
    if ( fY > T1 ) Y = pow((l_s+16.0) /116.0, 3.0);
    else                    Y = (l_s/903.3);
    
   // 
   
   if ( fY_temp > T1 ) fY =  sgn( Y)*pow(abs( Y) , ( 1./3. ) );
    else                    fY =  ( 7.787 * Y ) + ( 16. / 116. );
    

    
    float fX = a_s/500.0 + fY;
    float X;

    
    if ( fX > T2 )  X=pow(fX , 3.0 );
    else                 X = (fX - 16.0/116.0)/7.787;
    
 

    
   float fZ = fY - b_s/200.0;
    float Z;
   if ( fX > T2 )  Z=pow(fZ , 3.0 );
    else                 Z = (fZ - 16.0/116.0)/7.787;
    
  
    
    X = X * 0.950456;
    Z = Z * 1.088754;
    
    
/*
    var_X = X / 100. ;       //X from 0 to  95.047      (Observer = 2°, Illuminant = D65)
    var_Y = Y / 100. ;       //Y from 0 to 100.000
    var_Z = Z / 100. ;      //Z from 0 to 108.883
*/

     R = X *  3.240479 + Y * -1.537150 + Z * -0.498535;
     G = X * -0.969256 + Y *  1.875992 + Z *  0.041556;
     B = X *  0.055648 + Y * -0.204043 + Z *  1.057311;
/*
    if ( var_R > 0.0031308 ) var_R = 1.055 * pow(var_R , ( 1 / 2.4 ))  - 0.055;
    else                     var_R = 12.92 * var_R;
    if ( var_G > 0.0031308 ) var_G = 1.055 * pow(var_G , ( 1 / 2.4 ) )  - 0.055;
    else                     var_G = 12.92 * var_G;
    if ( var_B > 0.0031308 ) var_B = 1.055 * pow( var_B , ( 1 / 2.4 ) ) - 0.055;
    else                     var_B = 12.92 * var_B;
*/

}

MatrixXf  rgb2Lab(MatrixXf proj)
{    
   MatrixXf Lab;
   Lab.resize(proj.rows(),proj.cols());
       
    for(int i=0; i<proj.rows();i++)
            {    for(int j=0; j<proj.cols();j++)
          //  cout<<proj(i,j)<<" ";
             rgb2lab(proj(i,0),proj(i,1),proj(i,2),Lab(i,0),Lab(i,1),Lab(i,2));
             }
          
    return Lab;             
 }
float* getRBFActivations(Matrix<float, Dynamic, Dynamic,RowMajor>& Centers, float* betas, float* input)
  {
  
        float* sqrdDists = new float[Centers.cols()];
           for(int i=0; i<Centers.cols();i++)
             { 
              sqrdDists[i]=0.0;
             for(int j=0; j<Centers.rows();j++)
    
                sqrdDists[i]+= (Centers(j,i)-input[j])*(Centers(j,i)-input[j]);
                
                
            //  cout<<sqrdDists[i]<<" ";
              
              //cout<<endl;
               }
                    
 for(int i=0; i<Centers.cols();i++)
    sqrdDists[i] = exp(-1.0*betas[i]* sqrdDists[i]);     

 
   return   sqrdDists;          
           
  }
  float  evaluateFuncApproxRBFN(Matrix<float, Dynamic, Dynamic,RowMajor>& Centers, float* betas, float* Theta, bool normalize, float* input)
  {
  
  float result =0.0;
 float* phis= getRBFActivations(Centers, betas, input);
      

if (normalize)
 { 
 
 float sum =0.0;
  for(int i=0; i<Centers.cols(); i++) sum+=phis[i];
  for(int i=0; i<Centers.cols(); i++) phis[i] = phis[i] / sum;
   }
   
   result+=Theta[0];            
 for(int i=0; i<Centers.cols();i++) result+=phis[i]*Theta[i+1];
  
 
// cout<<result<<endl;
 
return result;

 
  }


float normVec(float* x,Matrix<float, Dynamic, Dynamic,RowMajor> &Centers, int row)

{
float result=0.0;
for(int i=0; i<Centers.rows(); i++) 
   result+=(x[i]-Centers(i,row))*(x[i]-Centers(i,row));

return result;
}
void ProjectToPCSpaceRBN(float* data,float* PCs,float* relativeOffset,int Qsize)

{
    int numC = 5;
    Matrix<float, Dynamic, Dynamic,RowMajor> b;
    b.resize(Qsize,3);
    for(int i=0; i<Qsize;i++)
    {
         b.data()[3*i+0]=data[3*i+0]-relativeOffset[i];
         b.data()[3*i+1]=data[3*i+1]-relativeOffset[i];
         b.data()[3*i+2]=data[3*i+2]-relativeOffset[i];
    }




    //cout<<b(0,0)<<endl<<b(0,1)<<endl<<b(0,2)<<endl;


     Matrix<float, Dynamic, Dynamic,RowMajor> A;
     A.resize(Qsize,numC);
    for(int i=0; i<Qsize;i++)
        for(int j=0; j<numC;j++)
        A.data()[numC*i+j]=PCs[numC*i+j];


    //  cout<<A(0,0)<<endl<<A(0,1)<<endl<<A(0,2)<<endl;


    REDSVD::RedSVD SVD(A,numC);

       MatrixXf Ut= SVD.matrixU().transpose();
       MatrixXf V = SVD.matrixV().transpose();

       MatrixXf Sinv;
       Sinv.resize(numC,numC);
       Sinv.setZero();
       for(int i=0; i<numC; i++) Sinv(i,i) = SVD.singularValues()[i]/( SVD.singularValues()[i]* SVD.singularValues()[i]);



          MatrixXf proj= V*Sinv*Ut*b;
   /*       
         proj(0,0)= 0.023970775;	proj(0,1)=0.024991909;	proj(0,2)=0.032395281;
 proj(1,0)=-0.053679522;	proj(1,1)=-0.026844800;	proj(1,2)=0.034950085;
 proj(2,0)=-0.11060172;	proj(2,1)=-0.10956725;	proj(2,2)=-0.091819309;
 proj(3,0)=0.13308536	;    proj(3,1)=0.12815318;	proj(3,2)=0.077729017;
 proj(4,0)=-0.072992906;	proj(4,1)=-0.077936530;	proj(4,2)=-0.069622166;
*/
    
  MatrixXf proj_Lab= rgb2Lab(proj);
          float* xnew = new float[proj_Lab.rows()];
          
        /*  float L,A_,B;
          rgb2lab(-0.05367952, -0.0268448,   0.03495008, L,A_,B);
          cout<<L<<endl;
          cout<<A_<<endl;
          cout<<B<<endl;
          */
         // proj_Lab(1,2)=-103.46116966;
           for(int i=0; i<proj_Lab.rows();i++){
        //   cout<<proj_Lab(i,2)<<endl;
           xnew[i]=proj_Lab(i,0)/100.0;
           
           }
           
           
     /*  for(int i=0; i<proj.rows();i++)
             { for(int j=0; j<proj.cols();j++)
                 cout<<proj_Lab(i,j)<<" ";
               
               cout<<endl;
               }*/
               
               DIR *dir;
struct dirent *ent;
if ((dir = opendir ("trained_RBFN_npz\\")) != NULL) {
  /* print all the files and directories within directory */
  while ((ent = readdir (dir)) != NULL) {
//    printf ("%s\n", ent->d_name);
    
    std::string str(ent->d_name);
    
    if(str.substr(str.find_last_of(".") + 1) == "npz") {
    printf ("%s\n", ent->d_name);
  }
  }
  closedir (dir);
} else {
  /* could not open directory */
  perror ("");

}
      
      cnpy::npz_t my_npz = cnpy::npz_load("trained_RBFN_npz\\RBFN_att_12_N_010_sigma_010_10-Jan-2016.npz");
       cnpy::NpyArray arr_mv1 = my_npz["Centers"];
       double* mv1 = reinterpret_cast<double*>(arr_mv1.data);
       Matrix<float, Dynamic, Dynamic,RowMajor> Centers;
       Centers.resize(arr_mv1.shape[1],arr_mv1.shape[0]);
    for(int i=0; i<arr_mv1.shape[0]*arr_mv1.shape[1];i++)
        Centers.data()[i]=mv1[i];
       
     
      arr_mv1 = my_npz["betas"];
      double* betasTemp = reinterpret_cast<double*>(arr_mv1.data);
      float* betas = new float[arr_mv1.shape[0]*arr_mv1.shape[1]];
      for(int i=0; i<arr_mv1.shape[0]*arr_mv1.shape[1];i++) betas[i]=betasTemp[i];
      
    //  for(int i=0; i<5;i++) cout<<betas[i]<<endl;
       
         arr_mv1 = my_npz["Theta"];
    
       double* ThetaTemp = reinterpret_cast<double*>(arr_mv1.data);
      float* Theta = new float[arr_mv1.shape[0]*arr_mv1.shape[1]];
      for(int i=0; i<arr_mv1.shape[0]*arr_mv1.shape[1];i++) Theta[i]=ThetaTemp[i];
      // for(int i=0; i<11;i++) cout<<Theta[i]<<endl;
       
      float ynew = evaluateFuncApproxRBFN(Centers, betas, Theta, true, xnew);
      
      float yobj=1.0;
      float alpha=0.05;
      
  
     

    while (abs(ynew-yobj)>EPS)
    {
 
    float* grad = new float[Centers.rows()];
   for(int i=0; i<Centers.rows();i++)
   {
   
   grad[i]=0.0;
   for(int j=0;j<Centers.cols(); j++)
   
   {
   
     grad[i] +=-2.0*betas[j]*Theta[j+1]*(xnew[i]-Centers(i,j))*exp(-1.0*betas[j]*normVec(xnew,Centers,j));
    
   }

   
   }

 float norm =0.0;
  for(int i=0; i<Centers.rows(); i++) norm+=grad[i]*grad[i];
  
  norm=sqrt(norm);
  
  
  for(int i=0; i<Centers.rows(); i++)
   { 
  grad[i] /= norm;  

  
   xnew[i] -= sgn(ynew-yobj)*alpha*grad[i];
  
   }
 
    ynew = evaluateFuncApproxRBFN(Centers, betas, Theta, true, xnew);
  
    
   // cout<<ynew<<endl;
}


         for(int i=0; i<proj.rows();i++)
             {
             //  cout<<xnew[i]*100.0<<endl;
                 lab2rgb(xnew[i]*100.0,proj_Lab(i,1),proj_Lab(i,2),proj(i,0),proj(i,1),proj(i,2));
                 
               //   for(int j=0; j<proj.cols();j++)
                 //    cout<<proj(i,j)<<" ";
                     
              // cout<<endl;
               }
         
         
         
      
/*
           for(int i=0; i<Centers.cols();i++)
             { for(int j=0; j<Centers.rows();j++)
                 cout<<Centers(j,i)<<" ";
               
               cout<<endl;
               }*/
               
               

      
   //  float yobj= 1.0;
  //  float ynew = evaluateFuncApproxRBFN(my_npz, true, xini);



          
          
        /*   cnpy::NpyArray my_npz = cnpy::npz_load("trained_RBFN_npz\\RBFN_att_01_N_010_sigma_010_10-Jan-2016.npz","sigma");
             cout<<" shape "<<my_npz.shape[0]<<" "<<my_npz.shape[1]<<endl;
             uint8_t* data =reinterpret_cast<uint8_t*>(my_npz.data);
             int suka=data[0];
             cout<<suka<<endl;*/
                   
         MatrixXf recon =  A*proj;
         
         for(int i=0; i<Qsize;i++)
             for(int j=0;j<3;j++)
                 data[i*3+j]= recon(i,j)+relativeOffset[i];

 free(mv1);
}
void ProjectToPCSpace(float* data,float* PCs,float* relativeOffset,int Qsize)

{
    int numC = 5;
    Matrix<float, Dynamic, Dynamic,RowMajor> b;
    b.resize(Qsize,3);
    for(int i=0; i<Qsize;i++)
    {
         b.data()[3*i+0]=data[3*i+0]-relativeOffset[i];
         b.data()[3*i+1]=data[3*i+1]-relativeOffset[i];
         b.data()[3*i+2]=data[3*i+2]-relativeOffset[i];
    }




    //cout<<b(0,0)<<endl<<b(0,1)<<endl<<b(0,2)<<endl;


     Matrix<float, Dynamic, Dynamic,RowMajor> A;
     A.resize(Qsize,numC);
    for(int i=0; i<Qsize;i++)
        for(int j=0; j<numC;j++)
        A.data()[numC*i+j]=PCs[numC*i+j];


    //  cout<<A(0,0)<<endl<<A(0,1)<<endl<<A(0,2)<<endl;


    REDSVD::RedSVD SVD(A,numC);

       MatrixXf Ut= SVD.matrixU().transpose();
       MatrixXf V = SVD.matrixV().transpose();

       MatrixXf Sinv;
       Sinv.resize(numC,numC);
       Sinv.setZero();
       for(int i=0; i<numC; i++) Sinv(i,i) = SVD.singularValues()[i]/( SVD.singularValues()[i]* SVD.singularValues()[i]);

         MatrixXf recon =  A*(V*Sinv*Ut*b);

         for(int i=0; i<Qsize;i++)
             for(int j=0;j<3;j++)
                 data[i*3+j]= recon(i,j)+relativeOffset[i];






}


int main(int argc, char *argv[])
{
    const char *filename = argv[1];


        // read brdf
        if (!read_brdf(filename))
        {
            fprintf(stderr, "Error reading %s\n", filename);
            exit(1);
        }
        int mask_size;
        int Qsize;





          cnpy::NpyArray my_npy = cnpy::npy_load("MaskMap.npy");
          
          cout<<"mask shape "<<my_npy.shape[0]<<" "<<my_npy.shape[1]<<endl;
          mask_size = my_npy.shape[0];
          bool* MaskMap = reinterpret_cast<bool*>(my_npy.data);


          my_npy = cnpy::npy_load("Median.npy");
           cout<<"median shape "<<my_npy.shape[0]<<" "<<my_npy.shape[1]<<endl;
          float* median =reinterpret_cast<float*>(my_npy.data);
          Qsize = my_npy.shape[0];

          my_npy = cnpy::npy_load("Q.npy");
           cout<<"ScaledEigenvectors shape "<<my_npy.shape[0]<<" "<<my_npy.shape[1]<<endl;
          float* Q  = reinterpret_cast<float*>(my_npy.data);
          
          my_npy = cnpy::npy_load("CosineMap.npy");
           cout<<"CosineMap "<<my_npy.shape[0]<<" "<<my_npy.shape[1]<<endl;
          double* CosineMap  = reinterpret_cast<double*>(my_npy.data);
          
      
          float* reshapedBRDF = reshape(CosineMap,MaskMap,mask_size,Qsize,brdfData);

          my_npy = cnpy::npy_load("RelativeOffset.npy");
           cout<<"RelativeOffset "<<my_npy.shape[0]<<" "<<my_npy.shape[1]<<endl;
          float* RelativeOffset  = reinterpret_cast<float*>(my_npy.data);
           
   
          

          MapBRDF(reshapedBRDF,median,Qsize);



           ProjectToPCSpaceRBN(reshapedBRDF,Q,RelativeOffset,Qsize);

          UnmapBRDF(CosineMap,reshapedBRDF,MaskMap,median,mask_size);
          //for(int i=0;i<12;i++)
         // cout<<brdfData[i]<<endl;


          if (! write_brdf("outputBRDF.binary"))
          {
              fprintf(stderr, "Error writing %s\n", "outputBRDF.binary");
              exit(1);
          }

        return 0;
        
   

}

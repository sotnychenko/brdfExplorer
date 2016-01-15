/*
Copyright Disney Enterprises, Inc. All rights reserved.

This license governs use of the accompanying software. If you use the software, you
accept this license. If you do not accept the license, do not use the software.

1. Definitions
The terms "reproduce," "reproduction," "derivative works," and "distribution" have
the same meaning here as under U.S. copyright law. A "contribution" is the original
software, or any additions or changes to the software. A "contributor" is any person
that distributes its contribution under this license. "Licensed patents" are a
contributor's patent claims that read directly on its contribution.

2. Grant of Rights
(A) Copyright Grant- Subject to the terms of this license, including the license
conditions and limitations in section 3, each contributor grants you a non-exclusive,
worldwide, royalty-free copyright license to reproduce its contribution, prepare
derivative works of its contribution, and distribute its contribution or any derivative
works that you create.
(B) Patent Grant- Subject to the terms of this license, including the license
conditions and limitations in section 3, each contributor grants you a non-exclusive,
worldwide, royalty-free license under its licensed patents to make, have made,
use, sell, offer for sale, import, and/or otherwise dispose of its contribution in the
software or derivative works of the contribution in the software.

3. Conditions and Limitations
(A) No Trademark License- This license does not grant you rights to use any
contributors' name, logo, or trademarks.
(B) If you bring a patent claim against any contributor over patents that you claim
are infringed by the software, your patent license from such contributor to the
software ends automatically.
(C) If you distribute any portion of the software, you must retain all copyright,
patent, trademark, and attribution notices that are present in the software.
(D) If you distribute any portion of the software in source code form, you may do
so only under this license by including a complete copy of this license with your
distribution. If you distribute any portion of the software in compiled or object code
form, you may only do so under a license that complies with this license.
(E) The software is licensed "as-is." You bear the risk of using it. The contributors
give no express warranties, guarantees or conditions. You may have additional
consumer rights under your local laws which this license cannot change.
To the extent permitted under your local laws, the contributors exclude the
implied warranties of merchantability, fitness for a particular purpose and non-
infringement.
*/

#include <cstdlib>
#include <string>
#include <fstream>
#include <string.h>
#include "BRDFMeasuredMERL.h"
#include "DGLShader.h"
#include "Paths.h"

#include "BRDFBase.h"
#include "projectToPCA.h"

#define BRDF_SAMPLING_RES_THETA_H       90
#define BRDF_SAMPLING_RES_THETA_D       90
#define BRDF_SAMPLING_RES_PHI_D         360


#include <eigen3/Eigen/Sparse>
#include"cnpy.h"

#include "redsvd.hpp"
#include "redsvdFile.hpp"



using namespace std;
using namespace Eigen;

#include <iostream>
#define EPS 0.01

using namespace std;

 float BRDFMeasuredMERL::sgn(float val)
  {
    if(val<0.0) return -1.0;
 
  return 1.0;

}

BRDFMeasuredMERL::BRDFMeasuredMERL()
                 : brdfData(NULL)
{
    std::string path = (getShaderTemplatesPath() + "measured.func");
  brdfParam=new brdfMERLparam;
    // read the shader
    std::ifstream ifs( path.c_str() );
    std::string temp;
    while( getline( ifs, temp ) )
        brdfFunction += (temp + "\n");
}

 BRDFMeasuredMERL::BRDFMeasuredMERL(float* data, string fileName,int numSamples)
 
 {
  name = fileName;
 numBRDFSamples =numSamples;
 brdfData=data;
 std::string path = (getShaderTemplatesPath() + "measured.func");

    // read the shader
    std::ifstream ifs( path.c_str() );
    std::string temp;
    while( getline( ifs, temp ) )
        brdfFunction += (temp + "\n");
 
 }
 void  BRDFMeasuredMERL::allocBRDF()
 {
  brdfData = new float[ numBRDFSamples * 3 ];
 
 }
void BRDFMeasuredMERL::setName(string n)
{name=n;}


BRDFMeasuredMERL::~BRDFMeasuredMERL()
{
    glBindBuffer(GL_TEXTURE_BUFFER_EXT, tbo);
    glDeleteBuffers( 1, &tbo);
}


std::string BRDFMeasuredMERL::getBRDFFunction()
{
    return brdfFunction;
}

bool BRDFMeasuredMERL::loadPCAprojection( const char* filename )

{



return true;
}
bool BRDFMeasuredMERL::loadMERLData( const char* filename )
{
    // the BRDF's name is just the filename
    name = std::string(filename);

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
    numBRDFSamples = dims[0] * dims[1] * dims[2];
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
    if (fread(brdf, sizeof(double), 3*numBRDFSamples, f) != size_t(3*numBRDFSamples)) {
        fprintf(stderr, "read error\n");
        fclose(f);
	return false;
    }
    fclose(f);

    // now transform it to RGBA floats
    brdfData = new float[ numBRDFSamples * 3 ];
    for( int i = 0; i < numBRDFSamples; i++ )
    {
            brdfData[i*3 + 0] = brdf[i*3 + 0];
            brdfData[i*3 + 1] = brdf[i*3 + 1];
            brdfData[i*3 + 2] = brdf[i*3 + 2];
    }

    // now we can dump the old data
    free( brdf );

    return true;
}



void BRDFMeasuredMERL::initGL()
{
    if( initializedGL )
        return;
 
    
    // create buffer object
    glGenBuffers(1, &tbo);
    glBindBuffer(GL_TEXTURE_BUFFER_EXT, tbo);
    
    
    
    // initialize buffer object
    unsigned int numBytes = numBRDFSamples * 3 * sizeof(float);
    //printf( "size = %d bytes (%f megs)\n", numBytes, float(numBytes) / 1048576.0f );
    glBufferData( GL_TEXTURE_BUFFER_EXT, numBytes, 0, GL_DYNAMIC_DRAW );
    
    //tex
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_BUFFER_EXT, tex);
    glTexBufferEXT(GL_TEXTURE_BUFFER_EXT, GL_INTENSITY32F_ARB, tbo);
    glBindBuffer(GL_TEXTURE_BUFFER_EXT, 0);
    
    
    
    glBindBuffer(GL_TEXTURE_BUFFER_EXT, tbo);
     float * p = (float*)glMapBuffer( GL_TEXTURE_BUFFER_EXT, GL_WRITE_ONLY );
    
   //
    memcpy( p, brdfData, numBytes );
    glUnmapBuffer(GL_TEXTURE_BUFFER_EXT);
    glBindBuffer(GL_TEXTURE_BUFFER_EXT, 0);
  
    
    delete[] brdfData;
    brdfData = NULL;
    
    initializedGL = true;
}


void BRDFMeasuredMERL::adjustShaderPreRender( DGLShader* shader )
{
    shader->setUniformTexture( "measuredData", tex, GL_TEXTURE_BUFFER_EXT );

    BRDFBase::adjustShaderPreRender( shader );
}


void  BRDFMeasuredMERL::reshapeFinal( )
 {
    

     int dims[3];
     dims[0]= BRDF_SAMPLING_RES_THETA_H;
     dims[1]= BRDF_SAMPLING_RES_THETA_D;
     dims[2]= BRDF_SAMPLING_RES_PHI_D / 2;
     int  numBRDFSamples=dims[0] * dims[1] * dims[2];



     float* brdf = new float[ numBRDFSamples * 3 ];

     for( int i = 0; i < numBRDFSamples; i++ )
     {

             brdf[i*3 + 0] = (brdfData[i*3 + 0])/(1.00/1500.0);
             brdf[i*3 + 1] = (brdfData[i*3 + 1])/(1.15/1500.0);
             brdf[i*3 + 2] = (brdfData[i*3 + 2])/(1.66/1500.0);
     }
    // cout<< brdf[0]<<endl<<brdf[1]<<endl<<brdf[2]<<endl<<brdf[3]<<endl<<brdf[4]<<endl<<brdf[5];
     float (&reshaped)[180][90][90][3] = *reinterpret_cast<float(*)[180][90][90][3]>(brdf);
     
     int ind=0;

     for(int i=0; i<3;i++)
              for(int j=0; j<90;j++)
                    for(int k=0; k<90;k++)
                          for(int l=0; l<180;l++)
                               brdfData[ind++]=reshaped[l][j][k][i];

    free(brdf);


  

 }

bool  BRDFMeasuredMERL::read_brdf( const char* filename )
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
     numBRDFSamples = dims[0] * dims[1] * dims[2];
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
float*   BRDFMeasuredMERL::reshape(double* CosineMap,bool* MaskMap,int mask_size,int Qsize, float* my_npy)
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

void   BRDFMeasuredMERL::MapBRDF(float* reshapedBRDF,float* median,int Qsize )
{
     for(int i=0; i<Qsize;i++)
     { reshapedBRDF[i*3+0]=  log((reshapedBRDF[i*3+0]+0.001)/(median[i]+0.001));
       reshapedBRDF[i*3+1]=  log((reshapedBRDF[i*3+1]+0.001)/(median[i]+0.001));
       reshapedBRDF[i*3+2]=  log((reshapedBRDF[i*3+2]+0.001)/(median[i]+0.001));
     }
}
void BRDFMeasuredMERL::rgb2labv2( float R, float G, float B, float & l_s, float &a_s, float &b_s )
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
void BRDFMeasuredMERL::lab2rgbv2( float l_s, float a_s, float b_s, float& R, float& G, float& B )
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
void  BRDFMeasuredMERL::rgb2lab( float var_R, float var_G, float var_B, float & l_s, float &a_s, float &b_s )
{


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
void  BRDFMeasuredMERL::lab2rgb( float l_s, float a_s, float b_s, float& R, float& G, float& B )
{
    float  T1 = 0.008856;
    float T2 = 0.206893;
    
   float fY = pow((l_s + 16.0)/116.0,1.0/3.0);
   
   float fY_temp=fY;
   float Y;
    if ( fY > T1 ) Y = pow((l_s+16.0) /116.0, 3.0);
    else                    Y = (l_s/903.3);
   
   
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
    
    

     R = X *  3.240479 + Y * -1.537150 + Z * -0.498535;
     G = X * -0.969256 + Y *  1.875992 + Z *  0.041556;
     B = X *  0.055648 + Y * -0.204043 + Z *  1.057311;

}

MatrixXf   BRDFMeasuredMERL::rgb2Lab(MatrixXf proj)
{    
   MatrixXf Lab;
   Lab.resize(proj.rows(),proj.cols());
       
    for(int i=0; i<proj.rows();i++)
            {    for(int j=0; j<proj.cols();j++)
          //  cout<<proj(i,j)<<" ";
             if(brdfParam->verOfColorSpace) rgb2labv2(proj(i,0),proj(i,1),proj(i,2),Lab(i,0),Lab(i,1),Lab(i,2));
             else rgb2lab(proj(i,0),proj(i,1),proj(i,2),Lab(i,0),Lab(i,1),Lab(i,2));
             }
          
    return Lab;             
 }
float*  BRDFMeasuredMERL::getRBFActivations(Matrix<float, Dynamic, Dynamic,RowMajor>& Centers, float* betas, float* input)
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
  float  BRDFMeasuredMERL::evaluateFuncApproxRBFN(Matrix<float, Dynamic, Dynamic,RowMajor>& Centers, float* betas, float* Theta, bool normalize, float* input)
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
  

return result;

 
  }


float BRDFMeasuredMERL::normVec(float* x,Matrix<float, Dynamic, Dynamic,RowMajor> &Centers, int row)

{
float result=0.0;
for(int i=0; i<Centers.rows(); i++) 
   result+=(x[i]-Centers(i,row))*(x[i]-Centers(i,row));

return result;
}
void BRDFMeasuredMERL::updateAttr(float* xnew)
{

            
      for(int i=0; i<brdfParam->npzFiles.size(); i++)
      {
    
       cnpy::NpyArray arr_mv1 = brdfParam->npzFiles.at(i)["Centers"];
       double* mv1 = reinterpret_cast<double*>(arr_mv1.data);
       Matrix<float, Dynamic, Dynamic,RowMajor> Centers;
       Centers.resize(arr_mv1.shape[1],arr_mv1.shape[0]);
      for(int i=0; i<arr_mv1.shape[0]*arr_mv1.shape[1];i++)
        Centers.data()[i]=mv1[i];
       
     
      arr_mv1 =  brdfParam->npzFiles.at(i)["betas"];
      double* betasTemp = reinterpret_cast<double*>(arr_mv1.data);
      float* betas = new float[arr_mv1.shape[0]*arr_mv1.shape[1]];
      for(int i=0; i<arr_mv1.shape[0]*arr_mv1.shape[1];i++) betas[i]=betasTemp[i];
      
  
       
         arr_mv1 =  brdfParam->npzFiles.at(i)["Theta"];
    
       double* ThetaTemp = reinterpret_cast<double*>(arr_mv1.data);
      float* Theta = new float[arr_mv1.shape[0]*arr_mv1.shape[1]];
      for(int i=0; i<arr_mv1.shape[0]*arr_mv1.shape[1];i++) Theta[i]=ThetaTemp[i];
    
       
      brdfParam->attrValues[i] = evaluateFuncApproxRBFN(Centers, betas, Theta, true, xnew);
      
     }

}
void  BRDFMeasuredMERL::UnmapBRDF(double* CosineMap,float* mappedData,bool* MaskMap,float* median,int mask_size)

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
         free(mappedData);

}

void  BRDFMeasuredMERL::ProjectToPCSpace(float* data,float* PCs,float* relativeOffset,int Qsize)

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

         MatrixXf project = V*Sinv*Ut*b;
         
         
          MatrixXf proj_Lab= rgb2Lab(project);
          float* xnew = new float[proj_Lab.rows()];
          
       
          for(int i=0; i<proj_Lab.rows();i++)    xnew[i]=proj_Lab(i,0)/100.0;
          
          updateAttr(xnew);
        
          
          for(int i=0; i<project.rows();i++)  
          if(brdfParam->verOfColorSpace)  lab2rgbv2(xnew[i]*100.0,proj_Lab(i,1),proj_Lab(i,2),project(i,0),project(i,1),project(i,2));     
                                                 else  lab2rgb(xnew[i]*100.0,proj_Lab(i,1),proj_Lab(i,2),project(i,0),project(i,1),project(i,2));
                 
           
  
       brdfParam->Q = A;
       brdfParam-> proj = project;
         
         

         MatrixXf recon =  A*project;
            

         for(int i=0; i<Qsize;i++)
             for(int j=0;j<3;j++)
                 data[i*3+j]= recon(i,j)+relativeOffset[i];
                 
                 
            delete xnew;


}
float*  BRDFMeasuredMERL::ProjectToPCSpaceShort()

{

    
         MatrixXf proj=  brdfParam->proj;
         
       MatrixXf proj_Lab= rgb2Lab(proj);
          float* xnew = new float[proj_Lab.rows()];
          
  
           for(int i=0; i<proj_Lab.rows();i++)    xnew[i]=proj_Lab(i,0)/100.0;
           

    
       cnpy::NpyArray arr_mv1 = brdfParam->npzFiles.at(brdfParam->idOfVal)["Centers"];
       double* mv1 = reinterpret_cast<double*>(arr_mv1.data);
       Matrix<float, Dynamic, Dynamic,RowMajor> Centers;
       Centers.resize(arr_mv1.shape[1],arr_mv1.shape[0]);
    for(int i=0; i<arr_mv1.shape[0]*arr_mv1.shape[1];i++)
        Centers.data()[i]=mv1[i];
       
     
      arr_mv1 =  brdfParam->npzFiles.at(brdfParam->idOfVal)["betas"];
      double* betasTemp = reinterpret_cast<double*>(arr_mv1.data);
      float* betas = new float[arr_mv1.shape[0]*arr_mv1.shape[1]];
      for(int i=0; i<arr_mv1.shape[0]*arr_mv1.shape[1];i++) betas[i]=betasTemp[i];

       
         arr_mv1 =  brdfParam->npzFiles.at(brdfParam->idOfVal)["Theta"];
    
       double* ThetaTemp = reinterpret_cast<double*>(arr_mv1.data);
      float* Theta = new float[arr_mv1.shape[0]*arr_mv1.shape[1]];
      for(int i=0; i<arr_mv1.shape[0]*arr_mv1.shape[1];i++) Theta[i]=ThetaTemp[i];
       
      float ynew = evaluateFuncApproxRBFN(Centers, betas, Theta, true, xnew);
      
      float yobj= brdfParam->newAttrVal;
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
  
  
}

     updateAttr(xnew);
     
         for(int i=0; i<proj.rows();i++)   
               if(brdfParam->verOfColorSpace)  lab2rgbv2(xnew[i]*100.0,proj_Lab(i,1),proj_Lab(i,2),proj(i,0),proj(i,1),proj(i,2));     
                                                 else  lab2rgb(xnew[i]*100.0,proj_Lab(i,1),proj_Lab(i,2),proj(i,0),proj(i,1),proj(i,2));
         
      
        
         MatrixXf recon =   brdfParam->Q*proj;
          brdfParam->proj=proj;
   

       float* data = new float[ brdfParam->Qsize * 3 ];
   
          brdfData = new float[ numBRDFSamples * 3 ];

         for(int i=0; i<brdfParam->Qsize;i++)
             for(int j=0;j<3;j++)
                 data[i*3+j]= recon(i,j)+brdfParam->RelativeOffset[i];

return data;
}



void BRDFMeasuredMERL::projectShort(int numBRDFsam,const char *filename)
{
  
         name = std::string(filename);
         numBRDFSamples = numBRDFsam;


          float* reshapedBRDF = ProjectToPCSpaceShort();


          UnmapBRDF(brdfParam->CosineMap,reshapedBRDF,brdfParam->MaskMap,brdfParam->median,brdfParam->mask_size);

           reshapeFinal();

        
}

void BRDFMeasuredMERL::project(const char *filename)
{
  
      name = std::string(filename);
      brdfParam->attrValues = new float[brdfParam->npzFiles.size()];
     
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
          
         brdfParam->MaskMap=MaskMap;
         brdfParam->mask_size=mask_size;
         
        
          my_npy = cnpy::npy_load("Median.npy");
           cout<<"median shape "<<my_npy.shape[0]<<" "<<my_npy.shape[1]<<endl;
          float* median =reinterpret_cast<float*>(my_npy.data);
          Qsize = my_npy.shape[0];
          brdfParam->median=median;
          brdfParam->Qsize=Qsize;
          

          my_npy = cnpy::npy_load("Q.npy");
           cout<<"ScaledEigenvectors shape "<<my_npy.shape[0]<<" "<<my_npy.shape[1]<<endl;

          float* Q  = reinterpret_cast<float*>(my_npy.data);
          
          my_npy = cnpy::npy_load("CosineMap.npy");
           cout<<"CosineMap "<<my_npy.shape[0]<<" "<<my_npy.shape[1]<<endl;
          double* CosineMap  = reinterpret_cast<double*>(my_npy.data);
          
          brdfParam->CosineMap=CosineMap;
          
          float* reshapedBRDF = reshape(CosineMap,MaskMap,mask_size,Qsize,brdfData);

          my_npy = cnpy::npy_load("RelativeOffset.npy");
           cout<<"RelativeOffset "<<my_npy.shape[0]<<" "<<my_npy.shape[1]<<endl;
          float* RelativeOffset  = reinterpret_cast<float*>(my_npy.data);
          
          brdfParam->RelativeOffset=RelativeOffset;

          MapBRDF(reshapedBRDF,median,Qsize);



           ProjectToPCSpace(reshapedBRDF,Q,RelativeOffset,Qsize);


          UnmapBRDF(CosineMap,reshapedBRDF,MaskMap,median,mask_size);
        //  for(int i=0;i<12;i++)
       //  cout<<brdfData[i]<<endl;


           reshapeFinal();

        
}



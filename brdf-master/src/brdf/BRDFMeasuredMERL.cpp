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

using namespace std;

BRDFMeasuredMERL::BRDFMeasuredMERL()
                 : brdfData(NULL)
{
    std::string path = (getShaderTemplatesPath() + "measured.func");

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
    cout<<"init gl"<<endl;
    
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

void  BRDFMeasuredMERL::ProjectToPCSpace(float* data,float* PCs,float* relativeOffset,int Qsize,float* var)

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
         
  
       brdfParam->Q = A;
       brdfParam-> proj = project;
         
         
         for(int i=0; i<project.cols();i++)
         for(int j=0; j<project.rows(); j++)  
           project(j,i)*=var[j];
        
        
         MatrixXf recon =  A*project;
         
         
         
       
        

         for(int i=0; i<Qsize;i++)
             for(int j=0;j<3;j++)
                 data[i*3+j]= recon(i,j)+relativeOffset[i];



}
float*  BRDFMeasuredMERL::ProjectToPCSpaceShort(float* var)

{

    
         MatrixXf proj=  brdfParam->proj;
         
         for(int i=0; i< brdfParam->proj.cols();i++)
         for(int j=0; j< brdfParam->proj.rows(); j++)  
               proj(j,i)*=var[j];
        
        
         MatrixXf recon =   brdfParam->Q*proj;
         
   

       float* data = new float[ brdfParam->Qsize * 3 ];
   
          brdfData = new float[ numBRDFSamples * 3 ];

         for(int i=0; i<brdfParam->Qsize;i++)
             for(int j=0;j<3;j++)
                 data[i*3+j]= recon(i,j)+brdfParam->RelativeOffset[i];

return data;
}



void BRDFMeasuredMERL::projectShort(int numBRDFsam,const char *filename,float* var)
{
  
         name = std::string(filename);
         numBRDFSamples = numBRDFsam;


          float* reshapedBRDF = ProjectToPCSpaceShort(var);


          UnmapBRDF(brdfParam->CosineMap,reshapedBRDF,brdfParam->MaskMap,brdfParam->median,brdfParam->mask_size);

           reshapeFinal();

        
}

void BRDFMeasuredMERL::project(const char *filename,float* var)
{
  
 name = std::string(filename);

        // read brdf
        if (!read_brdf(filename))
        {
            fprintf(stderr, "Error reading %s\n", filename);
            exit(1);
        }
        int mask_size;
        int Qsize;


        brdfParam=new brdfMERLparam;


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



           ProjectToPCSpace(reshapedBRDF,Q,RelativeOffset,Qsize,var);


          UnmapBRDF(CosineMap,reshapedBRDF,MaskMap,median,mask_size);
        //  for(int i=0;i<12;i++)
       //  cout<<brdfData[i]<<endl;


           reshapeFinal();

        
}



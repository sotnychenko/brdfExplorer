// my first program in C++
#include <iostream>
#include <cstdlib>
#include <eigen3/Eigen/Sparse>
#include"cnpy.h"

#include "redsvd.hpp"
#include "redsvdFile.hpp"

#define BRDF_SAMPLING_RES_THETA_H       90
#define BRDF_SAMPLING_RES_THETA_D       90
#define BRDF_SAMPLING_RES_PHI_D         360

using namespace std;
using namespace Eigen;

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



           ProjectToPCSpace(reshapedBRDF,Q,RelativeOffset,Qsize);

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

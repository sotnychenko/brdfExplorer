#ifndef PROJECTTOPCA_H
#define PROJECTTOPCA_H


class projectToPCA 
{


      
public:
        projectToPCA();
      ~projectToPCA();

void ProjectToPCSpace(float* data,float* PCs,float* relativeOffset,int Qsize);

float* project(const char *filename);

void UnmapBRDF(double* CosineMap,float* mappedData,bool* MaskMap,float* median,int mask_size);

void  MapBRDF(float* reshapedBRDF,float* median,int Qsize );

float*  reshape(double* CosineMap,bool* MaskMap,int mask_size,int Qsize, float* my_npy);

bool read_brdf( const char* filename );

 float* reshape( );
 bool loadMERLData( const char* filename );
 
  float* brdfData;

          
};


#endif // PROJECTTOPCA_H

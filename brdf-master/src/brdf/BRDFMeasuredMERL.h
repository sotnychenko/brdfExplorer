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

#ifndef  BRDF_MEASURED_MERL_H
#define BRDF_MEASURED_MERL_H

#include <string>
#include <GL/glew.h>
#include <GL/glut.h>
#include "BRDFBase.h"
#include "projectToPCA.h"
#include <eigen3/Eigen/Sparse>
using namespace Eigen;
using namespace std;
class BRDFMeasuredMERL : public BRDFBase
{
public:

    BRDFMeasuredMERL();
     BRDFMeasuredMERL(float* data, string fileName,int numSamples);
    virtual ~BRDFMeasuredMERL();

    bool loadMERLData( const char* filename );
    bool loadPCAprojection( const char* filename );
    void ProjectToPCSpace(float* data,float* PCs,float* relativeOffset,int Qsize);

void project(const char *filename);
float* ProjectToPCSpaceShort();

void UnmapBRDF(double* CosineMap,float* mappedData,bool* MaskMap,float* median,int mask_size);

void  MapBRDF(float* reshapedBRDF,float* median,int Qsize );
void projectShort(int numBRDFsam,const char *filename);

float*  reshape(double* CosineMap,bool* MaskMap,int mask_size,int Qsize, float* my_npy);

bool read_brdf( const char* filename );
void allocBRDF();
void reshapeFinal( );
void setName(string n);
void updateAttr(float* xnew);
float normVec(float* x,Matrix<float, Dynamic, Dynamic,RowMajor> &Centers, int row);
float evaluateFuncApproxRBFN(Matrix<float, Dynamic, Dynamic,RowMajor>& Centers, float* betas, float* Theta, bool normalize, float* input);
float* getRBFActivations(Matrix<float, Dynamic, Dynamic,RowMajor>& Centers, float* betas, float* input);
MatrixXf  rgb2Lab(MatrixXf proj);
void lab2rgb( float l_s, float a_s, float b_s, float& R, float& G, float& B );
void rgb2lab( float var_R, float var_G, float var_B, float & l_s, float &a_s, float &b_s );
void lab2rgbv2( float l_s, float a_s, float b_s, float& R, float& G, float& B );
void rgb2labv2( float var_R, float var_G, float var_B, float & l_s, float &a_s, float &b_s );
 float sgn(float val);
 
    float* brdfData;
    GLuint tbo;
    GLuint tex;
        int numBRDFSamples;
         brdfMERLparam* brdfParam;
   
       
protected:

    virtual void initGL();
    virtual std::string getBRDFFunction();

    virtual void adjustShaderPreRender( DGLShader* );

private:
    std::string brdfFunction;

    void createTBO();

    // IDs needed for the texture buffer object
 


 
        
};


#endif

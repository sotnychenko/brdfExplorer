#ifndef BRDF_MEASURED_MERL_H
#define BRDF_MEASURED_MERL_H

#include <string>
#include <GL/glew.h>
#include <GL/glut.h>
#include <eigen3/Eigen/Sparse>
#include "BRDFBase.h"
#include "projectToPCA.h"




#include "libqhullcpp/QhullFacetList.h"

using orgQhull::QhullFacetList;

using namespace Eigen;
using namespace std;
class BRDFMeasuredMERL : public BRDFBase {
public:
    BRDFMeasuredMERL();
    BRDFMeasuredMERL(float* data, string fileName, int numSamples);
    void project(const char* filename);
    void projectShort(int numBRDFsam, const char* filename);
    bool loadMERLData(const char* filename);

    brdfMERLparam* brdfParam;
    int numBRDFSamples;

    virtual ~BRDFMeasuredMERL();

    static float evaluateFuncApproxRBFN(Matrix<float, Dynamic, Dynamic, RowMajor>& Centers, float* betas, float* Theta, bool normalize, float* input);
    static MatrixXf rgb2Lab(MatrixXf proj);


protected:
    virtual void initGL();
    virtual std::string getBRDFFunction();

    virtual void adjustShaderPreRender(DGLShader*);

private:
    void ProjectToPCSpace(float* data, float* PCs, float* relativeOffset, int Qsize);

    float* ProjectToPCSpaceShort();

    void UnmapBRDF(double* CosineMap, float* mappedData, bool* MaskMap, float* median, int mask_size);

    void MapBRDF(float* reshapedBRDF, float* median, int Qsize);

    float* reshape(double* CosineMap, bool* MaskMap, int mask_size, int Qsize, float* my_npy);
    bool readBrdf(const char* filename);
    bool inhull(float* x,QhullFacetList& qlist,double tol);
    double dot(float* x,double* a,double* N);
    void reshapeFinal();
    void setName(string n);
    void updateAttr(float* xnew);
    float normVec(float* x, Matrix<float, Dynamic, Dynamic, RowMajor>& Centers, int row);
    static float* getRBFActivations(Matrix<float, Dynamic, Dynamic, RowMajor>& Centers, float* betas, float* input);
    bool inhull(float* xnew, Matrix<double, Dynamic, Dynamic, RowMajor>& hull, float tol );

    static void lab2rgb(float l_s, float a_s, float b_s, float& R, float& G, float& B);
    static void rgb2lab(float var_R, float var_G, float var_B, float& l_s, float& a_s, float& b_s);
    void lab2rgbv2(float l_s, float a_s, float b_s, float& R, float& G, float& B);
    void rgb2labv2(float var_R, float var_G, float var_B, float& l_s, float& a_s, float& b_s);
    static float sgn(float val);

    float* brdfData;
    GLuint tbo;
    GLuint tex;

    std::string brdfFunction;

    void createTBO();
};

#endif

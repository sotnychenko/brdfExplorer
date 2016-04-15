#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <string.h>
#include <limits>
#include <eigen3/Eigen/Sparse>

#include "BRDFMeasuredMERL.h"
#include "DGLShader.h"
#include "Paths.h"
#include "BRDFBase.h"
#include "projectToPCA.h"
#include "cnpy.h"
#include "redsvd.hpp"
#include "redsvdFile.hpp"

#include "libqhullcpp/RboxPoints.h"
#include "libqhullcpp/QhullError.h"
#include "libqhullcpp/QhullQh.h"
#include "libqhullcpp/QhullFacet.h"
#include "libqhullcpp/QhullFacetList.h"
#include "libqhullcpp/QhullLinkedList.h"
#include "libqhullcpp/QhullVertex.h"
#include "libqhullcpp/QhullVertexSet.h"
#include "libqhullcpp/Qhull.h"


using std::cerr;
using std::cin;
using std::cout;
using std::endl;

using orgQhull::Qhull;
using orgQhull::QhullPoint;
using orgQhull::QhullError;
using orgQhull::QhullFacet;
using orgQhull::QhullFacetList;
using orgQhull::QhullQh;
using orgQhull::RboxPoints;
using orgQhull::QhullVertex;
using orgQhull::QhullVertexSet;
using orgQhull::QhullLinkedList;

using namespace Eigen;

#define BRDF_SAMPLING_RES_THETA_H 90
#define BRDF_SAMPLING_RES_THETA_D 90
#define BRDF_SAMPLING_RES_PHI_D 360

#define EPS 0.01

bool custom_isnan(float var)
{
    volatile float d = var;
    return d != d;
}

float BRDFMeasuredMERL::sgn(float val)
{
    return (val < 0.0) ? -1.0 : 1.0;
}

BRDFMeasuredMERL::BRDFMeasuredMERL()
    : brdfData(NULL)
{
    std::string path = (getShaderTemplatesPath() + "measured.func");
    brdfParam = new brdfMERLparam;
    // read the shader
    std::ifstream ifs(path.c_str());
    std::string temp;
    while (getline(ifs, temp))
        brdfFunction += (temp + "\n");
}

BRDFMeasuredMERL::BRDFMeasuredMERL(float* data, string fileName, int numSamples)

{
    name = fileName;
    numBRDFSamples = numSamples;
    brdfData = data;
    std::string path = (getShaderTemplatesPath() + "measured.func");

    // read the shader
    std::ifstream ifs(path.c_str());
    std::string temp;
    while (getline(ifs, temp))
        brdfFunction += (temp + "\n");
}
void BRDFMeasuredMERL::setName(string n)
{
    name = n;
}

BRDFMeasuredMERL::~BRDFMeasuredMERL()
{
    glBindBuffer(GL_TEXTURE_BUFFER_EXT, tbo);
    glDeleteBuffers(1, &tbo);
}

std::string BRDFMeasuredMERL::getBRDFFunction()
{
    return brdfFunction;
}

bool BRDFMeasuredMERL::loadMERLData(const char* filename)
{
    // the BRDF's name is just the filename
    name = std::string(filename);

    // read in the MERL BRDF data
    FILE* f = fopen(filename, "rb");
    if (!f)
        return false;

    int dims[3];
    if (fread(dims, sizeof(int), 3, f) != 3) {
        fprintf(stderr, "read error\n");
        fclose(f);
        return false;
    }
    numBRDFSamples = dims[0] * dims[1] * dims[2];
    if (numBRDFSamples != BRDF_SAMPLING_RES_THETA_H * BRDF_SAMPLING_RES_THETA_D * BRDF_SAMPLING_RES_PHI_D / 2) {
        fprintf(stderr, "Dimensions don't match\n");
        fclose(f);
        return false;
    }

    // read the data
    double* brdf = (double*)malloc(sizeof(double) * 3 * numBRDFSamples);
    if (fread(brdf, sizeof(double), 3 * numBRDFSamples, f) != size_t(3 * numBRDFSamples)) {
        fprintf(stderr, "read error\n");
        fclose(f);
        return false;
    }
    fclose(f);

    // now transform it to RGBA floats
    brdfData = new float[numBRDFSamples * 3];
    for (int i = 0; i < numBRDFSamples; i++) {
        brdfData[i * 3 + 0] = brdf[i * 3 + 0];
        brdfData[i * 3 + 1] = brdf[i * 3 + 1];
        brdfData[i * 3 + 2] = brdf[i * 3 + 2];
    }

    // now we can dump the old data
    free(brdf);

    return true;
}

void BRDFMeasuredMERL::initGL()
{
    if (initializedGL)
        return;

    // create buffer object
    glGenBuffers(1, &tbo);
    glBindBuffer(GL_TEXTURE_BUFFER_EXT, tbo);

    // initialize buffer object
    unsigned int numBytes = numBRDFSamples * 3 * sizeof(float);
    //printf( "size = %d bytes (%f megs)\n", numBytes, float(numBytes) / 1048576.0f );
    glBufferData(GL_TEXTURE_BUFFER_EXT, numBytes, 0, GL_DYNAMIC_DRAW);

    //tex
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_BUFFER_EXT, tex);
    glTexBufferEXT(GL_TEXTURE_BUFFER_EXT, GL_INTENSITY32F_ARB, tbo);
    glBindBuffer(GL_TEXTURE_BUFFER_EXT, 0);

    glBindBuffer(GL_TEXTURE_BUFFER_EXT, tbo);
    float* p = (float*)glMapBuffer(GL_TEXTURE_BUFFER_EXT, GL_WRITE_ONLY);

    //
    memcpy(p, brdfData, numBytes);
    glUnmapBuffer(GL_TEXTURE_BUFFER_EXT);
    glBindBuffer(GL_TEXTURE_BUFFER_EXT, 0);

    delete[] brdfData;
    brdfData = NULL;

    initializedGL = true;
}

void BRDFMeasuredMERL::adjustShaderPreRender(DGLShader* shader)
{
    shader->setUniformTexture("measuredData", tex, GL_TEXTURE_BUFFER_EXT);

    BRDFBase::adjustShaderPreRender(shader);
}

void BRDFMeasuredMERL::reshapeFinal()
{

    int dims[3];
    dims[0] = BRDF_SAMPLING_RES_THETA_H;
    dims[1] = BRDF_SAMPLING_RES_THETA_D;
    dims[2] = BRDF_SAMPLING_RES_PHI_D / 2;
    int numBRDFSamples = dims[0] * dims[1] * dims[2];

    float* brdf = new float[numBRDFSamples * 3];

    for (int i = 0; i < numBRDFSamples; i++) {

        brdf[i * 3 + 0] = (brdfData[i * 3 + 0]) / (1.00 / 1500.0);
        brdf[i * 3 + 1] = (brdfData[i * 3 + 1]) / (1.15 / 1500.0);
        brdf[i * 3 + 2] = (brdfData[i * 3 + 2]) / (1.66 / 1500.0);
    }
    // cout<< brdf[0]<<endl<<brdf[1]<<endl<<brdf[2]<<endl<<brdf[3]<<endl<<brdf[4]<<endl<<brdf[5];
    float(&reshaped)[180][90][90][3] = *reinterpret_cast<float(*)[BRDF_SAMPLING_RES_PHI_D / 2][BRDF_SAMPLING_RES_THETA_D][BRDF_SAMPLING_RES_THETA_H][3]>(brdf);

    int ind = 0;

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < BRDF_SAMPLING_RES_THETA_H; j++)
            for (int k = 0; k < BRDF_SAMPLING_RES_THETA_D; k++)
                for (int l = 0; l < BRDF_SAMPLING_RES_PHI_D / 2; l++)
                    brdfData[ind++] = reshaped[l][j][k][i];

    free(brdf);
}

bool BRDFMeasuredMERL::readBrdf(const char* filename)
{

    // read in the MERL BRDF data
    FILE* f = fopen(filename, "rb");
    if (!f)
        return false;

    int dims[3];
    if (fread(dims, sizeof(int), 3, f) != 3) {
        fprintf(stderr, "read error\n");
        fclose(f);
        return false;
    }
    numBRDFSamples = dims[0] * dims[1] * dims[2];
    if (numBRDFSamples != BRDF_SAMPLING_RES_THETA_H * BRDF_SAMPLING_RES_THETA_D * BRDF_SAMPLING_RES_PHI_D / 2) {
        fprintf(stderr, "Dimensions don't match\n");
        fclose(f);
        return false;
    }

    // read the data
    double* brdf = (double*)malloc(sizeof(double) * 3 * numBRDFSamples);
    fread(brdf, sizeof(double), 3 * numBRDFSamples, f);

    fclose(f);

    // now transform it to RGBA floats

    brdfData = new float[numBRDFSamples * 3];
    double(&reshaped)[3][90][90][180] = *reinterpret_cast<double(*)[3][90][90][180]>(brdf);
    int ind;
    for (int i = 0; i < 3; i++) {
        ind = 0;
        for (int j = 0; j < 180; j++)
            for (int k = 0; k < 90; k++)
                for (int l = 0; l < 90; l++) {

                    brdfData[ind * 3 + i] = reshaped[i][k][l][j];
                    ind++;
                }
    }
    for (int i = 0; i < numBRDFSamples; i++) {
        brdfData[i * 3 + 0] *= 1.00 / 1500;
        brdfData[i * 3 + 1] *= 1.15 / 1500;
        brdfData[i * 3 + 2] *= 1.66 / 1500;
    }

    free(brdf);

    return true;
}
float* BRDFMeasuredMERL::reshape(double* CosineMap, bool* MaskMap, int mask_size, int Qsize, float* my_npy)
{
    float* result = new float[Qsize * 3];
    int ind = 0;
    for (int i = 0; i < mask_size; i++)
        if (MaskMap[i]) {
            result[ind * 3 + 0] = my_npy[i * 3 + 0] * CosineMap[ind];
            result[ind * 3 + 1] = my_npy[i * 3 + 1] * CosineMap[ind];
            result[ind * 3 + 2] = my_npy[i * 3 + 2] * CosineMap[ind];
            ind++;
        }

    return result;
}

void BRDFMeasuredMERL::MapBRDF(float* reshapedBRDF, float* median, int Qsize)
{
    for (int i = 0; i < Qsize; i++) {
        reshapedBRDF[i * 3 + 0] = log((reshapedBRDF[i * 3 + 0] + 0.001) / (median[i] + 0.001));
        reshapedBRDF[i * 3 + 1] = log((reshapedBRDF[i * 3 + 1] + 0.001) / (median[i] + 0.001));
        reshapedBRDF[i * 3 + 2] = log((reshapedBRDF[i * 3 + 2] + 0.001) / (median[i] + 0.001));
    }
}
void BRDFMeasuredMERL::rgb2labv2(float R, float G, float B, float& l_s, float& a_s, float& b_s)
{

    float var_R = R;
    float var_G = G;
    float var_B = B;

    //Observer. = 2°, Illuminant = D65
    float X = var_R * 0.4124 + var_G * 0.3576 + var_B * 0.1805;
    float Y = var_R * 0.2126 + var_G * 0.7152 + var_B * 0.0722;
    float Z = var_R * 0.0193 + var_G * 0.1192 + var_B * 0.9505;

    float var_X = X / 0.95047; //ref_X =  95.047   Observer= 2°, Illuminant= D65
    float var_Y = Y; //ref_Y = 100.000
    float var_Z = Z / 1.08883; //ref_Z = 108.883

    if (var_X > 0.008856)
        var_X = pow(var_X, (1.f / 3.f));
    else
        var_X = (7.787 * var_X) + (16. / 116.);
    if (var_Y > 0.008856)
        var_Y = pow(var_Y, (1.f / 3.f));
    else
        var_Y = (7.787 * var_Y) + (16. / 116.);
    if (var_Z > 0.008856)
        var_Z = pow(var_Z, (1.f / 3.f));
    else
        var_Z = (7.787 * var_Z) + (16. / 116.);

    l_s = (116. * var_Y) - 16.;
    a_s = 500. * (var_X - var_Y);
    b_s = 200. * (var_Y - var_Z);
}

//http://www.easyrgb.com/index.php?X=MATH&H=01#text1
void BRDFMeasuredMERL::lab2rgbv2(float l_s, float a_s, float b_s, float& R, float& G, float& B)
{
    float var_Y = (l_s + 16.) / 116.;
    float var_X = a_s / 500. + var_Y;
    float var_Z = var_Y - b_s / 200.;

    if (pow(var_Y, 3) > 0.008856)
        var_Y = pow(var_Y, 3);
    else
        var_Y = (var_Y - 16. / 116.) / 7.787;
    if (pow(var_X, 3) > 0.008856)
        var_X = pow(var_X, 3);
    else
        var_X = (var_X - 16. / 116.) / 7.787;
    if (pow(var_Z, 3) > 0.008856)
        var_Z = pow(var_Z, 3);
    else
        var_Z = (var_Z - 16. / 116.) / 7.787;

    float X = 95.047 * var_X; //ref_X =  95.047     Observer= 2°, Illuminant= D65
    float Y = 100.000 * var_Y; //ref_Y = 100.000
    float Z = 108.883 * var_Z; //ref_Z = 108.883

    var_X = X / 100.; //X from 0 to  95.047      (Observer = 2°, Illuminant = D65)
    var_Y = Y / 100.; //Y from 0 to 100.000
    var_Z = Z / 100.; //Z from 0 to 108.883

    R = var_X * 3.2406 + var_Y * -1.5372 + var_Z * -0.4986;
    G = var_X * -0.9689 + var_Y * 1.8758 + var_Z * 0.0415;
    B = var_X * 0.0557 + var_Y * -0.2040 + var_Z * 1.0570;
}
void BRDFMeasuredMERL::rgb2lab(float var_R, float var_G, float var_B, float& l_s, float& a_s, float& b_s)
{

    //Observer. = 2°, Illuminant = D65
    float X = var_R * 0.412453 + var_G * 0.357580 + var_B * 0.180423;
    float Y = var_R * 0.212671 + var_G * 0.715160 + var_B * 0.072169;
    float Z = var_R * 0.019334 + var_G * 0.119193 + var_B * 0.950227;

    float var_X = X / 0.95046; //ref_X =  95.047   Observer= 2°, Illuminant= D65
    float var_Y = Y; //ref_Y = 100.000
    float var_Z = Z / 1.088754; //ref_Z = 108.883
    float var_Yf = var_Y;
    if (var_X > 0.008856)
        var_X = sgn(var_X) * pow(abs(var_X), (1.f / 3.f));
    else
        var_X = (7.787 * var_X) + (16.f / 116.f);
    if (var_Y > 0.008856)
        var_Y = sgn(var_Y) * pow(abs(var_Y), (1.f / 3.f));
    else
        var_Y = (7.787 * var_Y) + (16.f / 116.f);
    if (var_Z > 0.008856)
        var_Z = sgn(var_Z) * pow(abs(var_Z), (1.f / 3.f));
    else
        var_Z = (7.787 * var_Z) + (16.f / 116.f);

    if (var_Yf > 0.008856)
        l_s = (116.0 * var_Y) - 16.0;
    else
        l_s = (903.3 * var_Yf);
    a_s = 500. * (var_X - var_Y);
    b_s = 200. * (var_Y - var_Z);
}

//http://www.easyrgb.com/index.php?X=MATH&H=01#text1
void BRDFMeasuredMERL::lab2rgb(float l_s, float a_s, float b_s, float& R, float& G, float& B)
{
    float T1 = 0.008856;
    float T2 = 0.206893;

    float fY = pow((l_s + 16.0) / 116.0, 1.0 / 3.0);

    float fY_temp = fY;
    float Y;
    if (fY > T1)
        Y = pow((l_s + 16.0) / 116.0, 3.0);
    else
        Y = (l_s / 903.3);

    if (fY_temp > T1)
        fY = sgn(Y) * pow(abs(Y), (1.f / 3.f));
    else
        fY = (7.787 * Y) + (16. / 116.);

    float fX = a_s / 500.0 + fY;
    float X;

    if (fX > T2)
        X = pow(fX, 3.0f);
    else
        X = (fX - 16.0 / 116.0) / 7.787;

    float fZ = fY - b_s / 200.0;
    float Z;
    if (fX > T2)
        Z = pow(fZ, 3.0f);
    else
        Z = (fZ - 16.0 / 116.0) / 7.787;

    X = X * 0.950456;
    Z = Z * 1.088754;

    R = X * 3.240479 + Y * -1.537150 + Z * -0.498535;
    G = X * -0.969256 + Y * 1.875992 + Z * 0.041556;
    B = X * 0.055648 + Y * -0.204043 + Z * 1.057311;
}

MatrixXf BRDFMeasuredMERL::rgb2Lab(MatrixXf proj)
{
    MatrixXf Lab;
    Lab.resize(proj.rows(), proj.cols());

    for (int i = 0; i < proj.rows(); i++) {
        for (int j = 0; j < proj.cols(); j++)
            //  cout<<proj(i,j)<<" ";
           // if (brdfParam->verOfColorSpace)
                rgb2labv2(proj(i, 0), proj(i, 1), proj(i, 2), Lab(i, 0), Lab(i, 1), Lab(i, 2));
          //  else
          //      rgb2lab(proj(i, 0), proj(i, 1), proj(i, 2), Lab(i, 0), Lab(i, 1), Lab(i, 2));
    }

    return Lab;
}
float* BRDFMeasuredMERL::getRBFActivations(Matrix<float, Dynamic, Dynamic, RowMajor>& Centers, float* betas, float* input)
{

    float* sqrdDists = new float[Centers.cols()];
    for (int i = 0; i < Centers.cols(); i++) {
        sqrdDists[i] = 0.0;
        for (int j = 0; j < Centers.rows(); j++)

            sqrdDists[i] += (Centers(j, i) - input[j]) * (Centers(j, i) - input[j]);

        //  cout<<sqrdDists[i]<<" ";

        //cout<<endl;
    }

    for (int i = 0; i < Centers.cols(); i++)
        sqrdDists[i] = exp(-1.0 * betas[i] * sqrdDists[i]);

    return sqrdDists;
}
float BRDFMeasuredMERL::evaluateFuncApproxRBFN(Matrix<float, Dynamic, Dynamic, RowMajor>& Centers, float* betas, float* Theta, bool normalize, float* input)
{

    float result = 0.0;
    float* phis = getRBFActivations(Centers, betas, input);

    if (normalize) {

        float sum = 0.0;
        for (int i = 0; i < Centers.cols(); i++)
            sum += phis[i];
        for (int i = 0; i < Centers.cols(); i++)
            phis[i] /= sum;
    }

    result += Theta[0];
    for (int i = 0; i < Centers.cols(); i++)
        result += phis[i] * Theta[i + 1];

    return result;
}

float BRDFMeasuredMERL::normVec(float* x, Matrix<float, Dynamic, Dynamic, RowMajor>& Centers, int row)

{
    float result = 0.0;
    for (int i = 0; i < Centers.rows(); i++)
        result += (x[i] - Centers(i, row)) * (x[i] - Centers(i, row));

    return result;
}
void BRDFMeasuredMERL::updateAttr(float* xnew)
{

    for (int i = 0; i < brdfParam->npzFiles.size(); i++) {

        cnpy::NpyArray arr_mv1 = brdfParam->npzFiles.at(i)["Centers"];
        double* mv1 = reinterpret_cast<double*>(arr_mv1.data);
        Matrix<float, Dynamic, Dynamic, RowMajor> Centers;
        Centers.resize(arr_mv1.shape[1], arr_mv1.shape[0]);
        for (int i = 0; i < arr_mv1.shape[0] * arr_mv1.shape[1]; i++)
            Centers.data()[i] = mv1[i];

        arr_mv1 = brdfParam->npzFiles.at(i)["betas"];
        double* betasTemp = reinterpret_cast<double*>(arr_mv1.data);
        float* betas = new float[arr_mv1.shape[0] * arr_mv1.shape[1]];
        for (int i = 0; i < arr_mv1.shape[0] * arr_mv1.shape[1]; i++)
            betas[i] = betasTemp[i];

        arr_mv1 = brdfParam->npzFiles.at(i)["Theta"];

        double* ThetaTemp = reinterpret_cast<double*>(arr_mv1.data);
        float* Theta = new float[arr_mv1.shape[0] * arr_mv1.shape[1]];
        for (int i = 0; i < arr_mv1.shape[0] * arr_mv1.shape[1]; i++)
            Theta[i] = ThetaTemp[i];

        brdfParam->attrValues[i] = evaluateFuncApproxRBFN(Centers, betas, Theta, true, xnew);
    }
}
void BRDFMeasuredMERL::UnmapBRDF(double* CosineMap, float* mappedData, bool* MaskMap, float* median, int mask_size)

{
    int ind = 0;
    for (int i = 0; i < mask_size; i++)
        if (MaskMap[i]) {
            brdfData[i * 3 + 0] = exp(mappedData[ind * 3 + 0]) * (median[ind] + 0.001) - 0.001;
            brdfData[i * 3 + 1] = exp(mappedData[ind * 3 + 1]) * (median[ind] + 0.001) - 0.001;
            brdfData[i * 3 + 2] = exp(mappedData[ind * 3 + 2]) * (median[ind] + 0.001) - 0.001;

            brdfData[i * 3 + 0] /= CosineMap[ind];
            brdfData[i * 3 + 1] /= CosineMap[ind];
            brdfData[i * 3 + 2] /= CosineMap[ind];

            ind++;
        }
        else
            brdfData[i * 3 + 0] = brdfData[i * 3 + 1] = brdfData[i * 3 + 2] = -1.0;

    free(mappedData);
}

void BRDFMeasuredMERL::ProjectToPCSpace(float* data, float* PCs, float* relativeOffset, int Qsize)

{
    int numC = 5;
    Matrix<float, Dynamic, Dynamic, RowMajor> b;
    b.resize(Qsize, 3);
    for (int i = 0; i < Qsize; i++) {
        b.data()[3 * i + 0] = data[3 * i + 0] - relativeOffset[i];
        b.data()[3 * i + 1] = data[3 * i + 1] - relativeOffset[i];
        b.data()[3 * i + 2] = data[3 * i + 2] - relativeOffset[i];
    }

    //cout<<b(0,0)<<endl<<b(0,1)<<endl<<b(0,2)<<endl;

    Matrix<float, Dynamic, Dynamic, RowMajor> A;
    A.resize(Qsize, numC);
    for (int i = 0; i < Qsize; i++)
        for (int j = 0; j < numC; j++)
            A.data()[numC * i + j] = PCs[numC * i + j];

    //  cout<<A(0,0)<<endl<<A(0,1)<<endl<<A(0,2)<<endl;

    REDSVD::RedSVD SVD(A, numC);

    MatrixXf Ut = SVD.matrixU().transpose();
    MatrixXf V = SVD.matrixV().transpose();

    MatrixXf Sinv;
    Sinv.resize(numC, numC);
    Sinv.setZero();
    for (int i = 0; i < numC; i++)
        Sinv(i, i) = SVD.singularValues()[i] / (SVD.singularValues()[i] * SVD.singularValues()[i]);

    MatrixXf project = V * Sinv * Ut * b;

    MatrixXf proj_Lab = rgb2Lab(project);
    float* xnew = new float[proj_Lab.rows()];

    for (int i = 0; i < proj_Lab.rows(); i++)
        xnew[i] = proj_Lab(i, 0) / 100.0;

    updateAttr(xnew);

    for (int i = 0; i < project.rows(); i++)
    //    if (brdfParam->verOfColorSpace)
           lab2rgbv2(xnew[i] * 100.0, proj_Lab(i, 1), proj_Lab(i, 2), project(i, 0), project(i, 1), project(i, 2));
     //   else
        //    lab2rgb(xnew[i] * 100.0, proj_Lab(i, 1), proj_Lab(i, 2), project(i, 0), project(i, 1), project(i, 2));

    brdfParam->Q = A;
    brdfParam->proj = project;

    MatrixXf recon = A * project;

    for (int i = 0; i < Qsize; i++)
        for (int j = 0; j < 3; j++)
            data[i * 3 + j] = recon(i, j) + relativeOffset[i];

    delete [] xnew;
}

double BRDFMeasuredMERL::dot(float* x, double* a, double* N)
{
    double result = 0.0;
    for(int i=0;i<5; i++)
        result+=(x[i]-a[i])*(-1.0)*N[i];

    return result;
}
float BRDFMeasuredMERL::evaluateBarLog(float* x,QhullFacetList& qlist,double tol)
{
    float result = tol;
    for( QhullLinkedList<QhullFacet>::iterator facet = qlist.begin(); facet!=qlist.end(); facet++ )
        result+= log(dot(x,facet->vertices().at(1).point().coordinates(),facet->hyperplane().coordinates()));


    return result;

}


float BRDFMeasuredMERL::evaluateBarDer(float* x,QhullFacetList& qlist, int grad)
{

    float der = 0.0;
    for( QhullLinkedList<QhullFacet>::iterator facet = qlist.begin(); facet!=qlist.end(); facet++ )
     {

        der+=(-1.0)*facet->hyperplane().coordinates()[grad]/(dot(x,facet->vertices().at(1).point().coordinates(),facet->hyperplane().coordinates())+EPS);
      }

    return der;

}

bool BRDFMeasuredMERL::inhull(float* x,QhullFacetList& qlist,double tol)
{
     bool inside = true;

    for( QhullLinkedList<QhullFacet>::iterator facet = qlist.begin(); facet!=qlist.end(); facet++ )
    {

       if(dot(x,facet->vertices().at(1).point().coordinates(),facet->hyperplane().coordinates())<-1.0*tol)
          { inside = false; break;}

    }

    return inside;
}

float* BRDFMeasuredMERL::ProjectToPCSpaceShort()

{

    MatrixXf proj = brdfParam->proj;
    MatrixXf proj_Lab = rgb2Lab(proj);
    float* xnew = new float[proj_Lab.rows()];

    for (int i = 0; i < proj_Lab.rows(); i++)
        xnew[i] = proj_Lab(i, 0) / 100.0;

    cnpy::NpyArray arr_mv1 = brdfParam->npzFiles.at(brdfParam->idOfVal)["Centers"];
    double* mv1 = reinterpret_cast<double*>(arr_mv1.data);
    Matrix<float, Dynamic, Dynamic, RowMajor> Centers;
    Centers.resize(arr_mv1.shape[1], arr_mv1.shape[0]);
    for (int i = 0; i < arr_mv1.shape[0] * arr_mv1.shape[1]; i++)
        Centers.data()[i] = mv1[i];




    arr_mv1 = brdfParam->npzFiles.at(brdfParam->idOfVal)["betas"];
    double* betasTemp = reinterpret_cast<double*>(arr_mv1.data);
    float* betas = new float[arr_mv1.shape[0] * arr_mv1.shape[1]];
    for (int i = 0; i < arr_mv1.shape[0] * arr_mv1.shape[1]; i++)
        betas[i] = betasTemp[i];

    arr_mv1 = brdfParam->npzFiles.at(brdfParam->idOfVal)["Theta"];

    double* ThetaTemp = reinterpret_cast<double*>(arr_mv1.data);
    float* Theta = new float[arr_mv1.shape[0] * arr_mv1.shape[1]];
    for (int i = 0; i < arr_mv1.shape[0] * arr_mv1.shape[1]; i++)
        Theta[i] = ThetaTemp[i];






    float ynew = evaluateFuncApproxRBFN(Centers, betas, Theta, true, xnew);

    if(brdfParam->paths.at(brdfParam->idOfVal).alpha.size()==0) {
        brdfParam->paths.at(brdfParam->idOfVal).alpha.resize(5);
        brdfParam->paths.at(brdfParam->idOfVal).alpha.at(0).push_back(xnew[0]);
        brdfParam->paths.at(brdfParam->idOfVal).alpha.at(1).push_back(xnew[1]);
        brdfParam->paths.at(brdfParam->idOfVal).alpha.at(2).push_back(xnew[2]);
        brdfParam->paths.at(brdfParam->idOfVal).alpha.at(3).push_back(xnew[3]);
        brdfParam->paths.at(brdfParam->idOfVal).alpha.at(4).push_back(xnew[4]);
    }

    float yobj = brdfParam->newAttrVal;
    float alpha = 0.01;

    cnpy::npz_t my_npz = cnpy::npz_load("hull.npz");
    arr_mv1 = my_npz["P"];
    double* PData = reinterpret_cast<double*>(arr_mv1.data);

 //   Matrix<double, Dynamic, Dynamic, RowMajor> hull;
  //  hull.resize(arr_mv1.shape[1], arr_mv1.shape[0]);
  //  for (int i = 0; i < arr_mv1.shape[0] * arr_mv1.shape[1]; i++)
  //      hull.data()[i] = PData[i];


   // double tol =1.e-12*(hull.cwiseAbs().mean());

    Qhull qhull;
    RboxPoints rbox;
    rbox.setDimension(5);
    rbox.append(arr_mv1.shape[1]*arr_mv1.shape[0],PData);

    qhull.runQhull(rbox, "");
    QhullFacetList qlist =qhull.facetList();


     float tol =0.001;
     float mu=150.0f;

   cout<<"start"<<endl;
  // cout<<xnew[0]<<";"<<xnew[1]<<";"<<xnew[2]<<";"<<xnew[3]<<";"<<xnew[4]<<";"<<endl;
if(abs(ynew - yobj) > EPS)
{
   if(!inhull(xnew,qlist,tol)){
       cout<<"was not in hull, changing to random"<<endl;
         //brdfParam->wasNotInHull = false;
    while(!inhull(xnew,qlist,tol))
    {

        xnew[0]= static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        xnew[1]= static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        xnew[2]= static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        xnew[3]= static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        xnew[4]= static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
        ynew =  evaluateFuncApproxRBFN(Centers, betas, Theta, true, xnew);

    }

   }

    while (abs(ynew - yobj) > EPS || mu>EPS) {


        float bar=evaluateBarLog (xnew,qlist,tol);


       /*if(!custom_isnan(bar))
       {

        brdfParam->xold[0]=xnew[0];
        brdfParam->xold[1]=xnew[1];
        brdfParam->xold[2]=xnew[2];
        brdfParam->xold[3]=xnew[3];
        brdfParam->xold[4]=xnew[4];
       }*/


        float* grad = new float[Centers.rows()];
        for (int i = 0; i < Centers.rows(); i++) {
            grad[i] = 0.0;
            for (int j = 0; j < Centers.cols(); j++) {
                grad[i] += -2.0 * betas[j] * Theta[j + 1] * (xnew[i] - Centers(i, j)) * exp(-1.0 * betas[j] * normVec(xnew, Centers, j)) - mu*(evaluateBarDer(xnew,qlist,i));
            }
        }

        float norm = 0.0;
        for (int i = 0; i < Centers.rows(); i++)
            norm += grad[i] * grad[i];

        norm = sqrt(norm);

        for (int i = 0; i < Centers.rows(); i++) {
            grad[i] /= norm;
            xnew[i] -= sgn(ynew - yobj) * alpha * grad[i];
        }


        ynew =  evaluateFuncApproxRBFN(Centers, betas, Theta, true, xnew)-mu*bar;
        brdfParam->paths.at(brdfParam->idOfVal).alpha.at(0).push_back(xnew[0]);
        brdfParam->paths.at(brdfParam->idOfVal).alpha.at(1).push_back(xnew[1]);
        brdfParam->paths.at(brdfParam->idOfVal).alpha.at(2).push_back(xnew[2]);
        brdfParam->paths.at(brdfParam->idOfVal).alpha.at(3).push_back(xnew[3]);
        brdfParam->paths.at(brdfParam->idOfVal).alpha.at(4).push_back(xnew[4]);

        mu*=0.46f;
        cout<<"ln(c(x))="<<bar<<endl;

       // if(!inhull(xnew,qlist,tol)){brdfParam->newAttrVal = ynew; cout<<"not in hull"<<endl; break; }

    }
     /*if(!inhull(xnew,qlist,tol)){
         cout<<"not in hull"<<endl;

          brdfParam->wasNotInHull = true;
          brdfParam->xold[0]=  xnew[0];
          brdfParam->xold[1]=  xnew[1];
          brdfParam->xold[2]= xnew[2];
          brdfParam->xold[3]= xnew[3];
          brdfParam->xold[4]= xnew[4];
          int pathSize = brdfParam->paths.at(brdfParam->idOfVal).alpha.at(0).size();
          xnew[0]=  brdfParam->paths.at(brdfParam->idOfVal).alpha.at(0).at(pathSize-3);
         xnew[1]= brdfParam->paths.at(brdfParam->idOfVal).alpha.at(0).at(pathSize-3);
          xnew[2]= brdfParam->paths.at(brdfParam->idOfVal).alpha.at(0).at(pathSize-3);
          xnew[3]= brdfParam->paths.at(brdfParam->idOfVal).alpha.at(0).at(pathSize-3);
          xnew[4]= brdfParam->paths.at(brdfParam->idOfVal).alpha.at(0).at(pathSize-3);


     }*/

    // if(inhull(xnew,qlist,tol)) cout<<"now in hull!"<<endl;

     cout<<"finished"<<endl;
   //  cout<<xnew[0]<<";"<<xnew[1]<<";"<<xnew[2]<<";"<<xnew[3]<<";"<<xnew[4]<<";"<<endl;

     ynew =  evaluateFuncApproxRBFN(Centers, betas, Theta, true, xnew);

    brdfParam->newAttrVal = ynew;

}

    updateAttr(xnew);

   for (int i = 0; i < proj.rows(); i++)
    //    if (brdfParam->verOfColorSpace)
            lab2rgbv2(xnew[i] * 100.0, proj_Lab(i, 1), proj_Lab(i, 2), proj(i, 0), proj(i, 1), proj(i, 2));
   //     else
       //     lab2rgb(xnew[i] * 100.0, proj_Lab(i, 1), proj_Lab(i, 2), proj(i, 0), proj(i, 1), proj(i, 2));

    MatrixXf recon = brdfParam->Q * proj;
    brdfParam->proj = proj;

    float* data = new float[brdfParam->Qsize * 3];

    brdfData = new float[numBRDFSamples * 3];

    for (int i = 0; i < brdfParam->Qsize; i++)
        for (int j = 0; j < 3; j++)
            data[i * 3 + j] = recon(i, j) + brdfParam->RelativeOffset[i];



    delete [] xnew;

    return data;
}

void BRDFMeasuredMERL::projectShort(int numBRDFsam, const char* filename)
{

    name = std::string(filename);
    numBRDFSamples = numBRDFsam;

    float* reshapedBRDF = ProjectToPCSpaceShort();

    UnmapBRDF(brdfParam->CosineMap, reshapedBRDF, brdfParam->MaskMap, brdfParam->median, brdfParam->mask_size);

    reshapeFinal();
}

void BRDFMeasuredMERL::project(const char* filename)
{

    name = std::string(filename);
    brdfParam->attrValues = new float[brdfParam->npzFiles.size()];
    brdfParam->paths.resize(brdfParam->npzFiles.size());
    brdfParam->idOfVal=0;
    brdfParam->wasNotInHull = false;
    brdfParam->xold = new float[5];

    // read brdf
    if (!readBrdf(filename)) {
        fprintf(stderr, "Error reading %s\n", filename);
        exit(1);
    }
    int mask_size;
    int Qsize;

    cnpy::NpyArray my_npy = cnpy::npy_load("MaskMap.npy");

    cout << "mask shape " << my_npy.shape[0] << " " << my_npy.shape[1] << endl;
    mask_size = my_npy.shape[0];
    bool* MaskMap = reinterpret_cast<bool*>(my_npy.data);

    brdfParam->MaskMap = MaskMap;
    brdfParam->mask_size = mask_size;

    my_npy = cnpy::npy_load("Median.npy");
    cout << "median shape " << my_npy.shape[0] << " " << my_npy.shape[1] << endl;
    float* median = reinterpret_cast<float*>(my_npy.data);
    Qsize = my_npy.shape[0];
    brdfParam->median = median;
    brdfParam->Qsize = Qsize;

    my_npy = cnpy::npy_load("Q.npy");
    cout << "ScaledEigenvectors shape " << my_npy.shape[0] << " " << my_npy.shape[1] << endl;

    float* Q = reinterpret_cast<float*>(my_npy.data);

    my_npy = cnpy::npy_load("CosineMap.npy");
    cout << "CosineMap " << my_npy.shape[0] << " " << my_npy.shape[1] << endl;
    double* CosineMap = reinterpret_cast<double*>(my_npy.data);

    brdfParam->CosineMap = CosineMap;

    float* reshapedBRDF = reshape(CosineMap, MaskMap, mask_size, Qsize, brdfData);

    my_npy = cnpy::npy_load("RelativeOffset.npy");
    cout << "RelativeOffset " << my_npy.shape[0] << " " << my_npy.shape[1] << endl;
    float* RelativeOffset = reinterpret_cast<float*>(my_npy.data);

    brdfParam->RelativeOffset = RelativeOffset;

    MapBRDF(reshapedBRDF, median, Qsize);

    ProjectToPCSpace(reshapedBRDF, Q, RelativeOffset, Qsize);

    UnmapBRDF(CosineMap, reshapedBRDF, MaskMap, median, mask_size);

    reshapeFinal();
}

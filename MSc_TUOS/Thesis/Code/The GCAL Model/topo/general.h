#include <morph/Visual.h>

#ifdef __OSX__
#include <OpenGL/gl3.h>
#else
#include <GL3/gl3.h>
#endif

#include <morph/tools.h>
#include <morph/VisualDataModel.h>
#include <morph/ColourMap.h>
#include <morph/HexGrid.h>
#include <morph/MathAlgo.h>
#include <morph/Vector.h>
#include <iostream>
#include <vector>
#include <array>


//! Helper function to save PNG images with a suitable name
void savePngs (const std::string& logpath, const std::string& name,
               unsigned int frameN, morph::Visual& v) {
    std::stringstream ff1;
    ff1 << logpath << "/" << name<< "_";
    ff1 << std::setw(5) << std::setfill('0') << frameN;
    ff1 << ".png";
    v.saveImage (ff1.str());
}

std::array<float, 4> GaussFit(std::vector<float> x,std::vector<float> y){

    // fit based on https://math.stackexchange.com/questions/1292889/parameters-estimation-for-gaussian-function-with-offset
    int n = x.size();

    float S = 0.0;
    float T = 0.0;

    float v0 = 0.0;
    float v1 = 0.0;
    float v2 = 0.0;
    float v3 = 0.0;

    float m00 = 0.0;
    float m01 = 0.0;
    float m02 = 0.0;
    float m03 = 0.0;
    float m10 = 0.0;
    float m11 = 0.0;
    float m12 = 0.0;
    float m13 = 0.0;
    float m20 = 0.0;
    float m21 = 0.0;
    float m22 = 0.0;
    float m23 = 0.0;
    float m30 = 0.0;
    float m31 = 0.0;
    float m32 = 0.0;
    float m33 = 0.0;


    float x1sq = x[0]*x[0];

    for(int k=1;k<n;k++){


        S=S+0.5*(y[k]+y[k-1])*(x[k]-x[k-1]);
        T=T+0.5*(x[k]*y[k]+x[k-1]*y[k-1])*(x[k]-x[k-1]);

        float dx = x[k]-x[0];
        float dy = y[k]-y[0];
        float dxsq = x[k]*x[k]-x1sq;

        v0 += S*dy;
        v1 += T*dy;
        v2 += dxsq*dy;
        v3 += dx*dy;

        m00 += S*S;
        m01 += S*T;
        m02 += S*dxsq;
        m03 += S*dx;
        m11 += T*T;
        m12 += T*dxsq;
        m13 += T*dx;
        m22 += dxsq*dxsq;
        m23 += dxsq*dx;
        m33 += dx*dx;
    }

    m10 = m01;
    m20 = m02;
    m21 = m12;
    m30 = m03;
    m31 = m13;
    m32 = m23;

    arma::Mat M = { {m00,m01,m02,m03},
                    {m10,m11,m12,m13},
                    {m20,m21,m22,m23},
                    {m30,m31,m32,m33} };
    arma::Col V = {v0,v1,v2,v3};


    arma::Mat<float> I = arma::inv(M);
    arma::Col C = I * V;

    float a = - C[0]/C[1];
    float b = -2.0/C[1];

    float sumy = 0.0;
    float sumt = 0.0;
    float sumtsq = 0.0;
    float sumty = 0.0;

    float t = 0.0;
    for(int k=0;k<n;k++){
        t = exp(-((x[k]-a)*(x[k]-a))/b);
        sumt += t;
        sumtsq += t*t;
        sumty += t*y[k];
        sumy += y[k];

    }

    arma::Mat P = { {sumtsq,sumt},
                    {sumt,(float)n} };
    arma::Col Q = {sumty,sumy};

    arma::Col<float> D = arma::inv(P) * Q;

    float c = D[0];
    float h = D[1];

    std::array<float, 4> coeffs = {a,b,c,h};
    return coeffs;

}

#include "opencv2/opencv.hpp"
#include <morph/HexGrid.h>
#include <morph/RD_Base.h>
#include <morph/Config.h>

template <class Flt>
class Network{

    /*
        High-level wrapper for specifying a network so that a simulation can be built by calling the methods (e.g., step/map) in a given order.
    */

    public:
        int time;
        float spatialScale;
        morph::Config conf;

    Network(morph::Config conf){
        time = 0;
        spatialScale = 1.0;
        this->conf = conf;
    };

    virtual void stepHidden(bool){
        // overwrite this in derived classes to expose the hidden transforms within a network (e.g., thalamic processes between sensors and cortex) to other objects that have a pointer to the network
    }

    ~Network(void){

    }

};



template <class Flt>
class Projection{

/*
    A projection class for connecting units on a source sheet to units on a destination sheet with topographically aligned weighted connections from a radius of units on the source sheet to each destination sheet unit.
*/

public:

    morph::HexGrid* hgSrc;
    morph::HexGrid* hgDst;
    Flt radius;                     // radius within which connections are made
    Flt strength;                   // strength of projection - multiplication after dot products
    Flt alpha;                      // learning rate
    unsigned int nSrc;              // number of units on source sheet
    unsigned int nDst;              // number of units on destination sheet

    std::vector<unsigned int> counts;                // number of connections in connection field for each unit
    std::vector<Flt> norms;                // 1./counts
    std::vector<Flt> alphas;                // learning rates for each unit may depend on e.g., the number of connections
    std::vector<std::vector<unsigned int> > srcId;            // identity of conneted units on the source sheet
    std::vector<std::vector<Flt> > weights;        // connection weights
    std::vector<std::vector<Flt> > distances;        // pre-compute distances between units in source and destination sheets
    std::vector<Flt> field;                // current activity patterns
    std::vector<Flt*> fSrc;                // pointers to the field elements on the source sheet
    std::vector<Flt*> fDst;                // pointers to the field elements on the destination sheet
    std::vector<double> weightPlot;            // for constructing activity plots
    bool normalizeAlphas;            // whether to normalize learning rate by individual unit connection density

    Projection(std::vector<Flt*> fSrc, std::vector<Flt*> fDst, morph::HexGrid* hgSrc, morph::HexGrid* hgDst, Flt radius, Flt strength, Flt alpha, Flt sigma, bool normalizeAlphas, Flt jitter){

        // Initialise the class with random weights (if sigma>0, the weights have a Gaussian pattern, else uniform random)

        this->fSrc = fSrc;
        this->fDst = fDst;
        this->hgSrc = hgSrc;
        this->hgDst = hgDst;
        this->radius = radius;
        this->strength = strength;
        this->alpha = alpha;
        this->normalizeAlphas = normalizeAlphas;

        nDst = hgDst->vhexen.size();
        nSrc = hgSrc->vhexen.size();

        field.resize(nDst);
        counts.resize(nDst);
        norms.resize(nDst);
        srcId.resize(nDst);
        weights.resize(nDst);
        alphas.resize(nDst);
        distances.resize(nDst);
        weightPlot.resize(nSrc);

        Flt radiusSquared = radius*radius;    // precompute for speed

        double OverTwoSigmaSquared = 1./(sigma*sigma*2.0);    // precompute normalisation constant

    // initialize connections for each destination sheet unit
    #pragma omp parallel for
        for(unsigned int i=0;i<nDst;i++){
            Flt jitterx = jitter * (morph::Tools::randDouble()-0.5);
            Flt jittery = jitter * (morph::Tools::randDouble()-0.5);
            for(unsigned int j=0;j<nSrc;j++){
                Flt dx = (hgSrc->vhexen[j]->x+jitterx-hgDst->vhexen[i]->x);
                Flt dy = (hgSrc->vhexen[j]->y+jittery-hgDst->vhexen[i]->y);
                Flt distSquared = dx*dx+dy*dy;
                if (distSquared<radiusSquared){
                    counts[i]++;
                    srcId[i].push_back(j);
                    Flt w = 1.0;
                    if(sigma>0.){
                        w = exp(-distSquared*OverTwoSigmaSquared);
                    }
                    weights[i].push_back(w);
                    distances[i].push_back(sqrt(distSquared));
                }
            }
            norms[i] = 1.0/(Flt)counts[i];
            alphas[i] = alpha;
            if(normalizeAlphas){
                alphas[i] *= norms[i];
            }

        }
    }



    void getWeightedSum(void){
    /*
        Dot product of each weight vector with the corresponding source sheet field values, multiplied by the strength of the projection
    */
#pragma omp parallel for
        for(unsigned int i=0;i<nDst;i++){
            field[i] = 0.;
            for(unsigned int j=0;j<counts[i];j++){
                field[i] += *fSrc[srcId[i][j]]*weights[i][j];
            }
            field[i] *= strength;
        }
    }


    void setFieldToWeights(int i){
        for(unsigned int j=0;j<counts[i];j++){
            *fSrc[srcId[i][j]] = weights[i][j];
        }
    }

    void learn(void){
    /*
     Hebbian adaptation of the weights
    */
    if(alpha>0.0){
    #pragma omp parallel for
        for(unsigned int i=0;i<nDst;i++){
            for(unsigned int j=0;j<counts[i];j++){
                weights[i][j] += *fSrc[srcId[i][j]] * *fDst[i] * alphas[i];
            }
        }
    }
    }

    std::vector<Flt> getWeights(void){
        std::vector<Flt> weightStore;
        for(unsigned int i=0;i<nDst;i++){
            for(unsigned int j=0;j<counts[i];j++){
                weightStore.push_back(weights[i][j]);
            }
        }
        return weightStore;
    }

    void setWeights(std::vector<Flt> weightStore){
        int k=0;
        for(unsigned int i=0;i<nDst;i++){
            for(unsigned int j=0;j<counts[i];j++){
                weights[i][j] = weightStore[k];
                k++;
            }
        }
    }

    void renormalize(void){

    #pragma omp parallel for
        for(unsigned int i=0;i<nDst;i++){
        Flt sumWeights = 0.0;
            for(unsigned int j=0;j<counts[i];j++){
                sumWeights += weights[i][j];
            }
            for(unsigned int j=0;j<counts[i];j++){
                weights[i][j] /= sumWeights;
            }
        }
    }

    void multiplyWeights(int i, double scale){
    #pragma omp parallel for
        for(unsigned int j=0;j<counts[i];j++){
            weights[i][j] *= scale;
        }
    }

    std::vector<double> getWeightPlot(int i){

        #pragma omp parallel for
        for(unsigned int j=0;j<weightPlot.size();j++){
            weightPlot[j] = 0.;
        }
        #pragma omp parallel for
        for(unsigned int j=0;j<counts[i];j++){
            weightPlot[srcId[i][j]] = weights[i][j];
        }
        return weightPlot;
    }

};


template <class Flt>
class RD_Sheet : public morph::RD_Base<Flt>
{
public:

    std::vector<Projection<Flt>> Projections;
    alignas(alignof(std::vector<Flt>)) std::vector<Flt> X;
    alignas(alignof(std::vector<Flt*>)) std::vector<Flt*> Xptr;
    std::vector<std::vector<int> > P;            // identity of projections to (potentially) joint normalize

    virtual void init (void) {
        this->stepCount = 0;
        this->zero_vector_variable (this->X);
    }

    virtual void allocate (void) {
        morph::RD_Base<Flt>::allocate();
        this->resize_vector_variable (this->X);
        for (unsigned int hi=0; hi<this->nhex; ++hi) {
            Xptr.push_back(&this->X[hi]);
        }
    }


    void addProjection(std::vector<Flt*> inXptr, morph::HexGrid* hgSrc, float radius, float strength, float alpha, float sigma, bool normalizeAlphas){
        Projections.push_back(Projection<Flt>(inXptr, this->Xptr, hgSrc, this->hg, radius, strength, alpha, sigma, normalizeAlphas, 0.0));
    }


    void addProjection(std::vector<Flt*> inXptr, morph::HexGrid* hgSrc, float radius, float strength, float alpha, float sigma, bool normalizeAlphas, float jitter){
        Projections.push_back(Projection<Flt>(inXptr, this->Xptr, hgSrc, this->hg, radius, strength, alpha, sigma, normalizeAlphas, jitter));
    }

    void zero_X (void) {
        #pragma omp parallel for
        for (unsigned int hi=0; hi<this->nhex; ++hi) {
            this->X[hi] = 0.;
        }
    }

    void clip_X(Flt lower, Flt upper){
        #pragma omp parallel for
        for (unsigned int hi=0; hi<this->nhex; ++hi) {
            if(this->X[hi]>upper){
                this->X[hi]=upper;
            }
            if(this->X[hi]<lower){
                this->X[hi]=lower;
            }

        }
    }

    void norm_X (void) {
        Flt maxVal = -1e9;
        #pragma omp parallel for
        for (unsigned int hi=0; hi<this->nhex; ++hi) {
            if(this->X[hi]>maxVal){ maxVal = this->X[hi]; }
        }
        maxVal = 1./maxVal;
        for (unsigned int hi=0; hi<this->nhex; ++hi) {
            this->X[hi]*=maxVal;
        }
    }

    void printMinMax (void) {
        Flt maxVal = -1e9;
        Flt minVal = +1e9;

        #pragma omp parallel for
        for (unsigned int hi=0; hi<this->nhex; ++hi) {
            if(this->X[hi]<minVal){ minVal = this->X[hi]; }
            if(this->X[hi]>maxVal){ maxVal = this->X[hi]; }
        }
        std::cout<<"Min: "<<minVal<<", Max: "<<maxVal<<std::endl;
    }

    void setNormalize(std::vector<int> proj){
        for(unsigned int p=0; p<proj.size();p++){
            for(unsigned int i=0; i<this->P.size();i++){
                for(unsigned int j=0; j<this->P[i].size();j++){
                    if(proj[p]==this->P[i][j]){
                        std::cout<<"Caution - projection may be mutiply joint normalized"<<std::endl;
                    }
                }
            }
        }
        this->P.push_back(proj);
    }

     void renormalize(void){

        for(unsigned int proj=0;proj<this->P.size();proj++){
        #pragma omp parallel for
            for(unsigned int i=0;i<this->nhex;i++){
                Flt sumWeights = 0.0;
                for(unsigned int p=0;p<this->P[proj].size();p++){
                    for(unsigned int j=0;j<this->Projections[this->P[proj][p]].counts[i];j++){
                        sumWeights += this->Projections[this->P[proj][p]].weights[i][j];
                    }
                }
                for(unsigned int p=0;p<this->P[proj].size();p++){
                    for(unsigned int j=0;j<this->Projections[this->P[proj][p]].counts[i];j++){
                        this->Projections[this->P[proj][p]].weights[i][j] /= sumWeights;
                    }
                }
            }
        }
     }

};

template <class Flt>
class NormByFirstProjection : public RD_Sheet<Flt>
{

public:
    Flt K;

    NormByFirstProjection(void){
        this->K = 1.0; // control offset
    }

    virtual void step (void) {

        this->stepCount++;
        this->zero_X();
        // indexing from 1 here assumes first projection is the divisive normalizer
        for(unsigned int i=1;i<this->Projections.size();i++){
            this->Projections[i].getWeightedSum();
        }
        // indexing from 1 here assumes first projection is the divisive normalizer
        for(unsigned int i=1;i<this->Projections.size();i++){
#pragma omp parallel for
            for (unsigned int hi=0; hi<this->nhex; ++hi) {
                this->X[hi] += this->Projections[i].field[hi];
            }
        }
#pragma omp parallel for
        for (unsigned int hi=0; hi<this->nhex; ++hi) {
            if(this->X[hi]<0){ this->X[hi]=0; }
        }
        // X is now the numerator

        this->Projections[0].getWeightedSum();
#pragma omp parallel for
        for (unsigned int hi=0; hi<this->nhex; ++hi) {
            if(this->Projections[0].field[hi]<0.0) { this->Projections[0].field[hi] = 0.0; }
        }
        // Projections[0].field is now the denominator

#pragma omp parallel for
        for(int i=0;i<this->nhex;i++){
            this->X[i] = fmax( this->X[i] / (this->K + this->Projections[0].field[i]) ,0.);
        }
    }

};


template <class Flt>
class CortexSOM : public RD_Sheet<Flt>
{
    public:

    Flt beta, lambda, mu, oneMinusBeta, thetaInit;

    alignas(alignof(std::vector<Flt>)) std::vector<Flt> Xavg;
    alignas(alignof(std::vector<Flt>)) std::vector<Flt> Theta;

    CortexSOM(void){

    }

    virtual void init (void) {
        oneMinusBeta = 1.-beta;
        this->stepCount = 0;
        this->zero_vector_variable (this->X);
        this->zero_vector_variable (this->Xavg);
        this->zero_vector_variable (this->Theta);
    }

    virtual void allocate (void) {
        morph::RD_Base<Flt>::allocate();
        this->resize_vector_variable (this->X);
        this->resize_vector_variable (this->Xavg);
        this->resize_vector_variable (this->Theta);
        for (unsigned int hi=0; hi<this->nhex; ++hi) {
            this->Xptr.push_back(&this->X[hi]);
        }
        for (unsigned int hi=0; hi<this->nhex; ++hi) {
            this->Xavg[hi] = mu;
            this->Theta[hi] = thetaInit;
        }
    }


    virtual void step (void) {
        this->stepCount++;

        for(unsigned int i=0;i<this->Projections.size();i++){
            this->Projections[i].getWeightedSum();
        }

        this->zero_X();

        for(unsigned int i=0;i<this->Projections.size();i++){
        #pragma omp parallel for
            for (unsigned int hi=0; hi<this->nhex; ++hi) {
                this->X[hi] += this->Projections[i].field[hi];
            }
        }

        #pragma omp parallel for
        for (unsigned int hi=0; hi<this->nhex; ++hi) {
            this->X[hi] = this->X[hi]-this->Theta[hi];
            if(this->X[hi]<0.0){
                this->X[hi] = 0.0;
            }
        }
    }

    virtual void step (std::vector<int> projectionIDs) {
        this->stepCount++;

        for(unsigned int i=0;i<projectionIDs.size();i++){
            this->Projections[projectionIDs[i]].getWeightedSum();
        }

        this->zero_X();

        for(unsigned int i=0;i<projectionIDs.size();i++){
        #pragma omp parallel for
            for (unsigned int hi=0; hi<this->nhex; ++hi) {
                this->X[hi] += this->Projections[projectionIDs[i]].field[hi];
            }
        }

        #pragma omp parallel for
        for (unsigned int hi=0; hi<this->nhex; ++hi) {
            this->X[hi] = this->X[hi]-this->Theta[hi];
            if(this->X[hi]<0.0){
                this->X[hi] = 0.0;
            }
        }
    }


    void homeostasis(void){

     #pragma omp parallel for
        for (unsigned int hi=0; hi<this->nhex; ++hi) {
            this->Xavg[hi] = oneMinusBeta*this->X[hi] + beta*this->Xavg[hi];
        }
     #pragma omp parallel for
        for (unsigned int hi=0; hi<this->nhex; ++hi) {
            this->Theta[hi] = this->Theta[hi] + lambda*(this->Xavg[hi]-mu);
        }
    }



};

template <class Flt>
class PatternGenerator_Sheet : public RD_Sheet<Flt>
{
public:
    PatternGenerator_Sheet(){ }

    virtual void step (void) {
        this->stepCount++;
    }


    void Gaussian(double x_center, double y_center, double theta, double sigmaA, double sigmaB, double amplitude){

        double cosTheta = cos(theta);
        double sinTheta = sin(theta);
        double overSigmaA = 1./sigmaA;
        double overSigmaB = 1./sigmaB;
        #pragma omp parallel for
        for (unsigned int hi=0; hi<this->nhex; ++hi) {
            Flt dx = this->hg->vhexen[hi]->x-x_center;
            Flt dy = this->hg->vhexen[hi]->y-y_center;
            this->X[hi] = amplitude*exp(-((dx*cosTheta-dy*sinTheta)*(dx*cosTheta-dy*sinTheta))*overSigmaA
                              -((dx*sinTheta+dy*cosTheta)*(dx*sinTheta+dy*cosTheta))*overSigmaB);
        }
    }


    void Gaussian(std::vector<double> x_center, std::vector<double> y_center, std::vector<double> theta, std::vector<double> sigmaA, std::vector<double> sigmaB, std::vector<double> amplitude){

        int n = x_center.size();

        for (unsigned int hi=0; hi<this->nhex; ++hi) {
            this->X[hi] = 0.0;
        }

        for(int i=0;i<n;i++){
            double cosTheta = cos(theta[i]);
            double sinTheta = sin(theta[i]);
            double overSigmaA = 1./sigmaA[i];
            double overSigmaB = 1./sigmaB[i];
            #pragma omp parallel for
            for (unsigned int hi=0; hi<this->nhex; ++hi) {
                Flt dx = this->hg->vhexen[hi]->x-x_center[i];
                Flt dy = this->hg->vhexen[hi]->y-y_center[i];
                double val = amplitude[i]*exp(-((dx*cosTheta-dy*sinTheta)*(dx*cosTheta-dy*sinTheta))*overSigmaA
                                  -((dx*sinTheta+dy*cosTheta)*(dx*sinTheta+dy*cosTheta))*overSigmaB);
                if(val > this->X[hi]){
                    this->X[hi] = val;
                }
            }
        }
    }


    void Grating(double theta, double phase, double width, double amplitude){

        double cosTheta = cos(theta);
        double sinTheta = sin(theta);

        #pragma omp parallel for
        for (unsigned int hi=0; hi<this->nhex; ++hi) {

            this->X[hi] = amplitude * sin( width * (this->hg->vhexen[hi]->x*sinTheta + this->hg->vhexen[hi]->y*cosTheta + phase) );

                //(dx*cosTheta-dy*sinTheta)

        }
    }


};




// This helper function is general-purpose and should really be moved into morphologica
std::vector<float> getPolyPixelVals(cv::Mat frame, std::vector<cv::Point> pp){
    cv::Point pts[4] = {pp[0],pp[1],pp[2],pp[3]};
    cv::Mat mask = cv::Mat::zeros(frame.rows, frame.cols, CV_8UC3);
    fillConvexPoly(mask, pts, 4, cv::Scalar(255,255,255) );
    cv::Mat result, resultGray;
    frame.copyTo(result,mask);
    //cvtColor(result,resultGray,CV_BGR2GRAY);
    cvtColor(result,resultGray,cv::COLOR_BGR2GRAY/*CV_BGR2GRAY*/);
    std::vector<cv::Point2i> positives;
    findNonZero(resultGray, positives);
    std::vector<float> polyPixelVals(positives.size());
    for(int j=0;j<positives.size();j++){
         cv::Scalar pixel = resultGray.at<uchar>(positives[j]);
        polyPixelVals[j] = (float)pixel.val[0]/255.;
    }
    return polyPixelVals;
}



class Square{
    public:
        int xid, yid;
        double x, y;
        double X;

    Square(int xid, int yid, double x, double y){
        this->xid = xid;
        this->yid = yid;
        this->x = x;
        this->y = y;
        X = 0.;
    }
};

class CartGrid{
    public:
        int n, nx, ny;
        std::vector<Square> vsquare;

    CartGrid(void){

    }

    CartGrid(int nx, int ny){
        init(ny, ny);
    }

    void init(int nx, int ny){

        this->nx = nx;
        this->ny = ny;
        n = nx*ny;

        double maxDim = (double)std::max(nx,ny);

        int k=0;
        for(int i=0;i<nx;i++){
            double xpos = ((double)i/((double)maxDim-1))-0.5;
            for(int j=0;j<ny;j++){
                double ypos = ((double)j/((double)maxDim-1))-0.5;
                vsquare.push_back(Square(i,j,xpos,ypos));
                k++;
            }
        }
    }
};


template <class Flt>
class HexCartSampler : public RD_Sheet<Flt>
{

public:
    CartGrid C;
    cv::VideoCapture cap;
    Flt radius, sigma;
    std::vector<std::vector<unsigned int> > srcId;
    std::vector<std::vector<Flt> > weights;
    std::vector<std::vector<Flt> > distances;
    std::vector<unsigned int> counts;
    std::vector<Flt> norms;
    Flt strength;
    std::vector<cv::Point> mask;
    unsigned int stepsize;
    int patternsWid;
    int nPatterns;

    std::vector<std::vector<Flt> > PreLoadedPatterns;

    HexCartSampler(void){

    }

    HexCartSampler(int nx, int ny, Flt radius, Flt sigma){
        init(nx,ny,radius,sigma);
    }

    void initProjection(int nx, int ny, Flt radius, Flt sigma){

        C.init(nx, ny);

        this->radius = radius;
        this->sigma = sigma;
        this->strength = 1.0;
        srcId.resize(this->nhex);
        counts.resize(this->nhex);
        weights.resize(this->nhex);
        distances.resize(this->nhex);
        norms.resize(this->nhex);

        Flt radiusSquared = radius*radius;    // precompute for speed
        Flt OverTwoSigmaSquared = 1./(2.0*sigma*sigma);
        #pragma omp parallel for
        for(unsigned int i=0;i<this->nhex;i++){
            for(unsigned int j=0;j<C.n;j++){
                Flt dx = (this->hg->vhexen[i]->x-C.vsquare[j].x);
                Flt dy = (this->hg->vhexen[i]->y-C.vsquare[j].y);
                Flt distSquared = dx*dx+dy*dy;
                if (distSquared<radiusSquared){
                    counts[i]++;
                    srcId[i].push_back(j);
                    Flt w = 1.0;
                    if(sigma>0.){
                        w = exp(-distSquared*OverTwoSigmaSquared);
                    }
                    weights[i].push_back(w);
                    distances[i].push_back(sqrt(distSquared));
                }
            }
            norms[i] = 1.0/(Flt)counts[i];
            for(unsigned int j=0;j<counts[i];j++){
                weights[i][j] *= norms[i];
            }
        }

    }

    int initCamera(int xoff, int yoff, int stepsize){
        this->stepsize = stepsize;
        mask.resize(4);
        mask[0] = cv::Point(xoff,yoff);
        mask[1] = cv::Point(xoff+stepsize*C.nx-1,yoff);
        mask[2] = cv::Point(xoff+stepsize*C.nx-1,yoff+stepsize*C.ny-1);
        mask[3] = cv::Point(xoff,yoff+stepsize*C.ny-1);
        return cap.open(0);
    }

    virtual void step (void) {
        this->zero_X();
        #pragma omp parallel for
        for(unsigned int i=0;i<this->nhex;i++){
            for(unsigned int j=0;j<counts[i];j++){
                this->X[i] += C.vsquare[srcId[i][j]].X*weights[i][j];
            }
            this->X[i] *= strength;
        }
    }

    void stepCamera(void){
        cv::Mat frame;
        cap >> frame;
        std::vector<float> img = getPolyPixelVals(frame,mask);
        std::vector<Flt> pat((img.size()/stepsize),0.);
        int iter = stepsize*C.ny*stepsize;
        int k=0;
        for(int i=0;i<C.nx;i++){
            int I = (C.nx-i-1)*stepsize;
            for(int j=0;j<C.ny;j++){
                C.vsquare[k].X = img[(C.ny-j-1)*iter+I];
                k++;
            }
        }
        step();
    }

    /*
    void preloadPatterns(std::string filename){
        std::stringstream fname; fname << filename;
        morph::HdfData data(fname.str(),1);
        std::vector<Flt> tmp;
        data.read_contained_vals ("P", tmp);
        int nPat = tmp.size()/C.n;
        PreLoadedPatterns.resize(nPat,std::vector<Flt>(C.n,0.));
        int k=0;
        for(int i=0;i<nPat;i++){
            for(int j=0;j<C.n;j++){
                PreLoadedPatterns[i][j] = tmp[k];
                k++;
            }
        }
    }
    */

    void preloadPatterns(std::string filename){
        std::stringstream fname; fname << filename;
        morph::HdfData data(fname.str(),1);
        std::vector<Flt> tmp;
        data.read_contained_vals ("P", tmp);
        std::vector<int> dims;
        data.read_contained_vals ("dims", dims);
        nPatterns = dims[2];
        patternsWid = dims[0]; // could/should check/deal with non-square images but not doing so now
        PreLoadedPatterns.resize(nPatterns,std::vector<Flt>(patternsWid*patternsWid,0.));
        int k=0;
        for(int i=0;i<nPatterns;i++){
            for(int j=0;j<patternsWid*patternsWid;j++){
                //if((tmp[k]<0.0) || (tmp[k]>1.0)){ std::cout<<"preloaded pattern with bad value: pixel "<<k<<", val: "<<tmp[k]<<std::endl; }
                PreLoadedPatterns[i][j] = tmp[k];
                k++;
            }
        }
    }

    /*
    void stepPreloaded(int p){

        for(int i=0;i<C.n;i++){
            C.vsquare[i].X = PreLoadedPatterns[p][i];
        }
        step();
    }
    */

    void stepPreloaded(int p){

        // Jitter if the sample width is less than the image width
        int offx = floor(morph::Tools::randDouble()*(patternsWid-C.nx));
        int offy = floor(morph::Tools::randDouble()*(patternsWid-C.ny));

        int k=0;
        for(int i=0;i<C.nx;i++){
            for(int j=0;j<C.ny;j++){
                int ind = (patternsWid + offy + i)*patternsWid + offx + j;
                C.vsquare[k].X = PreLoadedPatterns[p][ind];
                k++;
            }
        }
        step();
    }

    void stepPreloaded(int p, float offx, float offy){ // should supply xoff,yoff in sheet coords and convert inside this function really.
        int k=0;
        for(int i=0;i<C.nx;i++){
            for(int j=0;j<C.ny;j++){
                int ind = (patternsWid + offy + i)*patternsWid + offx + j;
                C.vsquare[k].X = PreLoadedPatterns[p][ind];
                k++;
            }
        }
        step();
    }

    void stepPreloaded(void){
        int p = floor(morph::Tools::randDouble()*PreLoadedPatterns.size());
        stepPreloaded(p);
    }
};



template <class Flt>
class CartHexSampler {

public:

    CartGrid C;
    std::vector<Flt*> fSrc;
    morph::HexGrid* hgSrc;
    Flt radius, sigma;
    std::vector<std::vector<unsigned int> > srcId;
    std::vector<std::vector<Flt> > weights;
    std::vector<std::vector<Flt> > distances;
    std::vector<unsigned int> counts;
    std::vector<Flt> norms;
    Flt strength;
    unsigned int stepsize;
    float rangex, rangey;

    void initProjection(int nx, int ny, float rangex, float rangey, std::vector<Flt*> fSrc, morph::HexGrid* hgSrc, Flt radius, Flt sigma){

        C.init(nx,ny);

        this->fSrc = fSrc;
        this->hgSrc = hgSrc;

        this->radius = radius;
        this->sigma = sigma;
        this->strength = 1.0;
        srcId.resize(C.n);
        counts.resize(C.n,0);
        weights.resize(C.n);
        distances.resize(C.n);
        norms.resize(C.n);

        Flt radiusSquared = radius*radius;    // precompute for speed
        Flt OverTwoSigmaSquared = 1./(2.0*sigma*sigma);
        #pragma omp parallel for
        for(unsigned int i=0;i<C.n;i++){
            for(unsigned int j=0;j<hgSrc->vhexen.size();j++){
                Flt dx = (hgSrc->vhexen[j]->x-C.vsquare[i].x*rangex);
                Flt dy = (hgSrc->vhexen[j]->y-C.vsquare[i].y*rangey);
                Flt distSquared = dx*dx+dy*dy;
                if (distSquared<radiusSquared){
                    counts[i]++;
                    srcId[i].push_back(j);
                    Flt w = 1.0;
                    if(sigma>0.){
                        w = exp(-distSquared*OverTwoSigmaSquared);
                    }
                    weights[i].push_back(w);
                    distances[i].push_back(sqrt(distSquared));
                }
            }
            float sumw = 0.;
            for(int j=0;j<weights[i].size();j++){
                sumw += weights[i][j];
            }
            for(int j=0;j<weights[i].size();j++){
                weights[i][j] /= sumw;
            }

        }
    }

    virtual void step (void) {
        #pragma omp parallel for
        for(unsigned int i=0;i<C.n;i++){
            C.vsquare[i].X = 0;
            for(unsigned int j=0;j<counts[i];j++){
                C.vsquare[i].X += *fSrc[srcId[i][j]]*weights[i][j];
            }
            C.vsquare[i].X *= strength;
        }
    }

    void stepOrientation (void) {

        /*
            Does vector average of pi-periodic values
        */

        #pragma omp parallel for
        for(unsigned int i=0;i<C.n;i++){
            float vx = 0.0;
            float vy = 0.0;
            for(unsigned int j=0;j<counts[i];j++){
                float val = *fSrc[srcId[i][j]]*2.0;
                vx += cos(val);
                vy += sin(val);
            }
            C.vsquare[i].X = 0.5*(atan2(vy,vx)+M_PI);
        }
    }

    //
    cv::Mat getDifferenceImage(std::vector<Flt> A, std::vector<Flt> B){

        std::vector<Flt> diff = A;
        for(int i=0;i<diff.size();i++){
            diff[i] -= B[i];
        }

        float maxV = -1e9;
        float minV = +1e9;
        for(int i=0;i<C.n;i++){
            if(maxV<fabs(diff[i])){
                maxV = fabs(diff[i]);
            } // THIS WAY ENSURE THAT ZERO DIFF ALWAYS MAPS TO VAL OF 0.5
            if(minV>diff[i]){ minV = diff[i]; }
        }
        float scale = 1./(2*maxV);

        cv::Mat I = cv::Mat::zeros(C.nx,C.ny,CV_32F);
        int k=0;
        for(int i=0;i<C.nx;i++){
            for(int j=C.ny; j-->0;){ // invert y axis of image
                I.at<float>(j,i) = (diff[k]+maxV)*scale;
                k++;
            }
        }
        return I;
    }

    //
    std::vector<std::vector<float> > fft(cv::Mat I, int nbins, int gaussBlur, bool view){

        if(view){
            cv::namedWindow("analysis", cv::WINDOW_NORMAL);
            cv::resizeWindow("analysis", I.cols*3,I.rows*3);
            cv::imshow("analysis", I);
        }

        // https://www.youtube.com/watch?v=MyX9UPu4cAA
        cv::Mat originalComplex[2] = {I, cv::Mat::zeros(I.size(), CV_32F)};

        I += cv::Scalar::all(-0.5);

        int gb = gaussBlur*2+1;

        cv::GaussianBlur(I, I, cv::Size(gb,gb), 0);

        cv::Mat dftReady;
        cv::merge(originalComplex, 2, dftReady);
        cv::Mat dftOfOriginal;
        cv::dft(dftReady, dftOfOriginal, cv::DFT_COMPLEX_OUTPUT);
        cv::Mat splitArray[2] = {cv::Mat::zeros(dftOfOriginal.size(), CV_32F), cv::Mat::zeros(dftOfOriginal.size(), CV_32F)};
        cv::split(dftOfOriginal, splitArray);
        cv::Mat dftMagnitude;
        cv::magnitude(splitArray[0], splitArray[1], dftMagnitude);
        dftMagnitude = dftMagnitude.mul(dftMagnitude);
        dftMagnitude += cv::Scalar::all(1);
        cv::log(dftMagnitude, dftMagnitude);
        cv::normalize(dftMagnitude, dftMagnitude, 0, 1, cv::NORM_MINMAX/*CV_MINMAX*/);
        cv::Mat magI = dftMagnitude;

        // crop the spectrum, if it has an odd number of rows or columns
        magI = magI(cv::Rect(0, 0, magI.cols & -2, magI.rows & -2));
        // rearrange the quadrants of Fourier image  so that the origin is at the image center
        int cx = magI.cols/2;
        int cy = magI.rows/2;
        cv::Mat q0(magI, cv::Rect(0, 0, cx, cy));   // Top-Left - Create a ROI per quadrant
        cv::Mat q1(magI, cv::Rect(cx, 0, cx, cy));  // Top-Right
        cv::Mat q2(magI, cv::Rect(0, cy, cx, cy));  // Bottom-Left
        cv::Mat q3(magI, cv::Rect(cx, cy, cx, cy)); // Bottom-Right
        cv::Mat tmp;                           // swap quadrants (Top-Left with Bottom-Right)
        q0.copyTo(tmp);
        q3.copyTo(q0);
        tmp.copyTo(q3);
        q1.copyTo(tmp);                    // swap quadrant (Top-Right with Bottom-Left)
        q2.copyTo(q1);
        tmp.copyTo(q2);

        if(view){
            cv::Mat magIinv = cv::Scalar::all(1) - magI;
            cv::namedWindow("FFT", cv::WINDOW_NORMAL);
            cv::resizeWindow("FFT", C.nx*3,C.ny*3);
            cv::imshow("FFT", magIinv);
        }


        // collect data against fft location
        std::vector<float> V(magI.rows*magI.cols,0.);
        std::vector<float> D(magI.rows*magI.cols,0.);

        // dx and dy - frequency at which image is sampled
        float dx = 1./((float)magI.rows);
        float dy = 1./((float)magI.cols);
        {
            int k=0;
            for(int i=0;i<magI.rows;i++){
                // xrangeval -- range of frequencies of fft x-axis (Nyquist freq)
                float xrangeval = ((float)i/((float)magI.rows-1)-0.5);
                // ... in units of dx
                float distx = xrangeval/dx;
                for(int j=0;j<magI.cols;j++){
                    // yrangeval -- range of frequencies of fft y-axis (Nyquist freq)
                    float yrangeval = ((float)j/((float)magI.cols-1)-0.5);
                    // ... in units of dy
                    float disty = yrangeval/dy;
                    // frequency at the distance of the fft point from the origin (max=sqrt(0.5))
                    D[k] = pow(distx*distx+disty*disty,0.5);
                    V[k] = magI.at<float>(j,i);
                    k++;
                }
            }
        }

        float binmax = 1.0/(2*dx); // sample up to Nyquist frequency

        std::vector<float> binS(nbins);
        std::vector<float> binE(nbins);
        std::vector<float> binM(nbins);
        std::vector<float> binV(nbins,0.);
        std::vector<float> binC(nbins,0);
        for(int i=0;i<nbins;i++){
            binS[i] = binmax*(float)i/(float)nbins;
            binE[i] = binmax*(float)(i+1)/(float)nbins;
            binM[i] = binE[i]+(binS[i]-binE[i])*0.5;
        }

        for(int i=0;i<V.size();i++){
            for(int j=0;j<nbins;j++){
                if((binS[j]<=D[i]) && (D[i]<binE[j])){
                    binV[j] += V[i];
                    binC[j] +=1.0;
                }
            }
        }
        for(int i=0;i<nbins;i++){
            binV[i] *= 1./(binC[i]);
        }

        std::vector<std::vector<float> > rtn;
        rtn.push_back(binM);
        rtn.push_back(binV);
        return rtn;
    }


};

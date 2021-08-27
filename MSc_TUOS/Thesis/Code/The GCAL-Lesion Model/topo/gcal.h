#include "topo.h"
#include "general.h"                // INTEGRATE THESE FUNCTIONS INTO MORPHOLOGICA
#include <morph/Scale.h>
#include <morph/Vector.h>
#include <morph/HdfData.h>
#include <morph/Config.h>

class gcal : public Network<FLT> {

public:

    HexCartSampler<FLT> HCM;
    PatternGenerator_Sheet<FLT> IN;
    NormByFirstProjection<FLT> LGN_ON, LGN_OFF;
    CortexSOM<FLT> CX;
    bool homeostasis;
    unsigned int settle, nGauss;
    float beta, lambda, mu, thetaInit, xRange, yRange, afferAlpha, excitAlpha, inhibAlpha;
    float afferStrength, excitStrength, inhibStrength, LGNstrength, lrscale, pscale, gainControlStrength, gainControlOffset;
    float sigmaA, sigmaB, afferRadius, excitRadius, inhibRadius, afferSigma, excitSigma, inhibSigma, LGNCenterSigma, LGNSurroundSigma;

    gcal(morph::Config conf) : Network(conf) {

        // GET PARAMS FROM JSON
        homeostasis = conf.getBool("homeostasis", true);
        settle = conf.getUInt("settle", 16);

        // homeostasis
        beta = conf.getFloat("beta", 0.991);
        lambda = conf.getFloat("lambda", 0.01);
        thetaInit = conf.getFloat("thetaInit", 0.00);
        mu = conf.getFloat("mu", 0.024);

        xRange = conf.getFloat("xRange", 2.0);
        yRange = conf.getFloat("yRange", 2.0);
        nGauss = conf.getUInt("nGauss", 2);

        // learning rates
        lrscale = conf.getFloat("lrscale", 1.0);
        afferAlpha = conf.getFloat("afferAlpha", 0.1) * lrscale;
        excitAlpha = conf.getFloat("excitAlpha", 0.0) * lrscale;
        inhibAlpha = conf.getFloat("inhibAlpha", 0.3) * lrscale;

        // projection strengths
        pscale = conf.getFloat("pscale", 1.0);
        afferStrength = conf.getFloat("afferStrength", 1.5) * pscale;
        excitStrength = conf.getFloat("excitStrength", 1.7) * pscale;
        inhibStrength = conf.getFloat("inhibStrength", -1.4) * pscale;
        LGNstrength = conf.getFloat("retinaLGNstrength", 2.0);    // note this is 14 in paper
        gainControlStrength = conf.getFloat("gainControlStrength", 2.0);
        gainControlOffset = conf.getFloat("gainControlOffset", 0.11);

        // spatial params
        spatialScale = conf.getFloat("scale", 0.5);
        //float spatialScale=0.4;
        sigmaA = conf.getFloat("sigmaA", 1.0) * spatialScale;
        sigmaB = conf.getFloat("sigmaB", 0.3) * spatialScale;
        afferRadius = conf.getFloat("afferRadius", 0.27) * spatialScale;
        excitRadius = conf.getFloat("excitRadius", 0.1) * spatialScale;
        inhibRadius = conf.getFloat("inhibRadius", 0.23) * spatialScale;
        afferSigma = conf.getFloat("afferSigma", 0.270) * spatialScale;
        excitSigma = conf.getFloat("excitSigma", 0.025) * spatialScale;
        inhibSigma = conf.getFloat("inhibSigma", 0.075) * spatialScale;
        LGNCenterSigma = conf.getFloat("LGNCenterSigma", 0.037) * spatialScale;
        LGNSurroundSigma = conf.getFloat("LGNSuroundSigma", 0.150) * spatialScale;
        float TCjitter = conf.getFloat("TCjitter", 0.0) * spatialScale;


        bool normAlphas = true;

        // Mapping Between Hexagonal and Cartesian Sheet
        HCM.svgpath = conf.getString("IN_svgpath", "boundaries/trialmod.svg");
        HCM.init();
        HCM.allocate();

        // INPUT SHEET
        IN.svgpath = conf.getString("IN_svgpath", "boundaries/trialmod.svg");
        IN.init();
        IN.allocate();

        // LGN ON CELLS
        LGN_ON.K = gainControlOffset;
        LGN_ON.svgpath = conf.getString("LGN_svgpath", "boundaries/trialmod.svg");
        LGN_ON.init();
        LGN_ON.allocate();

        LGN_ON.addProjection(LGN_ON.Xptr, LGN_ON.hg, afferRadius, gainControlStrength, 0.0, 0.125,
                             false); // self-conn should be first

        LGN_ON.addProjection(IN.Xptr, IN.hg, afferRadius, +LGNstrength, 0.0, LGNCenterSigma, false);
        LGN_ON.addProjection(IN.Xptr, IN.hg, afferRadius, -LGNstrength, 0.0, LGNSurroundSigma, false);

        for (unsigned int i = 0; i < LGN_ON.Projections.size(); i++) {
            LGN_ON.Projections[i].renormalize();
        }

        LGN_OFF.K = gainControlOffset;
        LGN_OFF.svgpath = conf.getString("LGN_svgpath", "boundaries/trialmod.svg");
        LGN_OFF.init();
        LGN_OFF.allocate();


        LGN_OFF.addProjection(LGN_OFF.Xptr, LGN_OFF.hg, afferRadius, gainControlStrength, 0.0, 0.125,
                              false); // self-conn should be first

        LGN_OFF.addProjection(IN.Xptr, IN.hg, afferRadius, -LGNstrength, 0.0, LGNCenterSigma, false);
        LGN_OFF.addProjection(IN.Xptr, IN.hg, afferRadius, +LGNstrength, 0.0, LGNSurroundSigma, false);

        for (unsigned int i = 0; i < LGN_OFF.Projections.size(); i++) {
            LGN_OFF.Projections[i].renormalize();
        }

        // CORTEX SHEET
        CX.beta = beta;
        CX.lambda = lambda;
        CX.mu = mu;
        CX.thetaInit = thetaInit;
        CX.svgpath = conf.getString("CX_svgpath", "boundaries/trialmod.svg");
        CX.init();
        CX.allocate();

        // Modified following five lines for Zheng Yuan
        float lesionX = conf.getFloat("lesionX", 0.0) * spatialScale;
        float lesionY = conf.getFloat("lesionY", 0.0) * spatialScale;
        float lesionRadius = conf.getFloat("lesionRadius", 0.27) * spatialScale;
        CX.addProjection(LGN_ON.Xptr, LGN_ON.hg, afferRadius, afferStrength * 0.5, afferAlpha, afferSigma, normAlphas,
                         TCjitter, lesionX, lesionY, lesionRadius);
        CX.addProjection(LGN_OFF.Xptr, LGN_OFF.hg, afferRadius, afferStrength * 0.5, afferAlpha, afferSigma, normAlphas,
                         TCjitter, lesionX, lesionY, lesionRadius);

        CX.addProjection(CX.Xptr, CX.hg, excitRadius, excitStrength, excitAlpha, excitSigma, normAlphas);
        CX.addProjection(CX.Xptr, CX.hg, inhibRadius, inhibStrength, inhibAlpha, inhibSigma, normAlphas);

        // SETUP FIELDS FOR JOINT NORMALIZATION
        std::vector<int> p1(2, 0);
        p1[1] = 1;
        CX.setNormalize(p1);
        CX.setNormalize(std::vector<int>(1, 2));
        CX.setNormalize(std::vector<int>(1, 3));
        CX.renormalize();

    }

    void stepHidden(bool learning) {
        LGN_ON.step();
        LGN_OFF.step();
    }

    void stepAfferent(unsigned type) {
        switch (type) {
            case (0): { // Gaussians
                std::vector<double> x(nGauss);
                std::vector<double> y(nGauss);
                std::vector<double> t(nGauss);
                std::vector<double> sA(nGauss);
                std::vector<double> sB(nGauss);
                std::vector<double> amp(nGauss);
                for (int i = 0; i < nGauss; i++) {
                    x[i] = (morph::Tools::randDouble() - 0.5) * xRange;
                    y[i] = (morph::Tools::randDouble() - 0.5) * yRange;
                    t[i] = morph::Tools::randDouble() * M_PI;
                    sA[i] = sigmaA;
                    sB[i] = sigmaB;
                    amp[i] = 1.0;
                }
                IN.Gaussian(x, y, t, sA, sB, amp);
            }
                break;
            case (1): { // Preloaded
                HCM.stepPreloaded();
                IN.X = HCM.X;
            }
                break;
            case (2): { // Camera input
                HCM.stepCamera();
                IN.X = HCM.X;
            }
                break;
            default: {
                for (int i = 0; i < HCM.C.n; i++) {
                    HCM.C.vsquare[i].X = morph::Tools::randDouble();
                }
                HCM.step();
                IN.X = HCM.X;
            }
        }
        //stepHidden(false);
    }

    void stepCortex(bool learning) {
        CX.zero_X();
        for (unsigned int j = 0; j < settle; j++) {
            CX.step();
        }
        if (learning) {
            for (unsigned int p = 0; p < CX.Projections.size(); p++) { CX.Projections[p].learn(); }
            CX.renormalize();
            if (homeostasis) { CX.homeostasis(); }
            time++;
        }
    }

    void save(std::string filename) {
        std::stringstream fname;
        fname << filename;
        morph::HdfData data(fname.str());
        std::vector<int> timetmp(1, time);
        data.add_contained_vals("time", timetmp);
        for (unsigned int p = 0; p < CX.Projections.size(); p++) {
            std::vector<FLT> proj = CX.Projections[p].getWeights();
            std::stringstream ss;
            ss << "proj_" << p;
            data.add_contained_vals(ss.str().c_str(), proj);
        }
    }

    void load(std::string filename) {
        std::stringstream fname;
        fname << filename;
        morph::HdfData data(fname.str(), 1);
        std::vector<int> timetmp;
        data.read_contained_vals("time", timetmp);
        time = timetmp[0];
        for (unsigned int p = 0; p < CX.Projections.size(); p++) {
            std::vector<FLT> proj;
            std::stringstream ss;
            ss << "proj_" << p;
            data.read_contained_vals(ss.str().c_str(), proj);
            CX.Projections[p].setWeights(proj);
        }
        std::cout << "Loaded weights and modified time to " << time << std::endl;
    }

    ~gcal(void) {}


};

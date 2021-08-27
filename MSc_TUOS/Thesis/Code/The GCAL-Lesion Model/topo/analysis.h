//#include "topo.h"
//#include "general.h"

#include <morph/ShapeAnalysis.h>

class orientationPinwheelDensity {

public:

    CartHexSampler<FLT> CHM;
    std::vector<std::vector<FLT> > orResponse;
    std::vector<std::vector<FLT> > orResponseSampled;
    morph::ShapeAnalysis<FLT> shapeAn;

    float ROIwid, ROIpinwheelCount, IsoORfrequency, IsoORcolumnSpacing, sampleRange, gratingWidth, amplitude;
    std::vector<FLT> orPref, orSel, intersects, binVals, histogram, sfPref;
    std::vector<std::vector<FLT> > IsoORcontours;
    int nBins, nPhase, polyOrder, sampwid, gaussBlur, nOr;
    float pinwheelDensity;
    std::vector<float> Freq;

    PatternGenerator_Sheet<FLT> *In;
    CortexSOM<FLT> *Out;
    Network<FLT> *Net;

    //new public var
    std::vector<float> max_X;
    std::vector<float> max_Y;
    std::vector<std::vector<float>> retinotopy_data;

    orientationPinwheelDensity(Network<FLT> *Net, PatternGenerator_Sheet<FLT> *In, CortexSOM<FLT> *Out) {

        this->In = In;
        this->Out = Out;
        this->Net = Net;

        polyOrder = Net->conf.getUInt("polyOrder", 4);
        sampleRange = Net->conf.getFloat("sampleRange", 1.0);
        ROIwid = Net->conf.getFloat("ROIwid", 0.4);
        sampwid = Net->conf.getUInt("sampwid", 100);
        gaussBlur = Net->conf.getUInt("gaussBlur", 1);
        nBins = Net->conf.getUInt("nBins", 50);
        nPhase = Net->conf.getUInt("nPhase", 8);


        const Json::Value freq = Net->conf.getArray("freq");
        for (int i = 0; i < freq.size(); i++) {
            Freq.push_back(freq[i].asFloat());
        }

        gratingWidth = 7.5;
        amplitude = 1.0;
        nOr = 4;

        CHM.initProjection(sampwid, sampwid, ROIwid, ROIwid, Out->Xptr, Out->hg, 0.05, 0.001);
        orPref.resize(Out->nhex, 0.);
        orSel.resize(Out->nhex, 0.);
        sfPref.resize(Out->nhex, 0.);
        IsoORcontours.resize(2, std::vector<FLT>(Out->nhex, 0.));
        binVals.resize(nBins, 0.);
        histogram.resize(nBins, 0.);
        intersects.resize(Out->nhex, 0.);
        orResponse.resize(nOr, std::vector<FLT>(Out->nhex, 0.));
        orResponseSampled.resize(nOr, std::vector<FLT>(CHM.C.n, 0));

        max_X.resize(Out->nhex, -1e9);
        max_Y.resize(Out->nhex, -1e9);
    }

    void updateORresponses(void) {

        // OR ANALYSIS STEP 1. RECORD (MAX) RESPONSE TO ORIENTATIONS 0, 45, 90, 135.

        int n = Out->nhex;

        double phaseInc = M_PI / (double) nPhase;
        std::vector<float> maxPhase(n, 0.);

        std::vector<int> aff(2, 0);
        aff[1] = 1;
        std::vector<float> theta(nOr);

        int nFreq = Freq.size();

        for (unsigned int i = 0; i < nOr; i++) {
            theta[i] = i * M_PI / (double) nOr;
        }

        std::vector<float> maxOverAll(n, -1e9);
        for (unsigned int i = 0; i < nOr; i++) {
            std::fill(maxPhase.begin(), maxPhase.end(), -1e9);
            for (unsigned int j = 0; j < nPhase; j++) {
                double phase = j * phaseInc;
                for (unsigned int k = 0; k < nFreq; k++) {
                    float sf = Freq[k] / Net->spatialScale;
                    In->Grating(theta[i], phase, sf, 1.0); //drawing on the retina
                    Net->stepHidden(false);
                    Out->zero_X();
                    Out->step(aff);
                    for (int l = 0; l < n; l++) {
                        if (maxPhase[l] < Out->X[l]) {
                            maxPhase[l] = Out->X[l];
                        }
                        if (maxOverAll[l] < Out->X[l]) {
                            maxOverAll[l] = Out->X[l];
                            sfPref[l] = sf;
                        }
                    }
                }
            }
            orResponse[i] = maxPhase; //storing response

            // tmp copy of max (over phases) response on cortex
            for (int k = 0; k < n; k++) {
                Out->X[k] = maxPhase[k];
            }
            // subsample the response and store a copy
            CHM.step();
            std::vector<float> r(CHM.C.n);
            for (int k = 0; k < CHM.C.n; k++) {
                r[k] = CHM.C.vsquare[k].X;
            }
            orResponseSampled[i] = r;
        }

        float maxSF = -1e9;
        float minSF = +1e9;
        for (int l = 0; l < n; l++) {
            if (sfPref[l] < minSF) {
                minSF = sfPref[l];
            }
            if (sfPref[l] > maxSF) {
                maxSF = sfPref[l];
            }
        }
        std::cout << "Min SF = " << minSF << ", Max SF = " << maxSF << std::endl;

    }

    void updateORpreferences(void) {

        // ANALYSIS STEP 2. MEASURE ORIENTATION PREFERENCE & SELECTIVITY

        // Get orientation preference and selectivity
        int nOr = 4;
        std::vector<float> theta(nOr);
        for (unsigned int i = 0; i < nOr; i++) {
            theta[i] = i * M_PI / (double) nOr;
        }
        for (int i = 0; i < Out->nhex; i++) {
            float vx = 0.;
            float vy = 0.;
            for (int j = 0; j < nOr; j++) {
                vx += orResponse[j][i] * cos(2.0 * theta[j]);
                vy += orResponse[j][i] * sin(2.0 * theta[j]);
            }
            orPref[i] = 0.5 * (atan2(vy, vx) + M_PI);
            orSel[i] = pow(vy * vy + vx * vx, 0.5);
        }
    }

    void updateIsoORcontoursDiffs(void) {

        // ANALYSIS STEP 3. COMPUTE ISO-ORIENTATION CONTOURS (From Difference Images)

        // Diff of response to 0 and 90 degree stimuli (on original hex grid)
        std::vector<float> df1 = orResponse[0];
        for (int i = 0; i < df1.size(); i++) {
            df1[i] -= orResponse[2][i];
        }

        // Diff of response to 45 and 135 degree stimuli (on original hex grid)
        std::vector<float> df2 = orResponse[1];
        for (int i = 0; i < df2.size(); i++) {
            df2[i] -= orResponse[3][i];
        }

        // Get zero-crossings of the two response difference maps

        IsoORcontours[0] = shapeAn.get_contour_map_flag_nonorm(Out->hg, df1, 0.0, 1);
        IsoORcontours[1] = shapeAn.get_contour_map_flag_nonorm(Out->hg, df2, 0.0, 1);

    }

    void updateIsoORcontoursPrefs(void) {

        // ANALYSIS STEP 3. COMPUTE ISO-ORIENTATION CONTOURS (From Preferences)

        std::vector<float> real(Out->nhex, 0.);
        std::vector<float> imag(Out->nhex, 0.);
        for (int i = 0; i < Out->nhex; i++) {
            real[i] = cos(orPref[i] * 2.0);
            imag[i] = sin(orPref[i] * 2.0);
        }

        // Get zero-crossings of the two response difference maps
        IsoORcontours[0] = shapeAn.get_contour_map_flag_nonorm(Out->hg, real, 0.0, 1);
        IsoORcontours[1] = shapeAn.get_contour_map_flag_nonorm(Out->hg, imag, 0.0, 1);

    }

    void updateROIpinwheelCount(void) {

        // ANALYSIS STEP 4. COUNT PINWHEELS WITHIN ROI

        intersects = IsoORcontours[0];
        for (int k = 0; k < Out->nhex; k++) {
            intersects[k] *= IsoORcontours[1][k];
        }

        // remove neighbouring intersects (these are fractures)
        int countSpurious = 0;
        for (int i = 0; i < Out->nhex; i++) {
            if (intersects[i] == 1) {

                bool remSelf = false;

                if (Out->hg->d_ne[i] != -1) {
                    if (intersects[Out->hg->d_ne[i]] == 1) {
                        intersects[Out->hg->d_ne[i]] = 0;
                        countSpurious++;
                        remSelf = true;
                    }
                }
                if (Out->hg->d_nne[i] != -1) {
                    if (intersects[Out->hg->d_nne[i]] == 1) {
                        intersects[Out->hg->d_nne[i]] = 0;
                        countSpurious++;
                        remSelf = true;
                    }
                }
                if (Out->hg->d_nnw[i] != -1) {
                    if (intersects[Out->hg->d_nnw[i]] == 1) {
                        intersects[Out->hg->d_nnw[i]] = 0;
                        countSpurious++;
                        remSelf = true;
                    }
                }
                if (Out->hg->d_nw[i] != -1) {
                    if (intersects[Out->hg->d_nw[i]] == 1) {
                        intersects[Out->hg->d_nw[i]] = 0;
                        countSpurious++;
                        remSelf = true;
                    }
                }
                if (Out->hg->d_nsw[i] != -1) {
                    if (intersects[Out->hg->d_nsw[i]] == 1) {
                        intersects[Out->hg->d_nsw[i]] = 0;
                        countSpurious++;
                        remSelf = true;
                    }
                }
                if (Out->hg->d_nse[i] != -1) {
                    if (intersects[Out->hg->d_nse[i]] == 1) {
                        intersects[Out->hg->d_nse[i]] = 0;
                        countSpurious++;
                        remSelf = true;
                    }
                }

                if (remSelf) {
                    countSpurious++;
                    intersects[i] = 0;
                }
            }
        }

        std::cout << "Spurious crossings removed : " << countSpurious << std::endl;

        // count within ROI
        float halfWid = ROIwid * 0.5;
        int count = 0;
        for (int k = 0; k < Out->nhex; k++) {
            if ((fabs(Out->hg->vhexen[k]->x) < halfWid) && (fabs(Out->hg->vhexen[k]->y) < halfWid)) {
                if (intersects[k]) {
                    count++;
                }
            }
        }
        ROIpinwheelCount = (float) count;

    }


    std::vector<float> updateIsoORfrequencyEstimate(bool showfft) {

        // ANALYSIS STEP 5. ESTIMATE ISO-ORIENTATION COLUMN SPACING

        binVals.resize(nBins, 0.);
        histogram.resize(nBins, 0.);

        // Get frequency histogram from response to 0-90 degrees
        cv::Mat I1 = CHM.getDifferenceImage(orResponseSampled[0], orResponseSampled[2]);
        std::vector<std::vector<float> > h1 = CHM.fft(I1, nBins, gaussBlur, showfft);

        // Get frequency histogram from response to 45-135 degrees
        cv::Mat I2 = CHM.getDifferenceImage(orResponseSampled[1], orResponseSampled[3]);
        std::vector<std::vector<float> > h2 = CHM.fft(I2, nBins, gaussBlur, showfft);

        // add together two histograms (maybe should be done before combining?)
        binVals = h2[0];      // get histogram bin mid-values
        histogram = h1[1];
        for (int i = 0; i < nBins; i++) {
            histogram[i] += h2[1][i];
            histogram[i] *= 0.5;
        }

        // sample portion of histogram to fit
        int nsamp = nBins * sampleRange;
        arma::vec xs(nsamp);
        arma::vec ys(nsamp);
        for (int i = 0; i < nsamp; i++) {
            xs[i] = binVals[i];
            ys[i] = histogram[i];
        }

        // do polynomial fit
        arma::vec cf = arma::polyfit(xs, ys, polyOrder);

        // make a high-resolution model for the data
        int fitres = 1000;
        arma::vec xfit(fitres);
        for (int i = 0; i < fitres; i++) {
            xfit[i] = binVals[nsamp - 1] * (float) i / ((float) fitres - 1);
        }
        arma::vec yfit = arma::polyval(cf, xfit);

        // get frequency at which high-res model peaks
        float maxVal = -1e9;
        float maxX = 0;
        for (int i = 0; i < fitres; i++) {
            if (yfit[i] > maxVal) {
                maxVal = yfit[i];
                maxX = xfit[i];
            }
        }

        IsoORfrequency = maxX; // units are cycles / ROI-width
        IsoORcolumnSpacing = ROIwid /
                             IsoORfrequency;  // spacing between iso-orientation columns in units of cortex sheet, e.g., to plot scale bar on maps

        // return coeffs in standard vector
        std::vector<float> coeffs(cf.size());
        for (int i = 0; i < cf.size(); i++) {
            coeffs[i] = (float) cf[i];
        }

        return coeffs;

    }

    void updatePinwheelDensity(void) {

        // ANALYSIS STEP 6. CALCULATE PINWHEEL DENSITY
        pinwheelDensity = ROIpinwheelCount / (IsoORfrequency * IsoORfrequency);
    }

    void printPinwheelDensity(void) {
        std::cout << "Pinwheel density: " << pinwheelDensity << std::endl;
    }

    void printMetricInfo(void) {
        std::cout << "Peak frequency = q = " << IsoORfrequency << " cycles/ROI_width." << std::endl;
        std::cout << "Wavelength = 1/q = " << 1. / IsoORfrequency << " ROI_widths." << std::endl;
        std::cout << "Column spacing = lambda = wavelen. * ROI_width = " << IsoORcolumnSpacing
                  << " cortex unit distance." << std::endl;
        std::cout << "Pinwheel count: " << ROIpinwheelCount << std::endl;
        std::cout << "Pinwheel density: " << pinwheelDensity << std::endl << std::endl;
    }

    void updateRetinotopicMapping() {
        float radius = 0.05;
        //float sigma = RetinotopyAnalysisSigma;
        std::vector<int> aff(2, 0);
        aff[1] = 1; //setting which projections are normalized together
        int n = Out->nhex;

        //std::vector<std::vector<float>> retinotopy_data;

        int nX = 100;
        std::vector<float> xLocations(nX, 0.0);
        for (int i = 0; i < nX; i++) {
            xLocations[i] = (float) i / (float) nX - 0.5;
        }

        int nY = 100;
        std::vector<float> yLocations(nY, 0.0);
        for (int i = 0; i < nY; i++) {
            yLocations[i] = (float) i / (float) nY - 0.5;
        }

        std::vector<float> max_Response(n, -1e9);

        for (int ix = 0; ix < nX; ix++) {
            for (int iy = 0; iy < nY; iy++) {
                //std::vector<float> max_Response(n,-1e9);
                float x = xLocations[ix];
                float y = yLocations[iy];

                In->CircleStimulus(x, y, radius); // Place a stimulus on retina at location (x, y)
                //In->GaussianStimulus(x, y, sigma);
                Net->stepHidden(false);
                Out->zero_X();
                Out->step(aff); // Get v1 responses (Out->X)

                for (int l = 0; l < n; l++) {
                    if (max_Response[l] < Out->X[l]) {
                        // if [default value] < [response of l_th unit(s) on v1]
                        max_Response[l] = Out->X[l];
                        max_X[l] = x;
                        max_Y[l] = y;
                    }
                }
            }
        }
    }

    void v1ResponseWithCircularStimuli(void) {
        float radius = 0.15;
        std::vector<int> aff(2, 0);
        aff[1] = 1; //setting which projections are normalized together
        float x = 0;
        float y = 0;

        In->CircleStimulus(x, y, radius); // Place a stimulus on retina at location (x, y)
        Net->stepHidden(false);
        Out->zero_X();
        Out->step(aff); // Get v1 responses (Out->X)

    }

    void v1ResponseWithSquareStimuli(void) {
        float half_width = 0.15 / 2;
        std::vector<int> aff(2, 0);
        aff[1] = 1; //setting which projections are normalized together
        float x = 0;
        float y = 0;

        In->SquareStimulus(x, y, half_width); // Place a stimulus on retina at location (x, y)
        Net->stepHidden(false);
        Out->zero_X();
        Out->step(aff); // Get v1 responses (Out->X)

    }
};
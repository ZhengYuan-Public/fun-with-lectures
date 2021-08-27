#include "../../topo/gcal.h"
#include "../../topo/analysis.h"

#include <morph/HdfData.h>
#include <morph/Config.h>
#include <morph/ColourMap.h>
#include <morph/Visual.h>
#include <morph/HexGridVisual.h>
#include <morph/GraphVisual.h>
#include <morph/VisualDataModel.h>
#include <morph/Scale.h>
#include <morph/Vector.h>

typedef morph::VisualDataModel<FLT>* VdmPtr;

using morph::Config;
using morph::Tools;
using morph::ColourMap;

//
class gcalV2 : public gcal {
public:
    float decayConstant, oneMinusDecay, noiseScale;
    CortexSOM<FLT> CX2;
    alignas(alignof(std::vector<FLT>)) std::vector<FLT> CXcpy;
    alignas(alignof(std::vector<FLT>)) std::vector<FLT> CX2cpy;
    float v1v2strength, v2exstrength, v2instrength, v2v1exstrength, v2v1instrength;
    float v1v2alpha, v2exalpha, v2inalpha, v2v1exalpha, v2v1inalpha;
    gcalV2(morph::Config conf) : gcal(conf){

        decayConstant = conf.getFloat ("decayConstant", 0.3);
        noiseScale = conf.getFloat ("noiseScale", 0.02);
        oneMinusDecay = 1.0-decayConstant;

        v1v2strength = conf.getFloat ("v1v2strength", afferStrength);
        v2exstrength = conf.getFloat ("v2exstrength", excitStrength);
        v2instrength = conf.getFloat ("v2instrength", inhibStrength);

        v2v1exstrength = conf.getFloat ("v2v1exstrength", excitStrength);
        v2v1instrength = conf.getFloat ("v2v1instrength", inhibStrength);


        v1v2alpha = conf.getFloat ("v1v2alpha", afferAlpha);
        v2exalpha = conf.getFloat ("v2exalpha", excitAlpha);
        v2inalpha = conf.getFloat ("v2inalpha", inhibAlpha);

        v2v1exalpha = conf.getFloat ("v2v1exalpha", excitAlpha);
        v2v1inalpha = conf.getFloat ("v2v1inalpha", inhibAlpha);


        // CORTEX SHEET
        CX2.beta = beta;
        CX2.lambda = lambda;
        CX2.mu = mu;
        CX2.thetaInit = thetaInit;
        CX2.svgpath = conf.getString ("CX_svgpath", "boundaries/trialmod.svg");
        CX2.init();
        CX2.allocate();

        CX2.addProjection(CX.Xptr, CX.hg, afferRadius, v1v2strength, v1v2alpha, 0.0, false);
        CX2.addProjection(CX2.Xptr, CX2.hg, excitRadius, v2exstrength, v2exalpha, excitSigma, true);
        CX2.addProjection(CX2.Xptr, CX2.hg, inhibRadius, v2instrength, v2inalpha, inhibSigma, true);

        CX2.setNormalize(std::vector<int>(1,0));
        CX2.setNormalize(std::vector<int>(1,1));
        CX2.setNormalize(std::vector<int>(1,2));
        CX2.renormalize();


        // feedback projections
        CX.addProjection(CX2.Xptr, CX2.hg, excitRadius, v2v1exstrength, v2v1exalpha, excitSigma, true);
        CX.addProjection(CX2.Xptr, CX2.hg, inhibRadius, v2v1instrength, v2v1inalpha, inhibSigma, true);

        CX.setNormalize(std::vector<int>(1,4));
        CX.setNormalize(std::vector<int>(1,5));
        CX.renormalize();

        CX.zero_vector_variable (this->CXcpy);
        CX2.zero_vector_variable (this->CX2cpy);
    }

    void copyPrev(void){
        CXcpy = CX.X;
        CX2cpy = CX2.X;
    }

    void decayX(void){
        for(int i=0;i<CX.nhex;i++){
            CX.X[i] = CX.X[i]*decayConstant + oneMinusDecay*CXcpy[i] + (morph::Tools::randDouble()-0.5)*noiseScale; // should be norm distr.
        }
        for(int i=0;i<CX2.nhex;i++){
            CX2.X[i] = CX2.X[i]*decayConstant + oneMinusDecay*CX2cpy[i] + (morph::Tools::randDouble()-0.5)*noiseScale;
        }
    }

    void stepHidden(bool learning){
        LGN_ON.step();
        LGN_OFF.step();
    }


    void save(std::string filename){
        std::stringstream fname; fname << filename;
        morph::HdfData data(fname.str());
        std::vector<int> timetmp(1,time);
        data.add_contained_vals ("time", timetmp);
        for(unsigned int p=0;p<CX.Projections.size();p++){
            std::vector<FLT> proj = CX.Projections[p].getWeights();
            std::stringstream ss; ss<<"proj_1_"<<p;
            data.add_contained_vals (ss.str().c_str(), proj);
        }
        for(unsigned int p=0;p<CX2.Projections.size();p++){
            std::vector<FLT> proj = CX2.Projections[p].getWeights();
            std::stringstream ss; ss<<"proj_2_"<<p;
            data.add_contained_vals (ss.str().c_str(), proj);
        }
    }

    void load(std::string filename){
        std::stringstream fname; fname << filename;
        morph::HdfData data(fname.str(),1);
        std::vector<int> timetmp;
        data.read_contained_vals ("time", timetmp);
        time = timetmp[0];
        for(unsigned int p=0;p<CX.Projections.size();p++){
            std::vector<FLT> proj;
            std::stringstream ss; ss<<"proj_1_"<<p;
            data.read_contained_vals (ss.str().c_str(), proj);
            CX.Projections[p].setWeights(proj);
        }
        for(unsigned int p=0;p<CX2.Projections.size();p++){
            std::vector<FLT> proj;
            std::stringstream ss; ss<<"proj_2_"<<p;
            data.read_contained_vals (ss.str().c_str(), proj);
            CX.Projections[p].setWeights(proj);
        }
        std::cout<<"Loaded weights and modified time to " << time << std::endl;
    }


};










int main(int argc, char **argv){

    if (argc < 5) { std::cerr << "\nUsage: ./test configfile logdir seed mode intype weightfile(optional)\n\n"; return -1; }

    std::srand(std::stoi(argv[3]));       // set seed
    int INTYPE = std::stoi(argv[4]); // 0,1,2 Gaussian,Loaded,Camera input

    std::string paramsfile (argv[1]);
    Config conf(paramsfile);
    if (!conf.ready) { std::cerr << "Error setting up JSON config: " << conf.emsg << std::endl; return 1; }

    std::string logpath = argv[2];
    std::ofstream logfile;
    morph::Tools::createDir (logpath);
    { std::stringstream ss; ss << logpath << "/log.txt"; logfile.open(ss.str());}
    logfile<<"Hello."<<std::endl;

    unsigned int nBlocks = conf.getUInt ("blocks", 100);
    unsigned int steps = conf.getUInt("steps", 100);


    // Creates the network
    gcalV2 Net(conf);


    // Input specific setup
    switch(INTYPE){
        case(0):{ // Gaussian patterns
        } break;
        case(1):{   // preload patterns
            int ncols = conf.getUInt("patternSampleCols", 100);
            int nrows = conf.getUInt("patternSampleRows", 100);
            Net.HCM.initProjection(ncols,nrows,0.01,20.);
            std::string filename = conf.getString ("patterns", "configs/testPatterns.h5");
            Net.HCM.preloadPatterns(filename);
        } break;

        case(2):{
            int ncols = conf.getUInt("cameraCols", 100);
            int nrows = conf.getUInt("cameraRows", 100);
            int stepsize = conf.getUInt("cameraSampleStep", 7);
            int xoff = conf.getUInt("cameraOffsetX", 100);
            int yoff = conf.getUInt("cameraOffsetY", 0);
            Net.HCM.initProjection(ncols,nrows,0.01,20.);
            if(!Net.HCM.initCamera(xoff, yoff, stepsize)){ return 0;}
        } break;
    }

    if(argc>5){
        std::cout<<"Using weight file: "<<argv[5]<<std::endl;
        Net.load(argv[5]);
    } else {
        std::cout<<"Using random weights"<<std::endl;
    }




    // PRE-COMPUTE ANALYSIS STUFF
    CartHexSampler<FLT> CHM;
    std::vector<std::vector<float> > orResponse, orResponse2;
    std::vector<std::vector<float> > orResponseSampled, orResponseSampled2;

    float ROIwid, gratingWidth, amplitude;
    std::vector<float> orPref, orSel, intersects, binVals, histogram, sfPref, phPref, phSel, sfSel;
    std::vector<float> orPref2, orSel2, sfPref2, phPref2, phSel2, sfSel2;
    std::vector<int> maxOr, maxPhase, maxFreq, maxOr2, maxPhase2, maxFreq2;
    int nPhase, sampwid, nOr, nFreq, nPlotUnits;
    std::vector<float> Freq;
    std::vector<int> plotUnits;


    const Json::Value freq = conf.getArray ("freq");
    for (int i=0; i<freq.size(); i++) {
        Freq.push_back(freq[i].asFloat());
    }
    nFreq = Freq.size();

    const Json::Value plotunits = conf.getArray ("plotUnits");
    for (int i=0; i<plotunits.size(); i++) {
        plotUnits.push_back(plotunits[i].asUInt());
    }
    nPlotUnits = plotUnits.size();

    ROIwid = conf.getFloat ("ROIwid", 0.4);
    sampwid = conf.getUInt ("sampwid", 100);
    nPhase = conf.getUInt ("nPhase", 8);

    gratingWidth = 7.5;
    amplitude = 1.0;
    nOr = 4;

    CHM.initProjection(sampwid,sampwid,ROIwid,ROIwid,Net.CX.Xptr,Net.CX.hg, 0.05, 0.001); // maybe not CX hexgrid but retina??

    orPref.resize(Net.CX.nhex,0.);
    orSel.resize(Net.CX.nhex,0.);
    phPref.resize(Net.CX.nhex,0.);
    phSel.resize(Net.CX.nhex,0.);
    sfPref.resize(Net.CX.nhex,0.);
    sfSel.resize(Net.CX.nhex,0.);

    orPref2.resize(Net.CX2.nhex,0.);
    orSel2.resize(Net.CX2.nhex,0.);
    phPref2.resize(Net.CX2.nhex,0.);
    phSel2.resize(Net.CX2.nhex,0.);
    sfPref2.resize(Net.CX2.nhex,0.);
    sfSel2.resize(Net.CX2.nhex,0.);

    maxOr.resize(Net.CX.nhex,0);
    maxPhase.resize(Net.CX.nhex,0);
    maxFreq.resize(Net.CX.nhex,0);

    maxOr2.resize(Net.CX2.nhex,0);
    maxPhase2.resize(Net.CX2.nhex,0);
    maxFreq2.resize(Net.CX2.nhex,0);

    // setup storage for all responses (CX 1)
    std::vector<std::vector<std::vector<std::vector<float> > > > R;
    R.resize(Net.CX.nhex);
    for(unsigned int h=0;h<Net.CX.nhex;h++){
        R[h].resize(nOr);
        for(unsigned int i=0;i<nOr;i++){
            R[h][i].resize(nFreq);
            for(unsigned int j=0;j<nFreq;j++){
                R[h][i][j].resize(nPhase,0.);
            }
        }
    }

    // setup storage for all responses (CX 2)
    std::vector<std::vector<std::vector<std::vector<float> > > > R2;
    R2.resize(Net.CX2.nhex);
    for(unsigned int h=0;h<Net.CX2.nhex;h++){
        R2[h].resize(nOr);
        for(unsigned int i=0;i<nOr;i++){
            R2[h][i].resize(nFreq);
            for(unsigned int j=0;j<nFreq;j++){
                R2[h][i][j].resize(nPhase,0.);
            }
        }
    }

    // pre-compute stim values
    std::vector<float> theta(nOr);
    for(unsigned int i=0;i<nOr;i++){
        theta[i] = i*M_PI/(double)nOr;
    }

    std::vector<float> spatFreq(nFreq);
    for(unsigned int j=0;j<nFreq;j++){
        spatFreq[j] = Freq[j]/Net.spatialScale;
    }

    std::vector<float> phase(nPhase);
    for(unsigned int k=0;k<nPhase;k++){
        phase[k] = k*2.0*M_PI/(double)nPhase; // 2*M_PI? (DOUBLE-CHECK FOR GCAL MODEL TOO!)
    }




    // SETUP PLOTS
    std::chrono::steady_clock::time_point lastrender = std::chrono::steady_clock::now();

    //const unsigned int plotevery = conf.getUInt ("plotevery", 1);
    const bool saveplots = conf.getBool ("saveplots", false);
    unsigned int framecount = 0;
    const unsigned int win_height = conf.getUInt ("win_height", 400);
    const unsigned int win_width = conf.getUInt ("win_width", win_height);

    morph::Visual v1 (win_width, win_height, "model");
    v1.backgroundWhite();
    v1.sceneLocked = conf.getBool ("sceneLocked", false);
    v1.scenetrans_stepsize = 0.1;
    v1.fov = 50;

    // plotting grids
    std::vector<unsigned int> grids1(5);
    std::vector<unsigned int> grids2(8);
    float grid1offx = -1.25f;
    float grid1offy = -0.5f;
    float grid2offx = +0.8f;
    float txtoff = -0.55f;

    // ADD PLOTS TO SCENE

    // general purpose objects
    morph::Scale<FLT> zscale; zscale.setParams (0.0f, 0.0f);
    morph::Scale<FLT> cscale; cscale.do_autoscale = true;
    morph::ColourMap<FLT> hsv(morph::ColourMapType::Fixed);

    // Retina display
    morph::HexGridVisual<FLT> hgvRetina (v1.shaderprog,v1.tshaderprog, Net.IN.hg,std::array<float,3>{grid1offx+0.0f,grid1offy-0.9f,0.0f}, &(Net.IN.X),zscale,cscale,morph::ColourMapType::Inferno);
    grids1[0] = v1.addVisualModel (&hgvRetina);
    v1.getVisualModel (grids1[0])->addLabel ("retina", {-0.15f, txtoff, 0.0f},
                                             morph::colour::black, morph::VisualFont::VeraSerif, 0.1, 56);
    hgvRetina.hexVisMode = morph::HexVisMode::Triangles;

    // LGN ON display
    morph::HexGridVisual<FLT> hgvLGNon (v1.shaderprog,v1.tshaderprog, Net.LGN_ON.hg,std::array<float,3>{grid1offx-0.6f,grid1offy+0.0f,0.0f}, &(Net.LGN_ON.X),zscale,cscale,morph::ColourMapType::Inferno);
    grids1[1] = v1.addVisualModel (&hgvLGNon);
    v1.getVisualModel (grids1[1])->addLabel ("LGN on", {-0.2f, txtoff, 0.0f},
                                             morph::colour::black, morph::VisualFont::VeraSerif, 0.1, 56);
    hgvLGNon.hexVisMode = morph::HexVisMode::Triangles;

    // LGN OFF display
    morph::HexGridVisual<FLT> hgvLGNoff (v1.shaderprog,v1.tshaderprog, Net.LGN_OFF.hg,std::array<float,3>{grid1offx+0.6f,grid1offy+0.0f,0.0f}, &(Net.LGN_OFF.X),zscale,cscale,morph::ColourMapType::Inferno);
    grids1[2] = v1.addVisualModel (&hgvLGNoff);
    v1.getVisualModel (grids1[2])->addLabel ("LGN off", {-0.2f, txtoff, 0.0f},
                                             morph::colour::black, morph::VisualFont::VeraSerif, 0.1, 56);
    hgvLGNoff.hexVisMode = morph::HexVisMode::Triangles;

    // Cortex V1 display
    morph::HexGridVisual<FLT> hgvV1 (v1.shaderprog,v1.tshaderprog, Net.CX.hg,std::array<float,3>{grid1offx+0.0f,grid1offy+0.9f,0.0f}, &(Net.CX.X),zscale,cscale,morph::ColourMapType::Inferno);
    grids1[3] = v1.addVisualModel (&hgvV1);
    v1.getVisualModel (grids1[3])->addLabel ("L4", {-0.05f, txtoff, 0.0f},
                                             morph::colour::black, morph::VisualFont::VeraSerif, 0.1, 56);
    hgvV1.hexVisMode = morph::HexVisMode::Triangles;

    // Cortex V2 display
    morph::HexGridVisual<FLT> hgvV2 (v1.shaderprog,v1.tshaderprog, Net.CX2.hg,std::array<float,3>{grid1offx+0.0f,grid1offy+2.0f,0.0f}, &(Net.CX2.X),zscale,cscale,morph::ColourMapType::Inferno);
    grids1[4] = v1.addVisualModel (&hgvV2);
    v1.getVisualModel (grids1[4])->addLabel ("L2/3", {-0.05f, txtoff, 0.0f},
                                             morph::colour::black, morph::VisualFont::VeraSerif, 0.1, 56);
    hgvV2.hexVisMode = morph::HexVisMode::Triangles;




    // Cortex map orientation preference and selectivity display
    morph::HexGridVisualManual<FLT> hgvORPrefSel(v1.shaderprog,v1.tshaderprog, Net.CX.hg,morph::Vector<float,3>{grid2offx,1.1f,0.0f},&(orPref),zscale,cscale,morph::ColourMapType::Rainbow);
    grids2[0] = v1.addVisualModel (&hgvORPrefSel);
    v1.getVisualModel (grids2[0])->addLabel ("L4 OR pref*sel", {-0.05f, txtoff, 0.0f},
                                             morph::colour::black, morph::VisualFont::VeraSerif, 0.1, 56);
    hgvORPrefSel.hexVisMode = morph::HexVisMode::Triangles;


    // Cortex map phase preference and selectivity display
    morph::HexGridVisualManual<FLT> hgvPHPrefSel(v1.shaderprog,v1.tshaderprog, Net.CX.hg,morph::Vector<float,3>{grid2offx,0.0f,0.0f},&(phPref),zscale,cscale,morph::ColourMapType::Rainbow);
    grids2[1] = v1.addVisualModel (&hgvPHPrefSel);
    v1.getVisualModel (grids2[1])->addLabel ("L4 phase pref*sel", {-0.05f, txtoff, 0.0f},
                                             morph::colour::black, morph::VisualFont::VeraSerif, 0.1, 56);
    hgvORPrefSel.hexVisMode = morph::HexVisMode::Triangles;


    // Cortex map spatial frequency preference display
    morph::HexGridVisual<FLT> hgvSFpref(v1.shaderprog,v1.tshaderprog, Net.CX.hg,morph::Vector<float,3>{grid2offx,-1.1f,0.0f},&(sfPref),zscale,cscale,morph::ColourMapType::Jet);
    grids2[2] = v1.addVisualModel (&hgvSFpref);
    v1.getVisualModel (grids2[2])->addLabel ("L4 SF pref", {-0.05f, txtoff, 0.0f},
                                             morph::colour::black, morph::VisualFont::VeraSerif, 0.1, 56);
    hgvSFpref.hexVisMode = morph::HexVisMode::Triangles;


     // Graph responses by phase
     std::vector<float> graphX = phase;
     float wid = 0.7;
     float hei = 0.7;

    //  graph responses for CX
    std::vector<std::vector<float> > graphY(nPlotUnits,std::vector<float>(nPhase,0.));
     morph::GraphVisual<float>* gvPhaseResponse = new morph::GraphVisual<float> (v1.shaderprog, v1.tshaderprog, morph::Vector<float>{grid2offx-wid*0.5f,2.1f-hei*0.5f,0.0f});
     gvPhaseResponse->xlabel="phase";
     gvPhaseResponse->ylabel="response";
     gvPhaseResponse->setsize(wid,hei);
     gvPhaseResponse->setlimits (0,2*M_PI,0,1.0); // plot up to nyquist (pixels / 2)

     for(int i=0;i<nPlotUnits;i++){
        float colVal = (float)i/((float)nPlotUnits-1);
         morph::DatasetStyle ds;
         ds.linewidth = 0.01;
         ds.linecolour = {colVal, 0.0, 1-colVal};
         ds.markerstyle = morph::markerstyle::circle;
         ds.markersize = 0.00;
         ds.markercolour = {0.0, 0.0, 0.0};
         ds.markergap = 0.0;
         gvPhaseResponse->setdata (graphX, graphY[i], ds);
    }
     gvPhaseResponse->finalize();
     grids2[6] = v1.addVisualModel (static_cast<morph::VisualModel*>(gvPhaseResponse));

    //  graph responses for CX2
    std::vector<std::vector<float> > graphY2(nPlotUnits,std::vector<float>(nPhase,0.));
     morph::GraphVisual<float>* gvPhaseResponse2 = new morph::GraphVisual<float> (v1.shaderprog, v1.tshaderprog, morph::Vector<float>{grid2offx+1.1f-wid*0.5f,2.1f-hei*0.5f,0.0f});
     gvPhaseResponse2->xlabel="phase";
     gvPhaseResponse2->ylabel="response";
     gvPhaseResponse2->setsize(wid,hei);
     gvPhaseResponse2->setlimits (0,2*M_PI,0,1.0); // plot up to nyquist (pixels / 2)

     for(int i=0;i<nPlotUnits;i++){
         float colVal = (float)i/((float)nPlotUnits-1);
         morph::DatasetStyle ds;
         ds.linewidth = 0.01;
         ds.linecolour = {colVal, 0.0, 1-colVal};
         ds.markerstyle = morph::markerstyle::circle;
         ds.markersize = 0.00;
         ds.markercolour = {0.0, 0.0, 0.0};
         ds.markergap = 0.0;
         gvPhaseResponse2->setdata (graphX, graphY2[i], ds);
    }
     gvPhaseResponse2->finalize();
     grids2[7] = v1.addVisualModel (static_cast<morph::VisualModel*>(gvPhaseResponse2));


    // V2
    // Cortex map orientation preference and selectivity display
    morph::HexGridVisualManual<FLT> hgvORPrefSel2(v1.shaderprog,v1.tshaderprog, Net.CX2.hg,morph::Vector<float,3>{grid2offx+1.1f,1.1f,0.0f},&(orPref2),zscale,cscale,morph::ColourMapType::Rainbow);
    grids2[3] = v1.addVisualModel (&hgvORPrefSel2);
    v1.getVisualModel (grids2[3])->addLabel ("L2/3 OR pref*sel", {-0.05f, txtoff, 0.0f},
                                             morph::colour::black, morph::VisualFont::VeraSerif, 0.1, 56);
    hgvORPrefSel2.hexVisMode = morph::HexVisMode::Triangles;


    // Cortex map phase preference and selectivity display
    morph::HexGridVisualManual<FLT> hgvPHPrefSel2(v1.shaderprog,v1.tshaderprog, Net.CX.hg,morph::Vector<float,3>{grid2offx+1.1f,0.0f,0.0f},&(phPref2),zscale,cscale,morph::ColourMapType::Rainbow);
    grids2[4] = v1.addVisualModel (&hgvPHPrefSel2);
    v1.getVisualModel (grids2[4])->addLabel ("L2/3 phase pref*sel", {-0.05f, txtoff, 0.0f},
                                             morph::colour::black, morph::VisualFont::VeraSerif, 0.1, 56);
    hgvORPrefSel.hexVisMode = morph::HexVisMode::Triangles;


    // Cortex map spatial frequency preference display
    morph::HexGridVisual<FLT> hgvSFpref2(v1.shaderprog,v1.tshaderprog, Net.CX.hg,morph::Vector<float,3>{grid2offx+1.1f,-1.1f,0.0f},&(sfPref2),zscale,cscale,morph::ColourMapType::Jet);
    grids2[5] = v1.addVisualModel (&hgvSFpref2);
    v1.getVisualModel (grids2[5])->addLabel ("L2/3 SF pref", {-0.05f, txtoff, 0.0f},
                                             morph::colour::black, morph::VisualFont::VeraSerif, 0.1, 56);
    hgvSFpref2.hexVisMode = morph::HexVisMode::Triangles;




    float speed = conf.getFloat ("speed", 1.0);
    int duration = 20;
    int stimIncs = 15;
    // RUN THE MODEL
    int p =0;
    float direction = 0.0;
    for(int b=0;b<nBlocks;b++){

        // UPDATE MODEL
        for(unsigned int i=0;i<steps;i++){

            p = floor(morph::Tools::randDouble() * Net.HCM.PreLoadedPatterns.size());
            direction = morph::Tools::randDouble()*M_PI*2.0;

            for(unsigned int j=0;j<stimIncs;j++){

                float dur = (float)j/(float)stimIncs;
                float offx = (Net.HCM.patternsWid-Net.HCM.C.nx)*0.5;
                float offy = (Net.HCM.patternsWid-Net.HCM.C.ny)*0.5;

                offx += speed*dur*cos(direction);
                offy += speed*dur*sin(direction);

                Net.HCM.stepPreloaded(p,(int)offx, (int)offy);
                Net.IN.X = Net.HCM.X;
                Net.IN.clip_X(0.0,1.0);

                Net.LGN_ON.step();
                Net.LGN_OFF.step();

                for(unsigned int k=0;k<duration;k++){

                    Net.copyPrev();
                    Net.CX.step();
                    Net.CX2.step();
                    Net.decayX();

                    if(1){

                        { // afferent display
                            VdmPtr avm = (VdmPtr)v1.getVisualModel (grids1[0]);
                            avm->updateData (&(Net.IN.X));
                            avm->clearAutoscaleColour();
                        }

                        { // LGN_ON display
                            VdmPtr avm = (VdmPtr)v1.getVisualModel (grids1[1]);
                            avm->updateData (&(Net.LGN_ON.X));
                            avm->clearAutoscaleColour();
                        }

                        { // LGN_OFF display
                            VdmPtr avm = (VdmPtr)v1.getVisualModel (grids1[2]);
                            avm->updateData (&(Net.LGN_OFF.X));
                            avm->clearAutoscaleColour();
                        }

                        { // Cortex V1 display
                            VdmPtr avm = (VdmPtr)v1.getVisualModel (grids1[3]);
                            avm->updateData (&(Net.CX.X));
                            avm->clearAutoscaleColour();
                        }

                        { // Cortex V2 display
                            VdmPtr avm = (VdmPtr)v1.getVisualModel (grids1[4]);
                            avm->updateData (&(Net.CX2.X));
                            avm->clearAutoscaleColour();
                        }

                        std::chrono::steady_clock::duration sincerender = std::chrono::steady_clock::now() - lastrender;
                        if (std::chrono::duration_cast<std::chrono::milliseconds>(sincerender).count() > 17) {
                            glfwPollEvents();
                            v1.render();
                            lastrender = std::chrono::steady_clock::now();
                        }

                    }
                } // ending duration

                Net.CX.Projections[0].learn();
                Net.CX.Projections[1].learn();
                Net.CX.Projections[3].learn();
                Net.CX.renormalize();
                Net.CX2.Projections[2].learn();
                Net.CX2.renormalize();
                //Net.CX2.Projections[2].renormalize();
                //Net.CX.homeostasis();

                std::cout<<"step : "<<Net.time<<std::endl;
                std::cout<<"Retina -- "; Net.IN.printMinMax();
                std::cout<<"LGN ON -- "; Net.LGN_ON.printMinMax();
                std::cout<<"LGN OFF -- "; Net.LGN_OFF.printMinMax();
                std::cout<<"Cortex V1 -- "; Net.CX.printMinMax();
                std::cout<<"Cortex V2 -- "; Net.CX2.printMinMax();
                std::cout<<std::endl;

            }

            Net.IN.zero_X();
            Net.LGN_ON.step();
            Net.LGN_OFF.step();

            for(unsigned int k=0;k<duration;k++){
                Net.copyPrev();
                Net.CX.step();
                Net.CX2.step();
                Net.decayX();
            }

            Net.time++;

        }


        // COLLECT RESPONSES TO RANGE OF TEST STIMULI
        for(unsigned int i=0;i<nOr;i++){
            for(unsigned int j=0;j<nFreq;j++){
                for(unsigned int k=0;k<nPhase;k++){
                    Net.IN.Grating(theta[i],phase[k],spatFreq[j],1.0);
                    Net.LGN_ON.step();
                    Net.LGN_OFF.step();
                    Net.CX.zero_X();
                    Net.CX2.zero_X();
                    Net.CX.zero_vector_variable (Net.CXcpy);
                    Net.CX2.zero_vector_variable (Net.CX2cpy);
                    for(unsigned int inc=0;inc<20;inc++){
                        Net.copyPrev();
                        Net.CX.step();
                        Net.CX2.step();
                        Net.decayX();
                    }
                    // store responses (CX 1)
                    for(unsigned int h=0;h<Net.CX.nhex;h++){
                        R[h][i][j][k] = Net.CX.X[h];
                    }
                    // store responses (CX 2)
                    for(unsigned int h=0;h<Net.CX2.nhex;h++){
                        R2[h][i][j][k] = Net.CX2.X[h];
                    }
                }
            }
        }

        // CALCULATE PREFERENCES

        // GET INDICES OF BEST STIMULI
        for(unsigned int h=0;h<Net.CX.nhex;h++){
            float maxVar = -1e9;
            for(unsigned int i=0;i<nOr;i++){
                for(unsigned int j=0;j<nFreq;j++){
                    for(unsigned int k=0;k<nPhase;k++){
                        if(maxVar<R[h][i][j][k]){
                            maxVar=R[h][i][j][k];
                            maxOr[h] = i;
                            maxFreq[h] = j;
                            maxPhase[h] = k;
                        }
                    }
                }
            }
        }

        // GET INDICES OF BEST STIMULI
        for(unsigned int h=0;h<Net.CX2.nhex;h++){
            float maxVar = -1e9;
            for(unsigned int i=0;i<nOr;i++){
                for(unsigned int j=0;j<nFreq;j++){
                    for(unsigned int k=0;k<nPhase;k++){
                        if(maxVar<R2[h][i][j][k]){
                            maxVar=R2[h][i][j][k];
                            maxOr2[h] = i;
                            maxFreq2[h] = j;
                            maxPhase2[h] = k;
                        }
                    }
                }
            }
        }

        // COMPUTE ORIENTATION PREFERENCE AND SELECTIVITY (CX 1)
        for(unsigned int h=0;h<Net.CX.nhex;h++){
            std::vector<float> resp(nOr,0.);
            for(unsigned int i=0;i<nOr;i++){
                float maxVar = -1e9;
                for(unsigned int j=0;j<nFreq;j++){
                    for(unsigned int k=0;k<nPhase;k++){
                        if(maxVar<R[h][i][j][k]){
                            maxVar=R[h][i][j][k];
                        }
                    }
                }
                resp[i]=maxVar;
            }
            float vx = 0.;
            float vy = 0.;
            for(int i=0;i<nOr;i++){
                vx += resp[i] * cos(2.0*theta[i]);
                vy += resp[i] * sin(2.0*theta[i]);
            }
            orPref[h] = 0.5*(atan2(vy,vx)+M_PI);
            orSel[h] = pow(vy*vy+vx*vx,0.5);
            float sumResp = 0.;
            for(int i=0;i<nOr;i++){
                sumResp+=resp[i];
            }
            orSel[h] /= sumResp;
        }

        // COMPUTE ORIENTATION PREFERENCE AND SELECTIVITY (CX 2)
        for(unsigned int h=0;h<Net.CX2.nhex;h++){
            std::vector<float> resp(nOr,0.);
            for(unsigned int i=0;i<nOr;i++){
                float maxVar = -1e9;
                for(unsigned int j=0;j<nFreq;j++){
                    for(unsigned int k=0;k<nPhase;k++){
                        if(maxVar<R2[h][i][j][k]){
                            maxVar=R2[h][i][j][k];
                        }
                    }
                }
                resp[i]=maxVar;
            }
            float vx = 0.;
            float vy = 0.;
            for(int i=0;i<nOr;i++){
                vx += resp[i] * cos(2.0*theta[i]);
                vy += resp[i] * sin(2.0*theta[i]);
            }
            orPref2[h] = 0.5*(atan2(vy,vx)+M_PI);
            orSel2[h] = pow(vy*vy+vx*vx,0.5);
            float sumResp = 0.;
            for(int i=0;i<nOr;i++){
                sumResp+=resp[i];
            }
            orSel2[h] /= sumResp;
        }

         // COMPUTE PHASE PREFERENCE AND SELECTIVITY (CX 1)
        for(unsigned int h=0;h<Net.CX.nhex;h++){
            std::vector<float> resp(nPhase,0.);
            for(unsigned int k=0;k<nPhase;k++){
                float maxVar = -1e9;
                for(unsigned int j=0;j<nFreq;j++){
                    for(unsigned int i=0;i<nOr;i++){
                        if(maxVar<R[h][i][j][k]){
                            maxVar=R[h][i][j][k];
                        }
                    }
                }
                resp[k]=maxVar;
            }
            float vx = 0.;
            float vy = 0.;
            for(int k=0;k<nPhase;k++){
                vx += resp[k] * cos(phase[k]);
                vy += resp[k] * sin(phase[k]);
            }
            phPref[h] = (atan2(vy,vx)+M_PI);
            phSel[h] = pow(vy*vy+vx*vx,0.5);
            float sumResp = 0.;
            for(int k=0;k<nPhase;k++){
                sumResp+=resp[k];
            }
            phSel[h] /= sumResp;
        }

        // COMPUTE PHASE PREFERENCE AND SELECTIVITY (CX 2)
        for(unsigned int h=0;h<Net.CX2.nhex;h++){
            std::vector<float> resp(nPhase,0.);
            for(unsigned int k=0;k<nPhase;k++){
                float maxVar = -1e9;
                for(unsigned int j=0;j<nFreq;j++){
                    for(unsigned int i=0;i<nOr;i++){
                        if(maxVar<R2[h][i][j][k]){
                            maxVar=R2[h][i][j][k];
                        }
                    }
                }
                resp[k]=maxVar;
            }
            float vx = 0.;
            float vy = 0.;
            for(int k=0;k<nPhase;k++){
                vx += resp[k] * cos(phase[k]);
                vy += resp[k] * sin(phase[k]);
            }
            phPref2[h] = (atan2(vy,vx)+M_PI);
            phSel2[h] = pow(vy*vy+vx*vx,0.5);
            float sumResp = 0.;
            for(int k=0;k<nPhase;k++){
                sumResp+=resp[k];
            }
            phSel2[h] /= sumResp;
        }

        // COMPUTE SF PREFERENCE AND SELECTIVITY (CX 1)
        for(unsigned int h=0;h<Net.CX.nhex;h++){
            std::vector<float> resp(nFreq,0.);
            for(unsigned int j=0;j<nFreq;j++){
                float maxVar = -1e9;
                for(unsigned int i=0;i<nOr;i++){
                    for(unsigned int k=0;k<nPhase;k++){
                        if(maxVar<R[h][i][j][k]){
                            maxVar=R[h][i][j][k];
                        }
                    }
                }
                resp[j]=maxVar;
            }
            float maxVar = -1e9;
            int maxInd = 0;
            float sumVar = 0.;
            for(int j=0;j<nFreq;j++){
                if(maxVar<resp[j]){
                    maxVar=resp[j];
                    maxInd = j;
                }
                sumVar += resp[j];
            }
            sfPref[h] = spatFreq[maxInd];
            sfSel[h] = maxVar / sumVar;
        }


        // COMPUTE SF PREFERENCE AND SELECTIVITY (CX 2)
        for(unsigned int h=0;h<Net.CX2.nhex;h++){
            std::vector<float> resp(nFreq,0.);
            for(unsigned int j=0;j<nFreq;j++){
                float maxVar = -1e9;
                for(unsigned int i=0;i<nOr;i++){
                    for(unsigned int k=0;k<nPhase;k++){
                        if(maxVar<R2[h][i][j][k]){
                            maxVar=R2[h][i][j][k];
                        }
                    }
                }
                resp[j]=maxVar;
            }
            float maxVar = -1e9;
            int maxInd = 0;
            float sumVar = 0.;
            for(int j=0;j<nFreq;j++){
                if(maxVar<resp[j]){
                    maxVar=resp[j];
                    maxInd = j;
                }
                sumVar += resp[j];
            }
            sfPref2[h] = spatFreq[maxInd];
            sfSel2[h] = maxVar / sumVar;
        }



        // UPDATE MAP DISPLAYS
        { // V1 OR pref display
            float maxSel = -1e9;
            float minSel = +1e9;
            for(int i=0;i<Net.CX.nhex;i++){
                if(maxSel<orSel[i]){ maxSel=orSel[i];}
                if(minSel>orSel[i]){ minSel=orSel[i];}
            }
            float rangeSel = 1./(maxSel-minSel);
            float overPi = 1./M_PI;

            for(int i=0;i<Net.CX.nhex;i++){
                float pref = orPref[i]*overPi;
                float sel = (orSel[i]-minSel)*rangeSel;
                std::array<float, 3> rgb2 = hsv.hsv2rgb(pref,1.0,sel);
                hgvORPrefSel.R[i] = rgb2[0];
                hgvORPrefSel.G[i] = rgb2[1];
                hgvORPrefSel.B[i] = rgb2[2];
            }
        }

        { // V2 OR pref display
            float maxSel = -1e9;
            float minSel = +1e9;
            for(int i=0;i<Net.CX2.nhex;i++){
                if(maxSel<orSel2[i]){ maxSel=orSel2[i];}
                if(minSel>orSel2[i]){ minSel=orSel2[i];}
            }
            float rangeSel = 1./(maxSel-minSel);
            float overPi = 1./M_PI;

            for(int i=0;i<Net.CX2.nhex;i++){
                float pref = orPref2[i]*overPi;
                float sel = (orSel2[i]-minSel)*rangeSel;
                std::array<float, 3> rgb2 = hsv.hsv2rgb(pref,1.0,sel);
                hgvORPrefSel2.R[i] = rgb2[0];
                hgvORPrefSel2.G[i] = rgb2[1];
                hgvORPrefSel2.B[i] = rgb2[2];
            }
        }

        { // V1 PH pref display
            float maxSel = -1e9;
            float minSel = +1e9;
            for(int i=0;i<Net.CX.nhex;i++){
                if(maxSel<phSel[i]){ maxSel=phSel[i];}
                if(minSel>phSel[i]){ minSel=phSel[i];}
            }
            float rangeSel = 1./(maxSel-minSel);
            float over2Pi = 0.5/M_PI;

            for(int i=0;i<Net.CX.nhex;i++){
                float pref = phPref[i]*over2Pi;
                float sel = (phSel[i]-minSel)*rangeSel;
                std::array<float, 3> rgb2 = hsv.hsv2rgb(pref,1.0,sel);
                hgvPHPrefSel.R[i] = rgb2[0];
                hgvPHPrefSel.G[i] = rgb2[1];
                hgvPHPrefSel.B[i] = rgb2[2];
            }
        }

        { // V2 PH pref display
            float maxSel = -1e9;
            float minSel = +1e9;
            for(int i=0;i<Net.CX2.nhex;i++){
                if(maxSel<phSel2[i]){ maxSel=phSel2[i];}
                if(minSel>phSel2[i]){ minSel=phSel2[i];}
            }
            float rangeSel = 1./(maxSel-minSel);
            float over2Pi = 0.5/M_PI;

            for(int i=0;i<Net.CX2.nhex;i++){
                float pref = phPref2[i]*over2Pi;
                float sel = (phSel2[i]-minSel)*rangeSel;
                std::array<float, 3> rgb2 = hsv.hsv2rgb(pref,1.0,sel);
                hgvPHPrefSel2.R[i] = rgb2[0];
                hgvPHPrefSel2.G[i] = rgb2[1];
                hgvPHPrefSel2.B[i] = rgb2[2];
            }
        }



        { // V1 OR preference map
            VdmPtr avm = (VdmPtr)v1.getVisualModel (grids2[0]);
            avm->updateData (&orPref);
            avm->clearAutoscaleColour();
        }

        { // V1 OR preference map
            VdmPtr avm = (VdmPtr)v1.getVisualModel (grids2[1]);
            avm->updateData (&phPref);
            avm->clearAutoscaleColour();
        }

        { // V1 Spatial Frequency Preference map
            VdmPtr avm = (VdmPtr)v1.getVisualModel (grids2[2]);
            avm->updateData (&sfPref);
            avm->clearAutoscaleColour();
        }


         { // Update histogram display
            for(int i=0;i<nPlotUnits;i++){
                for(int j=0;j<nPhase;j++){
                    int h = plotUnits[i];
                    graphY[i][j] = R[h][maxOr[h]][maxFreq[h]][j];
                }
                gvPhaseResponse->update (graphX, graphY[i], i);
            }
         }


         { // Update histogram display
            for(int i=0;i<nPlotUnits;i++){
                for(int j=0;j<nPhase;j++){
                    int h = plotUnits[i];
                    graphY2[i][j] = R2[h][maxOr2[h]][maxFreq2[h]][j];
                }
                gvPhaseResponse2->update (graphX, graphY2[i], i);
            }
         }


        { // V2 OR preference map

            VdmPtr avm = (VdmPtr)v1.getVisualModel (grids2[3]);
            avm->updateData (&orPref2);
            avm->clearAutoscaleColour();
        }

        { // V2 OR preference map

            VdmPtr avm = (VdmPtr)v1.getVisualModel (grids2[4]);
            avm->updateData (&phPref2);
            avm->clearAutoscaleColour();
        }

        { // V2 Spatial Frequency Preference map

            VdmPtr avm = (VdmPtr)v1.getVisualModel (grids2[5]);
            avm->updateData (&sfPref2);
            avm->clearAutoscaleColour();
        }

        std::chrono::steady_clock::duration sincerender = std::chrono::steady_clock::now() - lastrender;
        if (std::chrono::duration_cast<std::chrono::milliseconds>(sincerender).count() > 17) {
            glfwPollEvents();
            v1.render();
            lastrender = std::chrono::steady_clock::now();
        }

        if(saveplots){
            savePngs (logpath, "model", framecount, v1);
            framecount++;
        }

        // SAVE NETWORK WEIGHTS
        std::stringstream ss; ss << logpath << "/weights.h5";
        logfile<<"Weights saved at time: "<<Net.time<<std::endl;
        Net.save(ss.str());

    }


    return 0.;
}
















/*
 class orientationPinwheelDensity2 {

 public:

 CartHexSampler<FLT> CHM;
 std::vector<std::vector<FLT> > orResponse, orResponse2;
 std::vector<std::vector<FLT> > orResponseSampled, orResponseSampled2;
 morph::ShapeAnalysis<FLT> shapeAn;

 float ROIwid, ROIpinwheelCount, IsoORfrequency, IsoORcolumnSpacing, sampleRange, gratingWidth, amplitude;
 std::vector<FLT> orPref, orSel, intersects, binVals, histogram, sfPref, phPref, phSel;
 std::vector<FLT> orPref2, orSel2, sfPref2, phPref2, phSel2;
 std::vector<std::vector<FLT> > IsoORcontours;
 int nBins, nPhase, polyOrder, sampwid, gaussBlur, nOr, nFreq;
 float pinwheelDensity;
 std::vector<float> Freq;

 PatternGenerator_Sheet<FLT>* In;
 CortexSOM<FLT>* Out;
 CortexSOM<FLT>* Out2;
 gcalV2* Net;

 orientationPinwheelDensity2 (gcalV2* Net) {

 this->Net = Net;
 this-> In = &Net->IN;
 this-> Out = &Net->CX;
 this-> Out2 = &Net->CX2;

 polyOrder = Net->conf.getUInt ("polyOrder", 4);
 sampleRange = Net->conf.getFloat ("sampleRange", 1.0);
 ROIwid = Net->conf.getFloat ("ROIwid", 0.4);
 sampwid = Net->conf.getUInt ("sampwid", 100);
 gaussBlur = Net->conf.getUInt ("gaussBlur", 1);
 nBins = Net->conf.getUInt ("nBins", 50);
 nPhase = Net->conf.getUInt ("nPhase", 8);


 const Json::Value freq = Net->conf.getArray ("freq");
 for (int i=0; i<freq.size(); i++) {
 Freq.push_back(freq[i].asFloat());
 }
 nFreq = Freq.size();

 gratingWidth = 7.5;
 amplitude = 1.0;
 nOr = 4;

 CHM.initProjection(sampwid,sampwid,ROIwid,ROIwid,Out->Xptr,Out->hg, 0.05, 0.001);

 orPref.resize(Out->nhex,0.);
 orSel.resize(Out->nhex,0.);
 phPref.resize(Out->nhex,0.);
 phSel.resize(Out->nhex,0.);
 sfPref.resize(Out->nhex,0.);

 orPref2.resize(Out2->nhex,0.);
 orSel2.resize(Out2->nhex,0.);
 phPref2.resize(Out2->nhex,0.);
 phSel2.resize(Out2->nhex,0.);
 sfPref2.resize(Out2->nhex,0.);


 IsoORcontours.resize(2,std::vector<FLT>(Out->nhex,0.));
 binVals.resize(nBins,0.);
 histogram.resize(nBins,0.);
 intersects.resize(Out->nhex,0.);
 orResponse.resize(nOr,std::vector<FLT>(Out->nhex,0.));
 orResponseSampled.resize(nOr,std::vector<FLT>(CHM.C.n,0));

 orResponse2.resize(nOr,std::vector<FLT>(Out->nhex,0.));
 orResponseSampled2.resize(nOr,std::vector<FLT>(CHM.C.n,0));

 }

 void updateORresponses(void){

 // OR ANALYSIS STEP 1. RECORD (MAX) RESPONSE TO ORIENTATIONS 0, 45, 90, 135.
 int n = Out->nhex;
 int n2 = Out2->nhex;
 double phaseInc = M_PI/(double)nPhase;
 std::vector<float> maxPhase(n,0.);
 std::vector<float> maxPhase2(n2,0.);
 //std::vector<int> aff(2,0); aff[1]=1;
 std::vector<float> theta(nOr);
 int nFreq = Freq.size();
 for(unsigned int i=0;i<nOr;i++){
 theta[i] = i*M_PI/(double)nOr;
 }
 std::vector<float> maxOverAll(n,-1e9);
 std::vector<float> maxOverAll2(n2,-1e9);
 for(unsigned int i=0;i<nOr;i++){
 std::fill(maxPhase.begin(),maxPhase.end(),-1e9);
 std::fill(maxPhase2.begin(),maxPhase2.end(),-1e9);
 for(unsigned int j=0;j<nPhase;j++){

 double phase = j*phaseInc;

 for(unsigned int k=0;k<nFreq;k++){
 float sf = Freq[k]/Net->spatialScale;
 In->Grating(theta[i],phase,sf,1.0);
 Net->LGN_ON.step();
 Net->LGN_OFF.step();
 Net->CX.zero_X();
 Net->CX2.zero_X();
 Net->CX.zero_vector_variable (Net->CXcpy);
 Net->CX2.zero_vector_variable (Net->CX2cpy);
 for(unsigned int inc=0;inc<20;inc++){
 Net->copyPrev();
 Net->CX.step();
 Net->CX2.step();
 Net->decayX();
 }


 // CX
 for(int l=0;l<n;l++){
 if(maxPhase[l]<Out->X[l]){
 maxPhase[l] = Out->X[l];
 }
 if(maxOverAll[l]<Out->X[l]){
 maxOverAll[l] = Out->X[l];
 sfPref[l] = sf; // SF pref here is the SF of the best stim
 }
 }

 // CX2
 for(int l=0;l<n2;l++){
 if(maxPhase2[l]<Out2->X[l]){
 maxPhase2[l] = Out2->X[l];
 }
 if(maxOverAll2[l]<Out2->X[l]){
 maxOverAll2[l] = Out2->X[l];
 sfPref2[l] = sf;
 }
 }
 }
 }
 orResponse[i]=maxPhase;
 orResponse2[i]=maxPhase2;

 // tmp copy of max (over phases) response on cortex
 for(int k=0;k<n;k++){
 Out->X[k] = maxPhase[k];
 }

 for(int k=0;k<n2;k++){
 Out2->X[k] = maxPhase2[k];
 }

 // MAYBE THE REST JUST REQUIRED FOR PINWHEEL DENSITY +?
 // subsample the response and store a copy
 CHM.step();
 std::vector<float> r(CHM.C.n);
 for(int k=0;k<CHM.C.n;k++){
 r[k] = CHM.C.vsquare[k].X;
 }
 orResponseSampled[i] = r;

 }

 { // CX
 float maxSF = -1e9;
 float minSF = +1e9;
 for(int l=0;l<n;l++){
 if(sfPref[l]<minSF){
 minSF = sfPref[l];
 }
 if(sfPref[l]>maxSF){
 maxSF = sfPref[l];
 }
 }
 }

 { // CX2
 float maxSF = -1e9;
 float minSF = +1e9;
 for(int l=0;l<n;l++){
 if(sfPref2[l]<minSF){
 minSF = sfPref2[l];
 }
 if(sfPref2[l]>maxSF){
 maxSF = sfPref2[l];
 }
 }
 }

 }





 void updateORpreferences(void){

 // ANALYSIS STEP 2. MEASURE ORIENTATION PREFERENCE & SELECTIVITY

 // Get orientation preference and selectivity
 int nOr = 4;
 std::vector<float> theta(nOr);
 for(unsigned int i=0;i<nOr;i++){
 theta[i] = i*M_PI/(double)nOr;
 }
 for(int i=0;i<Out->nhex;i++){
 float vx = 0.;
 float vy = 0.;
 for(int j=0;j<nOr;j++){
 vx += orResponse[j][i] * cos(2.0*theta[j]);
 vy += orResponse[j][i] * sin(2.0*theta[j]);
 }
 orPref[i] = 0.5*(atan2(vy,vx)+M_PI);
 orSel[i] = pow(vy*vy+vx*vx,0.5);
 }


 // REPEAT FOR CX2 RESPONSES
 for(int i=0;i<Out2->nhex;i++){
 float vx = 0.;
 float vy = 0.;
 for(int j=0;j<nOr;j++){
 vx += orResponse2[j][i] * cos(2.0*theta[j]);
 vy += orResponse2[j][i] * sin(2.0*theta[j]);
 }
 orPref2[i] = 0.5*(atan2(vy,vx)+M_PI);
 orSel2[i] = pow(vy*vy+vx*vx,0.5);
 }


 }


 void updateIsoORcontoursDiffs(void){

 // ANALYSIS STEP 3. COMPUTE ISO-ORIENTATION CONTOURS (From Difference Images)

 // Diff of response to 0 and 90 degree stimuli (on original hex grid)
 std::vector<float> df1 = orResponse[0];
 for(int i=0;i<df1.size();i++){
 df1[i] -= orResponse[2][i];
 }

 // Diff of response to 45 and 135 degree stimuli (on original hex grid)
 std::vector<float> df2 = orResponse[1];
 for(int i=0;i<df2.size();i++){
 df2[i] -= orResponse[3][i];
 }

 // Get zero-crossings of the two response difference maps

 IsoORcontours[0] = shapeAn.get_contour_map_flag_nonorm(Out->hg, df1, 0.0, 1);
 IsoORcontours[1] = shapeAn.get_contour_map_flag_nonorm(Out->hg, df2, 0.0, 1);

 }

 void updateIsoORcontoursPrefs(void){

 // ANALYSIS STEP 3. COMPUTE ISO-ORIENTATION CONTOURS (From Preferences)

 std::vector<float> real(Out->nhex,0.);
 std::vector<float> imag(Out->nhex,0.);
 for(int i=0;i<Out->nhex;i++){
 real[i] = cos(orPref[i]*2.0);
 imag[i] = sin(orPref[i]*2.0);
 }

 // Get zero-crossings of the two response difference maps
 IsoORcontours[0] = shapeAn.get_contour_map_flag_nonorm(Out->hg, real, 0.0, 1);
 IsoORcontours[1] = shapeAn.get_contour_map_flag_nonorm(Out->hg, imag, 0.0, 1);

 }

 void updateROIpinwheelCount(void){

 // ANALYSIS STEP 4. COUNT PINWHEELS WITHIN ROI

 intersects = IsoORcontours[0];
 for(int k=0;k<Out->nhex;k++){
 intersects[k] *= IsoORcontours[1][k];
 }

 // remove neighbouring intersects (these are fractures)
 int countSpurious = 0;
 for(int i=0;i<Out->nhex;i++){
 if(intersects[i]==1){

 bool remSelf = false;

 if(Out->hg->d_ne[i] !=-1){ if(intersects[Out->hg->d_ne[i]] ==1){ intersects[Out->hg->d_ne[i]] =0; countSpurious++; remSelf = true;} }
 if(Out->hg->d_nne[i]!=-1){ if(intersects[Out->hg->d_nne[i]]==1){ intersects[Out->hg->d_nne[i]]=0; countSpurious++; remSelf = true;} }
 if(Out->hg->d_nnw[i]!=-1){ if(intersects[Out->hg->d_nnw[i]]==1){ intersects[Out->hg->d_nnw[i]]=0; countSpurious++; remSelf = true;} }
 if(Out->hg->d_nw[i] !=-1){ if(intersects[Out->hg->d_nw[i]] ==1){ intersects[Out->hg->d_nw[i]] =0; countSpurious++; remSelf = true;} }
 if(Out->hg->d_nsw[i]!=-1){ if(intersects[Out->hg->d_nsw[i]]==1){ intersects[Out->hg->d_nsw[i]]=0; countSpurious++; remSelf = true;} }
 if(Out->hg->d_nse[i]!=-1){ if(intersects[Out->hg->d_nse[i]]==1){ intersects[Out->hg->d_nse[i]]=0; countSpurious++; remSelf = true;} }

 if (remSelf) { countSpurious++; intersects[i] = 0; }
 }
 }

 std::cout<<"Spurious crossings removed : "<<countSpurious<<std::endl;

 // count within ROI
 float halfWid = ROIwid*0.5;
 int count=0;
 for(int k=0;k<Out->nhex;k++){
 if((fabs(Out->hg->vhexen[k]->x)<halfWid)&&(fabs(Out->hg->vhexen[k]->y)<halfWid)){
 if(intersects[k]){
 count++;
 }
 }
 }
 ROIpinwheelCount = (float)count;

 }


 std::vector<float> updateIsoORfrequencyEstimate(bool showfft){

 // ANALYSIS STEP 5. ESTIMATE ISO-ORIENTATION COLUMN SPACING

 binVals.resize(nBins,0.);
 histogram.resize(nBins,0.);

 // Get frequency histogram from response to 0-90 degrees
 cv::Mat I1 = CHM.getDifferenceImage(orResponseSampled[0],orResponseSampled[2]);
 std::vector<std::vector<float> > h1 = CHM.fft(I1, nBins, gaussBlur, showfft);

 // Get frequency histogram from response to 45-135 degrees
 cv::Mat I2 = CHM.getDifferenceImage(orResponseSampled[1],orResponseSampled[3]);
 std::vector<std::vector<float> > h2 = CHM.fft(I2, nBins, gaussBlur, showfft);

 // add together two histograms (maybe should be done before combining?)
 binVals = h2[0];      // get histogram bin mid-values
 histogram = h1[1];
 for(int i=0;i<nBins;i++){
 histogram[i] += h2[1][i];
 histogram[i] *= 0.5;
 }

 // sample portion of histogram to fit
 int nsamp = nBins*sampleRange;
 arma::vec xs(nsamp);
 arma::vec ys(nsamp);
 for(int i=0;i<nsamp;i++){
 xs[i] = binVals[i];
 ys[i] = histogram[i];
 }

 // do polynomial fit
 arma::vec cf = arma::polyfit(xs,ys,polyOrder);

 // make a high-resolution model for the data
 int fitres = 1000;
 arma::vec xfit(fitres);
 for(int i=0;i<fitres;i++){
 xfit[i] = binVals[nsamp-1]*(float)i/((float)fitres-1);
 }
 arma::vec yfit = arma::polyval(cf,xfit);

 // get frequency at which high-res model peaks
 float maxVal = -1e9;
 float maxX = 0;
 for(int i=0;i<fitres;i++){
 if(yfit[i]>maxVal){
 maxVal = yfit[i];
 maxX = xfit[i];
 }
 }

 IsoORfrequency = maxX; // units are cycles / ROI-width
 IsoORcolumnSpacing = ROIwid / IsoORfrequency;  // spacing between iso-orientation columns in units of cortex sheet, e.g., to plot scale bar on maps

 // return coeffs in standard vector
 std::vector<float> coeffs(cf.size());
 for(int i=0;i<cf.size();i++){
 coeffs[i] = (float)cf[i];
 }

 return coeffs;

 }

 void updatePinwheelDensity(void){

 // ANALYSIS STEP 6. CALCULATE PINWHEEL DENSITY
 pinwheelDensity = ROIpinwheelCount / (IsoORfrequency*IsoORfrequency);
 }

 void printPinwheelDensity(void){
 std::cout<<"Pinwheel density: "<<pinwheelDensity<<std::endl;
 }

 void printMetricInfo(void){
 std::cout<<"Peak frequency = q = "<<IsoORfrequency<<" cycles/ROI_width."<<std::endl;
 std::cout<<"Wavelength = 1/q = "<<1./IsoORfrequency<<" ROI_widths."<<std::endl;
 std::cout<<"Column spacing = lambda = wavelen. * ROI_width = "<<IsoORcolumnSpacing<<" cortex unit distance."<<std::endl;
 std::cout<<"Pinwheel count: "<<ROIpinwheelCount<<std::endl;
 std::cout<<"Pinwheel density: "<<pinwheelDensity<< std::endl<<std::endl;
 }


 };


 */

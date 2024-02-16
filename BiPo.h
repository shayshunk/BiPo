#include <array>
#include <fstream>
#include <iostream>
#include <string>

#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TLeaf.h"
#include "TTree.h"

// Invariables

// Variables

// Print flags
bool DETECTOR_VERBOSITY = 1;

// Utilities for parameters

enum Directions
{
    X = 0,
    Y,
    Z,
    DirectionSize
};

std::string AxisToString(int num)
{
    std::string name;

    switch (num)
    {
        case 0:
            name = "X";
            break;
        case 1:
            name = "Y";
            break;
        case 2:
            name = "Z";
            break;
        default:
            name = "X";
    }

    return name;
}

enum Signals
{
    Correlated = 0,
    Accidental,
    TotalDifference,
    SignalSize
};

std::string SignalToString(int num)
{
    std::string name;

    switch (num)
    {
        case 0:
            name = "Correlated";
            break;
        case 1:
            name = "Accidental";
            break;
        case 2:
            name = "Total Difference";
            break;
        default:
            name = "Total Difference";
    }

    return name;
}

enum Datasets
{
    Data = 0,
    DataUnbiased,
    DatasetSize
};

std::string DatasetToString(int num)
{
    std::string name;

    switch (num)
    {
        case 0:
            name = "Data";
            break;
        case 1:
            name = "Data Unbiased";
            break;
        default:
            name = "Data";
    }

    return name;
}

class BiPo
{
  public:
    // Variables
    double livetimeOff = 0, livetimeOn = 0;

    // Main functions
    BiPo();
    void ReadFileList();
    void SetUpHistograms();
    void SetBranchAddresses(std::shared_ptr<TTree> rootTree);
    void FillHistogram();
    void FillHistogramUnbiased(int signalSet);
    void CalculateUnbiasing();
    void SubtractBackgrounds();
    void FillOutputFile();

    // Inline functions
    inline void ResetLineNumber() { lineNumber = 0; }
    inline void ResetLineCounter() { lineCounter = 0; }
    inline void ResetIndex() { index = 0; }

  private:
    // Histogram to count IBDs
    std::array<std::array<std::array<TH1D, DirectionSize>, SignalSize>, DatasetSize> histogram;

    // File list
    std::array<std::string, 1740> files;

    // Values grabbed from ROOT tree
    double betaTime, deltaTime;
    float betaEnergy, betaPSD;
    float dx, dy, dz, displacement;
    float betaZ;
    int betaSegment;
    int alphaX, alphaY;
    int betaX, betaY;
    int multCluster, multClusterIoni;
    int dataSet;
    int direction;
    int lineNumber = 0, lineCounter = 0;
    int fillCount = 0;
    std::size_t index = 0;

    // Invariables
    static constexpr int xBins = 301;  // Number of bins for our histograms
    static constexpr int zBins = 801;
    static constexpr int totalDataLines = 1740;
    static constexpr float xHistogramMax = 150.5;  // Maximum value of bin for histograms
    static constexpr float zHistogramMax = 400.5;  // Maximum value of bin for histograms
    static constexpr float segmentWidth = 145.7;  // Distance between segment centers in mm
    static constexpr float atmosphericScaling = 1.000254;  // Atmosphering scaling coefficient

    char const* dataPath = "2019XList_RxOff.txt";  // Reactor off dataset
    char const* dataFileName = "/home/shay/Documents/PROSPECTData/BiPo_Data/%s/AD1_BiPo.root";

    // Cut values
    static constexpr float highAlphaEnergy = 1.0, lowAlphaEnergy = 0.72;  // Alpha energy cut
    static constexpr float highAlphaPSD = 0.34, lowAlphaPSD = 0.17;  // Alpha PSD cut
    static constexpr float highBetaEnergy = 4.0, lowBetaEnergy = 0;  // Beta energy cut
    static constexpr float highBetaPSD = 0.22, lowBetaPSD = 0.05;  // Beta PSD cut
    static constexpr float n2f = 1 / 12.0;  // Accidental scaling weight
    float tauBiPo = 0.1643 / log(2);  // BiPo lifetime
    float timeStart = 0.01, timeEnd = 3 * tauBiPo;  // Time window for BiPo
    float accTimeStart = 10 * tauBiPo;  // Accidental time window
    float accTimeEnd = accTimeStart + 12 * (timeEnd - timeStart);

    // Storing final counts
    std::array<std::array<float, DirectionSize>, DatasetSize> mean;
    std::array<std::array<float, DirectionSize>, DatasetSize> sigma;

    // Utility functions
    inline bool FiducialCut(int segment)
    {
        if (segment >= 140 || segment % 14 == 0 || (segment + 1) % 14 == 0 || segment == 25
            || segment == 26)
            return true;
        else
            return false;
    }

    // Setting up the leafs
    // Declaration of leaf types
    vector<int>* pseg;
    vector<double>* pt;
    vector<double>* pz;
    vector<double>* pPSD;
    vector<double>* pEtot;
    vector<int>* pmult_clust;
    vector<int>* pmult_clust_ioni;
    vector<int>* fseg;
    vector<double>* ft;
    vector<double>* fz;
    vector<double>* fPSD;
    vector<double>* fEtot;
    vector<int>* fmult_clust;
    vector<int>* fmult_clust_ioni;
    Int_t alphaSegment;
    Double_t alphaEnergy;
    Double_t alphaTime;
    Double_t alphaZ;
    Double_t alphaPSD;
    Int_t multCorrelated;
    Int_t multAccidental;

    // List of branches
    TBranch* b_pseg;
    TBranch* b_pt;
    TBranch* b_pz;
    TBranch* b_pPSD;
    TBranch* b_pEtot;
    TBranch* b_pmult_clust;
    TBranch* b_pmult_clust_ioni;
    TBranch* b_fseg;
    TBranch* b_ft;
    TBranch* b_fz;
    TBranch* b_fPSD;
    TBranch* b_fEtot;
    TBranch* b_fmult_clust;
    TBranch* b_fmult_clust_ioni;
    TBranch* b_aseg;
    TBranch* b_aE;
    TBranch* b_at;
    TBranch* b_az;
    TBranch* b_aPSD;
    TBranch* b_mult_prompt;
    TBranch* b_mult_far;
};

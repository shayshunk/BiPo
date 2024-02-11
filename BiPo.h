#include <array>
#include <fstream>
#include <iostream>
#include <string>

#include "TH1F.h"

// Invariables

// Variables

// Print flags
bool DETECTOR_VERBOSITY = 0;

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
    void FillHistogramUnbiased(int signalSet);
    void FillHistogram();
    void CalculateUnbiasing();
    void SubtractBackgrounds();
    void CalculateCovariances();
    void FillOutputFile();

    // Inline functions
    inline void ResetLineNumber() { lineNumber = 0; }
    inline void ResetLineCounter() { lineCounter = 0; }
    inline void ResetIndex() { index = 0; }

  private:
    // Histogram to count IBDs
    std::array<std::array<std::array<TH1F, DirectionSize>, SignalSize>, DatasetSize> histogram;

    // File list
    std::array<std::string, 4045> files;

    // Values grabbed from ROOT tree
    double Esmear;
    double nCaptTime;
    double xRx;
    double promptPosition, delayedPosition;
    int promptSegment, delayedSegment;
    int dataSet;
    int direction;
    int lineNumber = 0, lineCounter = 0;
    std::size_t index = 0;

    // Invariables
    static constexpr int xBins = 301;  // Number of bins for our histograms
    static constexpr int zBins = 801;
    static constexpr int totalDataLines = 4040;
    static constexpr int totalSimLines = 500;
    static constexpr float histogramMax = 150.5;  // Maximum value of bin for histograms
    static constexpr float segmentWidth = 145.7;  // Distance between segment centers in mm
    static constexpr float atmosphericScaling = 1.000254;  // Atmosphering scaling coefficient
    static constexpr char* dataPath = "2019XList_RxOff.txt";  // Reactor off dataset
    static constexpr char* dataFileName
        = "/home/shay/Documents/PROSPECTData/BiPo_Data/%s/AD1_BiPo.root";

    // Storing final counts
    std::array<std::array<float, DirectionSize>, DatasetSize> effectiveIBD;
    std::array<std::array<float, DirectionSize>, DatasetSize> totalIBD;
    std::array<std::array<float, DirectionSize>, DatasetSize> totalIBDError;
    std::array<std::array<float, DirectionSize>, DatasetSize> mean;
    std::array<std::array<float, DirectionSize>, DatasetSize> sigma;
};

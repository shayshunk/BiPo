#include "BiPo.h"

#include "DetectorConfig.h"
#include "Formatting.h"
#include "Timer.h"

using std::cout, std::string, std::ifstream, std::array, std::getline;

void FillDetectorConfig()
{
    // Filling values based on different periods

    int noSegments = 154;

    int excludeSize = excludeList.size();
    int counter = 0, tmp = 0;

    for (int j = 0; j < noSegments; j++)
    {
        // Filling live segments by checking exclude list
        tmp = excludeList[counter];

        if (j == tmp)
        {
            detectorConfig.push_back(0);
            if (counter + 1 < excludeSize)
            {
                counter += 1;
            }
        }
        else
        {
            detectorConfig.push_back(1);
        }
    }

    if (DETECTOR_VERBOSITY)
    {
        cout << "--------------------------------------------\n";
        cout << "Below is the detector configuration.\n";
        cout << "--------------------------------------------\n";

        for (int j = 140; j >= 0; j -= 14)
        {
            for (int k = 0; k < 14; k++)
            {
                if (detectorConfig[j + k])
                {
                    cout << "\u25A0 ";
                }
                else
                {
                    cout << "\u25A1 ";
                }
            }

            cout << '\n';
        }
        cout << "--------------------------------------------\n";
    }
}

bool CheckNeighbor(int segment, char direction)
{
    // Used for dead segment calculations

    bool neighbor = false;

    switch (direction)
    {
        case 'r':
            neighbor = detectorConfig[segment + 1];
            break;
        case 'l':
            neighbor = detectorConfig[segment - 1];
            break;
        case 'u':
            neighbor = detectorConfig[segment + 14];
            break;
        case 'd':
            neighbor = detectorConfig[segment - 14];
            break;
        default:
            cout << "That direction doesn't exist!\n";
            return false;
    }

    return neighbor;
}

BiPo::BiPo()
{
    for (int dataset = Data; dataset < DatasetSize; dataset++)  // Dataset
    {
        for (int signalSet = Correlated; signalSet < TotalDifference; signalSet++)  // Signal
        {
            for (int direction = X; direction < DirectionSize; direction++)
            {
                string data = DatasetToString(dataset);
                string signal = SignalToString(signalSet);
                string axis = AxisToString(direction);
                string histogramName = data + " " + signal + " " + axis;

                int bins;
                switch (direction)
                {
                    case X:
                        bins = xBins;
                        break;
                    case Y:
                        bins = xBins;
                        break;
                    case Z:
                        bins = zBins;
                        break;
                    default:
                        bins = xBins;
                }

                histogram[dataset][signalSet][direction]
                    = TH1F(histogramName.c_str(), data.c_str(), bins, -histogramMax, histogramMax);
            }
        }
    }

    ResetLineNumber();
}

void BiPo::ReadFileList() {}
void BiPo::SetUpHistograms() {}
void BiPo::FillHistogramUnbiased(int signalSet) {}
void BiPo::FillHistogram() {}
void BiPo::CalculateUnbiasing() {}
void BiPo::SubtractBackgrounds() {}
void BiPo::CalculateCovariances() {}
void BiPo::FillOutputFile() {}
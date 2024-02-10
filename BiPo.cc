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
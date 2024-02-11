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

void BiPo::ReadFileList()
{
    // Opening and checking file list
    ifstream file;
    file.open(dataPath, ifstream::in);

    if (!(file.is_open() && file.good()))
    {
        cout << "File list not found! Exiting.\n";
        cout << "Trying to find: " << dataPath << '\n';
        return;
    }

    while (file.good() && getline(file, files[lineNumber]))
    {
        lineNumber++;
    }
}

void BiPo::SetUpHistograms()
{
    while (index < files.size())
    {
        cout << "Reading file: " << lineCounter << "/" << totalDataLines << '\r';
        cout.flush();

        // Combining names into root file name
        TString rootFilename = Form(dataFileName, files[index].data());

        // Open the root file
        auto rootFile = std::make_unique<TFile>(rootFilename);

        // Grab rootTree and cast to unique pointer
        auto rootTree
            = std::shared_ptr<TTree>(static_cast<TTree*>(rootFile->Get("BiPoTreePlugin/BiPo")));

        long nEntries = rootTree->GetEntries();

        for (long i = 0; i < nEntries; i++)
        {
            rootTree->GetEntry(i);

            // Doing our own fiducial cut
            alphaSegment = rootTree->GetLeaf("aseg")->GetValue(0);

            if (FiducialCut(alphaSegment))
                continue;

            // Grabbing alpha values
            alphaEnergy = rootTree->GetLeaf("aE")->GetValue(0);
            alphaPSD = rootTree->GetLeaf("aPSD")->GetValue(0);
            alphaZ = rootTree->GetLeaf("az")->GetValue(0);

            // Applying alpha cuts
            if (abs(alphaZ) > 1000)
                continue;

            if (alphaEnergy < lowAlphaEnergy || alphaEnergy > highAlphaEnergy)
                continue;

            if (alphaPSD < lowAlphaPSD || alphaEnergy > highAlphaPSD)
                continue;

            multCorrelated = rootTree->GetLeaf("mult_prompt")->GetValue(0);

            FillHistogram(rootTree);
        }
        // Returns the next character in the input sequence, without extracting it: The character is
        // left as the next character to be extracted from the stream
        rootFile->Close();

        lineCounter++;
        index++;
    }
}
void BiPo::FillHistogram(std::shared_ptr<TTree> rootTree)
{
    for (int j = 0; j < multCorrelated; j++)
    {
        // Fiducial cut for beta
        betaSegment = rootTree->GetLeaf("pseg")->GetValue(j);

        if (FiducialCut(betaSegment))
            continue;

        // Grabbing beta values
        betaEnergy = rootTree->GetLeaf("pEtot")->GetValue(j);
        betaPSD = rootTree->GetLeaf("pPSD")->GetValue(0);
        betaZ = rootTree->GetLeaf("pz")->GetValue(0);
        multCluster = rootTree->GetLeaf("pmult_cluster")->GetValue(j);
        multClusterIoni = rootTree->GetLeaf("pmult_cluster_ioni")->GetValue(j);

        // Applying alpha cuts
        if (abs(betaZ) > 1000)
            continue;

        if (betaEnergy < lowBetaEnergy || betaEnergy > highBetaEnergy)
            continue;

        if (betaPSD < lowBetaPSD || betaEnergy > highBetaPSD)
            continue;

        if (multCluster != multClusterIoni)
            continue;

        // Alpha location
        alphaX = alphaSegment % 14;
        alphaY = alphaSegment / 14;

        // Beta location
        betaX = betaSegment % 14;
        betaY = betaSegment / 14;
        betaZ = rootTree->GetLeaf("bz")->GetValue(j);

        // Calculating prompt - delayed displacement
        float dx = 145.7 * (alphaX - betaX);
        float dy = 145.7 * (alphaY - betaY);
        float dz = alphaZ - betaZ;

        float distance = sqrt(dx * dx + dy * dy + dz * dz);

        if (distance > 550)
            continue;

        alphaTime = rootTree->GetLeaf("at")->GetValue(0);
        betaTime = rootTree->GetLeaf("pt")->GetValue(j);

        deltaTime = alphaTime - betaTime;
    }
}

void BiPo::FillHistogramUnbiased(int signalSet) {}
void BiPo::CalculateUnbiasing() {}
void BiPo::SubtractBackgrounds() {}
void BiPo::CalculateCovariances() {}
void BiPo::FillOutputFile() {}

int BiPoDirectionality()
{
    BiPo BiPoDirectionality;

    BiPoDirectionality.ReadFileList();
    BiPoDirectionality.SetUpHistograms();

    return 0;
}
#include "BiPo.h"
#include "DetectorConfig.h"
#include "Formatting.h"
#include "Timer.h"

using std::cout, std::string, std::ifstream, std::vector, std::array, std::getline;

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
                float histogramMax;

                switch (direction)
                {
                    case X:
                        bins = xBins;
                        histogramMax = xHistogramMax;
                        break;
                    case Y:
                        bins = xBins;
                        histogramMax = xHistogramMax;
                        break;
                    case Z:
                        bins = zBins;
                        histogramMax = zHistogramMax;
                        break;
                    default:
                        bins = xBins;
                        histogramMax = xHistogramMax;
                }

                histogram[dataset][signalSet][direction]
                    = TH1D(histogramName.c_str(), data.c_str(), bins, -histogramMax, histogramMax);
            }

            string data = DatasetToString(dataset);
            string signal = SignalToString(signalSet);
            string histogramName = "Multiplicity " + data + " " + signal;

            int bins = 9;

            multiplicity[dataset][signalSet] = TH1I(histogramName.c_str(), data.c_str(), bins, 1, bins);
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
        cout << "Reading file: " << lineCounter + 1 << "/" << totalDataLines << '\r';
        cout.flush();

        // Combining names into root file name
        TString rootFilename = Form(dataFileName, files[index].data());

        // Open the root file
        auto rootFile = std::make_unique<TFile>(rootFilename);

        // Grab rootTree and cast to unique pointer
        auto rootTree = std::shared_ptr<TTree>(static_cast<TTree*>(rootFile->Get("BiPoTreePlugin/BiPo")));

        SetBranchAddresses(rootTree);

        long nEntries = rootTree->GetEntries();

        for (long i = 0; i < nEntries; i++)
        {
            rootTree->GetEntry(i);

            // Doing our own fiducial cut
            if (FiducialCut(alphaSegment))
                continue;

            // Applying alpha cuts
            if (abs(alphaZ) > 1000)
                continue;

            if (alphaEnergy < lowAlphaEnergy || alphaEnergy > highAlphaEnergy)
                continue;

            if (alphaPSD < lowAlphaPSD || alphaPSD > highAlphaPSD)
                continue;

            FillHistogram();
        }

        // rootFile->Close();

        lineCounter++;
        index++;
    }
}

void BiPo::SetBranchAddresses(std::shared_ptr<TTree> rootTree)
{
    // Set object pointer
    pseg = 0;
    pt = 0;
    pz = 0;
    pPSD = 0;
    pEtot = 0;
    pmult_clust = 0;
    pmult_clust_ioni = 0;
    fseg = 0;
    ft = 0;
    fz = 0;
    fPSD = 0;
    fEtot = 0;
    fmult_clust = 0;
    fmult_clust_ioni = 0;

    // Set branch addresses and branch pointers
    if (!rootTree)
        return;

    // Prompt Window
    rootTree->SetBranchAddress("pseg", &pseg, &b_pseg);  // beta segment number
    rootTree->SetBranchAddress("pt", &pt, &b_pt);  // beta timing in us
    rootTree->SetBranchAddress("pz", &pz, &b_pz);  // beta position in Z position, given in mm
    rootTree->SetBranchAddress("pPSD", &pPSD, &b_pPSD);  // beta PSD
    rootTree->SetBranchAddress("pEtot", &pEtot, &b_pEtot);  // beta total energy in MeV
    rootTree->SetBranchAddress("pmult_clust", &pmult_clust, &b_pmult_clust);  // prompt cluster
                                                                              // multiplicity
    rootTree->SetBranchAddress("pmult_clust_ioni", &pmult_clust_ioni, &b_pmult_clust_ioni);  // prompt cluster multiplicity
                                                                                             // ionization?
    // Far Window
    rootTree->SetBranchAddress("fseg", &fseg, &b_fseg);  // beta segment number
    rootTree->SetBranchAddress("ft", &ft, &b_ft);  // beta timing in us
    rootTree->SetBranchAddress("fz", &fz, &b_fz);  // beta position in Z position, given in mm
    rootTree->SetBranchAddress("fPSD", &fPSD, &b_fPSD);  // beta PSD
    rootTree->SetBranchAddress("fEtot", &fEtot, &b_fEtot);  // beta total energy in MeV
    rootTree->SetBranchAddress("fmult_clust", &fmult_clust, &b_fmult_clust);  // prompt cluster
                                                                              // multiplicity
    rootTree->SetBranchAddress("fmult_clust_ioni", &fmult_clust_ioni, &b_fmult_clust_ioni);  // prompt cluster multiplicity
                                                                                             // ionization?
    // Alpha
    rootTree->SetBranchAddress("aseg", &alphaSegment, &b_aseg);  // alpha segment number
    rootTree->SetBranchAddress("aE", &alphaEnergy, &b_aE);  // alpha energy in MeV
    rootTree->SetBranchAddress("at", &alphaTime, &b_at);  // alpha timing in us
    rootTree->SetBranchAddress("az", &alphaZ, &b_az);
    rootTree->SetBranchAddress("aPSD", &alphaPSD, &b_aPSD);  // alpha PSD
    rootTree->SetBranchAddress("mult_prompt", &multCorrelated, &b_mult_prompt);  // prompt
                                                                                 // multiplicity for
                                                                                 // correlated
    rootTree->SetBranchAddress("mult_far", &multAccidental, &b_mult_far);  // prompt multiplicity
                                                                           // for accidentals
}

void BiPo::FillHistogram()
{
    for (int j = 0; j < multCorrelated; j++)
    {
        // Fiducial cut for beta
        betaSegment = pseg->at(j);

        if (FiducialCut(betaSegment))
            continue;

        // Grabbing beta values
        betaEnergy = pEtot->at(j);
        betaPSD = pPSD->at(j);
        betaZ = pz->at(j);
        multCluster = pmult_clust->at(j);
        multClusterIoni = pmult_clust_ioni->at(j);

        // Applying alpha cuts
        if (abs(betaZ) > 1000)
            continue;

        if (betaEnergy < lowBetaEnergy || betaEnergy > highBetaEnergy)
            continue;

        if (betaPSD < lowBetaPSD || betaPSD > highBetaPSD)
            continue;

        if (multCluster != multClusterIoni)
            continue;

        // Alpha location
        alphaX = alphaSegment % 14;
        alphaY = alphaSegment / 14;

        // Beta location
        betaX = betaSegment % 14;
        betaY = betaSegment / 14;
        betaZ = pz->at(j);

        // Calculating prompt - delayed displacement
        dx = 145.7 * (alphaX - betaX);
        dy = 145.7 * (alphaY - betaY);
        dz = alphaZ - betaZ;

        displacement = sqrt(dx * dx + dy * dy + dz * dz);

        if (displacement > 550)
            continue;

        betaTime = pt->at(j);

        deltaTime = alphaTime - betaTime;

        if (deltaTime > timeStart && deltaTime < timeEnd)
        {
            if (alphaSegment == betaSegment + 1 || alphaSegment == betaSegment - 1)
            {
                histogram[Data][Correlated][X].Fill(dx);
            }

            if (alphaSegment == betaSegment + 14 || alphaSegment == betaSegment - 14)
                histogram[Data][Correlated][Y].Fill(dy);

            if (alphaSegment == betaSegment)
            {
                histogram[Data][Correlated][X].Fill(0.0);
                histogram[Data][Correlated][Y].Fill(0.0);
                histogram[Data][Correlated][Z].Fill(dz);
                FillHistogramUnbiased(Correlated);
            }

            multiplicity[Data][Correlated].Fill(j + 1);
        }
    }

    for (int j = 0; j < multAccidental; j++)
    {
        // Fiducial cut for beta
        betaSegment = fseg->at(j);

        if (FiducialCut(betaSegment))
            continue;

        // Grabbing beta values
        betaEnergy = fEtot->at(j);
        betaPSD = fPSD->at(j);
        betaZ = fz->at(j);
        multCluster = fmult_clust->at(j);
        multClusterIoni = fmult_clust_ioni->at(j);

        // Applying alpha cuts
        if (abs(betaZ) > 1000)
            continue;

        if (betaEnergy < lowBetaEnergy || betaEnergy > highBetaEnergy)
            continue;

        if (betaPSD < lowBetaPSD || betaPSD > highBetaPSD)
            continue;

        if (multCluster != multClusterIoni)
            continue;

        // Alpha location
        alphaX = alphaSegment % 14;
        alphaY = alphaSegment / 14;

        // Beta location
        betaX = betaSegment % 14;
        betaY = betaSegment / 14;
        betaZ = fz->at(j);

        // Calculating prompt - delayed displacement
        dx = 145.7 * (alphaX - betaX);
        dy = 145.7 * (alphaY - betaY);
        dz = alphaZ - betaZ;

        displacement = sqrt(dx * dx + dy * dy + dz * dz);

        if (abs(dz) > 250)
            continue;

        if (displacement > 550)
            continue;

        betaTime = ft->at(j);

        deltaTime = betaTime - alphaTime;

        if (deltaTime > accTimeStart && deltaTime < accTimeEnd)
        {
            if (alphaSegment == betaSegment + 1 || alphaSegment == betaSegment - 1)
                histogram[Data][Accidental][X].Fill(dx, n2f);

            if (alphaSegment == betaSegment + 14 || alphaSegment == betaSegment - 14)
                histogram[Data][Accidental][Y].Fill(dy, n2f);

            if (alphaSegment == betaSegment)
            {
                histogram[Data][Accidental][X].Fill(0.0, n2f);
                histogram[Data][Accidental][Y].Fill(0.0, n2f);
                histogram[Data][Accidental][Z].Fill(dz, n2f);
                FillHistogramUnbiased(Accidental);
            }

            multiplicity[Data][Accidental].Fill(j + 1, n2f);
        }
    }
}

void BiPo::FillHistogramUnbiased(int signalSet)
{
    bool posDirectionX = false, negDirectionX = false;
    bool posDirectionY = false, negDirectionY = false;

    // Need to weight accidental datasets by deadtime correction factor
    double weight = (signalSet == Accidental) ? n2f : 1;

    // Check for live neighbors in different directions
    posDirectionX = CheckNeighbor(alphaSegment, 'r');
    negDirectionX = CheckNeighbor(alphaSegment, 'l');
    posDirectionY = CheckNeighbor(alphaSegment, 'u');
    negDirectionY = CheckNeighbor(alphaSegment, 'd');

    // Filling x axis
    if (posDirectionX && !negDirectionX)
        histogram[DataUnbiased][signalSet][X].Fill(segmentWidth, weight);
    else if (!posDirectionX && negDirectionX)
        histogram[DataUnbiased][signalSet][X].Fill(-segmentWidth, weight);
    else if (posDirectionX && negDirectionX)
        histogram[DataUnbiased][signalSet][X].Fill(0.0, weight);

    // Filling y axis
    if (posDirectionY && !negDirectionY)
        histogram[DataUnbiased][signalSet][Y].Fill(segmentWidth, weight);
    else if (!posDirectionY && negDirectionY)
        histogram[DataUnbiased][signalSet][Y].Fill(-segmentWidth, weight);
    else if (posDirectionY && negDirectionY)
        histogram[DataUnbiased][signalSet][Y].Fill(0.0, weight);

    // Filling z axis
    histogram[DataUnbiased][signalSet][Z].Fill(dz, weight);
}

void BiPo::SubtractBackgrounds()
{
    for (int dataset = Data; dataset < DatasetSize; dataset++)
    {
        for (int direction = X; direction < DirectionSize; direction++)
        {
            string histogramName;
            string data = DatasetToString(dataset);
            histogramName = data + " Total Difference " + AxisToString(direction);

            // Copying Correlated to start
            histogram[dataset][TotalDifference][direction] = TH1D(histogram[dataset][Correlated][direction]);

            histogram[dataset][TotalDifference][direction].SetNameTitle(histogramName.c_str(), data.c_str());

            histogram[dataset][TotalDifference][direction].Add(&histogram[dataset][Accidental][direction], -1.0);

            if (dataset == DataUnbiased || direction == Z)
                continue;

            mean[dataset][direction] = histogram[dataset][TotalDifference][direction].GetMean();

            float effEntries = histogram[dataset][TotalDifference][direction].GetEntries();
            sigma[dataset][direction] = histogram[dataset][TotalDifference][direction].GetStdDev() / sqrt(effEntries);
        }

        // Z is fit to a Guassian and only takes same segment inputs
        // Possible thanks to 1mm resolution in Z
        TF1 gaussian("Fit", "gaus", -250, 250);

        histogram[dataset][TotalDifference][Z].Fit("Fit", "RQ");

        float zMean = gaussian.GetParameter(1);
        float zError = gaussian.GetParError(1);

        mean[dataset][Z] = zMean;
        sigma[dataset][Z] = zError;

        // Deleting fit because we don't want the plot options stuck here
        delete histogram[dataset][TotalDifference][Z].GetListOfFunctions()->FindObject("Fit");

        // Multiplicity background subtraction
        string histogramName;
        string data = DatasetToString(dataset);
        histogramName = "Multiplicity " + data + " Total Difference";

        multiplicity[dataset][TotalDifference] = TH1I(multiplicity[dataset][Correlated]);
        multiplicity[dataset][TotalDifference].SetNameTitle(histogramName.c_str(), data.c_str());
        multiplicity[dataset][TotalDifference].Add(&multiplicity[dataset][Accidental], -1.0);
    }
}

void BiPo::CalculateUnbiasing()
{
    // Defining variables used in calculation. Check the error propagation technote for details on
    // the method
    double rPlus = 0, rMinus = 0;
    double p = 0, pError = 0;
    double nPlus = 0, nPlusPlus = 0, nMinus = 0, nMinusMinus = 0, nPlusMinus = 0;
    double nPlusError = 0, nPlusPlusError = 0, nMinusError = 0, nMinusMinusError = 0, nPlusMinusError = 0;

    for (int direction = X; direction < Z; direction++)
    {
        // Grabbing data from filled bins, rest should be empty
        nPlus = histogram[Data][TotalDifference][direction].GetBinContent(297);
        nPlusPlus = histogram[DataUnbiased][TotalDifference][direction].GetBinContent(297);
        nMinus = histogram[Data][TotalDifference][direction].GetBinContent(5);
        nMinusMinus = histogram[DataUnbiased][TotalDifference][direction].GetBinContent(5);
        nPlusMinus = histogram[DataUnbiased][TotalDifference][direction].GetBinContent(151);

        nPlusError = histogram[Data][TotalDifference][direction].GetBinError(297);
        nPlusPlusError = histogram[DataUnbiased][TotalDifference][direction].GetBinError(297);
        nMinusError = histogram[Data][TotalDifference][direction].GetBinError(5);
        nMinusMinusError = histogram[DataUnbiased][TotalDifference][direction].GetBinError(5);
        nPlusMinusError = histogram[DataUnbiased][TotalDifference][direction].GetBinError(151);

        rPlus = nPlus / (nPlusPlus + nPlusMinus);
        rMinus = nMinus / (nMinusMinus + nPlusMinus);

        p = segmentWidth * (rPlus - rMinus) / (rPlus + rMinus + 1);

        pError
            = segmentWidth
              * pow(1 / ((nMinus * (nPlusMinus + nPlusPlus) + (nMinusMinus + nPlusMinus) * (nPlus + nPlusMinus + nPlusPlus))), 2)
              * sqrt(pow((nMinusMinus + nPlusMinus) * (nPlusMinus + nPlusPlus), 2)
                         * (pow(nPlusError * (2 * nMinus + nMinusMinus + nPlusMinus), 2)
                            + pow(nMinusError * (2 * nPlus + nPlusPlus + nPlusMinus), 2))
                     + pow((nPlus * (nPlusMinus + nMinusMinus) * (2 * nMinus + nMinusMinus + nPlusMinus) * nPlusPlusError), 2)
                     + pow((nPlusMinusError
                            * (nPlus * pow((nMinusMinus + nPlusMinus), 2)
                               + nMinus * (2 * nMinusMinus * nPlus - 2 * nPlus * nPlusPlus - pow((nPlusMinus + nPlusPlus), 2)))),
                           2)
                     + pow((nMinus * (nPlusMinus + nPlusPlus) * (2 * nPlus + nPlusMinus + nPlusPlus) * nMinusMinusError), 2));

        mean[DataUnbiased][direction] = p;
        sigma[DataUnbiased][direction] = pError;

        if (NCOUNT_VERBOSITY)
        {
            cout << "N counts for: " << boldOn << "Data Unbiased " << AxisToString(direction) << '\n';
            cout << "N+: " << resetFormats << nPlus << '\n';
            cout << boldOn << "N-: " << resetFormats << nMinus << '\n';
            cout << boldOn << "N++: " << resetFormats << nPlusPlus << '\n';
            cout << boldOn << "N--: " << resetFormats << nMinusMinus << '\n';
            cout << boldOn << "N+-: " << resetFormats << nPlusMinus << '\n';
            cout << "--------------------------------------------\n";
        }
    }

    mean[DataUnbiased][Z] = mean[Data][Z];
    sigma[DataUnbiased][Z] = sigma[Data][Z];

    cout << "--------------------------------------------\n";
    cout << boldOn << cyanOn << "Calculated Means.\n" << resetFormats;
    cout << "--------------------------------------------\n";

    // Printing out values
    for (int dataset = Data; dataset < DatasetSize; dataset++)
    {
        cout << "Mean and sigma values for: " << boldOn << DatasetToString(dataset) << resetFormats << '\n';
        for (int direction = X; direction < DirectionSize; direction++)
        {
            cout << boldOn << "p" << AxisToString(direction) << ": " << resetFormats << mean[dataset][direction] << " ± "
                 << sigma[dataset][direction] << '\n';
        }
        cout << "--------------------------------------------\n";
    }
}

void BiPo::CalculateAngles()
{
    // Defining variables for readability of code
    double px, py, pz;
    double sigmaX, sigmaY, sigmaZ;

    for (int dataset = Data; dataset < DatasetSize; dataset++)
    {
        // Grabbing values from current dataset
        px = mean[dataset][X];
        py = mean[dataset][Y];
        pz = mean[dataset][Z];

        sigmaX = sigma[dataset][X];
        sigmaY = sigma[dataset][Y];
        sigmaZ = sigma[dataset][Z];

        // phi = arctan(y / x)
        float tanPhi = py / px;
        float phiTemp = atan(tanPhi) * 180.0 / pi;
        float tanPhiError = sqrt(pow((sigmaX * py) / pow(px, 2), 2) + pow(sigmaY / px, 2));
        float phiErrorTemp = (tanPhiError / (1 + pow(tanPhi, 2))) * 180.0 / pi;

        // theta = arctan(z / sqrt(x^2 + y^2))
        float tanTheta = pz / sqrt(pow(px, 2) + pow(py, 2));
        float thetaTemp = atan(tanTheta) * 180.0 / pi;
        float tanThetaError = sqrt((1 / (px * px + py * py))
                                   * (pow((px * pz * sigmaX / (px * px + py * py)), 2)
                                      + pow((py * pz * sigmaY / (px * px + py * py)), 2) + pow(sigmaZ, 2)));
        float thetaErrorTemp = (tanThetaError / (1 + pow(tanTheta, 2))) * 180.0 / pi;

        // Storing values
        phi[dataset] = phiTemp;
        phiError[dataset] = phiErrorTemp;
        theta[dataset] = thetaTemp;
        thetaError[dataset] = thetaErrorTemp;
    }
}

void BiPo::OffsetTheta()
{
    for (int dataset = Data; dataset < DatasetSize; dataset++)
    {
        theta[dataset] = 90 - theta[dataset];
    }
}

void BiPo::PrintAngles()
{
    cout << boldOn << cyanOn << "Final Angles!\n" << resetFormats;
    cout << "--------------------------------------------\n";

    for (int dataset = Data; dataset < DatasetSize; dataset++)
    {
        cout << "Angle values for: " << boldOn << DatasetToString(dataset) << resetFormats << '\n';
        cout << greenOn;
        cout << boldOn << underlineOn << "ϕ:" << resetFormats << greenOn << " " << phi[dataset] << "\u00B0 ± "
             << phiError[dataset] << "\u00B0.\n";
        cout << boldOn << underlineOn << "θ:" << resetFormats << greenOn << " " << theta[dataset] << "\u00B0 ± "
             << thetaError[dataset] << "\u00B0.\n"
             << resetFormats;
        cout << "--------------------------------------------\n";
    }
}

void BiPo::FillOutputFile()
{
    // Set up our output file
    TFile outputFile("BiPo.root", "recreate");

    outputFile.cd();

    for (int dataset = Data; dataset < DatasetSize; dataset++)
    {
        for (int signalSet = Correlated; signalSet < SignalSize; signalSet++)
        {
            for (int direction = X; direction < DirectionSize; direction++)
            {
                histogram[dataset][signalSet][direction].Write();
            }

            multiplicity[dataset][signalSet].Write();
        }
    }

    cout << boldOn << cyanOn << "Filled output file: " << resetFormats << blueOn << boldOn << "BiPo.root!\n" << resetFormats;
    cout << "--------------------------------------------\n";

    outputFile.Close();
}

int main(int argc, char* argv[])
{
    // Ignoring warnings
    gErrorIgnoreLevel = kError;

    // Using command line arguments for verbosity control
    for (int i = 1; i < argc; i++)
    {
        if (string(argv[i]) == "-D")
            DETECTOR_VERBOSITY = 1;
        else if (string(argv[i]) == "-N")
            NCOUNT_VERBOSITY = 1;
    }

    // Timing everything
    Timer timer;

    // Filling detector configuration
    FillDetectorConfig();

    // Setting up directionality class
    BiPo directionality;

    // Running analysis
    directionality.ReadFileList();
    directionality.SetUpHistograms();
    directionality.SubtractBackgrounds();
    directionality.CalculateUnbiasing();
    directionality.CalculateAngles();
    directionality.OffsetTheta();
    directionality.PrintAngles();
    directionality.FillOutputFile();

    return 0;
}
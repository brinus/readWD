// *********************************************************************************
// *                                                                               *
// *                   WaveDREAM Binary Data Analysis Tool                         *
// *                                                                               *
// *                                                                               *
// * Name:     readWD.cpp                                                          *
// * Author:   Patrick Schwendimann                                                *
// * E-Mail:   patrick.schwendimann@psi.ch                                         *
// * Version:  0.9.4 Beta (22.06.2020)                                             *
// *                                                                               *
// *                                                                               *
// *  Based on the old DRS4 binary data analysis programm orignially written for   *
// *  the ATAR experiment by Angela Papa <angela.papa@psi.ch>, later updated by    *
// *  Emanuele Ripiccini, Fred Gray and Giada Rutar and an example program         *
// *  written by Stefan Ritt <stefan.ritt@psi.ch> for the reading of a binary      *
// *  file written by the DRSOsc program.                                          *
// *                                                                               *
// *  Disclaimer: This code may not be fully suitable for whatever you may intend  *
// *              to use it for. No warranty is given. Use it always together with *
// *              some common sense.                                               *
// *              Should you lack sanity or reason ... you better don't use this   *
// *              code. ;)                                                         *
// *                                                                               *
// *  Compilation requires the root packages downloadable from the CERN website    *
// *                                                                               *
// *    https://root.cern.ch/downloading-root                                      *
// *                                                                               *
// *  and the curses library (should be installed by default)                      *
// *                                                                               *
// *                                                                               *
// *********************************************************************************

#define __READ_VERSION__ "0.9.4 Beta (22.06.2020)"

// -------------------------------- Includes ---------------------------------------

// Standard libraries
#include <curses.h>
#include <fstream>
#include <iostream>
#include <cstring>
#include <stdio.h>
#include <vector>

// Root libraries
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1F.h>
#include <TMath.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <TMatrixD.h>
#include <TVectorD.h>

// -------------------------------- Structs ----------------------------------------

/*! @struct FileHeader

    @brief File header.

    @details First line of binary data is the file header.
    It is made by a 4-byte variable: the fourth one is the version of the DRS.

    File header ||||
    -|-|-|-|
    'D' | 'R' | 'S' | '2' |
*/
struct FileHeader
{
    char tag[3];  ///< The tag of the file header.
    char version; ///< The version of the board.
};

/*! @struct EventHeader

    @brief Event header.

    @details Struct to store the event header. The binary file is encoded as in the next table.
    Serial number starts with 1, event date/time has 16-bit values.

    Event header ||||
    :--:|:---:|:---:|:---:|
    'E' | 'H' | 'D' | 'R' |
    Event serial number||||
    Year|| Mont ||
    Day|| Hour ||
    Minute|| Second ||
    Millisecond || Range ||
*/
struct EventHeader
{
    char tag[4];
    unsigned int serialNumber;
    unsigned short year;
    unsigned short month;
    unsigned short day;
    unsigned short hour;
    unsigned short min;
    unsigned short sec;
    unsigned short ms;
    unsigned short rangeCenter;
};

/*! @struct IntegrationWindow

    @brief Integration window.

    @details Simple struct to store infos on the integration window.
*/
struct IntegrationWindow
{
    int start; ///< Start of the integration window.
    int stop;  ///< Stop of the integration window.
};

/*! @struct Configuration

    @brief Configuration struct.

    @details Struct to store important variables.
*/
struct Configuration
{
    Configuration();
    bool firstOfRun;                                   ///< Setted true when is the first event of the run.
                                                       ///<
    short debug;                                       ///< Debug flag, passed as input `-d` or `--DEBUG` or setted from the configuration menu.
                                                       ///<
    short sigWF;                                       ///< To select wether the waveform must be positive or negative sign,
                                                       ///< setted from the configuration menu.
    bool subtractSine;                                 ///< Flag to subtract sine noise or not, setted from the configuration menu.
                                                       ///<
    float noiseFrequency;                              ///< The noise frequency, adjustable from the configuration menu.
                                                       ///<
    int runMode;                                       ///< 0 for DRS, 1 for WDB.
                                                       ///<
    float intRise;                                     ///< Number of Std.Dev. necessary to detect a signal.
                                                       ///<
    float intDecay;                                    ///< Number of Tau to set the @ref IntegrationWindow stop attribute.
                                                       ///<
    int run;                                           ///< Run number.
                                                       ///<
    unsigned int nSampleEvents;                        ///< Number of samples for integration window.
                                                       ///<
    unsigned int nSaveEvents;                          ///< Number of Events to save Waveforms.
                                                       ///<
    float cfFraction;                                  ///<
                                                       ///<
    float leThreshold;                                 ///<
                                                       ///<
    std::vector<unsigned int> nChannelsPerBoard;       ///< Number of channels for every board.
                                                       ///<
    std::vector<IntegrationWindow> integrationWindows; ///< Integration windows for every board.
                                                       ///<
    std::vector<std::vector<float *>> timeBinWidth;    ///< Width of time bin for every board.
                                                       ///<
    TFile *theFile;                                    ///< The ROOT path/file where the results will be stored.
                                                       ///<
};

/*! @brief @ref Configuration constructor

    @details Every new @ref Configuration variable is initialized as follows:
    @code{.cpp}
    Configuration::Configuration()
    {
        debug = 0;
        firstOfRun = true;
        sigWF = -1;
        subtractSine = false;
        noiseFrequency = 50.6e+6 * TMath::TwoPi();
        runMode = 1;
        cfFraction = 0.2;
        leThreshold = 0.05;
        run = 0;
        nSampleEvents = 2000;
        nSaveEvents = 0;
        theFile = 0;
        intRise = 6;
        intDecay = 3;
    }
    @endcode
*/
Configuration::Configuration()
{
    debug = 0;
    firstOfRun = true;
    sigWF = -1;
    subtractSine = false;
    noiseFrequency = 50.6e+6 * TMath::TwoPi();
    runMode = 1;
    cfFraction = 0.2;
    leThreshold = 0.05;
    run = 0;
    nSampleEvents = 2000;
    nSaveEvents = 0;
    theFile = 0;
    intRise = 6;
    intDecay = 3;
}

// ------------------------------- Global Variables --------------------------------

Configuration gCONFIG;
#define Debug          \
    if (gCONFIG.debug) \
    std::cout << "DEBUG L" << __LINE__ << ": "
#define DEBUG              \
    if (gCONFIG.debug > 1) \
    std::cout << "DEBUG L" << __LINE__ << ": "
#define SAMPLES_PER_WAVEFORM 1024

// -------------------------------- Declarations -----------------------------------

std::ifstream *Initialise(std::string);
int GetIntegrationBounds(std::ifstream *);
int ReadAnEvent(std::ifstream *, EventHeader &, std::vector<std::vector<float *>> &, std::vector<std::vector<unsigned short>> *);
int ReadFile(std::ifstream *);
int Config();
void PrintHelp();
void subtractSineNoise(float *, float *, float &, float &, char *);

/**
 * @brief Determines mean and standard deviation of pedestal
 *
 * @details Determine standard deviation of the first 100 samples and the second 100 samples.
 * The range with the smaller standard deviation will be used to determine the pedestal for the waveform.
 *
 * @param aWaveform
 * @param pedestal
 * @param stdv
 */
void getPedestal(float *, float &, float &);

/**
 * @brief Determines mean and standard deviation of pedestal
 *
 * @details Determine standard deviation of the first 100 samples and the second 100 samples.
 * The range with the smaller standard deviation will be used to determine the pedestal for the waveform.
 *
 * @param hSignal
 * @param pedestal
 * @param stdv
 */
void getPedestal(TH1F *, float &, float &);
float getTimeStamp(EventHeader);

// ------------------------------ Functions -----------------------------------------

void getPedestal(float *aWaveform, float &pedestal, float &stdv)
{
    pedestal = 0;
    stdv = 0;
    float vStdv1 = TMath::StdDev(100, aWaveform);
    float vStdv2 = TMath::StdDev(100, aWaveform + 100);
    if (vStdv1 < vStdv2)
    {
        pedestal = TMath::Mean(100, aWaveform);
        stdv = vStdv1;
    }
    else
    {
        pedestal = TMath::Mean(100, aWaveform + 100);
        stdv = vStdv2;
    }
}

void getPedestal(TH1F *hSignal, float &pedestal, float &stdv)
{
    pedestal = stdv = 0;
    int nBins = 100;
    float meanSquares1 = 0;
    float mean1 = 0;
    float stdv1 = 0;
    float meanSquares2 = 0;
    float mean2 = 0;
    float stdv2 = 0;
    float x = 0;
    for (int i = 0; i < nBins; ++i)
    {
        x = hSignal->GetBinContent(i + 1);
        mean1 += x;
        meanSquares1 += x * x;
        x = hSignal->GetBinContent(i + nBins + 1);
        mean2 += x;
        meanSquares2 += x * x;
    }
    mean1 /= nBins;
    mean2 /= nBins;
    stdv1 = sqrt(meanSquares1 / nBins - mean1 * mean1);
    stdv2 = sqrt(meanSquares2 / nBins - mean2 * mean2);
    if (stdv1 < stdv2)
    {
        pedestal = mean1;
        stdv = stdv1;
    }
    else
    {
        pedestal = mean2;
        stdv = stdv2;
    }
}

/**
 * @brief Subtract sine noise from waveform
 *
 * @details Subtract sine noise from waveform. To be further understood.
 *
 * @param wfData
 * @param wfTime
 * @param ampl
 * @param phi
 * @param name
 */
void subtractSineNoise(float *wfData, float *wfTime, float &ampl, float &phi, char *name = 0)
{
    float freq = gCONFIG.noiseFrequency;
    TMatrixD M(3, 3);
    TVectorD v(3);
    for (int i = 0; i < 200; ++i)
    {
        float co = TMath::Cos(freq * wfTime[i]);
        float si = TMath::Sin(freq * wfTime[i]);
        M[0][0] += co * co;
        M[0][1] += co * si;
        M[0][2] += co;
        M[1][1] += si * si;
        M[1][2] += si;
        M[2][2] += 1;
        v[0] += co * wfData[i];
        v[1] += si * wfData[i];
        v[2] += wfData[i];
    }
    M[1][0] = M[0][1];
    M[2][0] = M[0][2];
    M[2][1] = M[1][2];
    M.Invert();
    v = M * v;
    ampl = TMath::Sqrt(v[0] * v[0] + v[1] * v[1]);
    phi = TMath::ATan2(v[0], v[1]);
    if (name)
    {
        TGraph *aGraph = new TGraph(SAMPLES_PER_WAVEFORM, wfTime, wfData);
        TF1 *noise = new TF1("SineNoise", "[0] * sin([1] * x + [2]) + [3]");
        noise->SetNpx(1000);
        noise->SetParameters(ampl, freq, phi, v[2]);
        aGraph->GetListOfFunctions()->Add(noise);
        aGraph->SetName(name);
        aGraph->SetTitle(name);
        aGraph->Write("", TObject::kOverwrite);
    }
    for (int i = 0; i < SAMPLES_PER_WAVEFORM; ++i)
    {
        wfData[i] -= ampl * TMath::Sin(freq * wfTime[i] + phi);
    }
}

/**
 * @brief Get the timestamp
 *
 * @details Get the timestamp from the event header.
 *
 * @param eh The event header
 * @return Returns the timestamp
 */
float getTimeStamp(EventHeader eh)
{
    static int lastDay = 0;
    static int nDays = 0;
    if (eh.serialNumber == 0)
    {
        lastDay = eh.day;
        nDays = 0;
    }
    if (lastDay != eh.day)
    {
        lastDay = eh.day;
        ++nDays;
    }
    float timestamp = eh.ms / 1000. + eh.sec + eh.min * 60 + eh.hour * 3600 + nDays * 86400;
    return timestamp;
}

/**
 * @brief
 *
 * @param filename
 * @return std::ifstream*
 */
std::ifstream *Initialise(std::string filename)
{
    unsigned int nBoards = 1;
    unsigned int board = 0;
    float *times = 0;
    char word[5];
    word[4] = '\0';
    std::ifstream *file = new std::ifstream;
    file->open(filename, std::ios::in | std::ios::binary);

    // Reset the configuration
    gCONFIG.firstOfRun = true;
    gCONFIG.nChannelsPerBoard.resize(0);
    gCONFIG.integrationWindows.resize(0);
    gCONFIG.timeBinWidth.resize(0);

    // Return if unable to open the file
    if (!file->is_open())
    {
        std::cerr << "!! Could not open file " << filename << std::endl;
        delete file;
        return 0;
    }

    FileHeader fh;
    file->read((char *)&fh, sizeof(fh));

    // Check file header
    if (std::memcmp(fh.tag, "DRS", 3) != 0)
    {
        std::cerr << "!! No suitable file hedaer in file " << filename << std::endl;
        file->close();
        delete file;
        return 0;
    }
    if (fh.version - '0' < 8)
    {
        gCONFIG.runMode = 0; // DRS
        std::cout << "DRS Version " << fh.version << " : DRS Board" << std::endl;
    }
    else
    {
        gCONFIG.runMode = 1; // WDB
        std::cout << "DRS Version " << fh.version << " : WDB encoding" << std::endl;
    }

    // Skip the time 'TIME' header
    file->seekg(4, file->cur);
    file->read(word, 4);

    // Check board number
    if (word[0] != 'B')
    {
        std::cerr << "!! No board header found in file " << filename << std::endl;
        file->close();
        delete file;
        return 0;
    }
    std::cout << "Board 1: " << word[0] << word[1] << *(short *)(word + 2) << std::endl;
    gCONFIG.nChannelsPerBoard.push_back(0);
    gCONFIG.timeBinWidth.push_back({});
    while (word[0] == 'C' || word[0] == 'B')
    {
        file->read(word, 4);
        if (word[0] == 'C')
        {
            // Channel
            gCONFIG.nChannelsPerBoard.at(board) += 1;
            times = new float[SAMPLES_PER_WAVEFORM];
            file->read((char *)times, sizeof(float) * SAMPLES_PER_WAVEFORM);
            gCONFIG.timeBinWidth.at(board).push_back(times);
        }
        else if (word[0] == 'B')
        {
            // Board
            gCONFIG.nChannelsPerBoard.push_back(0);
            gCONFIG.timeBinWidth.push_back({});
            ++board;
            ++nBoards;
            std::cout << "Board " << nBoards << ": " << word[0] << word[1] << *(short *)(word + 2) << std::endl;
        }
    }

    // Check for event header
    if (strcmp(word, "EHDR") != 0)
    {
        std::cerr << "!! No event header found: Expected EHDR, got " << word << std::endl;
        Debug << "Cleaning up ... " << std::endl;
        file->close();
        for (board = 0; board < gCONFIG.timeBinWidth.size(); ++board)
        {
            for (float *arr : gCONFIG.timeBinWidth.at(board))
                delete[] arr;
        }
        delete file;
        return 0;
    }

    file->seekg(-4, file->cur);
    return file;
}

/**
 * @brief Get the Integration Bounds object
 *
 * @param file
 * @return int
 */
int GetIntegrationBounds(std::ifstream *file)
{
    int nChannelsTot = 0;
    char name[32];
    for (int nChannelsAtBoard : gCONFIG.nChannelsPerBoard)
    {
        nChannelsTot += nChannelsAtBoard;
    }

    // Create SampleSignal histograms
    std::vector<TH1F **> hSampleSig;
    for (int board = 0; board < gCONFIG.nChannelsPerBoard.size(); board++)
    {
        hSampleSig.push_back(new TH1F *[gCONFIG.nChannelsPerBoard.at(board)]);
        for (int ch = 0; ch < gCONFIG.nChannelsPerBoard.at(board); ++ch)
        {
            sprintf(name, "hSampleSigB%03dC%03d", board + 1, ch);
            hSampleSig.at(board)[ch] = new TH1F(name, name, SAMPLES_PER_WAVEFORM, -0.5, SAMPLES_PER_WAVEFORM - 0.5);

            // Attribute the histos with the root output file. Memory will be freed when closing the file
            hSampleSig.at(board)[ch]->SetDirectory(gCONFIG.theFile);
            Debug << "Created histogram " << name << std::endl;
        }
    }

    // Read some data and fill the histos
    EventHeader eh;
    std::vector<std::vector<float *>> wfData;
    for (int board = 0; board < gCONFIG.nChannelsPerBoard.size(); board++)
    {
        wfData.push_back({});
        for (int ch = 0; ch < gCONFIG.nChannelsPerBoard.at(board); ++ch)
            wfData.at(board).push_back(new float[SAMPLES_PER_WAVEFORM]);
    }
    float pedestal = 0;
    float pedStdv = 0;
    unsigned int nEvents = 0;
    int retVal = 0;
    int curPos = file->tellg();
    while (nEvents < gCONFIG.nSampleEvents && !file->eof())
    {
        retVal = ReadAnEvent(file, eh, wfData, 0);
        if (retVal)
            break;
        for (int board = 0; board < gCONFIG.nChannelsPerBoard.size(); board++)
        {
            for (int ch = 0; ch < gCONFIG.nChannelsPerBoard.at(board); ++ch)
            {
                getPedestal(wfData.at(board).at(ch), pedestal, pedStdv);
                for (int i = 0; i < SAMPLES_PER_WAVEFORM; ++i)
                {
                    hSampleSig.at(board)[ch]->Fill(i, wfData.at(board).at(ch)[i] - pedestal);
                }
            }
        }
    }
    if (retVal < 2)
    {
        float peak, stdv, tau;
        int peakBin;
        TH1F *hSample = 0;
        TF1 *gaus = new TF1("gaus", "gaus(0)");
        TF1 *decay = new TF1("decay", "[0] * exp( -(x - [1]) / [2])");
        TF1 *intWin = new TF1("intWindow", "[0]");
        intWin->SetLineColor(419);
        intWin->SetLineWidth(3);
        IntegrationWindow iw;
        std::cout << "Integration Windows:\nCh\tstart\tstop\n";
        // Fit the samples if reasonable
        for (int board = 0; board < gCONFIG.nChannelsPerBoard.size(); ++board)
        {
            for (int ch = 0; ch < gCONFIG.nChannelsPerBoard.at(board); ++ch)
            {
                hSample = hSampleSig.at(board)[ch];
                hSample->GetXaxis()->SetRange(2, 0.9 * SAMPLES_PER_WAVEFORM);
                peak = hSample->GetMinimum();
                peakBin = hSample->GetMinimumBin();
                hSample->GetXaxis()->SetRange(1, SAMPLES_PER_WAVEFORM);
                getPedestal(hSample, pedestal, pedStdv);
                Debug << "peak : " << peak << " bin: " << peakBin << std::endl;
                if (peak < pedestal - 5 * pedStdv)
                {
                    // There is at least a five sigma excess
                    // This channel probably collected a waveform
                    decay->FixParameter(0, peak);
                    decay->FixParameter(1, peakBin);
                    decay->SetParameter(2, 10);
                    decay->SetRange(peakBin, 1000);
                    hSample->Fit(decay, "RQ0");
                    hSample->GetFunction("decay")->ResetBit(TF1::kNotDraw);
                    gaus->FixParameter(0, peak);
                    gaus->FixParameter(1, peakBin);
                    gaus->SetParameter(2, 10);
                    gaus->SetRange(1, peakBin);
                    hSample->Fit(gaus, "RQ0+");
                    hSample->GetFunction("gaus")->ResetBit(TF1::kNotDraw);
                    stdv = fabs(gaus->GetParameter(2));
                    tau = fabs(decay->GetParameter(2));
                    iw.start = peakBin - gCONFIG.intRise * stdv;
                    iw.stop = peakBin + gCONFIG.intDecay * tau;
                    if (iw.start < 0)
                        iw.start = 0;
                    if (iw.stop > SAMPLES_PER_WAVEFORM - 1)
                        iw.stop = SAMPLES_PER_WAVEFORM - 1;
                }
                else
                {
                    // There is probably no signal - integrate over the whole time window
                    iw.start = 0;
                    iw.stop = SAMPLES_PER_WAVEFORM - 1;
                }
                intWin->SetRange(iw.start, iw.stop);
                intWin->FixParameter(0, pedestal);
                if (pedestal != 0)
                {
                    hSample->Fit(intWin, "RQ0+");
                    hSample->GetFunction("intWindow")->ResetBit(TF1::kNotDraw);
                }
                gCONFIG.integrationWindows.push_back(iw);
                std::printf("%d\t%d\t%d\n", ch, iw.start, iw.stop);
            }
        }

        // Deleting the functions
        gaus->Delete();
        decay->Delete();
        intWin->Delete();
    }

    // Clean up
    for (std::vector<float *> ch : wfData)
    {
        for (float *f : ch)
            delete[] f;
    }
    for (int board = 0; board < gCONFIG.nChannelsPerBoard.size(); ++board)
    {
        for (int ch = 0; ch < gCONFIG.nChannelsPerBoard.at(board); ++ch)
        {
            hSampleSig.at(board)[ch]->Write("", TObject::kOverwrite);
        }
        delete[] hSampleSig.at(board);
    }

    // Return 1 if eof was reached, 0 when nSampleEvents were read
    file->seekg(curPos);
    return (nEvents < gCONFIG.nSampleEvents);
}

/**
 * @brief Reads an event
 *
 * @details Reads an event in parallel from any channel available.
 *
 * @param file Input file, opened with Initialise() and the current position at the beginning of the EHDR
 * @param eh Event header at which the event will be read
 * @param wfData Waveform data at which the waveform data will be read
 * @param tCell Trigger cell
 * @return 0 on success, 1 for EOF, 2 otherwise
 */
int ReadAnEvent(std::ifstream *file, EventHeader &eh, std::vector<std::vector<float *>> &wfData, std::vector<std::vector<unsigned short>> *tCell)
{
    // Input: file - std::ifstream opened with Initialise() and the current position at the beginning of the EHDR.
    //        eh - to this variable the event header will be read
    //        wfData - to this variable the waveform data will be read. Memory already properly alocated!
    //        tCell - if valid address, the trigger cells will be stored to this address.
    // Return: 0 on success, 1 for eof, 2 otherwise
    //         the std::ifstream file will be at the very end of the event - next thing read should be the new EHDR.
    // Further: The details of the read data depends on gCONFIG.nChannelsPerBoard.

    char word[5];                                  // A default word to read from a file
    word[4] = '\0';                                // 4 bytes to read + 1 \0 character to terminate for printing purposes
    unsigned short voltages[SAMPLES_PER_WAVEFORM]; // Voltages for a single channel to read from the file.
    unsigned short index = 0;                      // Current channel number;

    // Get the last position in the file
    static long lastPos = 0;
    static long eventSize = 0;
    if (gCONFIG.firstOfRun || not lastPos || not eventSize)
    {
        long curPos = file->tellg();
        file->seekg(0, file->end);
        lastPos = file->tellg();
        file->seekg(curPos);
        gCONFIG.firstOfRun = false;
    }

    // Read event header first
    file->read((char *)&eh, sizeof(eh));
    if (std::memcmp(eh.tag, "EHDR", 4) != 0)
    {
        std::cerr << "!! No valid event header found. Found " << eh.tag << " instead" << std::endl;
        return 2;
    }
    DEBUG << "Reading event " << eh.serialNumber << std::endl;

    // Looping over all channels and all boards
    int totCh = 0;
    for (unsigned int board = 0; board < gCONFIG.nChannelsPerBoard.size(); ++board)
    {
        file->read(word, 4);
        if (std::memcmp(word, "B#", 2) != 0)
        {
            std::cerr << "No valid board header found. Found " << word << " instead" << std::endl;
            return 2;
        }
        DEBUG << " -found data for board " << *(short *)(word + 2) << std::endl;
        if (gCONFIG.runMode == 0)
        {
            // DRS encoding
            // Read trigger cell
            file->read(word, 4);
            if (std::memcmp(word, "T#", 2) != 0)
            {
                std::cerr << "No valid trigger cell found. Found " << word << " instead" << std::endl;
                return 2;
            }
            DEBUG << "Trigger cell: " << *(unsigned short *)(word + 2) << std::endl;
            if (tCell)
            {
                for (unsigned int i = 0; i < tCell->at(board).size(); ++i)
                    tCell->at(board).at(i) = *(unsigned short *)(word + 2);
            }
        }
        file->read(word, 4);
        while (word[0] == 'C')
        {
            ++totCh;
            char iCh[10];
            std::memcpy(iCh, word + 1, 3);
            iCh[3] = 0;
            index = std::stoi(iCh);
            DEBUG << "Found data for channel " << word << std::endl;

            // Read scaler - and ignore it for now
            file->read(word, 4);
            if (gCONFIG.runMode == 1)
            {
                // Read trigger cell
                file->read(word, 4);
                if (std::memcmp(word, "T#", 2) != 0)
                {
                    std::cerr << "No valid trigger cell found. Found " << word << " instead" << std::endl;
                    return 2;
                }
                DEBUG << "Trigger cell: " << *(unsigned short *)(word + 2) << std::endl;
                if (tCell)
                    tCell->at(board).at(index) = *(unsigned short *)(word + 2);
            }

            // Read voltages from file
            int curPos = file->tellg();
            file->read((char *)voltages, sizeof(voltages));
            std::cout << index << std::endl;
            for (int bin = 0; bin < SAMPLES_PER_WAVEFORM; ++bin)
            {
                wfData.at(board).at(index)[bin] = (voltages[bin] / 65536. + eh.rangeCenter / 1000. - 0.5);
            }
            if (!file->eof())
            {
                file->read(word, 4);
            }
        }
    }

    for (int board = 0; board < gCONFIG.nChannelsPerBoard.size(); ++board)
    {
        for (int ch = 0; ch < gCONFIG.nChannelsPerBoard.at(board); ++ch)
        {
            if (wfData.at(board).at(ch)[0] == 0)
            {
                for (int bin = 0; bin < SAMPLES_PER_WAVEFORM; ++bin)
                    wfData.at(board).at(ch)[bin] = -1.;
            }
        }
    }

    // if (std::memcmp(word, "EHDR", 4) != 0)
    // {
    //     std::cerr << "No valid channel header or event header found. Found " << word << " instead" << std::endl;
    //     return 2;
    // }
    eventSize = 0;
    if (gCONFIG.runMode == 1)
    {
        // WDB event encoding
        eventSize += sizeof(EventHeader);                  // 1 Event header
        eventSize += 4 * gCONFIG.nChannelsPerBoard.size(); // 4 byte board header per board

        // 4 Byte channel header, 4 byte scaler, 4 byte Trigger header + data
        eventSize += (12 + sizeof(voltages)) * totCh;
    }
    else
    {
        // DRS event encoding
        eventSize += sizeof(EventHeader);                  // 1 Event header
        eventSize += 8 * gCONFIG.nChannelsPerBoard.size(); // 4 byte board header per board + trigger header

        // 4 Byte channel header, 4 byte scaler + data
        eventSize += (8 + sizeof(voltages)) * totCh;
    }

    DEBUG << "File status: " << file->tellg() << "/" << lastPos << "  Event Size: " << eventSize << std::endl;

    // Assert that the file is not at its end
    if (lastPos - file->tellg() < eventSize)
    {
        DEBUG << "File is at its end: " << file->tellg() << "/" << lastPos << "   Event Size: " << eventSize << std::endl;
        return 1;
    }
    file->seekg(-4, file->cur);
    return 0;
}

/**
 * @brief Read the whole file
 *
 * @details Read the whole .dat file and write the values to a tree. The tree gets saved to
 * the file gCONFIG.theFile.
 *
 * @param file Opened .dat file. The current position is at the first event header.
 * @return 0 on success, 1 otherwise
 */
int ReadFile(std::ifstream *file)
{
    int retVal = 0; ///< The value to be returned

    // Initialise the tree
    //  1. Number of channels and leaflist description
    int nChannelsTot = 0;
    std::string leaflistDescription = "ch0";
    for (int nChannels : gCONFIG.nChannelsPerBoard)
    {
        nChannelsTot += nChannels;
    }
    for (int ch = 1; ch < nChannelsTot; ++ch)
    {
        leaflistDescription += ":ch";
        leaflistDescription += std::to_string(ch);
    }
    Debug << leaflistDescription << std::endl;

    // 2. Variables to fill the tree from
    int fRun = gCONFIG.run;
    int fEvent = 0;
    float fTimestamp = 0;
    float *fTime = new float[nChannelsTot];
    float *fTimeLE = new float[nChannelsTot];
    float *fAmplitude = new float[nChannelsTot];
    float *fArea = new float[nChannelsTot];
    float *fArea2 = new float[nChannelsTot];
    float *fPed = new float[nChannelsTot];
    float *fStdv = new float[nChannelsTot];
    float *fThreshold = new float[nChannelsTot];
    float *fStdvRaw = new float[nChannelsTot];
    float *fSineAmplitude = new float[nChannelsTot];
    float *fSinePhase = new float[nChannelsTot];

    // 3. TTree itself
    TTree *theTree = new TTree("T", "WaveDREAM Tree");
    theTree->SetDirectory(gCONFIG.theFile);
    theTree->Branch("runnumber", &fRun);
    theTree->Branch("event", &fEvent);
    theTree->Branch("timestamp", &fTimestamp);
    theTree->Branch("time", fTime, leaflistDescription.c_str());
    theTree->Branch("timeLE", fTimeLE, leaflistDescription.c_str());
    theTree->Branch("amplitude", fAmplitude, leaflistDescription.c_str());
    theTree->Branch("area", fArea, leaflistDescription.c_str());
    theTree->Branch("area2", fArea2, leaflistDescription.c_str());
    theTree->Branch("ped", fPed, leaflistDescription.c_str());
    theTree->Branch("Stdv", fStdv, leaflistDescription.c_str());
    theTree->Branch("Threshold", fThreshold, leaflistDescription.c_str());

    if (gCONFIG.subtractSine)
    {
        theTree->Branch("StdvRaw", fStdvRaw, leaflistDescription.c_str());
        theTree->Branch("SineAmplitude", fSineAmplitude, leaflistDescription.c_str());
        theTree->Branch("SinePhase", fSinePhase, leaflistDescription.c_str());
    }
    Debug << "TTree initialised " << std::endl;
    if (gCONFIG.debug)
        theTree->Print();

    // Start the mainloop
    Debug << "Mainloop over the file started " << std::endl;
    EventHeader eh;
    std::vector<std::vector<float *>> wfData;
    std::vector<std::vector<float *>> wfTime;
    std::vector<std::vector<unsigned short>> tCell;
    for (int board = 0; board < gCONFIG.nChannelsPerBoard.size(); board++)
    {
        wfData.push_back({});
        wfTime.push_back({});
        tCell.push_back({});
        for (int ch = 0; ch < gCONFIG.nChannelsPerBoard.at(board); ++ch)
        {
            wfData.at(board).push_back(new float[SAMPLES_PER_WAVEFORM]);
            wfTime.at(board).push_back(new float[SAMPLES_PER_WAVEFORM]);
            tCell.at(board).push_back(0);
        }
    }

    float t1, t2, dt;
    int index = 0;

    // Variables to be later filled into the tree
    float timeCF, timeLE, ampl, area, area2, ped, stdv, stdvRaw, sinAmpl, phi, thr;
    bool isSignal = false;
    int minBin = 0;
    float firstMin = 0;

    while (not file->eof())
    {
        retVal = ReadAnEvent(file, eh, wfData, &tCell);
        if (retVal)
            break;
        if (eh.serialNumber % 100 == 0)
        {
            std::cout << "Processing Event " << eh.serialNumber << std::endl;
        };

        // Calculate times for each channel
        for (int board = 0; board < gCONFIG.nChannelsPerBoard.size(); ++board)
        {
            for (int ch = 0; ch < gCONFIG.nChannelsPerBoard.at(board); ++ch)
            {
                int bin = tCell.at(board).at(ch);
                double time = 0;
                for (int i = 0; i < SAMPLES_PER_WAVEFORM; ++i)
                {
                    wfTime.at(board).at(ch)[i] = time;
                    time += gCONFIG.timeBinWidth.at(board).at(ch)[bin];
                    ++bin;
                    bin %= 1024;
                }
            }
        }

        // Align cell #0 of all channels
        index = 0;
        for (int board = 0; board < gCONFIG.nChannelsPerBoard.size(); ++board)
        {
            t1 = wfTime.at(board).at(index)[(1024 - tCell.at(board)[index]) % 1024];
            ++index;
            for (int ch = 1; ch < gCONFIG.nChannelsPerBoard.at(board); ++ch)
            {
                t2 = wfTime.at(board).at(index)[(1024 - tCell.at(board)[index]) % 1024];
                dt = t1 - t2;
                for (int i = 0; i < SAMPLES_PER_WAVEFORM; ++i)
                {
                    wfTime.at(board).at(index)[i] += dt;
                }
                ++index;
            }
        }

        // Analysis of individual channels
        int nCh = 0;
        for (int board = 0; board < gCONFIG.nChannelsPerBoard.size(); ++board)
        {
            for (int ch = 0; ch < gCONFIG.nChannelsPerBoard.at(board); ++ch)
            {
                float *aWF = wfData.at(board).at(ch);
                float *aWFT = wfTime.at(board).at(ch);

                // Flip waveform if needed -> want a negative one
                if (gCONFIG.sigWF == 1)
                {
                    for (int i = 0; i < SAMPLES_PER_WAVEFORM; ++i)
                        aWF[i] *= -1;
                    Debug << "Positive Waveform - flipping" << std::endl;
                }

                // Subtract sine noise if activated
                if (gCONFIG.subtractSine && abs(aWF[0]) != 1.)
                {
                    getPedestal(aWF, ped, stdvRaw);
                    if (eh.serialNumber < gCONFIG.nSaveEvents)
                    {
                        char name[32];
                        sprintf(name, "event%03d_ch%03d_raw", eh.serialNumber, ch);
                        subtractSineNoise(aWF, aWFT, sinAmpl, phi, name);
                    }
                    else
                    {
                        subtractSineNoise(aWF, aWFT, sinAmpl, phi);
                    }
                }

                // Get the pedestal
                getPedestal(aWF, ped, stdv);

                // Check for a signal waveform - aka 5 sigma excess
                // May needs a somewhat less restrictive cut
                isSignal = false;
                IntegrationWindow iw = gCONFIG.integrationWindows.at(nCh);
                for (int i = iw.start; i < iw.stop && not isSignal; ++i)
                {
                    if ((ped - aWF[i]) > 5 * stdv)
                        isSignal = true;
                }

                // Fill a graph if needed:
                if (eh.serialNumber < gCONFIG.nSaveEvents && abs(ped) != 1.)
                {
                    TGraph *aGraph = new TGraph(SAMPLES_PER_WAVEFORM, aWFT, aWF);
                    TF1 *pedestal = new TF1("pedestal", "[0]", aWFT[iw.start], aWFT[iw.stop]);
                    pedestal->SetParameter(0, ped);
                    aGraph->GetListOfFunctions()->Add(pedestal);
                    TString name = TString::Format("event%03d_ch%03d", eh.serialNumber, ch);
                    aGraph->SetName(name);
                    aGraph->SetTitle(name);
                    aGraph->Write("", TObject::kOverwrite);
                }
                std::cout << "ped:" << ped << "stdv: " << stdv << std::endl;
                if (isSignal && abs(ped) != 1.)
                {
                    // Find peak
                    ampl = firstMin = 0;
                    for (int bin = iw.start + 1; bin < iw.stop - 2; ++bin)
                    {
                        if (ampl > aWF[bin] - ped)
                        {
                            // The maximal amplitude in the intergation window
                            ampl = aWF[bin] - ped;
                        }
                        if (firstMin > -gCONFIG.leThreshold && aWF[bin] < aWF[bin - 1] && aWF[bin] < aWF[bin + 1] && aWF[bin] < aWF[bin + 2])
                        {
                            firstMin = aWF[bin] - ped;
                            minBin = bin;
                        }
                    }

                    // CF + LE time extraction
                    // start at minBin where the first (local) min has been found.
                    // Then move to earlier times to find the closest crossing.
                    timeLE = -1e+5;
                    timeCF = -1e+5;
                    for (int bin = minBin; bin > iw.start; --bin)
                    {
                        // DEBUG << aWF[bin] - ped << "  " << -gCONFIG.leThreshold << " " << gCONFIG.cfFraction * firstMin << std::endl;
                        // LE
                        if (timeLE < 0 && aWF[bin] - ped > -gCONFIG.leThreshold)
                        {
                            std::cout << "gCONFIG.leThreshold: " << gCONFIG.leThreshold << std::endl;
                            timeLE = (-gCONFIG.leThreshold - aWF[bin] + ped) / (aWF[bin + 1] - aWF[bin]) * (aWFT[bin + 1] - aWFT[bin]) + aWFT[bin];
                            std::cout << "ciao ciao" << std::endl;
                        }

                        // Constant Fraction
                        if (timeCF < 0 && aWF[bin] - ped > gCONFIG.cfFraction * firstMin)
                        {
                            timeCF = (gCONFIG.cfFraction * firstMin - aWF[bin] + ped) / (aWF[bin + 1] - aWF[bin]) * (aWFT[bin + 1] - aWFT[bin]) + aWFT[bin];
                        }
                        if (timeCF > 0 && timeLE > 0)
                            break;
                    }

                    // Charge Integration
                    area = 0;
                    area2 = 0;
                    for (int bin = iw.start; bin < iw.stop; ++bin)
                    {
                        area += ped - aWF[bin];
                        area2 += (ped - aWF[bin]) * gCONFIG.timeBinWidth.at(board).at(ch)[(tCell.at(board).at(ch) + bin) % 1024];
                    }
                }
                else
                {
                    timeCF = -1e5;
                    timeLE = -1e5;
                    ampl = 0;
                    area = 0;
                    area2 = 0;
                    thr = 0;
                }
                DEBUG << ch << "/" << nChannelsTot << std::endl;
                fEvent = eh.serialNumber;
                fTimestamp = getTimeStamp(eh);
                fTime[nCh] = timeCF;
                fTimeLE[nCh] = timeLE;
                fAmplitude[nCh] = -ampl; // prefer a positive amplitude in the end;
                fArea[nCh] = area;
                fArea2[nCh] = area2;
                fPed[nCh] = ped;
                fStdv[nCh] = stdv;
                fThreshold[nCh] = thr;
                fStdvRaw[nCh] = stdvRaw;
                fSineAmplitude[nCh] = sinAmpl;
                fSinePhase[nCh] = phi;
                nCh++;
            }
            theTree->Fill();
        }
    }

    // Clean up at the end of the file
    Debug << "Finished reading - clean up " << std::endl;
    theTree->Write("", TObject::kOverwrite);
    delete[] fTime;
    delete[] fTimeLE;
    delete[] fAmplitude;
    delete[] fArea;
    delete[] fArea2;
    delete[] fPed;
    delete[] fStdv;
    delete[] fThreshold;
    delete[] fStdvRaw;
    delete[] fSineAmplitude;
    delete[] fSinePhase;
    for (int board = 0; board < gCONFIG.nChannelsPerBoard.size(); ++board)
    {
        for (float *f : wfData.at(board))
            delete[] f;
    }
    return 0;
}

/**
 * @brief
 *
 * @return int
 */
int Config()
{
    int ret = 0;
    enum position
    {
        posNSAMPLE,
        posNSAVE,
        posSIGN,
        posCF,
        posLE,
        posRAISE,
        posDECAY,
        posNOISE,
        posFREQ,
        posDEBUG,
        posCONTINUE,
        posCANCEL
    };
    int pos = posCONTINUE;
    int key = 0;
    const int labelWidth = 45;

    // Initialise the screen
    initscr();
    nl();
    noecho();
    cbreak();
    keypad(stdscr, TRUE);
    while (key != 0xA || pos < posCONTINUE)
    {
        clear();
        curs_set(0);

        // Print some default text
        mvaddstr(1, 4, "***********************************************************************");
        mvaddstr(2, 4, "*                                                                     *");
        mvaddstr(3, 4, "*           Welcome to the WaveDREAM data analysis tool               *");
        mvprintw(4, 4, "*                  Version: %-30s            *\n", __READ_VERSION__);
        mvaddstr(5, 4, "*                      Configuration Screen                           *");
        mvaddstr(6, 4, "***********************************************************************");

        mvaddstr(8 + posNSAMPLE, 4, "nSamples for Integration Windows:");
        mvaddstr(8 + posNSAVE, 4, "Number of Events to save Waveforms:");
        mvaddstr(8 + posSIGN, 4, "Waveform type:");
        mvaddstr(8 + posCF, 4, "CF Fraction:");
        mvaddstr(8 + posLE, 4, "LE Fraction:");
        mvaddstr(8 + posRAISE, 4, "Integration rise (stdv):");
        mvaddstr(8 + posDECAY, 4, "Integration decay (tau):");
        mvaddstr(8 + posNOISE, 4, "Sine Noise Subtraction:");
        mvaddstr(8 + posFREQ, 4, "Sine Noise Frequency (MHz):");
        mvaddstr(8 + posDEBUG, 4, "Debug Mode:");

        // Print the values. The one at the current position (pos) is supposed to stand out
        if (pos == posNSAMPLE)
            standout();
        mvprintw(8 + posNSAMPLE, labelWidth, "%d", gCONFIG.nSampleEvents);
        standend();
        if (pos == posNSAVE)
            standout();
        mvprintw(8 + posNSAVE, labelWidth, "%d", gCONFIG.nSaveEvents);
        standend();
        if (pos == posSIGN)
            standout();
        mvprintw(8 + posSIGN, labelWidth, "%+d - %s", gCONFIG.sigWF, gCONFIG.sigWF == 1 ? "positive" : "negative");
        standend();
        if (pos == posCF)
            standout();
        mvprintw(8 + posCF, labelWidth, "%1.2f", gCONFIG.cfFraction);
        standend();
        if (pos == posLE)
            standout();
        mvprintw(8 + posLE, labelWidth, "%1.2f", gCONFIG.leThreshold);
        standend();
        if (pos == posRAISE)
            standout();
        mvprintw(8 + posRAISE, labelWidth, "%1.2f", gCONFIG.intRise);
        standend();
        if (pos == posDECAY)
            standout();
        mvprintw(8 + posDECAY, labelWidth, "%1.2f", gCONFIG.intDecay);
        standend();
        if (pos == posNOISE)
            standout();
        mvprintw(8 + posNOISE, labelWidth, "%s", gCONFIG.subtractSine ? "true" : "false");
        standend();
        if (pos == posFREQ)
            standout();
        mvprintw(8 + posFREQ, labelWidth, "%1.2f", gCONFIG.noiseFrequency * 1e-6 / TMath::TwoPi());
        standend();
        if (pos == posDEBUG)
            standout();
        if (gCONFIG.debug == 0)
        {
            mvaddstr(8 + posDEBUG, labelWidth, "0 - No debug output");
        }
        else if (gCONFIG.debug == 1)
        {
            mvaddstr(8 + posDEBUG, labelWidth, "1 - Normal debug output");
        }
        else
        {
            mvaddstr(8 + posDEBUG, labelWidth, "2 - Extensive debug output");
        }
        standend();
        if (pos == posCONTINUE)
            standout();
        mvaddstr(8 + posCANCEL, 4, "continue");
        standend();
        if (pos == posCANCEL)
            standout();
        mvaddstr(8 + posCANCEL, 20, "abort");
        standend();

        // Get it interacting, registering key strokes
        key = getch();
        // 0xA = Enter
        while (key != 0xA && key != KEY_DOWN && key != KEY_UP && key != KEY_LEFT && key != KEY_RIGHT)
            key = getch();
        if (key == KEY_DOWN && pos < posCONTINUE)
        {
            // Down arrow pressed, move to the next entry
            pos += 1;
        }
        else if (key == KEY_UP && pos > 0)
        {
            // Up arrow pressed, move to the previous entry
            pos -= 1;
        }
        else if (key == KEY_LEFT)
        {
            // Left arrow pressed, toggle entries where possible
            if (pos == posSIGN)
            {
                gCONFIG.sigWF = -1;
            }
            else if (pos == posDEBUG && gCONFIG.debug > 0)
            {
                gCONFIG.debug -= 1;
            }
            else if (pos == posNOISE)
            {
                gCONFIG.subtractSine = false;
            }
            else if (pos == posCANCEL)
            {
                pos = posCONTINUE;
            }
        }
        else if (key == KEY_RIGHT)
        {
            // Right arrow pressed, toggle entries where possible
            if (pos == posSIGN)
            {
                gCONFIG.sigWF = 1;
            }
            else if (pos == posNOISE)
            {
                gCONFIG.subtractSine = true;
            }
            else if (pos == posDEBUG && gCONFIG.debug < 2)
            {
                gCONFIG.debug += 1;
            }
            else if (pos == posCONTINUE)
            {
                pos = posCANCEL;
            }
        }
        else if (key == 0xA && pos < posCONTINUE && pos != posDEBUG && pos != posSIGN && pos != posNOISE)
        {
            // Enter key pressed - Get to entry mode
            move(8 + pos, labelWidth);
            clrtoeol();  // clear what is currently there
            curs_set(1); // enable cursor
            std::vector<int> nr;
            while ((key = getch()) != 0xA && key != 0x1B)
            {
                // 0xA -> Enter, 0x1B -> Esc
                int n = key - '0';
                if (n > -1 && n < 10)
                {
                    // A number was typed, add it to memory and to the screen
                    nr.push_back(n);
                    addch(key);
                }
                else if (key == '.')
                {
                    // A decimal point. Code as -1 in the vector
                    bool hasDecimal = false;
                    for (unsigned int i = 0; i < nr.size() && not hasDecimal; ++i)
                        hasDecimal = (nr.at(i) == -1);
                    if (not hasDecimal)
                    {
                        nr.push_back(-1);
                        addch(key);
                    }
                }
                else if (key == 127)
                {
                    // Backspace - drop the last entry from memory, move a step back and delete last char printed
                    int x, y;
                    getyx(stdscr, y, x);
                    move(y, x - 1);
                    clrtoeol();
                    nr.pop_back();
                }
            }
            if (key == 0xA)
            {
                // if terminated with the Enter key, set the new value
                float newVal = 0;
                int exponent = 0;
                for (unsigned int i = 0; i < nr.size(); ++i)
                {
                    int digit = nr.at(nr.size() - 1 - i);
                    if (digit != -1)
                    {
                        newVal += nr.at(nr.size() - 1 - i) * pow(10, exponent);
                        ++exponent;
                    }
                    else if (digit == -1)
                    {
                        newVal *= pow(10, -exponent);
                        exponent = 0;
                    }
                }
                if (pos == posNSAMPLE)
                {
                    gCONFIG.nSampleEvents = newVal;
                }
                else if (pos == posNSAVE)
                {
                    gCONFIG.nSaveEvents = newVal;
                }
                else if (pos == posCF)
                {
                    gCONFIG.cfFraction = newVal;
                }
                else if (pos == posLE)
                {
                    gCONFIG.leThreshold = newVal;
                }
                else if (pos == posRAISE)
                {
                    gCONFIG.intRise = newVal;
                }
                else if (pos == posDECAY)
                {
                    gCONFIG.intDecay = newVal;
                }
                else if (pos == posFREQ)
                {
                    gCONFIG.noiseFrequency = newVal * 1e+6 * TMath::TwoPi();
                }
            }
        }
    }

    // Restore the usual terminal configuration
    nocbreak();
    endwin();
    if (pos == posCONTINUE)
    {
        // Enter was hit, while "continue" was active
        ret = 0;
    }
    else if (pos == posCANCEL)
    {
        // Enter was hit, while "abort" was active
        ret = 1;
    }
    return ret;
}

/**
 * @brief
 *
 */
void PrintHelp()
{
    // Print a short help text to the screen
    std::cout << " Use to read data from binary WaveDREAM output files to root files\n\n";
    std::cout << " ./readWD [args] [-o {outputPath}] file1.dat [file2.dat ... fileN.dat]\n\n";
    std::cout << " Reads the files file1.dat ... fileN.dat and creates a root file for each\n\n";
    std::cout << "   -o {outputPath}        set the path where the rootfiles are created\n";
    std::cout << "                          default is ./\n\n";
    std::cout << " Possible [args] are : \n";
    std::cout << "   -d, --DEBUG            print debug information\n";
    std::cout << "   -f, --ForceOverwrite   enforce the overwriting of existing rootfiles\n";
    std::cout << "   -h, --help             print this help text and quit\n";
    std::cout << "   -p, --pos_wf           use for positive waveforms\n";
    std::cout << "   -c, --config           enter configuration interface, requires user input\n";
    std::cout << "   -s, --subtractNoise    subtract sine noise from waveforms\n";
    std::cout << std::endl;
}

/**
 * @brief Main function
 *
 * @param argc
 * @param argv
 * @return int
 */
int main(int argc, const char **argv)
{
    // Main function called at start of the program.
    // Deal with user input and call other functions accordingly
    std::printf("***********************************************************************\n");
    std::printf("*                                                                     *\n");
    std::printf("*           Welcome to the WaveDREAM data analysis tool               *\n");
    std::printf("*                  Version: %-32s          *\n", __READ_VERSION__);
    std::printf("*                                                                     *\n");
    std::printf("***********************************************************************\n");
    std::cout << std::endl;

    // Some configuration flags that only take effect in the main function
    bool overwrite = false;
    bool enterConfig = false;
    std::vector<std::string> fileList;
    std::string outputPath("./");

    // Reading the command line arguments
    for (int i = 1; i < argc; ++i)
    {
        if (strcmp(argv[i], "-d") == 0 || strcmp(argv[i], "--DEBUG") == 0)
        {
            gCONFIG.debug = 1;
        }
        else if (strcmp(argv[i], "-D") == 0)
        {
            gCONFIG.debug = 2;
        }
        else if (strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--config") == 0)
        {
            enterConfig = true;
        }
        else if (strcmp(argv[i], "-p") == 0 || strcmp(argv[i], "--pos_wf") == 0)
        {
            gCONFIG.sigWF = 1;
        }
        else if (strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "--ForceOverwrite") == 0)
        {
            overwrite = true;
        }
        else if (strcmp(argv[i], "-s") == 0 || strcmp(argv[i], "--subtractNoise") == 0)
        {
            gCONFIG.subtractSine = true;
        }
        else if (strcmp(argv[i], "-o") == 0)
        {
            outputPath = argv[++i];
        }
        else if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0)
        {
            PrintHelp();
            return 1;
        }
        else if (argv[i][0] == '-')
        {
            std::cerr << "Unkown flag: " << argv[i] << "\n\n";
            PrintHelp();
            return 1;
        }
        else
        {
            // Anything that does not start on '-' is taken as a file to read
            fileList.push_back(argv[i]);
        }
    }
    if (enterConfig)
    {
        if (Config())
            return 1;
    }

    // Assert that at least one file is given
    if (fileList.size() == 0)
    {
        std::cerr << "!! No input file is given\n"
                  << std::endl;
        PrintHelp();
        return 1;
    }

    // Return value for various functions is stored here.
    // If different from 0, something went wrong and the file seems faulty
    int returnValue = 0;
    for (unsigned int i = 0; i < fileList.size(); ++i)
    {
        // Determining file paths and names
        std::string inFileName(fileList.at(i));
        std::string outFileName(fileList.at(i));
        std::cout << std::endl;
        std::cout << "==========================================================\n";
        std::cout << "  Reading from file: " << inFileName << std::endl;

        // Keep only file name for the output file name
        size_t pos = outFileName.rfind('/');
        if (pos != std::string::npos)
            outFileName.replace(0, pos + 1, "");

        // Replace the .dat by a .root
        pos = outFileName.rfind('.');
        if (pos != std::string::npos)
        {
            outFileName.replace(pos + 1, outFileName.length(), "root");
        }
        else
        {
            // The file actually did not have any extension - try by adding .root extension
            outFileName += ".root";
        }

        // Extract runnumber:
        gCONFIG.run = 0;
        int c = outFileName[outFileName.length() - 6] - '0';
        int j = 0;
        while (c > -1 && c < 10)
        {
            // It is a numeric char
            gCONFIG.run += c * pow(10, j);
            c = outFileName[outFileName.length() - 6 - ++j] - '0';
        }

        // The output file is to be located at the output path
        outFileName = outputPath + outFileName;
        std::cout << "  Writing to file: " << outFileName << std::endl;
        std::cout << "  Run number : " << gCONFIG.run << std::endl;
        std::cout << "==========================================================" << std::endl;

        // Check accessibility of the output file:
        if (!overwrite && !gSystem->AccessPathName(outFileName.c_str()))
        {
            std::cout << "The file " << outFileName << " already exists.\n";
            std::cout << "Use the -f flag to force overwriting existing files." << std::endl;
            continue; // The file already exists and may not be overwritten.
        }

        // Initialise the input file stream
        std::ifstream *inFile = Initialise(inFileName);
        if (not inFile)
        {
            // Something went wrong
            continue;
        }

        // Initialise the root file
        gCONFIG.theFile = TFile::Open(outFileName.c_str(), "recreate");
        if (!gCONFIG.theFile)
        {
            // For some reason, the file could not be opened.
            std::cerr << "Could not open file " << outFileName.c_str() << std::endl;
            Debug << "... cleaning up." << std::endl;
            inFile->close();
            delete inFile;
            for (unsigned int board = 0; board < gCONFIG.timeBinWidth.size(); ++board)
            {
                for (float *arr : gCONFIG.timeBinWidth.at(board))
                    delete[] arr;
            }
            continue;
        }
        returnValue = GetIntegrationBounds(inFile);
        if (returnValue == 2)
            return 2;
        if (gCONFIG.debug)
        {
            if (returnValue)
            {
                Debug << "EOF reached" << std::endl;
            }
            else
            {
                Debug << gCONFIG.nSampleEvents << " read" << std::endl;
            }
        }
        returnValue = ReadFile(inFile);
        gCONFIG.theFile->Close();
        inFile->close();
        for (unsigned int board = 0; board < gCONFIG.timeBinWidth.size(); ++board)
        {
            for (float *arr : gCONFIG.timeBinWidth.at(board))
                delete[] arr;
        }
        delete inFile;
    }
    return 0;
}
#include <Bela.h>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iomanip>

#include "DebouncedButton.h"

class BeatDetector {
    public:
        // Create default constructor
        BeatDetector() {}

        BeatDetector(BelaContext *context) : context(context){
            // Create a button object for each physical button
            debouncedButton1 = DebouncedButton(P8_07, 0.1*context->audioSampleRate);
            debouncedButton2 = DebouncedButton(P8_08, 0.1*context->audioSampleRate);
            // Use the first second after the activation button is pressed for
            // system calibration
            calibrationTime = 1*context->audioSampleRate;
            // Declare hop and window size in samples for analysis
            hopSize = int(0.01*context->audioSampleRate);
            winSize = int(0.02*context->audioSampleRate);
            overlapFactor = ceil(winSize/hopSize);
            // Create circular buffer of size to store windows for processing
            // in real-time
            audioBufferSize = winSize+((overlapFactor-1)*hopSize);
            audioBuffer.resize(audioBufferSize, 0);

            // Create a write pointer at position 0 in buffer
            bufferWritePtr = 0;
            // Resize buffer to allow for a read pointer for each overlapping
            // window
            bufferReadPtrs.resize(overlapFactor);
            bufferReadPtrsPos.resize(overlapFactor);

            //energyFile.open("/ShannonEnergy.csv");
            //peaksFile.open("/Peaks.csv");

            // bufferReadPtrs stores the index to read in relation to the
            // circular buffer.
            // bufferReadPtrsPos stores the index position in relation to the
            // window currently being read

            // Equivelants to Numpy style vector operations (np.arange and *=)
            // taken from:
            // https://codereview.stackexchange.com/questions/77546/multiply-vector-elements-by-a-scalar-value-using-stl-and-templates
            // http://stackoverflow.com/questions/17694579/use-stdfill-to-populate-vector-with-increasing-numbers
            // Fill vector with values from 0-overlapFactor
            std::iota(std::begin(bufferReadPtrs), std::end(bufferReadPtrs), 0);
            // Multiply all values by -hopSize
            std::transform(bufferReadPtrs.begin(), bufferReadPtrs.end(), bufferReadPtrs.begin(), std::bind1st(std::multiplies<int>(), -hopSize));
            // Fill vector with values from 0-overlapFactor
            std::iota(std::begin(bufferReadPtrsPos), std::end(bufferReadPtrsPos), 0);
            // Multiply all values by -hopSize
            std::transform(bufferReadPtrsPos.begin(), bufferReadPtrsPos.end(), bufferReadPtrsPos.begin(), std::bind1st(std::multiplies<int>(), -hopSize));

            // Store 5 seconds of previous normalisaed shannon energy values
            normShanEngySize = round((1.0*context->audioSampleRate) / hopSize);
            normalisedShannonEnergy.resize(normShanEngySize);
            shannonEnergyPeaks.resize(normShanEngySize, false);

            nn = 100;
            shannonEnergy.resize(nn, 0);
        }


        bool calculateAverageShannonEnergy() {
            bool newVal = false;
            // for each read pointer
            for(int i=0; i<bufferReadPtrs.size(); i++) {
                // Get pointer positions in relation to circular buffer and to
                // the current window being read by each pointer
                int ptr = bufferReadPtrs[i];
                int pos = bufferReadPtrsPos[i];
                // If the current read pointer has reached the end of a
                // window...
                if(pos == winSize-1) {
                    float sigAccum = 0.0;
                    float powSamp;
                    // For window index up to the current pointer index...
                    for(int i=(ptr-(winSize-1))%audioBufferSize; i != ptr; i++, i%=audioBufferSize){
                        // Perform samplewise operations for shannon energy
                        powSamp = pow(audioBuffer[i], 2.0);
                        sigAccum += powSamp * log(powSamp);
                    }
                    sigAccum = (-1.0/winSize)*sigAccum;
                    // If a Nan value is created due to a log(0), set the
                    // output to 0
                    if(isnan(sigAccum)) {
                        sigAccum = 0;
                    }
                    // Calculate consecutive mean values by accumulating n
                    // previous values and dividing by n
                    shanEngyMeanAccum -= shannonEnergy[normShanEngyPtr];
                    // Store current shannon energy value to buffer
                    shannonEnergy[normShanEngyPtr] = sigAccum;
                    shanEngyMeanAccum += sigAccum;
                    shanEngyMean = shanEngyMeanAccum / shannonEnergy.size();
                    // Calculate standard deviation of previous n values
                    float stdAccum = 0.0;
                    for(int j=0; j<shannonEnergy.size(); j++) {
                        stdAccum += pow(std::abs(shannonEnergy[j] - shanEngyMean), 2.0);
                    }
                    float stdMean = stdAccum / shannonEnergy.size();
                    shanEngyStd = sqrt(stdMean);
                    //rt_printf("%f\n", shanEngyStd);

                    // Normalise shannon energy using mean and standard
                    // deviation
                    float normShanEngyVal = (sigAccum - shanEngyMean)/shanEngyStd;
                    // If output value is Nan, set to 0
                    if(isnan(normShanEngyVal)) {
                        normShanEngyVal = 0;
                    }
                    // Save final value to shannon energy buffer
                    normalisedShannonEnergy[normShanEngyPtr] = normShanEngyVal;
                    //energyFile << std::fixed << std::setprecision(8) << normalisedShannonEnergy[normShanEngyPtr] << std::endl;


                    // Signal that a new shannon energy value has been
                    // calculated
                    newVal = true;
                }

            }
            return newVal;
        }


        float processSample(const float &sample, int n) {
            // Check if button has been pressed...
            if(debouncedButton1.getCurrentVal(context, n)) {
                sampleCounter = 0;
            }
            // If in the calibration phase of processing, run calibration and
            // return 0
            if(sampleCounter < calibrationTime){
                calibrateSystem(sample);
                //rt_printf("%d\n", sampleCounter);
                sampleCounter += 1;
                return 0.0;
            }


            /*
            if(sampleCounter > calibrationTime*17) {
                //rt_printf("YUP");
                //energyFile.close();
                //peaksFile.close();
                exit(0);
            }
            */

            audioBuffer[bufferWritePtr] = sample;
            // Calculate new values for shannon energy and get status of
            // whether a new value has been generated
            bool newVal = calculateAverageShannonEnergy();

            // Every time a new shannon energy value is calculated...
            if(newVal) {
                // Calculate if it is a peak above the threshold
                if(findNewPeaks() && (sampleCounter > 0.05 * context->audioSampleRate)){
                    rt_printf("BEEP %d\n", sampleCounter);
                    sampleCounter = 0;
                }

            }
            // Increment pointers and wrap around their respective
            // container sizes
            shanEngyPtr += 1;
            shanEngyPtr %= nn;
            normShanEngyPtr += 1;
            normShanEngyPtr %= normShanEngySize;
            // Increment the buffer write pointer and wrap to keep within size
            // of buffer
            bufferWritePtr += 1;
            bufferWritePtr %= audioBufferSize;

            transform(bufferReadPtrs.begin(), bufferReadPtrs.end(), bufferReadPtrs.begin(), bind2nd(std::plus<int>(), 1.0));
            transform(bufferReadPtrs.begin(), bufferReadPtrs.end(), bufferReadPtrs.begin(), bind2nd(std::modulus<int>(), audioBufferSize));
            transform(bufferReadPtrsPos.begin(), bufferReadPtrsPos.end(), bufferReadPtrsPos.begin(), bind2nd(std::plus<int>(), 1.0));
            transform(bufferReadPtrsPos.begin(), bufferReadPtrsPos.end(), bufferReadPtrsPos.begin(), bind2nd(std::modulus<int>(), audioBufferSize));
            // Increment the counter for the number of processed samples since
            // the button was pressed. This will eventually overflow but won't
            // be used by that point unless calibration is crazy long so it has
            // been left as is.
            sampleCounter += 1;
            return sample;
        }


        bool findNewPeaks() {
            if(normalisedShannonEnergy[normShanEngyPtr] > threshold) {
                // If a peak hasn't already been found for the current
                // overshoot...
                if(!peakFound) {
                    //////////////////////////////////////////////////////////////////////
                    // Search for a peak since the Pa values first raised above the
                    // threshold
                    //////////////////////////////////////////////////////////////////////
                    //peaksFile << int(1) << std::endl;

                    float diff = normalisedShannonEnergy[normShanEngyPtr]-normalisedShannonEnergy[(normShanEngyPtr-1)%normShanEngySize];
                    //float prevDiff = normalisedShannonEnergy[(normShanEngyPtr-1)%normalisedShannonEnergy.size()]-normalisedShannonEnergy[(normShanEngyPtr-2)%normalisedShannonEnergy.size()];
                    //rt_printf("%f %f\n", diff, prevDiff);
                    // Simplified plateu case handeling. The start of the
                    // plateu is classified as the peak
                    if(diff == 0) {
                        diff = -1.0;
                    }
                    if(prevDiff == 0) {
                        prevDiff = -1.0;
                    }
                    // Check if last frame is a peak
                    if((prevDiff < 0.0) && (diff > 0.0)) {
                        peakFound = true;
                        shannonEnergyPeaks[normShanEngyPtr] = true;
                    }
                    else {
                        peakFound = false;
                        shannonEnergyPeaks[normShanEngyPtr] = false;
                        prevDiff = diff;
                    }
                }
                else {
                    shannonEnergyPeaks[normShanEngyPtr] = false;
                    //peaksFile << int(0) << std::endl;
                }
            }
            else {
                peakFound = false;
                shannonEnergyPeaks[normShanEngyPtr] = false;
                //peaksFile << int(0) << std::endl;
                prevDiff = 1.0;
            }

            //peaksFile << int(shannonEnergyPeaks[normShanEngyPtr]) << std::endl;
            return peakFound;
        }
    private:
        void calibrateSystem(const float &sample) {
            if(sample > calibratedMax) {
                calibratedMax = sample;
            }
        }
        // Counts how many samples have elapsed since the last trigger
        unsigned int sampleCounter = 0;
        DebouncedButton debouncedButton1;
        DebouncedButton debouncedButton2;
        BelaContext* context;

        //////////////////////////////////////////////////////////////////////
        // Calibration variables
        //////////////////////////////////////////////////////////////////////
        // Sets the amount of time to be used for calibration
        int calibrationTime = 0;

        // Store th maximum value found in the calibration phase
        float calibratedMax = 0.0;

        //////////////////////////////////////////////////////////////////////
        // Sample buffering variables
        //////////////////////////////////////////////////////////////////////

        int hopSize;
        int winSize;
        int overlapFactor;
        // Buffer for storing samples until full windows are accumulated
        std::vector<float> audioBuffer;
        int audioBufferSize = 0;
        // Vector of read pointers for buffer
        std::vector<int> bufferReadPtrs;
        std::vector<int> bufferReadPtrsPos;
        int bufferWritePtr = 0;

        //////////////////////////////////////////////////////////////////////
        // Shannon energy processing variables
        //////////////////////////////////////////////////////////////////////
        // Stores the number of shannon energy values to be used in the moving
        // mean and standard deviation calculations
        int nn;
        // Vector for storing the un-normalised shannon energy values
        std::vector<float> shannonEnergy;
        int shanEngyPtr = 0;

        // Vector to store normalised shannon energy
        int normShanEngySize = 0;
        std::vector<float> normalisedShannonEnergy;
        std::vector<bool> shannonEnergyPeaks;
        int normShanEngyPtr = 0;
        float shanEngyMean = 0;
        float shanEngyMeanAccum = 0;
        float shanEngyStd = 0;

        // Test file output stream
        //std::ofstream energyFile;
        //std::ofstream peaksFile;


        //////////////////////////////////////////////////////////////////////
        // Peak finding variables
        //////////////////////////////////////////////////////////////////////
        // Threshold that the shannon energy must pass above to trigger peak
        // finding
        float threshold = 0.5;
        // Stores whether a peak has already been found for the current
        // threshold overshoot
        bool peakFound = false;
        float prevDiff = 1.0;
        
};

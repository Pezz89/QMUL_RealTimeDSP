#include <Bela.h>
#include <vector>
#include <cmath>
#include <algorithm>

#include "DebouncedButton.h"
#include "LED.h"

class BeatDetector {
    public:
        // Create default constructor
        BeatDetector() {}

        BeatDetector(BelaContext *context) : context(context){
            // Create a button object for each physical button
            debouncedButton1 = DebouncedButton(P8_07, 0.1*context->audioSampleRate);
            // Use the first second after the activation button is pressed for
            // system calibration
            calibrationTime = 1*context->audioSampleRate;
            // Declare hop and window size in samples for analysis
            hopSize = int(0.01*context->audioSampleRate);
            winSize = int(0.02*context->audioSampleRate);
            overlapFactor = ceil(winSize/hopSize);
            // Create circular buffer of size to store windows for processing
            // in real-time
            audioBuffer.resize(winSize+((overlapFactor-1)*hopSize), 0);

            // Create a write pointer at position 0 in buffer
            bufferWritePtr = 0;
            // Resize buffer to allow for a read pointer for each overlapping
            // window
            bufferReadPtrs.resize(overlapFactor);
            bufferReadPtrsPos.resize(overlapFactor);


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
            normShanEngySize = round((5.0*context->audioSampleRate) / hopSize);
            normalisedShannonEnergy.resize(normShanEngySize);
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
                    for(int i=ptr-(winSize-1); i != ptr; i++, i%=audioBuffer.size()){
                        powSamp = pow(audioBuffer[i], 2);
                        sigAccum += powSamp * log(powSamp);
                    }
                    rt_printf("sigAc:%f\n", sigAccum);
                    sigAccum = (-1/winSize)*sigAccum;
                    // If a Nan value is created due to a log(0), set the
                    // output to 0
                    if(isnan(sigAccum)) {
                        sigAccum = 0;
                    }
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
                rt_printf("%d\n", sampleCounter);
                sampleCounter += 1;
                return 0.0;
            }

            audioBuffer[bufferWritePtr] = sample;
            calculateAverageShannonEnergy();
            // Increment the buffer write pointer and wrap to keep within size
            // of buffer
            bufferWritePtr += 1;
            bufferWritePtr %= audioBuffer.size();

            transform(bufferReadPtrs.begin(), bufferReadPtrs.end(), bufferReadPtrs.begin(), bind2nd(std::plus<int>(), 1.0));
            transform(bufferReadPtrs.begin(), bufferReadPtrs.end(), bufferReadPtrs.begin(), bind2nd(std::modulus<int>(), audioBuffer.size()));
            transform(bufferReadPtrsPos.begin(), bufferReadPtrsPos.end(), bufferReadPtrsPos.begin(), bind2nd(std::plus<int>(), 1.0));
            transform(bufferReadPtrsPos.begin(), bufferReadPtrsPos.end(), bufferReadPtrsPos.begin(), bind2nd(std::modulus<int>(), audioBuffer.size()));
            // Increment the counter for the number of processed samples since
            // the button was pressed. This will eventually overflow but won't
            // be used by that point unless calibration is crazy long so it has
            // been left as is.
            sampleCounter += 1;
            return sample;
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
        // Vector of read pointers for buffer
        std::vector<int> bufferReadPtrs;
        std::vector<int> bufferReadPtrsPos;
        int bufferWritePtr = 0;

        //////////////////////////////////////////////////////////////////////
        // Shannon energy processing variables
        //////////////////////////////////////////////////////////////////////
        // Stores the number of shannon energy values to be used in the moving
        // mean and standard deviation calculations
        int n = 100;
        // Vector for storing the un-normalised shannon energy values
        std::vector<float> shannonEnergy;

        // Vector to store normalised shannon energy
        int normShanEngySize = 0;
        std::vector<float> normalisedShannonEnergy;
        int normShanEngyPtr = 0;
};

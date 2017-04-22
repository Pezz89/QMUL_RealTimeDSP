#include <Bela.h>
#include <vector>
#include <cmath>

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
            bufferReadPtrs.resize(overlapFactor)

            // Equivelants to Numpy style vector operations (np.arange and *=)
            // taken from: 
            // https://codereview.stackexchange.com/questions/77546/multiply-vector-elements-by-a-scalar-value-using-stl-and-templates
            // http://stackoverflow.com/questions/17694579/use-stdfill-to-populate-vector-with-increasing-numbers

            // Fill vector with values from 0-overlapFactor
            std::iota(std::begin(bufferReadPtrs), std::end(bufferReadPtrs), 0);
            // Multiply all values by -hopSize
            std::transform(bufferReadPtrs.begin(), bufferReadPtrs.end(), bufferReadPtrs.begin(), std::bind1st (std::multiplies<T>(), -hopSize));

            bufferReadPtrsPos = np.arange(self.overlapFactor, dtype=int)
            bufferReadPtrsPos *= -self.hopSize
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
        int bufferWritePtr = 0;
};

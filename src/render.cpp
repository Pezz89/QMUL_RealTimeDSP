
/*
 * assignment2_drums
 *
 * Second assignment for ECS732 RTDSP, to create a sequencer-based
 * drum machine which plays sampled drum sounds in loops.
 *
 * This code runs on BeagleBone Black with the Bela/BeagleRT environment.
 *
 * 2016, 2017
 * Becky Stewart
 *
 * 2015
 * Andrew McPherson and Victor Zappi
 * Queen Mary University of London
 */


#include <Bela.h>
#include <cmath>
#include "drums.h"
#include <array>
#include <vector>
#include <algorithm>
#include "filter.h"

/* Variables which are given to you: */

/* Drum samples are pre-loaded in these buffers. Length of each
 * buffer is given in gDrumSampleBufferLengths.
 */
extern float *gDrumSampleBuffers[NUMBER_OF_DRUMS];
extern int gDrumSampleBufferLengths[NUMBER_OF_DRUMS];

int gIsPlaying = 0;         /* Whether we should play or not. Implement this in Step 4b. */

/* Read pointer into the current drum sample buffer.
 *
 * TODO (step 3): you will replace this with two arrays, one
 * holding each read pointer, the other saying which buffer
 * each read pointer corresponds to.
 */
std::array<int, 16> gReadPointers;
std::array<int, 16> gDrumBufferForReadPointer;

/* Patterns indicate which drum(s) should play on which beat.
 * Each element of gPatterns is an array, whose length is given
 * by gPatternLengths.
 */
extern int *gPatterns[NUMBER_OF_PATTERNS];
extern int gPatternLengths[NUMBER_OF_PATTERNS];

/* These variables indicate which pattern we're playing, and
 * where within the pattern we currently are. Used in Step 4c.
 */
int gCurrentPattern = 0;
int gCurrentIndexInPattern = 0;

/* Triggers from buttons (step 2 etc.). Read these here and
 * do something if they are nonzero (resetting them when done). */
int gTriggerButton1 = 0;
int gPrevTriggerButton1 = 0;
int gTriggerButton2 = 0;
int gPrevTriggerButton2 = 0;

/* This variable holds the interval between events in **milliseconds**
 * To use it (Step 4a), you will need to work out how many samples
 * it corresponds to.
 */
float gEventIntervalMilliseconds = 0.05;

/* This variable indicates whether samples should be triggered or
 * not. It is used in Step 4b, and should be set in gpio.cpp.
 */
extern int gIsPlaying;

/* This indicates whether we should play the samples backwards.
 */
int gPlaysBackwards = 0;

/* For bonus step only: these variables help implement a fill
 * (temporary pattern) which is triggered by tapping the board.
 */
int gShouldPlayFill = 0;
bool gStartOfFill = false;
int gPreviousPattern = 0;

int gAudioFramesPerAnalogFrame;

int gSampleCounter = 0;

class Accelerometer {
    public:
        Accelerometer() {}
        Accelerometer(BelaContext *context, int pin1, int pin2, int pin3) : pinX(pin1), pinY(pin2), pinZ(pin3) {
            calibrationSamplesX.resize(context->audioSampleRate*0.5);
            calibrationSamplesY.resize(context->audioSampleRate*0.5);
            calibrationSamplesZ.resize(context->audioSampleRate*0.5);
            itX = calibrationSamplesX.begin();
            itY = calibrationSamplesY.begin();
            itZ = calibrationSamplesZ.begin();

            spikeFilter = Filter(900.0, context->audioSampleRate*0.5, true);
        }

        float readX(BelaContext *context, int n) {
            // On even audio samples: read analog input and return x value
            float val = analogRead(context, n/gAudioFramesPerAnalogFrame, pinX);
            if(!needsCalibrating)
                val -= averageX;
            return val;
        }

        float readY(BelaContext *context, int n) {
            // On even audio samples: read analog input and return y value
            float val = analogRead(context, n/gAudioFramesPerAnalogFrame, pinY);
            if(!needsCalibrating)
                val -= averageY;
            return val;
        }

        float readZ(BelaContext *context, int n) {
            float val = (5.0 / 3.3) * analogRead(context, n/gAudioFramesPerAnalogFrame, pinZ);
            // On even audio samples: read analog input and return z value
            if(!needsCalibrating)
                val = map(val, 1-averageZ, averageZ, -1.0, 1.0);
            return val;
        }

        void calibrate(BelaContext *context, int n) {
            if(needsCalibrating) {
                if(itX != calibrationSamplesX.end()) {
                    *itX = readX(context, n);
                    *itY = readY(context, n);
                    *itZ = readZ(context, n);
                    itX++;
                    itY++;
                    itZ++;
                }
                else
                {
                    averageX = accumulate(calibrationSamplesX.begin(), calibrationSamplesX.end(), 0.0)/calibrationSamplesX.size();
                    averageY = accumulate(calibrationSamplesY.begin(), calibrationSamplesY.end(), 0.0)/calibrationSamplesY.size();
                    averageZ = accumulate(calibrationSamplesZ.begin(), calibrationSamplesZ.end(), 0.0)/calibrationSamplesZ.size();
                    needsCalibrating = false;
                }
            }
        }

        enum orientation {
            upright = 1,
            left,
            right,
            front,
            back,
            upsidedown
        };

        int calculateOrientation(BelaContext *context, int n) {
            float x = readX(context, n);
            float y = readY(context, n);
            float z = readZ(context, n);

            int xOrient = hysterisisThreshold(x, -0.2, -0.1, hysts[0]) + hysterisisThreshold(x, 0.1, 0.2, hysts[1]);
            int yOrient = hysterisisThreshold(y, -0.2, -0.1, hysts[2]) + hysterisisThreshold(y, 0.1, 0.2, hysts[3]);
            int zOrient = hysterisisThreshold(z, -0.5, -0.3, hysts[4]) + hysterisisThreshold(z, 0.3, 0.5, hysts[5]);


            float filteredZ = spikeFilter.applyFilter(z);
            if(filteredZ > 0.3) {
                gShouldPlayFill = 1;
                gStartOfFill = true;
            }

            if(xOrient == 1 && yOrient == 1 && zOrient == 2) {
                prevOrient = upright;
                return upright;
            }
            if(xOrient == 0 && yOrient == 1 && zOrient == 1) {
                prevOrient = left;
                return left;
            }
            if(xOrient == 2 && yOrient == 1 && zOrient == 1) {
                prevOrient = right;
                return right;
            }
            if(xOrient == 1 && yOrient == 0 && zOrient == 1) {
                prevOrient = front;
                return front;
            }
            if(xOrient == 1 && yOrient == 2 && zOrient == 1) {
                prevOrient = back;
                return back;
            }
            if(xOrient == 1 && yOrient == 1 && zOrient == 0) {
                prevOrient = upsidedown;
                return upsidedown;
            }
            else {
                return prevOrient;
            }
        }

        bool hysterisisThreshold(float val, float low, float high, bool& hyst) {
            switch(hyst)
            {
                case false:
                    if(val < high)
                        return 0;
                    else
                        hyst = true;
                        return 1;
                    break;
                case true:
                    if(val > low)
                        return 1;
                    else
                        hyst = false;
                        return 0;
                    break;
            }
        }


        bool curentlyCalibrating() {
            return needsCalibrating;
        }

    private:
        int pinX, pinY, pinZ = 0;

        bool hysts[9] = {false};

        int prevOrient = upright;

        bool needsCalibrating = true;

        Filter spikeFilter;

        float averageX = 0;
        float averageY = 0;
        float averageZ = 0;
        std::vector<float> calibrationSamplesX;
        std::vector<float>::iterator itX;
        std::vector<float> calibrationSamplesY;
        std::vector<float>::iterator itY;
        std::vector<float> calibrationSamplesZ;
        std::vector<float>::iterator itZ;
};

Accelerometer gAccelerometer;

class LED {
    public:
        LED() {}
        LED(int pin, int time) : pin(pin), interval(time){}

        void trigger() {
            active = true;
            timer = 0;
        }

        void onIfActive(BelaContext *context, int n) {
            if(timer > interval) {
                active = false;
                timer = 0;
            }
            if(active) {
                digitalWrite(context, n, P8_07, 1);
                timer++;
            }
            else {
                digitalWrite(context, n, P8_07, 0);
            }
        }

    private:
        int pin;
        bool active = false;
        int timer = 0;
        int interval;
};

LED gLED;

class DebouncedButton {
    public:
        DebouncedButton() {}
        DebouncedButton(int pin, int debounceTime, bool defaultToggle=false) : counter(0), pin(pin), debounceTime(debounceTime), toggleBool(defaultToggle) {}

        // Reads value if enough time has past since the last positive
        // button reading
        bool getCurrentVal(BelaContext *context, int n) {
            // Inverted to account for the way that the buttons were wired to
            // produce 0s on press
            bool val = !digitalRead(context, n, pin);
            // If the button was previously off but is now on then start the
            // counter
            if(val == true && prevDBVal == false) {
                prevDBVal = val;
                counter++;
            }
            // If the button is now off but the debounce time has not elapsed
            // then it is assumed to be a bounce and remains on.
            else if(val == false && prevDBVal == true && counter < debounceTime) {
                val = true;
                counter++;
            }
            // If it is still on then just increment the counter.
            else if (val == true && prevDBVal == true) {
                counter++;
            }
            // if the debounce time has elapsed at it is no longer on then set
            // the button to the off state and reset the counter.
            else {
                prevDBVal = val;
                counter = 0;
            }
            return val;
        }

        // Use this function to convert buttons to toggles. These will return
        // 0 or 1 constatntly until the button is pressed again.
        bool getCurrentToggle(BelaContext *context, int n) {
            // Get the current debounced button value
            bool val = this->getCurrentVal(context, n);

            // if it has been pressed then toggle the toggle.
            if(prevTogVal == true && val == false) {
                toggleBool = !toggleBool;
            }
            prevTogVal = val;

            // return inverted toggle value. Somwhere the button value has
            // become inverted accidentally... Doesn't cause any problems and
            // buttons all work fine though...
            return !toggleBool;
        }
    private:
        int counter;
        int pin;
        int debounceTime;
        bool prevDBVal;
        bool prevTogVal;
        bool toggleBool;
};


// Create two global button variables, one for each physical button
DebouncedButton gDebouncedButton1;
DebouncedButton gDebouncedButton2;

// Boolean for storing whether bonus features are enabled or not.
bool gBonus = false;

bool setup(BelaContext *context, void *userData)
{
    // Set all read pointers to -1, signaling that they are currently
    // unused
    gDrumBufferForReadPointer.fill(-1);
    // Set all read pointers to the starting index of their buffers.
    gReadPointers.fill(0);

    // Code from one of the examples... Havent had time to work out whether
    // it's neccesary or not but doesn't do any harm.
    if(context->analogFrames == 0 || context->analogFrames > context->audioFrames) {
        rt_printf("Error: this example needs analog enabled, with 4 channels\n");
        return false;
    }

    // Useful calculations
    gAudioFramesPerAnalogFrame = context->audioFrames / context->analogFrames;

    // Activate digital pins for buttons and LED. Should have been put in their
    // classes really
    pinMode(context, 0, P8_07, OUTPUT);
    pinMode(context, 0, P8_08, INPUT);
    pinMode(context, 0, P8_09, INPUT);

    // Create an instance of the LED object for handeling LED related output
    gLED = LED(P8_07, 0.05*context->audioSampleRate);

    // Create instance of the buttons for reading input later
    gDebouncedButton1 = DebouncedButton(P8_08, 0.1*context->audioSampleRate);
    gDebouncedButton2 = DebouncedButton(P8_09, 0.1*context->audioSampleRate);

    // Create n acceerometer instance for reading in later
    gAccelerometer = Accelerometer(context, 7, 6, 5);

    return true;
}

// render() is called regularly at the highest priority by the audio engine.
// Input and output are given from the audio hardware and the other
// ADCs and DACs (if available). If only audio is available, numMatrixFrames
// will be 0.

void render(BelaContext *context, void *userData)
{
    // For every sample in the current block...
    for(unsigned int n=0; n<context->digitalFrames; n++){
        // For every sample of the analog inputs in the current block...
        if(!(n % gAudioFramesPerAnalogFrame)) {
            // Run calibration for the first half a second of the program
            gAccelerometer.calibrate(context, n);
            // When calibration has finished, begin processing analog inputs
            if(!gAccelerometer.curentlyCalibrating()) {
                // Caluclate the orientation of the drum machine as one of 6
                // orientations
                int orientation = gAccelerometer.calculateOrientation(context, n);

                // If bonus features are active...
                if(gBonus) {
                    // Map tempo to the Y axis of the accelerometer
                    float y = gAccelerometer.readY(context, n);
                    gEventIntervalMilliseconds = map(y, -0.2, 0.2, 0.05, 0.6);
                    // for all other orientations, change the drum pattern or
                    // play samples backwards as normal.
                    if(orientation < 4) {
                        gCurrentPattern = orientation-1;
                        gPlaysBackwards = false;
                    }
                    else if(orientation == 6) {
                        gPlaysBackwards = true;
                    }

                }
                else {
                    // Map the input of the potentiometer to the interval time
                    // in milliseconds. This allows for tempo adjustments using
                    // the pot.
                    gEventIntervalMilliseconds = map(analogRead(context, n/gAudioFramesPerAnalogFrame, 4), 0.0, 0.84, 0.05, 1.0);
                    // For each of the 5 orientations, set the current pattern
                    // to that orientation's index
                    if(orientation < 6) {
                        gCurrentPattern = orientation-1;
                        gPlaysBackwards = false;
                    }
                    // If the board is upside down (orientation no. 6) set the
                    // reverse samples flag
                    else {
                        gPlaysBackwards = true;
                    }
                }
            }
        }

        //Calculate the integer interval in samples
        int interval = round(gEventIntervalMilliseconds*context->audioSampleRate);
        // Set Bonus 
        if(gDebouncedButton1.getCurrentToggle(context, n))
            gBonus = true;
        else
            gBonus = false;
        // If the toggle is on then increment the counter (play samples)
        if(gDebouncedButton2.getCurrentToggle(context, n))
            gSampleCounter++;
        // When enough samples have elapsed, play the next beat.
        if(gSampleCounter >= interval) {
            gSampleCounter = 0;
            startNextEvent();
            // Also trigger an LED
            gLED.trigger();
        }
        // If the LED has been triggered, then light up for a specific amount
        // of time.
        gLED.onIfActive(context, n);


        float out = 0;
        // For each read pointer that isn't currently active
        for(int i=0; i<gDrumBufferForReadPointer.size(); i++) {
            int currentBuff = gDrumBufferForReadPointer[i];
            int currentPtr = gReadPointers[i];
            if(gDrumBufferForReadPointer[i] > -1) {
                if(gReadPointers[i] < gDrumSampleBufferLengths[currentBuff]) {
                    // If samples need to be played in reverse...
                    if(gPlaysBackwards){
                        // Read samples backwards from the end of the buffer.
                        out += gDrumSampleBuffers[currentBuff][gDrumSampleBufferLengths[currentBuff] - currentPtr];
                    } else {
                        // Read samples fowards from the current buffer using
                        // it's allocated pointer
                        out += gDrumSampleBuffers[currentBuff][currentPtr];
                    }
                    gReadPointers[i]++;
                }
                else
                {
                    // If the buffer has reached the end/start of the buffer,
                    // then reset pointer and mark it as free to use
                    gReadPointers[i] = 0;
                    gDrumBufferForReadPointer[i] = -1;
                }
            }
        }
        // write out audio
        audioWrite(context, n, 0, out);
        audioWrite(context, n, 1, out);
    }
}

/* Start playing a particular drum sound given by drumIndex.
 */
void startPlayingDrum(int drumIndex) {
    for(int i=0; i<gReadPointers.size(); i++) {
        if(gDrumBufferForReadPointer[i] == -1) {
            gDrumBufferForReadPointer[i] = drumIndex;
            gReadPointers[i] = 0;
            break;
        }
    }
}

/* Start playing the next event in the pattern */
void startNextEvent() {
    /* TODO in Step 4 */
    if(gShouldPlayFill) {
        gPreviousPattern = gCurrentPattern;
        gCurrentPattern = FILL_PATTERN;
    }
    if(gStartOfFill) {
        gCurrentIndexInPattern = 0;
        gStartOfFill = false;
    }
    const int currentPatterLength = gPatternLengths[gCurrentPattern];
    for(int i = 0; i < sizeof(gDrumSampleBuffers) * sizeof(float*); i++){
        if(eventContainsDrum(gPatterns[gCurrentPattern][gCurrentIndexInPattern], i))
            startPlayingDrum(i);
    }
    if(gCurrentIndexInPattern == currentPatterLength-1 && gCurrentPattern == FILL_PATTERN) {
        gCurrentPattern = gPreviousPattern;
        gShouldPlayFill = false;
    }
    gCurrentIndexInPattern = (gCurrentIndexInPattern+1+currentPatterLength)%currentPatterLength;
}

/* Returns whether the given event contains the given drum sound */
int eventContainsDrum(int event, int drum) {
    if(event & (1 << drum))
        return 1;
    return 0;
}

// cleanup_render() is called once at the end, after the audio has stopped.
// Release any resources that were allocated in initialise_render().

void cleanup(BelaContext *context, void *userData)
{

}

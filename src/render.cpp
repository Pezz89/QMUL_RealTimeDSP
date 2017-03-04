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
int gEventIntervalMilliseconds = 250;

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
int gPreviousPattern = 0;

int gAudioFramesPerAnalogFrame = 0;

int gSampleCounter = 0;

class DebouncedButton {
    public:
        DebouncedButton() {}
        DebouncedButton(int pin, int debounceTime) : counter(0), pin(pin), debounceTime(debounceTime) {}

        bool getCurrentVal(BelaContext *context, int n) {
            bool val = !digitalRead(context, n, pin);
            if(val == true && prevDBVal == false) {
                prevDBVal = val;
                counter++;
            }
            else if(val == false && prevDBVal == true && counter < debounceTime) {
                val = true;
                counter++;
            }
            else if (val == true && prevDBVal == true) {
                counter++;
            }
            else {
                prevDBVal = val;
                counter = 0;
            }
            return val;
        }

        bool getCurrentToggle(BelaContext *context, int n) {
            bool val = this->getCurrentVal(context, n);

            if(prevTogVal == false && val == true) {
                toggleBool = !toggleBool;
            }
            prevTogVal = val;

            return toggleBool;
        }
    private:
        int counter;
        int pin;
        int debounceTime;
        bool prevDBVal;
        bool prevTogVal;
        bool toggleBool;

        //bool toggleVal;
};


DebouncedButton gDebouncedButton1;
DebouncedButton gDebouncedButton2;

/* TODO: Declare any further global variables you need here */

// setup() is called once before the audio rendering starts.
// Use it to perform any initialisation and allocation which is dependent
// on the period size or sample rate.
//
// userData holds an opaque pointer to a data structure that was passed
// in from the call to initAudio().
//
// Return true on success; returning false halts the program.

bool setup(BelaContext *context, void *userData)
{
    gDrumBufferForReadPointer.fill(-1);
    gReadPointers.fill(0);

    if(context->analogFrames == 0 || context->analogFrames > context->audioFrames) {
        rt_printf("Error: this example needs analog enabled, with 4 channels\n");
        return false;
    }

    // Useful calculations
    gAudioFramesPerAnalogFrame = context->audioFrames / context->analogFrames;

    pinMode(context, 0, P8_07, OUTPUT);
    pinMode(context, 0, P8_08, INPUT);
    pinMode(context, 0, P8_09, INPUT);

    gDebouncedButton1 = DebouncedButton(P8_08, 1.02*context->audioSampleRate);
    gDebouncedButton2 = DebouncedButton(P8_09, 1.02*context->audioSampleRate);
    return true;
}

// render() is called regularly at the highest priority by the audio engine.
// Input and output are given from the audio hardware and the other
// ADCs and DACs (if available). If only audio is available, numMatrixFrames
// will be 0.

void render(BelaContext *context, void *userData)
{
    for(unsigned int n=0; n<context->digitalFrames; n++){
        bool button1 = gDebouncedButton1.getCurrentVal(context, n);
        bool button2 = gDebouncedButton2.getCurrentVal(context, n);
        //False value means the button is pressed due to the use of a pull up
        //resistor
        if(button1 == true){
            gTriggerButton1 = 1;
        }
        else
        {
            gTriggerButton1 = 0;
        }
        if(gTriggerButton1 && gTriggerButton1 != gPrevTriggerButton1) {
            startPlayingDrum(3);
            startPlayingDrum(0);
        }
        gPrevTriggerButton1 = gTriggerButton1;

        if(button2 == true){
            gTriggerButton2 = 1;
        }
        else
        {
            gTriggerButton2 = 0;
        }
        if(gTriggerButton2 && gTriggerButton2 != gPrevTriggerButton2) {
            startPlayingDrum(4);
        }
        gPrevTriggerButton2 = gTriggerButton2;

        int interval = round(1*context->audioSampleRate);
        // Increment counter every sample

        /*
        gSampleCounter++;
        if(gSampleCounter == interval) {
            gSampleCounter = 0;
            startNextEvent();
        }
        */

        float out = 0;
        for(int i=0; i<gDrumBufferForReadPointer.size(); i++) {
            int currentBuff = gDrumBufferForReadPointer[i];
            int currentPtr = gReadPointers[i];
            if(gDrumBufferForReadPointer[i] > -1) {
                if(gReadPointers[i] < gDrumSampleBufferLengths[currentBuff]) {
                    out += gDrumSampleBuffers[currentBuff][currentPtr];
                    gReadPointers[i]++;
                }
                else 
                {
                    gReadPointers[i] = 0;
                    gDrumBufferForReadPointer[i] = -1;
                }
            }
        }
        audioWrite(context, n, 0, out);
        audioWrite(context, n, 1, out);
    }
    /* TODO: your audio processing code goes here! */

    /* Step 2: use gReadPointer to play a drum sound */

    /* Step 3: use multiple read pointers to play multiple drums */

    /* Step 4: count samples and decide when to trigger the next event */
}

/* Start playing a particular drum sound given by drumIndex.
 */
void startPlayingDrum(int drumIndex) {
    /* TODO in Steps 3a and 3b */

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
    startPlayingDrum(0);
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

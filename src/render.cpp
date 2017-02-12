/*
 * assignment1_crossover
 * RTDSP 2017
 *
 * First assignment for ECS732 RTDSP, to implement a 2-way audio crossover
 * using the BeagleBone Black.
 *
 * Andrew McPherson and Victor Zappi
 * Modified by Becky Stewart
 * Queen Mary, University of London
 */

#include <Bela.h>
#include <cmath>
#include <Utilities.h>
#include <memory>
#include <Scope.h>

#include "filter.h"
#include "userOptions.h"

/* TASK: declare any global variables you need here */

// setup() is called once before the audio rendering starts.
// Use it to perform any initialisation and allocation which is dependent
// on the period size or sample rate.
//
// userData holds an opaque pointer to a data structure that was passed
// in from the call to initAudio().
//
// Return true on success; returning false halts the program.

// Create filter objects using unique pointers for automatic memory
// deallocation
std::unique_ptr<Filter> gLowPass;
std::unique_ptr<Filter> gHighPass;

// Instantiate the scope
Scope scope;
bool setup(BelaContext *context, void *userData)
{
    // Set the Bela IDE scope to read two outputs
    scope.setup(2, context->audioSampleRate);
    // Set default values for when command line arguments are not provided by
    // the user
    float crossoverFrequency = 1000.0;
    bool linkwitzRiley = false;
    // Retrieve a parameter passed in from the initAudio() call
    if(userData != 0)
        crossoverFrequency = (*(UserOpts *)userData).frequency;
        linkwitzRiley = (*(UserOpts *)userData).linkwitzRiley;

    // Create a low-pass and high-pass filter object using the crossover
    // frequency and filter design specified by CLI arguments
    gLowPass.reset(new Filter(crossoverFrequency, context->audioSampleRate, false, linkwitzRiley));
    gHighPass.reset(new Filter(crossoverFrequency, context->audioSampleRate, true, linkwitzRiley));

    return true;
}

// render() is called regularly at the highest priority by the audio engine.
// Input and output are given from the audio hardware and the other
// ADCs and DACs (if available). If only audio is available, numMatrixFrames
// will be 0.

void render(BelaContext *context, void *userData)
{
    /* TASK:
     * Mix the two input channels together.
     *
     * Apply a lowpass filter and a highpass filter, sending the
     * lowpass output to the left channel and the highpass output
     * to the right channel.
     */
    for(unsigned int n = 0; n < context->audioFrames; n++) {
        // Read audio inputs
        float leftIn = audioRead(context,n,0);
        float rightIn = audioRead(context,n,1);

        // Convert input to mono
        float monoSamp = (leftIn + rightIn) * 0.5;

        // Apply filter to current sample for both channels.
        // filtering/buffering is handeled within the Filter objects.
        float leftOut = gLowPass->applyFilter(monoSamp);
        float rightOut = gHighPass->applyFilter(monoSamp);
        
        // Plot the output of the filters to the IDE scope
        scope.log(leftOut, rightOut);
        // Write the sample into the output buffer
        audioWrite(context, n, 0, leftOut);
        audioWrite(context, n, 1, rightOut);
    }
}

// cleanup_render() is called once at the end, after the audio has stopped.
// Release any resources that were allocated in initialise_render().

void cleanup(BelaContext *context, void *userData)
{
    // Cleanup wasn't necessary through the use of unique pointers.
}

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

#include "BeatDetector.h"

// Create a beat detector object for processing PCG data
std::unique_ptr<BeatDetector> gBeatDetector;

bool setup(BelaContext *context, void *userData)
{
    gBeatDetector = std::make_unique<BeatDetector>(context);
    return true;
}

void render(BelaContext *context, void *userData)
{
    for(unsigned int n = 0; n < context->audioFrames; n++) {
        // Read audio inputs
        float leftIn = audioRead(context,n,0);
        float rightIn = audioRead(context,n,1);

        // Convert input to mono
        float monoSamp = (leftIn + rightIn) * 0.5;

        float out = gBeatDetector->processSample(monoSamp, n);

        // Write the sample into the output buffer
        audioWrite(context, n, 0, out);
        audioWrite(context, n, 1, out);
    }
}

void cleanup(BelaContext *context, void *userData)
{
    // Cleanup wasn't necessary through the use of unique pointers.
}

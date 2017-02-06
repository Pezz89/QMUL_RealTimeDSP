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

/* TASK: declare any global variables you need here */

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
	float crossoverFrequency;
	// Retrieve a parameter passed in from the initAudio() call
	if(userData != 0)
		crossoverFrequency = *(float *)userData;

	/* TASK:
	 * Calculate the filter coefficients based on the given
	 * crossover frequency.
	 *
	 * Initialise any previous state (clearing buffers etc.)
	 * to prepare for calls to render()
	 */

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
}

// cleanup_render() is called once at the end, after the audio has stopped.
// Release any resources that were allocated in initialise_render().

void cleanup(BelaContext *context, void *userData)
{
	/* TASK:
	 * If you allocate any memory, be sure to release it here.
	 * You may or may not need anything in this function, depending
	 * on your implementation.
	 */
}

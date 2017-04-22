/*
 * assignment1_crossover
 * RTDSP 2017
 *
 * First assignment for ECS732 RTDSP, to implement a 2-way audio crossover
 * using the BeagleBone Black.
 *
 * Andrew McPherson and Victor Zappi
 * Edited by Becky Stewart
 * Queen Mary, University of London
 */

#include <unistd.h>
#include <iostream>
#include <cstdlib>
#include <libgen.h>
#include <signal.h>
#include <getopt.h>
#include <Bela.h>

#include "userOptions.h"

using namespace std;

// Handle Ctrl-C by requesting that the audio rendering stop
void interrupt_handler(int var)
{
    gShouldStop = true;
}

// Print usage information
void usage(const char * processName)
{
    cerr << "Usage: " << processName << " [options]" << endl;

    Bela_usage();

    cerr << "   --help [-h]:                Print this menu\n";
}

int main(int argc, char *argv[])
{
    BelaInitSettings settings;  // Standard audio settings
    
    // Initialize default values for command line options
    UserOpts uOpts = {false};

    struct option customOptions[] =
    {
        {"help", 0, NULL, 'h'},
        {"frequency", 1, NULL, 'f'},
        {"linkwitzriley", 0, NULL, 'l'},
        {NULL, 0, NULL, 0}
    };

    // Set default settings
    Bela_defaultSettings(&settings);

    // Parse command-line arguments
    while (1) {
        int c;
        if ((c = Bela_getopt_long(argc, argv, "hf:", customOptions, &settings)) < 0)
                break;
        switch (c) {
        case 'h':
                usage(basename(argv[0]));
                exit(0);
        case 'f':
                break;
        case 'l':
                break;
        case '?':
        default:
                usage(basename(argv[0]));
                exit(1);
        }
    }

    // Initialise the PRU audio device
    if(Bela_initAudio(&settings, &uOpts) != 0) {
        cout << "Error: unable to initialise audio" << endl;
        return -1;
    }

    // Start the audio device running
    if(Bela_startAudio()) {
        cout << "Error: unable to start real-time audio" << endl;
        return -1;
    }

    // Set up interrupt handler to catch Control-C
    signal(SIGINT, interrupt_handler);
    signal(SIGTERM, interrupt_handler);

    // Run until told to stop
    while(!gShouldStop) {
        usleep(100000);
    }

    // Stop the audio device
    Bela_stopAudio();

    // Clean up any resources allocated for audio
    Bela_cleanupAudio();

    // All done!
    return 0;
}

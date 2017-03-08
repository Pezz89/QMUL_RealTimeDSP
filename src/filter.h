#include <Bela.h>
#include <cmath>
#include <array>
#include <vector>

class Filter {
    public:
        Filter() {}
        Filter(const float &crossoverFrequency, const float &fs, const bool &highpass=false, const bool &linkwitzRiley=false) {
            // Filter class constructor is used for the calculation of filter
            // coefficients and delay line memory allocation.

            // Calculate ratio between cutoff frequency and sampling rate
            const double wc = crossoverFrequency/fs;

            // Deifine Q as the square root of 2
            const double q = sqrt(2.0);

            // Warp the frequency to convert from continuous to discrete time cutoff
            const double wd1 = 1.0 / tan(M_PI*wc);

            // Calculate coefficients from equation and store in a vector
            numerator.push_back(1.0 / (1.0 + q*wd1 + pow(wd1, 2)));
            numerator.push_back(2 * numerator[0]);
            numerator.push_back(numerator[0]);
            denominator.push_back(1.0);
            denominator.push_back(2.0 * (pow(wd1, 2) - 1.0) * numerator[0]);
            denominator.push_back(-(1.0 - q * wd1 + pow(wd1, 2)) * numerator[0]);
            // If the filter is a high pass filter, convert numerator
            // coefficients to reflect this
            if(highpass) {
                numerator[0] = numerator[0] * pow(wd1, 2);
                numerator[1] = -numerator[1] * pow(wd1, 2);
                numerator[2] = numerator[2] * pow(wd1, 2);
            }
            // If the filter is using the Linkwitz-Riley filter structure,
            // convolve the numerator and denominator generated for the 2nd
            // order butterworth filter with themselves. This creates the 5
            // coefficients of 2 cascaded 2nd order butterworth filters needed
            // for this filter structure.
            if(linkwitzRiley) {
                numerator = convolve(numerator, numerator);
                denominator = convolve(denominator, denominator);
                //
                // Print coefficients to the console
                rt_printf("\nNumerator:\t\t%f %f %f %f %f\nDenominator:\t\t%f %f %f %f %f\nCrossover Frequency:\t%f\nHighpass:\t\t%s\nFilter Type:\t\t%s\n\n",
                        numerator[0],
                        numerator[1],
                        numerator[2],
                        numerator[3],
                        numerator[4],
                        denominator[0],
                        denominator[1],
                        denominator[2],
                        denominator[3],
                        denominator[4],
                        crossoverFrequency,
                        highpass ? "true" : "false",
                        linkwitzRiley ? "4th-Order Linkwitz-Riley" : "2nd-Order Butterworth");
            } else {
                // Print coefficients to the console
                rt_printf("\nNumerator:\t\t%f %f %f\nDenominator:\t\t%f %f %f\nCrossover Frequency:\t%f\nHighpass:\t\t%s\nFilter Type:\t\t%s\n\n",
                        numerator[0],
                        numerator[1],
                        numerator[2],
                        denominator[0],
                        denominator[1],
                        denominator[2],
                        crossoverFrequency,
                        highpass ? "true" : "false",
                        linkwitzRiley ? "4th-Order Linkwitz-Riley" : "2nd-Order Butterworth");
            }
            // Allocate memory for delay line based on the number of
            // coefficients generated. Initialize vectors with values of 0.
            inputDelayBuf.assign(int(numerator.size()), 0.0);
            outputDelayBuf.assign(int(denominator.size()), 0.0);
            // Store the delay size of delay buffers
            inputDelaySize = inputDelayBuf.size();
            outputDelaySize = outputDelayBuf.size();
        }

        double applyFilter(const float &x0) {
            // Increment the write pointer of the delay buffer storing input
            // samples
            ++inputDelayBufWritePtr;
            // Wrap values to withink size of buffer. Prevents an integer
            // overflow
            inputDelayBufWritePtr = (inputDelayBufWritePtr+inputDelaySize)%inputDelaySize;

            // Increment the write pointer of the delay buffer storing output
            // samples
            ++outputDelayBufWritePtr;
            // Wrap values to withink size of buffer. Prevents an integer
            // overflow
            outputDelayBufWritePtr = (outputDelayBufWritePtr+outputDelaySize)%outputDelaySize;

            // Set the current value of the input delay buffer to the value of
            // the sample provided to the function
            inputDelayBuf[(inputDelayBufWritePtr+inputDelaySize)%inputDelaySize] = x0;

            // Initialize a variable to store an output value
            double y = 0;
            // Accumulate each sample in the input delay buffer, multiplied by
            // it's corresponding coefficient
            for(unsigned int i = 0; i < numerator.size(); i++) {
                y += inputDelayBuf[(inputDelayBufWritePtr-i+inputDelaySize)%inputDelaySize] * numerator[i];
            }
            // decumulate each sample in the output delay buffer (aside from
            // the current index), multiplied by it's corresponding coefficient
            for(unsigned int i = 1; i < denominator.size(); i++) {
                y += outputDelayBuf[(outputDelayBufWritePtr-i+outputDelaySize)%outputDelaySize] * denominator[i];
            }
            // Scale by first coefficient in the denominator (always 1 in
            // current implementation, so added only for generalization of the
            // method for future use)
            //y /= denominator[0];

            // Store the calculated output sample in the output sample delay
            // buffer
            outputDelayBuf[(outputDelayBufWritePtr+outputDelaySize)%outputDelaySize] = y;

            return y;
        }

    private:
        // Vector used for dynamic allocation of numerators and denominators
        // based on filter order
        // doubles used for maximum precision. Probably not necessary...
        std::vector<double> numerator;
        std::vector<double> denominator;
        // Vector used for dynamic allocation of delay line size based on
        // filter type
        std::vector<float> inputDelayBuf;
        std::vector<float> outputDelayBuf;
        std::vector<float>::size_type inputDelaySize, outputDelaySize;
        unsigned int inputDelayBufWritePtr, outputDelayBufWritePtr = 0;

        // Convolution function adapted from: http://stackoverflow.com/questions/24518989/how-to-perform-1-dimensional-valid-convolution
        template<typename T>
        std::vector<T> convolve(std::vector<T> const &f, std::vector<T> const &g) {
            // Calculate the size of input vectors
            const int nf = f.size();
            const int ng = g.size();
            // Calculate the size of output vector as the combined size of both
            // input vectors, minus 1
            const int n  = nf + ng - 1;
            // Initialize vector of the same input type as input vectors
            // Allocate memory for all elements of the output to be calculated
            std::vector<T> out(n, T());
            // For each output element...
            for(auto i(0); i < n; ++i) {
                // Calculate minimum and maximum indexes to iterate over each
                // vector
                const int jmn = (i >= ng - 1)? i - (ng - 1) : 0;
                const int jmx = (i <  nf - 1)? i : nf - 1;
                // Accumulate the multiplication of elements in both vectors,
                // based on the indexes calculated, to give the output value at
                // the current output index
                for(auto j(jmn); j <= jmx; ++j) {
                    out[i] += (f[j] * g[i - j]);
                }
            }
            return out;
        }
};

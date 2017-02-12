#include <Bela.h>
#include <cmath>
#include <array>
#include <vector>

class Filter {
    public:
        Filter(float crossoverFrequency, float fs, bool highpass=false, bool linkwitzRiley=false) {

            // Calculate ratio between cutoff frequency and sampling rate
            double wc = crossoverFrequency/fs;

            // Deifine Q as the square root of 2
            double q = sqrt(2.0);

            // Warp the frequency to convert from continuous to discrete time cutoff
            double wd1 = 1.0 / tan(M_PI*wc);
 
            // Calculate coefficients from equation
            numerator.push_back(1.0 / (1.0 + q*wd1 + pow(wd1, 2)));
            numerator.push_back(2 * numerator[0]);
            numerator.push_back(numerator[0]);
            denominator.push_back(1.0);
            denominator.push_back(-2.0 * (pow(wd1, 2) - 1.0) * numerator[0]);
            denominator.push_back((1.0 - q * wd1 + pow(wd1, 2)) * numerator[0]);
            if(highpass) {
                numerator[0] = numerator[0] * pow(wd1, 2);
                numerator[1] = -numerator[1] * pow(wd1, 2);
                numerator[2] = numerator[2] * pow(wd1, 2);
            }
            if(linkwitzRiley) {
                //rt_printf("Num size: %d\n", numerator.size());
                //rt_printf("Den size: %d\n", denominator.size());
                numerator = convolve(numerator, numerator);
                denominator = convolve(denominator, denominator);
                rt_printf("Num size: %d\n", numerator.size());
                rt_printf("Den size: %d\n", denominator.size());
                rt_printf("%f %f %f %f %f\n", numerator[0], numerator[1], numerator[2], numerator[3], numerator[4]);
                rt_printf("%f %f %f %f %f\n", denominator[0], denominator[1], denominator[2], denominator[3], denominator[4]);
            }
            else {
            rt_printf("Numerator: %f %f %f\nDenominator: %f %f %f\nHighpass: %s\n", 
                    numerator[0], 
                    numerator[1], 
                    numerator[2], 
                    denominator[0], 
                    denominator[1],
                    denominator[2],
                    highpass ? "true":"false");
            }
            inputDelayBuf.assign(int(numerator.size()), 0.0);
            outputDelayBuf.assign(int(denominator.size()), 0.0);
            // Initialize delay buffer to be two samples long for the
            // second-order Butterworth filter
            inputDelaySize = inputDelayBuf.size();
            outputDelaySize = outputDelayBuf.size();
        }
        
        float applyFilter(float x0) {
            ++inputDelayBufWritePtr;
            inputDelayBufWritePtr = (inputDelayBufWritePtr+inputDelaySize)%inputDelaySize;

            ++outputDelayBufWritePtr;
            outputDelayBufWritePtr = (outputDelayBufWritePtr+outputDelaySize)%outputDelaySize;

            inputDelayBuf[(inputDelayBufWritePtr+inputDelaySize)%inputDelaySize] = x0;

            float y = 0;
            for(unsigned int i = 0; i < inputDelaySize; i++) {
                y += inputDelayBuf[(inputDelayBufWritePtr-i+inputDelaySize)%inputDelaySize] * numerator[i];
            }
            for(unsigned int i = 1; i < outputDelaySize; i++) {
                y -= outputDelayBuf[(outputDelayBufWritePtr-i+outputDelaySize)%outputDelaySize] * denominator[i];
            }
            y /= denominator[0];

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
          int const nf = f.size();
          int const ng = g.size();
          int const n  = nf + ng - 1;
          std::vector<T> out(n, T());
          for(auto i(0); i < n; ++i) {
            int const jmn = (i >= ng - 1)? i - (ng - 1) : 0;
            int const jmx = (i <  nf - 1)? i            : nf - 1;
            for(auto j(jmn); j <= jmx; ++j) {
              out[i] += (f[j] * g[i - j]);
            }
          }
          return out; 
        }
};

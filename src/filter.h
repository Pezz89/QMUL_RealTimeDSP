#include <Bela.h>
#include <cmath>
#include <array>
#include <vector>

class Filter {
    public:
        Filter(float crossoverFrequency, float fs, bool highpass=false, bool linkwitzRiley=false) : 
            inputDelayBuf(3, 0), outputDelayBuf(3, 0) {
            // Initialize delay buffer to be two samples long for the
            // second-order Butterworth filter
            inputDelaySize = inputDelayBuf.size();
            outputDelaySize = outputDelayBuf.size();

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
            //denominator.push_back(1.0);
            denominator.push_back(-2.0 * (pow(wd1, 2) - 1.0) * numerator[0]);
            denominator.push_back((1.0 - q * wd1 + pow(wd1, 2)) * numerator[0]);
            if(highpass) {
                numerator[0] = numerator[0] * pow(wd1, 2);
                numerator[1] = -numerator[1] * pow(wd1, 2);
                numerator[2] = numerator[2] * pow(wd1, 2);
            }
            if(linkwitzRiley) {
            }
            rt_printf("Numerator: %f %f %f\nDenominator: %f %f\nHighpass: %s\n", 
                    numerator[0], 
                    numerator[1], 
                    numerator[2], 
                    denominator[0], 
                    denominator[1],
                    highpass ? "true":"false");
            
        }
        
        float applyFilter(float x0) {
            float y = 0;

                

            ++inputDelayBufWritePtr;
            inputDelayBufWritePtr = (inputDelayBufWritePtr+inputDelaySize)%inputDelaySize;

            ++outputDelayBufWritePtr;
            outputDelayBufWritePtr = (outputDelayBufWritePtr+outputDelaySize)%outputDelaySize;

            inputDelayBuf[(inputDelayBufWritePtr+inputDelaySize)%inputDelaySize] = x0;

            y = x0 * numerator[0] 
                + inputDelayBuf[(inputDelayBufWritePtr-1+inputDelaySize)%inputDelaySize] * numerator[1] 
                + inputDelayBuf[(inputDelayBufWritePtr-2+inputDelaySize)%inputDelaySize] * numerator[2]
                - outputDelayBuf[(outputDelayBufWritePtr-1+outputDelaySize)%outputDelaySize] * denominator[0]
                - outputDelayBuf[(outputDelayBufWritePtr-2+outputDelaySize)%outputDelaySize] * denominator[1];

            //rt_printf("%d\n", y);

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

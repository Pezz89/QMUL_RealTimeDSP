#include <Bela.h>

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

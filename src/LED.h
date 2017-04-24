
class LED {
    public:
        LED() {}
        LED(int pin, int time) : pin(pin), interval(time){}

        void trigger() {
            // Trigger the LED to be illuminated and reset the timer
            active = true;
            timer = 0;
        }

        void onIfActive(BelaContext *context, int n) {
            // If the LED has been illuminated for the required duration, turn
            // it off and reset the timer
            if(timer > interval) {
                active = false;
                timer = 0;
            }
            // If the LED is on then write this out to the digital pin and
            // increment the times
            if(active) {
                digitalWrite(context, n, pin, 1);
                timer++;
            }
            // Otherwise, turn the LED off
            else {
                digitalWrite(context, n, pin, 0);
            }
        }

    private:
        int pin;
        bool active = false;
        int timer = 0;
        int interval;
};

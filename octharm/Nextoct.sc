Nextoct : UGen {
    *kr { arg freq = 440.0, candidate = 60.0;
        ^this.multiNew('control', freq, candidate)
    }
}
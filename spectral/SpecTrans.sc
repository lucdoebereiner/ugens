PV_MagLagUD : PV_ChainUGen {
	*new { arg buffer, lagU = 0.0, lagD = 0.0;
		^this.multiNew('control', buffer, lagU, lagD)
	}
}


PV_PollMags : PV_ChainUGen {
	*new { arg buffer, poll = 0.0;
		^this.multiNew('control', buffer, poll)
	}
}


PV_MagMinusOct : UGen
{
	*new { arg bufferA, bufferB;
		^this.multiNew('control', bufferA, bufferB)
	}
}

PV_FreqDiffs : UGen
{
	*new { arg bufferA, bufferB, bufferC;
		^this.multiNew('control', bufferA, bufferB, bufferC)
	}
}



PV_MagPeaksDecay : PV_ChainUGen {
	*new { arg buffer, minPeak = 0 , maxPeak = 0.5, minDecay = 0, maxDecay = 0.3, upTime = 0, sustain = 20, maxNSustain = 5, minDistFac = 1.2;
		^this.multiNew('control', buffer, minPeak, maxPeak, minDecay, maxDecay, upTime, sustain, maxNSustain, minDistFac)
	}
}


MagPeaksFreqs : MultiOutUGen {

	*kr { arg buffer, n, minPeak = 0.001;
		^this.multiNew('control', buffer, n, minPeak)
	}


	init { arg ... theInputs;
		inputs = theInputs;
		^this.initOutputs(inputs[1]*2, 'control');
	}
	
}


WrapOct : UGen {
    *kr { arg freq = 440.0, minFreq = 20.0, maxFreq = 12000.0;
        ^this.multiNew('control', freq, minFreq, maxFreq)
    }
}




OnOffChange {

	*kr { arg in, level = 1, onTimes = [1], offTime = 4;

		var changed = Changed.kr(in, 1);

		var thisOnTime = TChoose.kr(changed, onTimes);

		var trigAmp = Trig1.kr(changed, thisOnTime+offTime);
		
		var amp = EnvGen.kr(Env([0.000001,1,1,0.000001],[0.4,0.1,0.5],[3,1,-3]),
			trigAmp,
			levelScale: Latch.kr(level,trigAmp),
			timeScale: Latch.kr(thisOnTime,trigAmp));

		var freq = Gate.kr(in, amp < 0.000001);

		^[freq, amp];
	}
}




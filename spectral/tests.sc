s.boot;
s.reboot;
s.quit;

s.options.numOutputBusChannels = 6;

// test output
{ Out.ar(1, SinOsc.ar(800) * 0.01) }.play

b = Buffer.read(s, "/home/luc/Work/plastik/bass.wav")

~c = Buffer.read(s, "/home/luc/Work/plastik/cello_transitions.wav")

~c = Buffer.read(s, "/home/luc/Work/plastik/multi.wav")

~c = Buffer.read(s, "/home/luc/Work/plastik/my1.wav")

~c.play

(
{
	var maxed, lagged, subtr, lagShort, input = PlayBuf.ar(1,~c.bufnum,loop:1);
	//	var input = SoundIn.ar(0);
	var chain = FFT(LocalBuf(2048), input);
	var whitened = PV_Whiten(PV_Copy(chain, LocalBuf(2048)), LocalBuf(2048));
	var contrast = PV_MagMul(chain,whitened);
	//	chain = PV_Diffuser(contrast, 1);
	chain = contrast;
	
	lagged = PV_MagLagUD(PV_Copy(chain, LocalBuf(2048)), 0.2, 2.5);

	//	maxed = PV_LocalMax(PV_Copy(chain, LocalBuf(2048)), 0.5);
	maxed = PV_MagAbove(PV_Copy(chain, LocalBuf(2048)), 0.25);

	// new local maxima
	subtr = PV_MagMinus(maxed, lagged);

	lagShort = PV_MagLagUD(PV_Copy(subtr, LocalBuf(2048)), 0, 0.8);
	

	[IFFT.ar(lagShort), input] * 0.5;
}.play
)


(
{
	var subtr, lagged, input = SinOsc.ar(MouseX.kr(100,5000,1));
	var chain = FFT(LocalBuf(2048), input, 0.5);
	lagged = PV_MagLagUD(PV_Copy(chain, LocalBuf(2048)), 0.25, 1.5);
	subtr = PV_MagMinus(chain, lagged);
	IFFT.ar(chain).dup * 0.1;
}.play
)

(
{
	var contrast, whitened, past, new, peaks,rndPhase,lagged,
	input = PlayBuf.ar(1,~c.bufnum,loop:1);
	//	input = Mix.ar(SinOsc.ar([300,4000])) * 0.3;
	var chain = FFT(LocalBuf(4096), input, 0.5);
	//var rnd = FFT(LocalBuf(4096), WhiteNoise.ar(0.6), 0.5);
	whitened = PV_Whiten(PV_Copy(chain, LocalBuf(4096)), LocalBuf(4096));
	contrast = PV_MagMul(chain,whitened);
	past = PV_MagLagUD(PV_Copy(contrast, LocalBuf(4096)), 0.5, 0.3);
	new = PV_MagMinus(PV_Copy(contrast, LocalBuf(4096)), past);
	peaks = PV_MagPeaksDecay(new, 0.1, 1, 0.25, 3, 0, 100, 4, 2.2);
	//	lagged = PV_MagLagUD(peaks, 0, 1);
	[IFFT.ar(new), input * 0.1] ;
}.play
)


// todo
// peakedness merken und als stack dauer verwenden
// magnitude scaling / amplitude of it

//44100 / 4096

x = { SinOsc.ar(350) * 0.1 }.play
x.free

(
{
	var fftA, fftB, whitenedA, whitenedB, inputA, inputB, contrastA, contrastB,
	longTermA, longTermB, newA, newB, possible4A, possible4B;
	//	inputA = PlayBuf.ar(1,~c.bufnum,loop:1);
	inputA = PMOsc.ar(200, 120, MouseX.kr(0, 10, 0)) * 0.5;
	inputB = PMOsc.ar(1220, 190, MouseX.kr(0, 10, 0)) * 0.5;
	//DelayL.ar(inputA, 10, MouseY.kr(0,10).poll);
	fftA = FFT(LocalBuf(2048), inputA);
	fftB = FFT(LocalBuf(2048), inputB);
	whitenedA = PV_Whiten(PV_Copy(fftA, LocalBuf(2048)), LocalBuf(2048));
	whitenedB = PV_Whiten(PV_Copy(fftB, LocalBuf(2048)), LocalBuf(2048));
	contrastA = PV_MagMul(fftA,whitenedA);
	contrastB = PV_MagMul(fftB,whitenedB);
	
	longTermA = PV_MagLagUD(PV_Copy(contrastA, LocalBuf(2048)), 0.2, 2.5);
	longTermB = PV_MagLagUD(PV_Copy(contrastB, LocalBuf(2048)), 0.2, 2.5);

	newA = PV_MagMinus(PV_Copy(contrastA, LocalBuf(2048)), longTermA);
	newB = PV_MagMinus(PV_Copy(contrastB, LocalBuf(2048)), longTermB);

	possible4A = PV_MagMinus(newB, PV_Copy(longTermA, LocalBuf(2048)));
	possible4B = PV_MagMinus(newA, PV_Copy(longTermB, LocalBuf(2048)));
	
	possible4A = PV_MagLagUD(PV_Copy(possible4A, LocalBuf(2048)), 0, 0.8);
	possible4B = PV_MagLagUD(PV_Copy(possible4B, LocalBuf(2048)), 0, 0.8);
	

	[IFFT.ar(possible4A), IFFT.ar(possible4B)] * 0.5;
}.play
)


// sensibleres system, function auf differenzen
// eigenverhalten mit kopplung




(
x = { arg speed = 1, n = 4;
	var contrast, contrast2, whitened, whitened2, past, past2, peaks2, new, new2,
	peaks,rndPhase,lagged, chain, chain2, input2, subtract, subtractIfft,
	input = HPF.ar(PlayBuf.ar(1,~c.bufnum,loop:1),100) * 1.5;
	input2 = DelayC.ar(input, 10, 5);
	chain = FFT(LocalBuf(4096), input, 0.5);
	chain2 = FFT(LocalBuf(4096), input2, 0.5);
	// whitened = PV_Whiten(PV_Copy(chain, LocalBuf(4096)), LocalBuf(4096));
	// whitened2 = PV_Whiten(PV_Copy(chain2, LocalBuf(4096)), LocalBuf(4096));
	// contrast = PV_MagMul(chain,whitened);
	// contrast2 = PV_MagMul(chain2,whitened2);
	past = PV_MagLagUD(PV_Copy(chain, LocalBuf(4096)), 0.5, 0.3);
	past2 = PV_MagLagUD(PV_Copy(chain2, LocalBuf(4096)), 0.5, 0.3);
	new = PV_MagMinus(PV_Copy(chain, LocalBuf(4096)), past);
	new2 = PV_MagMinus(PV_Copy(chain2, LocalBuf(4096)), past2);
	
	peaks = PV_MagPeaksDecay(new, 0.1, 1, 0.25*speed, 3*speed, 0, 120*speed, n, 1.8);
	peaks2 = PV_MagPeaksDecay(new2, 0.1, 1, 0.25*speed, 3*speed, 0, 120*speed, n, 1.8);

	subtract = PV_MagMinus(PV_Copy(peaks, LocalBuf(4096)), PV_Copy(peaks2, LocalBuf(4096)));
	// fade in new ones
	subtract = PV_MagLagUD(PV_Copy(subtract, LocalBuf(4096)), 0.4*speed, 0);
	subtractIfft = 	Limiter.ar(IFFT.ar(subtract)*3,0.9);
	subtractIfft = Compander.ar(subtractIfft, subtractIfft, 0.1, 1, 0.2) * 2;
	subtractIfft = HPF.ar(subtractIfft, 100);
	
	subtractIfft.dup; // + ([input,input2] * 0.5);
	//	IFFT.ar(subtract)*2;
}.play
)

x.set(\speed, 2);

s.meter;



// live input


(
x = { arg speed = 1, n = 4;
	var contrast, contrast2, whitened, whitened2, past, past2, peaks2, new, new2,
	peaks,rndPhase,lagged, chain, chain2, input2, subtract, subtractIfft,
	input = HPF.ar(SoundIn.ar(0),100) * 1.5;
	input2 = DelayC.ar(input, 10, 5);
	chain = FFT(LocalBuf(4096), input, 0.5);
	chain2 = FFT(LocalBuf(4096), input2, 0.5);
	// whitened = PV_Whiten(PV_Copy(chain, LocalBuf(4096)), LocalBuf(4096));
	// whitened2 = PV_Whiten(PV_Copy(chain2, LocalBuf(4096)), LocalBuf(4096));
	// contrast = PV_MagMul(chain,whitened);
	// contrast2 = PV_MagMul(chain2,whitened2);
	past = PV_MagLagUD(PV_Copy(chain, LocalBuf(4096)), 0.5, 0.3);
	past2 = PV_MagLagUD(PV_Copy(chain2, LocalBuf(4096)), 0.5, 0.3);
	new = PV_MagMinus(PV_Copy(chain, LocalBuf(4096)), past);
	new2 = PV_MagMinus(PV_Copy(chain2, LocalBuf(4096)), past2);
	
	// is uptime working with sustain?
	peaks = PV_MagPeaksDecay(new, 0.1, 1, 0.25*speed, 3*speed, 0, 100*speed, n, 1.8);
	peaks2 = PV_MagPeaksDecay(new2, 0.1, 1, 0.25*speed, 3*speed, 0, 100*speed, n, 1.8);

	subtract = PV_MagMinus(PV_Copy(peaks, LocalBuf(4096)), PV_Copy(peaks2, LocalBuf(4096)));
	// fade in new ones
	subtract = PV_MagLagUD(PV_Copy(subtract, LocalBuf(4096)), 0.3*speed, 0);

	subtractIfft = 	Limiter.ar(IFFT.ar(subtract)*3,0.9);
	subtractIfft = Compander.ar(subtractIfft, subtractIfft, 0.1, 1, 0.1) * 3;
	subtractIfft = HPF.ar(subtractIfft, 100);
	
	subtractIfft.dup; // + ([input,input2] * 0.5);
	//	IFFT.ar(subtract)*2;
}.play
)

s.meter

x.set(\speed, 1.5, \n, 7);

x.set(\speed, 0.1, \n, 7);

s.meter

/// octaves
/// possible techniques and transitions
/// time delta change



(
{
	var input, input2, chain, chain2, subtract, rndPhase;
	input = SinOsc.ar(300)*0.5;
	input2 = SinOsc.ar(300);
	chain = FFT(LocalBuf(4096), input, 1);
	chain2 = FFT(LocalBuf(4096), input2, 1);
	subtract = PV_MagMinusOct(chain,chain2);
	subtract = PV_CopyPhase(subtract, FFT(LocalBuf(4096), DC.ar(0)));
	(Limiter.ar(IFFT.ar(subtract),1) * 0.01).dup
}.play
)







(
x = { arg speed = 1, n = 2;
	var contrast, contrast2, whitened, whitened2, past, past2, peaks2, new, new2,
	peaks,rndPhase,lagged, chain, chain2, chainMul1, chainMul2, input2, subtract1, subtract2, subtractIfft,
	input = HPF.ar(PlayBuf.ar(1,~c.bufnum,loop:1),50)*1.5;
	input = Compander.ar(input,input, 0.1, 1, 0.1)*1.5;
	input2 = DelayC.ar(input, 5, 5);
	chain = FFT(LocalBuf(4096*2), HPF.ar(input,250), 0.5);
	chain2 = FFT(LocalBuf(4096*2), HPF.ar(input2,250), 0.5);
	chainMul1 = FFT(LocalBuf(4096*2), DC.ar(0));
	chainMul2 = FFT(LocalBuf(4096*2), DC.ar(0));
	whitened = PV_Whiten(PV_Copy(chain, LocalBuf(4096*2)), LocalBuf(4096*2));
	whitened2 = PV_Whiten(PV_Copy(chain2, LocalBuf(4096*2)), LocalBuf(4096*2));
	contrast = PV_MagMul(chain,whitened);
	contrast2 = PV_MagMul(chain2,whitened2);
	past = PV_MagLagUD(PV_Copy(chain, LocalBuf(4096*2)), 0.5, 0.3);
	past2 = PV_MagLagUD(PV_Copy(chain2, LocalBuf(4096*2)), 0.5, 0.3);
	new = PV_MagMinusOct(PV_Copy(contrast, LocalBuf(4096*2)), past);
	new2 = PV_MagMinusOct(PV_Copy(contrast2, LocalBuf(4096*2)), past2);
	
	peaks = PV_MagPeaksDecay(new, 0.1, 1, 0.4*speed, 3*speed, 0, 120*speed, n, 1.8);
	peaks2 = PV_MagPeaksDecay(new2, 0.1, 1, 0.4*speed, 3*speed, 0, 120*speed, n, 1.8);

	subtract1 = PV_MagMinusOct(PV_Copy(peaks, LocalBuf(4096*2)), PV_Copy(peaks2, LocalBuf(4096*2)));
	subtract2 = PV_MagMinusOct(PV_Copy(peaks2, LocalBuf(4096*2)), PV_Copy(peaks, LocalBuf(4096*2)));

	// fade in new ones
	subtract1 = PV_MagLagUD(PV_Copy(subtract1, LocalBuf(4096*2)), 0.2*speed, 0.2);
	subtract2 = PV_MagLagUD(PV_Copy(subtract2, LocalBuf(4096*2)), 0.2*speed, 0.2);


	chainMul1 = PV_FreqDiffs(chainMul1, subtract1, subtract2);
	chainMul2 = PV_FreqDiffs(chainMul2, subtract2, subtract1);

	chainMul1 = PV_MagLagUD(PV_Copy(chainMul1, LocalBuf(4096*2)), 0.2*speed, 0.2);
	chainMul2 = PV_MagLagUD(PV_Copy(chainMul2, LocalBuf(4096*2)), 0.2*speed, 0.2);
	
	subtractIfft = 	Limiter.ar([IFFT.ar(subtract1),IFFT.ar(subtract2),IFFT.ar(chainMul1),IFFT.ar(chainMul2)]*2,0.9)*2;
	subtractIfft = Compander.ar(subtractIfft, subtractIfft, 0.05, 1, 0.1) * 1.5;
	subtractIfft = HPF.ar(subtractIfft, 50);

	//	subtractIfft[2].poll;
	
	Out.ar(0, [subtractIfft[2] * 0.2, subtractIfft[3] * 0.2]);

	Out.ar(2, input*0.15); 

	Out.ar(3, [subtractIfft[0] * 0.2, subtractIfft[1] * 0.2]); 

}.play
)


s.reboot



~x = Buffer.alloc(s, 4096);

(
{
	var input, input2, chain, chain2, chain3, diffs;
	input = SinOsc.ar(11300)*1;
	input2 = SinOsc.ar(11900)*1;
	chain = FFT(LocalBuf(4096), input, 0.5);
	chain2 = FFT(LocalBuf(4096), input2, 0.5);
	chain3 = FFT(LocalBuf(4096), DC.ar(0));
	diffs = PV_FreqDiffs(chain3, chain, chain2);
	//	diffs = PV_CopyPhase(diffs, FFT(LocalBuf(4096), DC.ar(0)));
	(Limiter.ar(IFFT.ar(chain3).poll,0.05) * 0.01);
}.play
)

(
{
	var input, chain, analysis;
	//input = PlayBuf.ar(1, ~c.bufnum)*2;
	input = SinOsc.ar(MouseY.kr(90,700,1).poll);
	chain = FFT(LocalBuf(4096*2), input, 0.5);
	analysis = MagPeaksFreqs.kr(chain,3,0.1);
	analysis = OnOffChange.kr(analysis,[15,13,10],1.5,2);
	analysis.poll;
	DC.ar(0);
	//	[ input, SinOsc.ar(WrapOct.kranalysis[0],0,analysis[1]).sum * 0.2]
}.play
)

6989.16 / 1299

44100/((4096*2)-1) * 1299

{ SinOsc.ar(700) * 0.3 }.play

(
{
	var in = LFNoise0.kr(7).exprange(100,400);
	var onTime = 1, offTime = 4, ampLag = 3;
	
	var trigOn = Trig1.kr(Changed.kr(in),onTime);
	var trigOnFreq = Trig1.kr(trigOn,onTime+(ampLag*1.5)).poll;
	var delTrig = 1-Trig1.kr(Changed.kr(TDelay.kr(trigOn,onTime)),offTime+onTime);
	
	var gated = Gate.kr(in,1-trigOnFreq);
	SinOsc.ar(gated) * 0.05 * (trigOn*delTrig).lag3(ampLag);
}.play
)


s.reboot;
(
{
	var ch = OnOffChange.kr(LFNoise0.kr(20!5).exprange(150,300),Rand(1,5!5),Rand(1,8!5),Rand(0.5,3!5));
	Splay.ar(SinOsc.ar(ch[0]) * ch[1]) * 0.2;
}.play
)


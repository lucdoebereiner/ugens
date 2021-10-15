#include "SC_fftlib.h"
#include "FFT_UGens.h"


#define PV_GET_BUF1 \
	float fbufnum = ZIN0(0); \
	uint32 ibufnum = (uint32)fbufnum; \
	World *world = unit->mWorld; \
	SndBuf *buf; \
	if (ibufnum >= world->mNumSndBufs) { \
		int localBufNum = ibufnum - world->mNumSndBufs; \
		Graph *parent = unit->mParent; \
		if(localBufNum <= parent->localBufNum) { \
			buf = parent->mLocalSndBufs + localBufNum; \
		} else { \
			buf = world->mSndBufs; \
		} \
	} else { \
		buf = world->mSndBufs + ibufnum; \
	} \
	LOCK_SNDBUF(buf); \
	int numbins = (buf->samples - 2) >> 1;	\


#define PV_GET_BUF3 \
	float fbufnum1 = ZIN0(0); \
	float fbufnum2 = ZIN0(1); \
	float fbufnum3 = ZIN0(2);					\
	if (fbufnum1 < 0.f || fbufnum2 < 0.f || fbufnum3 < 0.f) { ZOUT0(0) = -1.f; return; } \
	ZOUT0(0) = fbufnum1; \
	uint32 ibufnum1 = (int)fbufnum1; \
	uint32 ibufnum2 = (int)fbufnum2; \
	uint32 ibufnum3 = (int)fbufnum3; \
	World *world = unit->mWorld; \
	SndBuf *buf1; \
	SndBuf *buf2; \
	SndBuf *buf3; \
	if (ibufnum1 >= world->mNumSndBufs) { \
		int localBufNum = ibufnum1 - world->mNumSndBufs; \
		Graph *parent = unit->mParent; \
		if(localBufNum <= parent->localBufNum) { \
			buf1 = parent->mLocalSndBufs + localBufNum; \
		} else { \
			buf1 = world->mSndBufs; \
		} \
	} else { \
		buf1 = world->mSndBufs + ibufnum1; \
	} \
	if (ibufnum2 >= world->mNumSndBufs) { \
		int localBufNum = ibufnum2 - world->mNumSndBufs; \
		Graph *parent = unit->mParent; \
		if(localBufNum <= parent->localBufNum) { \
			buf2 = parent->mLocalSndBufs + localBufNum; \
		} else { \
			buf2 = world->mSndBufs; \
		} \
	} else { \
		buf2 = world->mSndBufs + ibufnum2; \
	} \
	if (ibufnum3 >= world->mNumSndBufs) { \
		int localBufNum = ibufnum3 - world->mNumSndBufs; \
		Graph *parent = unit->mParent; \
		if(localBufNum <= parent->localBufNum) { \
			buf3 = parent->mLocalSndBufs + localBufNum; \
		} else { \
			buf3 = world->mSndBufs; \
		} \
	} else { \
		buf3 = world->mSndBufs + ibufnum3; \
	} \
	LOCK_SNDBUF2(buf1, buf2); \
	LOCK_SNDBUF(buf3); \
	if (buf1->samples != buf2->samples) return; \
	int numbins = (buf1->samples - 2) >> 1;




struct PV_PollMags : public PV_Unit {
  double last_poll;
};



struct PV_MagLagUD : public PV_Unit {
  double* m_y1s;
  float m_lagu, m_lagd;
  double m_b1u, m_b1d, m_y1;
};


struct PV_MagPeaksDecay : public PV_Unit {
  double* slopes;
  double* m_y1s;
  int* sustain_counters;
  int n_sustain;
  double* m_lagds;
  double* rnd_phases;
  double lagd_decay;
};


struct MagPeaksFreqs : public Unit {
  int nFreqs;
  float binFreqStep;
  double* slopes;
  double* mags;
  int* lastFreqIdxs;
};


struct WrapOct : public Unit {
  
};


extern "C" {
  void PV_MagLagUD_Ctor(PV_MagLagUD *unit);
  void PV_MagLagUD_next(PV_MagLagUD *unit, int inNumSamples);

  void PV_MagPeaksDecay_Ctor(PV_MagPeaksDecay *unit);
  void PV_MagPeaksDecay_next(PV_MagPeaksDecay *unit, int inNumSamples);

  void PV_MagMinusOct_Ctor(PV_Unit *unit);
  void PV_MagMinusOct_next(PV_Unit *unit, int inNumSamples);

  void PV_FreqDiffs_Ctor(PV_Unit *unit);
  void PV_FreqDiffs_next(PV_Unit *unit, int inNumSamples);

  void MagPeaksFreqs_Ctor(MagPeaksFreqs *unit);
  void MagPeaksFreqs_next(MagPeaksFreqs *unit, int inNumSamples);
  void MagPeaksFreqs_Dtor(MagPeaksFreqs *unit);

  void WrapOct_next_k(WrapOct *unit, int inNumSamples);
  void WrapOct_Ctor(WrapOct* unit);

  void PV_PollMags_Ctor(PV_PollMags *unit);
  void PV_PollMags_next(PV_PollMags *unit, int inNumSamples);

}



inline int closesDiffFreqIdx(int idx1, int idx2, int numbins, int sr) {

  double binFreqStep = (double) sr/ (double) numbins;

  double fr1 = idx1 * binFreqStep;
  double fr2 = idx2 * binFreqStep;

  double diffFreqs[6] = { fr1 - fr2, fr2 - fr1, (2 * fr2) - fr1, (2 * fr1) - fr2, (2 * fr1) - (2 * fr2), (2 * fr2) - (2 * fr1) };

  int smallestDiffIdx = 0;
  int smallestDiff = abs(fr1 - abs(diffFreqs[0]));

  for (int i = 1; i < 6; i++) {
    double diff = abs(fr1 - abs(diffFreqs[i]));
    if ((diff < smallestDiff) && (abs(diffFreqs[i]) < 10000)) {
      smallestDiff = diff;
      smallestDiffIdx = i;
    }
  };
  
  // Print("idx1: %d, step: %f, fr1: %f, fr2: %f, smallestDiff: %f, diffidx: %d\n", idx1, binFreqStep, fr1, fr2, smallestDiff, smallestDiffIdx);
  
  return  round(abs(diffFreqs[smallestDiffIdx]) / binFreqStep);
  
}



void PV_FreqDiffs_next(PV_Unit *unit, int inNumSamples)
{
	
  PV_GET_BUF3


  
	SCPolarBuf *out = ToPolarApx(buf1);
	SCPolarBuf *in1 = ToPolarApx(buf2);
	SCPolarBuf *in2 = ToPolarApx(buf3);

	float thresh = 0.0001;
	int in2Idx = 4;
	int diffFreqIdx = 0;


	// set all to zero first
	for (int i=0; i<numbins; ++i) {
	  out->bin[i].phase = in1->bin[i].phase;
	}
	
	for (int i=0; i<numbins; ++i) {
	  
	  float thisMag = in1->bin[i].mag;
	  bool in2MagFound = false;

	  
	  if ((thisMag > thresh) && (i > 4) && (i < numbins - 2)) {

	    
	    // find next relevant bin in in2
	    while (!in2MagFound && (in2Idx < numbins)) {
	      if (in2->bin[in2Idx].mag > thresh) {
		in2MagFound = true;
	      } else {
		in2Idx++;
	      }
	    };

	    if  (in2MagFound) {
	      diffFreqIdx = closesDiffFreqIdx(i, in2Idx, numbins, 44100);
	      out->bin[diffFreqIdx].mag = thisMag;
	    } 

	    in2Idx++;
	  
	  };

	}

}



void PV_MagMinusOct_next(PV_Unit *unit, int inNumSamples)
{
	PV_GET_BUF2

	SCPolarBuf *p = ToPolarApx(buf1);
	SCPolarBuf *q = ToPolarApx(buf2);

	
	for (int i=0; i<numbins; ++i) {

       	  // subtract same freqs
	  p->bin[i].mag = sc_max(p->bin[i].mag - q->bin[i].mag, 0.f);

	  // subtract 4 octaves up and down
	  for (int k=1; k<4; ++k) {

	    int octFac = k * 2;
	    int lowerOctIdx = i / octFac;
	    int upperOctIdx = i * octFac;

	    if (lowerOctIdx > 0) {
	      p->bin[lowerOctIdx].mag = sc_max(p->bin[lowerOctIdx].mag - q->bin[i].mag, 0.f);
	    };

	    if ((upperOctIdx > 0) && (upperOctIdx < numbins)) {
	      p->bin[upperOctIdx].mag = sc_max(p->bin[upperOctIdx].mag - q->bin[i].mag, 0.f);
	    };

	  };

	}
}

void PV_MagMinusOct_Ctor(PV_Unit *unit)
{
	SETCALC(PV_MagMinusOct_next);
	ZOUT0(0) = ZIN0(0);
}


void PV_FreqDiffs_Ctor(PV_Unit *unit)
{
  SETCALC(PV_FreqDiffs_next);
  ZOUT0(0) = ZIN0(0);
}






void PV_MagLagUD_next(PV_MagLagUD *unit, int inNumSamples) {

  PV_GET_BUF

  SCPolarBuf *p = ToPolarApx(buf);

  double b1u = unit->m_b1u;
  double b1d = unit->m_b1d;
  
  float lagu= ZIN0(1);
  float lagd= ZIN0(2);

  double *y1s = unit->m_y1s;
  double y1;

  // Update lag parameters
  if ( !((lagu == unit->m_lagu) && (lagd == unit->m_lagd)) ) {

    unit->m_b1u = lagu == 0.f ? 0.f : exp(log001 / (lagu * numbins));
    unit->m_lagu = lagu;
    unit->m_b1d = lagd == 0.f ? 0.f : exp(log001 / (lagd * numbins));
    unit->m_lagd = lagd;
  }
  
  for (int i=0; i<numbins; ++i) {
    float y0 = p->bin[i].mag;
    y1 = y1s[i];
    
    if ( y0 > y1 )
      p->bin[i].mag = y1 = y0 + b1u * (y1 - y0);
    else
      p->bin[i].mag = y1 = y0 + b1d * (y1 - y0);

    y1s[i] = zapgremlins(y1);
  }
  
}



void PV_MagLagUD_Ctor(PV_MagLagUD *unit) {
  SETCALC(PV_MagLagUD_next);

  ZOUT0(0) = ZIN0(0);

  PV_GET_BUF
    
  unit->m_y1s = (double*)RTAlloc(unit->mWorld, numbins * sizeof(double));
    
}


double mean(double *array, int n) {
  double sum = 0;
  for (int i = 0; i<n; ++i) {
    sum += array[i];
  }
  return sum / (double) n;
}


inline int sgn(double x) {
    if (x==0)
        return 0;
    else
        return (x>0) ? 1 : -1;
}


inline double linterp(double x, double a, double b) {
	return a + x * (b - a);
}

// for testing purposes, get the idx of max element in array
int maxIdx(double* array, int n) {
  double max = 0.0;
  int idx = -1;
  for (int i = 0; i<n; ++i) {
    if (array[i] > max) {
      max = array[i];
      idx = i;
    }
  }
  return idx;
}

inline bool distance_to_next_sustained(int i, int* sustain_counters, int n, int lowerLimit, int upperLimit) {
  int distUp = 0;
  int distLow = 0;
  int upper, lower;
  //Print("i: %d, n: %d\n", i, n);
  while (true) {
    upper = sc_clip(i+distUp, 0, n-1);
    lower = sc_clip(i+distLow, 0, n-1);
    // we are outside the limits, so min dist is ok
    if ((distUp >= upperLimit) && (distLow <= lowerLimit)) {
      return true;
    };
    // we found a sustained bin
    if ((sustain_counters[upper] > 0) || (sustain_counters[lower] > 0)) {
      // Print("i %d, lowerLimit: %d, upperLimit %d, distLow: %d, distUp: %d\n", i, lowerLimit, upperLimit, distLow, distUp);
      return false;
    };
    if (distUp < upperLimit) {
      distUp++;
    };
    if (distLow > lowerLimit) {
      distLow--;
    };    
    // we checked the entire array, should never happen
    if (((i+distUp) >= n) && ((i+distLow) < 0)) {
      return false;
    };
  }
}


int find_lowest_sustained(double* mags, int* counters, int n) {
  int lowestIdx = 0;
  double lowestMag = 0.0;
  for (int i = 0; i<n; i++) {
    if (counters[i] > 0) {
      if (mags[i] < lowestMag) {
	lowestMag = mags[i];
	lowestIdx = i;
      }
    }
  }
  return lowestIdx;
}

void PV_MagPeaksDecay_next(PV_MagPeaksDecay *unit, int inNumSamples) {
  
  PV_GET_BUF

  SCPolarBuf *p = ToPolarApx(buf);
  
  double *slopes = unit->slopes;

  double binFreqStep = (double) 44100 / (double) numbins;
  
  float minPeak = ZIN0(1);
  float maxPeak = ZIN0(2);
  float minDec = ZIN0(3);
  float maxDec = ZIN0(4);
  float upTime = ZIN0(5);
  float sustainC = ZIN0(6);
  int max_n_sustain = ZIN0(7);
  float minDistFac = ZIN0(8); // factor
  
  double *y1s = unit->m_y1s;
  double y1;

  int meanWindowSize = 10;
  double meanArray[meanWindowSize] = { 0. };

  // validate input
  max_n_sustain = sc_clip(max_n_sustain, 1, numbins);
  
  // create mean
  for (int i=0; i<numbins; ++i) {
    meanArray[i % meanWindowSize] = p->bin[i].mag;
    slopes[i] = mean(meanArray, meanWindowSize);
  }

  // derivative
  double last = 0;
  for (int i=0; i<numbins; ++i) {
    double prev = slopes[i];
    slopes[i] = prev - last;
    last = prev;
  }

  // velocity at zero crossing
  double x1, x2, x3;
  x1 = x2 = x3 = std::numeric_limits<double>::quiet_NaN();

  for (int i=0; i<numbins; ++i) {
    double ret = 0.;
    double inval = slopes[i];
    if ( !(isnan(x1) || isnan(x2) || isnan(x3)) ) {
      if (sgn(x2) != sgn(x1)) {
	ret = ((x2 - x3) + (x1 - x2) + (inval - x1)) / 3.0;
      } else {
	ret = 0.;
      }
    };
    x3 = x2;
    x2 = x1;
    x1 = inval;
    slopes[i] = ret;    

  }


  // double maxMag = 0.0;
  // int idx = -1;
  // for (int i = 1; i<numbins; ++i) {
  //   double curMag = p->bin[i].mag;
  //   if (curMag > maxMag) {
  //     maxMag = curMag;
  //     idx = i;
  //   }
  // }
  // Print("maxMag %f, idx: %d\n", maxMag, idx);
      
  // lag 
  double* lagds = unit->m_lagds;
  int* sustain_counters = unit->sustain_counters;
  double lagd_decay = unit->lagd_decay;
  double* rnd_phases = unit->rnd_phases;
    
  for (int i=0; i<numbins; ++i) {

    double range = maxPeak-minPeak;
    double lagu = upTime;
    double cur_slope = 0, lagd, b1u, b1d;
    double binFreq = binFreqStep * i;
    
    double y0 = p->bin[i].mag;
    y1 = y1s[i];

    
    // the actual peak velocity of current bin is shifted by meanWindowsSize/2 + 1
    // without + 1 seems to hit the max
    int velocityOffset = meanWindowSize/2; // + 1;

    if (sustain_counters[i] > 0) {
      if (y0 > y1) {
	y1s[i] = y0;
	sustain_counters[i] = sustainC;
	//	Print("i: %d\n", i);
      } else {
	p->bin[i].mag = y1;
	sustain_counters[i] = sustain_counters[i] - 1;
	if (sustain_counters[i] == 0) {
	  unit->n_sustain--;
	}
      }
    } else {
    
      // get current slope and lagd
      if (( i + velocityOffset ) > (numbins-1)) {
	lagd = 0.0;
      } else {
	cur_slope = sc_abs(sc_min(slopes[i+velocityOffset], 0));
	if (cur_slope < minPeak) {
	  lagd = 0;
	} else if (cur_slope > maxPeak) {
	  lagd = maxDec;
	} else {
	  lagd = linterp( sc_clip((cur_slope-minPeak) / range, 0.0, 1.0), minDec, maxDec);
	}
      }

      
      // ensure distance to next sustained bin
      double maxFreqLimit = minDistFac > 1.0 ? binFreq * minDistFac: binFreq / minDistFac;
      double minFreqLimit = minDistFac > 1.0 ? binFreq / minDistFac: binFreq *  minDistFac;

      int upperLimit = (int) (maxFreqLimit - binFreq) / binFreqStep;
      int lowerLimit = (int) (minFreqLimit - binFreq) / binFreqStep; // will be negative

      // if ((i % 500) == 0) {
      // 	Print("i: %d, minDistFac %f, minFreqLimit: %f, lowerLimit: %d, maxFreqLimit: %f, upperLimit: %d\n", i, minDistFac, minFreqLimit, lowerLimit, maxFreqLimit, upperLimit);
      // }
      
      bool distance_to_next = false; 
      bool possibleSlope = (cur_slope > minPeak);
      bool possibleSustain = (unit->n_sustain < max_n_sustain);
      // if ((i % 500) == 0) {
      // 	Print("distance: %d\n", distance_to_next);
      // }
      
      // adjust mag according to peakedness and remove non-peaked
      if (possibleSlope) {
	distance_to_next = distance_to_next_sustained(i, sustain_counters, numbins, lowerLimit, upperLimit);
      }

      // see if one of the sustained is smaller than the candidate
      if (distance_to_next && possibleSlope && !possibleSustain) {
	int smallest = find_lowest_sustained(y1s, sustain_counters, numbins);
	if (smallest > 0) { // 0 is returned if no smallest is found
	  if (y1s[smallest] < y0) {
	    sustain_counters[smallest] = 0;
	    y1s[smallest] = 0;
	    possibleSustain = true;
	  }
	}
      }

      
      if (distance_to_next && possibleSlope && possibleSustain) {
	p->bin[i].mag =
	  sc_clip(p->bin[i].mag * linterp( sc_clip((cur_slope-minPeak) / range, 0.0, 1.0), 0.5, 1.5), 70, 200);
	sustain_counters[i] = sustainC;
	y1s[i] = p->bin[i].mag;
	//	y1s[i] = y0; // ?
	unit->n_sustain++;
      }
      else {
	p->bin[i].mag = 0;
	y0 = 0;
      }
    
      b1u = lagu == 0.f ? 0.f : exp(log001 / (lagu * numbins));
      b1d = lagd == 0.f ? 0.f : exp(log001 / (lagd * numbins));
    
      // lag of lag times
      double last_decay = lagds[i];
      if ( last_decay > b1d ) {
	b1d = b1d + lagd_decay * (last_decay - b1d);
      };
    
      lagds[i] = zapgremlins(b1d);

      // set phase to static rnd value
      p->bin[i].phase = rnd_phases[i];
      
      if ( y0 > y1 ) 
  	p->bin[i].mag = y1 = y0 + b1u * (y1 - y0);
      else
  	p->bin[i].mag = y1 = y0 + b1d * (y1 - y0);
          
      y1s[i] = zapgremlins(y1);

    }
    
  }

}



void PV_MagPeaksDecay_Ctor(PV_MagPeaksDecay *unit) {
  SETCALC(PV_MagPeaksDecay_next);

  PV_GET_BUF
  RGen& rgen = *unit->mParent->mRGen;
    
  unit->slopes = (double*)RTAlloc(unit->mWorld, numbins * sizeof(double));
  unit->m_y1s = (double*)RTAlloc(unit->mWorld, numbins * sizeof(double));
  unit->m_lagds = (double*)RTAlloc(unit->mWorld, numbins * sizeof(double));
  unit->sustain_counters = (int*)RTAlloc(unit->mWorld, numbins * sizeof(double));
  unit->rnd_phases = (double*)RTAlloc(unit->mWorld, numbins * sizeof(double));
  unit->n_sustain = 0;
  for (int i = 0; i < numbins; ++i) {
    unit->m_lagds[i] = 0;
    unit->m_y1s[i] = 0;
    unit->sustain_counters[i] = 0;
    unit->rnd_phases[i] = rgen.frand() * twopi;
  }
  unit->lagd_decay = exp(log001 / (8 * numbins));

  ZOUT0(0) = ZIN0(0);
}




void MagPeaksFreqs_next(MagPeaksFreqs *unit, int inNumSamples) {
  
  PV_GET_BUF1

  SCPolarBuf *p = ToPolarApx(buf);
  
  double *slopes = unit->slopes;



  float minPeak = ZIN0(2);
  
  float binFreqStep = unit->binFreqStep;

  int *lastFreqIdxs = unit->lastFreqIdxs;
  double *mags = unit->mags;
  int nFreqs = unit->nFreqs;
  //  int maxSlopeIdxs[nFreqs] = { 0 };
  int *maxSlopeIdxs = lastFreqIdxs;
  
  int meanWindowSize = 10;
  double meanArray[meanWindowSize] = { 0. };
  int velocityOffset = meanWindowSize/2; // + 1;
  
//  float** out = unit->mOutBuf;
  
  // create mean if analysis is there
  if (!(fbufnum < 0.f)) {

      for (int i=0; i<numbins; ++i) {
	meanArray[i % meanWindowSize] = p->bin[i].mag;
	slopes[i] = mean(meanArray, meanWindowSize);

	// if (p->bin[i].mag > 20) {
	//   Print("i: %d, mag: %f, freq: %f\n", i, p->bin[i].mag, binFreqStep * i);
	// }
      }

      // derivative
      double last = 0;
      for (int i=0; i<numbins; ++i) {
	double prev = slopes[i];
	slopes[i] = prev - last;
	last = prev;
      }

      // velocity at zero crossing
      double x1, x2, x3;
      x1 = x2 = x3 = std::numeric_limits<double>::quiet_NaN();

      for (int i=0; i<numbins; ++i) {
	double ret = 0.;
	double inval = slopes[i];
	if ( !(isnan(x1) || isnan(x2) || isnan(x3)) ) {
	  if (sgn(x2) != sgn(x1)) {
	    ret = ((x2 - x3) + (x1 - x2) + (inval - x1)) / 3.0;
	  } else {
	    ret = 0.;
	  }
	};
	x3 = x2;
	x2 = x1;
	x1 = inval;
	slopes[i] = ret;    
      }

      // get the biggest n slope indices
      for (int i=0; (i + velocityOffset) <numbins; ++i) {
	float curSlope = sc_abs(sc_min(slopes[i+velocityOffset], 0));

	float smallestMax = sc_abs(sc_min(slopes[maxSlopeIdxs[0]+velocityOffset], 0));;
	int smallestMaxIdx = 0;
	bool alreadyMember = false;

	for (int k = 1; k < nFreqs; ++k) {
	  if (sc_abs((maxSlopeIdxs[k] - i)) < 3) {
	    alreadyMember = true;
	  }
	}

	for (int k = 1; k < nFreqs; ++k) {
	  float maxCandidate = sc_abs(sc_min(slopes[maxSlopeIdxs[k]+velocityOffset], 0));
	  if (maxCandidate < smallestMax) {
	    smallestMax = maxCandidate;
	    smallestMaxIdx = k;
	  }
	}
	
	  
	if ((curSlope >= minPeak) && (curSlope > smallestMax) && !alreadyMember && (p->bin[i].mag > 10)) {
	  //	  Print("curslope: %f, smallestmax: %f\n", curSlope, smallestMax);
	    maxSlopeIdxs[smallestMaxIdx] = i;
	    mags[smallestMaxIdx] = p->bin[i].mag;
	  }
      }
      
    }

  for (int k=0; k < inNumSamples; ++k) {
    for (int i=0; i<(nFreqs*2); ++i) {
    
      float thisFreq = ((float) maxSlopeIdxs[i]) * binFreqStep;

      // if (sc_abs(sc_min(slopes[maxSlopeIdxs[i]+velocityOffset], 0)) > minPeak) {
      // 	lastFreqIdxs[i] = maxSlopeIdxs[i];
      // 	thisFreq = ((float) maxSlopeIdxs[i]) * binFreqStep;
      // }

      // Print("i: %d, k: %d,lastIdx: %d, maxSlopeIdx: %d, minPeak: %f, step: %f, thisFreq: %f\n", i, k, lastFreqIdxs[i], maxSlopeIdxs[i], minPeak, binFreqStep, thisFreq);

      if (i < nFreqs) {
	OUT(i)[k] = thisFreq;
      } else {
	OUT(i)[k] = mags[i-nFreqs];
      }
    }
  }

}



void MagPeaksFreqs_Ctor(MagPeaksFreqs *unit) {
  SETCALC(MagPeaksFreqs_next);

  PV_GET_BUF1

  unit->nFreqs = IN0(1);     
  unit->slopes = (double*)RTAlloc(unit->mWorld, numbins * sizeof(double));
  unit->mags = (double*)RTAlloc(unit->mWorld,  unit->nFreqs * sizeof(double));
  unit->lastFreqIdxs = (int*)RTAlloc(unit->mWorld, unit->nFreqs * sizeof(int));
  unit->binFreqStep = (double) 44100 / (((double) numbins + 1) * 2);
  for (int i = 0; i < unit->nFreqs; ++i) {
    unit->lastFreqIdxs[i] = 0;
    unit->mags[i] = 0;
  }

  MagPeaksFreqs_next(unit, 1);
}



void MagPeaksFreqs_Dtor(MagPeaksFreqs *unit) {

  RTFree(unit->mWorld, unit->slopes);
  RTFree(unit->mWorld, unit->lastFreqIdxs);
  RTFree(unit->mWorld, unit->mags);

}



void WrapOct_Ctor(WrapOct* unit)
{
  SETCALC(WrapOct_next_k);
  WrapOct_next_k(unit, 1);
}


void WrapOct_next_k(WrapOct *unit, int inNumSamples) {
  float *out = OUT(0);
  float inFreq = IN0(0);
  float minFreq = IN0(1);
  float maxFreq = IN0(2);
  
  
  for (int i=0; i < inNumSamples; ++i) {
    while ((inFreq > maxFreq) || (inFreq < minFreq))  {
      if (inFreq > maxFreq) {
	inFreq = inFreq * 0.5f;
      } else {
	inFreq = inFreq * 2.f;
      }
    }
    out[i] = inFreq;
  }
}


//// PollMags


void PV_PollMags_next(PV_PollMags *unit, int inNumSamples) {

  PV_GET_BUF

  SCPolarBuf *p = ToPolarApx(buf);

  float this_poll = ZIN0(1);

  double last_poll = unit->last_poll;

  // Update lag parameters
  if ((last_poll <= 0.0) && (this_poll > 0)) { 
  
    for (int i=0; i<numbins; ++i) {
      float m = p->bin[i].mag;
      Print("i: %d, mag: %.12f\n", i, m);      
    };
  };
  
}

void PV_PollMags_Ctor(PV_PollMags *unit) {
  SETCALC(PV_PollMags_next);

  ZOUT0(0) = ZIN0(0);

  PV_GET_BUF
    
    unit->last_poll = 0.0;
    
}



#define DefinePVUnit(name)						\
	(*ft->fDefineUnit)(#name, sizeof(PV_Unit), (UnitCtorFunc)&name##_Ctor, 0, 0);


InterfaceTable *ft;

PluginLoad(SpecTrans)
{
  ft = inTable;
  
  DefinePVUnit(PV_MagLagUD);
  DefinePVUnit(PV_MagPeaksDecay);
  DefinePVUnit(PV_MagMinusOct);
  DefinePVUnit(PV_FreqDiffs);
  DefinePVUnit(PV_PollMags);
  DefineDtorUnit(MagPeaksFreqs);
  DefineSimpleUnit(WrapOct)
}


// TODO
// write destructors

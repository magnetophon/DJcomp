
declare name "DJcomp";
declare version "0.1";
declare author "Bart Brouns bart@magnetophon.nl";
declare license "AGPLv3";

import("stdfaust.lib");


process =
  DJcomp;

DJcomp(l,r) =
  (gain_computer_lin(l,r)
   <:(_*l,_*r,_)
  )~(_,_);

gain_computer_lin(l,r,FBl,FBr) =
  level
  : gain_computer_non_smoothed_db(1,threshold,knee)
  : ba.db2linear
    // : smootherCascade(N, ma.EPSILON, rel)
  : si.onePoleSwitching(sdRel,0)
    // : si.onePoleSwitching(rel*0.25,0)
    // : si.onePoleSwitching(rel*0.25,0)
with {
  level =
    max(abs(l),abs(r))
    : ba.linear2db;
  FBlevel =
    // hslider("FBlevel", 0, -70, 70, 0.1);
    max(abs(FBl),abs(FBr))
    : ba.linear2db;
  sdRel =
    rel/
    select2(checkbox("slowdown")
           , 1
           , (
             (((FBlevel-threshold)*-1/sdKnee:min(1):max(ma.EPSILON)))
             : si.onePoleSwitching(rel,SDatt)
               // :pow(hslider("power", 1, 0.1, 5, 0.1))
               // :hbargraph("slow down", 0, 1)
           ))
  ;
  SDatt =
    select2(FBlevel<threshold
           , 0
           ,0.42*0.001);
  // , hslider("sd att", 0.42, 0, 200, 0.01)*0.001);
  N=4;
  T = ma.T;
  PI = ma.PI;
  TWOPI = 2.0 * PI;
  TWOPIT = TWOPI * T;
  /* Cascaded one-pole smoothers with attack and release times. */
  smoother(N, att, rel, x) = loop ~ _
  with {
    loop(fb) = coeff * fb + (1.0 - coeff) * x
    with {
    cutoffCorrection = 1.0 / sqrt(pow(2.0, 1.0 / N) - 1.0);
    coeff = ba.if(x > fb, attCoeff, relCoeff);
    TWOPITC = TWOPIT * cutoffCorrection;
    attCoeff = exp(-TWOPITC / att);
    relCoeff = exp(-TWOPITC / rel);
  };
  };
  smootherCascade(N, att, rel, x) = x : seq(i, N, smoother(N, att, rel));
};


// strength goes from  0 to 1, the rest is in dB
gain_computer_non_smoothed_db(strength,thresh,knee,level) =
  select3((level>(thresh-(knee/2)))+(level>(thresh+(knee/2))),
          0,
          ((level-thresh+(knee/2)) : pow(2)/(2*max(ma.EPSILON,knee))),
          (level-thresh))
  : max(0)*-strength;

//*****************************************************************************
//                 GUI
//*****************************************************************************

threshold = hslider("[03]threshold", 0, -70, 0, 0.1) : si.smoo ;
knee = hslider("knee", 0, 0, 30, 0.1) : si.smoo;
sdKnee = hslider("slow down knee", ma.EPSILON, ma.EPSILON, 30, 0.1) : si.smoo;
// threshold = hslider("threshold", 0, -70.0, 0.0, 0.01) : si.smoo: ba.db2linear ;
// knee = hslider("knee", 0.5, 0.0, 1.0, 0.001) : si.smoo;
pre_gain = hslider("pre-gain", 0, -30.0, 30.0, 0.01) : si.smoo: ba.db2linear ;
// pre_gain = hslider("pre-gain", ba.linear2db(ma.PI/2), -30.0, 30.0, 0.01) : ba.db2linear : si.smoo;

rel = hslider("release", 80, 0, 500, 0.1)*0.001;

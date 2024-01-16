declare name "parallelLim";
declare version "0.1";
declare author "Bart Brouns bart@magnetophon.nl";
declare license "AGPLv3";

import("stdfaust.lib");

process(l,r) =
  DJcomp(l,r);
// (l,r):co.FFcompressor_N_chan(1,threshold,0.01,0.04,knee,1,1,_,2);

DJcomp(l,r) =
  max(abs(l),abs(r)):gain_computer_mono
                     <:(_*pre_gain*l,_*pre_gain*r,_) ;

soft_clipper(pre_gain, threshold, knee, x) =
  pregain
  : clip
  : nonlin
  : postgain
with {
  pregain = (x / threshold) * pre_gain;
  abs_pregain = abs(pregain);
  sign = -1, 1 : select2(x >= 0.0);
  clip(in) =
    in:max(-1):min(1);
  nonlin(in) =
    in,
    ((1.0 - knee
      + sin( 1 - (1 -abs_pregain )/knee)
      * knee) * sign)
    : select2((abs_pregain <= breakpoint) & (abs_pregain >= (1.0 - knee)));
  postgain(in) = in * threshold;
  breakpoint = 1.0 + knee *(ma.PI/2-1);
};

gain_computer_mono =
  // :parMean
  // ba.linear2db
  _*pre_gain
  <: parGainComputer
  : ba.parallelMin(N)
  : si.onePoleSwitching(hslider("[14]post rel", 0.013, 0, 0.1, 0.001),0)
  : hbargraph("[13]GR", -24, 0)
  : ba.db2linear
with {
  parGainComputer =
    (
      // pregains,
      strengths,thresholds,attacks,releases,knees,prePosts, si.bus(N))
    : ro.interleave(N,7)
    : par(i, N, peak_compression_gain_mono_db:hgroup("[12]meters", vbargraph("%i", -24, 0)))
  ;
  // peak_compression_gain_mono_db(strength,thresh,att,rel,knee,prePost)

  peak_compression_gain_mono_db(strength,thresh,att,rel,knee,prePost) =
    abs : ba.bypass1(prePost,si.onePoleSwitching(att,rel)) : ba.linear2db : gain_computer(strength,thresh,knee): ba.db2linear : ba.bypass1((prePost !=1),si.onePoleSwitching(rel,att): ba.linear2db)
  with {
    gain_computer(strength,thresh,knee,level) =
      select3((level>(thresh-(knee/2)))+(level>(thresh+(knee/2))),
              0,
              ((level-thresh+(knee/2)) : pow(2)/(2*max(ma.EPSILON,knee))),
              (level-thresh))
      : max(0)*-strength;
  };
  pregains = par(i, N, pre_gain);
  strengths =
    shapedArray(
      hslider("[01]strength",0, 0, 1, 0.001)
     ,1
     , hslider("[02]strength shape", 0, -1, 1, 0.001 ),N)
    : ro.cross(N)
  ;
  thresholds =
    // par(i, N, threshold);
    shapedArray(0,hslider("[03]end threshold", 0, -30, 0, 0.1), hslider("[04]threshold shape", 0, -1, 1, 0.001) ,N);
  attacks =
    shapedArray(0,hslider("[05]end attack", 0.1, 0, 0.25, 0.001), hslider("[06]attack shape", 0, -1, 1, 0.001) ,N);
  releases=
    shapedArray(hslider("[07]start rel", 0.1, 0, 0.5, 0.001),hslider("[08]end rel", 0.1, 0, 0.5, 0.001), hslider("[09]rel shape", 0, -1, 1, 0.001) ,N);
  knees =
    shapedArray(ma.EPSILON,hslider("[10]end knee", 30, 0, 70, 0.1), hslider("[11]knee shape", 0, -1, 1, 0.001) ,N);
  prePosts = par(i, N, 1);
  //par(i, N, checkbox("[14]prepost"));
  // par(i, N+1, 0.5);
  // gain_computer(pre_gain, threshold, strength, knee, x) =
  // (x*pre_gain-soft_clipper(pre_gain, threshold, knee, x))
  // : max(0)*-strength;
  gain_computer(thresh,strength,knee,level) =
    select3((level>(thresh-(knee/2)))+(level>(thresh+(knee/2))),
            0,
            ((level-thresh+(knee/2)) : pow(2)/(2*max(ma.EPSILON,knee))),
            (level-thresh))
    : max(0)*-strength;
}
;

shapedArray(bottom,top,shape,0) =   0:! ;
shapedArray(bottom,top, shape ,nrElements) =
  par(i,nrElements,
      (i/(nrElements-1))
      :shaper(shape)
       *(top-bottom)
       +bottom
     )
with {
  // https://www.desmos.com/calculator/pn4myus6x4
  shaper(s,x) = (x-x*s)/(s-x*2*s+1);
};

parMean =
  sequentialOperatorParOut(N,+)
  : par(i, N+1, _/(1<<i)) ;

parRMS =
  pow(2):sequentialOperatorParOut(N,+)
  : par(i, N+1, _/(1<<i):sqrt) ;

sequentialOperatorParOut(N,op) =
  seq(i, N, operator(i))
with {
  operator(i) = si.bus(i), (_<: _ , op(_,_@(pow2(i))));
  pow2(i) = 1<<i;
  // same as:
  // pow2(i) = int(pow(2,i));
  // but in the block diagram, it will be displayed as a number, instead of a formula
};


//*****************************************************************************
//                 GUI
//*****************************************************************************

threshold = hslider("[03]threshold", 0, -70, 0, 0.01) : si.smoo ;
knee = hslider("knee", 0, 0, 30, 0.001) : si.smoo;
// threshold = hslider("threshold", 0, -70.0, 0.0, 0.01) : si.smoo: ba.db2linear ;
// knee = hslider("knee", 0.5, 0.0, 1.0, 0.001) : si.smoo;
pre_gain = hslider("pre-gain", 0, -30.0, 30.0, 0.01) : si.smoo: ba.db2linear ;
// pre_gain = hslider("pre-gain", ba.linear2db(ma.PI/2), -30.0, 30.0, 0.01) : ba.db2linear : si.smoo;


//*****************************************************************************
//                 CONSTANTS
//*****************************************************************************
// number of parallel means
// each mean doubles in size
// so the biggest mean is 2^N
// N=15 is 32768 samples, so 682 ms at 48k
N = 8;

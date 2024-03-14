
declare name "DJcomp";
declare version "0.1";
declare author "Bart Brouns";
declare license "AGPL-3.0-only";
declare copyright "2024, Bart Brouns";

import("stdfaust.lib");

// TODO:
// Auto Makeup: compensates for the drop in output level caused by gain reduction. It
// determines the theoretical compensation at a given threshold / ratio setting – for instance
// -20dB and 4:1 gives 5dB – and applies half of that value (2.5dB) to achieve the same perceived
// loudness.�
//  Also use attack + release time?
//
// Fast release after a transient,
// Slow release after longer periods of gain reduction.
// Aditionally adjust the release time depending on the current amount of gain reduction

top = hslider("top", 1, 0.001, 100, 0.001);
bottom = hslider("bottom", 1, 0.001, 100, 0.001);
speedFactor = hslider("speedF", 0, 0, 1, 0.001);

process =
  DJcomp;

linLogInf =
  (it.interpolate_linear(speedFactor,bottom,bottom/ma.EPSILON)
   : hbargraph("[0]linear", 0, 100))
,
  (interpolate_logarithmic(speedFactor,bottom,
                           bottom/ma.EPSILON
                          )
   : hbargraph("[1]logarithmic", 0, 100)
  )
,
  (fade_to_inf(speedFactor,bottom)
   : hbargraph("[2]to inf", 0, 100))
;


fade_to_inf(dv,v0) =
  v0/max(1-dv,ma.EPSILON);

interpolate_logarithmic(dv,v0,v1) =
  pow(((v1/v0)),dv)*v0
;
DJcomp =
  compressor_N_chan(strength,thresh,attack,release,knee,1,2);

compressor_N_chan(strength,thresh,att,rel,knee,link,N) =
  par(i, N, _*inputGain)
  <: si.bus(N*2)
  : (
  si.bus(N)
  // par(i, N, _)
  ,(compression_gain_N_chan(strength,thresh,att,rel,knee,link,N)
    <: si.bus(N*2))
)
    // <: si.bus(N*2)
  :
  (
    ro.interleave(N,2)
    : par(i,N, *)
  )
, si.bus(N)
;

compression_gain_N_chan(strength,thresh,att,rel,knee,link,1) =
  abs:ba.linear2db
  : compression_gain_mono(strength,thresh,att,rel,knee);

compression_gain_N_chan(strength,thresh,att,rel,knee,link,N) =
  par(i, N, abs:ba.linear2db)
  <: (si.bus(N),(ba.parallelMax(N) <: si.bus(N)))
  : ro.interleave(N,2)
  : par(i,N,(it.interpolate_linear(link))
       )
  : par(i,N,compression_gain_mono(strength,thresh,att,rel,knee)) ;

compression_gain_mono(strength,thresh,att,rel,knee) =
  loop~(_,_)
       : (_,!)
       : ba.db2linear
with {
  loop(prevGain,prevRef) =
    gain,ref
  with {
  gain =
    gain_computer(strength,thresh,knee)
    : smootherARorder(maxOrder, orderRel,orderAtt, adaptiveRel, att)
    : hbargraph("GR[unit:dB]", -24, 0);
  adaptiveRel =
    fade_to_inf(1-dv,rel) ;

  ref =
    (prevGain+transitionRange)
    : min(0)
      // : smootherOrder(maxOrder,refOrder,refRel,0)
    : smootherOrder(1,1,refRel,0)
    : hbargraph("ref[unit:dB]", -24, 0)
  ;
  refRel =
    interpolate_logarithmic(dv, slowRelease,slowRelease/ma.EPSILON) ;
  dv= (fastGR
       :min(0)
        / (transitionRange
          )
       *-1):min(1)
            // : hbargraph("dv", 0, 1)
  ;
  fastGR =
    (prevGain-prevRef):min(0)
                       // :hbargraph("fast GR[unit:dB]", -24, 0)
  ;
};
};


gain_computer(strength,thresh,knee,level) =
  select3((level>(thresh-(knee/2)))+(level>(thresh+(knee/2))),
          0,
          ((level-thresh+(knee/2)) : pow(2)/(2*max(ma.EPSILON,knee))),
          (level-thresh))
  : max(0)*-strength;


///////////////////////////////////////////////////////////////////////////////
//                               smoothers                                   //
///////////////////////////////////////////////////////////////////////////////

// fixed order
smoother(order, att, rel, xx) =
  smootherOrder(order,order, att, rel, xx);

smootherOrder(maxOrder,order, att, rel, xx) =
  smootherARorder(maxOrder,order, order, att, rel, xx);

smootherARorder(maxOrder,orderAtt, orderRel, att, rel, xx) =
  xx : seq(i, maxOrder, loop(i) ~ _)
with {
  loop(i,fb, x) = coeff(i) * fb + (1.0 - coeff(i)) * x
  with {
  cutoffCorrection(order) = 1.0 / sqrt(pow(2.0, 1.0 / order) - 1.0);
  coeff(i) =
    ba.if(x > fb, attCoeff(i), relCoeff(i) );
  attCoeff(i) =
    exp(-TWOPIT * cutoffCorrection(orderAtt) / max(ma.EPSILON, att))
    * (i<orderAtt);
  relCoeff(i) =
    exp(-TWOPIT * cutoffCorrection(orderRel) / max(ma.EPSILON, rel))
    * (i<orderRel);
  TWOPIT = 2 * ma.PI * ma.T;
};
};

///////////////////////////////////////////////////////////////////////////////
//                                    GUI                                   //
///////////////////////////////////////////////////////////////////////////////

inputGain = hslider("[01]input gain[unit:dB]", 0, -24, 24, 0.1):ba.db2linear:si.smoo;
strength = hslider("[02]strength[unit:%]", 100, 0, 100, 1) * 0.01;
thresh = hslider("[03]thresh[unit:dB]",-1,-30,0,0.1);
attack = hslider("[04]attack[unit:ms] [scale:log]",9, 1000/48000, maxAttack*1000,0.1)*0.001;
orderAtt =
  4;
// hslider("[05]attack order", 4, 1, maxOrder, 1);
release = hslider("[06]fast release[unit:ms] [scale:log]",60,0.1,maxRelease*1000,1)*0.001;
transitionRange = hslider("[07]release transition range[unit:dB]",9,0,30,0.1);
slowRelease = hslider("[08]slow release[unit:ms] [scale:log]",4000,50,10000,50)*0.001;
orderRel =
  4;
// hslider("[09]release order", 4, 1, maxOrder, 1);
knee = hslider("[10]knee[unit:dB]",1,0,72,0.1);

maxOrder = 4;
refOrder =
  hslider("[11]ref release order", 1, 1, maxOrder, 1);
// 100 ms
maxAttack = 0.1;
// 2 sec
maxRelease = 2;


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
    // <:  (ba.db2linear,_)
    // : autoSmoother
    // : ba.db2linear
    : smootherARorder(maxOrder, orderRelR,orderAttR, adaptiveRel, att)
      // : ba.linear2db
    : hbargraph("GR[unit:dB]", -24, 0);
  adaptiveRel =
    fade_to_inf(1-dv,rel) ;

  ref =
    (prevGain+slowKnee)
    : min(0)
      // : ba.bypass1(linLog,ba.db2linear)
    : smootherOrder(maxOrder,refOrder,refRel,0)
      // : ba.bypass1(linLog,ba.linear2db)
    : hbargraph("ref[unit:dB]", -24, 0)
  ;
  linLog = checkbox("linLog");
  refRel =
    interpolate_logarithmic(dv, relR,relR/ma.EPSILON) ;
  dv= (fastGR
       :min(0)
        / (slowKnee
          )
       *-1):min(1)
      : hbargraph("dv", 0, 1)
  ;
  fastGR =
    (prevGain-prevRef):min(0)
                       // :hbargraph("fast GR[unit:dB]", -24, 0)
  ;
  autoSmoother(lin,db) =
    lin
    : smoother(1,autoRelease(db),0)
    : smoother(4,rel,att)
  ;

  autoRelease(GR) =
    rel *
    ( ((GR+dead)
       :min(0)
        *-1
       :smootherARorder(maxOrder,orderAttR, orderRelR, attR, relR)
        * factor
      ):hbargraph("usf", 0, 10)
    );
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

inputGain = hslider("[01]input gain", 0, -24, 24, 0.1):ba.db2linear:si.smoo;
strength = hslider("[02]strength", 100, 0, 100, 1) * 0.01;
thresh = hslider("[03]thresh",-1,-30,0,0.1);
attack = hslider("[04]attack[unit:ms] [scale:log]",9, 1000/48000, maxAttack*1000,0.1)*0.001;
release = hslider("[06]release[unit:ms] [scale:log]",60,0.1,maxRelease*1000,1)*0.001;
knee = hslider("[09]knee",1,0,72,0.1);

maxOrder = 8;
attR = hslider("[04]slow attack[unit:ms] [scale:log]",700, 10, 3000,10)*0.001;
orderAttR =
  // 1;
  hslider("[05]slow attack order", 4, 1, maxOrder, 1);
relR = hslider("[06]ref release[unit:s] [scale:log]",4,1,50,1);
orderRelR =
  // 1;
  hslider("[07]slow release order", 4, 1, maxOrder, 1);
factor = hslider("[08]factor", 0.1, 0, 10, 0.1);
slowKnee = hslider("[09]slow knee",12,0,72,0.1);
refOrder =
  // 1;
  hslider("[10]ref release order", 1, 1, maxOrder, 1);
dead = hslider("dead", 0, 0, 24, 1);
// 100 ms
maxAttack = 0.1;
// 2 sec
maxRelease = 2;

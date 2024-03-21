
declare name "DJcomp";
declare version "0.1";
declare author "Bart Brouns";
declare license "AGPL-3.0-only";
declare copyright "2024, Bart Brouns";

import("stdfaust.lib");


enableGRout = 0;
// enable gain reduction outputs

process =
  DJcomp;

DJcomp =
  compressor_N_chan(strength,thresh,attack,fastRelease,knee,1,2) ;

compressor_N_chan(strength,thresh,att,rel,knee,link,N) =
  par(i, N, _*inputGain)
  <: si.bus(N*2)
  : (
  si.bus(N)
  // par(i, N, _)
 ,(compression_gain_N_chan(strength,thresh,att,rel,knee,link,N)
   <: si.bus(N*2)
  )
)
  : (
    ro.interleave(N,2)
    : par(i,N, *)
  )
, si.bus(N)
  : postProc(N,enableGRout)
;

compression_gain_N_chan(strength,thresh,att,rel,knee,link,1) =
  abs:ba.linear2db
  : compression_gain_mono(strength,thresh,att,rel,knee);

compression_gain_N_chan(strength,thresh,att,rel,knee,1,N) =
  par(i, N, abs)
  : ba.parallelMax(N)
  : ba.linear2db
  : compression_gain_mono(strength,thresh,att,rel,knee);

compression_gain_N_chan(strength,thresh,att,rel,knee,link,N) =
  par(i, N, abs:ba.linear2db)
  <: (si.bus(N),(ba.parallelMax(N) <: si.bus(N)))
  : ro.interleave(N,2)
  : par(i,N,(it.interpolate_linear(link))
       )
  : par(i,N,compression_gain_mono(strength,thresh,att,rel,knee)) ;


compression_gain_mono(strength,thresh,att,rel,knee,level) =
  loop~(_,_)
       : (_,!)
       : ba.db2linear
with {
  loop(prevGain,prevRef) =
    gain,ref
  with {
  gain =
    (  gain_computer(1,thresh,knee,level)
       : ba.db2linear
       : smootherARorder(maxOrder, orderRel,orderAtt, adaptiveRel, att)
         // , ((level-limThres):max(0)*-1: ba.db2linear)
     , (gain_computer(1,limThres,limKnee,level): ba.db2linear)
    ):min
    : smootherARorder(maxOrder, orderRelLim,4, releaseLim, 0)
    : ba.linear2db
      * strength
    : hbargraph("GR[unit:dB]", -24, 0);

  adaptiveRel =
    fade_to_inf(
      // shaper(adaptShape,
      1-dv
      // )
     ,rel) ;

  ref =
    (prevGain+transitionRange)
    : min(0)
    : ba.db2linear
      // : smootherOrder(maxOrder,refOrder,refRel,0)
    : smootherOrder(1,1,refRel,0)
    : ba.linear2db
      // : hbargraph("ref[unit:dB]", -24, 0)
  ;
  refRel =
    interpolate_logarithmic(
      // shaper(refShape,
      dv
      // )
    , slowRelease,slowRelease/ma.EPSILON) ;
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

postProc(N,0) = si.bus(N),par(i, N, !);
postProc(N,1) = si.bus(N*2);

///////////////////////////////////////////////////////////////////////////////
//                               smoothers                                   //
///////////////////////////////////////////////////////////////////////////////

// smoother adapted from Dario Sanfilippo
// https://github.com/dariosanfilippo/limiterStereo/blob/da1c38cc393f08b5dd79e56ffd4e6256af07a708/limiterStereo.dsp#L90-L101
//
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
//                              Utilities                                    //
///////////////////////////////////////////////////////////////////////////////

fade_to_inf(dv,v0) =
  v0/max(1-dv,ma.EPSILON);

interpolate_logarithmic(dv,v0,v1) =
  pow(((v1/v0)),dv)*v0 ;

// https://www.desmos.com/calculator/pn4myus6x4
shaper(s,x) = (x-x*s)/(s-x*2*s+1);
///////////////////////////////////////////////////////////////////////////////
//                                    GUI                                   //
///////////////////////////////////////////////////////////////////////////////

oneKnob = hslider("one knob", 0, 0, 1, 0.01):si.smoo;

inputGainSlider = hslider("[01]input gain[unit:dB]", 0, -24, 24, 0.1):si.smoo;
inputGain = (inputGainSlider + it.remap(0, 0.5, -3, 6,oneKnob:min(0.5))):ba.db2linear;
strength = it.remap(0, 0.5, 0, 1,oneKnob:min(0.5));
// hslider("[02]strength[unit:%]", 100, 0, 100, 1) * 0.01;

thresh =
  // limThres + hslider("[03]thresh offset[unit:dB]",0,-12,12,0.1);
  limThres + it.remap(0.5, 1, 0, 6,oneKnob:max(0.5));
attack =
  0.009;
// hslider("[04]attack[unit:ms] [scale:log]",9, 1000/48000, maxAttack*1000,0.1)*0.001;
orderAtt =
  4;
// hslider("[05]attack order", 4, 1, maxOrder, 1);
fastRelease =
  0.42;
// hslider("[06]fast release[unit:ms] [scale:log]",420,0.1,maxRelease*1000,1)*0.001;
transitionRange =
  9;
// hslider("[07]release transition range[unit:dB]",9,0,30,0.1);
slowRelease =
  // hslider("[08]slow release[unit:ms] [scale:log]",2000,50,10000,50)*0.001;
  it.remap(0.5, 1, 2, 0.5,oneKnob:max(0.5));
orderRel =
  1;
// hslider("[09]release order", 1, 1, maxOrder, 1);
knee =
  // hslider("[10]knee[unit:dB]",1,0,72,0.1);
  it.remap(0.5, 1, 9, 15,oneKnob:max(0.5));
maxOrder = 4;
refOrder =
  hslider("[11]ref release order", 1, 1, maxOrder, 1);
adaptShape = hslider("[12]adapt shape", 0, -1, 1, 0.001);
refShape = hslider("[13]ref shape", 0, -1, 1, 0.001);

// give it some headroom.
// when built into a mixer, there should be a limiter on the master.
// if that is done, this can come back up to -1dB
limThres = hslider("[14]thresh[unit:dB]",-6,-30,30,0.1);
releaseLim =
  // hslider("[15]release limiter[unit:ms] [scale:log]",60,5,500,1)*0.001;
  it.remap(0.5, 1, 0.06, 0.12,oneKnob:max(0.5));

orderRelLim =
  4;
// hslider("[16]lim release order", 4, 1, maxOrder, 1);
limKnee =
  // hslider("[17]lim knee[unit:dB]",1,0,72,0.1);
  it.remap(0.5, 1, 12, 0,oneKnob:max(0.5));
// 100 ms
maxAttack = 0.1;
// 2 sec
maxRelease = 2;

// https://www.desmos.com/calculator/qcjwfaaqc5

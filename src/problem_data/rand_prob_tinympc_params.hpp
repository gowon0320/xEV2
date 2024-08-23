#include "types.hpp"
#include "Arduino.h"

#pragma once

const PROGMEM tinytype Adyn_data[2*2] = {
-1.1500000, -0.9980000,
0.4000000, -1.2920000
};

const PROGMEM tinytype Bdyn_data[2*1] = {
0.0000000,
0.0002500
};


const PROGMEM tinytype Q_data[2] = {100000.0000000, 10000000.0000000};

const PROGMEM tinytype Qf_data[2] = {100000.0000000, 10000000.0000000};

const PROGMEM tinytype R_data[1] = {0.0100000};

const PROGMEM tinytype umin[4] = {
  -30000.0,	
  -30000.0,	
  -30000.0,	
  -30000.0 
};

const PROGMEM tinytype umax[4] = {
  30000.0,	
  30000.0,	
  30000.0,	
  30000.0	
};

const PROGMEM tinytype xmin[4] = {
  -10000.0,	
  -10000.0,	
  -10000.0,	
  -10000.0,	
};

const PROGMEM tinytype xmax[4] = {
  10000.0,	
  10000.0,	
  10000.0,	
  10000.0,	
};

const PROGMEM tinytype rho_value = 3;

const PROGMEM tinytype Kinf_data[1*2] = {
147.1, -4679
};

const PROGMEM tinytype Pinf_data[2*2] = {
9.934e+07, 8.033e+07,
8.033e+07, 2.641e+08
};

const PROGMEM tinytype Quu_inv_data[1*1] = {
0.04442
};

const PROGMEM tinytype AmBKt_data[2*2] = {
-1.15,  0.3632,
-0.998, -0.1223
};

const PROGMEM tinytype coeff_d2p_data[4*4] = {
  -8.146955332577477e-14,	9.96876192704832e-14,	2.1733493910797153e-13,	-3.061439990403869e-14,	
  -6.510764150036152e-14,	7.967584925161475e-14,	1.7257723028407668e-13,	-2.492624162631074e-14,	
  -6.855627177060342e-14,	8.348183255790786e-14,	1.8213555663670888e-13,	-2.588294162331195e-14,	
  1.8038522064944829e-13,	-2.152167333235866e-13,	-4.731701142013378e-13,	6.641171987342709e-14,	
};


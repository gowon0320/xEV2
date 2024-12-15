#include "types.hpp"
#include "Arduino.h"

#pragma once

const PROGMEM tinytype Adyn_data[2 * 2] = {
    0.9809, -0.0096,
    0.0261, 0.9783};

const PROGMEM tinytype Bdyn_data[2 * 1] = {
    -1.204e-08,
    2.473e-06};

const PROGMEM tinytype Q_data[2] = {4.5e+04, 4.5e+06};

const PROGMEM tinytype Qf_data[2] = {4.5e+04, 4.5e+06};

const PROGMEM tinytype R_data[1] = {5.01};

const PROGMEM tinytype umin[1] = {
    -800.0};

const PROGMEM tinytype umax[1] = {
    800.0};

const PROGMEM tinytype xmin[2] = {
    -0.1745,
    -0.5};

const PROGMEM tinytype xmax[2] = {
    0.1745,
    0.5};

const PROGMEM tinytype rho_value = 5;

const PROGMEM tinytype Kinf_data[1 * 2] = {
    5.263e-08, 2.415e-06};

const PROGMEM tinytype Pinf_data[2 * 2] = {
    4.501e+04, 0.08058,
    0.08058, 4.5e+06};

const PROGMEM tinytype Quu_inv_data[1 * 1] = {
    0.1996};

const PROGMEM tinytype AmBKt_data[2 * 2] = {
    0.9809, 0.0261,
    -0.0096, 0.9783};

const PROGMEM tinytype coeff_d2p_data[2 * 1] = {
    -0.2899,
    -10.89};

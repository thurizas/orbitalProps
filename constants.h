#ifndef _constants_h_
#define _constants_h_

#include "common.h"

#include <cstdint>

const double PI = 3.14159265359;
const double DEG2RAD = PI / 180.0f;       // 1 degree = 0.01745329252 radians
const double EPSILON = 1E-6;              // sensitivity for numerical methods

const double  AU = 1.495978707E8;         // 1AU is 1.495 x 10^8 km
const double  AU2M = 1000 * AU;
const double  M2AU = 1 / AU2M;            // conversion from meters to AU

const double SIDEREALDAY = 23.934469444;  // 1 Sidereal day is 23h 56m 4.09s
const double DAY2SEC = 8.6164089984E4;    // 86,164 seconds per sidereal day

const double G = 6.67430E-11;             // G is 6.67430 x 10-11 Nm^2/Kg^2
const double K = 0.01720209895;           // Gausian gravitational constant (units: L^{1.5}T^{-1}M^{-0.5}_

const double RSOL = 6.96340E5;            // radius of the sun in kilometers
const double MSOL = 1.988400E30;          // mass of sun in kilograms

const double REARTH = 6371;               // radius of the earth in kilometers
const double MEARTH = 5.9722E24;          // mass of earth in kilograms


static const char epoch[] = "J2000.0";    //  reference epoch from NASA

// see: https://ssd.jpl.nasa.gov/planets/approx_pos.html, table 1
//      earth data from: https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
// data for epoch J2000, valid for [1800AD, 2050AD]
//                          n          a            e            L               i            w_bar               W           m         o
orbitalPropT planets[] = { {"mercury",  0.38709927f, 0.20563593f,  252.25032350f, 7.00497902f,  77.45779628f,  48.33076593f,  0.0553f,  0.03f},
                           {"venus"  ,  0.72333566f, 0.00677672f,  181.97909950f, 3.39467605f, 131.60246718f,  76.67984255f,  0.815f ,  2.64f},
                           {"earth"  ,  1.00000261f, 0.01671123f,  100.46457166f,        0.0f,    102.94719f,    -11.26064f,   1.00f , 23.44f},
                           {"mars"   ,  1.52371034f, 0.09339410f,   -4.55343205f, 1.84969142f, -23.94362959f,  49.55953891f,  0.1075f, 25.19f},
                           {"jupiter",  5.20288700f, 0.04838624f,   34.39644051f, 1.30439695f,  14.72847983f, 100.47390909f, 317.9f   ,  3.13f},
                           {"saturn" ,  9.53667594f, 0.05386179f,   49.95424423f, 2.48599187f,  92.59887831f, 113.66242448f, 95.2f   , 26.73f},
                           {"uranus" , 19.18916464f, 0.04725744f,  313.23810451f, 0.77263783f, 170.95427630f,  74.01692503f, 14.6f   , 82.23f},
                           {"neptune", 30.06992276f, 0.008959048f, -55.12002969f, 1.77004347f,  44.96476227f, 131.78422574f, 17.2f   , 28.32f},
};
// note argument of periapsis, w, given by w = w_bar-W
//      mean anomoly, M, given by M = L - w_bar 

static const char* months[] = { "January", "Feburary", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December" };


#endif

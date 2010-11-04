
%module libnova

%{
#include "libnova.h"
%}


/* libnova public header files */

%include aberration.h
%include angular_separation.h
%include apparent_position.h
%include asteroid.h
%include comet.h
%include dynamical_time.h
%include earth.h
%include elliptic_motion.h
%include hyperbolic_motion.h
%include julian_day.h
%include jupiter.h
%include ln_types.h
%include lunar.h
%include mars.h
%include mercury.h
%include neptune.h
%include nutation.h
%include parabolic_motion.h
%include pluto.h
%include precession.h
%include proper_motion.h
%include refraction.h
%include rise_set.h
%include saturn.h
%include sidereal_time.h
%include solar.h
%include transform.h
%include uranus.h
%include utility.h
%include venus.h
%include vsop87.h



/* provide SWIG pointer library for time_t pointers */

typedef long int time_t;

%include cpointer.i
%pointer_functions(time_t, time_t) 









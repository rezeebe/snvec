#===============================================================#
#
# ReadMe.txt
#
# snvec: Pre-computed Precession-Tilt solutions and C code.
# 
# *** snvec comes with ABSOLUTELY NO WARRANTY ***		 
# *** Use at your own risk. DO NOT DISTRIBUTE  ***
#
# When using snvec, cite as:
#
# A deep-time dating tool for paleo-applications utilizing obliquity 
# and precession cycles: The role of dynamical ellipticity and tidal 
# dissipation, Richard E. Zeebe and Lucas J. Lourens, Paleoceanography 
# and Paleoclimatology, 2022.
#
# Richard E. Zeebe
# School of Ocean and Earth 
# Science and Technology
# University of Hawaii at Manoa
# 1000 Pope Road, MSB 504
# Honolulu, HI 96822, USA
# email: zeebe@soest.hawaii.edu
#
#===============================================================#   

NOTE: The provided, pre-computed Precession-Tilt (PT) solutions can 
be used without installing/using the C code package. The code may be 
used to generate solutions with parameter values not included in 
the pre-computed PT solutions.

(1) FILES/FOLDERS

ReadMe.txt This file.
snvec-$v.tar.gz  snvec C code package.
OS/ZB18a         Orbital Solution ZB18a (Zeebe and Lourens, Science, 2019)
ZB18a/asc        Precession-Tilt solutions for ZB18a. ASCII  format
ZB18a/bin        Precession-Tilt solutions for ZB18a. Binary format
ZB18a/asc.tar.gz Precession-Tilt solutions for ZB18a. ASCII  format (all zipped)
ZB18a/bin.tar.gz Precession-Tilt solutions for ZB18a. Binary format (all zipped)

ASCII/BIN files. PT.Dex.xxxxTdy.yyyy.ext is the PT solution for 
ZB18a using dynamical ellipticity Ed = x.xxxx and tidal dissipation 
Td = y.yyyy.

ASCII file content (columns):
time                   [kyr] in the past
obliquity              [radian]
precession             [radian] from ECLIPJ2000
climatic precession    [] e*sin(\varpi*)

where \varpi* is the longitude of perihelion relative to the moving equinox.

CAUTION: The output time intervals are not equally spaced. If required for 
e.g.,spectral analysis, the user needs to interpolate to equally spaced
time intervals.

BIN file content (blocks): {M,time[days],obliquity,precession,climatic precession}
where M = length(time). See below for a MATLAB example to read the binary output.

(2) Code: Installation & Run

On most linux/unix/mac systems, all that should be required is a 
C compiler such as gcc (free). Possible free options for Windows users 
include installing Dev-C++ or making use of projects such as Cygwin 
or MinGW.

On linux/unix/mac: unzip snvec-$v.tar.gz. Do not change
the relative location of the fun folder. The program snvec requires 
the Orbital Solution (OS) OS/ZB18a/ems-plan3.dat. The assumed
default location is /dat/PrecTilt/OS/ZB18a/. Either create and 
save ems-plan3.dat in that folder or change the arguments of the
executable (see below).

compile (example: gcc):
$ gcc -o snvec.x snvec-3.7.5.c -lm

run (example):
$ ./snvec.x -1.e3 1.0 0.0 /dat/PrecTilt/OS/ZB18a ems-plan3.dat

argument list:              Defaults
 [1] tend (kyr, negative)   -1.e3
 [2] Ed                      1.0000
 [3] Td                      0.0000
 [4] dir  OrbitSoln          /dat/PrecTilt/OS/ZB18a
 [5] file OrbitSoln          ems-plan3.dat

running the following command uses the above defaults:
$ ./snvec.x

For the content of the output files out.dat and out.bin, see
ASCII/BIN file content above.

NOTE: By default, effects of tidal dissipation on obliquity
(secular trend) were not included. However, the C code includes 
this option, following Quinn et al. (1991). To activate, change
the macro UEPSDOT to EPSDOT and recompile.

(3) MATLAB example: read bin output

ky2d = 365250.;
fid  = fopen('out.bin','r');
M    = fread(fid,1,'int');
t    = fread(fid,M,'double')/ky2d; % days => kyr
epl  = fread(fid,M,'double');
phi  = fread(fid,M,'double');
cp   = fread(fid,M,'double');
fclose(fid);


/* Vocal Tract Simulation Code (C) 2025 Jyothiraditya Nellakra
 * See end of file for software licence. */

// Constants & variables defined by <main_male.cc> and <main_female.cc>.
// See my paper for explanations of these constants.
// Many values taken from Stevens [1].

extern const long double M1, M2; // Masses.
extern const long double k1, kC; // Spring constants.
extern const long double k2_start, k2_stop, k2_step; // Progressively iterated.
extern const long double t_start, t_stop, t_step; // Time range simulated.
extern const long double x1_0, x2_0, y2_0; // Initial positions.
extern const long double Px1_0, Px2_0, Py2_0; // Initial momenta.
extern const long double x1_eq, x2_eq, y2_eq, l_eq; // Equilibrium positions.

/* Defining Constants & Declaring Variables */

// Masses.
const long double M1 = 0.1; // g/cm
const long double M2 = 0.02; // g/cm

// Spring constants.
const long double k1 = 33000; // dyne/cm^2
const long double kC = 20000; // dyne/cm^2

// We will step through different values of k2 in the video so that we can find
// the frame with the most accurate value of k2.

const long double k2_start = 10000;  // dyne/cm^2
const long double k2_stop = 30000; // dyne/cm^2
const long double k2_step = (k2_stop - k2_start) / (10.0 * 24.0); // 10 seconds, 24 fps.

// Time range simulated.
const long double t_start = 0; // s
const long double t_stop = 0.05; // s
const long double t_step = 1.0 / 100000000.0; // 0.1 GHz

// Initial position coordinates
const long double x1_0 = 0.13; // cm
const long double x2_0 = 0.23; // cm
const long double y2_0 = 0.38; // cm

// Initial momentum coordinates
const long double Px1_0 = 0; // cm/s
const long double Px2_0 = 0; // cm/s
const long double Py2_0 = 0; // cm/s

// Equilibrium position coordinates for C_1 and C_2.
const long double x1_eq = 0.19; // cm
const long double x2_eq = 0.12; // cm
const long double y2_eq = 0.36; // cm
const long double l_eq = 0.37; // cm; equilibrium length for C_C.

/* References */

// 1: K. N. Stevens, Acoustic phonetics (MIT Press, July 24, 2000), 628 pp.

/* License for this software. */

/* This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program. If not, see <https://www.gnu.org/licenses/>. */

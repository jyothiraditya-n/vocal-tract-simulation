#! /usr/bin/python3

# Vocal Tract Simulation Code (C) 2025 Jyothiraditya Nellakra
# See end of file for software licence.

# Imports to get the code to work.
from scipy.fft import fft, fftfreq
from math import sin, sqrt, tan, pi
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys, os

# We need at least one argument on the command line: the CSV file we are going
# convert into a graph.
if len(sys.argv) != 2:
    print(f"usage: {sys.argv[0]} <Path>\n")
    exit(1)

# Create the output file name from the input.
input_file = sys.argv[1]
output_file = os.path.splitext(input_file)[0] + ".png"

# Read the input CSV.
data_set = pd.read_csv(input_file, sep='\t', header = None)
data_frame = pd.DataFrame(data_set)

# Get the values we want.
values = np.array(data_frame.values)

# Split the simulation data into the relevant variables.
k2 = values[0, 0]
ts = values[:, 1]
x1s = values[:, 2]
x2s = values[:, 3]
y2s = values[:, 4]

# Code for creating a graph of the Fourier transform of x_2.

N = len(ts) # Number of entries
T = 0.05 / N # Time elapsed between entries, must be copied from the C code
    # whenever the value in the C code changes.

fs = fftfreq(N, T)[10:1000] # Frequency values; 10 Hz--1 kHz limit.
x2fs = np.log(np.abs(fft(x2s)))[10:1000] # ln(|amplitude|) for each f.
# x1fs = fft(x1s)[10:1000]
# y2fs = fft(y2s)[10:1000]

# Remove all values that are less than 0.01 by setting them to 0.01.
x2fs = np.clip(x2fs, 0.01, np.max(x2fs))

# The following is code for simulating the upper vocal tract as a linear
# acoustic filter.

# everything in CGS units. See my paper for explanations of these constants.
c = 34000 # speed of sound
mu = 0.00019 # viscosity of air
rho = 0.001225 # density of air
lambda_var = 0.000055
c_p = 0.24

def my_tan(x):
    # Remove asymptotes to produce cleaner graphs.
    return max(min(tan(x), 100), -100)

def my_cot(x):
    return 1 / my_tan(x)

# Cross-sectional areas of the vocal tract at 0.5 cm intervals from the glottis
# to the mouth, pulled from Table 2.33-1(A) from Fant [1]. These are the
# standard vowels of Russian.

l_div = 0.5 # cm distance per sample point along glottal path.

vowel_areas = { # cross-sectional areas.
    'ɒ': [
        5, 5, 5, 5, 6.5, 8, 8, 8, 8, 8, 8, 8, 8, 6.5, 5, 4, 3.2, 1.6, 2.6, 2.6,
        2, 1.6, 1.3, 1, 0.65, 0.65, 0.65, 1, 1.6, 2.6, 4, 1, 1.3, 1.6, 2.6
    ],
    'o': [
        3.2, 3.2, 3.2, 3.2, 6.5, 13, 13, 16, 13, 10.5, 10.5, 8, 8, 6.5, 6.5, 5,
        5, 4, 3.2, 2, 1.6, 2.6, 1.3, 0.65, 0.65, 1, 1, 1.3, 1.6, 2, 3.2, 4, 5,
        5, 1.3, 1.3, 1.6, 2.6
    ],
    'u': [
        0.65, 0.65, 0.32, 0.32, 2, 5, 10.5, 13, 13, 16, 13, 10.5, 10.5, 8, 8,
        6.5, 6.5, 5, 5, 4, 3.2, 2, 1.6, 2.6, 1.3, 0.65, 0.65, 1, 1, 1.3, 1.6,
        2, 3.2, 4, 5, 5, 1.3, 1.3, 1.6, 2.6
    ],
    'ɨ': [
        6.5, 6.5, 2, 6.5, 8, 8, 8, 5, 3.2, 2.6, 2, 2, 1.6, 1.3, 1, 1, 1.3, 1.6,
        2.6, 2, 4, 5, 6.5, 6.5, 8, 10.5, 10.5, 10.5, 10.5, 10.5, 13, 13, 10.5,
        10.5, 3.2, 3.2, 3.2, 3.2, 3.2
    ],
    'i': [
        4, 4, 3.2, 1.6, 1.3, 1, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 1.3,
        2.6, 4, 6.5, 8, 8, 10.5, 10.5, 10.5, 10.5, 10.5, 10.5, 10.5, 10.5,
        10.5, 8, 8, 2, 2, 2.6, 3.2
    ],
    'e': [
        8, 8, 5, 5, 4, 2.6, 2, 2.6, 2.6, 3.2, 4, 4, 4, 5, 5, 6.5, 8, 6.5, 8,
        10.5, 10.5, 10.5, 10.5, 10.5, 8, 8, 6.5, 6.5, 6.5, 6.5, 1.3, 1.6, 2.0,
        2.6
    ]
}

# Helper function to generate the transfer function representing the upper
# vocal tract during the production of a given vowel.
def construct_T_f(vowel):
    # Get the appropriate list of cross-sectional areas.
    areas = vowel_areas[vowel]

    # We want to split the list of cross-sectional areas into sub sections
    # based on the points of constriction in the vocal tract.

    areas_collections = [] # collections of lists of cross-sectional areas.
    running = [areas[0]] # current list of cross-sectional areas.

    for area in areas[1:]:
        if area >= running[-1]:
            running.append(area)
        else:
            areas_collections.append(running)
            running = [area]

    areas_collections.append(running)

    # Our "sections" of the acoustic pipe model, where we want to know if we
    # are going to use a tan or a cot function (modelling an open or a closed
    # pipe, respectively), and we need to know the cross-sectional area A of
    # the pipe and the length l to compute the impedance of this section of
    # pipe. See my paper for more details.

    trig_area_lengths = []

    for i in areas_collections:
        area = i[0]
        length = 1

        for j in i[1:]:
            if j == area:
                length += 1
            else:
                break

        trig_area_lengths.append((my_cot, area, length * l_div))
        i = i[length:]

        if len(i) != 0:
            trig_area_lengths.append((my_tan, sum(i)/len(i), len(i) * l_div))

    # Return a transfer function built on this set of acoustic pipes. See my
    # paper for details on how we mathematically derive this.
    def abs_T(f):
        return abs(sum([ # Note: this is a Python list comprehension.
            1 / (
                (1 / (
                    (1j * rho * c * (1/A_i) * trig(2 * pi * f * l_i / c))
                    + (
                        l_i * (2 * pi * sqrt(A_i / pi) / (A_i ** 2))
                        * sqrt(pi * f * rho * mu)
                    )
                )) + (1 / (
                    (
                        rho * (c ** 2) / (0.4 * (2 * pi * sqrt(A_i / pi))))
                        * sqrt((c_p * rho) / (pi * lambda_var * f)
                    )
                ))
            )

            for trig, A_i, l_i in trig_area_lengths
        ]))

    return abs_T

# Construct the transfer functions for all of the vowels.
vowel_T_fs = {key: construct_T_f(key) for key in vowel_areas}
vowel_x2fs = {key: np.array(x2fs) for key in vowel_areas}

# Apply the transfer functions to get amplitude values for all frequencies.
for index, f in enumerate(fs):
    for vowel in vowel_x2fs:
        vowel_x2fs[vowel][index] *= vowel_T_fs[vowel](f)

for vowel in vowel_x2fs:
    vowel_x2fs[vowel] = 10000 * np.log(vowel_x2fs[vowel])
    vowel_x2fs[vowel] = np.clip(vowel_x2fs[vowel], 0, np.max(vowel_x2fs[vowel]))

# The following are different output functions for the kind of graphs we want
# to produce. Scroll down to the bottom to change which one is currently
# selected.

# Code for creating a graph of the movement of each coordinate over time.
def generate_vocal_fold_motion_graphs():
        plt.scatter(ts, x1s, c=['#880000'], s = 0.1, alpha = 1)
        plt.scatter(ts, x2s, c=['#008800'], s = 0.1, alpha = 1)
        plt.scatter(ts, y2s, c=['#000088'], s = 0.1, alpha = 1)
        plt.xlabel('$t$ / s')
        plt.ylabel('$x_1$ (red), $x_2$ (green), $y_2$ (blue) / cm')
        plt.title(f"Vocal Fold Simulation. $k_2 = {k2}$ dynes/cm$^2$")
        plt.savefig(output_file)

# Code for creating a graph of the Fourier transform of x_2.
def generate_vocal_fold_motion_spectrum():
        plt.stackplot(fs, 2.0/N * x2fs)
        plt.xlabel('$f$ / Hz')
        plt.ylabel('$\\ln\\left[\\left|\\hat{x}_2(f)\\right|\\right]$')
        plt.title(f"Vocal Fold Simulation. $k_2 = {k2}$ dynes/cm$^2$")
        plt.savefig(output_file)

# Code for generating a vowel spectrum based on the acoustic filtering provided
# by the supralaryngeal vocal tract.
def generate_vowel_spectrum(vowel):
    plt.step(fs, 2.0/N * vowel_x2fs[vowel])
    plt.xlabel('$f$ / Hz')
    plt.ylabel('$\\ln\\left[\\left|\\hat{x}_2(f)\\right||T(f)|\\right] * 10^4$')
    plt.title(f"Vocal Fold Simulation. Vowel /{vowel}/. "
              + f"$k_2 = {k2}$ dynes/cm$^2$")
    plt.savefig(output_file)

# The function we want to run.
generate_vowel_spectrum('u')

# References

# 1 G. Fant, Acoustic theory of speech production: with calculations based on
#   x-ray studies of russian articulations, D A C S R Series (De Gruyter,
#   Incorporated, 1971).

# Licence for this software.

# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.

# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.

# You should have received a copy of the GNU General Public License along with
# this program. If not, see <https://www.gnu.org/licenses/>.

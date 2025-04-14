/* Vocal Tract Simulation Code (C) 2025 Jyothiraditya Nellakra
 * See end of file for software licence. */

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <list>
#include <sstream>
#include <string>
#include <syncstream>
#include <thread>
#include <vector>

/* Defining Constants & Declaring Variables */

// Constants & variables defined by <main_male.cc> and <main_female.cc>.
// See my paper for explanations of these constants.

extern const long double M1, M2; // Masses.
extern const long double k1, kC; // Spring constants.
extern const long double k2_start, k2_stop, k2_step; // Progressively iterated.
extern const long double t_start, t_stop, t_step; // Time range simulated.
extern const long double x1_0, x2_0, y2_0; // Initial positions.
extern const long double Px1_0, Px2_0, Py2_0; // Initial momenta.
extern const long double x1_eq, x2_eq, y2_eq, l_eq; // Equilibrium positions.

// Variables used as we progress in the simulation.
long double k2; // dyne/cm^2
long double t; // s

// How many cycles to wait before writing down the computation values.
size_t log_freq = 100; // output frames per simulated frames.

long double x1, x2, y2; // position coordinates
long double Px1, Px2, Py2; // momenta

// All of the files that we will generate will go in the "output" folder.
std::string folder_name{"output"};
std::vector<std::string> csv_files;

/* Helper functions */

// Define the change over time functions that we found from the Hamiltonian.
long double dx1_dt() { return Px1 / M1; }
long double dx2_dt() { return Px2 / M2; }
long double dy2_dt() { return Py2 / M2; }

long double dPx1_dt() {
	return -(
		2 * k1 * (x1 - x1_eq)
		+ (2 * kC * (sqrt(pow(x1 - x2, 2) + pow(y2, 2)) - l_eq) * (x1 - x2))
		/ sqrt(pow(x1 - x2, 2) + pow(y2, 2))
	);
}

long double dPx2_dt() {
	return -(
		2 * k2 * (x2 - x2_eq)
		- (2 * kC * (sqrt(pow(x1 - x2, 2) + pow(y2, 2)) - l_eq) * (x1 - x2))
		/ sqrt(pow(x1 - x2, 2) + pow(y2, 2))
	);
}

long double dPy2_dt() {
	return -(
		2 * k2 * (y2 - y2_eq)
		+ (2 * kC * (sqrt(pow(x1 - x2, 2) + pow(y2, 2)) - l_eq) * (y2))
		/ sqrt(pow(x1 - x2, 2) + pow(y2, 2))
	);
}

std::string simulate(size_t frame) {
	std::ostringstream ss;
	ss << folder_name << "/" << std::setfill('0') << std::setw(5)
		<< frame << ".csv";

	std::string file_name{ss.str()};
	std::ofstream csv_file{file_name};

	x1 = x1_0; x2 = x2_0; y2 = y2_0;
	Px1 = Px1_0; Px2 = Px2_0; Py2 = Py2_0;

	csv_file << std::setprecision(3) << k2 << "\t" << t << "\t" << x1
		<< "\t" << x2 << "\t" << y2 << "\n";

	size_t log = 0;

	for(t = t_start; t < t_stop; t += t_step) {
		long double dx1 = dx1_dt() * t_step;
		long double dx2 = dx2_dt() * t_step;
		long double dy2 = dy2_dt() * t_step;

		long double dPx1 = dPx1_dt() * t_step;
		long double dPx2 = dPx2_dt() * t_step;
		long double dPy2 = dPy2_dt() * t_step;

		x1 += dx1; x2 += dx2; y2 += dy2;
		Px1 += dPx1; Px2 += dPx2; Py2 += dPy2;

		log++;
		if(log % log_freq) continue;

		csv_file << k2 << "\t" << t << "\t" << x1
		<< "\t" << x2 << "\t" << y2 << "\n";
	}

	return file_name;
}

void call_python(std::string file_name) {
	std::string prog{"./main.py "};
	prog.append(file_name);

	std::osyncstream{std::cout} << "+ " << prog << "\n";
	std::system(prog.c_str());
}

int main(int argc, char **argv) {
	std::ostringstream ss;
	for(int i = 1; i < argc; i++) ss << argv[i] << "_";

	if(argc > 1) folder_name = ss.str();
	std::filesystem::create_directory(folder_name);

	std::list<std::thread> threads;
	size_t frame = 0;

	for(k2 = k2_start; k2 < k2_stop; k2 += k2_step) {
		std::string file_name = simulate(frame);
		csv_files.push_back(file_name);

		if(threads.size() >= std::thread::hardware_concurrency()) {
			threads.begin() -> join();
			threads.pop_front();
		}

		threads.push_back(
			std::thread{call_python, *csv_files.rbegin()}
		);

		frame++;
	}

	for(auto& i: threads) i.join();

	std::string prog{"./main.sh "};
	prog.append(folder_name);

	std::cout << "+ " << prog << "\n";
	std::system(prog.c_str());
}

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

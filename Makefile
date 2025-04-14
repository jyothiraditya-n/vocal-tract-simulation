# Makefile Automated Build System

# Vocal Tract Simulation Code (C) 2025 Jyothiraditya Nellakra
# See end of file for software licence.

.DEFAULT_GOAL := main
.PHONY: main run clean

# Goals

main.o: main.cc
	g++ -c main.cc -O3 -o main.o --std=c++20

main_male.o: main_male.cc
	g++ -c main_male.cc -O3 -o main_male.o --std=c++20

main_female.o: main_female.cc
	g++ -c main_female.cc -O3 -o main_female.o --std=c++20

main_male: main_male.o main.o
	g++ main_male.o main.o -O3 -o main_male --std=c++20

main_female: main_female.o main.o
	g++ main_female.o main.o -O3 -o main_female --std=c++20

# Commands

main: main_male main_female

run:	main
	./main_male
	mv output/ male/
	mv output.mp4 male.mp4
	./main_female
	mv output/ female/
	mv output.mp4 female.mp4

clean:
	-rm output/ output.mp4 male/ male.mp4 female/ female.mp4
	-rm *.o main main_male main_female

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

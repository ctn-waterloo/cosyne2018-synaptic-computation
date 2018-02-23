#   This file is part of NEF Synaptic Computation
#   (c) Andreas Stöckel 2017, 2018
#
#   NEF Synaptic Computation is free software: you can redistribute it and/or
#   modify it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   NEF Synaptic Computation is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

LIBPATH ?= two_compartment_lif.so

all: $(LIBPATH)

$(LIBPATH): two_compartment_lif.cpp Makefile
	g++ -std=c++11 -Wall -O3 -c -fPIC two_compartment_lif.cpp -lm -o $(LIBPATH).o -fno-exceptions -fno-rtti
	g++ -std=c++11 -shared -s -Wl,--as-needed,-soname,two_compartment_lif.so -lm -o $(LIBPATH) $(LIBPATH).o

clean:
	rm -f $(LIBPATH).o $(LIBPATH)


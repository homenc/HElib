# Copyright (C) 2012,2013 IBM Corp.
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
# 
CC = g++
LD = g++
AR = ar
#
# CFLAGS = -g -O2 -Wfatal-errors -Wshadow -Wall -std=c++11 -I/usr/local/include -I.. -DUSE_ALT_CRT -DDEBUG_PRINTOUT -DUSE_ZZX_POLY
CFLAGS = -O2 -Wfatal-errors -Wshadow -Wall -std=c++11 -I/usr/local/include -I..
GMP=-lgmp 
LDLIBS = -L/usr/local/lib -lntl $(GMP) -lm

all: Test_AES_x

%.o: %.cpp
	$(CC) $(CFLAGS) -c $<

Test_AES_x: Test_AES.cpp simpleAES.o homAES.o ../fhe.a
	$(CC) $(CFLAGS) -o $@ $< simpleAES.o homAES.o ../fhe.a $(LDLIBS)

clean:
	rm -f *.o *_x *_x.exe *.a core.*
	rm -rf *.dSYM

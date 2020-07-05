# Copyright (C) 2012-2017 IBM Corp.
#
# This program is Licensed under the Apache License, Version 2.0
# (the "License"); you may not use this file except in compliance
# with the License. You may obtain a copy of the License at
#   http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License. See accompanying LICENSE file.
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

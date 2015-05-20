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
#
INCDIR  = ../include/
#
TARGET = libfhe.a
PREFIX = /usr/local/
#
CFLAGS = -g -O2 -Wfatal-errors -Wshadow -Wall -I$(PREFIX)/include/ -I$(INCDIR)

# useful flags:
#   -std=c++11
#   -std=c++0x
#   -DUSE_ALT_CRT  tells helib to use an alterntive to the default DoubleCRT
#                  representation...experimental...can be faster or slower
#                  depending on the operation mixture
#
#   -DNO_HALF_SIZE_PRIME  tells helib to not use the half size prime
#                         in the prime chain
#
#   -DEVALMAP_CACHED=0  tells helib to cache certain constants as ZZX's
#   -DEVALMAP_CACHED=1  tells helib to cache certain constants as DoubleCRT's
#                       these flags only affect bootstrapping
#
#   -DFHE_THREADS  tells helib to enable generic multithreading capabilities;
#                  must be used with a thread-enabled NTL and the -pthread
#                  flag should be passed to gcc
#
#   -DFHE_DCRT_THREADS  tells helib to use a multithreading strategy at the
#                       DoubleCRT level; requires -DFHE_THREADS (see above)
#
#   -DFHE_BOOT_THREADS  tells helib to use a multithreading strategy for
#                       bootstrapping; requires -DFHE_THREADS (see above)

#  If you get compilation errors, you may need to add -std=c++11 or -std=c++0x
# CFLAGS = -O2 -I$(INCDIR)

$(info HElib requires NTL version 9.0.1 or higher, see http://shoup.net/ntl)
$(info If you get compilation errors, try to add/remove -std=c++11 in Makefile)
$(info )

LD = g++
AR = ar
ARFLAGS=ruv
GMP=-lgmp
LDLIBS = -L/usr/local/lib -lntl $(GMP) -lm

_HEADER = EncryptedArray.h FHE.h Ctxt.h CModulus.h FHEContext.h PAlgebra.h DoubleCRT.h NumbTh.h bluestein.h IndexSet.h timing.h IndexMap.h replicate.h hypercube.h matching.h powerful.h permutations.h polyEval.h multicore.h
SRC = KeySwitching.cpp EncryptedArray.cpp FHE.cpp Ctxt.cpp CModulus.cpp FHEContext.cpp PAlgebra.cpp DoubleCRT.cpp NumbTh.cpp bluestein.cpp IndexSet.cpp timing.cpp replicate.cpp hypercube.cpp matching.cpp powerful.cpp BenesNetwork.cpp permutations.cpp PermNetwork.cpp OptimizePermutations.cpp eqtesting.cpp polyEval.cpp extractDigits.cpp EvalMap.cpp OldEvalMap.cpp recryption.cpp debugging.cpp
OBJ = NumbTh.o timing.o bluestein.o PAlgebra.o CModulus.o FHEContext.o IndexSet.o DoubleCRT.o FHE.o KeySwitching.o Ctxt.o EncryptedArray.o replicate.o hypercube.o matching.o powerful.o BenesNetwork.o permutations.o PermNetwork.o OptimizePermutations.o eqtesting.o polyEval.o extractDigits.o EvalMap.o OldEvalMap.o recryption.o debugging.o
TESTPROGS = Test_General_x Test_PAlgebra_x Test_IO_x Test_Replicate_x Test_LinPoly_x Test_matmul_x Test_Powerful_x Test_Permutations_x Test_Timing_x Test_PolyEval_x Test_extractDigits_x Test_EvalMap_x Test_bootstrapping_x
CHECKPROGS = Test_General_x Test_LinPoly_x Test_Permutations_x Test_PolyEval_x Test_Replicate_x Test_EvalMap_x Test_extractDigits_x Test_bootstrapping_x

HEADER = $(patsubst %,$(INCDIR)/FHE/%, $(_HEADER))


all: $(TARGET)

check: $(CHECKPROGS)
	./Test_General_x R=1 m=255
	./Test_General_x R=1 k=10 p=2 r=1 d=1
	./Test_General_x R=1 k=10 p=2 r=2 d=1
	./Test_General_x R=1 k=10 p=2 r=1 d=2
	./Test_General_x R=1 k=10 p=3 r=1 d=1
	./Test_General_x R=2 k=10 p=3 r=2 d=1 c=3
	./Test_General_x R=1 k=10 p=3 r=1 d=2 c=3
	./Test_LinPoly_x
	./Test_Permutations_x
	./Test_PolyEval_x
	./Test_Replicate_x m=1247
	./Test_EvalMap_x mvec="[7 3 221]" gens="[3979 3095 3760]" ords="[6 2 -8]"
	./Test_extractDigits_x m=2047 p=5
	./Test_bootstrapping_x

test: $(TESTPROGS)

obj: $(OBJ)

DoubleCRT.o: DoubleCRT.cpp AltCRT.cpp $(HEADER)
	$(CC) $(CFLAGS) -c DoubleCRT.cpp

%.o: %.cpp $(HEADER)
	$(CC) $(CFLAGS) -c $<

$(TARGET): $(OBJ)
	$(AR) $(ARFLAGS) $(TARGET) $(OBJ)

./%_x: %.cpp $(TARGET)
	$(CC) $(CFLAGS) -o $@ $< $(TARGET) $(LDLIBS)

install: $(TARGET)
	cp -r $(INCDIR) -t $(PREFIX)
	cp $(TARGET) $(PREFIX)/lib

uninstall:
	rm -r $(PREFIX)/include/FHE
	rm $(PREFIX)/lib/$(TARGET)

clean:
	rm -f *.o *_x *_x.exe *.a core.*
	rm -rf *.dSYM

info:
	: HElib require NTL version 9.0.1 or higher
	: Compilation flags are 'CFLAGS=$(CFLAGS)'
	: If errors occur, try adding/removing '-std=c++11' in Makefile
	:

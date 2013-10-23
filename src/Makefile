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
#CFLAGS = -g -Wall -I/usr/local/include
CFLAGS = -g -O2 -Wfatal-errors -Wshadow -Wall -I/usr/local/include 
#CFLAGS= -g -O2 -Wfatal-errors -Wshadow -Wall -I/usr/local/include -DUSE_ALT_CRT
#CFLAGS = -pg -Wall -I/usr/local/include
LD = g++


AR = ar
ARFLAGS=ruv

GMP=-lgmp 
LDLIBS = -L/usr/local/lib -lntl $(GMP) -lm

HEADER = EncryptedArray.h FHE.h Ctxt.h CModulus.h FHEContext.h PAlgebra.h SingleCRT.h DoubleCRT.h NumbTh.h bluestein.h IndexSet.h timing.h IndexMap.h replicate.h hypercube.h matching.h powerful.h permutations.h

SRC = KeySwitching.cpp EncryptedArray.cpp FHE.cpp Ctxt.cpp CModulus.cpp FHEContext.cpp PAlgebra.cpp SingleCRT.cpp DoubleCRT.cpp NumbTh.cpp PAlgebraMod.cpp bluestein.cpp IndexSet.cpp timing.cpp replicate.cpp hypercube.cpp matching.cpp powerful.cpp BenesNetwork.cpp permutations.cpp PermNetwork.cpp OptimizePermutaions.cpp

#OBJ = EncryptedArray.o FHE.o Ctxt.o CModulus.o FHEContext.o PAlgebra.o SingleCRT.o DoubleCRT.o NumbTh.o bluestein.o IndexSet.o timing.o KeySwitching.o PAlgebraMod.o

OBJ = NumbTh.o timing.o bluestein.o PAlgebra.o  CModulus.o FHEContext.o IndexSet.o DoubleCRT.o SingleCRT.o FHE.o KeySwitching.o Ctxt.o EncryptedArray.o replicate.o hypercube.o matching.o powerful.o BenesNetwork.o permutations.o PermNetwork.o OptimizePermutaions.o

#TESTPROGS = Test_PAlgebra_x Test_DoubleCRT_x Test_CModulus_x Test_FHE_x Test_Arrays_x
TESTPROGS = Test_General_x Test_Replicate_x Test_LinPoly_x Test_matmul_x Test_Powerful_x Test_Permutations_x


all: fhe.a

check: Test_General_x
	./Test_General_x R=1 m=255
	./Test_General_x R=1 k=10 p=2 r=1 d=1
	./Test_General_x R=1 k=10 p=2 r=2 d=1
	./Test_General_x R=1 k=10 p=2 r=1 d=2
	./Test_General_x R=1 k=10 p=3 r=1 d=1
	./Test_General_x R=2 k=10 p=3 r=2 d=1 c=3
	./Test_General_x R=1 k=10 p=3 r=1 d=2 c=3

test: $(TESTPROGS)

obj: $(OBJ)

%.o: %.cpp $(HEADER)
	$(CC) $(CFLAGS) -c $<

fhe.a: $(OBJ)
	$(AR) $(ARFLAGS) fhe.a $(OBJ)

./%_x: %.cpp fhe.a
	$(CC) $(CFLAGS) -o $@ $< fhe.a $(LDLIBS)

clean:
	rm -f *.o *_x *_x.exe *.a core.*
	rm -rf *.dSYM

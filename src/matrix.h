/* Copyright (C) 2012,2013 IBM Corp.
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */
#ifndef FHE_matrix_H_
#define FHE_matrix_H_
/**
 * @file matrix.h
 * @brief some matrix / linear algenra stuff
 */

#include "EncryptedArray.h"

void free_mat_mul(const EncryptedArray& ea, Ctxt& ctxt, const PlaintextBlockMatrixBaseInterface& mat);
void free_mat_mul(const EncryptedArray& ea, Ctxt& ctxt, const PlaintextMatrixBaseInterface& mat);
void free_mat_mul_dense(const EncryptedArray& ea, Ctxt& ctxt, const PlaintextMatrixBaseInterface& mat);

void free_compMat(const EncryptedArray& ea, 
                 CachedPtxtMatrix& cmat, const PlaintextMatrixBaseInterface& mat);

void free_compMat(const EncryptedArray& ea, 
                 CachedDCRTPtxtMatrix& cmat, const PlaintextMatrixBaseInterface& mat);

#endif /* ifdef FHE_matrix_H_ */


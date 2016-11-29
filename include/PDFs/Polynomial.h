/*
 * Polynomial.h
 *
 *  Created on: 22/11/2016
 *      Author: hschrein
 */

#pragma once

#include <hydra/Function.h>
#include <hydra/Types.h>
#include <hydra/Parameter.h>
#include <initializer_list>
#include <PDFs/BasePDF.h>
#include <thrust/functional.h>

namespace hydra {

namespace pdfs{

template<size_t N>
struct Polynomial : public BasePDF<
             Polynomial<N>,   // Curiously Recurring Tempate Pattern (CRTP) - this is the current class
             GReal_t, // The type
             N        // The number of parameters
             > {


	GUInt_t  fPosition;
    
    // Standard constructor
	Polynomial(
            std::array<Parameter, N> vals,
            GUInt_t position=0 ):
        fPosition(position),
		BasePDF<Polynomial<N>,GReal_t,N>(vals) {
        }

    // Copy constructor
	__host__ __device__
	inline Polynomial(Polynomial<N> const& other):
	    BasePDF<Polynomial<N>,GReal_t,N>(other),
	    fPosition(other.fPosition) {}


    // Assignment operator
	__host__ __device__
	inline Polynomial& operator=(Polynomial<N> const& other) {
		if(this == &other) return *this;

		BasePDF<Polynomial<N>,GReal_t,N>::operator=(other);
		this->fPosition = other.fPosition;

		return *this;
	}

    // The evaluate method
	template<typename T>
	__host__ __device__
	inline GReal_t Evaluate(T* x, T* p=0) {
        T xval = x[fPosition];

        GReal_t value = 0;
        for(size_t i = 1; i<=N; i++)
            value += pow(xval, i) * this->fParams[i];
        
        thrust::maximum<GReal_t> mx;
        return mx(value,0);
	}
};
}
}

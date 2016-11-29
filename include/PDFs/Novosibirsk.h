/*
 * NovosibirskPdf.h
 *
 *  Created on: 22/11/2016
 *      Author: hschrein
 */

#pragma once

#include <hydra/Types.h>
#include <hydra/Parameter.h>
#include "PDFs/BasePDF.h"
#include <initializer_list>

namespace hydra {

namespace pdfs{

class Novosibirsk:public BasePDF<
             Novosibirsk,   // Curiously Recurring Tempate Pattern (CRTP) - this is the current class
             GReal_t, // The type
             3        // The number of parameters
             > {


	GUInt_t  fPosition;

public:
    // Standard constructor
	Novosibirsk(
            Parameter const& mean,
            Parameter const& sigma,
            Parameter const& tail,
            GUInt_t position=0 ):
		BasePDF<Novosibirsk,GReal_t,3>({mean, sigma, tail}),
		fPosition(position) {}

    // Copy constructor
	__host__ __device__
	inline Novosibirsk(Novosibirsk const& other):
	    BasePDF<Novosibirsk,GReal_t,3>(other),
	    fPosition(other.fPosition) {}


    // Assignment operator
	__host__ __device__
	inline Novosibirsk& operator=( Novosibirsk const& other) {
		if(this == &other) return *this;

		BasePDF<Novosibirsk,GReal_t,3>::operator=(other);
		this->fPosition = other.fPosition;

		return *this;
	}

    // The evaluate method
	template<typename T>
	__host__ __device__
	inline GReal_t Evaluate(T* x, T* p=0) {
        T xval = x[fPosition];
        Parameter& Mean  = fParams[0];
        Parameter& Sigma = fParams[1];
        Parameter& Tail  = fParams[2];
        
        if(fabs(Tail) < 1.e-7)
            return exp(.5 * (xval-Mean)/Sigma * (xval-Mean)/Sigma);

        T arg = 1. - (xval-Mean) * Tail / Sigma;

        if(arg < 1.e-7)
            return 0.;

        T log_val = log(arg);
        T xi = 2*sqrt(log(4.));

        T width_zero = (2. / xi) * asinh(Tail * xi * .5 );
        T width_zero2 = width_zero * width_zero;
        T exponent = ( -0.5 / (width_zero2) * log_val * log_val) - ( width_zero2 * 0.5 );
        return exp(exponent);
	}
};
}
}

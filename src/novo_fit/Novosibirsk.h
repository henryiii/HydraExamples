/*
 * NovosibirskPdf.h
 *
 *  Created on: 22/11/2016
 *      Author: hschrein
 */

#pragma ONCE

#include <hydra/Function.h>
#include <hydra/Types.h>
#include <hydra/Parameter.h>
#include <initializer_list>

namespace hydra {

namespace pdfs{

class Novosibirsk:public BaseFunctor<
             Novosibirsk,   // Curiously Recurring Tempate Pattern (CRTP) - this is the current class
             GReal_t, // The type
             3        // The number of parameters
             > {


	GUInt_t  fPosition;
	Parameter fMean;
	Parameter fSigma;
    Parameter fTail;
public:
    // Standard constructor
	Novosibirsk(
            Parameter const& mean,
            Parameter const& sigma,
            Parameter const& tail,
            GUInt_t position=0 ):
		BaseFunctor<Novosibirsk,GReal_t,3>(),
		fPosition(position),
		fMean(mean),
		fSigma(sigma),
        fTail(tail) {
            RegistryParameters({&fMean, &fSigma, &fTail});
        }

    // Copy constructor
	__host__ __device__
	inline Novosibirsk(Novosibirsk const& other):
	    BaseFunctor<Novosibirsk,GReal_t,3>(other),
	    fPosition(other.fPosition),
	    fMean(other.fMean),
	    fSigma(other.fSigma),
        fTail(other.fTail) {
            RegistryParameters({&(this->fMean), &(this->fSigma), &(this->fTail)});
        }


    // Assignment operator
	__host__ __device__
	inline Novosibirsk& operator=( Novosibirsk const& other) {
		if(this == &other) return *this;

		BaseFunctor<Novosibirsk,GReal_t,3>::operator=(other);
		this->fMean = other.fMean;
		this->fSigma = other.fSigma;
		this->fTail = other.fTail;
		this->fPosition = other.fPosition;
		this->RegistryParameters({&(this->fMean), &(this->fSigma), &(this->fTail)});

		return *this;
	}

    // The evaluate method
	template<typename T>
	__host__ __device__
	inline GReal_t Evaluate(T* x, T* p=0) {
        T xval = x[fPosition];
        
        if(fabs(fTail) < 1.e-7)
            return exp(.5 * (xval-fMean)/fSigma * (xval-fMean)/fSigma);

        T arg = 1. - (xval-fMean) * fTail / fSigma;

        if(arg < 1.e-7)
            return 0.;

        T log_val = log(arg);
        T xi = 2*sqrt(log(4.));

        T width_zero = (2. / xi) * asinh(fTail * xi * .5 );
        T width_zero2 = width_zero * width_zero;
        T exponent = ( -0.5 / (width_zero2) * log_val * log_val) - ( width_zero2 * 0.5 );
        return exp(exponent);
	}

};

}
}

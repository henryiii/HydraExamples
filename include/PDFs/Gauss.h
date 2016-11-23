/*
 * Gauss.h
 *
 *  Created on: 26/08/2016
 *      Author: augalves
 */

#pragma ONCE

#include <hydra/Function.h>
#include <hydra/Types.h>
#include <hydra/Parameter.h>
#include <initializer_list>

namespace hydra {

namespace pdfs{

class Gauss:public BaseFunctor<
             Gauss,   // Curiously Recurring Tempate Pattern (CRTP) - this is the current class
             GReal_t, // The type
             2        // The number of parameters
             > {


	GUInt_t  fPosition;
	Parameter fM;
	Parameter fS;
public:
    // Standard constructor
	Gauss(
            Parameter const& mean,
            Parameter const& sigma,
            GUInt_t position=0 ):
		BaseFunctor<Gauss,GReal_t,2>(),
		fPosition(position),
		fM(mean),
		fS(sigma) {
            RegistryParameters({&fM, &fS});
        }

    // Copy constructor
	__host__ __device__
	inline Gauss(Gauss const& other):
	    BaseFunctor<Gauss,GReal_t,2>(other),
	    fPosition(other.fPosition),
	    fM(other.fM),
	    fS(other.fS) {
            RegistryParameters({&(this->fM), &(this->fS)});
        }


    // Assignment operator
	__host__ __device__
	inline Gauss& operator=( Gauss const& other) {
		if(this == &other) return *this;

		BaseFunctor<Gauss,GReal_t,2>::operator=(other);
		this->fM = other.fM;
		this->fS = other.fS;
		this->fPosition = other.fPosition;
		this->RegistryParameters({&(this->fM), &(this->fS)});

		return *this;
	}

    // The evaluate method
	template<typename T>
	__host__ __device__
	inline GReal_t Evaluate(T* x, T* p=0) {
		return exp(-((x[fPosition] - fM) * (x[fPosition] - fM)) / (2 * fS * fS));
	}

};

}
}

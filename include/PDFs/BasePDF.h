/*
 * BasePDF.h
 *
 *  Created on: 29/11/2016
 *      Author: hschrein
 */

#pragma once

#include <hydra/Function.h>
#include <hydra/Types.h>
#include <hydra/Parameter.h>
#include <initializer_list>


namespace hydra {

namespace pdfs {


template<typename Functor, typename ReturnType, size_t NPARAM>
struct BasePDF : public BaseFunctor<Functor, ReturnType, NPARAM> {

    std::array<Parameter, NPARAM> fParams;

    void RegistryParametersArray() {
        #ifndef __CUDA_ARCH__
		int i=0;
		for(Parameter &var: fParams) {
			this->fParameters[i]=&var;
			i++;
		}
		this->fParamResgistered=1;
        #endif
    }

    // Standard constructor
	BasePDF(std::array<Parameter, NPARAM> par):
		BaseFunctor<Functor, ReturnType, NPARAM>(), fParams(par) {
            this->RegistryParametersArray();
        }

    // Copy constructor
	__host__ __device__
	inline BasePDF(BasePDF const& other):
	    BaseFunctor<Functor, ReturnType, NPARAM>(other),
	    fParams(other.fParams) {
            this->RegistryParametersArray();
        }


    // Assignment operator
	__host__ __device__
	inline BasePDF& operator=( BasePDF const& other) {
		if(this == &other) return *this;

		BaseFunctor<Functor, ReturnType, NPARAM>::operator=(other);
		this->fParams = other.fParams;
        this->RegistryParametersArray();

		return *this;
	}

};
}
}

#include <string>
#include <algorithm>
#include <hydra/Parameter.h>
#include "PDFs/Novosibirsk.h"

namespace HyFit {

class RealParam : public hydra::Parameter {
    static UserParameters upar;
    std::string var_name;

    public:
        template<typename... T>
        RealParam(std::string name, T... args):
            var_name(name), Parameter(var_name, std::forward<T>(args...)) {
            
            upar.AddParameter(this);
        }

}


class RealVar {
    static unsigned int s_key = 0;
    hydra::GReal_t m_min;
    hydra::GReal_t m_max;
    hydra::GUInt_t m_key;
    std::string m_name;
public:
    RealVar(std::string name, hydra::GReal_t min, hydra::GReal_t max) :
        m_name(name), m_key(s_key++), m_min(min), m_max(max) {}

    RealVar(std::string name, unsigned int key, hydra::GReal_t min, hydra::GReal_t max) :
        m_name(name), m_key(key), m_min(min), m_max(max) {}

    hydra::GUInt_t getKey() const {return m_key;}
    hydra::GUInt_t getMin() const {return m_min;}
    hydra::GUInt_t getMax() const {return m_max;}
}



class NovoPDF : public hydra::PDFs::Novosibirsk {
    std::string m_name;
    
    NovoPDF(std::string name, RealVar &x, RealParam &mean, RealParam &sigma, RealParam &tail):
                m_name(name), Novosibirsk(mean, sigma, tail, x.getKey())  {}

}




}

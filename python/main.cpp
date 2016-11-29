#include <hydra/Parameter.h>
#include <PDFs/Novosibirsk.h>
#include <pybind11/pybind11.h>
#include <sstream>

namespace py = pybind11;
using namespace hydra;
using namespace pybind11::literals;

struct Variable : public Parameter {
    std::string m_name_string;

    template<typename... T>
    Variable(const std::string& input, T&&... args) : m_name_string(input), Parameter(m_name_string, std::forward<T>(args)...) {}

    void SetName(const std::string& name) {m_name_string = name; return Parameter::SetName(m_name_string);}
    std::string GetName() const {return m_name_string;}
};

PYBIND11_PLUGIN(hydra_example) {
    py::module m("hydra_example", "An example of pybind");

    py::class_<Variable>(m, "Variable")
        .def(py::init<std::string const&, GReal_t, GReal_t, GReal_t, GReal_t>(),"name"_a, "value"_a, "error"_a, "downlim"_a, "uplim"_a)
        .def(py::init<std::string const&, GReal_t, GReal_t>(),"name"_a, "value"_a, "error"_a)
        .def(py::init<std::string const&, GReal_t>(),"name"_a, "value"_a)
        .def("__repr__", [](const Variable &p){std::stringstream os; os << "<" << p << ">"; return os.str();})
        .def_property("name", &Variable::GetName, &Variable::SetName)
        ;

    py::class_<pdfs::Novosibirsk>(m, "Novosibirsk")
        .def(py::init<Parameter const&, Parameter const&, Parameter const&, GUInt_t>(), "mean"_a, "sigma"_a, "tail"_a, "position"_a)
        .def("__repr__", [](const pdfs::Novosibirsk& n){return std::string("<Novosibirsk>");})
        ;

    return m.ptr();
}




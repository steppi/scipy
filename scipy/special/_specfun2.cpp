#include "special/specfun.h"
#include "ufunc.h"

extern const char *bei_doc;
extern const char *beip_doc;
extern const char *ber_doc;
extern const char *berp_doc;
extern const char *exp1_doc;
extern const char *expi_doc;

static PyModuleDef _specfun2_def = {
    PyModuleDef_HEAD_INIT,
    .m_name = "_specfun2",
    .m_size = -1,
};

PyMODINIT_FUNC PyInit__specfun2() {
    SpecFun_Initialize();

    PyObject *specfun2 = PyModule_Create(&_specfun2_def);
    if (specfun2 == nullptr) {
        return nullptr;
    }

    PyObject *bei = SpecFun_UFunc<special::bei>("bei", bei_doc);
    PyModule_AddObjectRef(specfun2, "bei", bei);

    PyObject *beip = SpecFun_UFunc<special::beip>("beip", beip_doc);
    PyModule_AddObjectRef(specfun2, "beip", beip);

    PyObject *ber = SpecFun_UFunc<special::ber>("ber", ber_doc);
    PyModule_AddObjectRef(specfun2, "ber", ber);

    PyObject *berp = SpecFun_UFunc<special::berp>("berp", berp_doc);
    PyModule_AddObjectRef(specfun2, "berp", berp);

    PyObject *exp1 = SpecFun_UFunc<special::exp1, special::cexp1>("exp1", exp1_doc);
    PyModule_AddObjectRef(specfun2, "exp1", exp1);

    PyObject *expi = SpecFun_UFunc<special::expi, special::cexpi>("expi", expi_doc);
    PyModule_AddObjectRef(specfun2, "expi", expi);

    return specfun2;
}

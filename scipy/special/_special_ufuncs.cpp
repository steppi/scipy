#include <cmath>

#include "_special.h"
#include "amos_wrappers.h"
#include "special/beta.h"
#include "special/gamma.h"
#include "special/lgamma.h"
#include "special/specfun.h"
#include "special/trig.h"
#include "special/zeta.h"
#include "ufunc.h"

// This is needed by sf_error, it is defined in "_ufuncs_extra_code_common.pxi" for "_generate_pyx.py".
// It claims to exist to call "call PyUFunc_getfperr in a context where PyUFunc_API array is initialized", but we are
// already in such a context
extern "C" int wrap_PyUFunc_getfperr() { return PyUFunc_getfperr(); }

extern const char *_cosine_cdf_doc;
extern const char *_cosine_invcdf_doc;
extern const char *_cospi_doc;
extern const char *_sinpi_doc;
extern const char *_zeta_doc;
extern const char *airy_doc;
extern const char *airye_doc;
extern const char *bei_doc;
extern const char *beip_doc;
extern const char *ber_doc;
extern const char *berp_doc;
extern const char *beta_doc;
extern const char *betaln_doc;
extern const char *exp1_doc;
extern const char *expi_doc;
extern const char *gamma_doc;
extern const char *gammaln_doc;
extern const char *hankel1_doc;
extern const char *hankel1e_doc;
extern const char *hankel2_doc;
extern const char *hankel2e_doc;
extern const char *it2i0k0_doc;
extern const char *it2j0y0_doc;
extern const char *it2struve0_doc;
extern const char *itairy_doc;
extern const char *iti0k0_doc;
extern const char *itj0y0_doc;
extern const char *itmodstruve0_doc;
extern const char *itstruve0_doc;
extern const char *kei_doc;
extern const char *keip_doc;
extern const char *kelvin_doc;
extern const char *ker_doc;
extern const char *kerp_doc;
extern const char *mathieu_a_doc;
extern const char *mathieu_b_doc;
extern const char *mathieu_cem_doc;
extern const char *mathieu_modcem1_doc;
extern const char *mathieu_modcem2_doc;
extern const char *mathieu_modsem1_doc;
extern const char *mathieu_modsem2_doc;
extern const char *mathieu_sem_doc;
extern const char *modfresnelm_doc;
extern const char *modfresnelp_doc;

static PyModuleDef _special_ufuncs_def = {
    PyModuleDef_HEAD_INIT,
    .m_name = "_special_ufuncs",
    .m_size = -1,
};

PyMODINIT_FUNC PyInit__special_ufuncs() {
    SpecFun_Initialize();

    PyObject *_special_ufuncs = PyModule_Create(&_special_ufuncs_def);
    if (_special_ufuncs == nullptr) {
        return nullptr;
    }

    PyObject *_cospi = SpecFun_UFunc<special::cephes::cospi<float>, special::cephes::cospi<double>,
                                     special::cospi<float>, special::cospi<double>>("_cospi", _cospi_doc);
    PyModule_AddObjectRef(_special_ufuncs, "_cospi", _cospi);

    PyObject *_sinpi = SpecFun_UFunc<special::cephes::sinpi<float>, special::cephes::sinpi<double>,
                                     special::sinpi<float>, special::sinpi<double>>("_sinpi", _sinpi_doc);
    PyModule_AddObjectRef(_special_ufuncs, "_sinpi", _sinpi);

    PyObject *_zeta = SpecFun_UFunc<special::zeta<float>, special::zeta<double>>("_zeta", _zeta_doc);
    PyModule_AddObjectRef(_special_ufuncs, "_zeta", _zeta);

    PyObject *bei = SpecFun_UFunc<special::bei<float>, special::bei<double>>("bei", bei_doc);
    PyModule_AddObjectRef(_special_ufuncs, "bei", bei);

    PyObject *beip = SpecFun_UFunc<special::beip<float>, special::beip<double>>("beip", beip_doc);
    PyModule_AddObjectRef(_special_ufuncs, "beip", beip);

    PyObject *ber = SpecFun_UFunc<special::ber<float>, special::ber<double>>("ber", ber_doc);
    PyModule_AddObjectRef(_special_ufuncs, "ber", ber);

    PyObject *berp = SpecFun_UFunc<special::berp<float>, special::berp<double>>("berp", berp_doc);
    PyModule_AddObjectRef(_special_ufuncs, "berp", berp);

    PyObject *beta = SpecFun_UFunc<special::beta<float>, special::beta<double>>("beta", beta_doc);
    PyModule_AddObjectRef(_special_ufuncs, "beta", beta);

    PyObject *betaln = SpecFun_UFunc<special::lbeta<float>, special::lbeta<double>>("betaln", betaln_doc);
    PyModule_AddObjectRef(_special_ufuncs, "betaln", betaln);

    PyObject *exp1 = SpecFun_UFunc<special::exp1<float>, special::exp1<double>, special::cexp1>("exp1", exp1_doc);
    PyModule_AddObjectRef(_special_ufuncs, "exp1", exp1);

    PyObject *expi = SpecFun_UFunc<special::expi<float>, special::expi<double>, special::cexpi>("expi", expi_doc);
    PyModule_AddObjectRef(_special_ufuncs, "expi", expi);

    PyObject *gamma = SpecFun_UFunc<special::gamma<float>, special::gamma<double>, cgamma>("gamma", gamma_doc);
    PyModule_AddObjectRef(_special_ufuncs, "gamma", gamma);

    PyObject *gammaln = SpecFun_UFunc<special::lgam<float>, special::lgam<double>>("gammaln", gammaln_doc);
    PyModule_AddObjectRef(_special_ufuncs, "gammaln", gammaln);

    PyObject *it2i0k0 = SpecFun_UFunc<special::it2i0k0<float>, special::it2i0k0<double>>("it2i0k0", it2i0k0_doc, 2);
    PyModule_AddObjectRef(_special_ufuncs, "it2i0k0", it2i0k0);

    PyObject *it2j0y0 = SpecFun_UFunc<special::it2j0y0<float>, special::it2j0y0<double>>("it2j0y0", it2j0y0_doc, 2);
    PyModule_AddObjectRef(_special_ufuncs, "it2j0y0", it2j0y0);

    PyObject *it2struve0 =
        SpecFun_UFunc<special::it2struve0<float>, special::it2struve0<double>>("it2struve0", it2struve0_doc);
    PyModule_AddObjectRef(_special_ufuncs, "it2struve0", it2struve0);

    PyObject *itairy = SpecFun_UFunc<special::itairy<float>, special::itairy<double>>("itairy", itairy_doc, 4);
    PyModule_AddObjectRef(_special_ufuncs, "itairy", itairy);

    PyObject *iti0k0 = SpecFun_UFunc<special::it1i0k0<float>, special::it1i0k0<double>>("iti0k0", iti0k0_doc, 2);
    PyModule_AddObjectRef(_special_ufuncs, "iti0k0", iti0k0);

    PyObject *itj0y0 = SpecFun_UFunc<special::it1j0y0<float>, special::it1j0y0<double>>("itj0y0", itj0y0_doc, 2);
    PyModule_AddObjectRef(_special_ufuncs, "itj0y0", itj0y0);

    PyObject *itmodstruve0 =
        SpecFun_UFunc<special::itmodstruve0<float>, special::itmodstruve0<double>>("itmodstruve0", itmodstruve0_doc);
    PyModule_AddObjectRef(_special_ufuncs, "itmodstruve0", itmodstruve0);

    PyObject *itstruve0 =
        SpecFun_UFunc<special::itstruve0<float>, special::itstruve0<double>>("itstruve0", itstruve0_doc);
    PyModule_AddObjectRef(_special_ufuncs, "itstruve0", itstruve0);

    PyObject *kei = SpecFun_UFunc<special::kei<float>, special::kei<double>>("kei", kei_doc);
    PyModule_AddObjectRef(_special_ufuncs, "kei", kei);

    PyObject *keip = SpecFun_UFunc<special::keip<float>, special::keip<double>>("keip", keip_doc);
    PyModule_AddObjectRef(_special_ufuncs, "keip", keip);

    PyObject *kelvin = SpecFun_UFunc<special::kelvin<float>, special::kelvin<double>>("kelvin", kelvin_doc, 4);
    PyModule_AddObjectRef(_special_ufuncs, "kelvin", kelvin);

    PyObject *ker = SpecFun_UFunc<special::ker<float>, special::ker<double>>("ker", ker_doc);
    PyModule_AddObjectRef(_special_ufuncs, "ker", ker);

    PyObject *kerp = SpecFun_UFunc<special::kerp<float>, special::kerp<double>>("kerp", kerp_doc);
    PyModule_AddObjectRef(_special_ufuncs, "kerp", kerp);

    PyObject *mathieu_a = SpecFun_UFunc<special::cem_cva<float>, special::cem_cva<double>>("mathieu_a", mathieu_a_doc);
    PyModule_AddObjectRef(_special_ufuncs, "mathieu_a", mathieu_a);

    PyObject *mathieu_b = SpecFun_UFunc<special::sem_cva<float>, special::sem_cva<double>>("mathieu_b", mathieu_b_doc);
    PyModule_AddObjectRef(_special_ufuncs, "mathieu_b", mathieu_b);

    PyObject *mathieu_cem = SpecFun_UFunc<special::cem<float>, special::cem<double>>("mathieu_cem", mathieu_cem_doc, 2);
    PyModule_AddObjectRef(_special_ufuncs, "mathieu_cem", mathieu_cem);

    PyObject *mathieu_modcem1 =
        SpecFun_UFunc<special::mcm1<float>, special::mcm1<double>>("mathieu_modcem1", mathieu_modcem1_doc, 2);
    PyModule_AddObjectRef(_special_ufuncs, "mathieu_modcem1", mathieu_modcem1);

    PyObject *mathieu_modcem2 =
        SpecFun_UFunc<special::mcm2<float>, special::mcm2<double>>("mathieu_modcem2", mathieu_modcem2_doc, 2);
    PyModule_AddObjectRef(_special_ufuncs, "mathieu_modcem2", mathieu_modcem2);

    PyObject *mathieu_modsem1 =
        SpecFun_UFunc<special::msm1<float>, special::msm1<double>>("mathieu_modsem1", mathieu_modsem1_doc, 2);
    PyModule_AddObjectRef(_special_ufuncs, "mathieu_modsem1", mathieu_modsem1);

    PyObject *mathieu_modsem2 =
        SpecFun_UFunc<special::msm2<float>, special::msm2<double>>("mathieu_modsem2", mathieu_modsem2_doc, 2);
    PyModule_AddObjectRef(_special_ufuncs, "mathieu_modsem2", mathieu_modsem2);

    PyObject *mathieu_sem = SpecFun_UFunc<special::sem<float>, special::sem<double>>("mathieu_sem", mathieu_sem_doc, 2);
    PyModule_AddObjectRef(_special_ufuncs, "mathieu_sem", mathieu_sem);

    PyObject *modfresnelm =
        SpecFun_UFunc<special::modified_fresnel_minus<float>, special::modified_fresnel_minus<double>>(
            "modfresnelm", modfresnelm_doc, 2);
    PyModule_AddObjectRef(_special_ufuncs, "modfresnelm", modfresnelm);

    PyObject *modfresnelp =
        SpecFun_UFunc<special::modified_fresnel_plus<float>, special::modified_fresnel_plus<double>>(
            "modfresnelp", modfresnelp_doc, 2);
    PyModule_AddObjectRef(_special_ufuncs, "modfresnelp", modfresnelp);

    return _special_ufuncs;
}

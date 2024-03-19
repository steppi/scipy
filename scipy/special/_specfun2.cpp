#include <cmath>

extern "C" {

#include "_cosine.h"
#include "amos_wrappers.h"
}

#include "_special.h"
#include "special/cephes/beta.h"
#include "special/cephes/gamma.h"
#include "special/cephes/zeta.h"
#include "special/lgamma.h"
#include "special/specfun.h"
#include "special/trig.h"
#include "ufunc.h"

// This is needed by sf_error, it is defined in "_ufuncs_extra_code_common.pxi" for "_generate_pyx.py".
// It claims to exist to call "call PyUFunc_getfperr in a context where PyUFunc_API array is initialized", but we are
// already in such a context
extern "C" int wrap_PyUFunc_getfperr() { return PyUFunc_getfperr(); }

extern const char *_cosine_cdf_doc;
extern const char *_cosine_invcdf_doc;
extern const char *airy_doc;
extern const char *airye_doc;
extern const char *bei_doc;
extern const char *beip_doc;
extern const char *ber_doc;
extern const char *berp_doc;
extern const char *exp1_doc;
extern const char *expi_doc;
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

    PyObject *_cosine_cdf = SpecFun_UFunc<cosine_cdf>("_cosine_cdf", _cosine_cdf_doc);
    PyModule_AddObjectRef(specfun2, "_cosine_cdf", _cosine_cdf);

    PyObject *_cosine_invcdf = SpecFun_UFunc<cosine_invcdf>("_cosine_invcdf", _cosine_invcdf_doc);
    PyModule_AddObjectRef(specfun2, "_cosine_invcdf", _cosine_invcdf);

    PyObject *airy = SpecFun_UFunc<airy_wrap_v, cairy_wrap_v>("airy", airy_doc, 4);
    PyModule_AddObjectRef(specfun2, "airy", airy);

    PyObject *airye = SpecFun_UFunc<cairy_wrap_e_real_v, cairy_wrap_e_v>("airye", airye_doc, 4);
    PyModule_AddObjectRef(specfun2, "airye", airye);

    PyObject *bei = SpecFun_UFunc<special::bei<float>, special::bei<double>>("bei", bei_doc);
    PyModule_AddObjectRef(specfun2, "bei", bei);

    PyObject *beip = SpecFun_UFunc<special::beip<float>, special::beip<double>>("beip", beip_doc);
    PyModule_AddObjectRef(specfun2, "beip", beip);

    PyObject *ber = SpecFun_UFunc<special::ber<float>, special::ber<double>>("ber", ber_doc);
    PyModule_AddObjectRef(specfun2, "ber", ber);

    PyObject *berp = SpecFun_UFunc<special::berp<float>, special::berp<double>>("berp", berp_doc);
    PyModule_AddObjectRef(specfun2, "berp", berp);

    PyObject *beta = SpecFun_UFunc<special::cephes::beta>("beta", nullptr);
    PyModule_AddObjectRef(specfun2, "beta", beta);

    PyObject *betaln = SpecFun_UFunc<special::cephes::lbeta>("betaln", nullptr);
    PyModule_AddObjectRef(specfun2, "betaln", betaln);

    PyObject *cospi = SpecFun_UFunc<special::cephes::cospi, special::cospi>("cospi", nullptr);
    PyModule_AddObjectRef(specfun2, "cospi", cospi);

    PyObject *exp1 = SpecFun_UFunc<special::exp1<float>, special::exp1<double>, special::cexp1>("exp1", exp1_doc);
    PyModule_AddObjectRef(specfun2, "exp1", exp1);

    PyObject *expi = SpecFun_UFunc<special::expi<float>, special::expi<double>, special::cexpi>("expi", expi_doc);
    PyModule_AddObjectRef(specfun2, "expi", expi);

    PyObject *gamma = SpecFun_UFunc<special::cephes::Gamma, cgamma>("gamma", nullptr);
    PyModule_AddObjectRef(specfun2, "gamma", gamma);

    PyObject *gammaln = SpecFun_UFunc<special::lgam<float>, special::lgam<double>>("gammaln", nullptr);
    PyModule_AddObjectRef(specfun2, "gammaln", gammaln);

    PyObject *hankel1 = SpecFun_UFunc<cbesh_wrap1>("hankel1", hankel1_doc);
    PyModule_AddObjectRef(specfun2, "hankel1", hankel1);

    PyObject *hankel1e = SpecFun_UFunc<cbesh_wrap1_e>("hankel1e", hankel1e_doc);
    PyModule_AddObjectRef(specfun2, "hankel1e", hankel1e);

    PyObject *hankel2 = SpecFun_UFunc<cbesh_wrap2>("hankel2", hankel2_doc);
    PyModule_AddObjectRef(specfun2, "hankel2", hankel2);

    PyObject *hankel2e = SpecFun_UFunc<cbesh_wrap2_e>("hankel2e", hankel2e_doc);
    PyModule_AddObjectRef(specfun2, "hankel2e", hankel2e);

    PyObject *it2i0k0 = SpecFun_UFunc<special::it2i0k0<float>, special::it2i0k0<double>>("it2i0k0", it2i0k0_doc, 2);
    PyModule_AddObjectRef(specfun2, "it2i0k0", it2i0k0);

    PyObject *it2j0y0 = SpecFun_UFunc<special::it2j0y0<float>, special::it2j0y0<double>>("it2j0y0", it2j0y0_doc, 2);
    PyModule_AddObjectRef(specfun2, "it2j0y0", it2j0y0);

    PyObject *it2struve0 =
        SpecFun_UFunc<special::it2struve0<float>, special::it2struve0<double>>("it2struve0", it2struve0_doc);
    PyModule_AddObjectRef(specfun2, "it2struve0", it2struve0);

    PyObject *itairy = SpecFun_UFunc<special::itairy<float>, special::itairy<double>>("itairy", itairy_doc, 4);
    PyModule_AddObjectRef(specfun2, "itairy", itairy);

    PyObject *iti0k0 = SpecFun_UFunc<special::it1i0k0<float>, special::it1i0k0<double>>("iti0k0", iti0k0_doc, 2);
    PyModule_AddObjectRef(specfun2, "iti0k0", iti0k0);

    PyObject *itj0y0 = SpecFun_UFunc<special::it1j0y0<float>, special::it1j0y0<double>>("itj0y0", itj0y0_doc, 2);
    PyModule_AddObjectRef(specfun2, "itj0y0", itj0y0);

    PyObject *itmodstruve0 =
        SpecFun_UFunc<special::itmodstruve0<float>, special::itmodstruve0<double>>("itmodstruve0", itmodstruve0_doc);
    PyModule_AddObjectRef(specfun2, "itmodstruve0", itmodstruve0);

    PyObject *itstruve0 =
        SpecFun_UFunc<special::itstruve0<float>, special::itstruve0<double>>("itstruve0", itstruve0_doc);
    PyModule_AddObjectRef(specfun2, "itstruve0", itstruve0);

    PyObject *kei = SpecFun_UFunc<special::kei<float>, special::kei<double>>("kei", kei_doc);
    PyModule_AddObjectRef(specfun2, "kei", kei);

    PyObject *keip = SpecFun_UFunc<special::keip<float>, special::keip<double>>("keip", keip_doc);
    PyModule_AddObjectRef(specfun2, "keip", keip);

    PyObject *kelvin = SpecFun_UFunc<special::kelvin<float>, special::kelvin<double>>("kelvin", kelvin_doc, 4);
    PyModule_AddObjectRef(specfun2, "kelvin", kelvin);

    PyObject *ker = SpecFun_UFunc<special::ker<float>, special::ker<double>>("ker", ker_doc);
    PyModule_AddObjectRef(specfun2, "ker", ker);

    PyObject *kerp = SpecFun_UFunc<special::kerp<float>, special::kerp<double>>("kerp", kerp_doc);
    PyModule_AddObjectRef(specfun2, "kerp", kerp);

    PyObject *mathieu_a = SpecFun_UFunc<special::cem_cva<float>, special::cem_cva<double>>("mathieu_a", mathieu_a_doc);
    PyModule_AddObjectRef(specfun2, "mathieu_a", mathieu_a);

    PyObject *mathieu_b = SpecFun_UFunc<special::sem_cva<float>, special::sem_cva<double>>("mathieu_b", mathieu_b_doc);
    PyModule_AddObjectRef(specfun2, "mathieu_b", mathieu_b);

    PyObject *mathieu_cem = SpecFun_UFunc<special::cem<float>, special::cem<double>>("mathieu_cem", mathieu_cem_doc, 2);
    PyModule_AddObjectRef(specfun2, "mathieu_cem", mathieu_cem);

    PyObject *mathieu_modcem1 = SpecFun_UFunc<special::mcm1>("mathieu_modcem1", mathieu_modcem1_doc, 2);
    PyModule_AddObjectRef(specfun2, "mathieu_modcem1", mathieu_modcem1);

    PyObject *mathieu_modcem2 = SpecFun_UFunc<special::mcm2>("mathieu_modcem2", mathieu_modcem2_doc, 2);
    PyModule_AddObjectRef(specfun2, "mathieu_modcem2", mathieu_modcem2);

    PyObject *mathieu_modsem1 = SpecFun_UFunc<special::msm1>("mathieu_modsem1", mathieu_modsem1_doc, 2);
    PyModule_AddObjectRef(specfun2, "mathieu_modsem1", mathieu_modsem1);

    PyObject *mathieu_modsem2 = SpecFun_UFunc<special::msm2>("mathieu_modsem2", mathieu_modsem2_doc, 2);
    PyModule_AddObjectRef(specfun2, "mathieu_modsem2", mathieu_modsem2);

    PyObject *mathieu_sem = SpecFun_UFunc<special::sem<float>, special::sem<double>>("mathieu_sem", mathieu_sem_doc, 2);
    PyModule_AddObjectRef(specfun2, "mathieu_sem", mathieu_sem);

    PyObject *modfresnelm =
        SpecFun_UFunc<special::modified_fresnel_minus<float>, special::modified_fresnel_minus<double>>(
            "modfresnelm", modfresnelm_doc, 2);
    PyModule_AddObjectRef(specfun2, "modfresnelm", modfresnelm);

    PyObject *modfresnelp =
        SpecFun_UFunc<special::modified_fresnel_plus<float>, special::modified_fresnel_plus<double>>(
            "modfresnelp", modfresnelp_doc, 2);
    PyModule_AddObjectRef(specfun2, "modfresnelp", modfresnelp);

    PyObject *sinpi = SpecFun_UFunc<special::cephes::sinpi, special::sinpi>("sinpi", nullptr);
    PyModule_AddObjectRef(specfun2, "sinpi", sinpi);

    PyObject *_zeta = SpecFun_UFunc<special::cephes::zeta>("_zeta", nullptr);
    PyModule_AddObjectRef(specfun2, "_zeta", _zeta);

    return specfun2;
}

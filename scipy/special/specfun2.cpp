#include <array>
#include <map>
#include <tuple>

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_API_VERSION
#include <numpy/arrayobject.h>
#include <numpy/npy_3kcompat.h>
#include <numpy/ufuncobject.h>

#include "special/specfun.h"

// Just initializes everything needed, can also go in a common header
bool SpecFun_Initialize() {
    Py_Initialize();

    import_array();
    if (PyErr_Occurred()) {
        // import array failed
        return false;
    }

    import_umath();
    if (PyErr_Occurred()) {
        // import umath failed
        return false;
    }

    return true;
}

#if (PY_VERSION_HEX < 0x030a00f0)
// this is in Python >=3.10
int PyModule_AddObjectRef(PyObject *module, const char *name, PyObject *value) {
    Py_INCREF(value);
    return PyModule_AddObject(module, name, value);
}
#endif

template <auto F>
struct ufunc_traits;

template <typename T>
struct npy_type;

template <>
struct npy_type<double> {
    static constexpr int value = NPY_FLOAT64;
};

template <>
struct npy_type<std::complex<double>> {
    static constexpr int value = NPY_COMPLEX128;
};

template <typename T>
struct npy_type;

template <typename Res, typename Arg0, Res (*F)(Arg0)>
struct ufunc_traits<F> {
    static constexpr int nin = 1;
    static constexpr int nout = 1;

    static constexpr std::array<char, 2> type = {npy_type<Arg0>::value, npy_type<Res>::value};

    static void func(char **args, const npy_intp *dimensions, const npy_intp *steps, void *data) {
        for (npy_intp i = 0; i < dimensions[0]; ++i) {
            *reinterpret_cast<Res *>(args[1]) = F(*reinterpret_cast<Arg0 *>(args[0]));

            args[0] += steps[0];
            args[1] += steps[1];
        }
    }
};

template <int N>
struct SpecFun_UFuncEntry {
    PyUFuncGenericFunction func[N];
    std::array<std::array<char, 2>, N> types;
};

// This function now generates a ufunc
template <auto... F>
PyObject *SpecFun_UFunc(const char *name, const char *doc) {
    using entry_type = SpecFun_UFuncEntry<sizeof...(F)>;

    static std::map<std::string, entry_type> m;

    auto &e = m[name];
    e = entry_type{{ufunc_traits<F>::func...}, {ufunc_traits<F>::type...}};

    return PyUFunc_FromFuncAndData(e.func, nullptr, reinterpret_cast<char *>(e.types.data()), sizeof...(F), 1, 1,
                                   PyUFunc_None, name, doc, 0);
}

static PyModuleDef specfun2_def = {
    PyModuleDef_HEAD_INIT,
    .m_name = "specfun2",
    .m_size = -1,
};

PyMODINIT_FUNC PyInit_specfun2() {
    SpecFun_Initialize();

    PyObject *specfun2 = PyModule_Create(&specfun2_def);
    if (specfun2 == nullptr) {
        return nullptr;
    }

    PyObject *expi = SpecFun_UFunc<special::expi, special::cexpi>("expi", "docstring goes here");
    PyModule_AddObjectRef(specfun2, "expi", expi);

    // More ufuncs would follow ...

    return specfun2;
}
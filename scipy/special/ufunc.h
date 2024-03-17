#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <numpy/arrayobject.h>
#include <numpy/npy_3kcompat.h>
#include <numpy/ufuncobject.h>

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

template <auto F>
struct ufunc_traits;

template <typename Res, typename Arg0, Res (*F)(Arg0)>
struct ufunc_traits<F> {
    static constexpr int nin = 1;
    static constexpr int nout = 1;

    static constexpr char type[2] = {npy_type<Arg0>::value, npy_type<Res>::value};

    static void func(char **args, const npy_intp *dimensions, const npy_intp *steps, void *data) {
        for (npy_intp i = 0; i < dimensions[0]; ++i) {
            *reinterpret_cast<Res *>(args[1]) = F(*reinterpret_cast<Arg0 *>(args[0]));

            args[0] += steps[0];
            args[1] += steps[1];
        }
    }
};

inline PyObject *SpecFun_UFunc(const char *name, const char *doc, std::initializer_list<PyUFuncGenericFunction> funcs,
                               std::initializer_list<const char *> types, int nin, int nout) {
    int ntypes = types.size();

    PyUFuncGenericFunction *funcs2 = new PyUFuncGenericFunction[ntypes];
    std::copy(funcs.begin(), funcs.end(), funcs2);

    char *types2 = new char[ntypes * (nin + nout)];
    for (auto it = types.begin(); it != types.end(); ++it) {
        memcpy(types2 + (it - types.begin()) * (nin + nout), *it, nin + nout);
    }

    return PyUFunc_FromFuncAndData(funcs2, nullptr, types2, ntypes, nin, nout, PyUFunc_None, name, doc, 0);
}

// This function now generates a ufunc
template <auto F0, auto... F>
PyObject *SpecFun_UFunc(const char *name, const char *doc) {
    return SpecFun_UFunc(name, doc, {ufunc_traits<F0>::func, ufunc_traits<F>::func...},
                         {ufunc_traits<F0>::type, ufunc_traits<F>::type...}, ufunc_traits<F0>::nin,
                         ufunc_traits<F0>::nout);
}
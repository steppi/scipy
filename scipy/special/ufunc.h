#include <utility>
#include <vector>

#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <numpy/arrayobject.h>
#include <numpy/npy_3kcompat.h>
#include <numpy/ufuncobject.h>

void sf_error_check_fpe(const char *func_name) {
    int status = PyUFunc_getfperr();
    if (status & UFUNC_FPE_DIVIDEBYZERO) {
        special::set_error(func_name, SF_ERROR_SINGULAR, "floating point division by zero");
    }
    if (status & UFUNC_FPE_UNDERFLOW) {
        special::set_error(func_name, SF_ERROR_UNDERFLOW, "floating point underflow");
    }
    if (status & UFUNC_FPE_OVERFLOW) {
        special::set_error(func_name, SF_ERROR_OVERFLOW, "floating point overflow");
    }
    if (status & UFUNC_FPE_INVALID) {
        special::set_error(func_name, SF_ERROR_DOMAIN, "floating point invalid value");
    }
}

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
struct arity_of;

template <typename Res, typename... Args, Res (*F)(Args...)>
struct arity_of<F> {
    static constexpr size_t value = sizeof...(Args);
};

template <auto F>
constexpr size_t arity_of_v = arity_of<F>::value;

template <auto F>
struct has_result;

template <typename Res, typename... Args, Res (*F)(Args...)>
struct has_result<F> {
    static constexpr bool value = true;
};

template <typename... Args, void (*F)(Args...)>
struct has_result<F> {
    static constexpr bool value = false;
};

template <auto F>
constexpr bool has_result_v = has_result<F>::value;

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

template <auto F, typename I = std::make_index_sequence<arity_of_v<F>>>
struct ufunc_traits;

template <typename Res, typename... Args, Res (*F)(Args...), size_t... I>
struct ufunc_traits<F, std::index_sequence<I...>> {
    static constexpr char type[sizeof...(Args) + 1] = {npy_type<Args>::value..., npy_type<Res>::value};

    static void func(char **args, const npy_intp *dimensions, const npy_intp *steps, void *data) {
        for (npy_intp i = 0; i < dimensions[0]; ++i) {
            *reinterpret_cast<Res *>(args[sizeof...(Args)]) = F(*reinterpret_cast<Args *>(args[I])...);

            for (npy_uintp j = 0; j < sizeof...(Args); ++j) {
                args[j] += steps[j];
            }
            args[sizeof...(Args)] += steps[sizeof...(Args)]; // output
        }

        const char *func_name = static_cast<char *>(data);
        sf_error_check_fpe(func_name);
    }
};

template <typename... Args, void (*F)(Args...), size_t... I>
struct ufunc_traits<F, std::index_sequence<I...>> {
    static constexpr char type[sizeof...(Args)] = {npy_type<Args>::value...};

    static void func(char **args, const npy_intp *dimensions, const npy_intp *steps, void *data) {
        for (npy_intp i = 0; i < dimensions[0]; ++i) {
            F(*reinterpret_cast<Args *>(args[I])...);

            for (npy_uintp j = 0; j < sizeof...(Args); ++j) {
                args[j] += steps[j];
            }
        }

        const char *func_name = static_cast<char *>(data);
        sf_error_check_fpe(func_name);
    }
};

template <int NTypes, int NInAndNOut>
struct SpecFun_UFuncFuncAndData {
    static constexpr int ntypes = NTypes;
    static constexpr int nin_and_nout = NInAndNOut;

    PyUFuncGenericFunction func[ntypes];
    char types[ntypes * nin_and_nout];
    void *data[ntypes];

    SpecFun_UFuncFuncAndData(std::initializer_list<PyUFuncGenericFunction> func,
                             std::initializer_list<const char *> types, const char *name) {
        std::copy(func.begin(), func.end(), this->func);

        for (auto it = types.begin(); it != types.end(); ++it) {
            std::copy(*it, *it + nin_and_nout, this->types + (it - types.begin()) * nin_and_nout);
        }

        std::fill_n(data, ntypes, const_cast<char *>(name));
    }
};

// This function now generates a ufunc
template <auto F0, auto... F>
PyObject *SpecFun_UFunc(const char *name, const char *doc, int nout) {
    static_assert(((arity_of_v<F0> == arity_of_v<F>) && ... && true),
                  "all functions must have the same number of arguments");
    static_assert(((has_result_v<F0> == has_result_v<F>) && ... && true),
                  "all functions must be void if any function is");

    using func_and_data_type = SpecFun_UFuncFuncAndData<sizeof...(F) + 1, arity_of_v<F0> + has_result_v<F0>>;

    static std::vector<func_and_data_type> entries;
    entries.push_back(
        {{ufunc_traits<F0>::func, ufunc_traits<F>::func...}, {ufunc_traits<F0>::type, ufunc_traits<F>::type...}, name});

    func_and_data_type &func_and_data = entries.back();
    return PyUFunc_FromFuncAndData(func_and_data.func, func_and_data.data, func_and_data.types,
                                   func_and_data_type::ntypes, func_and_data_type::nin_and_nout - nout, nout,
                                   PyUFunc_None, name, doc, 0);
}

template <auto F0, auto... F>
PyObject *SpecFun_UFunc(const char *name, const char *doc) {
    return SpecFun_UFunc<F0, F...>(name, doc, has_result_v<F0>);
}

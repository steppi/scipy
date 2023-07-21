#define PY_MAJOR_VERSION 3
#undef ENABLE_PYTHON_MODULE
#include <pythonic/core.hpp>
#include <pythonic/python/core.hpp>
#include <pythonic/types/bool.hpp>
#include <pythonic/types/int.hpp>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <pythonic/include/types/float64.hpp>
#include <pythonic/types/float64.hpp>
#include <pythonic/include/builtins/abs.hpp>
#include <pythonic/include/builtins/range.hpp>
#include <pythonic/include/builtins/round.hpp>
#include <pythonic/include/numpy/exp.hpp>
#include <pythonic/include/numpy/inf.hpp>
#include <pythonic/include/numpy/log.hpp>
#include <pythonic/include/numpy/square.hpp>
#include <pythonic/include/operator_/add.hpp>
#include <pythonic/include/operator_/div.hpp>
#include <pythonic/include/operator_/iadd.hpp>
#include <pythonic/include/operator_/lt.hpp>
#include <pythonic/include/operator_/mul.hpp>
#include <pythonic/include/operator_/neg.hpp>
#include <pythonic/builtins/abs.hpp>
#include <pythonic/builtins/range.hpp>
#include <pythonic/builtins/round.hpp>
#include <pythonic/numpy/exp.hpp>
#include <pythonic/numpy/inf.hpp>
#include <pythonic/numpy/log.hpp>
#include <pythonic/numpy/square.hpp>
#include <pythonic/operator_/add.hpp>
#include <pythonic/operator_/div.hpp>
#include <pythonic/operator_/iadd.hpp>
#include <pythonic/operator_/lt.hpp>
#include <pythonic/operator_/mul.hpp>
#include <pythonic/operator_/neg.hpp>
namespace 
{
  namespace __pythran__exp1
  {
    struct _exp1
    {
      typedef void callable;
      typedef void pure;
      template <typename argument_type0 >
      struct type
      {
        typedef double __type0;
        typedef typename pythonic::lazy<__type0>::type __type1;
        typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::log{})>::type>::type __type2;
        typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type3;
        typedef __type3 __type4;
        typedef decltype(std::declval<__type2>()(std::declval<__type4>())) __type5;
        typedef decltype(pythonic::operator_::mul(std::declval<__type0>(), std::declval<__type5>())) __type6;
        typedef typename pythonic::assignable<__type0>::type __type8;
        typedef __type8 __type9;
        typedef decltype(pythonic::operator_::neg(std::declval<__type9>())) __type10;
        typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::range{})>::type>::type __type11;
        typedef long __type12;
        typedef decltype(std::declval<__type11>()(std::declval<__type12>(), std::declval<__type12>())) __type13;
        typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type13>::type::iterator>::value_type>::type __type14;
        typedef __type14 __type15;
        typedef decltype(pythonic::operator_::mul(std::declval<__type10>(), std::declval<__type15>())) __type16;
        typedef decltype(pythonic::operator_::mul(std::declval<__type16>(), std::declval<__type4>())) __type18;
        typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::square{})>::type>::type __type19;
        typedef decltype(pythonic::operator_::add(std::declval<__type15>(), std::declval<__type12>())) __type21;
        typedef decltype(std::declval<__type19>()(std::declval<__type21>())) __type22;
        typedef decltype(pythonic::operator_::div(std::declval<__type18>(), std::declval<__type22>())) __type23;
        typedef typename pythonic::assignable<__type23>::type __type24;
        typedef typename __combined<__type8,__type24>::type __type25;
        typedef __type25 __type26;
        typedef typename __combined<__type8,__type26>::type __type27;
        typedef __type27 __type28;
        typedef decltype(pythonic::operator_::mul(std::declval<__type4>(), std::declval<__type28>())) __type29;
        typedef decltype(pythonic::operator_::add(std::declval<__type6>(), std::declval<__type29>())) __type30;
        typedef typename pythonic::lazy<__type30>::type __type31;
        typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::exp{})>::type>::type __type32;
        typedef decltype(pythonic::operator_::neg(std::declval<__type4>())) __type34;
        typedef decltype(std::declval<__type32>()(std::declval<__type34>())) __type35;
        typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::round{})>::type>::type __type37;
        typedef decltype(pythonic::operator_::div(std::declval<__type0>(), std::declval<__type4>())) __type39;
        typedef decltype(std::declval<__type37>()(std::declval<__type39>())) __type40;
        typedef decltype(pythonic::operator_::add(std::declval<__type12>(), std::declval<__type40>())) __type41;
        typedef typename pythonic::lazy<__type41>::type __type42;
        typedef __type42 __type43;
        typedef decltype(std::declval<__type11>()(std::declval<__type43>(), std::declval<__type12>(), std::declval<__type12>())) __type44;
        typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type44>::type::iterator>::value_type>::type __type45;
        typedef __type45 __type46;
        typedef __type1 __type49;
        typedef decltype(pythonic::operator_::add(std::declval<__type4>(), std::declval<__type49>())) __type50;
        typedef decltype(pythonic::operator_::div(std::declval<__type46>(), std::declval<__type50>())) __type51;
        typedef decltype(pythonic::operator_::add(std::declval<__type0>(), std::declval<__type51>())) __type52;
        typedef decltype(pythonic::operator_::div(std::declval<__type46>(), std::declval<__type52>())) __type53;
        typedef typename pythonic::lazy<__type53>::type __type54;
        typedef typename __combined<__type1,__type54>::type __type55;
        typedef __type55 __type56;
        typedef decltype(pythonic::operator_::add(std::declval<__type4>(), std::declval<__type56>())) __type57;
        typedef decltype(pythonic::operator_::div(std::declval<__type0>(), std::declval<__type57>())) __type58;
        typedef decltype(pythonic::operator_::mul(std::declval<__type35>(), std::declval<__type58>())) __type59;
        typedef typename pythonic::lazy<__type59>::type __type60;
        typedef typename __combined<__type1,__type31,__type60>::type __type61;
        typedef __type61 __type62;
        typedef typename pythonic::returnable<__type62>::type __type63;
        typedef __type63 result_type;
      }  
      ;
      template <typename argument_type0 >
      inline
      typename type<argument_type0>::result_type operator()(argument_type0 x) const
      ;
    }  ;
    template <typename argument_type0 >
    inline
    typename _exp1::type<argument_type0>::result_type _exp1::operator()(argument_type0 x) const
    {
      typedef double __type0;
      typedef typename pythonic::assignable<__type0>::type __type1;
      typedef __type1 __type2;
      typedef decltype(pythonic::operator_::neg(std::declval<__type2>())) __type3;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::range{})>::type>::type __type4;
      typedef long __type5;
      typedef decltype(std::declval<__type4>()(std::declval<__type5>(), std::declval<__type5>())) __type6;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type6>::type::iterator>::value_type>::type __type7;
      typedef __type7 __type8;
      typedef decltype(pythonic::operator_::mul(std::declval<__type3>(), std::declval<__type8>())) __type9;
      typedef typename std::remove_cv<typename std::remove_reference<argument_type0>::type>::type __type10;
      typedef __type10 __type11;
      typedef decltype(pythonic::operator_::mul(std::declval<__type9>(), std::declval<__type11>())) __type12;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::square{})>::type>::type __type13;
      typedef decltype(pythonic::operator_::add(std::declval<__type8>(), std::declval<__type5>())) __type15;
      typedef decltype(std::declval<__type13>()(std::declval<__type15>())) __type16;
      typedef decltype(pythonic::operator_::div(std::declval<__type12>(), std::declval<__type16>())) __type17;
      typedef typename pythonic::assignable<__type17>::type __type18;
      typedef typename __combined<__type1,__type18>::type __type19;
      typedef __type19 __type20;
      typedef typename __combined<__type1,__type20>::type __type21;
      typedef typename pythonic::assignable<__type21>::type __type22;
      typedef typename pythonic::assignable<__type19>::type __type23;
      typedef typename pythonic::lazy<__type0>::type __type24;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::builtins::functor::round{})>::type>::type __type25;
      typedef decltype(pythonic::operator_::div(std::declval<__type0>(), std::declval<__type11>())) __type27;
      typedef decltype(std::declval<__type25>()(std::declval<__type27>())) __type28;
      typedef decltype(pythonic::operator_::add(std::declval<__type5>(), std::declval<__type28>())) __type29;
      typedef typename pythonic::lazy<__type29>::type __type30;
      typedef __type30 __type31;
      typedef decltype(std::declval<__type4>()(std::declval<__type31>(), std::declval<__type5>(), std::declval<__type5>())) __type32;
      typedef typename std::remove_cv<typename std::iterator_traits<typename std::remove_reference<__type32>::type::iterator>::value_type>::type __type33;
      typedef __type33 __type34;
      typedef __type24 __type37;
      typedef decltype(pythonic::operator_::add(std::declval<__type11>(), std::declval<__type37>())) __type38;
      typedef decltype(pythonic::operator_::div(std::declval<__type34>(), std::declval<__type38>())) __type39;
      typedef decltype(pythonic::operator_::add(std::declval<__type0>(), std::declval<__type39>())) __type40;
      typedef decltype(pythonic::operator_::div(std::declval<__type34>(), std::declval<__type40>())) __type41;
      typedef typename pythonic::lazy<__type41>::type __type42;
      typedef typename __combined<__type24,__type42>::type __type43;
      typedef typename pythonic::lazy<__type43>::type __type44;
      typedef typename pythonic::lazy<__type30>::type __type45;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::log{})>::type>::type __type46;
      typedef decltype(std::declval<__type46>()(std::declval<__type11>())) __type48;
      typedef decltype(pythonic::operator_::mul(std::declval<__type0>(), std::declval<__type48>())) __type49;
      typedef __type21 __type51;
      typedef decltype(pythonic::operator_::mul(std::declval<__type11>(), std::declval<__type51>())) __type52;
      typedef decltype(pythonic::operator_::add(std::declval<__type49>(), std::declval<__type52>())) __type53;
      typedef typename pythonic::lazy<__type53>::type __type54;
      typedef typename std::remove_cv<typename std::remove_reference<decltype(pythonic::numpy::functor::exp{})>::type>::type __type55;
      typedef decltype(pythonic::operator_::neg(std::declval<__type11>())) __type57;
      typedef decltype(std::declval<__type55>()(std::declval<__type57>())) __type58;
      typedef __type43 __type60;
      typedef decltype(pythonic::operator_::add(std::declval<__type11>(), std::declval<__type60>())) __type61;
      typedef decltype(pythonic::operator_::div(std::declval<__type0>(), std::declval<__type61>())) __type62;
      typedef decltype(pythonic::operator_::mul(std::declval<__type58>(), std::declval<__type62>())) __type63;
      typedef typename pythonic::lazy<__type63>::type __type64;
      typedef typename __combined<__type24,__type54,__type64>::type __type65;
      typedef typename pythonic::lazy<__type65>::type __type66;
      __type66 e1;
      if (pythonic::operator_::lt(x, 0.0))
      {
        e1 = +pythonic::numpy::inf;
      }
      else
      {
        {
          __type44 t0;
          __type45 m;
          if (pythonic::operator_::lt(x, 1.0))
          {
            __type22 e1_ = 1.0;
            __type23 r = 1.0;
            {
              long  __target140436893154464 = 26L;
              for (long  k=1L; k < __target140436893154464; k += 1L)
              {
                r = pythonic::operator_::div(pythonic::operator_::mul(pythonic::operator_::mul(pythonic::operator_::neg(r), k), x), pythonic::numpy::functor::square{}(pythonic::operator_::add(k, 1L)));
                e1_ += r;
                if (pythonic::operator_::lt(pythonic::builtins::functor::abs{}(r), pythonic::operator_::mul(pythonic::builtins::functor::abs{}(e1_), 1e-15)))
                {
                  break;
                }
              }
            }
            e1 = pythonic::operator_::add(pythonic::operator_::mul(-0.5772156649015328, pythonic::numpy::functor::log{}(x)), pythonic::operator_::mul(x, e1_));
          }
          else
          {
            m = pythonic::operator_::add(20L, pythonic::builtins::functor::round{}(pythonic::operator_::div(80.0, x)));
            t0 = 0.0;
            {
              long  __target140436894009568 = 0L;
              for (long  k_=m; k_ > __target140436894009568; k_ += -1L)
              {
                t0 = pythonic::operator_::div(k_, pythonic::operator_::add(1.0, pythonic::operator_::div(k_, pythonic::operator_::add(x, t0))));
              }
            }
            e1 = pythonic::operator_::mul(pythonic::numpy::functor::exp{}(pythonic::operator_::neg(x)), pythonic::operator_::div(1.0, pythonic::operator_::add(x, t0)));
          }
        }
      }
      return e1;
    }
  }
}
#include <pythonic/python/exception_handler.hpp>
inline
typename __pythran__exp1::_exp1::type<double>::result_type _exp1(double x) 
{
  return __pythran__exp1::_exp1()(x);
}
#ifdef ENABLE_PYTHON_MODULE

static PyMethodDef Methods[] = {

    {NULL, NULL, 0, NULL}
};


#if PY_MAJOR_VERSION >= 3
  static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_exp1",            /* m_name */
    "",         /* m_doc */
    -1,                  /* m_size */
    Methods,             /* m_methods */
    NULL,                /* m_reload */
    NULL,                /* m_traverse */
    NULL,                /* m_clear */
    NULL,                /* m_free */
  };
#define PYTHRAN_RETURN return theModule
#define PYTHRAN_MODULE_INIT(s) PyInit_##s
#else
#define PYTHRAN_RETURN return
#define PYTHRAN_MODULE_INIT(s) init##s
#endif
PyMODINIT_FUNC
PYTHRAN_MODULE_INIT(_exp1)(void)
#ifndef _WIN32
__attribute__ ((visibility("default")))
#if defined(GNUC) && !defined(__clang__)
__attribute__ ((externally_visible))
#endif
#endif
;
PyMODINIT_FUNC
PYTHRAN_MODULE_INIT(_exp1)(void) {
    import_array()
    #if PY_MAJOR_VERSION >= 3
    PyObject* theModule = PyModule_Create(&moduledef);
    #else
    PyObject* theModule = Py_InitModule3("_exp1",
                                         Methods,
                                         ""
    );
    #endif
    if(! theModule)
        PYTHRAN_RETURN;
    PyObject * theDoc = Py_BuildValue("(sss)",
                                      "0.13.1",
                                      "2023-05-31 16:31:30.944558",
                                      "764f86d2d75091d442b06d4c987469b4aee1481ded370700fbaebc7bae9feda4");
    if(! theDoc)
        PYTHRAN_RETURN;
    PyModule_AddObject(theModule,
                       "__pythran__",
                       theDoc);


PyModule_AddObject(theModule, "_exp1",
                   PyCapsule_New((void*)&_exp1, "_exp1(float64)", NULL)
);
    PYTHRAN_RETURN;
}

#endif
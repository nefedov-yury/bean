#ifndef _Shared_DLLDEFINES_H_
#define _Shared_DLLDEFINES_H_

// This header defines BeanUserShared_EXPORT macro for use in
// a Windows system; in other operating systems it is empty

#if defined (_WIN32)
// CMake automatically adds <libname>_EXPORTS as a definition for
// source files compiled in a shared library called <libname>
// See about __declspec(dllexport / dllimport) here:
// https://learn.microsoft.com/en-us/cpp/build/
// + exporting-from-a-dll-using-declspec-dllexport
// + importing-into-an-application-using-declspec-dllimport
#if defined(BeanUser_EXPORTS)
#define  BeanUserShared_EXPORT __declspec(dllexport)
#else
#define  BeanUserShared_EXPORT __declspec(dllimport)
#endif
#else
#define BeanUserShared_EXPORT
#endif

#endif

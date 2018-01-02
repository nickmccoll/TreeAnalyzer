// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME UtilitiesdIHiggsCombineImportdIsrcdIdOdOdIobjdITreeAnalyzerDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "Utilities/HiggsCombineImport/src/../interface/Accumulators.h"
#include "Utilities/HiggsCombineImport/src/../interface/FastTemplate.h"
#include "Utilities/HiggsCombineImport/src/../interface/SimpleCacheSentry.h"
#include "Utilities/HiggsCombineImport/src/../interface/VerticalInterpHistPdf.h"
#include "Utilities/HiggsCombineImport/src/../interface/utils.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_SimpleCacheSentry(void *p = 0);
   static void *newArray_SimpleCacheSentry(Long_t size, void *p);
   static void delete_SimpleCacheSentry(void *p);
   static void deleteArray_SimpleCacheSentry(void *p);
   static void destruct_SimpleCacheSentry(void *p);
   static void streamer_SimpleCacheSentry(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::SimpleCacheSentry*)
   {
      ::SimpleCacheSentry *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::SimpleCacheSentry >(0);
      static ::ROOT::TGenericClassInfo 
         instance("SimpleCacheSentry", ::SimpleCacheSentry::Class_Version(), "Utilities/HiggsCombineImport/src/../interface/../interface/SimpleCacheSentry.h", 8,
                  typeid(::SimpleCacheSentry), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::SimpleCacheSentry::Dictionary, isa_proxy, 16,
                  sizeof(::SimpleCacheSentry) );
      instance.SetNew(&new_SimpleCacheSentry);
      instance.SetNewArray(&newArray_SimpleCacheSentry);
      instance.SetDelete(&delete_SimpleCacheSentry);
      instance.SetDeleteArray(&deleteArray_SimpleCacheSentry);
      instance.SetDestructor(&destruct_SimpleCacheSentry);
      instance.SetStreamerFunc(&streamer_SimpleCacheSentry);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::SimpleCacheSentry*)
   {
      return GenerateInitInstanceLocal((::SimpleCacheSentry*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::SimpleCacheSentry*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_VerticalInterpHistPdf(void *p = 0);
   static void *newArray_VerticalInterpHistPdf(Long_t size, void *p);
   static void delete_VerticalInterpHistPdf(void *p);
   static void deleteArray_VerticalInterpHistPdf(void *p);
   static void destruct_VerticalInterpHistPdf(void *p);
   static void streamer_VerticalInterpHistPdf(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::VerticalInterpHistPdf*)
   {
      ::VerticalInterpHistPdf *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::VerticalInterpHistPdf >(0);
      static ::ROOT::TGenericClassInfo 
         instance("VerticalInterpHistPdf", ::VerticalInterpHistPdf::Class_Version(), "Utilities/HiggsCombineImport/src/../interface/VerticalInterpHistPdf.h", 20,
                  typeid(::VerticalInterpHistPdf), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::VerticalInterpHistPdf::Dictionary, isa_proxy, 16,
                  sizeof(::VerticalInterpHistPdf) );
      instance.SetNew(&new_VerticalInterpHistPdf);
      instance.SetNewArray(&newArray_VerticalInterpHistPdf);
      instance.SetDelete(&delete_VerticalInterpHistPdf);
      instance.SetDeleteArray(&deleteArray_VerticalInterpHistPdf);
      instance.SetDestructor(&destruct_VerticalInterpHistPdf);
      instance.SetStreamerFunc(&streamer_VerticalInterpHistPdf);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::VerticalInterpHistPdf*)
   {
      return GenerateInitInstanceLocal((::VerticalInterpHistPdf*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::VerticalInterpHistPdf*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_FastVerticalInterpHistPdfBase(void *p);
   static void deleteArray_FastVerticalInterpHistPdfBase(void *p);
   static void destruct_FastVerticalInterpHistPdfBase(void *p);
   static void streamer_FastVerticalInterpHistPdfBase(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FastVerticalInterpHistPdfBase*)
   {
      ::FastVerticalInterpHistPdfBase *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::FastVerticalInterpHistPdfBase >(0);
      static ::ROOT::TGenericClassInfo 
         instance("FastVerticalInterpHistPdfBase", ::FastVerticalInterpHistPdfBase::Class_Version(), "Utilities/HiggsCombineImport/src/../interface/VerticalInterpHistPdf.h", 70,
                  typeid(::FastVerticalInterpHistPdfBase), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::FastVerticalInterpHistPdfBase::Dictionary, isa_proxy, 16,
                  sizeof(::FastVerticalInterpHistPdfBase) );
      instance.SetDelete(&delete_FastVerticalInterpHistPdfBase);
      instance.SetDeleteArray(&deleteArray_FastVerticalInterpHistPdfBase);
      instance.SetDestructor(&destruct_FastVerticalInterpHistPdfBase);
      instance.SetStreamerFunc(&streamer_FastVerticalInterpHistPdfBase);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FastVerticalInterpHistPdfBase*)
   {
      return GenerateInitInstanceLocal((::FastVerticalInterpHistPdfBase*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdfBase*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_FastVerticalInterpHistPdf(void *p = 0);
   static void *newArray_FastVerticalInterpHistPdf(Long_t size, void *p);
   static void delete_FastVerticalInterpHistPdf(void *p);
   static void deleteArray_FastVerticalInterpHistPdf(void *p);
   static void destruct_FastVerticalInterpHistPdf(void *p);
   static void streamer_FastVerticalInterpHistPdf(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FastVerticalInterpHistPdf*)
   {
      ::FastVerticalInterpHistPdf *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::FastVerticalInterpHistPdf >(0);
      static ::ROOT::TGenericClassInfo 
         instance("FastVerticalInterpHistPdf", ::FastVerticalInterpHistPdf::Class_Version(), "Utilities/HiggsCombineImport/src/../interface/VerticalInterpHistPdf.h", 128,
                  typeid(::FastVerticalInterpHistPdf), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::FastVerticalInterpHistPdf::Dictionary, isa_proxy, 16,
                  sizeof(::FastVerticalInterpHistPdf) );
      instance.SetNew(&new_FastVerticalInterpHistPdf);
      instance.SetNewArray(&newArray_FastVerticalInterpHistPdf);
      instance.SetDelete(&delete_FastVerticalInterpHistPdf);
      instance.SetDeleteArray(&deleteArray_FastVerticalInterpHistPdf);
      instance.SetDestructor(&destruct_FastVerticalInterpHistPdf);
      instance.SetStreamerFunc(&streamer_FastVerticalInterpHistPdf);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FastVerticalInterpHistPdf*)
   {
      return GenerateInitInstanceLocal((::FastVerticalInterpHistPdf*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdf*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_FastVerticalInterpHistPdf2D(void *p = 0);
   static void *newArray_FastVerticalInterpHistPdf2D(Long_t size, void *p);
   static void delete_FastVerticalInterpHistPdf2D(void *p);
   static void deleteArray_FastVerticalInterpHistPdf2D(void *p);
   static void destruct_FastVerticalInterpHistPdf2D(void *p);
   static void streamer_FastVerticalInterpHistPdf2D(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FastVerticalInterpHistPdf2D*)
   {
      ::FastVerticalInterpHistPdf2D *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::FastVerticalInterpHistPdf2D >(0);
      static ::ROOT::TGenericClassInfo 
         instance("FastVerticalInterpHistPdf2D", ::FastVerticalInterpHistPdf2D::Class_Version(), "Utilities/HiggsCombineImport/src/../interface/VerticalInterpHistPdf.h", 183,
                  typeid(::FastVerticalInterpHistPdf2D), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::FastVerticalInterpHistPdf2D::Dictionary, isa_proxy, 16,
                  sizeof(::FastVerticalInterpHistPdf2D) );
      instance.SetNew(&new_FastVerticalInterpHistPdf2D);
      instance.SetNewArray(&newArray_FastVerticalInterpHistPdf2D);
      instance.SetDelete(&delete_FastVerticalInterpHistPdf2D);
      instance.SetDeleteArray(&deleteArray_FastVerticalInterpHistPdf2D);
      instance.SetDestructor(&destruct_FastVerticalInterpHistPdf2D);
      instance.SetStreamerFunc(&streamer_FastVerticalInterpHistPdf2D);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FastVerticalInterpHistPdf2D*)
   {
      return GenerateInitInstanceLocal((::FastVerticalInterpHistPdf2D*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdf2D*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_FastVerticalInterpHistPdf2Base(void *p);
   static void deleteArray_FastVerticalInterpHistPdf2Base(void *p);
   static void destruct_FastVerticalInterpHistPdf2Base(void *p);
   static void streamer_FastVerticalInterpHistPdf2Base(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FastVerticalInterpHistPdf2Base*)
   {
      ::FastVerticalInterpHistPdf2Base *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::FastVerticalInterpHistPdf2Base >(0);
      static ::ROOT::TGenericClassInfo 
         instance("FastVerticalInterpHistPdf2Base", ::FastVerticalInterpHistPdf2Base::Class_Version(), "Utilities/HiggsCombineImport/src/../interface/VerticalInterpHistPdf.h", 233,
                  typeid(::FastVerticalInterpHistPdf2Base), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::FastVerticalInterpHistPdf2Base::Dictionary, isa_proxy, 16,
                  sizeof(::FastVerticalInterpHistPdf2Base) );
      instance.SetDelete(&delete_FastVerticalInterpHistPdf2Base);
      instance.SetDeleteArray(&deleteArray_FastVerticalInterpHistPdf2Base);
      instance.SetDestructor(&destruct_FastVerticalInterpHistPdf2Base);
      instance.SetStreamerFunc(&streamer_FastVerticalInterpHistPdf2Base);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FastVerticalInterpHistPdf2Base*)
   {
      return GenerateInitInstanceLocal((::FastVerticalInterpHistPdf2Base*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdf2Base*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_FastVerticalInterpHistPdf2(void *p = 0);
   static void *newArray_FastVerticalInterpHistPdf2(Long_t size, void *p);
   static void delete_FastVerticalInterpHistPdf2(void *p);
   static void deleteArray_FastVerticalInterpHistPdf2(void *p);
   static void destruct_FastVerticalInterpHistPdf2(void *p);
   static void streamer_FastVerticalInterpHistPdf2(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FastVerticalInterpHistPdf2*)
   {
      ::FastVerticalInterpHistPdf2 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::FastVerticalInterpHistPdf2 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("FastVerticalInterpHistPdf2", ::FastVerticalInterpHistPdf2::Class_Version(), "Utilities/HiggsCombineImport/src/../interface/VerticalInterpHistPdf.h", 297,
                  typeid(::FastVerticalInterpHistPdf2), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::FastVerticalInterpHistPdf2::Dictionary, isa_proxy, 16,
                  sizeof(::FastVerticalInterpHistPdf2) );
      instance.SetNew(&new_FastVerticalInterpHistPdf2);
      instance.SetNewArray(&newArray_FastVerticalInterpHistPdf2);
      instance.SetDelete(&delete_FastVerticalInterpHistPdf2);
      instance.SetDeleteArray(&deleteArray_FastVerticalInterpHistPdf2);
      instance.SetDestructor(&destruct_FastVerticalInterpHistPdf2);
      instance.SetStreamerFunc(&streamer_FastVerticalInterpHistPdf2);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FastVerticalInterpHistPdf2*)
   {
      return GenerateInitInstanceLocal((::FastVerticalInterpHistPdf2*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdf2*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_FastVerticalInterpHistPdf2D2(void *p = 0);
   static void *newArray_FastVerticalInterpHistPdf2D2(Long_t size, void *p);
   static void delete_FastVerticalInterpHistPdf2D2(void *p);
   static void deleteArray_FastVerticalInterpHistPdf2D2(void *p);
   static void destruct_FastVerticalInterpHistPdf2D2(void *p);
   static void streamer_FastVerticalInterpHistPdf2D2(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FastVerticalInterpHistPdf2D2*)
   {
      ::FastVerticalInterpHistPdf2D2 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::FastVerticalInterpHistPdf2D2 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("FastVerticalInterpHistPdf2D2", ::FastVerticalInterpHistPdf2D2::Class_Version(), "Utilities/HiggsCombineImport/src/../interface/VerticalInterpHistPdf.h", 351,
                  typeid(::FastVerticalInterpHistPdf2D2), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::FastVerticalInterpHistPdf2D2::Dictionary, isa_proxy, 16,
                  sizeof(::FastVerticalInterpHistPdf2D2) );
      instance.SetNew(&new_FastVerticalInterpHistPdf2D2);
      instance.SetNewArray(&newArray_FastVerticalInterpHistPdf2D2);
      instance.SetDelete(&delete_FastVerticalInterpHistPdf2D2);
      instance.SetDeleteArray(&deleteArray_FastVerticalInterpHistPdf2D2);
      instance.SetDestructor(&destruct_FastVerticalInterpHistPdf2D2);
      instance.SetStreamerFunc(&streamer_FastVerticalInterpHistPdf2D2);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FastVerticalInterpHistPdf2D2*)
   {
      return GenerateInitInstanceLocal((::FastVerticalInterpHistPdf2D2*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdf2D2*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_FastVerticalInterpHistPdf3D(void *p = 0);
   static void *newArray_FastVerticalInterpHistPdf3D(Long_t size, void *p);
   static void delete_FastVerticalInterpHistPdf3D(void *p);
   static void deleteArray_FastVerticalInterpHistPdf3D(void *p);
   static void destruct_FastVerticalInterpHistPdf3D(void *p);
   static void streamer_FastVerticalInterpHistPdf3D(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::FastVerticalInterpHistPdf3D*)
   {
      ::FastVerticalInterpHistPdf3D *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::FastVerticalInterpHistPdf3D >(0);
      static ::ROOT::TGenericClassInfo 
         instance("FastVerticalInterpHistPdf3D", ::FastVerticalInterpHistPdf3D::Class_Version(), "Utilities/HiggsCombineImport/src/../interface/VerticalInterpHistPdf.h", 393,
                  typeid(::FastVerticalInterpHistPdf3D), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::FastVerticalInterpHistPdf3D::Dictionary, isa_proxy, 16,
                  sizeof(::FastVerticalInterpHistPdf3D) );
      instance.SetNew(&new_FastVerticalInterpHistPdf3D);
      instance.SetNewArray(&newArray_FastVerticalInterpHistPdf3D);
      instance.SetDelete(&delete_FastVerticalInterpHistPdf3D);
      instance.SetDeleteArray(&deleteArray_FastVerticalInterpHistPdf3D);
      instance.SetDestructor(&destruct_FastVerticalInterpHistPdf3D);
      instance.SetStreamerFunc(&streamer_FastVerticalInterpHistPdf3D);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::FastVerticalInterpHistPdf3D*)
   {
      return GenerateInitInstanceLocal((::FastVerticalInterpHistPdf3D*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdf3D*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr SimpleCacheSentry::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *SimpleCacheSentry::Class_Name()
{
   return "SimpleCacheSentry";
}

//______________________________________________________________________________
const char *SimpleCacheSentry::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SimpleCacheSentry*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int SimpleCacheSentry::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::SimpleCacheSentry*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *SimpleCacheSentry::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SimpleCacheSentry*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *SimpleCacheSentry::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::SimpleCacheSentry*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr VerticalInterpHistPdf::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *VerticalInterpHistPdf::Class_Name()
{
   return "VerticalInterpHistPdf";
}

//______________________________________________________________________________
const char *VerticalInterpHistPdf::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::VerticalInterpHistPdf*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int VerticalInterpHistPdf::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::VerticalInterpHistPdf*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *VerticalInterpHistPdf::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::VerticalInterpHistPdf*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *VerticalInterpHistPdf::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::VerticalInterpHistPdf*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr FastVerticalInterpHistPdfBase::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *FastVerticalInterpHistPdfBase::Class_Name()
{
   return "FastVerticalInterpHistPdfBase";
}

//______________________________________________________________________________
const char *FastVerticalInterpHistPdfBase::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdfBase*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int FastVerticalInterpHistPdfBase::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdfBase*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *FastVerticalInterpHistPdfBase::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdfBase*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *FastVerticalInterpHistPdfBase::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdfBase*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr FastVerticalInterpHistPdf::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *FastVerticalInterpHistPdf::Class_Name()
{
   return "FastVerticalInterpHistPdf";
}

//______________________________________________________________________________
const char *FastVerticalInterpHistPdf::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdf*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int FastVerticalInterpHistPdf::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdf*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *FastVerticalInterpHistPdf::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdf*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *FastVerticalInterpHistPdf::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdf*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr FastVerticalInterpHistPdf2D::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *FastVerticalInterpHistPdf2D::Class_Name()
{
   return "FastVerticalInterpHistPdf2D";
}

//______________________________________________________________________________
const char *FastVerticalInterpHistPdf2D::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdf2D*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int FastVerticalInterpHistPdf2D::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdf2D*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *FastVerticalInterpHistPdf2D::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdf2D*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *FastVerticalInterpHistPdf2D::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdf2D*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr FastVerticalInterpHistPdf2Base::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *FastVerticalInterpHistPdf2Base::Class_Name()
{
   return "FastVerticalInterpHistPdf2Base";
}

//______________________________________________________________________________
const char *FastVerticalInterpHistPdf2Base::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdf2Base*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int FastVerticalInterpHistPdf2Base::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdf2Base*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *FastVerticalInterpHistPdf2Base::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdf2Base*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *FastVerticalInterpHistPdf2Base::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdf2Base*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr FastVerticalInterpHistPdf2::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *FastVerticalInterpHistPdf2::Class_Name()
{
   return "FastVerticalInterpHistPdf2";
}

//______________________________________________________________________________
const char *FastVerticalInterpHistPdf2::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdf2*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int FastVerticalInterpHistPdf2::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdf2*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *FastVerticalInterpHistPdf2::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdf2*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *FastVerticalInterpHistPdf2::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdf2*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr FastVerticalInterpHistPdf2D2::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *FastVerticalInterpHistPdf2D2::Class_Name()
{
   return "FastVerticalInterpHistPdf2D2";
}

//______________________________________________________________________________
const char *FastVerticalInterpHistPdf2D2::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdf2D2*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int FastVerticalInterpHistPdf2D2::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdf2D2*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *FastVerticalInterpHistPdf2D2::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdf2D2*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *FastVerticalInterpHistPdf2D2::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdf2D2*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr FastVerticalInterpHistPdf3D::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *FastVerticalInterpHistPdf3D::Class_Name()
{
   return "FastVerticalInterpHistPdf3D";
}

//______________________________________________________________________________
const char *FastVerticalInterpHistPdf3D::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdf3D*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int FastVerticalInterpHistPdf3D::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdf3D*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *FastVerticalInterpHistPdf3D::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdf3D*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *FastVerticalInterpHistPdf3D::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::FastVerticalInterpHistPdf3D*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void SimpleCacheSentry::Streamer(TBuffer &R__b)
{
   // Stream an object of class SimpleCacheSentry.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      RooAbsArg::Streamer(R__b);
      _deps.Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, SimpleCacheSentry::IsA());
   } else {
      R__c = R__b.WriteVersion(SimpleCacheSentry::IsA(), kTRUE);
      RooAbsArg::Streamer(R__b);
      _deps.Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_SimpleCacheSentry(void *p) {
      return  p ? new(p) ::SimpleCacheSentry : new ::SimpleCacheSentry;
   }
   static void *newArray_SimpleCacheSentry(Long_t nElements, void *p) {
      return p ? new(p) ::SimpleCacheSentry[nElements] : new ::SimpleCacheSentry[nElements];
   }
   // Wrapper around operator delete
   static void delete_SimpleCacheSentry(void *p) {
      delete ((::SimpleCacheSentry*)p);
   }
   static void deleteArray_SimpleCacheSentry(void *p) {
      delete [] ((::SimpleCacheSentry*)p);
   }
   static void destruct_SimpleCacheSentry(void *p) {
      typedef ::SimpleCacheSentry current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_SimpleCacheSentry(TBuffer &buf, void *obj) {
      ((::SimpleCacheSentry*)obj)->::SimpleCacheSentry::Streamer(buf);
   }
} // end of namespace ROOT for class ::SimpleCacheSentry

//______________________________________________________________________________
void VerticalInterpHistPdf::Streamer(TBuffer &R__b)
{
   // Stream an object of class VerticalInterpHistPdf.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      RooAbsPdf::Streamer(R__b);
      _x.Streamer(R__b);
      _funcList.Streamer(R__b);
      _coefList.Streamer(R__b);
      R__b >> _smoothRegion;
      R__b >> _smoothAlgo;
      R__b.CheckByteCount(R__s, R__c, VerticalInterpHistPdf::IsA());
   } else {
      R__c = R__b.WriteVersion(VerticalInterpHistPdf::IsA(), kTRUE);
      RooAbsPdf::Streamer(R__b);
      _x.Streamer(R__b);
      _funcList.Streamer(R__b);
      _coefList.Streamer(R__b);
      R__b << _smoothRegion;
      R__b << _smoothAlgo;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_VerticalInterpHistPdf(void *p) {
      return  p ? new(p) ::VerticalInterpHistPdf : new ::VerticalInterpHistPdf;
   }
   static void *newArray_VerticalInterpHistPdf(Long_t nElements, void *p) {
      return p ? new(p) ::VerticalInterpHistPdf[nElements] : new ::VerticalInterpHistPdf[nElements];
   }
   // Wrapper around operator delete
   static void delete_VerticalInterpHistPdf(void *p) {
      delete ((::VerticalInterpHistPdf*)p);
   }
   static void deleteArray_VerticalInterpHistPdf(void *p) {
      delete [] ((::VerticalInterpHistPdf*)p);
   }
   static void destruct_VerticalInterpHistPdf(void *p) {
      typedef ::VerticalInterpHistPdf current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_VerticalInterpHistPdf(TBuffer &buf, void *obj) {
      ((::VerticalInterpHistPdf*)obj)->::VerticalInterpHistPdf::Streamer(buf);
   }
} // end of namespace ROOT for class ::VerticalInterpHistPdf

//______________________________________________________________________________
void FastVerticalInterpHistPdfBase::Streamer(TBuffer &R__b)
{
   // Stream an object of class FastVerticalInterpHistPdfBase.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      RooAbsPdf::Streamer(R__b);
      _x.Streamer(R__b);
      _funcList.Streamer(R__b);
      _coefList.Streamer(R__b);
      R__b >> _smoothRegion;
      R__b >> _smoothAlgo;
      R__b.CheckByteCount(R__s, R__c, FastVerticalInterpHistPdfBase::IsA());
   } else {
      R__c = R__b.WriteVersion(FastVerticalInterpHistPdfBase::IsA(), kTRUE);
      RooAbsPdf::Streamer(R__b);
      _x.Streamer(R__b);
      _funcList.Streamer(R__b);
      _coefList.Streamer(R__b);
      R__b << _smoothRegion;
      R__b << _smoothAlgo;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_FastVerticalInterpHistPdfBase(void *p) {
      delete ((::FastVerticalInterpHistPdfBase*)p);
   }
   static void deleteArray_FastVerticalInterpHistPdfBase(void *p) {
      delete [] ((::FastVerticalInterpHistPdfBase*)p);
   }
   static void destruct_FastVerticalInterpHistPdfBase(void *p) {
      typedef ::FastVerticalInterpHistPdfBase current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_FastVerticalInterpHistPdfBase(TBuffer &buf, void *obj) {
      ((::FastVerticalInterpHistPdfBase*)obj)->::FastVerticalInterpHistPdfBase::Streamer(buf);
   }
} // end of namespace ROOT for class ::FastVerticalInterpHistPdfBase

//______________________________________________________________________________
void FastVerticalInterpHistPdf::Streamer(TBuffer &R__b)
{
   // Stream an object of class FastVerticalInterpHistPdf.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      FastVerticalInterpHistPdfBase::Streamer(R__b);
      _x.Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, FastVerticalInterpHistPdf::IsA());
   } else {
      R__c = R__b.WriteVersion(FastVerticalInterpHistPdf::IsA(), kTRUE);
      FastVerticalInterpHistPdfBase::Streamer(R__b);
      _x.Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_FastVerticalInterpHistPdf(void *p) {
      return  p ? new(p) ::FastVerticalInterpHistPdf : new ::FastVerticalInterpHistPdf;
   }
   static void *newArray_FastVerticalInterpHistPdf(Long_t nElements, void *p) {
      return p ? new(p) ::FastVerticalInterpHistPdf[nElements] : new ::FastVerticalInterpHistPdf[nElements];
   }
   // Wrapper around operator delete
   static void delete_FastVerticalInterpHistPdf(void *p) {
      delete ((::FastVerticalInterpHistPdf*)p);
   }
   static void deleteArray_FastVerticalInterpHistPdf(void *p) {
      delete [] ((::FastVerticalInterpHistPdf*)p);
   }
   static void destruct_FastVerticalInterpHistPdf(void *p) {
      typedef ::FastVerticalInterpHistPdf current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_FastVerticalInterpHistPdf(TBuffer &buf, void *obj) {
      ((::FastVerticalInterpHistPdf*)obj)->::FastVerticalInterpHistPdf::Streamer(buf);
   }
} // end of namespace ROOT for class ::FastVerticalInterpHistPdf

//______________________________________________________________________________
void FastVerticalInterpHistPdf2D::Streamer(TBuffer &R__b)
{
   // Stream an object of class FastVerticalInterpHistPdf2D.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      FastVerticalInterpHistPdfBase::Streamer(R__b);
      _x.Streamer(R__b);
      _y.Streamer(R__b);
      R__b >> _conditional;
      R__b.CheckByteCount(R__s, R__c, FastVerticalInterpHistPdf2D::IsA());
   } else {
      R__c = R__b.WriteVersion(FastVerticalInterpHistPdf2D::IsA(), kTRUE);
      FastVerticalInterpHistPdfBase::Streamer(R__b);
      _x.Streamer(R__b);
      _y.Streamer(R__b);
      R__b << _conditional;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_FastVerticalInterpHistPdf2D(void *p) {
      return  p ? new(p) ::FastVerticalInterpHistPdf2D : new ::FastVerticalInterpHistPdf2D;
   }
   static void *newArray_FastVerticalInterpHistPdf2D(Long_t nElements, void *p) {
      return p ? new(p) ::FastVerticalInterpHistPdf2D[nElements] : new ::FastVerticalInterpHistPdf2D[nElements];
   }
   // Wrapper around operator delete
   static void delete_FastVerticalInterpHistPdf2D(void *p) {
      delete ((::FastVerticalInterpHistPdf2D*)p);
   }
   static void deleteArray_FastVerticalInterpHistPdf2D(void *p) {
      delete [] ((::FastVerticalInterpHistPdf2D*)p);
   }
   static void destruct_FastVerticalInterpHistPdf2D(void *p) {
      typedef ::FastVerticalInterpHistPdf2D current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_FastVerticalInterpHistPdf2D(TBuffer &buf, void *obj) {
      ((::FastVerticalInterpHistPdf2D*)obj)->::FastVerticalInterpHistPdf2D::Streamer(buf);
   }
} // end of namespace ROOT for class ::FastVerticalInterpHistPdf2D

//______________________________________________________________________________
void FastVerticalInterpHistPdf2Base::Streamer(TBuffer &R__b)
{
   // Stream an object of class FastVerticalInterpHistPdf2Base.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      RooAbsPdf::Streamer(R__b);
      _coefList.Streamer(R__b);
      R__b >> _smoothRegion;
      R__b >> _smoothAlgo;
      {
         vector<FastVerticalInterpHistPdf2Base::Morph> &R__stl =  _morphs;
         R__stl.clear();
         TClass *R__tcl1 = TBuffer::GetClass(typeid(struct FastVerticalInterpHistPdfBase::Morph));
         if (R__tcl1==0) {
            Error("_morphs streamer","Missing the TClass object for struct FastVerticalInterpHistPdfBase::Morph!");
            return;
         }
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            FastVerticalInterpHistPdfBase::Morph R__t;
            R__b.StreamObject(&R__t,R__tcl1);
            R__stl.push_back(R__t);
         }
      }
      R__b.CheckByteCount(R__s, R__c, FastVerticalInterpHistPdf2Base::IsA());
   } else {
      R__c = R__b.WriteVersion(FastVerticalInterpHistPdf2Base::IsA(), kTRUE);
      RooAbsPdf::Streamer(R__b);
      _coefList.Streamer(R__b);
      R__b << _smoothRegion;
      R__b << _smoothAlgo;
      {
         vector<FastVerticalInterpHistPdf2Base::Morph> &R__stl =  _morphs;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
         TClass *R__tcl1 = TBuffer::GetClass(typeid(struct FastVerticalInterpHistPdfBase::Morph));
         if (R__tcl1==0) {
            Error("_morphs streamer","Missing the TClass object for struct FastVerticalInterpHistPdfBase::Morph!");
            return;
         }
            vector<FastVerticalInterpHistPdf2Base::Morph>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b.StreamObject((FastVerticalInterpHistPdfBase::Morph*)&(*R__k),R__tcl1);
            }
         }
      }
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_FastVerticalInterpHistPdf2Base(void *p) {
      delete ((::FastVerticalInterpHistPdf2Base*)p);
   }
   static void deleteArray_FastVerticalInterpHistPdf2Base(void *p) {
      delete [] ((::FastVerticalInterpHistPdf2Base*)p);
   }
   static void destruct_FastVerticalInterpHistPdf2Base(void *p) {
      typedef ::FastVerticalInterpHistPdf2Base current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_FastVerticalInterpHistPdf2Base(TBuffer &buf, void *obj) {
      ((::FastVerticalInterpHistPdf2Base*)obj)->::FastVerticalInterpHistPdf2Base::Streamer(buf);
   }
} // end of namespace ROOT for class ::FastVerticalInterpHistPdf2Base

//______________________________________________________________________________
void FastVerticalInterpHistPdf2::Streamer(TBuffer &R__b)
{
   // Stream an object of class FastVerticalInterpHistPdf2.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      FastVerticalInterpHistPdf2Base::Streamer(R__b);
      _x.Streamer(R__b);
      R__b.StreamObject(&(_cacheNominal),typeid(_cacheNominal));
      R__b.StreamObject(&(_cacheNominalLog),typeid(_cacheNominalLog));
      R__b.CheckByteCount(R__s, R__c, FastVerticalInterpHistPdf2::IsA());
   } else {
      R__c = R__b.WriteVersion(FastVerticalInterpHistPdf2::IsA(), kTRUE);
      FastVerticalInterpHistPdf2Base::Streamer(R__b);
      _x.Streamer(R__b);
      R__b.StreamObject(&(_cacheNominal),typeid(_cacheNominal));
      R__b.StreamObject(&(_cacheNominalLog),typeid(_cacheNominalLog));
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_FastVerticalInterpHistPdf2(void *p) {
      return  p ? new(p) ::FastVerticalInterpHistPdf2 : new ::FastVerticalInterpHistPdf2;
   }
   static void *newArray_FastVerticalInterpHistPdf2(Long_t nElements, void *p) {
      return p ? new(p) ::FastVerticalInterpHistPdf2[nElements] : new ::FastVerticalInterpHistPdf2[nElements];
   }
   // Wrapper around operator delete
   static void delete_FastVerticalInterpHistPdf2(void *p) {
      delete ((::FastVerticalInterpHistPdf2*)p);
   }
   static void deleteArray_FastVerticalInterpHistPdf2(void *p) {
      delete [] ((::FastVerticalInterpHistPdf2*)p);
   }
   static void destruct_FastVerticalInterpHistPdf2(void *p) {
      typedef ::FastVerticalInterpHistPdf2 current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_FastVerticalInterpHistPdf2(TBuffer &buf, void *obj) {
      ((::FastVerticalInterpHistPdf2*)obj)->::FastVerticalInterpHistPdf2::Streamer(buf);
   }
} // end of namespace ROOT for class ::FastVerticalInterpHistPdf2

//______________________________________________________________________________
void FastVerticalInterpHistPdf2D2::Streamer(TBuffer &R__b)
{
   // Stream an object of class FastVerticalInterpHistPdf2D2.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      FastVerticalInterpHistPdf2Base::Streamer(R__b);
      _x.Streamer(R__b);
      _y.Streamer(R__b);
      R__b >> _conditional;
      R__b.StreamObject(&(_cacheNominal),typeid(_cacheNominal));
      R__b.StreamObject(&(_cacheNominalLog),typeid(_cacheNominalLog));
      R__b.CheckByteCount(R__s, R__c, FastVerticalInterpHistPdf2D2::IsA());
   } else {
      R__c = R__b.WriteVersion(FastVerticalInterpHistPdf2D2::IsA(), kTRUE);
      FastVerticalInterpHistPdf2Base::Streamer(R__b);
      _x.Streamer(R__b);
      _y.Streamer(R__b);
      R__b << _conditional;
      R__b.StreamObject(&(_cacheNominal),typeid(_cacheNominal));
      R__b.StreamObject(&(_cacheNominalLog),typeid(_cacheNominalLog));
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_FastVerticalInterpHistPdf2D2(void *p) {
      return  p ? new(p) ::FastVerticalInterpHistPdf2D2 : new ::FastVerticalInterpHistPdf2D2;
   }
   static void *newArray_FastVerticalInterpHistPdf2D2(Long_t nElements, void *p) {
      return p ? new(p) ::FastVerticalInterpHistPdf2D2[nElements] : new ::FastVerticalInterpHistPdf2D2[nElements];
   }
   // Wrapper around operator delete
   static void delete_FastVerticalInterpHistPdf2D2(void *p) {
      delete ((::FastVerticalInterpHistPdf2D2*)p);
   }
   static void deleteArray_FastVerticalInterpHistPdf2D2(void *p) {
      delete [] ((::FastVerticalInterpHistPdf2D2*)p);
   }
   static void destruct_FastVerticalInterpHistPdf2D2(void *p) {
      typedef ::FastVerticalInterpHistPdf2D2 current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_FastVerticalInterpHistPdf2D2(TBuffer &buf, void *obj) {
      ((::FastVerticalInterpHistPdf2D2*)obj)->::FastVerticalInterpHistPdf2D2::Streamer(buf);
   }
} // end of namespace ROOT for class ::FastVerticalInterpHistPdf2D2

//______________________________________________________________________________
void FastVerticalInterpHistPdf3D::Streamer(TBuffer &R__b)
{
   // Stream an object of class FastVerticalInterpHistPdf3D.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      FastVerticalInterpHistPdfBase::Streamer(R__b);
      _x.Streamer(R__b);
      _y.Streamer(R__b);
      _z.Streamer(R__b);
      R__b >> _conditional;
      R__b.CheckByteCount(R__s, R__c, FastVerticalInterpHistPdf3D::IsA());
   } else {
      R__c = R__b.WriteVersion(FastVerticalInterpHistPdf3D::IsA(), kTRUE);
      FastVerticalInterpHistPdfBase::Streamer(R__b);
      _x.Streamer(R__b);
      _y.Streamer(R__b);
      _z.Streamer(R__b);
      R__b << _conditional;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_FastVerticalInterpHistPdf3D(void *p) {
      return  p ? new(p) ::FastVerticalInterpHistPdf3D : new ::FastVerticalInterpHistPdf3D;
   }
   static void *newArray_FastVerticalInterpHistPdf3D(Long_t nElements, void *p) {
      return p ? new(p) ::FastVerticalInterpHistPdf3D[nElements] : new ::FastVerticalInterpHistPdf3D[nElements];
   }
   // Wrapper around operator delete
   static void delete_FastVerticalInterpHistPdf3D(void *p) {
      delete ((::FastVerticalInterpHistPdf3D*)p);
   }
   static void deleteArray_FastVerticalInterpHistPdf3D(void *p) {
      delete [] ((::FastVerticalInterpHistPdf3D*)p);
   }
   static void destruct_FastVerticalInterpHistPdf3D(void *p) {
      typedef ::FastVerticalInterpHistPdf3D current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_FastVerticalInterpHistPdf3D(TBuffer &buf, void *obj) {
      ((::FastVerticalInterpHistPdf3D*)obj)->::FastVerticalInterpHistPdf3D::Streamer(buf);
   }
} // end of namespace ROOT for class ::FastVerticalInterpHistPdf3D

namespace ROOT {
   static TClass *vectorlEFastVerticalInterpHistPdfBasecLcLMorphgR_Dictionary();
   static void vectorlEFastVerticalInterpHistPdfBasecLcLMorphgR_TClassManip(TClass*);
   static void *new_vectorlEFastVerticalInterpHistPdfBasecLcLMorphgR(void *p = 0);
   static void *newArray_vectorlEFastVerticalInterpHistPdfBasecLcLMorphgR(Long_t size, void *p);
   static void delete_vectorlEFastVerticalInterpHistPdfBasecLcLMorphgR(void *p);
   static void deleteArray_vectorlEFastVerticalInterpHistPdfBasecLcLMorphgR(void *p);
   static void destruct_vectorlEFastVerticalInterpHistPdfBasecLcLMorphgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<FastVerticalInterpHistPdfBase::Morph>*)
   {
      vector<FastVerticalInterpHistPdfBase::Morph> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<FastVerticalInterpHistPdfBase::Morph>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<FastVerticalInterpHistPdfBase::Morph>", -2, "vector", 450,
                  typeid(vector<FastVerticalInterpHistPdfBase::Morph>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEFastVerticalInterpHistPdfBasecLcLMorphgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<FastVerticalInterpHistPdfBase::Morph>) );
      instance.SetNew(&new_vectorlEFastVerticalInterpHistPdfBasecLcLMorphgR);
      instance.SetNewArray(&newArray_vectorlEFastVerticalInterpHistPdfBasecLcLMorphgR);
      instance.SetDelete(&delete_vectorlEFastVerticalInterpHistPdfBasecLcLMorphgR);
      instance.SetDeleteArray(&deleteArray_vectorlEFastVerticalInterpHistPdfBasecLcLMorphgR);
      instance.SetDestructor(&destruct_vectorlEFastVerticalInterpHistPdfBasecLcLMorphgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<FastVerticalInterpHistPdfBase::Morph> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<FastVerticalInterpHistPdfBase::Morph>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEFastVerticalInterpHistPdfBasecLcLMorphgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<FastVerticalInterpHistPdfBase::Morph>*)0x0)->GetClass();
      vectorlEFastVerticalInterpHistPdfBasecLcLMorphgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEFastVerticalInterpHistPdfBasecLcLMorphgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEFastVerticalInterpHistPdfBasecLcLMorphgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<FastVerticalInterpHistPdfBase::Morph> : new vector<FastVerticalInterpHistPdfBase::Morph>;
   }
   static void *newArray_vectorlEFastVerticalInterpHistPdfBasecLcLMorphgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<FastVerticalInterpHistPdfBase::Morph>[nElements] : new vector<FastVerticalInterpHistPdfBase::Morph>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEFastVerticalInterpHistPdfBasecLcLMorphgR(void *p) {
      delete ((vector<FastVerticalInterpHistPdfBase::Morph>*)p);
   }
   static void deleteArray_vectorlEFastVerticalInterpHistPdfBasecLcLMorphgR(void *p) {
      delete [] ((vector<FastVerticalInterpHistPdfBase::Morph>*)p);
   }
   static void destruct_vectorlEFastVerticalInterpHistPdfBasecLcLMorphgR(void *p) {
      typedef vector<FastVerticalInterpHistPdfBase::Morph> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<FastVerticalInterpHistPdfBase::Morph>

namespace {
  void TriggerDictionaryInitialization_TreeAnalyzerDict_Impl() {
    static const char* headers[] = {
"Utilities/HiggsCombineImport/src/../interface/Accumulators.h",
"Utilities/HiggsCombineImport/src/../interface/FastTemplate.h",
"Utilities/HiggsCombineImport/src/../interface/SimpleCacheSentry.h",
"Utilities/HiggsCombineImport/src/../interface/VerticalInterpHistPdf.h",
"Utilities/HiggsCombineImport/src/../interface/utils.h",
0
    };
    static const char* includePaths[] = {
"../..",
"/Users/nmccoll/Dropbox/Work/programs/root/include",
"/Users/nmccoll/Dropbox/Work/GitRepositories/TreeAnalyzer/framework/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "TreeAnalyzerDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$Utilities/HiggsCombineImport/src/../interface/SimpleCacheSentry.h")))  SimpleCacheSentry;
class __attribute__((annotate("$clingAutoload$Utilities/HiggsCombineImport/src/../interface/VerticalInterpHistPdf.h")))  VerticalInterpHistPdf;
class __attribute__((annotate("$clingAutoload$Utilities/HiggsCombineImport/src/../interface/VerticalInterpHistPdf.h")))  FastVerticalInterpHistPdfBase;
class __attribute__((annotate("$clingAutoload$Utilities/HiggsCombineImport/src/../interface/VerticalInterpHistPdf.h")))  FastVerticalInterpHistPdf;
class __attribute__((annotate("$clingAutoload$Utilities/HiggsCombineImport/src/../interface/VerticalInterpHistPdf.h")))  FastVerticalInterpHistPdf2D;
class __attribute__((annotate("$clingAutoload$Utilities/HiggsCombineImport/src/../interface/VerticalInterpHistPdf.h")))  FastVerticalInterpHistPdf2Base;
class __attribute__((annotate("$clingAutoload$Utilities/HiggsCombineImport/src/../interface/VerticalInterpHistPdf.h")))  FastVerticalInterpHistPdf2;
class __attribute__((annotate("$clingAutoload$Utilities/HiggsCombineImport/src/../interface/VerticalInterpHistPdf.h")))  FastVerticalInterpHistPdf2D2;
class __attribute__((annotate("$clingAutoload$Utilities/HiggsCombineImport/src/../interface/VerticalInterpHistPdf.h")))  FastVerticalInterpHistPdf3D;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "TreeAnalyzerDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "Utilities/HiggsCombineImport/src/../interface/Accumulators.h"
#include "Utilities/HiggsCombineImport/src/../interface/FastTemplate.h"
#include "Utilities/HiggsCombineImport/src/../interface/SimpleCacheSentry.h"
#include "Utilities/HiggsCombineImport/src/../interface/VerticalInterpHistPdf.h"
#include "Utilities/HiggsCombineImport/src/../interface/utils.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"FastVerticalInterpHistPdf", payloadCode, "@",
"FastVerticalInterpHistPdf2", payloadCode, "@",
"FastVerticalInterpHistPdf2Base", payloadCode, "@",
"FastVerticalInterpHistPdf2D", payloadCode, "@",
"FastVerticalInterpHistPdf2D2", payloadCode, "@",
"FastVerticalInterpHistPdf3D", payloadCode, "@",
"FastVerticalInterpHistPdfBase", payloadCode, "@",
"SimpleCacheSentry", payloadCode, "@",
"VerticalInterpHistPdf", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("TreeAnalyzerDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_TreeAnalyzerDict_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_TreeAnalyzerDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_TreeAnalyzerDict() {
  TriggerDictionaryInitialization_TreeAnalyzerDict_Impl();
}

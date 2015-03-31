#ifndef __ConvertBuffer_h
#define __ConvertBuffer_h

#include "TypeInfo.h"

#include <cassert>

template <typename T>
inline void flipEndianness(T &val)
{
  int mid = sizeof(T)/2;
  char *data = reinterpret_cast<char*>(&val);
  for (int i = 0; i < mid; ++i)
    {
    std::swap(data[i], data[sizeof(T) - i - 1]);
    }
}

template <typename SrcT, typename DstT>
void convertBuffer(const SrcT *in, size_t nelms, DstT *out)
{
  for (size_t i = 0; i < nelms; ++i)
    {
    SrcT val = *in++;
    flipEndianness(val);
    *out++ = static_cast<DstT>(val);
    }
}

template <typename T>
void convertBufferWithTypeInfo(const char *in, const TypeInfo &ti,
                               size_t nelms, T *out)
{
  switch (ti.getId())
    {
    case TypeInfo::ID_CHAR:
      convertBuffer(in, nelms, out);
      break;
    case TypeInfo::ID_SHORT:
      convertBuffer(reinterpret_cast<const short*>(in), nelms, out);
      break;
    case TypeInfo::ID_USHORT:
      convertBuffer(reinterpret_cast<const unsigned short*>(in), nelms, out);
      break;
    case TypeInfo::ID_INT:
      convertBuffer(reinterpret_cast<const int*>(in), nelms, out);
      break;
    case TypeInfo::ID_FLOAT:
      convertBuffer(reinterpret_cast<const float*>(in), nelms, out);
      break;
    case TypeInfo::ID_DOUBLE:
      convertBuffer(reinterpret_cast<const double*>(in), nelms, out);
      break;
    default:
      assert(false);
      break;
    }
}

#endif


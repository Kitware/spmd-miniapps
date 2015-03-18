#ifndef smp_h
#define smp_h

#if !defined(USE_TBB_BACKEND) && !defined(USE_OMP_BACKEND)
# error "Please choose a valid smp backend (tbb, omp)"
#endif

#if defined(USE_TBB_BACKEND)

#include <tbb/blocked_range.h>
#include <tbb/blocked_range2d.h>
#include <tbb/blocked_range3d.h>
#include <tbb/enumerable_thread_specific.h>
#include <tbb/parallel_for.h>

namespace smp
{

typedef tbb::blocked_range<unsigned> Range;
typedef tbb::blocked_range2d<unsigned, unsigned> Range2D;
typedef tbb::blocked_range3d<unsigned, unsigned, unsigned> Range3D;

template <typename T>
class ThreadLocalStorage
{
public:
  typedef typename tbb::enumerable_thread_specific<T>::iterator iterator;
  typedef typename tbb::enumerable_thread_specific<T>::const_iterator
    const_iterator;

  ThreadLocalStorage() {}
  explicit ThreadLocalStorage(const T &init) : internal(init) {}

  const T& local() const
  {
    return internal.local();
  }

  T& local()
  {
    return internal.local();
  }

  iterator begin()
  {
    return internal.begin();
  }

  iterator end()
  {
    return internal.end();
  }

  const_iterator begin() const
  {
    return internal.begin();
  }

  const_iterator end() const
  {
    return internal.end();
  }

  size_t size() const
  {
    return internal.size();
  }

private:
  tbb::enumerable_thread_specific<T> internal;
};

template <typename RangeT, typename FuncT>
void parallel_for(const RangeT &range, const FuncT &func)
{
  tbb::parallel_for(range, func);
}

}

#endif

#if defined(USE_OMP_BACKEND)

#include <boost/thread/tss.hpp>

#include <omp.h>

#include <algorithm>
#include <deque>

namespace smp
{

class Range
{
public:
  Range() : from(0), to(0), gsize(1) {}
  Range(unsigned from, unsigned to, unsigned gsize = 1)
    : from(from), to(to), gsize(gsize) {}

  unsigned begin() const
  {
    return from;
  }

  unsigned end() const
  {
    return to;
  }

  unsigned grain() const
  {
    return gsize;
  }

private:
  unsigned from, to, gsize;
};

class Range2D
{
public:
  Range2D() : r(0, 0), c(0, 0) {}
  Range2D(unsigned rfrom, unsigned rto, unsigned rgrain,
          unsigned cfrom, unsigned cto, unsigned cgrain)
    : r(rfrom, rto, rgrain), c(cfrom, cto, cgrain) {}
  Range2D(unsigned rfrom, unsigned rto, unsigned cfrom, unsigned cto)
    : r(rfrom, rto), c(cfrom, cto) {}

  const Range& rows() const
  {
    return r;
  }

  const Range& cols() const
  {
    return c;
  }

private:
  Range r, c;
};

class Range3D
{
public:
  Range3D() : p(0, 0), r(0, 0), c(0, 0) {}
  Range3D(unsigned pfrom, unsigned pto, unsigned pgrain,
          unsigned rfrom, unsigned rto, unsigned rgrain,
          unsigned cfrom, unsigned cto, unsigned cgrain)
    : p(pfrom, pto, pgrain), r(rfrom, rto, rgrain), c(cfrom, cto, cgrain) {}
  Range3D(unsigned pfrom, unsigned pto, unsigned rfrom, unsigned rto,
          unsigned cfrom, unsigned cto)
    : p(pfrom, pto), r(rfrom, rto), c(cfrom, cto) {}

  const Range& pages() const
  {
    return p;
  }

  const Range& rows() const
  {
    return r;
  }

  const Range& cols() const
  {
    return c;
  }

private:
  Range p, r, c;
};

namespace detail
{

template <typename T>
class AccessContainerIteratorWrapper
{
private:
  typedef AccessContainerIteratorWrapper<T> self;
  typedef typename std::deque<T*>::iterator Wrapee;

public:
  explicit AccessContainerIteratorWrapper(const Wrapee &it) : internal(it) {}

  self& operator++()
    {
    ++internal;
    return *this;
    }

  self operator++(int)
    {
    self copy = *this;
    ++internal;
    return copy;
    }

  bool operator==(const self& other) const
    {
    return internal == other.internal;
    }

  bool operator!=(const self& other) const
    {
    return internal != other.internal;
    }

  T& operator*() const
    {
    return **internal;
    }

  T* operator->() const
    {
    return *internal;
    }

private:
  Wrapee internal;
};

template <typename T>
class AccessContainerConstIteratorWrapper
{
private:
  typedef AccessContainerConstIteratorWrapper<T> self;
  typedef typename std::deque<T*>::const_iterator Wrapee;

public:
  explicit AccessContainerConstIteratorWrapper(const Wrapee &it)
    : internal(it) {}

  self& operator++()
    {
    ++internal;
    return *this;
    }

  self operator++(int)
    {
    self copy = *this;
    ++internal;
    return copy;
    }

  bool operator==(const self& other) const
    {
    return internal == other.internal;
    }

  bool operator!=(const self& other) const
    {
    return internal != other.internal;
    }

  const T& operator*() const
    {
    return **internal;
    }

  const T* operator->() const
    {
    return *internal;
    }

private:
  Wrapee internal;
};

}

template <typename T>
class ThreadLocalStorage
{
public:
  typedef detail::AccessContainerIteratorWrapper<T> iterator;
  typedef detail::AccessContainerConstIteratorWrapper<T> const_iterator;

  ThreadLocalStorage()
  {
    omp_init_lock(&accessLock);
  }

  explicit ThreadLocalStorage(const T &init) : exemplar(init)
  {
    omp_init_lock(&accessLock);
  }

  ~ThreadLocalStorage()
  {
    omp_destroy_lock(&accessLock);
  }

  const T& local() const
  {
    return *getLocal();
  }

  T& local()
  {
    return *getLocal();
  }

  iterator begin()
  {
    return detail::AccessContainerIteratorWrapper<T>(access.begin());
  }

  iterator end()
  {
    return detail::AccessContainerIteratorWrapper<T>(access.end());
  }

  const_iterator begin() const
  {
    return detail::AccessContainerConstIteratorWrapper<T>(access.begin());
  }

  const_iterator end() const
  {
    return detail::AccessContainerConstIteratorWrapper<T>(access.end());
  }

  size_t size() const
  {
    return access.size();
  }

private:
  T* getLocal()
  {
    if (!internal.get())
      {
      T *inst = new T(exemplar);
      internal.reset(inst);

      omp_set_lock(&accessLock);
      access.push_back(inst);
      omp_unset_lock(&accessLock);
      }
    return internal.get();
  }

  T exemplar;
  boost::thread_specific_ptr<T> internal;
  std::deque<T*> access;
  omp_lock_t accessLock;
};

template <typename FuncT>
void parallel_for(const Range &range, const FuncT &func)
{
  unsigned nblocks = (range.end() - range.begin() + range.grain() - 1) /
                     range.grain();

# pragma omp parallel for schedule(runtime)
  for (unsigned i = 0; i < nblocks; ++i)
    {
    unsigned from = i * range.grain();
    unsigned to = std::min(from + range.grain(), range.end());
    func(Range(from, to));
    }
}

template <typename FuncT>
void parallel_for(const Range2D &range, const FuncT &func)
{
  unsigned numBlockRows = (range.rows().end() - range.rows().begin() +
                          range.rows().grain() - 1) / range.rows().grain();
  unsigned numBlockCols = (range.cols().end() - range.cols().begin() +
                          range.cols().grain() - 1) / range.cols().grain();
  unsigned nblocks = numBlockRows * numBlockCols;

# pragma omp parallel for schedule(runtime)
  for (unsigned i = 0; i < nblocks; ++i)
    {
    unsigned blockRowIdx = i/numBlockCols;
    unsigned blockColIdx = i%numBlockCols;

    unsigned rfrom = blockRowIdx * range.rows().grain();
    unsigned rto = std::min(rfrom + range.rows().grain(), range.rows().end());
    unsigned cfrom = blockColIdx * range.cols().grain();
    unsigned cto = std::min(cfrom + range.cols().grain(), range.cols().end());

    func(Range2D(rfrom, rto, cfrom, cto));
   }
}

template <typename FuncT>
void parallel_for(const Range3D &range, const FuncT &func)
{
  unsigned numBlockPages = (range.pages().end() - range.pages().begin() +
                           range.pages().grain() - 1) / range.pages().grain();
  unsigned numBlockRows = (range.rows().end() - range.rows().begin() +
                          range.rows().grain() - 1) / range.rows().grain();
  unsigned numBlockCols = (range.cols().end() - range.cols().begin() +
                          range.cols().grain() - 1) / range.cols().grain();
  unsigned nblocks = numBlockPages * numBlockRows * numBlockCols;
  unsigned nblocksPerPage = numBlockRows * numBlockCols;

# pragma omp parallel for schedule(runtime)
  for (unsigned i = 0; i < nblocks; ++i)
    {
    unsigned blockPageIdx = i/nblocksPerPage;
    unsigned blockRowIdx = (i%nblocksPerPage)/numBlockCols;
    unsigned blockColIdx = (i%nblocksPerPage)%numBlockCols;

    unsigned pfrom = blockPageIdx * range.pages().grain();
    unsigned pto = std::min(pfrom+range.pages().grain(), range.pages().end());
    unsigned rfrom = blockRowIdx * range.rows().grain();
    unsigned rto = std::min(rfrom + range.rows().grain(), range.rows().end());
    unsigned cfrom = blockColIdx * range.cols().grain();
    unsigned cto = std::min(cfrom + range.cols().grain(), range.cols().end());

    func(Range3D(pfrom, pto, rfrom, rto, cfrom, cto));
   }
}

}

#endif

#endif // smp_h

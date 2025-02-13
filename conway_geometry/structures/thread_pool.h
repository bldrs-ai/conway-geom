#pragma once

#include <optional>
#include <functional>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>

namespace conway {

#if defined( __EMSCRIPTEN_PTHREADS__ )

class ThreadPool {
public:

  ThreadPool() {

    for ( size_t where = 1, end = std::thread::hardware_concurrency(); where < end; ++where ) {

      threads_.emplace_back( [&]() {

        worker();

      } );
    }
  }

  ~ThreadPool() {

    exit_ = true;

    condition_.notify_all();

    for ( std::thread& thread : threads_ ) {

      thread.join();
    }
  }


  /** Simple parallel for, not re-entrant */
  template < typename FunctionType > 
  void parallel_for( size_t start, size_t end, FunctionType function, size_t threadStride = 2 ) {

    if ( ( end - start ) <= threadStride ) {

      for ( size_t where = start; where < end; ++where ) {

        function( where );
      }

      return;
    }

    std::function< void ( size_t ) > wrapper( function );

    {
      std::lock_guard< std::mutex > guard( lock_ );

      complete_     = start;
      threadStride_ = threadStride;
      end_          = end;
      counter_      = start;

      iteration_ = &wrapper;
    }

    condition_.notify_all();

    work();

    {
      std::unique_lock<std::mutex> latch( lock_ );

      join_.wait( latch, [&]() { return complete_.load() == end; } );
    }

    iteration_ = nullptr;
  }

  static ThreadPool& instance() { 

    if ( !instance_.has_value() ) {

      instance_.emplace();
    }

    return instance_.value();
  }

private:

  void worker() {

    while ( !exit_.load() ) {

      std::unique_lock<std::mutex> latch( lock_ );
      
      condition_.wait( latch, [&]() { return exit_ || counter_.load() < end_.load(); } );

      if ( !exit_ ) {

        work();
      }
    }
  }

  void work() {

    size_t end    = end_.load();
    size_t stride = threadStride_.load();

    while ( counter_.load() < end ) {

      size_t cursor = counter_.fetch_add( stride );

      size_t cursorEnd = std::min( cursor + stride, end );

      while ( cursor < cursorEnd ) {

        (*iteration_)( cursor );

        ++cursor;

        if ( ( ++complete_ ) >= end ) {

          join_.notify_one();
        }
      }
    }
  }
  
  std::atomic< size_t > threadStride_ = 0;
  std::atomic< size_t > end_ = 0;
  std::atomic< size_t > counter_ = 0;
  std::atomic< size_t > complete_ = 0;

  std::atomic< bool > exit_ = false;

  std::mutex lock_;
  std::condition_variable condition_;
  std::condition_variable join_;

  std::vector< std::thread > threads_;

  std::function< void ( size_t ) >* iteration_ = nullptr;

  static std::optional< ThreadPool > instance_;

};

#else

class ThreadPool {
  public:
  
    ThreadPool() {
    }
  
    ~ThreadPool() {
    }
  
  
    /** Simple parallel for, not re-entrant */
    template < typename FunctionType > 
    void parallel_for( size_t start, size_t end, FunctionType function, [[maybe_unused]]size_t threadStride = 1 ) {
  
      for ( size_t where = start; where < end; where += increment ) {

        function( where );
      }
    }
  
    static ThreadPool& instance() { 
  
      if ( !instance_.has_value() ) {
  
        instance_.emplace();
      }
  
      return instance_.value();
    }
  
  private:
    
    static std::optional< ThreadPool > instance_;
  
  };
  

#endif

}
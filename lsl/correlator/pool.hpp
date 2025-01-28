#pragma once

#include <vector>
#include <mutex>
#include <tuple>
#include <stdexcept>
#include <fftw3.h>
#include <complex>
#include <unordered_map>
#include <functional>
#include <memory>

/*
 Thread-safe memory pool for 64-bit aligned allocations.
*/

class Aligned64BufferPool {
private:
    struct Buffer {
        void* ptr;
        size_t size_bytes;
        bool in_use;
    };
    
    std::vector<Buffer> buffers;
    std::mutex mtx;
    std::string context_name;

public:
    explicit Aligned64BufferPool(const std::string& context) 
        : context_name(context) {}

    ~Aligned64BufferPool() {
        std::lock_guard<std::mutex> lock(mtx);
        for (auto& buf : buffers) {
            aligned64_free(buf.ptr);
        }
    }

    // Prevent copying and moving
    Aligned64BufferPool(const Aligned64BufferPool&) = delete;
    Aligned64BufferPool& operator=(const Aligned64BufferPool&) = delete;
    Aligned64BufferPool(Aligned64BufferPool&&) = delete;
    Aligned64BufferPool& operator=(Aligned64BufferPool&&) = delete;

    template<typename T>
    T* acquire(size_t count) {
        size_t bytes = count * sizeof(T);
        std::lock_guard<std::mutex> lock(mtx);
        
        // Look for existing buffer of correct size
        for (auto& buf : buffers) {
            if (buf.size_bytes == bytes && !buf.in_use) {
                buf.in_use = true;
                return static_cast<T*>(buf.ptr);
            }
        }

        // Allocate new buffer
        T* newPtr = static_cast<T*>(aligned64_malloc(bytes));
        if (!newPtr) {
            throw std::bad_alloc();
        }
        buffers.push_back({newPtr, bytes, true});
        return newPtr;
    }

    template<typename T>
    void release(T* ptr) {
        std::lock_guard<std::mutex> lock(mtx);
        
        for (auto& buf : buffers) {
            if (buf.ptr == ptr) {
                buf.in_use = false;
                return;
            }
        }
        throw std::runtime_error("Attempted to release unmanaged aligned buffer");
    }

    std::tuple<size_t, size_t> get_stats() {
        std::lock_guard<std::mutex> lock(mtx);
        size_t total = buffers.size();
        size_t in_use = std::count_if(buffers.begin(), buffers.end(),
                                     [](const Buffer& b) { return b.in_use; });
        return {total, in_use};
    }
};

inline Aligned64BufferPool& get_aligned64_buffer_pool(const std::string& context) {
    static std::unordered_map<std::string, std::unique_ptr<Aligned64BufferPool>> caches;
    static std::mutex cache_mtx;
    
    std::lock_guard<std::mutex> lock(cache_mtx);
    auto it = caches.find(context);
    if (it == caches.end()) {
        caches[context] = std::make_unique<Aligned64BufferPool>(context);
        return *caches[context];
    }
    return *it->second;
}


/*
 Thread-safe memory pool for FFTW single-precision allocations.
*/

class FFTWBufferPool {
private:
    struct Buffer {
        void* ptr;
        size_t size_bytes;
        bool in_use;
    };
    
    std::vector<Buffer> buffers;
    std::mutex mtx;
    std::string context_name;

public:
    explicit FFTWBufferPool(const std::string& context) 
        : context_name(context) {}
    ~FFTWBufferPool() {
        std::lock_guard<std::mutex> lock(mtx);
        for (auto& buf : buffers) {
            fftwf_free(buf.ptr);
        }
    }

    FFTWBufferPool(const FFTWBufferPool&) = delete;
    FFTWBufferPool& operator=(const FFTWBufferPool&) = delete;
    FFTWBufferPool(FFTWBufferPool&&) = delete;
    FFTWBufferPool& operator=(FFTWBufferPool&&) = delete;

    template<typename T>
    T* acquire(size_t count) {
        size_t bytes = count * sizeof(T);
        std::lock_guard<std::mutex> lock(mtx);
        
        for (auto& buf : buffers) {
            if (buf.size_bytes == bytes && !buf.in_use) {
                buf.in_use = true;
                return static_cast<T*>(buf.ptr);
            }
        }

        T* newPtr = static_cast<T*>(fftwf_malloc(bytes));
        if (!newPtr) {
            throw std::bad_alloc();
        }

        buffers.push_back({newPtr, bytes, true});
        return newPtr;
    }

    template<typename T>
    void release(T* ptr) {
        std::lock_guard<std::mutex> lock(mtx);
        
        for (auto& buf : buffers) {
            if (buf.ptr == ptr) {
                buf.in_use = false;
                return;
            }
        }
        throw std::runtime_error("Attempted to release unmanaged buffer");
    }

    std::tuple<size_t, size_t> get_stats() {
        std::lock_guard<std::mutex> lock(mtx);
        size_t total = buffers.size();
        size_t in_use = std::count_if(buffers.begin(), buffers.end(),
                                     [](const Buffer& b) { return b.in_use; });
        return {total, in_use};
    }
};

inline FFTWBufferPool& get_fftw_buffer_pool(const std::string& context) {
    static std::unordered_map<std::string, std::unique_ptr<FFTWBufferPool>> caches;
    static std::mutex cache_mtx;

    std::lock_guard<std::mutex> lock(cache_mtx);
    auto it = caches.find(context);
    if (it == caches.end()) {
        caches[context] = std::make_unique<FFTWBufferPool>(context);
        return *caches[context];
    }
    return *it->second;
}

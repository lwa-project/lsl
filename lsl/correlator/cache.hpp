#pragma once

#include <vector>
#include <mutex>
#include <memory>
#include <string>
#include <functional>
#include <unordered_map>
#include <fftw3.h>

/*
 FFTW plan cache for specific function contexts.
*/

class FFTWPlanCache {
private:
    enum class PlanType {
        MANY_R2C,       // Real to complex many transform
        MANY_C2C,       // Complex to complex many transform
        C2C_2D,         // 2D complex to complex
        R2R_1D          // 1D real to real (for REDFT01)
    };

    struct PlanEntry {
        PlanType type;
        std::vector<int> dims;    // Transform dimensions
        int ntransforms;          // Number of transforms (for many)
        int istride;              // Input stride
        int ostride;             // Output stride
        int idist;               // Input distance
        int odist;               // Output distance
        void* in;                // Input buffer
        void* out;               // Output buffer
        int direction;           // FFTW_FORWARD or FFTW_BACKWARD
        int r2r_kind;           // FFTW_REDFT01 etc for R2R
        fftwf_plan plan;        // The FFTW plan
        bool in_use;
    };
    
    std::vector<PlanEntry> plans;
    std::mutex mtx;
    std::string context_name;

public:
    explicit FFTWPlanCache(const std::string& context) 
        : context_name(context) {}

    ~FFTWPlanCache() {
        std::lock_guard<std::mutex> lock(mtx);
        for (auto& entry : plans) {
            fftwf_destroy_plan(entry.plan);
        }
    }

    // No copy or move
    FFTWPlanCache(const FFTWPlanCache&) = delete;
    FFTWPlanCache& operator=(const FFTWPlanCache&) = delete;
    FFTWPlanCache(FFTWPlanCache&&) = delete;
    FFTWPlanCache& operator=(FFTWPlanCache&&) = delete;

    // Real-to-complex many transform - matches fftwf_plan_many_dft_r2c
    fftwf_plan plan_many_dft_r2c(int rank, const int* n, int howmany,
                                float* in, const int* inembed, int istride, int idist,
                                fftwf_complex* out, const int* onembed, int ostride, int odist,
                                unsigned flags) {
        std::lock_guard<std::mutex> lock(mtx);
        
        for (auto& entry : plans) {
            if (entry.type == PlanType::MANY_R2C &&
                entry.dims[0] == n[0] &&
                entry.ntransforms == howmany &&
                entry.istride == istride && entry.ostride == ostride &&
                entry.idist == idist && entry.odist == odist &&
                entry.in == in && entry.out == out && !entry.in_use) {
                entry.in_use = true;
                return entry.plan;
            }
        }

        fftwf_plan plan = fftwf_plan_many_dft_r2c(rank, n, howmany,
                                                in, inembed, istride, idist,
                                                out, onembed, ostride, odist,
                                                flags);
        if (!plan) {
            throw std::runtime_error("Failed to create FFTW plan");
        }

        plans.push_back({
            PlanType::MANY_R2C,
            {n[0]}, howmany, istride, ostride, idist, odist,
            in, out, FFTW_FORWARD, 0, plan, true
        });
        return plan;
    }

    // Complex-to-complex many transform - matches fftwf_plan_many_dft
    fftwf_plan plan_many_dft(int rank, const int* n, int howmany,
                           fftwf_complex* in, const int* inembed, int istride, int idist,
                           fftwf_complex* out, const int* onembed, int ostride, int odist,
                           int sign, unsigned flags) {
        std::lock_guard<std::mutex> lock(mtx);
        
        for (auto& entry : plans) {
            if (entry.type == PlanType::MANY_C2C &&
                entry.dims[0] == n[0] &&
                entry.ntransforms == howmany &&
                entry.istride == istride && entry.ostride == ostride &&
                entry.idist == idist && entry.odist == odist &&
                entry.in == in && entry.out == out &&
                entry.direction == sign && !entry.in_use) {
                entry.in_use = true;
                return entry.plan;
            }
        }

        fftwf_plan plan = fftwf_plan_many_dft(rank, n, howmany,
                                            in, inembed, istride, idist,
                                            out, onembed, ostride, odist,
                                            sign, flags);
        if (!plan) {
            throw std::runtime_error("Failed to create FFTW plan");
        }

        plans.push_back({
            PlanType::MANY_C2C,
            {n[0]}, howmany, istride, ostride, idist, odist,
            in, out, sign, 0, plan, true
        });
        return plan;
    }

    // 2D complex-to-complex transform - matches fftwf_plan_dft_2d
    fftwf_plan plan_dft_2d(int nx, int ny,
                         fftwf_complex* in, fftwf_complex* out,
                         int sign, unsigned flags) {
        std::lock_guard<std::mutex> lock(mtx);
        
        for (auto& entry : plans) {
            if (entry.type == PlanType::C2C_2D &&
                entry.dims[0] == nx && entry.dims[1] == ny &&
                entry.in == in && entry.out == out &&
                entry.direction == sign && !entry.in_use) {
                entry.in_use = true;
                return entry.plan;
            }
        }

        fftwf_plan plan = fftwf_plan_dft_2d(nx, ny, in, out, sign, flags);
        if (!plan) {
            throw std::runtime_error("Failed to create FFTW plan");
        }

        plans.push_back({
            PlanType::C2C_2D,
            {nx, ny}, 0, 0, 0, 0, 0,
            in, out, sign, 0, plan, true
        });
        return plan;
    }

    // 1D real-to-real transform - matches fftwf_plan_r2r_1d
    fftwf_plan plan_r2r_1d(int n, float* in, float* out,
                         fftwf_r2r_kind kind, unsigned flags) {
        std::lock_guard<std::mutex> lock(mtx);
        
        for (auto& entry : plans) {
            if (entry.type == PlanType::R2R_1D &&
                entry.dims[0] == n &&
                entry.in == in && entry.out == out &&
                entry.r2r_kind == kind && !entry.in_use) {
                entry.in_use = true;
                return entry.plan;
            }
        }

        fftwf_plan plan = fftwf_plan_r2r_1d(n, in, out, kind, flags);
        if (!plan) {
            throw std::runtime_error("Failed to create FFTW plan");
        }

        plans.push_back({
            PlanType::R2R_1D,
            {n}, 0, 0, 0, 0, 0,
            in, out, 0, kind, plan, true
        });
        return plan;
    }

    void release_plan(fftwf_plan plan) {
        std::lock_guard<std::mutex> lock(mtx);
        
        for (auto& entry : plans) {
            if (entry.plan == plan) {
                entry.in_use = false;
                return;
            }
        }
        throw std::runtime_error("Attempted to release unmanaged plan");
    }

    std::tuple<size_t, size_t> get_stats() {
        std::lock_guard<std::mutex> lock(mtx);
        size_t total = plans.size();
        size_t in_use = std::count_if(plans.begin(), plans.end(),
                                     [](const PlanEntry& p) { return p.in_use; });
        return {total, in_use};
    }
};

inline FFTWPlanCache& get_fftw_plan_cache(const std::string& context) {
    static std::unordered_map<std::string, std::unique_ptr<FFTWPlanCache>> caches;
    static std::mutex cache_mtx;

    std::lock_guard<std::mutex> lock(cache_mtx);
    auto it = caches.find(context);
    if (it == caches.end()) {
        caches[context] = std::make_unique<FFTWPlanCache>(context);
        return *caches[context];
    }
    return *it->second;
}

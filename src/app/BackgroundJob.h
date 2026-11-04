#pragma once

// A long-running computation executed on a worker thread while the UI keeps
// rendering. The UI thread starts it, then polls once per frame; the result
// (or the error message) is delivered exactly once, on the UI thread.
//
// The work function must not touch viewer state: give it copies of
// everything it needs (expression text, settings, output paths).

#include <chrono>
#include <exception>
#include <future>
#include <optional>
#include <string>
#include <utility>

namespace iso_app {

template <typename Result>
class BackgroundJob
{
public:
    BackgroundJob() = default;
    BackgroundJob(const BackgroundJob&) = delete;
    BackgroundJob& operator=(const BackgroundJob&) = delete;
    ~BackgroundJob() { wait(); }

    bool running() const { return future_.valid(); }

    // Ignored while a previous job is still running.
    template <typename Work>
    void start(Work&& work)
    {
        if (!running())
            future_ = std::async(std::launch::async, std::forward<Work>(work));
    }

    // If the job has finished, calls on_success(Result) or on_failure(message).
    template <typename OnSuccess, typename OnFailure>
    void poll(OnSuccess&& on_success, OnFailure&& on_failure)
    {
        if (!running() || future_.wait_for(std::chrono::seconds(0)) != std::future_status::ready)
            return;
        std::optional<Result> result;
        std::string error;
        try
        {
            result.emplace(future_.get());
        }
        catch (const std::exception& e)
        {
            error = e.what();
        }
        catch (...)
        {
            error = "unknown error";
        }
        if (result)
            on_success(std::move(*result));
        else
            on_failure(error);
    }

    void wait()
    {
        if (future_.valid())
            future_.wait();
    }

private:
    std::future<Result> future_;
};

} // namespace iso_app

#include "NanoLog.h"
#include "RuntimeLogger.h"

/**
 * This file implements the public API to NanoLog
 */
namespace NanoLog {
    using namespace NanoLogInternal;

    void preallocate() {
        RuntimeLogger::preallocate();
    }

    void setLogFile(const char *filename) {
        RuntimeLogger::setLogFile(filename);
    }

    LogLevel getLogLevel() {
        return RuntimeLogger::getLogLevel();
    }

    void setLogLevel(LogLevel logLevel) {
        RuntimeLogger::setLogLevel(logLevel);
    }

    void poll() {
      RuntimeLogger::poll();
    }

    void setThreadName(const char* name) {
      RuntimeLogger::setThreadName(name);
    }

    void setLogCB(LogCBFn cb, LogLevel maxCBLogLevel, uint64_t minCBPeriodInSec) {
      RuntimeLogger::setLogCB(cb, maxCBLogLevel, minCBPeriodInSec);
    }
};

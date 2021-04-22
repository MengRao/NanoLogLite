/* Copyright (c) 2016-2020 Stanford University
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR(S) DISCLAIM ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL AUTHORS BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */

#ifndef RUNTIME_NANOLOG_H
#define RUNTIME_NANOLOG_H

#include <cassert>

#include <mutex>
#include <string>
#include <thread>
#include <vector>
#include <unistd.h>
#include <sys/syscall.h>
#include "NanoLog.h"

namespace NanoLogInternal {
using namespace NanoLog;

struct TwoNibbles
{
  uint8_t first : 4;
  uint8_t second : 4;
};

/**
 * Given an unsigned integer and a char array, find the fewest number of
 * bytes needed to represent the integer, copy that many bytes into the
 * char array, and bump the char array pointer.
 *
 * \param[in/out] buffer
 *      char array pointer used to store the compressed value and bump
 * \param val
 *      Unsigned integer to pack into the buffer
 *
 * \return
 *      Special 4-bit value indicating how the primitive was packed
 */
template<typename T>
inline typename std::enable_if<std::is_integral<T>::value && !std::is_signed<T>::value, int>::type pack(char** buffer,
                                                                                                        T val) {
  // Binary search for the smallest container. It is also worth noting that
  // with -O3, the compiler will strip out extraneous if-statements based on T
  // For example, if T is uint16_t, it would only leave the t < 1U<<8 check

  // TODO(syang0) Is this too costly vs. a simple for loop?
  int numBytes;
  if (val < (1UL << 8)) {
    numBytes = 1;
  }
  else if (val < (1UL << 16)) {
    numBytes = 2;
  }
  else if (val < (1UL << 24)) {
    numBytes = 3;
  }
  else if (val < (1UL << 32)) {
    numBytes = 4;
  }
  else if (val < (1UL << 40)) {
    numBytes = 5;
  }
  else if (val < (1UL << 48)) {
    numBytes = 6;
  }
  else if (val < (1UL << 56)) {
    numBytes = 7;
  }
  else {
    numBytes = 8;
  }

  // Although we store the entire value here, we take advantage of the fact
  // that x86-64 is little-endian (storing the least significant bits first)
  // and lop off the rest by only partially incrementing the buffer pointer
  std::memcpy(*buffer, &val, sizeof(T));
  *buffer += numBytes;

  return numBytes;
}

/**
 * Below are a series of pack functions that take in a signed integer,
 * test to see if the value will be smaller if negated, and then invoke
 * the unsigned version of the pack() function above.
 *
 * \param[in/out] buffer
 *      char array to copy the value into and bump
 * \param val
 *      Unsigned integer to pack into the buffer
 *
 * \return
 *      Special 4-bit value indicating how the primitive was packed
 */
inline int pack(char** buffer, int32_t val) {
  if (val >= 0 || val <= int32_t(-(1 << 24)))
    return pack<uint32_t>(buffer, static_cast<uint32_t>(val));
  else
    return 8 + pack<uint32_t>(buffer, static_cast<uint32_t>(-val));
}

inline int pack(char** buffer, int64_t val) {
  if (val >= 0 || val <= int64_t(-(1LL << 56)))
    return pack<uint64_t>(buffer, static_cast<uint64_t>(val));
  else
    return 8 + pack<uint64_t>(buffer, static_cast<uint64_t>(-val));
}

// TODO(syang0) we should measure the performance of doing it this way
// vs taking both the negated and non-negated versions and encoding the smaller
inline int pack(char** buffer, long long int val) {
  if (val >= 0 || val <= int64_t(-(1LL << 56)))
    return pack<uint64_t>(buffer, static_cast<uint64_t>(val));
  else
    return 8 + pack<uint64_t>(buffer, static_cast<uint64_t>(-val));
}

/**
 * Pointer specialization for the pack template that will copy the value
 * without compression.
 *
 * \param[in/out] buffer
 *      char array to copy the integer into and bump
 * \param val
 *      Unsigned integer to pack into the buffer
 *
 * \return - Special 4-bit value indicating how the primitive was packed
 */
template<typename T>
inline int pack(char** in, T* pointer) {
  return pack<uint64_t>(in, reinterpret_cast<uint64_t>(pointer));
}

/**
 * Floating point specialization for the pack template that will copy the value
 * without compression.
 *
 * \param[in/out] buffer
 *      char array to copy the float into and bump
 * \param val
 *      Unsigned integer to pack into the buffer
 *
 * \return - Special 4-bit value indicating how the primitive was packed
 */
template<typename T>
inline typename std::enable_if<std::is_floating_point<T>::value, int>::type pack(char** buffer, T val) {
  std::memcpy(*buffer, &val, sizeof(T));
  *buffer += sizeof(T);
  return sizeof(T);
}

/**
 * Below are various unpack functions that will take in a data array pointer
 * and the special pack code, return the value originally pack()-ed and bump
 * the pointer to "consume" the pack()-ed value.
 *
 * \param in
 *      data array pointer to read the data back from and increment.
 * \param packResult
 *      special 4-bit code returned from pack()
 *
 * \return
 *      original full-width value before compression
 */

template<typename T>
inline typename std::enable_if<!std::is_floating_point<T>::value && !std::is_pointer<T>::value, T>::type
unpack(const char** in, uint8_t packResult) {
  int64_t packed = 0;

  if (packResult <= 8) {
    memcpy(&packed, (*in), packResult);
    (*in) += packResult;
    return static_cast<T>(packed);
  }

  int bytes = packResult == 0 ? 16 : (0x0f & (packResult - 8));
  memcpy(&packed, (*in), bytes);
  (*in) += bytes;

  return static_cast<T>(-packed);
}

template<typename T>
inline typename std::enable_if<std::is_pointer<T>::value, T>::type unpack(const char** in, uint8_t packNibble) {
  return (T*)(unpack<uint64_t>(in, packNibble));
}

template<typename T>
inline typename std::enable_if<std::is_floating_point<T>::value, T>::type unpack(const char** in, uint8_t packNibble) {
  if (packNibble == sizeof(float)) {
    float result;
    std::memcpy(&result, *in, sizeof(float));
    *in += sizeof(float);
    return result;
  }

  // Double case
  T result;
  int bytes = packNibble == 0 ? 16 : packNibble;
  std::memcpy(&result, (*in), bytes);
  (*in) += bytes;

  return result;
}

/**
 * Given a stream of nibbles, return the total number of bytes used to represent
 * the values encoded with the nibbles.
 *
 * \param nibbles
 *      The start of the nibbles
 * \param numNibbles
 *      Number of nibbles to process
 *
 * \return
 *      The number of bytes encoded used to encode the values
 */
inline static uint32_t getSizeOfPackedValues(const TwoNibbles* nibbles, uint32_t numNibbles) {
  uint32_t size = 0;
  for (uint32_t i = 0; i < numNibbles / 2; ++i) {
    size += nibbles[i].first + nibbles[i].second;
    if (nibbles[i].first == 0) size += 16;
    if (nibbles[i].first > 0x8) size -= 8;
    if (nibbles[i].second == 0) size += 16;
    if (nibbles[i].second > 0x8) size -= 8;
  }

  if (numNibbles & 0x1) {
    size += nibbles[numNibbles / 2].first;
    if (nibbles[numNibbles / 2].first == 0) size += 16;
    if (nibbles[numNibbles / 2].first > 0x8) size -= 8;
  }

  return size;
}

/**
 * Describes the type of parameter that would be passed into a printf-like
 * function.
 *
 * These types are optimized to store enough information to determine
 * (a) whether a 'const char*' parameter indicates string (%s) or not (%p)
 * (b) if a string parameter (%s) needs to be truncated due to precision
 * (c) whether a parameter is a dynamic precision/width specifier
 */
enum ParamType : int32_t
{
  // Indicates that there is a problem with the parameter
  INVALID = -6,

  // Indicates a dynamic width (i.e. the '*' in  %*.d)
  DYNAMIC_WIDTH = -5,

  // Indicates dynamic precision (i.e. the '*' in %.*d)
  DYNAMIC_PRECISION = -4,

  // Indicates that the parameter is not a string type (i.e. %d, %lf)
  NON_STRING = -3,

  // Indicates the parameter is a string and has a dynamic precision
  // (i.e. '%.*s' )
  STRING_WITH_DYNAMIC_PRECISION = -2,

  // Indicates a string with no precision specified (i.e. '%s' )
  STRING_WITH_NO_PRECISION = -1,

  // All non-negative values indicate a string with a precision equal to its
  // enum value casted as an int32_t
  STRING = 0
};

// Default, uninitialized value for log identifiers associated with log
// invocation sites.
static constexpr int UNASSIGNED_LOGID = -1;

// https://github.com/MengRao/SPSC_Queue
template<uint32_t Bytes>
class SPSCVarQueueOPT
{
public:
  struct MsgHeader
  {
    // size of this msg, including header itself
    // auto set by lib, can be read by user
    uint16_t size;
    uint16_t msg_type;
    // userdata can be used by caller, e.g. save timestamp or other stuff
    // we assume that user_msg is 8 types alligned so there'll be 4 bytes padding anyway, otherwise we can choose to
    // eliminate userdata
    uint32_t userdata;
  };
  static constexpr uint32_t BLK_CNT = Bytes / sizeof(MsgHeader);

  MsgHeader* alloc(uint16_t size_) {
    size = size_ + sizeof(MsgHeader);
    uint32_t blk_sz = (size + sizeof(MsgHeader) - 1) / sizeof(MsgHeader);
    if (blk_sz >= free_write_cnt) {
      asm volatile("" : "=m"(read_idx) : :); // force read memory
      uint32_t read_idx_cache = read_idx;
      if (read_idx_cache <= write_idx) {
        free_write_cnt = BLK_CNT - write_idx;
        if (blk_sz >= free_write_cnt && read_idx_cache != 0) { // wrap around
          blk[0].size = 0;
          asm volatile("" : : "m"(blk) :); // memory fence
          blk[write_idx].size = 1;
          write_idx = 0;
          free_write_cnt = read_idx_cache;
        }
      }
      else {
        free_write_cnt = read_idx_cache - write_idx;
      }
      if (free_write_cnt <= blk_sz) {
        return nullptr;
      }
    }
    return &blk[write_idx];
  }

  void push() {
    uint32_t blk_sz = (size + sizeof(MsgHeader) - 1) / sizeof(MsgHeader);
    blk[write_idx + blk_sz].size = 0;

    asm volatile("" : : "m"(blk) :); // memory fence
    blk[write_idx].size = size;
    write_idx += blk_sz;
    free_write_cnt -= blk_sz;
  }

  template<typename Writer>
  bool tryPush(uint16_t size, Writer writer) {
    MsgHeader* header = alloc(size);
    if (!header) return false;
    writer(header);
    push();
    return true;
  }

  MsgHeader* front() {
    asm volatile("" : "=m"(blk) : :); // force read memory
    uint16_t size = blk[read_idx].size;
    if (size == 1) { // wrap around
      read_idx = 0;
      size = blk[0].size;
    }
    if (size == 0) return nullptr;
    return &blk[read_idx];
  }

  void pop() {
    asm volatile("" : "=m"(blk) : "m"(read_idx) :); // memory fence
    uint32_t blk_sz = (blk[read_idx].size + sizeof(MsgHeader) - 1) / sizeof(MsgHeader);
    read_idx += blk_sz;
    asm volatile("" : : "m"(read_idx) :); // force write memory
  }

  template<typename Reader>
  bool tryPop(Reader reader) {
    MsgHeader* header = front();
    if (!header) return false;
    reader(header);
    pop();
    return true;
  }

private:
  alignas(64) MsgHeader blk[BLK_CNT] = {};

  alignas(128) uint32_t write_idx = 0;
  uint32_t free_write_cnt = BLK_CNT;
  uint16_t size;

  alignas(128) uint32_t read_idx = 0;
};

// https://github.com/MengRao/tscns
class TSCNS
{
public:
  // If you haven't calibrated tsc_ghz on this machine, set tsc_ghz as 0.0 and it will auto wait 10 ms and calibrate.
  // Of course you can calibrate again later(e.g. after system init is done) and the longer you wait the more precise
  // tsc_ghz calibrate can get. It's a good idea that user waits as long as possible(more than 1 min) once, and save the
  // resultant tsc_ghz returned from calibrate() somewhere(e.g. config file) on this machine for future use. Or you can
  // cheat, see README and cheat.cc for details.
  //
  // If you have calibrated/cheated before on this machine as above, set tsc_ghz and skip calibration.
  //
  // One more thing: you can re-init and calibrate TSCNS at later times if you want to re-sync with
  // system time in case of NTP or manual time changes.
  double init(double tsc_ghz = 0.0) {
    syncTime(base_tsc, base_ns);
    if (tsc_ghz > 0) {
      tsc_ghz_inv = 1.0 / tsc_ghz;
      adjustOffset();
      return tsc_ghz;
    }
    else {
      return calibrate();
    }
  }

  double calibrate(uint64_t min_wait_ns = 10000000) {
    uint64_t delayed_tsc, delayed_ns;
    do {
      syncTime(delayed_tsc, delayed_ns);
    } while ((delayed_ns - base_ns) < min_wait_ns);
    tsc_ghz_inv = (double)(int64_t)(delayed_ns - base_ns) / (int64_t)(delayed_tsc - base_tsc);
    adjustOffset();
    return 1.0 / tsc_ghz_inv;
  }

  static uint64_t rdtsc() { return __builtin_ia32_rdtsc(); }

  uint64_t tsc2ns(uint64_t tsc) const { return ns_offset + (int64_t)((int64_t)tsc * tsc_ghz_inv); }

  uint64_t rdns() const { return tsc2ns(rdtsc()); }

  // If you want cross-platform, use std::chrono as below which incurs one more function call:
  // return std::chrono::high_resolution_clock::now().time_since_epoch().count();
  static uint64_t rdsysns() {
    timespec ts;
    ::clock_gettime(CLOCK_REALTIME, &ts);
    return ts.tv_sec * 1000000000 + ts.tv_nsec;
  }

  // For checking purposes, see test.cc
  uint64_t rdoffset() const { return ns_offset; }

private:
  // Linux kernel sync time by finding the first try with tsc diff < 50000
  // We do better: we find the try with the mininum tsc diff
  void syncTime(uint64_t& tsc, uint64_t& ns) {
    const int N = 10;
    uint64_t tscs[N + 1];
    uint64_t nses[N + 1];

    tscs[0] = rdtsc();
    for (int i = 1; i <= N; i++) {
      nses[i] = rdsysns();
      tscs[i] = rdtsc();
    }

    int best = 1;
    for (int i = 2; i <= N; i++) {
      if (tscs[i] - tscs[i - 1] < tscs[best] - tscs[best - 1]) best = i;
    }
    tsc = (tscs[best] + tscs[best - 1]) >> 1;
    ns = nses[best];
  }

  void adjustOffset() { ns_offset = base_ns - (int64_t)((int64_t)base_tsc * tsc_ghz_inv); }

  alignas(64) double tsc_ghz_inv = 1.0; // make sure tsc_ghz_inv and ns_offset are on the same cache line
  uint64_t ns_offset = 0;
  uint64_t base_tsc = 0;
  uint64_t base_ns = 0;
};

extern TSCNS tscns;
extern uint64_t midnight_ns;

/**
 * Stores the static log information associated with a log invocation site
 * (i.e. filename/line/fmtString combination).
 */
struct StaticLogInfo
{

  // Function signature of the compression function used in the
  // non-preprocessor version of NanoLog
  typedef void (*CompressionFn)(int, const ParamType*, char**, char**);

  // Constructor
  constexpr StaticLogInfo(CompressionFn compress, const char* filename, const uint32_t lineNum, LogLevel severity,
                          const char* fmtString, const int numNibbles, const ParamType* paramTypes)
    : compressionFunction(compress)
    , filename(filename)
    , lineNum(lineNum)
    , severity(severity)
    , numPrintFragments(0)
    , formatString(fmtString)
    , numNibbles(numNibbles)
    , paramTypes(paramTypes)
    , fragments(nullptr)
    , last_ts(0) {}

  bool createFragments(char** microCode);

  // Stores the compression function to be used on the log's dynamic arguments
  CompressionFn compressionFunction;

  // File where the log invocation is invoked
  const char* filename;

  // Line number in the file for the invocation
  const uint32_t lineNum;

  // LogLevel severity associated with the log invocation
  const LogLevel severity;

  uint8_t numPrintFragments;

  // printf format string associated with the log invocation
  const char* formatString;

  // Number of nibbles needed to compress all the non-string log arguments
  const int numNibbles;

  // Mapping of parameter index (i.e. order in which it appears in the
  // argument list starting at 0) to parameter type as inferred from the
  // printf log message invocation
  const ParamType* paramTypes;

  void* fragments;

  uint64_t last_ts;
};

/**
 * RuntimeLogger provides runtime support to the C++ code generated by the
 * Preprocessor component.
 * Its main responsibilities are to manage fast thread-local storage to stage
 * uncompressed log messages and manage a background thread to compress the
 * log messages to an output file.
 */
class RuntimeLogger
{
public:
  using MsgHeader = SPSCVarQueueOPT<1 << 20>::MsgHeader;
  /**
   * See function below.
   */
  inline void registerInvocationSite_internal(int& logId, StaticLogInfo info) {
    // TODO(syang0) Make this into a spin lock
    std::lock_guard<std::mutex> lock(nanoLogSingleton.registrationMutex);

    if (logId != UNASSIGNED_LOGID) return;

    logId = static_cast<int32_t>(invocationSites.size());
    invocationSites.push_back(info);
  }

  /**
   * Assigns a globally unique identifier to static log information and
   * stages it for persistence to disk.
   *
   * \param info
   *      Static log info to associate and persist
   *
   * \param[in/out] logId
   *       Unique log identifier to be assigned. A value other than -1
   *       indicates that the id has already been assigned and this
   *       function becomes a no-op.
   */
  static inline void registerInvocationSite(StaticLogInfo info, int& logId) {
    nanoLogSingleton.registerInvocationSite_internal(logId, info);
  }

  /**
   * Allocate thread-local space for the generated C++ code to store an
   * uncompressed log message, but do not make it available for compression
   * yet. The caller should invoke finishAlloc() to make the space visible
   * to the compression thread and this function shall not be invoked
   * again until the corresponding finishAlloc() is invoked first.
   *
   * Note this will block of the buffer is full.
   *
   * \param nbytes
   *      number of bytes to allocate in the
   *
   * \return
   *      pointer to the allocated space
   */
  static inline MsgHeader* reserveAlloc(size_t nbytes) {
    if (stagingBuffer == nullptr) nanoLogSingleton.ensureStagingBufferAllocated();

    // NOLINTNEXTLINE(clang-analyzer-core.CallAndMessage)
    return stagingBuffer->reserveProducerSpace(nbytes);
  }

  /**
   * Complement to reserveAlloc, makes the bytes previously
   * reserveAlloc()-ed visible to the compression/output thread.
   *
   * \param nbytes
   *      Number of bytes to make visible
   */
  static inline void finishAlloc() { stagingBuffer->finishReservation(); }

  static void preallocate();
  static void setLogFile(const char* filename);
  static void setLogLevel(LogLevel logLevel);
  static void poll() { nanoLogSingleton.poll_(); }

  static inline LogLevel getLogLevel() { return nanoLogSingleton.currentLogLevel; }

  static void setThreadName(const char* name) {
    nanoLogSingleton.ensureStagingBufferAllocated();
    stagingBuffer->setName(name);
  }

  static void setLogCB(LogCBFn cb, LogLevel maxCBLogLevel, uint64_t minCBPeriodInSec) {
    nanoLogSingleton.setLogCB_(cb, maxCBLogLevel, minCBPeriodInSec);
  }

private:
  // Forward Declarations
  class StagingBuffer;
  class StagingBufferDestroyer;

  // Storage for staging uncompressed log statements for compression
  static __thread StagingBuffer* stagingBuffer;

  // Destroys the __thread StagingBuffer upon its own destruction, which
  // is synchronized with thread death
  static thread_local StagingBufferDestroyer sbc;

  // Singleton RuntimeLogger that manages the thread-local structures and
  // background output thread.
  static RuntimeLogger nanoLogSingleton;

  RuntimeLogger();

  ~RuntimeLogger();

  void compressionThreadMain();

  void setLogFile_internal(const char* filename);

  void poll_();

  void setLogCB_(LogCBFn cb, LogLevel maxCBLogLevel_, uint64_t minCBPeriodInSec) {
    logCB = cb;
    maxCBLogLevel = maxCBLogLevel_;
    minCBPeriod = minCBPeriodInSec * 1000000000;
  }

  /**
   * Allocates thread-local structures if they weren't already allocated.
   * This is used by the generated C++ code to ensure it has space to
   * log uncompressed messages to and by the user if they wish to
   * preallocate the data structures on thread creation.
   */
  inline void ensureStagingBufferAllocated() {
    if (stagingBuffer == nullptr) {
      std::unique_lock<std::mutex> guard(bufferMutex);
      // Unlocked for the expensive StagingBuffer allocation
      guard.unlock();
      stagingBuffer = new StagingBuffer();
      guard.lock();

      threadBuffers.push_back(stagingBuffer);
    }
  }

  // Globally the thread-local stagingBuffers
  std::vector<StagingBuffer*> threadBuffers;

  // Protects reads and writes to threadBuffers
  std::mutex bufferMutex;

  FILE* outputFp;

  // Minimum log level that RuntimeLogger will accept. Anything lower will
  // be dropped.
  LogLevel currentLogLevel;

  // Used to control access to invocationSites
  std::mutex registrationMutex;

  // Maps unique identifiers to log invocation sites encountered thus far
  // by the non-preprocessor version of NanoLog
  std::vector<StaticLogInfo> invocationSites;
  std::vector<StaticLogInfo> bgInvocationSites;
  char* rawFragments;

  char* rawFragmentEnd;

  LogCBFn logCB;

  LogLevel maxCBLogLevel;

  uint64_t minCBPeriod;

  /**
   * Implements a circular FIFO producer/consumer byte queue that is used
   * to hold the dynamic information of a NanoLog log statement (producer)
   * as it waits for compression via the NanoLog background thread
   * (consumer). There exists a StagingBuffer for every thread that uses
   * the NanoLog system.
   */
  class StagingBuffer
  {
  public:
    /**
     * Attempt to reserve contiguous space for the producer without
     * making it visible to the consumer. The caller should invoke
     * finishReservation() before invoking reserveProducerSpace()
     * again to make the bytes reserved visible to the consumer.
     *
     * This mechanism is in place to allow the producer to initialize
     * the contents of the reservation before exposing it to the
     * consumer. This function will block behind the consumer if
     * there's not enough space.
     *
     * \param nbytes
     *      Number of bytes to allocate
     *
     * \return
     *      Pointer to at least nbytes of contiguous space
     */
    inline MsgHeader* reserveProducerSpace(size_t nbytes) { return varq.alloc(nbytes); }

    inline void finishReservation() { varq.push(); }

    MsgHeader* peek() { return varq.front(); }

    inline void consume() { varq.pop(); }

    void setName(const char* name_) { strncpy(name, name_, sizeof(name) - 1); }

    const char* getName() { return name; }

    /**
     * Returns true if it's safe for the compression thread to delete
     * the StagingBuffer and remove it from the global vector.
     *
     * \return
     *      true if its safe to delete the StagingBuffer
     */
    bool checkCanDelete() { return shouldDeallocate; }

    StagingBuffer()
      : varq()
      , shouldDeallocate(false) {
      // Empty function, but causes the C++ runtime to instantiate the
      // sbc thread_local (see documentation in function).
      sbc.stagingBufferCreated();

      uint32_t tid = static_cast<pid_t>(::syscall(SYS_gettid));
      snprintf(name, sizeof(name), "%d", tid);
    }

    ~StagingBuffer() {}

  private:
    SPSCVarQueueOPT<1 << 20> varq;

    bool shouldDeallocate;

    char name[16] = {0};

    friend RuntimeLogger;
    friend StagingBufferDestroyer;
  };

  // This class is intended to be instantiated as a C++ thread_local to
  // synchronize marking the thread local stagingBuffer for deletion with
  // thread death.
  //
  // The reason why this class exists rather than wrapping the stagingBuffer
  // in a unique_ptr or declaring the stagingBuffer itself to be thread_local
  // is because of performance. Dereferencing the former costs 10 ns and the
  // latter allocates large amounts of resources for every thread that is
  // created, which is wasteful for threads that do not use the RuntimeLogger.
  class StagingBufferDestroyer
  {
  public:
    // TODO(syang0) I wonder if it'll be better if stagingBuffer was
    // actually a thread_local wrapper with dereference operators
    // implemented.

    explicit StagingBufferDestroyer() {}

    // Weird C++ hack; C++ thread_local are instantiated upon first use
    // thus the StagingBuffer has to invoke this function in order
    // to instantiate this object.
    void stagingBufferCreated() {}

    virtual ~StagingBufferDestroyer() {
      if (stagingBuffer != nullptr) {
        stagingBuffer->shouldDeallocate = true;
        stagingBuffer = nullptr;
      }
    }
  };

}; // RuntimeLogger
}; // Namespace NanoLogInternal

#endif /* RUNTIME_NANOLOG_H */

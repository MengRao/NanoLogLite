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

#include <fcntl.h>
#include <iosfwd>
#include <iostream>
#include <locale>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <unistd.h>
#include <regex>

#include "RuntimeLogger.h"

namespace NanoLogInternal {

// Define the static members of RuntimeLogger here
__thread RuntimeLogger::StagingBuffer* RuntimeLogger::stagingBuffer = nullptr;
thread_local RuntimeLogger::StagingBufferDestroyer RuntimeLogger::sbc;
RuntimeLogger RuntimeLogger::nanoLogSingleton;

TSCNS tscns;
uint64_t midnight_ns;

class TimeIniter
{
public:
  TimeIniter() {
    tscns.init();
    time_t rawtime = tscns.rdns() / 1000000000;
    struct tm* timeinfo = localtime(&rawtime);
    timeinfo->tm_sec = timeinfo->tm_min = timeinfo->tm_hour = 0;
    midnight_ns = mktime(timeinfo) * 1000000000;
  }
};

static TimeIniter _;

/**
 * Describes how to interpret the dynamic log stream and partially
 * reconstruct the original log message.
 */
struct PrintFragment
{
  // The type of the argument to pull from the dynamic buffer to the
  // partial format string (formatFragment)
  uint8_t argType : 5;

  // Indicates that the fragment requires a dynamic width/precision
  // argument in addition to one required by the format specifier.
  bool hasDynamicWidth : 1;
  bool hasDynamicPrecision : 1;

  // TODO(syang0) is this necessary? The format fragment is null-terminated
  // Length of the format fragment
  uint16_t fragmentLength;

  // A fragment of the original LOG statement that contains at most
  // one format specifier.
  char formatFragment[];
};

/**
 * These enums help encode LOG parameter types in the dynamic paramter
 * stream. These enums should match types generated by the preprocessor
 * and have a 1:1 correspondence to the actual type without the additional
 * underscores and "_t" at the end.
 */
enum FormatType : uint8_t
{
  NONE,

  unsigned_char_t,
  unsigned_short_int_t,
  unsigned_int_t,
  unsigned_long_int_t,
  unsigned_long_long_int_t,
  uintmax_t_t,
  size_t_t,
  wint_t_t,

  signed_char_t,
  short_int_t,
  int_t,
  long_int_t,
  long_long_int_t,
  intmax_t_t,
  ptrdiff_t_t,

  double_t,
  long_double_t,
  const_void_ptr_t,
  const_char_ptr_t,
  const_wchar_t_ptr_t,

  MAX_FORMAT_TYPE
};

class Nibbler
{
private:
  // Position in the nibble stream
  const TwoNibbles* nibblePosition;

  // Indicates whether whether to use the first nibble or second
  bool onFirstNibble;

  // Number of nibbles in this stream
  int numNibbles;

  // Position in the stream marking the next packed value
  const char* currPackedValue;

  // End of the last valid packed value
  const char* endOfValues;

public:
  /**
   * Nibbler Constructor
   *
   * \param nibbleStart
   *      Data stream consisting of the Nibbles followed by pack()ed values.
   * \param numNibbles
   *      Number of nibbles in the data stream
   */
  Nibbler(const char* nibbleStart, int numNibbles)
    : nibblePosition(reinterpret_cast<const TwoNibbles*>(nibbleStart))
    , onFirstNibble(true)
    , numNibbles(numNibbles)
    , currPackedValue(nibbleStart + (numNibbles + 1) / 2)
    , endOfValues(nullptr) {
    endOfValues = nibbleStart + (numNibbles + 1) / 2 + getSizeOfPackedValues(nibblePosition, numNibbles);
  }

  /**
   * Returns the next pack()-ed value in the stream
   *
   * \tparam T
   *      Type of the value in the stream
   * \return
   *      Next pack()-ed value in the stream
   */
  template<typename T>
  T getNext() {
    assert(currPackedValue < endOfValues);

    uint8_t nibble = (onFirstNibble) ? nibblePosition->first : nibblePosition->second;

    T ret = unpack<T>(&currPackedValue, nibble);

    if (!onFirstNibble) ++nibblePosition;

    onFirstNibble = !onFirstNibble;

    return ret;
  }

  /**
   * Returns a pointer to to the first byte beyond the last pack()-ed value
   *
   * \return
   *      Pointer to the first byte beyond last pack()-ed value.
   */
  const char* getEndOfPackedArguments() { return endOfValues; }
};

static FormatType getFormatType(std::string length, char specifier) {

  // Signed Integers
  if (specifier == 'd' || specifier == 'i') {
    if (length.empty()) return int_t;

    if (length.size() == 2) {
      if (length[0] == 'h') return signed_char_t;
      if (length[0] == 'l') return long_long_int_t;
    }

    switch (length[0]) {
      case 'h': return short_int_t;
      case 'l': return long_int_t;
      case 'j': return intmax_t_t;
      case 'z': return size_t_t;
      case 't': return ptrdiff_t_t;
      default: break;
    }
  }

  // Unsigned integers
  if (specifier == 'u' || specifier == 'o' || specifier == 'x' || specifier == 'X') {
    if (length.empty()) return unsigned_int_t;

    if (length.size() == 2) {
      if (length[0] == 'h') return unsigned_char_t;
      if (length[0] == 'l') return unsigned_long_long_int_t;
    }

    switch (length[0]) {
      case 'h': return unsigned_short_int_t;
      case 'l': return unsigned_long_int_t;
      case 'j': return uintmax_t_t;
      case 'z': return size_t_t;
      case 't': return ptrdiff_t_t;
      default: break;
    }
  }

  // Strings
  if (specifier == 's') {
    if (length.empty()) return const_char_ptr_t;
    if (length[0] == 'l') return const_wchar_t_ptr_t;
  }

  // Pointer
  if (specifier == 'p') {
    if (length.empty()) return const_void_ptr_t;
  }

  // Floating points
  if (specifier == 'f' || specifier == 'F' || specifier == 'e' || specifier == 'E' || specifier == 'g' ||
      specifier == 'G' || specifier == 'a' || specifier == 'A') {
    if (length.size() == 1 && length[0] == 'L')
      return long_double_t;
    else
      return double_t;
  }

  if (specifier == 'c') {
    if (length.empty()) return int_t;
    if (length[0] == 'l') return wint_t_t;
  }

  fprintf(stderr, "Attempt to decode format specifier failed: %s%c\r\n", length.c_str(), specifier);
  return MAX_FORMAT_TYPE;
}

bool StaticLogInfo::createFragments(char** microCode) {
  char* microCodeStartingPos = *microCode;

  fragments = microCodeStartingPos;
  size_t formatStringLength = strlen(formatString) + 1;
  std::regex regex("^%"
                   "([-+ #0]+)?"            // Flags (Position 1)
                   "([\\d]+|\\*)?"          // Width (Position 2)
                   "(\\.(\\d+|\\*))?"       // Precision (Position 4; 3 includes '.')
                   "(hh|h|l|ll|j|z|Z|t|L)?" // Length (Position 5)
                   "([diuoxXfFeEgGaAcspn])" // Specifier (Position 6)
  );

  size_t i = 0;
  std::cmatch match;
  int consecutivePercents = 0;
  size_t startOfNextFragment = 0;
  PrintFragment* pf = nullptr;

  // The key idea here is to split up the format string in to fragments (i.e.
  // PrintFragments) such that there is at most one specifier per fragment.
  // This then allows the decompressor later to consume one argument at a
  // time and print the fragment (vs. buffering all the arguments first).

  while (i < formatStringLength) {
    char c = formatString[i];

    // Skip the next character if there's an escape
    if (c == '\\') {
      i += 2;
      continue;
    }

    if (c != '%') {
      ++i;
      consecutivePercents = 0;
      continue;
    }

    // If there's an even number of '%'s, then it's a comment
    if (++consecutivePercents % 2 == 0 || !std::regex_search(formatString + i, match, regex)) {
      ++i;
      continue;
    }

    // Advance the pointer to the end of the specifier & reset the % counter
    consecutivePercents = 0;
    i += match.length();

    // At this point we found a match, let's start analyzing it
    pf = reinterpret_cast<PrintFragment*>(*microCode);
    *microCode += sizeof(PrintFragment);

    std::string width = match[2].str();
    std::string precision = match[4].str();
    std::string length = match[5].str();
    char specifier = match[6].str()[0];

    FormatType type = getFormatType(length, specifier);
    if (type == MAX_FORMAT_TYPE) {
      fprintf(stderr, "Error: Couldn't process this: %s\r\n", match.str().c_str());
      *microCode = microCodeStartingPos;
      return false;
    }

    pf->argType = 0x1F & type;
    pf->hasDynamicWidth = (width.empty()) ? false : width[0] == '*';
    pf->hasDynamicPrecision = (precision.empty()) ? false : precision[0] == '*';

    // Tricky tricky: We null-terminate the fragment by copying 1
    // extra byte and then setting it to NULL
    pf->fragmentLength = static_cast<uint16_t>(i - startOfNextFragment + 1);
    memcpy(*microCode, formatString + startOfNextFragment, pf->fragmentLength);
    *microCode += pf->fragmentLength;
    *(*microCode - 1) = '\0';

    startOfNextFragment = i;
    ++numPrintFragments;
  }

  // If we didn't encounter any specifiers, make one for a basic string
  if (pf == nullptr) {
    pf = reinterpret_cast<PrintFragment*>(*microCode);
    *microCode += sizeof(PrintFragment);
    numPrintFragments = 1;

    pf->argType = FormatType::NONE;
    pf->hasDynamicWidth = pf->hasDynamicPrecision = false;
    pf->fragmentLength = static_cast<uint16_t>(formatStringLength);
    memcpy(*microCode, formatString, formatStringLength);
    *microCode += formatStringLength;
  }
  else {
    // Extend the last fragment to include the rest of the string
    size_t endingLength = formatStringLength - startOfNextFragment;
    memcpy(pf->formatFragment + pf->fragmentLength - 1, // -1 to erase \0
           formatString + startOfNextFragment, endingLength);
    pf->fragmentLength = static_cast<uint16_t>(pf->fragmentLength - 1 + endingLength);
    *microCode += endingLength;
  }
  return true;
}

// RuntimeLogger constructor
RuntimeLogger::RuntimeLogger()
  : threadBuffers()
  , bufferMutex()
  , outputFp(nullptr)
  , currentLogLevel(INFO)
  , registrationMutex()
  , rawFragments(nullptr)
  , rawFragmentEnd(nullptr)
  , logCB(nullptr)
  , maxCBLogLevel(WARNING)
  , minCBPeriod(0) {

  rawFragmentEnd = rawFragments = static_cast<char*>(malloc(1024 * 1024));
}

// RuntimeLogger destructor
RuntimeLogger::~RuntimeLogger() {

  if (rawFragments) {
    free(rawFragments);
    rawFragments = nullptr;
  }

  if (outputFp) {
    fclose(outputFp);
    outputFp = nullptr;
  }
}

// See documentation in NanoLog.h
void RuntimeLogger::preallocate() {
  nanoLogSingleton.ensureStagingBufferAllocated();
  // I wonder if it'll be a good idea to update minFreeSpace as well since
  // the user is already willing to invoke this up front cost.
}

template<typename T>
static inline void printSingleArg(char** output, const char* formatString, T arg, int width = -1, int precision = -1) {

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-security"
#pragma GCC diagnostic ignored "-Wformat-nonliteral"

  if (width < 0 && precision < 0) {
    *output += sprintf(*output, formatString, arg);
  }
  else if (width >= 0 && precision < 0)
    *output += sprintf(*output, formatString, width, arg);
  else if (width >= 0 && precision >= 0)
    *output += sprintf(*output, formatString, width, precision, arg);
  else
    *output += sprintf(*output, formatString, precision, arg);

#pragma GCC diagnostic pop
}

void RuntimeLogger::poll_() {
  char compress_buf[1024 * 1024];
  char output_buf[1024 * 1024];

  {
    std::unique_lock<std::mutex> lock(registrationMutex);
    for (size_t i = bgInvocationSites.size(); i < invocationSites.size(); i++) {
      invocationSites[i].createFragments(&rawFragmentEnd);
      bgInvocationSites.push_back(invocationSites[i]);
    }
  }

  {
    std::unique_lock<std::mutex> lock(bufferMutex);
    for (size_t i = 0; i < threadBuffers.size(); i++) {
      StagingBuffer* sb = threadBuffers[i];
      while (true) {
        auto header = sb->peek();
        if (!header) {
          // If there's no work, check if we're supposed to delete
          // the stagingBuffer
          if (sb->checkCanDelete()) {
            delete sb;

            std::swap(threadBuffers[i], threadBuffers.back());
            threadBuffers.pop_back();
            i--;
          }
          break;
        }

        if (header->userdata >= bgInvocationSites.size()) break;
        StaticLogInfo& info = bgInvocationSites[header->userdata];
        header++;
        uint64_t timestamp = *(uint64_t*)(header);
        header++;
        char* argData = (char*)header;
        char* writePos = compress_buf;
        info.compressionFunction(info.numNibbles, info.paramTypes, &argData, &writePos);

        Nibbler nb(compress_buf, info.numNibbles);
        const char* nextStringArg = nb.getEndOfPackedArguments();
        char* output = output_buf;

        // Print out the actual log message, piece by piece
        PrintFragment* pf = reinterpret_cast<PrintFragment*>(info.fragments);
        for (int i = 0; i < info.numPrintFragments; ++i) {
          const wchar_t* wstrArg;

          int width = -1;
          if (pf->hasDynamicWidth) width = nb.getNext<int>();

          int precision = -1;
          if (pf->hasDynamicPrecision) precision = nb.getNext<int>();

          switch (pf->argType) {
            case NONE:

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wformat-security"
#pragma GCC diagnostic ignored "-Wformat-nonliteral"
              output += sprintf(output, pf->formatFragment);
#pragma GCC diagnostic pop
              break;

            case unsigned_char_t:
              printSingleArg(&output, pf->formatFragment, nb.getNext<unsigned char>(), width, precision);
              break;

            case unsigned_short_int_t:
              printSingleArg(&output, pf->formatFragment, nb.getNext<unsigned short int>(), width, precision);
              break;

            case unsigned_int_t:
              printSingleArg(&output, pf->formatFragment, nb.getNext<unsigned int>(), width, precision);
              break;

            case unsigned_long_int_t:
              printSingleArg(&output, pf->formatFragment, nb.getNext<unsigned long int>(), width, precision);
              break;

            case unsigned_long_long_int_t:
              printSingleArg(&output, pf->formatFragment, nb.getNext<unsigned long long int>(), width, precision);
              break;

            case uintmax_t_t:
              printSingleArg(&output, pf->formatFragment, nb.getNext<uintmax_t>(), width, precision);
              break;

            case size_t_t: printSingleArg(&output, pf->formatFragment, nb.getNext<size_t>(), width, precision); break;

            case wint_t_t: printSingleArg(&output, pf->formatFragment, nb.getNext<wint_t>(), width, precision); break;

            case signed_char_t:
              printSingleArg(&output, pf->formatFragment, nb.getNext<signed char>(), width, precision);
              break;

            case short_int_t:
              printSingleArg(&output, pf->formatFragment, nb.getNext<short int>(), width, precision);
              break;

            case int_t: printSingleArg(&output, pf->formatFragment, nb.getNext<int>(), width, precision); break;

            case long_int_t:
              printSingleArg(&output, pf->formatFragment, nb.getNext<long int>(), width, precision);
              break;

            case long_long_int_t:
              printSingleArg(&output, pf->formatFragment, nb.getNext<long long int>(), width, precision);
              break;

            case intmax_t_t:
              printSingleArg(&output, pf->formatFragment, nb.getNext<intmax_t>(), width, precision);
              break;

            case ptrdiff_t_t:
              printSingleArg(&output, pf->formatFragment, nb.getNext<ptrdiff_t>(), width, precision);
              break;

            case double_t: printSingleArg(&output, pf->formatFragment, nb.getNext<double>(), width, precision); break;

            case long_double_t:
              printSingleArg(&output, pf->formatFragment, nb.getNext<long double>(), width, precision);
              break;

            case const_void_ptr_t:
              printSingleArg(&output, pf->formatFragment, nb.getNext<const void*>(), width, precision);
              break;

            // The next two are strings, so handle it accordingly.
            case const_char_ptr_t:
              printSingleArg(&output, pf->formatFragment, nextStringArg, width, precision);

              nextStringArg += strlen(nextStringArg) + 1; // +1 for NULL
              break;

            case const_wchar_t_ptr_t:

              /**
               * I've occasionally encountered the following assertion:
               * __wcsrtombs: Assertion `data.__outbuf[-1] == '\0'' failed
               *
               * I don't know why this occurs, but it appears to be caused
               * by a wcslen() call deep inside printf returning the wrong
               * value when called on the dynamic buffer. I've found
               * that copying the data into a stack buffer first fixes
               * the problem (not implemented here).
               *
               * I don't know why this is the case and I'm inclined to
               * believe it's a problem with the library because I've...
               *  (1) verified byte by byte that the copied wchar_t
               *      strings and the surrounding bytes are exactly the
               *      same in the dynamic and stack allocated buffers.
               *  (2) copied the wcslen code from the glib sources into
               *      to this file and calling it on the dynamic buffer
               *      works, but the public wcslen() API still returns an
               *      incorrect value.
               *  (3) verified that no corruption occurs in the buffers
               *      before and after the wcslen() returns the wrong val
               *
               * If wide character support becomes important and this
               * assertion keeps erroring out, I would work around it
               * by copying the wide string into a stack buffer before
               * passing it to printf.
               */
              wstrArg = reinterpret_cast<const wchar_t*>(nextStringArg);
              printSingleArg(&output, pf->formatFragment, wstrArg, width, precision);
              // +1 for NULL
              nextStringArg += (wcslen(wstrArg) + 1) * sizeof(wchar_t);
              break;

            case MAX_FORMAT_TYPE:
            default: output += sprintf(output, "Error: Corrupt log header in header file\r\n"); exit(-1);
          }

          pf =
            reinterpret_cast<PrintFragment*>(reinterpret_cast<char*>(pf) + pf->fragmentLength + sizeof(PrintFragment));
        }

        static const char* logLevelNames[] = {"NAN", "ERR", "WRN", "INF", "DBG"};
        const char* logLevel = logLevelNames[info.severity];
        uint64_t ns = tscns.tsc2ns(timestamp);
        uint64_t t = (ns - midnight_ns) / 1000;
        uint32_t us = t % 1000000;
        t /= 1000000;
        uint32_t s = t % 60;
        t /= 60;
        uint32_t m = t % 60;
        t /= 60;
        uint32_t h = t % 24;
        if (outputFp) {
          fprintf(outputFp, "%02d:%02d:%02d.%06d %s:%d %s[%s]: %s\n", h, m, s, us, info.filename, info.lineNum,
                  logLevel, sb->getName(), output_buf);
        }

        if (logCB && info.severity <= maxCBLogLevel && info.last_ts + minCBPeriod < ns) {
          info.last_ts = ns;
          logCB(ns, info.severity, output_buf, output - output_buf);
        }

        sb->consume();
      }
    }
  }
  if (outputFp) {
    fflush(outputFp);
  }
}

void RuntimeLogger::setLogFile_internal(const char* filename) {
  FILE* newFp = fopen(filename, "a");
  if (!newFp) {
    std::string err = "Unable to open file new log file: '";
    err.append(filename);
    err.append("': ");
    err.append(strerror(errno));
    throw std::ios_base::failure(err);
  }

  if (outputFp) fclose(outputFp);
  outputFp = newFp;
}

/**
 * Set where the NanoLog should output its compressed log. If a previous
 * log file was specified, NanoLog will attempt to sync() the remaining log
 * entries before swapping files. For best practices, the output file shall
 * be set before the first invocation to log by the main thread as this
 * function is *not* thread safe.
 *
 * By default, the NanoLog will output to /tmp/compressedLog
 *
 * \param filename
 *      File for NanoLog to output the compress log
 *
 * \throw is_base::failure
 *      if the file cannot be opened or crated
 */
void RuntimeLogger::setLogFile(const char* filename) {
  nanoLogSingleton.setLogFile_internal(filename);
}

/**
 * Sets the minimum log level new NANO_LOG messages will have to meet before
 * they are saved. Anything lower will be dropped.
 *
 * \param logLevel
 *      LogLevel enum that specifies the minimum log level.
 */
void RuntimeLogger::setLogLevel(LogLevel logLevel) {
  if (logLevel < 0)
    logLevel = static_cast<LogLevel>(0);
  else if (logLevel >= NUM_LOG_LEVELS)
    logLevel = static_cast<LogLevel>(NUM_LOG_LEVELS - 1);
  nanoLogSingleton.currentLogLevel = logLevel;
}

}; // namespace NanoLogInternal

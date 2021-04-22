# NanoLogLite
NanoLogLite is a revised version of [NanoLog](https://github.com/PlatformLab/NanoLog) which is easier to use without performance compromise.

## The major changes are:
* NanoLogLite writes directly in human readable format instead of writing to a binary file and needing an additional decoder for reading, this behaves like other log libraries.
* User can now register a log msg callback function. This is useful in circumstances where warning/error msgs need to be published out for alerting.
* NanoLogLite won't create a default log file before user calling `NanoLog::setLogFile`.
* Log msg header is improved: NanoLogLite by default sets a meaningful thread id, and also allows user to customize the thread name; Time string is shortened, and folder part is removed from the file path.
* NanoLogLite won't create background thread internally, it requires user to poll it periodically. The idea is that user should have the ability to manage the threads in their program.  This would not make it harder to use as the user can simply create a thread himself to poll if he doesn't care.
* Log timestamp precision is much improved, it's almost synchronized with system time.
* Some shortcut macros are added for easier writing logging code: `logd`, `logi`, `logw`, `loge` for logging DEBUG/INFO/WARNING/ERROR msgs respectively.

## Examples of the new features:
```c++
// easier to write log msg
logi("Simple log message with 1 parameters, %d", i);

// set a meaningful name for current thread
NanoLog::setThreadName("main"); 

// set a callback/hook to capture msgs with WARNING or above log level
void logcb(uint64_t ns, NanoLog::LogLevel level, const char* msg, size_t msg_len) {
  printf("logcb, ns: %ld, msg: %s\n", ns, msg);
}
NanoLog::setLogCB(logcb, NanoLog::WARNING, 0);

// user need to poll it
NanoLog::poll(); 
```

![image](https://user-images.githubusercontent.com/11496526/115710553-1019f500-a3a5-11eb-8688-74d9bac60fa0.png)

## Test
Check [log_test.cc](https://github.com/MengRao/NanoLogLite/blob/main/test/log_test.cc) for example usage.

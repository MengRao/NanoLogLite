add_executable(log_test log_test.cc)
target_link_libraries(log_test nanologlite)
install(TARGETS log_test)

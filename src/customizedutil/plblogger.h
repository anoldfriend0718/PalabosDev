#pragma once

#include "io/parallelIO.h"
#include "plog/Severity.h"
#include <memory>
#include <plog/Init.h>
#include <plog/Log.h>
#include <plog/Severity.h>

namespace plb {

namespace util {

template <class Formatter>
class PlbLoggerAppender
    : public plog::IAppender // All appenders MUST inherit IAppender interface.
{
public:
  void
  write(const plog::Record &record) override // This is a method from IAppender
                                             // that MUST be implemented.
  {
    plog::util::nstring str = Formatter::format(
        record); // Use the formatter to get a string from a record.

    pcout << str;

    *logfile << str;
  }
  static std::shared_ptr<plb_ofstream> logfile;
};

template <class Formatter>
std::shared_ptr<plb_ofstream> PlbLoggerAppender<Formatter>::logfile =
    std::make_shared<plb_ofstream>("plb.log");

plog::Severity getSeverity(std::string severity) {
  if (severity == "verbose") {
    return plog::verbose;
  } else if (severity == "debug") {
    return plog::debug;
  } else if (severity == "info") {
    return plog::info;
  } else if (severity == "warning") {
    return plog::warning;
  } else if (severity == "error") {
    return plog::error;
  } else if (severity == "fatal") {
    return plog::verbose;
  } else if (severity == "none") {
    return plog::none;
  } else {
    return plog::info;
  }
}

} // namespace util

} // namespace plb
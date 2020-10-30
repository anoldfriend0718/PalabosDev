#pragma once

#include "io/parallelIO.h"
#include "plog/Severity.h"
#include <memory>
#include <plog/Init.h>
#include <plog/Log.h>

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

} // namespace util
} // namespace plb
//
// Created by Luis Yanes (EI) on 2019-06-26.
//

#include <sdglib/utilities/OutputLog.hpp>
#include "Journal.hpp"

JournalOperation::JournalOperation(const std::string &name, const std::string &tool, const std::string &detail) : name(name), timestamp(time(nullptr)), tool(tool), detail(detail) {}

void JournalOperation::addEntry(const std::string &detail) {
    entries.emplace_back(detail);
}

void JournalOperation::status() const {
    sdglib::OutputLog(false) << "Operation: " << name << " applied on " << ctime(&timestamp);
    sdglib::OutputLog(false) << "Tool: " << tool << std::endl;
    sdglib::OutputLog(false) << "Details: " << detail << std::endl << std::endl;
}

bool JournalOperation::operator==(const JournalOperation &o) const {
    return std::tie(name, tool, detail, timestamp, entries) == std::tie(o.name, o.tool, o.detail, o.timestamp, o.entries);
}

std::ostream &operator<<(std::ostream &os, const JournalOperation &opj) {
    os << "Operation: " << opj.name << " applied on " << ctime(&opj.timestamp);
    os << "Tool: " << opj.tool << std::endl;
    os << "Details: " << opj.detail << std::endl << std::endl;
    return os;
}

JournalEntry::JournalEntry(const std::string &detail) : detail(detail) {}

bool JournalEntry::operator==(const JournalEntry &o) const { return detail == o.detail;}

std::ostream &operator<<(std::ostream &os, const JournalEntry &je) {
    os << je.detail << std::endl;
    return os;
}

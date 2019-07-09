//
// Created by Luis Yanes (EI) on 2019-06-26.
//

#include <sdglib/logger/OutputLog.hpp>
#include "OperationJournal.hpp"

OperationJournal::OperationJournal(const std::string &name, const std::string &tool, const std::string &detail) : name(name), timestamp(time(nullptr)), tool(tool), detail(detail) {}

void OperationJournal::addEntry(const std::string &detail) {
    entries.emplace_back(detail);
}

void OperationJournal::status() const {
    sdglib::OutputLog(false) << "Operation: " << name << " applied on " << ctime(&timestamp);
    sdglib::OutputLog(false) << "Tool: " << tool << std::endl;
    sdglib::OutputLog(false) << "Details: " << detail << std::endl << std::endl;
}

bool OperationJournal::operator==(const OperationJournal &o) const {
    return std::tie(name, tool, detail, timestamp, entries) == std::tie(o.name, o.detail, o.detail, o.timestamp, o.entries);
}

JournalEntry::JournalEntry(const std::string &detail) : detail(detail) {}

bool JournalEntry::operator==(const JournalEntry &o) const { return detail == o.detail;}

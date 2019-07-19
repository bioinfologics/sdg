//
// Created by Luis Yanes (EI) on 2019-06-26.
//

#pragma once

#include <string>
#include <vector>
#include <tuple>

class JournalEntry {
public:
    JournalEntry()=default;
    explicit JournalEntry(const std::string &detail);
    std::string detail="";

    friend std::ostream& operator<<(std::ostream &os, const JournalEntry &je);

        bool operator==(const JournalEntry& o) const;
};

class JournalOperation {
public:
    JournalOperation() = default;
    JournalOperation(const std::string &name, const std::string &tool, const std::string &detail);

    friend std::ostream& operator<<(std::ostream &os, const JournalOperation &opj);

    void addEntry(const std::string &detail);

    void status() const;

    std::string name={};
    std::time_t timestamp=time(nullptr);
    std::string tool={};
    std::string detail={};
    std::vector<JournalEntry> entries = {};

    bool operator==(const JournalOperation& o) const;
};
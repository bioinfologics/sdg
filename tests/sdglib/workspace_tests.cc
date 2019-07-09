//
// Created by Luis Yanes (EI) on 2019-06-27.
//

#include <catch.hpp>
#include <sdglib/workspace/WorkSpace.hpp>

TEST_CASE("Test log/journal"){
    WorkSpace ws;
    ws.add_operation("Op1", "tool1", "Changed nothing");
    ws.add_operation("Op2", "tool2", "Changed nothing again");
    ws.status();

    REQUIRE(true);
}
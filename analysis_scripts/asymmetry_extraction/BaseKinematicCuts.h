#pragma once
#include <TTreeReader.h>

class BaseKinematicCuts {
public:
    BaseKinematicCuts(TTreeReader& reader) : reader(reader) {}
    virtual ~BaseKinematicCuts() {}

    virtual bool applyCuts(int currentFits, bool isMC) = 0; // Pure virtual function

protected:
    TTreeReader& reader;
};
#include <array>
#include <fstream>
#include <iostream>
#include <string>

// Invariables

// Variables

// Print flags
bool DETECTOR_VERBOSITY = 0;

// Utilities for parameters

enum Directions
{
    X = 0,
    Y,
    Z,
    DirectionSize
};

std::string AxisToString(int num)
{
    std::string name;

    switch (num)
    {
        case 0:
            name = "X";
            break;
        case 1:
            name = "Y";
            break;
        case 2:
            name = "Z";
            break;
        default:
            name = "X";
    }

    return name;
}

enum Signals
{
    Correlated = 0,
    Accidental,
    TotalDifference,
    SignalSize
};

std::string SignalToString(int num)
{
    std::string name;

    switch (num)
    {
        case 0:
            name = "Correlated";
            break;
        case 1:
            name = "Accidental";
            break;
        case 2:
            name = "Total Difference";
            break;
        default:
            name = "Total Difference";
    }

    return name;
}

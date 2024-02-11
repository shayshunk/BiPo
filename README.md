
<p align="center">
    <img src="PROSPECT.png" width="200">
</p>

---

<h1 align="center">
    <br>
    BiPo Directionality Study
    <br>
</h1>

## Calculating the average direction of the ${}^{214}$ Bi $\to {}^{214}$ Po decay in the PROSPECT detector. 

<h2>
    Summary
</h2>

The analysis code borrows heavily from Don Jones' work on BiPo in PROSPECT, available [here](https://github.com/jonesdc76/BiPoAnalysis). For the directionality analysis, I follow roughly the same method as the antineutrino directionality analysis. My repository for the neutrino directionality is [here](https://github.com/shayshunk/NeutrinoDirectionality). The process is to track the average displacement between prompt and delayed events in $x$, $y$, and $z$. The missing segments are accounted for with a modified method that is described in the paper. The purpose is to verify our directionality methods using BiPo (${}^{214}$ Bi $\to {}^{214}$ Po, seen in the ${}^{222}$ Rn decay chain) and extract any systematic uncertainties we might have. 

<h2>
    Requirements
</h2>

* ROOT 5+
* C++ compliler (I'm using g++ in the instructions)
* Analyzed and calibrated PROSPECT data, using the P2x BiPo plugin 

<h2>
    How To Use (Debian)
</h2>

```bash
# Clone this repository
git clone https://github.com/shayshunk/BiPo

# Go into the repository
cd BiPo

# Install dependencies (quick way below but not recommended)
# Go here to find the CERN recommended way: https://root.cern/install/
sudo snap install root-framework

# Edit the header NeutrinoDirectionality.h to point to your data files
# Line 26-29 contain my local paths but yours will be different

# Run the code
# I'm using g++ because it's straightforward but you can configure another compiler
# You can also run in macro mode but it'll be slower
g++ BiPo.cc -o BiPo `root-config --cflags --glibs`
./BiPo

# Make the plots
# Just use macro mode, it's fast enough that there's no time lost
# Compiling changes the plot aspect ratio for some reason
root BiPoPlots.cc
```

<h2>
Customization Options
</h2>

A few command line options are available for customizing the output of the code. The follwing variables in `BiPo.h` let you define how verbose you want the output of the values calculated to be.

```C++
DETECTOR_VERBOSITY 
BIPOCOUNT_VERBOSITY 
MEAN_VERBOSITY 
```
Each option is accessed by adding a flag after the executable while running in the terminal. Example: `./BiPo -D -B -M`. The options are are as follows:

 * `-D` will set `DETECTOR_VERBOSITY` to true. It will print the detector configuration used for the modified method. We are using the PRD configuration here because data splitting has not been applied to BiPo.
 * `-B` sets `IBD_COUNT VERBOSITY` to true. It will print the total and effective BiPo counts in each direction for each dataset. Effective IBDs are calculated through Poisson statistics. 
 * `-M` sets `MEAN_VERBOSITY` to true. It will print the *p* components and respective errors that are used to extract systematic uncertainty.

The other option is contained in `Formatting.h`. I added a few quick functions that return a certain formatting (bold/underline) or color for more aesthetically pleasing output. These only work on Linux terminals. If working on another platform or the output simply looks jumbled or unpleasant, turn off the special formatting on line 4 by setting it to 0.

```C++
#define COLORFUL_FORMATTING 1
```

---

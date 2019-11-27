# Generator macros for O2 simulations

Macros to test external generators in o2-sim.

Usage:
`o2-sim -n 10 -g extgen --extGenFile PATH_TO_MACRO/macroname.C`

## AMPT (amptgen.C)
* Needs 'ampt.dat' output file from AMPT generator
    * 'ampt.dat' needs to be in the same directory where o2-sim is run
    * Something wrong with the macro: crashes when number of event given to o2-sim is the same as number of event generated using AMPT
        * use one less event when running o2-sim (in '-n' argument), then it works

## Pythia8 (py8hiflow.C)
* Uses ROOT class TPythia8 so Pythia8 needs to be enabled in ROOT (check using `root-config --has-pythia8`)
* Can be enabled by adding following options to cmake command in aliBuild recipe for ROOT (*PATH_TO_ALICE_DIR/alidist/root.sh*):
    * `-DPYTHIA8_DIR=PATH_TO_PYTHIA8_INSTALLATION_DIR`
    * `-DPYTHIA8_INCLUDE_DIR=PATH_TO_PYTHIA8_INSTALLATION_DIR/include`
    * `-DPYTHIA8_LIBRARY=PATH_TO_PYTHIA8_INSTALLATION_DIR/lib/libpythia8.so`
    * `-Dpythia8=ON`
* Adding `pythia8` to the `FEATURES` variable checks that the feature is enabled
* Build O2 or ROOT again with these options

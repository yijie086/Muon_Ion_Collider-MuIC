# How to change and modify the Marlin Processors

## Files to be changed

### Calorimeter Processor

**Path:** /opt/ilcsoft/muonc/DDMarlinPandora/v00-13-MC

**File:** DDCaloDigi.h, DDCaloDigi.cc, DDCaloDigi_BIB.h, DDCaloDigi_BIB.cc

### Tracker Processor

**Path:** /opt/ilcsoft/muonc/MarlinTrkProcessors/v02-13-MC/source/Digitisers

**File:** DDPlanarDigiProcessor.h, DDPlanarDigiProcessor.cc

## Make Marlin with CMake

```
cd build
```

```
cmake -C $ILCSOFT/ILCSoft.cmake ..
```

```
make install
```

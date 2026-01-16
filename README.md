# StandaloneEEC Build Instructions

## Overview
StandaloneEEC is a C++ project for Energy-Energy Correlators (EEC) with multiple benchmarks, unit tests, and variants.

## Prerequisites
- CMake 3.16+
- boost
- C++17 compiler (g++ or clang)
- ROOT (with MathMore, GenVector, Minuit2 components)
- Eigen3
- libnpy (included as submodule)

## Configuration

1. **Install dependencies** (example for conda):
   ```bash
   conda install -c conda-forge root eigen boost
   ```

2. **Generate build files**:
   ```bash
   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
   ```

   For custom ROOT/Eigen paths:
   ```bash
   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
     -DROOT_DIR=/path/to/root \
     -DEigen3_DIR=/path/to/eigen3
   ```

## Building

**Standard build**:
```bash
cmake --build build
```

**Parallel build** (faster):
```bash
cmake --build build -j
```

**With verbose output**:
```bash
cmake --build build -j --verbose
```

## Cleaning

**Remove all build artifacts**:
```bash
rm -rf build
```

Then reconfigure and rebuild as above.

## Outputs

- **Libraries**: `lib/EEC.so`, `lib/EEC_byhand.so`, `lib/SimonTools.so`, `lib/Matching.so`
- **Executables**: `bin/res*`, `bin/CAres*`, `bin/proj_unittest`, `bin/test_matching_v2`, `bin/testshower`

## Running Tests

```bash
./bin/res3_unittest
./bin/res4_unittest
./bin/CAres3_unittest
./bin/CAres4_unittest
./bin/proj_unittest
```

## Notes

- Variants ending in `_checkbyhand` are compiled with `CHECK_BY_HAND` flag
- Release builds include `-O3 -mtune=native -flto` optimizations
- All warnings treated as errors (`-Werror`)

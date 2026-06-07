# LiteMesh: An Efficient Quality-Aware Algorithm for Feature-Preserving Triangle Mesh Simplification

homepage of LiteMesh: [link]

paper link of LiteMesh: [link]

## Dependency
- CGAL
- Eigen3
- C++ 17
- vcpkg
- CMake
 
## Build
We use vcpkg to install dependencies and use CMake to build the project. Please make sure they are installed on your computer.
```
git clone https://github.com/ochenglino/LiteMesh.git && cd LiteMesh

# install dependencies by vcpkg
$<VCPKG_ROOT>\vcpkg.exe install --triplet x64-windows

# configuration
cmake -B build -S . -DCMAKE_TOOLCHAIN_FILE="$<VCPKG_ROOT>\scripts\buildsystems\vcpkg.cmake" -A x64

# compilation
cmake --build build --config Release
```
The output LiteMesh.exe will be in the folder ``` /build/src/Release```

## Run 
After the project is successfully built, you can easily use it from the command line, as follows:
```
LiteCAD.exe <INPUT> <OUTPUT> <SIZE> <FACTOR_OF_DECAY> <FACTOR_FOR_PARTICLE_SYSTEM>
# example: LiteMesh.exe ../../../data/bunny.obj ../../../data 1000 0.95 0.3 
```
It is worth noting that we recommend setting the decay factor and particle system factor to 0.95 and 0.3 respectively if you want to achieve the desired mesh simplification effect.

## Citation
If this code contributes to academic work, please cite as:
```

```


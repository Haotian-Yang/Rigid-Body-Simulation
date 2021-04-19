## Introduction

This is my implementation of rigid body simulation assignment in [CSC417/CSC2549-Physics-based Animation](https://github.com/dilevin/CSC417-physics-based-animation).
![Fun with interactive rigid bodies](images/rb_example.gif)

## Build & Execution
```
git submodule update --init --recursive
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make 
./a5-cloth-simulation
```
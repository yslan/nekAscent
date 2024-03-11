# nekAscent

An CXX/Fortran interface to link Ascent for Nek5000

Original implementation by Andres Sewell and Victor Mateevitsi    

The code is initially copied from https://github.com/siramok/nekIBM-ascent
It's striped down and refactored for Nek5000.

## Build
It requires MPI and c++ compiler. To run examples, it needs a mpi-fortran compiler.    

```
git clone https://github.com/yslan/nekAscent.git
cd nekAscent

make ASCENT_DIR=~/ascent/install/ascent-develop CXX=mpic++ FC=mpif77 UNDERSCORE=1
```

## Environment Variables (WIP)

`NEKASCENT_VERBOSE_LEVEL`=0, 1 (default), 2, 3

At runtime, it might also need the following path. 
The version for vtkm and conduit depends on the version of Ascent.
```
export LD_LIBRARY_PATH+=${ASCENT_DIR}/../ascent-develop/lib:
export LD_LIBRARY_PATH+=${ASCENT_DIR}/../vtk-m-v2.1.0/lib:
export LD_LIBRARY_PATH+=${ASCENT_DIR}/../conduit-v0.8.8/lib:
```

## Run Example

```
cd build/examples/fortran/

# folder for images (path can be cahnged)
mkdir -p img

# sample of Ascent's yaml files
cp ../../../action_files/ascent_actions.yaml .
cp ../../../action_files/trigger_render.yaml .

# run the case
mpiexec -n 4 ./eddy2d |tee log
```

(TODO) img to mp4


## Portability (WIP)

- docker
  ```
  ASCENT_DIR="/ascent/install
  ```

## See Also (WIP)
- nekrs (PR)
- nek5000 (PR)
- nekrsAscentExamples
- nek5kAscentExamples

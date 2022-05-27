# Mesh flow Project
The mesh flow project for mesh simplification.

## To download and build

```
git clone --recurse-submodules https://github.com/csyzzkdcz/mesh_flow.git
cd mesh_flow
mkdir build
cd build
cmake ..
make -j4

```

## To run
Suppose you are under build folder, then
```
./bin/mesh-flow_bin ../testdata/0_vd.obj ../testdata/0_mc.obj
```
and click the Iso-flow button.

## issue
these codes doesn't work on Apple with M1 chips (the sign_distance was broken on M1 according to https://github.com/libigl/libigl/issues/1686#issuecomment-956293758), I wrappered a subversion of signed distance of libigl which doesn't use any winding number function

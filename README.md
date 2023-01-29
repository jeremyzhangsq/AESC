# AESC

## Environment
- System: Ubuntu 18.04
- GCC: 7.5.0
- CMake: 2.8.12
- Python 3.x
- Boost: 1.65.1

## Running Commands
Constructing xxx
```
$ cd datasets
$ python calEigen.py -f facebook/ -w 128
```
AESC approximation
```
$ cd aesc
$ cmake .
$ make clean all -j16
# run TGT with epsilon 0.05 and omega 128 on FB
$ ./aesc -f ../datasets/ -g facebook -a tgt -e 0.05 -w 128       
# run TGT+ with epsilon 0.05 and omega 128 on FB
$ ./aesc -f ../datasets/ -g facebook -a tgt+ -e 0.05 -w 128
```
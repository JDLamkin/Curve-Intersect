# Curve-Intersect

Research done by Jordan Lamkin to use the [SLV library](http://www-polsys.lip6.fr/˜elias/soft.html) for 2D curve intersection. This is meant to be
1. A test to determine how much of an improvement SLV gives over the previous state-of-the-art, Sturm sequences
2. A stepping stone on the way to creating a full 3D surface toolset for intersections and boolean operations

See [Efficient Polynomial Root Isolation Applied to Computational Geometry](https://oaktrust.library.tamu.edu/bitstream/handle/1969.1/186132/JordanLamkin_EH-CSCE_honors_thesis.pdf?sequence=1&isAllowed=y) for more details and for results.

## Usage

This software is packaged in two parts: The original root isolation library (`slv`) and the curve intersection library (`cci`). The single top-level Makefile will generate a test suite for each library. Both libraries require [GMP](https://gmplib.org/) and each Makefile has variables `GMP_INC` and `GMP_LIB` to point to the installation thereof.

### Quick Run
After installing GMP, the following commands should work and give a basic overview of how to use `slv` and `cci`:
```bash
make
./slv -h
./slv -f slv-v0.5/test.dat -i 20 -p
./cci < cci-v0.0/cciTest.dat
```

### SLV

SLV isolates roots of a given integer polynomial, input as a file in the form `[degree] [coef x^0] [coef x^1] ...`. For example, if `test.dat` contains
```
5 -120 600 -600 200 -25 1
```
then it represents the polynomial -120 + 600x - 600x<sup>2</sup> + 200x<sup>3</sup> - 25x<sup>4</sup> + x<sup>5</sup>.

According to SLV's own help message:
```
The syntax is:
slv [-h] [-f file] [-i prec] [-p]
Details:
-h : prints help message
-f fname : read the coeffs from the fname
-i prec : the output intervals have width 2^(-prec)
-p : print roots (Default: no print)
```

So running `./slv -f slv-v0.5/test.dat -i 20 -p` should yeild something like
```
Isolating test.dat
Finish isolating
/*------------------------------------------------------------*/
Statistics:
#degree : 5
           Neg     Pos   Total
#bound :    -1;     6
#roots :     0     5     5

#nodes :     0     10    10
#depth :     0     5     5
#trans :     0     5     5
#homo :      0     9     9
#pos_h_1 :   0     2;     2
#pos_h_2 :   0     3     3
#half_h :    0     1     1
/*------------------------------------------------------------*/
Refinement process...
#roots = 5
[
 -1; (276363/2^20, 276364/2^20)
  1 (1482060/2^20, 1482061/2^20)
 -1; (3771125/2^20, 3771126/2^20)
  1 (7430010/2^20, 7430011/2^20)
 -1; (13254840/2^20, 13254841/2^20)
]

Solving time: 0.000100 s
Refinement time: 0.000103 s
Total time: 0.000203 s
```

### CCI

CCI finds intersection points between curves defined by 2 variable polynomials. It accepts input of two rectangular polynomials each in the form
```
[x degree] [y degree]

[coef x^0 y^0] [coef x^1 y^0] ...
[coef x^0 y^1] [coef x^1 y^1] ...
       .              .       .  
       .              .        . 
       .              .         .
```

So if the input to `./cci` is
```
2 2

3   16  -16
24  0   0
-16 0   0

4 4

-18 12  9   -8  2
12  0   -8  0   0
9   -8  4   0   0
-8  0   0   0   0
2   0   0   0   0
```

then the output should be something like
```
-16x^{2} + 16x - 16y^{2} + 24y + 3 = 0
2x^{4} - 8x^{3} + 4x^{2}y^{2} - 8x^{2}y + 9x^{2} - 8xy^{2} + 12x + 2y^{4} - 8y^{3} + 9y^{2} + 12y - 18 = 0
max(2^{0} abs(x - -1/2), 2^{1} abs(y - 5/4)) <= 0.5
(-851/2048, 2359/2048)
max(2^{0} abs(x - 1/2), 2^{1} abs(y - 7/4)) <= 0.5
(1247/2048, 3571/2048)
max(2^{32} abs(x - 10508069981/8589934592), 2^{32} abs(y - 510812781/8589934592)) <= 0.5
(10508069981/8589934592, 510812781/8589934592)
max(2^{2} abs(x - 11/8), 2^{1} abs(y - 3/4)) <= 0.5
(3067/2048, 1389/2048)
0.00310600000000000 0.00841000000000000
```

representing the two input polynomials -16x<sup>2</sup> + 16x - 16y<sup>2</sup> + 24y + 3 and 2x<sup>4</sup> - 8x<sup>3</sup> + 4x<sup>2</sup>y<sup>2</sup> - 8x<sup>2</sup>y + 9x<sup>2</sup> - 8xy<sup>2</sup> + 12x + 2y<sup>4</sup> - 8y<sup>3</sup> + 9y<sup>2</sup> + 12y - 18, and their intersection regions and estimated points. This output format is meant to be easily pastable into [Desmos](https://www.desmos.com/calculator).

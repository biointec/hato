# Hato

Heuristic AlgoriThms and Optimality algorithms for the design of search schemes for approximate pattern matching. 
Hato was introduced in our paper  "Automated design of efficient search schemes for lossless approximate pattern matching"[^Renders2023].

**Generating optimal search schemes!**

This program implements a search scheme solver that can design search schemes either greedily or optimally under our objective function. 
It also provides functionality to calculate the expected number of enumerated substrings during the approximate pattern matching of a random pattern of length `m` to a random text of length `n` over an alphabet of size `sigma` using the Hamming distance and the provided/calculated search scheme[^Kucherov2015].

Hato makes use of an ILP formulation that is solved using CPLEX. Both this ILP formulation and a greedy heuristic are described in our aforementioned paper [^Renders2023]. The ILP formulation is inspired by the work by Kianfar et al. [^Kianfar2017].

## Usage

To use the program, run the executable with the following command-line arguments:
```
$ ./hato <mode> -k <number> [-s <path>] [-c <critical>] [-g] [-n <number>] [-m <number>] [-sigma <number>]
```

Modes:

  *  solver: Adapt the search scheme using either the greedy or optimal algorithm.
  *  expNodes: Report the number of expected nodes using the Hamming distance.

Arguments:

   * -k <number>: The number of allowed errors.
   * -s <path>: Path to the search scheme. (optional, if not set a pigeonhole principle based scheme will be used)
   * -c <critical>: The critical search index (only if k is even and p = k + 1).
  *  -g: Flag to indicate to solve using only the greedy algorithm (optional).
  *  -n <number>: Length of the reference text for the expected number of nodes. (optional, default: 3000000000).
  *  -m <number>: Length of the pattern for the expected number of nodes calculation(optional, default: 150).
  * -sigma <number>: Size of the alphabet for the expected number of nodes calculation(optional, default: 4).
 
### Search Scheme Notation
 
 The search scheme file  must adhere to the following format:

- Each line represents a single search definition, with three lists enclosed in curly braces. The lists are separated by whitespace.
- Each list should contain a comma-separated list of integers.
- The three sets must represent the following:
  1. The permutation (`pi`) of the parts.
  2. The lower bounds (`L`) for each part.
  3. The upper bounds (`U`) for each part.
- The sets should not have any spaces within the curly braces.

Example of a valid search definition line:

```
{0,1,2,3} {0,0,0,1} {0,1,3,3}
```



#### Example
Consider file `search_scheme.txt` for `k = 2` with the following content:

```
{0,1,2,3} {0,0,0,2} {0,1,2,2}
{1,2,3,0} {0,0,1,1} {0,1,1,2}
{2,3,1,0} {0,0,0,0} {0,0,2,2}
```
This file adheres to the format.

## Dependencies

    * C++ compiler with C++11 support
    * CMake 3.0 or higher
    * IBM ILOG CPLEX Optimization Studio 22.11 or higher

Before building the program, make sure to set the CPLEX_ROOT_DIR and CPLEX_DIR variables in the CMakeLists.txt file to the appropriate paths of your CPLEX installation.

## Building the program

To build the program, follow these steps:

    1. Create a build directory and navigate to this directory: mkdir build && cd build
    2. Configure the project: cmake ..
    3. Build the project: make

The executable `hato` will be generated.
 
 ## Generating Search Schemes
 
 To generate the minU search schemes[^Renders2023], use the following command:
 ```
$ ./hato sover -k <number> 
```
In order to generate the co-optimal search schemes for even `k`, one needs to add the `-c` flag together with the index of the requested critical search.
 
To generate the greedy adaptations of the pigeonhole principle use:
 
  ```
$ ./hato sover -k <number> -g
```
If you want to greedily adapt another scheme, or use this scheme as the basis for the ILP, you can use the `-s` flag to let the program know what search scheme to use.

[^Renders2023]: Renders, L., Depuydt, L., Rahmann, S., Fostier, J. (2023). "Automated design of efficient
search schemes for lossless approximate pattern matching" submitted to RECOMB.
[^Kucherov2015]: Kucherov, G., Salikhov, K., Tsur, D. (2015). "Approximate String Matching Using a Bidirectional Index." In A. S. Kulikov, S. O. Kuznetsov, & P. Pevzner (Eds.), *Combinatorial Pattern Matching* (pp. 222-231). Springer International Publishing.
[^Kianfar2017]: Kianfar, K., Pockrandt, C., Torkamandi, B., Luo, H., & Reinert, K. (2017). "Optimum Search Schemes for approximate string matching using bidirectional FM-index." arXiv preprint arXiv:1711.02035.


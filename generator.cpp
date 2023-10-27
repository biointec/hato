/******************************************************************************
 *  Hato: Heuristic AlgoriThms and Optimality algorithms for the design of    *
 * search schemes                                                             *
 *  Copyright (C) 2023      - Luca Renders <luca.renders@ugent.be> and        *
 *                            Jan Fostier <jan.fostier@ugent.be> and          *
 *                            Sven Rahmann <rahmann@cs.uni-saarland.de>       *
 *                                                                            *
 *  This program is free software: you can redistribute it and/or modify      *
 *  it under the terms of the GNU Affero General Public License as            *
 *  published by the Free Software Foundation, either version 3 of the        *
 *  License, or (at your option) any later version.                           *
 *                                                                            *
 *  This program is distributed in the hope that it will be useful,           *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU Affero General Public License for more details.                       *
 *                                                                            *
 * You should have received a copy of the GNU Affero General Public License   *
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.     *
 ******************************************************************************/

#include "generator.h"
#include <iostream>

using namespace std;

uint64_t nChoosek(uint16_t n, uint16_t k) {
    assert(n > k);
    double res = 1;
    for (uint16_t i = 1; i <= k; i++) {
        res *= (n - k + i) * 1.0 / i;
    }

    return (uint64_t)(res + 0.01);
}

bool successor(vector<uint16_t>& A, uint16_t m, uint16_t n) {
    assert(m < A.size());
    int8_t i = m - 1;

    while (A[i] == n - m + i) {
        i--;
        if (i < 0) {
            // undefined successor
            return false;
        }
    }

    for (uint16_t j = m; j-- > i;) {
        A[j] = A[i] + 1 + j - i;
    }
    return true;
}

vector<ErrorDistribution> generateErrorDistributions(const uint16_t k,
                                                     const uint16_t p) {

    // Create first p-subset
    vector<uint16_t> A(p);
    iota(A.begin(), A.end(), 0); // iota fill vector with increasing values

    // add k + p to A to handle edge case for final part
    // this removes the need for an if-clause in the loop where a distribution
    // is created from the p-subset
    A.emplace_back(k + p);

    // create return vector
    vector<ErrorDistribution> distributions;
    distributions.reserve(nChoosek(p + k, k));

    cout << distributions.capacity() << endl;

    uint64_t index = 0;

    // start generating all subsets and distributions
    while (true) {
        // generate the error distribution from the p-subset
        vector<uint16_t> distribution(p);

        for (uint16_t i = 0; i < p; i++) {

            distribution[i] = A[i + 1] - A[i] - 1;
        }

        distributions.emplace_back(distribution, index++);

        // go to next subset
        if (!successor(A, p, p + k)) {
            // no successor so end function
            return distributions;
        }
    }
}

ostream& operator<<(ostream& os, const Search& obj) {
    os << "{";
    uint32_t numParts = obj.getNumParts();
    for (uint32_t i = 0; i < numParts; i++) {
        os << obj.getPart(i) << ((i == numParts - 1) ? "}" : ",");
    }
    os << " {";
    for (uint32_t i = 0; i < numParts; i++) {
        os << obj.getLowerBound(i) << ((i == numParts - 1) ? "}" : ",");
    }
    os << " {";
    for (uint32_t i = 0; i < numParts; i++) {
        os << obj.getUpperBound(i) << ((i == numParts - 1) ? "}" : ",");
    }
    return os;
}
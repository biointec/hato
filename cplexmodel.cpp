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
#include "cplexmodel.h"

using namespace std;

void initializeBoolVarMatrix(uint16_t d1, uint16_t d2, const string& name,
                             IloEnv& env, BoolVarMatrix& matrix) {

    for (uint16_t i = 0; i < d1; i++) {
        matrix[i] = IloBoolVarArray(env, d2);
        for (uint16_t j = 0; j < d2; j++) {
            matrix[i][j] = IloBoolVar(
                env, (name + "_" + to_string(i) + "_" + to_string(j)).c_str());
        }
    }
}

void initializeIntVarMatrix(uint16_t d1, uint16_t d2, uint16_t max,
                            const string& name, IloEnv& env,
                            IntVarMatrix& matrix) {

    for (uint16_t i = 0; i < d1; i++) {
        matrix[i] = IloIntVarArray(env, d2);
        for (uint16_t j = 0; j < d2; j++) {
            matrix[i][j] = IloIntVar(
                env, 0, max,
                (name + "_" + to_string(i) + "_" + to_string(j)).c_str());
        }
    }
}

void initializeIntVar3Matrix(uint16_t d1, uint16_t d2, uint16_t d3,
                             uint16_t max, const string& name, IloEnv& env,
                             IntVar3Matrix& matrix) {
    for (int i = 0; i < d1; i++) {
        matrix[i] = IntVarMatrix(env, d2);
        for (int j = 0; j < d2; j++) {
            matrix[i][j] = IloIntVarArray(env, d3);

            for (int k = 0; k < d3; k++) {
                matrix[i][j][k] = IloIntVar(env, 0, max,
                                            (name + "_" + to_string(i) + "_" +
                                             to_string(j) + "_" + to_string(k))
                                                .c_str());
            }
        }
    }
}

void initializeBoolVar3Matrix(uint16_t d1, uint16_t d2, uint16_t d3,
                              const string& name, IloEnv& env,
                              BoolVar3Matrix& matrix) {
    for (int i = 0; i < d1; i++) {
        matrix[i] = BoolVarMatrix(env, d2);
        for (int j = 0; j < d2; j++) {
            matrix[i][j] = IloBoolVarArray(env, d3);

            for (int k = 0; k < d3; k++) {
                matrix[i][j][k] =
                    IloBoolVar(env, (name + "_" + to_string(i) + "_" +
                                     to_string(j) + "_" + to_string(k))
                                        .c_str());
            }
        }
    }
}
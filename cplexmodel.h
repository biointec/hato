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
#ifndef CPLEXMODEL_H
#define CPLEXMODEL_H

#include "generator.h"
#include <cmath>
#include <ilconcert/ilolinear.h>
#include <ilcplex/ilocplex.h>

// ============================================================================
// VARIABLE INTIALIZATION HELPERS
// ============================================================================

typedef IloArray<IloIntVarArray> IntVarMatrix;
typedef IloArray<IntVarMatrix> IntVar3Matrix;

typedef IloArray<IloBoolVarArray> BoolVarMatrix;
typedef IloArray<BoolVarMatrix> BoolVar3Matrix;

/**
 * @brief Initializes a two-dimensional matrix of integer variables.
 *
 * This function initializes a two-dimensional matrix of integer variables with
 * the specified dimensions and maximum value. The minimum value is zero. Each
 * variable in the matrix is assigned a unique name based on the provided name
 * prefix and its indices.
 *
 * @param d1 The number of rows in the matrix.
 * @param d2 The number of columns in the matrix.
 * @param max The maximum value allowed for each variable in the matrix.
 * @param name The prefix for variable names.
 * @param env The CPLEX environment.
 * @param matrix The matrix of integer variables to be initialized.
 */
void initializeIntVarMatrix(uint16_t d1, uint16_t d2, uint16_t max,
                            const std::string& name, IloEnv& env,
                            IntVarMatrix& matrix);
/**
 * @brief Initializes a three-dimensional matrix of Integer variables.
 *
 * This function initializes a three-dimensional matrix of Integer variables
 * with the specified dimensions and maximum value. The minimum value is 0. Each
 * variable in the matrix is assigned a unique name based on the provided name
 * prefix and its indices.
 *
 * @param d1 The number of elements in the first dimension.
 * @param d2 The number of elements in the second dimension.
 * @param d3 The number of elements in the third dimension.
 * @param max The maximum value allowed for each variable in the matrix.
 * @param name The prefix for variable names.
 * @param env The CPLEX environment.
 * @param matrix The matrix of integer variables to be initialized.
 */
void initializeIntVar3Matrix(uint16_t d1, uint16_t d2, uint16_t d3,
                             uint16_t max, const std::string& name, IloEnv& env,
                             IntVar3Matrix& matrix);

/**
 * @brief Initializes a two-dimensional matrix of boolean variables.
 *
 * This function initializes a two-dimensional matrix of boolean variables with
 * the specified dimensions. Each variable in the matrix is
 * assigned a unique name based on the provided name prefix and its indices.
 *
 * @param d1 The number of rows in the matrix.
 * @param d2 The number of columns in the matrix.
 * @param name The prefix for variable names.
 * @param env The CPLEX environment.
 * @param matrix The matrix of boolean variables to be initialized.
 */
void initializeBoolVarMatrix(uint16_t d1, uint16_t d2, const std::string& name,
                             IloEnv& env, BoolVarMatrix& matrix);

/**
 * @brief Initializes a three-dimensional matrix of boolean variables.
 *
 * This function initializes a three-dimensional matrix of boolean variables
 * with the specified dimensions. Each variable in the matrix is assigned a
 * unique name based on the provided name prefix and its indices.
 *
 * @param d1 The number of elements in the first dimension.
 * @param d2 The number of elements in the second dimension.
 * @param d3 The number of elements in the third dimension.
 * @param name The prefix for variable names.
 * @param env The CPLEX environment.
 * @param matrix The matrix of boolean variables to be initialized.
 */
void initializeBoolVar3Matrix(uint16_t d1, uint16_t d2, uint16_t d3,
                              const std::string& name, IloEnv& env,
                              BoolVar3Matrix& matrix);

// ============================================================================
// CLASS CPLEXMODEL
// ============================================================================
/**
 * @brief CPLEX-based model for solving the scheme adaptation problem.
 *
 * The CPLEXModel class provides an implementation of the scheme adaptation
 * problem using the CPLEX optimization solver. It constructs and solves an
 * optimization model to adapt the scheme by optimizing the coverage of
 * distributions using searches. The class initializes the model, sets up the
 * constraints and objectives, and provides functions to solve the model and
 * retrieve the adapted scheme.
 */
class CPLEXModel {

  private:
    IloEnv& env;       // CPLEX environment
    IloModel model;    // CPLEX model
    IloCplex cplex;    // CPLEX solver
    IloRangeArray con; // Array of constraints

    uint16_t k; // Maximal number of allowed errors
    uint16_t S; // the number of searchs
    uint64_t Q; // the number of error distributions
    uint16_t p; // the number of parts for the searches

    BoolVarMatrix lambda; // coverage per distribution per search
    IntVarMatrix U;       // upper bounds per search
    IntVarMatrix L;       // lower bounds per search

    IntVar3Matrix cumPerErrors;

    BoolVarMatrix y; // Binary encoding of pi-matrix

    AdaptableScheme& scheme;    // The current scheme
    uint16_t criticalSearchIdx; //  Index of the critical searc

    /**
     * @brief Initializes the lambda variable matrix and adds it to the
     * environment.
     *
     * This function initializes a matrix of binary decision variables called
     * lambda. Each element in the matrix represents a binary decision variable
     * lambda[q][s], where q is the index of a distribution and s is the index
     * of a search.
     */
    void initializeLambda() {
        initializeBoolVarMatrix(Q, S, "lambda", env, lambda);
    }

    /**
     * @brief Initializes the U variable matrix and adds it to the environment.
     *
     * This function initializes a matrix of integer decision variables called
     * U. Each element in the matrix represents an integer decision variable
     * U[s][i], where s is the index of a search and i is the index of a part.
     * The minimum value is 0 and maximum value is k.
     */
    void initializeUpper() {
        initializeIntVarMatrix(scheme.getSearches().size(),
                               scheme.getNumParts(), k, "U", env, U);
    }

    /**
     * @brief Initializes the L variable matrix and adds it to the environment.
     *
     * This function initializes a matrix of integer decision variables called
     * L. Each element in the matrix represents an integer decision variable
     * L[s][i], where s is the index of a search and i is the index of a part.
     * The minimum value is 0 and maximum value is k.
     */
    void initializeLower() {
        initializeIntVarMatrix(scheme.getSearches().size(),
                               scheme.getNumParts(), k, "L", env, L);
    }

    /**
     * @brief Initializes the cum errors encountered variables
     */
    void initializeCumPerError() {
        initializeIntVar3Matrix(
            scheme.getDistributions().size(), scheme.getSearches().size(),
            scheme.getNumParts(), k, "cumPerErrors", env, cumPerErrors);
    }

    /**
     * @brief Initializes the y variable matrix and adds it to the
     * environment.
     *
     */
    void initializeY() {
        initializeBoolVarMatrix(S, p, "y", env, y);
    }

    /**
     * @brief Creates the objective function for minimizing U.
     *
     * This function creates the objective function expression for
     * minimizing the value of U. The objective function is defined as the
     * sum over both searches and parts of U, where for each part i the
     * value is mutliplied with (p-i) squared.
     *
     * @return The objective expression for minimizing U.
     */
    IloExpr createObjectiveMinU() const {
        IloExpr obj(env);
        for (uint16_t s = 0; s < S; s++) {
            for (uint16_t i = 0; i < p; i++) {
                auto weight = getWeightsU(i);
                obj += weight * U[s][i];
            }
        }
        return obj;
    }

    IloExpr createObjectiveMinUMaxL() const {
        IloExpr obj(env);
        for (uint16_t s = 0; s < S; s++) {

            obj += getWeightsU(0) * U[s][0];
            for (uint16_t i = 1; i < p; i++) {

                obj += getWeightsU(i) * (U[s][i] - L[s][i - 1]);
            }
        }
        return obj;
    }

    double getWeightsU(uint16_t i) const {
        return pow(k + 1, p - i - 1);
    }

    /**
     * @brief Creates the objective function for maximizing L.
     *
     * This function creates the objective function expression for
     * minimizing the value of U. The objective function is defined as the
     * sum over both searches and parts of U, where for each part i the
     * value is mutliplied with (p-i) squared.
     *
     * @return The objective expression for maximizing L.
     */
    IloExpr createObjectiveMaxL() const {
        IloExpr obj(env);
        for (uint16_t s = 0; s < S; s++) {
            for (uint16_t i = 1; i < p; i++) {
                obj += pow(p - i, 1) * L[s][i];
            }
        }
        return -obj;
    }

    /**
     * @brief Creates the objective function for minimizing coverage.
     *
     * This function creates the objective function expression for
     * minimizing the coverage. The objective function is defined as the sum
     * of the coverage for each error distribution.
     *
     * @return The objective expression for minimizing coverage.
     */
    IloExpr createObjectiveMinCov() const {
        IloExpr obj(env);
        for (uint16_t q = 0; q < Q; q++) {
            for (uint16_t s = 0; s < S; s++) {
                obj += lambda[q][s];
            }
        }
        return obj;
    }

    /**
     * @brief Sets the minimal cover constraints.
     *
     * This function sets the minimal cover constraints for the optimization
     * problem. The constraints ensure that each distribution is covered by at
     * least one search.
     */
    void setMinimalCoverconstraints() {
        // create minimal cover constraint
        // each distribution must be covered by at least 1 search
        IloExpr totalCov(env);

        for (uint64_t q = 0; q < Q; q++) {
            IloExpr coverage(env);
            for (uint16_t s = 0; s < S; s++) {
                coverage += lambda[q][s];
            }
            con.add(coverage >= 1);
            coverage.end();
        }
    }

    /**
     * @brief Sets the constraints on the y variables.
     *
     * This function sets the constraints on the binary representations (y
     * variables) of the pi-strings. The function ensures that the orderings are
     * distinct and creates an ordering of the searches to remove symmetry from
     * the model. It also links y to the cumulative permutation errors that help
     * in setting the coverage variables correctly.
     */
    void setYConstraints() {
        // enforce distrinct permutations and create ordering
        for (uint16_t s = 0; s < S - 1; s++) {

            if (S == p) {
                // number of zeros in y (index of start part) must be equal to s
                IloExpr numOnes(env);
                for (uint16_t i = 0; i < p - 1; i++) {
                    numOnes += y[s][i];
                }
                con.add((p - 1) - numOnes == s);

            } else {
                IloExpr sumVarNext(env);
                IloExpr sumVarCur(env);

                for (uint16_t i = 0; i < p - 1; i++) {
                    sumVarNext += y[s + 1][i] * pow(2, i);
                    sumVarCur += y[s][i] * pow(2, i);
                }
                con.add(sumVarNext - sumVarCur <= -1);

                sumVarNext.end();
                sumVarCur.end();
            }
        }

        // link y variables to cumulative permutation errors

        // Create the r' variables
        BoolVar3Matrix rPrime(env, S);

        initializeBoolVar3Matrix(S, p, p, "r'", env, rPrime);

        for (uint16_t s = 0; s < S; s++) {
            for (uint16_t i = 0; i < p; i++) {
                IloExpr numOnes(env);
                for (uint16_t h = i; h < p - 1; h++) {
                    numOnes += y[s][h];
                }

                IloExpr sumRPrimeTimesJ(env);
                IloExpr sumRPrimeBinaries(env);

                for (uint16_t j = 0; j < p; j++) {
                    sumRPrimeTimesJ += j * rPrime[s][i][j];
                    sumRPrimeBinaries += rPrime[s][i][j];
                }
                con.add((p - 1) - numOnes - sumRPrimeTimesJ == 0);
                con.add(sumRPrimeBinaries == 1);

                sumRPrimeTimesJ.end();
                sumRPrimeBinaries.end();

                for (uint16_t j = 0; j < i; j++) {
                    con.add(rPrime[s][i][j] == 0);
                }
            }
        }

        // linking r' to the errors in the distribution
        const auto& distributions = scheme.getDistributions();
        for (uint16_t q = 0; q < Q; q++) {
            const auto& distribution = distributions[q];
            std::vector<uint16_t> E(p, distribution[0]);

            for (uint16_t i = 1; i < p; i++) {
                E[i] = E[i - 1] + distribution[i];
            }

            for (uint16_t s = 0; s < S; s++) {
                for (uint16_t i = 0; i < p; i++) {

                    IloExpr toRightBorder(env);

                    for (uint16_t j = i; j < p; j++) {
                        toRightBorder += rPrime[s][i][j] * E[j];
                    }

                    IloExpr toLeftBorder(env);
                    for (uint16_t j = 0; j < p - i - 1; j++) {
                        toLeftBorder += rPrime[s][i][j + i + 1] * E[j];
                    }

                    con.add(cumPerErrors[q][s][i] -
                                (toRightBorder - toLeftBorder) ==
                            0);
                }
            }
        }
    }

    /**
     * @brief Sets some extra constraints on the y variables in the case where S
     = p = k +1.
     */
    void setYConstraintsSymmetryRemoving() {
        const uint16_t& S = scheme.getSearches().size();
        const uint16_t& p = scheme.getNumParts();

        // add constraint that first search is 11111... (all extensions to
        // the right)
        if (S == p && p == k + 1) {
            for (uint16_t i = 0; i < p - 1; i++) {
                con.add(y[0][i] == 1);
            }

            for (uint16_t i = 0; i < p - 1; i++) {
                // the last search must be all zero's (all extensions to the
                // left)
                con.add(y[S - 1][i] == 0);
            }
        }
    }

    /**
     * @brief Sets the lower bound constraints.
     *
     * This function sets the lower bound constraints for the optimization
     * problem. The constraints enforce that the lower bounds on each part
     * of each search are consistent and do not violate any upper bounds.
     */
    void setLowerBoundConstraints() {
        const uint16_t& S = scheme.getSearches().size();
        const uint16_t& p = scheme.getNumParts();

        for (uint16_t s = 0; s < S; s++) {

            for (uint16_t i = 1; i < p; i++) {
                // lower bound cannot be bigger than upper bound
                con.add(L[s][i] - U[s][i] <= 0);

                // lower bound can only increase
                con.add(L[s][i] - L[s][i - 1] >= 0);
            }
        }
    }

    /**
     * @brief Sets the upper bound constraints.
     *
     * This function sets the upper bound constraints for the optimization
     * problem. The constraints enforce that the upper bounds on each part
     * of each search are consistent and satisfy specific conditions based
     * on k, p, and s. In particular this function ensures that all
     * U-strings start with 01, except in the case where k is even and p = k
     * +1, then the search at criticalSearchIdx has a U-string that starts with
     * 02.
     */
    void setUpperBoundConstraints() {
        const uint16_t& S = scheme.getSearches().size();
        const uint16_t& p = scheme.getNumParts();

        for (uint16_t s = 0; s < S; s++) {
            // add constraints that bounds cannot decrease
            for (uint16_t i = 1; i < p; i++) {
                con.add(U[s][i] - U[s][i - 1] >= 0);
            }
        }
    }

    /**
     * @brief adds extra constraints for the cases where p == k + 1 or k + 2
     * retrieved from info about the critical search.
     *
     */
    void setExtraUpperAndLowerBoundConstraints() {

        if (p == k + 1 || p == k + 2) {

            // upper bound
            for (uint16_t s = 0; s < S; s++) {
                // first part must be 0
                con.add(U[s][0] == 0);
                con.add(L[s][0] == 0);

                // second part must be 1, unless k is even, p = k+ 1 and s =
                // criticalSearchIdx
                if (k % 2 == 0 && p == k + 1 && s == criticalSearchIdx) {
                    con.add(U[s][1] == 2);
                } else {
                    con.add(U[s][1] <= 1);
                }

                // add constraint that if S= p then upperbound in the end
                // must be k! because of ..11011.. error distribution
                if (S <= p) {
                    con.add(U[s][p - 1] == k);
                }
            }

            // sum of lower bounds on last part must be sum 0 + 1 + 2 +... + k
            IloExpr sumLastL(env);
            for (uint16_t s = 0; s < S; s++) {
                // lower bound on first part must equal to zero
                con.add(L[s][0] == 0);

                sumLastL += L[s][p - 1];
            }
            if (k == p + 1 && S == p) {
                con.add(sumLastL == (k * (k + 1)) / 2);
            }
        }
    }

    /**
     * @brief Sets the constraints for covered distributions in each search.
     *
     * This function sets the constraints that force lambda to 1 if a search
     * covers a distribution and forces lambda to 0 if the search does not cover
     * the distribution. The constraints are added to the optimization problem
     * to ensure that lambda is appropriately restricted based on the cumulative
     * number of errors in distributions and the upper and lower bounds.
     */
    void setCoveredInSearchConstraints() {

        // muL[q][s][i] indicates that the cumulative number of errors in
        // distribution q up untill index i is not smaller than L[i] for search
        // s muU[q][s][i] indicates that the cumulative number of errors in
        // distribution q up untill index i is not greater than L[i] for search
        // s
        BoolVar3Matrix muL(env, Q);
        BoolVar3Matrix muU(env, Q);
        initializeBoolVar3Matrix(Q, S, p, "µL", env, muL);
        initializeBoolVar3Matrix(Q, S, p, "µU", env, muU);

        for (uint16_t q = 0; q < Q; q++) {
            for (uint16_t s = 0; s < S; s++) {

                IloExpr lowerOkay(env);
                IloExpr upperOkay(env);

                for (int i = 0; i < p; i++) {

                    // add the lower and upper indicator for this part
                    lowerOkay += muL[q][s][i];
                    upperOkay += muU[q][s][i];

                    // for lambda to zero if muL or muU is not okay
                    con.add(lambda[q][s] - muL[q][s][i] <= 0);
                    con.add(lambda[q][s] - muU[q][s][i] <= 0);

                    // add constraint that l_i = 1 => L[i] <=
                    // errorsEncountered
                    con.add(L[s][i] - cumPerErrors[q][s][i] -
                                k * (1 - muL[q][s][i]) <=
                            0);
                    // add constraint that L[i] <= errorsEncountered => l_1
                    // = 1
                    con.add(L[s][i] - cumPerErrors[q][s][i] +
                                (k + 1) * muL[q][s][i] >=
                            1);
                    // add constraint that u_i = 1 => errorsEncountered <=
                    // U[i]
                    con.add(cumPerErrors[q][s][i] - U[s][i] -
                                k * (1 - muU[q][s][i]) <=
                            0);
                    // add constraint that errorsEncountered <= U [i] => ui
                    // = 1
                    con.add(cumPerErrors[q][s][i] - U[s][i] +
                                (k + 1) * muU[q][s][i] >=
                            1);
                }

                // 2p - 1 >= lowerOkay + upperOkay - lambda forces
                // lambda to 1 if covered
                con.add(lowerOkay + upperOkay - lambda[q][s] <= 2 * p - 1);

                lowerOkay.end();
                upperOkay.end();
            }
        }
    }

    /**
     * @brief Sets the warm start for the optimization model.
     *
     * This function sets the initial variable values for the optimization
     * model using a warm start. It retrieves the upper bounds, lower
     * bounds, and part assignments from the current searches in the scheme
     * and sets the corresponding variable values in the model. This
     * function can only be executed if the model has been extracted.
     *
     * @param lower boolean value that indicates whether the lower bounds of
     * the warm start should be considered
     */
    void setWarmStart(bool lower = true) {
        // warm start
        IloNumVarArray startVar(env);
        IloNumArray startVal(env);
        const auto& searches = scheme.getSearches();
        const auto& S = searches.size();
        const auto& p = scheme.getNumParts();
        for (uint16_t s = 0; s < S; s++) {
            // set U, L and pi values

            const auto firstPart = searches[s].getPart(0);
            for (uint16_t i = 0; i < p; i++) {

                startVar.add(U[s][i]);
                startVal.add(searches[s].getUpperBound(i));

                if (lower) {
                    startVar.add(L[s][i]);
                    startVal.add(searches[s].getLowerBound(i));
                }

                // create the binary encoding of the search (y variables)
                if (i > 0) {
                    startVar.add(y[s][i - 1]);
                    startVal.add((searches[s].getPart(i) > firstPart) ? 1 : 0);
                }
            }
        }

        cplex.addMIPStart(startVar, startVal);
        startVal.end();
        startVar.end();
    }

    /**
     * @brief Greedily solves the scheme adaptation problem.
     *
     * This function applies a greedy algorithm to adapt the scheme by
     * iteratively updating the searches to improve the coverage. After
     * adaptation, it prints the adapted scheme, coverage statistics, and
     * the maximum coverage achieved.
     */
    void greedySolve() {
        // greedily adapt scheme
        scheme.adapt();

        std::cout << "Greedy solution: " << std::endl;
        std::cout << scheme << std::endl;
        std::cout << "Coverage: " << std::endl;
        scheme.printCoverages();
        std::cout << "max cov: " << scheme.getMaxCoverage() << std::endl;
    }

    /**
     * @brief Creates the optimal searches found by the model
     *
     * Using the encodings in the model which was solved to optimality the
     * optimal searches are created.
     *
     * @returns the optimal searches
     */
    std::vector<AdaptableSearch> getOptSearches() const {
        const auto& S = scheme.getSearches().size();
        const auto& p = scheme.getNumParts();
        // Create the optimal searches
        std::vector<AdaptableSearch> optSearches;
        optSearches.reserve(S);

        for (uint16_t s = 0; s < S; s++) {
            std::vector<uint16_t> optU(p);
            std::vector<uint16_t> optL(p);
            std::vector<uint16_t> optPi(p, -1);

            // create permutation from binary encoding
            uint16_t start = 0;
            for (uint16_t i = 0; i < p - 1; i++) {

                start += (fabs(cplex.getValue(y[s][i])) <
                          std::numeric_limits<float>::epsilon())
                             ? 1
                             : 0;
                optPi[0] = start;
            }

            uint16_t leftMost = start;
            uint16_t rightMost = start;

            for (int i = 0; i < p; i++) {

                optU[i] = round(cplex.getValue(U[s][i]));
                optL[i] = round(cplex.getValue(L[s][i]));

                if (i > 0) {

                    if (fabs(cplex.getValue(y[s][i - 1])) <
                        std::numeric_limits<float>::epsilon()) {
                        // left
                        optPi[i] = --leftMost;
                    } else {
                        optPi[i] = ++rightMost;
                    }
                }
            }

            optSearches.push_back(Search::makeSearch(optPi, optL, optU, s));
        }
        return optSearches;
    }
    /**
     * @brief Solves the scheme adaptation problem using the CPLEX
     * optimization model.
     *
     * This function constructs and solves the CPLEX optimization model for
     * the scheme adaptation problem. It adds the constraints and objectives
     * to the model, sets a warm start with the current solution if the
     * model has not been initialized yet, solves the model, and retrieves
     * the optimal solution if found. The adapted searches are updated in
     * the scheme based on the optimal solution.
     *
     * @param init Indicates whether the model has already been extracted
     * and solved using extra constraints. Default value is false.
     */
    void solveCplexModel(bool init = false) {

        if (!init) {

            // add all constraints to the model
            model.add(con);

            // add the objectives to the model
            model.add(IloMinimize(env, createObjectiveMinU() +
                                           createObjectiveMaxL() +
                                           createObjectiveMinCov() -
                                           scheme.getDistributions().size()));
            // extract the model
            cplex.extract(model);
        } else {
            // keep the sum of the values for the upper bounds

            uint64_t sumU = 0;
            const uint16_t& S = scheme.getSearches().size();
            const uint16_t& p = scheme.getNumParts();
            for (uint16_t s = 0; s < S; s++) {
                for (uint16_t i = 0; i < p; i++) {
                    auto weight = getWeightsU(i);
                    sumU += weight * scheme.getSearches()[s].getUpperBound(i);
                }
            }
            con.add(createObjectiveMinU() == sumU);

            cplex.deleteMIPStarts(0, cplex.getNMIPStarts());
        }

        // give a warm start with the current solution
        setWarmStart();

        setParameters();

        // solve the model
        cplex.solve();

        // Check if the solution is optimal
        if (cplex.getStatus() == IloAlgorithm::Optimal) {
            // Get the objective value
            double objectiveValue = cplex.getObjValue();
            std::cout << "Objective value: " << objectiveValue << std::endl;

            auto optSearches = getOptSearches();

            // set the scheme to the newly found solution
            scheme.setSearches(optSearches);
        } else {
            std::cout << "No optimal solution found." << std::endl;
            std::cout << "Objective value: " << cplex.getObjValue()
                      << std::endl;
        }

        std::cout << cplex.getStatus() << std::endl;
        cplex.end();
    }

    /**
     * @brief Solves the scheme adaptation problem using the CPLEX
     * optimization model without considering lower bounds.
     *
     * This function constructs and solves the CPLEX optimization model for
     * the scheme adaptation problem without considering the lower bounds.
     * It adds the constraints and objectives to the model, sets a warm
     * start with the current solution, solves the model, and retrieves the
     * optimal solution if found. The adapted searches are updated in the
     * scheme based on the optimal solution.
     */
    void solveCplexModelWithoutL() {

        const auto& S = scheme.getSearches().size();
        const auto& p = scheme.getNumParts();

        // add all constraints to the model
        model.add(con);

        // add the objectives to the model
        model.add(IloMinimize(env, createObjectiveMinU() +
                                       createObjectiveMaxL() +
                                       createObjectiveMinCov() -
                                       scheme.getDistributions().size()));

        // add constraint that all lower bounds are zero
        IloConstraintArray lbZero(env);
        for (uint16_t s = 0; s < S; s++) {
            for (uint16_t i = 0; i < p; i++) {
                // lower bound must be 0
                lbZero.add(L[s][i] == 0);
            }
        }
        model.add(lbZero);

        // extract the model
        cplex.extract(model);

        // give a warm start with current solution, without considering
        // lower bounds
        setWarmStart(false);

        setParameters();

        // solve the model
        cplex.solve();

        // Check if the solution is optimal
        if (cplex.getStatus() == IloAlgorithm::Optimal) {
            // Get the objective value
            double objectiveValue = cplex.getObjValue();
            std::cout << "Objective value: " << objectiveValue << std::endl;

            auto optSearches = getOptSearches();

            for (const auto& s : optSearches) {
                std::cout << s << std::endl;
            }

            // set the scheme to the newly found solution
            scheme.setSearches(optSearches);
        } else {
            std::cout << "No optimal solution found." << std::endl;
            std::cout << "Objective value: " << cplex.getObjValue()
                      << std::endl;
        }

        std::cout << cplex.getStatus() << std::endl;

        // remove the constraints that set lower bound to zero
        model.remove(lbZero);
    }

    /**
     * @brief Sets the parameters for the CPLEX solver.
     *
     * This function sets the parameters for the CPLEX solver, such as
     * display level, emphasis on BEST bound, and the backtrack strategy.
     */
    void setParameters() {
        cplex.setParam(IloCplex::Param::MIP::Display, 4);
        cplex.setParam(IloCplex::Param::Emphasis::MIP,
                       CPX_MIPEMPHASIS_BESTBOUND);
        cplex.setParam(IloCplex::Param::MIP::Strategy::Backtrack, 0.1);
        cplex.setParam(IloCplex::Param::Read::Scale, 1);
        cplex.setParam(IloCplex::Param::MIP::Strategy::HeuristicFreq, -1);
        cplex.setParam(IloCplex::Param::MIP::Strategy::VariableSelect,
                       CPX_VARSEL_PSEUDOREDUCED);
    }

  public:
    /**
     * @brief Constructs a CPLEXModel object.
     *
     * This constructor initializes a CPLEXModel object with the given
     * environment, scheme, and parameters. It sets up the optimization
     * model, constraints, and variables based on the provided scheme.
     *
     * @param env The CPLEX environment.
     * @param k The parameter value.
     * @param scheme The adaptable scheme to be optimized.
     * @param criticalSearchIdx The index of the critical search.
     */
    CPLEXModel(IloEnv& env, uint16_t k, AdaptableScheme& scheme,
               uint16_t criticalSearchIdx)
        : env(env), model(env), cplex(model), con(env), k(k),
          S(scheme.getSearches().size()), Q(scheme.getDistributions().size()),
          p(scheme.getNumParts()), lambda(env, Q), U(env, S), L(env, S),
          cumPerErrors(env, Q), y(env, S), scheme(scheme),
          criticalSearchIdx(criticalSearchIdx) {
        // initialize the binary and integer variables
        initializeLambda();
        initializeUpper();
        initializeLower();
        initializeCumPerError();
        initializeY();

        // Set the constraints
        setMinimalCoverconstraints();

        setYConstraints();
        setYConstraintsSymmetryRemoving();

        setLowerBoundConstraints();
        setUpperBoundConstraints();

        setExtraUpperAndLowerBoundConstraints();

        setCoveredInSearchConstraints();
    }

    /**
     * @brief Solves the model.
     *
     * This function solves the scheme adaptation problem by performing a
     * series of steps. It first sorts the scheme on the order of the first
     * part. If the parameter k is even and the number of parts in the
     * scheme is k+1, it sets U[1] to 2 for the critical search and
     * recalculates the coverage. Otherwise, it recalculates the initial
     * coverage after sorting. Then, it applies a greedy algorithm to adapt
     * the scheme and solves the CPLEX optimization model using the greedy
     * solution as a warm start.
     */
    void solveModel() {
        // sort scheme on order
        scheme.sortOnFirstPart();
        if (k % 2 == 0 && scheme.getNumParts() == k + 1) {

            // set U[1] == 2 for the critical search
            scheme.setUAtIndextoTwoAllOthers1(criticalSearchIdx);
        }
        // recalculate coverage
        scheme.recalculateInitialCoverage();

        greedySolve();

        bool step1 = false;

        std::cout << "Step 1: Solve without considering lower bound."
                  << std::endl;
        step1 = true;
        solveCplexModelWithoutL();
        greedySolve();

        std::cout << "Step 2: Solve entire problem." << std::endl;
        solveCplexModel(step1);
    }

    /**
     * @brief Returns the adaptable scheme.
     *
     * This function returns a reference to the adaptable scheme being
     * optimized.
     *
     * @return A reference to the adaptable scheme.
     */
    AdaptableScheme& getScheme() {
        return scheme;
    }
};

#endif
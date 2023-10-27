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
#ifndef GREEDY_GENERATOR_H
#define GREEDY_GENERATOR_H

#include <algorithm>
#include <cassert>
#include <cinttypes>
#include <iostream>
#include <numeric>
#include <vector>

// ============================================================================
// CLASS ERROR PATTERN
// ============================================================================
/**
 * @brief ErrorDistribution class represents an error distribution.
 *
 * The ErrorDistribution class is derived from std::vector<uint16_t> and
 * represents an error distribution. It stores the distribution as a vector of
 * uint16_t values and includes an index to identify the distribution.
 */
class ErrorDistribution : public std::vector<uint16_t> {
  private:
    uint64_t index; // the unique index of the distribution

  public:
    /**
     * @brief Constructor for creating an ErrorDistribution object.
     *
     * This constructor initializes an ErrorDistribution object with the
     * specified distribution and index.
     *
     * @param distribution The error distribution represented as a vector of
     * uint16_t values.
     * @param index The index to identify the error distribution.
     */
    ErrorDistribution(const std::vector<uint16_t>& distribution, uint64_t index)
        : std::vector<uint16_t>(distribution), index(index) {
    }

    /**
     * @brief Get the index of the error distribution.
     *
     * This function returns the index of the error distribution.
     *
     * @return The index of the error distribution.
     */
    uint64_t getIndex() const {
        return index;
    }

    /**
     * @brief Overloaded output stream operator for printing the error
     * distribution.
     *
     * This friend function overloads the output stream operator to print the
     * error distribution in the format: "<index>: <element1>, <element2>, ...,
     * <elementN>".
     *
     * @param os The output stream to write the error distribution.
     * @param distribution The ErrorDistribution object to be printed.
     * @return The output stream after writing the error distribution.
     */
    friend std::ostream& operator<<(std::ostream& os,
                                    const ErrorDistribution& distribution) {
        os << distribution.getIndex() << ": ";

        for (uint16_t i = 0; i < distribution.size() - 1; i++) {
            os << distribution[i] << ",";
        }

        os << distribution.back();
        return os;
    }
};

/**
 * @brief Calculate the number of combinations (n choose k).
 *
 * This function calculates the number of possibilities to choose k distinct
 * elements out of a set of n distinct elements using the combination formula.
 *
 * @param n The number of elements to choose from.
 * @param k The number of elements to choose.
 * @returns The number of combinations (n choose k).
 */
uint64_t nChoosek(uint16_t n, uint16_t k);

/**
 * @brief Generate the lexicographic m-subset successor of A.
 *
 * This function generates the lexicographic m-subset successor of a given
 * subset A of a set with n elements. The first m elements of A will be replaced
 * by its successor if it exists. If no successor exists, the function returns
 * false.
 *
 * @param A The subset of {0,..,n-1} with m elements. At the end of the
 * function, the first m elements of A will be overwritten by its successor if
 * it exists.
 * @param m The number of elements in the subsets.
 * @param n The number of symbols in the superset.
 * @returns True if the successor exists and A is updated, false otherwise.
 */
bool successor(std::vector<uint16_t>& A, uint16_t m, uint16_t n);

/**
 * @brief Generate all error distributions with at most k errors divided over p
 * parts.
 *
 * This function generates all error distributions where at most k errors are
 * divided over p parts. It generates all possible combinations of indices for
 * the given number of parts and errors.
 *
 * @param k The maximum number of errors.
 * @param p The number of parts.
 * @returns A vector containing all error distributions with p parts and at most
 * k errors.
 */
std::vector<ErrorDistribution> generateErrorDistributions(const uint16_t k,
                                                          const uint16_t p);

// ============================================================================
// CLASS SEARCH
// ============================================================================

/**
 * @brief Search class represents a search configuration.
 *
 * The Search class represents a search configuration with order, lower bounds,
 * upper bounds, and an index. It provides functions to create search
 * configurations, as well as perform sanity checks and comparisons.
 */
class Search {
  protected:
    std::vector<uint16_t> L;     // the vector with lower bounds
    std::vector<uint16_t> U;     // the vector with upper bounds
    std::vector<uint16_t> order; // the vector with the order of the parts
    uint16_t sIdx;               // the index of the search
  private:
    /**
     * @brief Private constructor for creating a Search object.
     *
     * This constructor initializes a Search object with the specified order,
     * lower bounds, upper bounds, and index.
     *
     * @param order The order of the search.
     * @param lowerBounds The lower bounds of the search.
     * @param upperBounds The upper bounds of the search.
     * @param idx The index of the search.
     */
    Search(std::vector<uint16_t>& order, std::vector<uint16_t>& lowerBounds,
           std::vector<uint16_t>& upperBounds, uint16_t idx)
        : L(lowerBounds), U(upperBounds), order(order), sIdx(idx) {
    }

    /**
     * @brief Checks if the connectivity property is satisfied.
     *
     * This function checks if the connectivity property is satisfied by
     * examining the order of the parts. The connectivity property requires that
     * a part p can only be processed at index i (i > 1) if part p-1 or part p+1
     * has been processed at one of the indexes lower than i.
     *
     * @returns True if the connectivity property is satisfied, false otherwise.
     */
    bool connectivitySatisfied() const {
        // Initialize the highest and lowest seen parts as the first part in the
        // order
        uint16_t highestSeen = order[0];
        uint16_t lowestSeen = order[0];

        // Iterate through the remaining parts in the order
        for (uint16_t i = 1; i < order.size(); i++) {
            // Check if the current part is the next consecutive part after the
            // highest seen part
            if (order[i] == highestSeen + 1) {
                highestSeen++; // Update the highest seen part
            }
            // Check if the current part is the previous consecutive part before
            // the lowest seen part
            else if (order[i] == lowestSeen - 1) {
                lowestSeen--; // Update the lowest seen part
            }
            // If the current part is not adjacent to the highest or lowest seen
            // parts, the connectivity property is violated
            else {
                return false;
            }
        }
        return true; // All parts satisfy the connectivity property
    }

    /**
     * @brief Checks if the search is zero-based.
     *
     * This function checks if the search is zero-based, meaning that the order
     * must contain zero.
     *
     * @returns True if the search is zero-based, false otherwise.
     */
    bool zeroBased() const {
        return *std::min_element(order.begin(), order.end()) == 0;
    }

    /**
     * @brief Performs a sanity check on the search configuration.
     *
     * This function performs a sanity check on the search configuration to
     * ensure that it satisfies the necessary conditions. It checks if the
     * search is zero-based, if the connectivity property is satisfied, and if
     * the bounds are valid.
     *
     * @throws std::runtime_error If any of the conditions are violated.
     */
    void sanityCheck() const {
        if (!zeroBased()) {
            throw std::runtime_error("Search is not zero-based!");
        }
        if (!connectivitySatisfied()) {
            for (const auto& i : order) {
                std::cout << i << " ";
            }
            std::cout << std::endl;
            throw std::runtime_error("Connectivity property violated!");
        }
        if (!validBounds()) {
            throw std::runtime_error("Bounds are not valid!");
        }
    }

  public:
    /**
     * @brief Static function to construct a search.
     *
     * This static function creates a Search object based on the specified
     * order, lower bounds, upper bounds, and index. It performs a sanity check
     * on the search configuration to ensure its validity.
     *
     * @param order The order of the search.
     * @param lowerBounds The lower bounds of the search.
     * @param upperBounds The upper bounds of the search.
     * @param sIdx The index of the search.
     * @returns The constructed Search object.
     *
     * @throws std::runtime_error If the sizes of the vectors are not equal or
     * if the sanity check fails.
     */
    static Search makeSearch(std::vector<uint16_t> order,
                             std::vector<uint16_t> lowerBounds,
                             std::vector<uint16_t> upperBounds, uint16_t sIdx) {
        // check correctness of sizes
        if (order.size() != lowerBounds.size() ||
            order.size() != upperBounds.size()) {
            throw std::runtime_error("Could not create search, the sizes of "
                                     "all vectors are not equal");
        }

        const Search s = Search(order, lowerBounds, upperBounds, sIdx);

        // do a sanity check
        s.sanityCheck();

        // sanity check passed
        return s;
    }

    /**
     * @brief Get the lower bound for the specified part.
     *
     * This function returns the lower bound for the part at the specified
     * index.
     *
     * @param idx The index of the part.
     * @returns The lower bound for the specified part.
     */
    uint16_t getLowerBound(uint16_t idx) const {
        assert(idx < L.size());
        return L[idx];
    }

    /**
     * @brief Get the upper bound for the specified part.
     *
     * This function returns the upper bound for the part at the specified
     * index.
     *
     * @param idx The index of the part.
     * @returns The upper bound for the specified part.
     */
    uint16_t getUpperBound(uint16_t idx) const {
        assert(idx < U.size());
        return U[idx];
    }
    /**
     * @brief Get the part at the specified index.
     *
     * This function returns the part at the specified index in the search
     * configuration.
     *
     * @param idx The index of the part.
     * @returns The part at the specified index.
     */
    uint16_t getPart(uint16_t idx) const {
        assert(idx < order.size());
        return order[idx];
    }

    /**
     * @brief Get the number of parts in this search configuration.
     *
     * This function returns the number of parts in the search configuration.
     *
     * @returns The number of parts.
     */
    uint16_t getNumParts() const {
        return order.size();
    }

    /**
     * @brief Get the index of the search configuration.
     *
     * This function returns the index of the search configuration.
     *
     * @returns The index of the search configuration.
     */
    uint16_t getIndex() const {
        return sIdx;
    }

    /**
     * @brief Checks if the upper and lower bounds are valid.
     *
     * This function checks if the upper and lower bounds are valid, i.e., if
     * they are non-decreasing and L[i] <= U[i] for each index i.
     *
     * @returns True if the bounds are valid, false otherwise.
     */
    bool validBounds() const {
        // check if bounds at index 0 make sense
        if (L[0] > U[0]) {
            return false;
        }
        for (uint16_t i = 1; i < order.size(); i++) {
            if (L[i] > U[i] || L[i] < L[i - 1] || U[i] < U[i - 1]) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Compare two search configurations.
     *
     * This function compares two search configurations and determines their
     * ordering based on their hardness. It first compares the upper bounds,
     * then the lower bounds, and finally the indices. A harder search is
     * considered less than an easier search.
     *
     * @param rhs The Search object to compare with.
     * @returns True if the current search configuration is less than the rhs
     * configuration, false otherwise.
     *
     * @throws std::out_of_range If the number of parts in the configurations is
     * not equal.
     */
    bool operator<(const Search& rhs) const {
        if (rhs.getNumParts() != getNumParts()) {
            throw std::out_of_range(
                "The number of parts in the compared searches are not equal!");
        }

        // A) compare U-string
        for (uint16_t i = 0; i < getNumParts(); i++) {
            if (U[i] != rhs.U[i]) {
                return U[i] > rhs.U[i];
            }
        }

        // B) compare L-string
        for (uint16_t i = 0; i < getNumParts(); i++) {
            if (L[i] != rhs.L[i]) {
                return L[i] < rhs.L[i];
            }
        }

        // C) sort on reverse index
        return sIdx > rhs.sIdx;
    }
};

// ============================================================================
// CLASS ADAPTABLE SEARCH
// ============================================================================

/**
 * @brief Represents an adaptable search.
 *
 * The AdaptableSearch class is a derived class of the Search class and provides
 * additional functionality to adapt it. It includes methods for modifying U
 * (upper bound), L (lower bound), order, and the index (sIdx) of the search.
 */
class AdaptableSearch : public Search {
  public:
    /**
     * Decreases the value of U[i] by 1.
     * @param i The index of U to be decreased.
     */
    void decreaseU(uint16_t i) {
        U[i]--;
    }

    /**
     * Increases the value of U[i] by 1.
     * @param i The index of U to be increased.
     */
    void increaseU(uint16_t i) {
        U[i]++;
    }

    /**
     * Decreases the value of L[i] by 1.
     * @param i The index of L to be decreased.
     */
    void decreaseL(uint16_t i) {
        L[i]--;
    }

    /**
     * Increases the value of L[i] by 1.
     * @param i The index of L to be increased.
     */
    void increaseL(uint16_t i) {
        L[i]++;
    }

    /**
     * Sets the value of L[index].
     * @param index The index of L to be set.
     * @param value The value to set for L[index].
     */
    void setL(uint16_t index, uint16_t value) {
        L[index] = value;
    }

    /**
     * Sets the value of U[index].
     * @param index The index of U to be set.
     * @param value The value to set for U[index].
     */
    void setU(uint16_t index, uint16_t value) {
        U[index] = value;
    }

    /**
     * Sets the value of order[index].
     * @param index The index of order to be set.
     * @param value The value to set for order[index].
     */
    void setPi(uint16_t index, uint16_t value) {
        order[index] = value;
    }

    /**
     * Constructs an AdaptableSearch object from a given Search object.
     * @param s The Search object to be used for construction.
     */
    AdaptableSearch(const Search& s) : Search(s) {
    }

    /**
     * Sets the value of sIdx.
     * @param i The value to set for sIdx.
     */
    void setIdx(uint16_t i) {
        sIdx = i;
    }
};

/**
 * Operator overloading. Outputs the search to the output stream
 * @param os, the output stream
 * @param obj, the search to print
 */
std::ostream& operator<<(std::ostream& os, const Search& obj);

// ============================================================================
// CLASS ADAPTABLE SCHEME
// ============================================================================

/**
 * @brief Represents an adaptable scheme.
 *
 * The AdaptableScheme class represents a scheme that consists of adaptable
 * searches. It provides methods for adapting the scheme to minimize the upper
 * bounds (U-string) and maximize the lower bounds (L-string) of the searches.
 * The scheme maintains a set of error distributions and tracks the coverage of
 * each distribution by the searches. The scheme ensures that each distribution
 * is covered by at least one search and keeps track of critical distributions
 * that are covered by only one search.
 */
class AdaptableScheme {
  private:
    const std::vector<ErrorDistribution>
        distributions; // the error distributions for k errors and p parts
    std::vector<AdaptableSearch> searches; // the searches

    std::vector<std::vector<bool>>
        coverageMatrix; // matrix that for each distribution index and each
                        // search index tell whether this search covers this
                        // distribution

    std::vector<std::vector<uint32_t>>
        criticalDistributions; // For each search the list of critical
                               // distribution indexes is kept. A distribution
                               // is considered critical if it is not covered by
                               // any other search.

    /**
     * @brief Initializes the coverageMatrix and criticalDistributions.
     *
     * This function sets up the initial coverage information and identifies
     * critical distributions. It populates the coverageMatrix matrix based on
     * the coverage of each distribution by the searches. It also determines the
     * critical distributions that are only covered by a single search.
     */
    void setInitialCoverage() {

        uint16_t nSearches = searches.size(); // the number of searches

        // Resize the vectors
        criticalDistributions.clear();
        criticalDistributions.resize(nSearches);
        coverageMatrix.resize(distributions.size());

        for (uint32_t pIdx = 0; pIdx < distributions.size(); pIdx++) {
            coverageMatrix[pIdx].resize(nSearches);

            // Set coverage for each search
            for (const auto& s : searches) {
                coverageMatrix[pIdx][s.getIndex()] = covers(s, pIdx);
            }

            uint16_t coveringSearches = numberOfCoveringSearches(pIdx);

            if (coveringSearches == 0) {
                const auto& problem = distributions[pIdx];
                std::cout << problem << std::endl;
                /* std::stringstream ss();
                ss << problem; */
                throw std ::runtime_error(
                    "Not all distributions are covered\nDistribution: " +
                    std::to_string(pIdx) + " is not covered");
            }
            if (coveringSearches == 1) {
                // add this error distribution to the list of critical
                // distributions for its covering search
                addToCriticalDistributions(pIdx);
            }
        }
    }

    /**
     * @brief Adds a distribution to the list of critical distributions for a
     * search.
     *
     * This function is called when a distribution becomes critical for a
     * search. It adds the distribution with the specified index to the list of
     * critical distributions for the only search that covers it. This function
     * can only be called when the number of covering searches is 1 (and hence
     * the distribution is critical for its covering search).
     *
     * @param pIdx The index of the distribution to be added to the critical
     * distributions.
     */
    void addToCriticalDistributions(uint32_t pIdx) {
        criticalDistributions[idxOfCoveringSearch(pIdx)].emplace_back(pIdx);
    }

    /**
     * @brief Calculates the number of covering searches for a distribution.
     *
     * This function determines the number of searches that cover the
     * distribution with the specified index.
     *
     * @param pIdx The index of the distribution for which the number of
     * covering searches is calculated.
     * @returns The number of searches that cover the distribution.
     */
    uint16_t numberOfCoveringSearches(uint32_t pIdx) const {
        return std::accumulate(coverageMatrix[pIdx].begin(),
                               coverageMatrix[pIdx].end(), 0);
    }

    /**
     * @brief Finds the index of the solely covering search for a distribution.
     *
     * This function returns the index of the search that solely covers the
     * distribution with the specified index. It can only be called when the
     * number of covering searches for the distribution equals one.
     *
     * @param pIdx The index of the distribution for which the index of the
     * covering search is returned.
     * @returns The index of the solely covering search for the distribution.
     */
    uint16_t idxOfCoveringSearch(uint32_t pIdx) const {
        assert(numberOfCoveringSearches(pIdx) == 1);

        auto it = std::find(coverageMatrix[pIdx].begin(),
                            coverageMatrix[pIdx].end(), true);
        return std::distance(coverageMatrix[pIdx].begin(), it);
    }

    /**
     * @brief Checks if a search covers a distribution.
     *
     * This function checks if the distribution with the specified index is
     * covered by the given search. It compares the cumulative number of errors
     * in the search's parts with the bounds of the search to determine if the
     * distribution is covered.
     *
     * @param s The search to check for coverage.
     * @param pIdx The index of the distribution to be checked.
     * @return True if the distribution is covered by the search, False
     * otherwise.
     */
    bool covers(const Search& s, uint32_t pIdx) const {
        uint16_t totalErrors = 0;

        // Iterate over the parts of the search
        for (uint16_t i = 0; i < s.getNumParts(); i++) {
            uint16_t part = s.getPart(i);
            // add the number of errors
            totalErrors += distributions[pIdx][part];

            // check if the cumulative number of errors up untill this part are
            // between the bounds of the search
            if (totalErrors > s.getUpperBound(i) ||
                totalErrors < s.getLowerBound(i)) {
                // one of the bounds is violated
                return false;
            }
        }

        return true;
    }

    /**
     * @brief Sorts a vector containing the indexes of the searches based on the
     * hardness of the search.
     *
     * This function sorts the vector of search indexes (`order`) based on the
     * hardness of each search. The index of the hardest search will come first
     * in the vector, followed by the indexes of the remaining searches in
     * descending order of their hardness.
     *
     * @param order A vector containing the indexes of the searches. This vector
     * will be modified during the sorting process.
     */
    void sortSIdxOrder(std::vector<uint16_t>& order) const {
        std::sort(order.begin(), order.end(),
                  [this](const uint16_t lhs, const uint16_t rhs) {
                      return searches[lhs] < searches[rhs];
                  });
    }

    /**
     * @brief Checks if the adapted search at index `sIdx` still covers all its
     * critical distributions.
     *
     * This helper function checks if the search at the specified index (`sIdx`)
     * still covers all the critical distributions assigned to it. It iterates
     * over the critical distributions list of `sIdx` and checks if each
     * distribution is still covered by the search. If any distribution is not
     * covered, it returns false; otherwise, it returns true.
     *
     * @param sIdx The index of the search to check.
     * @return A boolean value indicating whether the search still covers all
     * its critical distributions.
     */
    bool coversCritical(uint16_t sIdx) {
        const auto& s = searches[sIdx];
        for (const auto& pIdx : criticalDistributions[sIdx]) {
            if (!covers(s, pIdx)) {
                return false;
            }
        }
        return true;
    }

    /**
     * @brief Updates the critical distributions after the search at index
     * `sIdx` has been adapted.
     *
     * This helper function updates the critical distributions list and coverage
     * matrix after the search at index `sIdx` has been adapted. It iterates
     * over the distributions that were previously covered by `sIdx` and checks
     * if any of them have become critical for another search. If a distribution
     * is no longer covered by `sIdx`, its coverage in the coverage matrix is
     * set to false, and if the distribution has become critical, it is added to
     * the critical distributions list.
     *
     * @param sIdx The index of the search that was adapted.
     */
    void updateCritical(uint16_t sIdx) {
        for (uint32_t pIdx = 0; pIdx < distributions.size(); pIdx++) {
            // Check if this search used to cover the distribution
            if (coverageMatrix[pIdx][sIdx]) {
                // Check if search does not cover distribution anymore
                if (!covers(searches[sIdx], pIdx)) {
                    // set coverage to false
                    coverageMatrix[pIdx][sIdx] = false;
                    // Check if distribution has become critical
                    if (numberOfCoveringSearches(pIdx) == 1) {
                        addToCriticalDistributions(pIdx);
                    }
                }
            }
        }
    }

  public:
    /**
     * @brief Constructor for creating an AdaptableScheme object.
     *
     * This constructor initializes an AdaptableScheme object with the specified
     * searches and error distributions.
     *
     * @param searches The searches that will be adapted.
     * @param distributions The error distributions for this scheme.
     */
    AdaptableScheme(std::vector<Search>& searches,
                    const std::vector<ErrorDistribution> distributions)
        : distributions(distributions) {

        for (auto& s : searches) {
            this->searches.emplace_back(s);
        }

        setInitialCoverage();
    }

    /**
     * @brief Constructor for creating an AdaptableScheme object.
     *
     * This constructor initializes an AdaptableScheme object with the specified
     * adaptable searches and error distributions.
     *
     * @param searches The adaptable searches that will be adapted.
     * @param distributions The error distributions for this scheme.
     */
    AdaptableScheme(std::vector<AdaptableSearch>& searches,
                    const std::vector<ErrorDistribution>& distributions)
        : distributions(distributions), searches(searches) {

        setInitialCoverage();
    }

    /**
     * @brief Adapts the scheme by minimizing the U-string and maximizing the
     * L-string.
     *
     * This function greedily adapts the scheme to minimize the U-string and
     * maximize the L-string. It first minimizes the U-string by iterating over
     * the parts. At each index i, the searches are sorted based on their
     * hardness, and starting with the hardest search, the U[i] value is
     * decreased until either the bounds are violated or not all critical
     * distributions are covered anymore. The final violating decrease is
     * undone, and the algorithm switches to the next hardest search after
     * updating the critical distributions. Once all indexes i are finished, the
     * algorithm switches to maximizing the lower bounds. This is done
     * analogously, but now the parts are processed from right to left, and L[i]
     * is increased until the bounds are violated or not all critical
     * distributions are covered.
     */
    void adapt() {

        std::vector<uint16_t> sIdxOrder(searches.size());
        std::iota(sIdxOrder.begin(), sIdxOrder.end(), 0);

        uint16_t numParts = getNumParts();

        // A) Minimize the upper bounds
        for (uint16_t i = 0; i < numParts; i++) {

            sortSIdxOrder(sIdxOrder);
            for (const auto& sIdx : sIdxOrder) {

                auto& s = searches[sIdx];

                while (true) {
                    // A.1) Check if 0 is reached
                    if (s.getUpperBound(i) == 0) {
                        break;
                    }
                    // A.2) Decrease value at position i
                    s.decreaseU(i);
                    // A.3) check validity of search scheme
                    if (!s.validBounds() || !coversCritical(sIdx)) {
                        // Undo final decrease
                        s.increaseU(i);
                        break;
                    }
                }
                // Update the critical error distributions for the other
                // searches
                updateCritical(sIdx);
            }
        }
        // B) Maximize the lower bounds
        for (uint16_t i = numParts; i-- > 0;) {

            sortSIdxOrder(sIdxOrder);
            for (const auto& sIdx : sIdxOrder) {

                auto& s = searches[sIdx];

                while (true) {
                    // B.1) INCREASE value at position i
                    s.increaseL(i);
                    // B.2) check validity of search scheme
                    if (!s.validBounds() || !coversCritical(sIdx)) {
                        // Undo final increase
                        s.decreaseL(i);
                        break;
                    }
                }

                // Update the critical error distributions
                updateCritical(sIdx);
            }
        }
    }

    /**
     * @brief Output the search scheme to an output stream.
     *
     * This friend function outputs the search scheme to the specified output
     * stream. It prints the searches one by one.
     *
     * @param os The output stream to output the searches to.
     * @param scheme The AdaptableScheme object to be output.
     * @return The output stream after writing the search scheme.
     */
    friend std::ostream& operator<<(std::ostream& os,
                                    const AdaptableScheme& scheme) {

        os << "Searches: "
           << "\n";

        for (const auto& s : scheme.searches) {

            os << s << "\n";
        }

        return os;
    }

    /**
     * @brief Calculates the number of distributions per coverage.
     *
     * This function calculates the number of distributions for each coverage
     * level and returns a vector of size `searches.size()`. Each element of the
     * vector represents the number of distributions that have a coverage value
     * equal to the index + 1. The coverage is determined based on the number of
     * searches covering each distribution.
     *
     * @return A vector where each element represents the number of
     * distributions for each coverage level.
     */
    std::vector<uint64_t> getCoverages() const {
        std::vector<uint64_t> errorDistributionsPerNSearches(searches.size(),
                                                             0);
        for (uint64_t pidx = 0; pidx < distributions.size(); pidx++) {
            errorDistributionsPerNSearches[numberOfCoveringSearches(pidx) -
                                           1]++;
        }
        return errorDistributionsPerNSearches;
    }

    /**
     * @brief Get the vector of searches in the scheme.
     *
     * This function returns a constant reference to the vector of searches in
     * the scheme.
     *
     * @return A constant reference to the vector of searches.
     */
    const std::vector<AdaptableSearch>& getSearches() const {
        return searches;
    }

    /**
     * @brief Print the coverages of the search scheme.
     *
     * This function prints the number of error distributions covered by each
     * number of searches. It calculates the coverages based on the number of
     * searches that cover each distribution.
     */
    void printCoverages() const {

        std::vector<uint64_t> errorDistributionsPerNSearches = getCoverages();

        for (size_t i = 0; i < searches.size(); i++) {
            std::cout << errorDistributionsPerNSearches[i] << "x " << i + 1
                      << " covering searches" << std::endl;
        }
    }

    /**
     * @brief Get the maximal number of covering searches for this scheme.
     *
     * This function returns the number of searches that cover the distribution
     * with the most covering searches.
     *
     * @return The number of searches that cover the distribution with the most
     * covering searches.
     */
    int getMaxCoverage() const {
        std::vector<uint64_t> errorDistributionsPerNSearches = getCoverages();
        for (size_t i = searches.size(); i-- > 0;) {
            if (errorDistributionsPerNSearches[i] > 0) {
                return i + 1;
            }
        }

        return 0;
    }

    /**
     * @brief Get the vector of error distributions in the scheme.
     *
     * This function returns a constant reference to the vector of error
     * distributions in the scheme.
     *
     * @return A constant reference to the vector of error distributions.
     */
    const std::vector<ErrorDistribution>& getDistributions() const {
        return distributions;
    }

    /**
     * @brief Get the number of parts in the scheme.
     *
     * This function returns the number of parts in the scheme. It retrieves the
     * number of parts from the first search in the vector of searches.
     *
     * @return The number of parts in the scheme.
     */
    uint16_t getNumParts() const {
        return searches[0].getNumParts();
    }

    /**
     * @brief Helper function to sort the searches based on the first processed
     * part.
     *
     * This function sorts the searches in the scheme so that the search that
     * processes part i as its first part comes before all searches for which
     * their first processed part is greater than part i.
     */
    void sortOnFirstPart() {
        std::sort(searches.begin(), searches.end(),
                  [](const AdaptableSearch& lhs, const AdaptableSearch& rhs) {
                      for (uint16_t i = 0; i < lhs.getNumParts(); i++) {
                          if (lhs.getPart(i) != rhs.getPart(i)) {
                              return lhs.getPart(i) < rhs.getPart(i);
                          }
                      }

                      // order is equal return based on id
                      return lhs.getIndex() < rhs.getIndex();
                  });

        // reset idxes
        for (size_t s = 0; s < searches.size(); s++) {
            searches[s].setIdx(s);
        }
    }

    /**
     * @brief Set the U[1] value of the search the specified index to 2 and set
     * all other U[1] values to 1.
     *
     * This function sets the U[1] value of the search at the specified index to
     * 2 and sets all other U[1] values to 1. It is intended to be called when
     * the scheme is sorted on the first part.
     *
     * @param index The index of the search for which the U[1] value needs to
     * be 2.
     */
    void setUAtIndextoTwoAllOthers1(uint16_t index) {

        if (index >= searches.size()) {
            throw std::runtime_error("Critical search index (" +
                                     std::to_string(index) + ") is too high");
        }

        for (size_t i = 0; i < searches.size(); i++) {
            auto value = (index == i) ? 2 : 1;
            searches[i].setU(1, value);
        }
    }

    /**
     * @brief Recalculate the initial coverage of the search scheme.
     *
     * This function recalculates the initial coverage of the search scheme. It
     * is used to update the coverage after externally modifying the searches.
     */
    void recalculateInitialCoverage() {
        setInitialCoverage();
    }

    /**
     * @brief Set the searches in the scheme.
     *
     * This function sets the searches in the scheme to the specified vector of
     * adaptable searches. It also updates the indexes of the searches
     * accordingly and calculates the initial coverages.
     *
     * @param nSearches The new vector of adaptable searches.
     */
    void setSearches(std::vector<AdaptableSearch>& nSearches) {
        searches = nSearches;
        auto idx = 0;
        for (auto& s : searches) {
            s.setIdx(idx++);
        }
        setInitialCoverage();
    }
};

#endif
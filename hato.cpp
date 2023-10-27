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

#include <cstdint>
#include <iostream>
#include <string>

using namespace std;

// ============================================================================
// SCHEME EXTRACTION FROM FOLDER
// ============================================================================
/**
 * Extracts a vector of unsigned 16-bit integers from a string representation.
 * The string should be in the format of a vector, enclosed in (curly) brackets
 * and separated by commas.
 *
 * @param vectorString The string representation of the vector.
 * @param vector A reference to a vector to store the extracted values.
 * @throws std::runtime_error If the input string is not a valid vector.
 */
void getVector(const string& vectorString, vector<uint16_t>& vector) {

    if (vectorString.size() < 2) {
        throw runtime_error(vectorString +
                            " is not a valid vector for a search");
    }
    string bracketsRemoved = vectorString.substr(1, vectorString.size() - 2);

    stringstream ss(bracketsRemoved);
    string token;
    while (getline(ss, token, ',')) {
        vector.emplace_back(stoull(token));
    }
}

/**
 * Creates a Search object based on the input line and index.
 * The line should contain three vectors: order, lower bound, and upper bound.
 *
 * @param line The input line containing the vectors.
 * @param idx The index of the Search object.
 * @return A Search object constructed from the input line.
 * @throws std::runtime_error If the input line does not contain three vectors,
 * if any vector is empty, or if one of the vectors is not valid.
 */
Search makeSearch(const string& line, uint16_t idx) {
    stringstream ss;
    ss << line;

    string order, lower, upper;

    ss >> order >> lower >> upper;

    if (order.empty() || lower.empty() || upper.empty()) {
        throw runtime_error("A search should have 3 vectors: order, "
                            "lower bound and upper bound!");
    }

    vector<uint16_t> orderV;
    getVector(order, orderV);

    vector<uint16_t> L;
    getVector(lower, L);

    vector<uint16_t> U;
    getVector(upper, U);

    return Search::makeSearch(orderV, L, U, idx);
}

/**
 * Parses search schemes from a file and constructs a vector of Search objects.
 * Each search consists of three vectors: order, lower bound, and upper
 * bound. The scheme is read from a file in the specified folder based
 * on the number of errors.
 *
 * @param folder The folder containing the search scheme file.
 * @param k The number of errors.
 * @return A vector of Search objects representing the search scheme.
 * @throws std::runtime_error If the search scheme file is not found, has an
 * invalid format, or is empty of if the individual searches do not have the
 * same number of parts.
 */
vector<Search> getSearches(const string& folder, uint16_t k) {
    ifstream ifs(folder + +"/" + to_string(k) + "/searches.txt");
    if (!ifs) {
        throw runtime_error("No search scheme found in folder " + folder +
                            " with " + to_string(k) + " errors!");
    }

    vector<Search> v;
    uint16_t idx = 0;
    string line;
    while (getline(ifs, line)) {
        v.push_back(makeSearch(line, idx++));
    }

    if (v.size() == 0) {
        throw runtime_error("Empty search scheme");
    }

    uint16_t p = v.front().getNumParts();
    for (const auto& s : v) {
        if (p != s.getNumParts()) {
            throw runtime_error("Not all searches have same number of parts!");
        }
    }

    return v;
}

// ============================================================================
// SCHEME GENERATION
// ============================================================================

/**
 * @brief Generates pigeonhole searches based on a parameter k.
 *
 * This function generates a vector of Search objects representing pigeonhole
 * searches.
 *
 * @param k An integer parameter.
 * @return A vector of Search objects representing the generated pigeonhole
 * searches.
 */
vector<Search> generatePigeonholeSearches(uint16_t k,
                                          uint16_t criticalSearchIdx) {
    vector<Search> searches;
    searches.reserve(k + 1);

    bool switchOrder = (k % 2 == 0 && criticalSearchIdx < k / 2);

    for (uint16_t s = 0; s < k + 1; s++) {
        // create order
        vector<uint16_t> order;
        // order from s to k (inclusive) and then from s - 1 to 0
        for (uint16_t i = s; i <= k; i++) {
            order.emplace_back(i);
        }
        for (uint16_t i = s; i-- > 0;) {
            order.emplace_back(i);
        }

        if (switchOrder) {
            for (uint16_t i = 0; i < k + 1; i++) {
                order[i] = k - order[i];
            }
        }
        // create upper: 0kkkk
        vector<uint16_t> upper = {0};
        for (uint16_t i = 1; i < k + 1; i++) {
            upper.emplace_back(k);
        }
        // create lower: 00000
        vector<uint16_t> lower(k + 1, 0);

        searches.emplace_back(Search::makeSearch(order, lower, upper, k - s));
    }
    return searches;
}

// ============================================================================
// ADAPTABLE SCHEME GENERATION
// ============================================================================
/**
 * @brief Creates an adaptable scheme based on a vector of searches and a
 * parameter k.
 *
 * This function generates the error distributions and creates an adaptable
 * scheme using the provided searches and error distributions.
 *
 * @param searches A vector of Search objects representing the searches.
 * @param k An integer parameter.
 * @return An AdaptableScheme object representing the created scheme.
 */
AdaptableScheme createScheme(vector<Search> searches, uint16_t k) {

    uint16_t p = searches.front().getNumParts();

    // Generate the error distributions
    cout << "Generating error distributions" << endl;
    const auto& eDistributions = generateErrorDistributions(k, p);

    cout << eDistributions.size() << " error distributions generated" << endl;

    // Create the adaptable scheme
    AdaptableScheme scheme(searches, eDistributions);

    return scheme;
}

/**
 * @brief Creates an adaptable scheme based on a parameter k and a folder.
 *
 * This function creates an adaptable scheme using either pigeonhole searches
 * or searches loaded from a path. If the provided  path is empty, it uses the
 * pigeonhole principle to generate searches.
 *
 * @param k An integer parameter.
 * @param path A string representing the paht from which to load searches
 * (can be empty string).
 * @param criticalSearch An integer parameter representing the critical search
 * (optional).
 * @return An AdaptableScheme object representing the created scheme.
 */
AdaptableScheme createScheme(uint16_t k, const string& path,
                             int criticalSearch = 0) {
    vector<Search> searches;
    // get the searches
    if (path.empty()) {
        cout << "No search scheme specified, using pigeonhole principle"
             << endl;
        searches = generatePigeonholeSearches(k, criticalSearch);
    } else {
        searches = getSearches(path, k);
    }

    return createScheme(searches, k);
}

// ============================================================================
// SOLVERS
// ============================================================================

/** @brief Adapts greedily for a given number of errors and search scheme.
 *
 * This function applies the greedy adaptation algorithm to a specified number
 * of errors and search scheme. It updates the scheme by greedily selecting the
 * best adaptations. It also prints the resulting scheme and the coverages.
 *
 * @param scheme The adaptable scheme to start from. (will be updated)
 * @param k The number of errors.
 * @param criticalSearchIdx The index of the critical search (optional, default:
 * 0). This value is ignored if k is odd or if the number of parts is not equal
 * to k
 * + 1.
 */
void solveGreedy(AdaptableScheme& scheme, int k, int criticalSearchIdx = 0) {
    // sort scheme on order
    scheme.sortOnFirstPart();
    if (k % 2 == 0 && scheme.getNumParts() == k + 1) {

        // set U[1] == 2 for the critical search
        scheme.setUAtIndextoTwoAllOthers1(criticalSearchIdx);
    }
    // recalculate coverage
    scheme.recalculateInitialCoverage();

    // greedily adapt the scheme
    scheme.adapt();

    cout << scheme << endl;
    scheme.printCoverages();
}
/**
 * @brief Solves the CPLEX model for a given number of errors and search scheme.
 *
 * This function solves the CPLEX model for a specified number of errors and
 * search scheme. It invokes the CPLEX solver to find the optimal solution. It
 * also prints the resulting scheme and the coverages.
 *
 * @param scheme the scheme to start from (will be updated)
 * @param k The number of errors.
 * @param criticalSearch The index of the critical search (optional, default:
 * 0).
 */
void solveMILP(AdaptableScheme& scheme, int k, int criticalSearch = 0) {

    // Solve the cplex model
    IloEnv env;
    try {
        CPLEXModel cplex(env, k, scheme, criticalSearch);
        cplex.solveModel();

        vector<AdaptableSearch> optSearches = cplex.getScheme().getSearches();
        scheme.setSearches(optSearches);

        cout << scheme << endl;
        scheme.printCoverages();

    }

    catch (IloException& e) {
        cerr << "Error: " << e << endl;
        e.end();
    } catch (runtime_error& e) {
        cerr << "Error: " << e.what() << endl;

    } catch (...) {
        cerr << "Unknown exception caught!" << endl;
    }
}

// ============================================================================
// FUNCTIONS RELATED TO EXPECTED NUMBER OF NODES
// ============================================================================

/**
 * @brief Calculates the number of nodes at depth l in the tree where d errors
 * have been encountered.
 *
 * @param l The depth in the tree
 * @param Ux The upper bounds vector
 * @param Lx The lower bounds vector
 * @param sigma The size of the alphabet
 * @param d The number of errors encountered before this level
 * @return The number of encountered errors in levels 0 to l (inclusive)
 */
uint64_t nodesld(uint16_t l, const vector<uint16_t>& Ux,
                 const vector<uint16_t>& Lx, uint16_t sigma, uint16_t d) {

    if (l >= 1 && Lx[l - 1] <= d && d <= Ux[l - 1]) {
        return nodesld(l - 1, Ux, Lx, sigma, d) +
               (sigma - 1) * nodesld(l - 1, Ux, Lx, sigma, d - 1);
    } else if (l == 0 && d == 0) {
        return 1;
    }
    return 0;
}

/**
 * @brief Calculates the total number of nodes for a given depth in the tree
 *
 * @param l The depth in the tree.
 * @param Ux The upper bounds vector
 * @param Lx The lower bounds vector
 * @param sigma The size of the alphabet
 * @return The total number of nodes number of nodes for a given depth in the
 * tree
 */
uint64_t nodesl(uint16_t l, const vector<uint16_t>& Ux,
                const vector<uint16_t>& Lx, uint16_t sigma) {

    uint64_t total = 0;
    for (uint16_t d = Lx[l - 1]; d <= Ux[l - 1]; d++) {
        auto ld = nodesld(l, Ux, Lx, sigma, d);
        total += ld;
    }
    return total;
}

/**
 * @brief Calculates the logarithm of a number in a given base.
 *
 * @param base The base of the logarithm
 * @param a The number for which logarithm is calculated
 * @return The logarithm of the number in the given base
 */
double logb(uint64_t base, uint64_t a) {
    return log2(a) / log2(base);
}

/**
 * @brief Calculates the expected number of nodes for a given scheme.
 *
 * @param scheme The adaptable scheme
 * @param sigma The size of the alphabet
 * @param n The length of the reference text
 * @param m The length of the pattern
 * @param pLengths The lengths of the parts
 * @return The expected number of nodes
 */
double expectedNodes(const AdaptableScheme& scheme, const uint16_t sigma,
                     const uint64_t n, const uint16_t m,
                     const vector<uint16_t> pLengths) {

    double total = 0.0; // the total number of expected nodes
    const auto& searches = scheme.getSearches();
    uint16_t nParts = scheme.getNumParts(); // the number of parts in the scheme

    // assertions
    assert(pLengths.size() == nParts);
    assert(accumulate(pLengths.begin(), pLengths.end(), 0) == m);

    for (const auto& s : searches) {
        // create Lx and Ux
        vector<uint16_t> Lx(m, 0);
        vector<uint16_t> Ux(m, 0);

        uint16_t start = 0;

        for (uint16_t i = 0; i < nParts; i++) {
            uint16_t p = s.getPart(i);
            uint16_t pLength = pLengths[p];
            const uint16_t end = start + pLength;

            fill(Lx.begin() + start, Lx.begin() + end, s.getLowerBound(i));
            fill(Ux.begin() + start, Ux.begin() + end, s.getUpperBound(i));

            start = end;
        }

        uint16_t u = ceil(logb(sigma, n)) + 50;
        for (uint16_t l = 1; l < u; l++) {
            double factor = (1 - exp(-(n * 1.0 / pow(sigma, l))));
            if (factor < std::numeric_limits<double>::epsilon()) {
                break;
            }
            uint64_t nodes = nodesl(l, Ux, Lx, sigma);

            total += nodes * factor;
        }
    }
    return total;
}

/**
 * @brief Prints the expected number of nodes for a given adaptable scheme.
 *
 * @param scheme The adaptable scheme
 * @param k The value of k
 * @param m The length of the pattern
 * @param n The length of the reference text
 * @param sigma The size of the alphabet
 */
void printExpectedNodes(const AdaptableScheme& scheme, uint16_t k, uint16_t m,
                        uint64_t n, uint16_t sigma) {
    const auto& p = scheme.getNumParts();
    vector<uint16_t> partitioning(p, m / p);
    partitioning.back() =
        m - accumulate(partitioning.begin(), partitioning.end() - 1, 0);

    cout << "Expected nodes for n = " << n << ", m = " << m
         << ", sigma = " << sigma << ": " << endl;
    cout << expectedNodes(scheme, sigma, n, m, partitioning) << endl;
}

// ============================================================================
// MAIN
// ============================================================================

void printHelp() {
    cout << "Usage: program <mode> -k <number> [-s "
            "<path>] [-c <criticalSearchIndex>] [-g] [-n <number>] [-m "
            "<number>] [-sigma "
            "<number>]"
         << endl;
    cout << "Modes:" << endl;
    cout << "  solver     : Adapt the search scheme using either greedy or "
            "optimal algorithm"
         << endl;
    cout << "  expNodes   : Report the number of expected nodes using hamming "
            "distance"
         << endl;
    cout << "Arguments:" << endl;
    cout << "  -k <number> : The number of allowed errors" << endl;
    cout << "  -s <path>   : Path to the search scheme" << endl;
    cout << "  -c <criticalSearchIndex>  : The critical search index (only if "
            "k is even and p = "
            "k + 1). Must be even!"
         << endl;
    cout << "  -g          : Solve using only the greedy algorithm (optional)"
         << endl;
    cout << "  n           : Length of the reference text (optional, default: "
            "3000000000)"
         << endl;
    cout << "  m           : Length of the pattern (optional, default: 150)"
         << endl;
    cout << "  sigma       : Size of the alphabet (optional, default: 4)"
         << endl;
}

/**
 * @brief The main function of the program.
 *
 * This function serves as the entry point of the program. It reads the
 * command-line arguments, extracts the values for the mode, the number of
 * errors (k), search scheme folder path,  critical search index, and the greedy
 * flag. In mode "solver" the search scheme is adapted either greedily or
 * optimally. In mode "expNodes" the number of expected nodes using hamming
 * distance is reported.
 *
 * @param argc The number of command-line arguments.
 * @param argv An array of command-line argument strings.
 * @return The exit status of the program.
 */
int main(int argc, char* argv[]) {

    if (argc < 2) {
        printHelp();
        return 0;
    }

    string mode = argv[1]; // Mode: solver or expNodes

    if (mode == "help") {
        printHelp();
        return 0;
    }
    if (argc < 4) {
        cerr << "Not enough arguments. Usage: program <mode> -k <number> [-s "
                "<path>] [-c <critical>] [-g] [-n <number>] [-m <number>] "
                "[-sigma "
                "<number>]"
             << endl;
        printHelp();
        return 1;
    }

    uint16_t k = 0; // the number of allowed errors
    uint16_t critical =
        0;       // the critical search index (only if k even and p = k + 1)
    string path; // path to the search scheme
    bool greedy = false;     // flag to solve using only greedy algorithm
    uint64_t n = 3000000000; // Length of the reference text
    uint16_t m = 150;        // Length of the pattern
    uint16_t sigma = 4;      // Size of the alphabet

    for (int i = 2; i < argc; i++) {
        if (string(argv[i]) == "-k") {
            try {
                k = stoi(argv[i + 1]);
            } catch (...) {
                cerr << "Invalid value for -k argument. Must be a number."
                     << endl;
                printHelp();
                return 1;
            }
        } else if (string(argv[i]) == "-s") {
            try {
                path = argv[i + 1];
            } catch (...) {
                cerr << "-s parameter takes an argument!" << endl;
                printHelp();
                return 1;
            }

        } else if (string(argv[i]) == "-c") {
            try {
                critical = stoi(argv[i + 1]);
                if (critical % 2 == 1) {
                    cerr << "Invalid value for -c argument. Must be a even."
                         << endl;
                    printHelp();
                    return 1;
                }
            } catch (...) {
                cerr << "Invalid value for -c argument. Must be a number."
                     << endl;
                printHelp();
                return 1;
            }
        } else if (string(argv[i]) == "-g") {
            greedy = true;
        } else if (i + 2 < argc) {
            if (string(argv[i]) == "-n") {
                try {
                    n = stoul(argv[i + 1]);
                } catch (...) {
                    cerr << "Invalid value for -n argument. Must be a positive "
                            "integer."
                         << endl;
                    printHelp();
                    return 1;
                }
            } else if (string(argv[i]) == "-m") {
                try {
                    m = stoul(argv[i + 1]);
                } catch (...) {
                    cerr << "Invalid value for -m argument. Must be a positive "
                            "integer."
                         << endl;
                    printHelp();
                    return 1;
                }
            } else if (string(argv[i]) == "-sigma") {
                try {
                    sigma = stoul(argv[i + 1]);
                } catch (...) {
                    cerr << "Invalid value for -sigma argument. Must be a "
                            "positive integer."
                         << endl;
                    printHelp();
                    return 1;
                }
            }
        }
    }
    if (k == 0) {
        cerr << "Invalid value for -k argument. Must be a positive number."
             << endl;
        return 1;
    } else if (k > 13) {
        cerr << "Invalid value for -k argument. Must be no greater than 13."
             << endl;
        return 1;
    }

    auto scheme = createScheme(k, path, critical);

    if (mode == "solver") { // Mode 1: Solver
        if (greedy) {
            solveGreedy(scheme, k, critical);
        } else {
            solveMILP(scheme, k, critical);
        }
        printExpectedNodes(scheme, k, m, n, sigma);

    } else if (mode == "expNodes") { // Mode 2: ExpNodes
        if (path.empty()) {
            cerr << "Missing path argument (-s) for expNodes mode." << endl;
            printHelp();
            return 1;
        }
        printExpectedNodes(scheme, k, m, n, sigma);

    } else {
        cerr << "Invalid mode. Available modes: solver, expNodes." << endl;
        return 1;
    }

    return 0;
}

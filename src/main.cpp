#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

struct Atom {
    int index{}; // 1-based index
    std::string element;
};

struct Bond {
    int from{}; // 1-based atom index
    int to{};   // 1-based atom index
    int order{};
};

struct MoleculeGraph {
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
    std::vector<std::vector<int>> adjacency; // adjacency list of atom indices (1-based)
};

struct RepeatedFragment {
    std::string canonical;
    std::vector<int> edge_indices;
    int multiplicity{};
};

int atomic_number(const std::string &symbol) {
    static const std::unordered_map<std::string, int> kAtomicNumbers = {
        {"H", 1},   {"He", 2},  {"Li", 3},  {"Be", 4},  {"B", 5},   {"C", 6},
        {"N", 7},   {"O", 8},   {"F", 9},   {"Ne", 10}, {"Na", 11}, {"Mg", 12},
        {"Al", 13}, {"Si", 14}, {"P", 15}, {"S", 16},  {"Cl", 17}, {"Ar", 18},
        {"K", 19},  {"Ca", 20}, {"Sc", 21}, {"Ti", 22}, {"V", 23},  {"Cr", 24},
        {"Mn", 25}, {"Fe", 26}, {"Co", 27}, {"Ni", 28}, {"Cu", 29}, {"Zn", 30},
        {"Ga", 31}, {"Ge", 32}, {"As", 33}, {"Se", 34}, {"Br", 35}, {"Kr", 36},
        {"Rb", 37}, {"Sr", 38}, {"Y", 39},  {"Zr", 40}, {"Nb", 41}, {"Mo", 42},
        {"Tc", 43}, {"Ru", 44}, {"Rh", 45}, {"Pd", 46}, {"Ag", 47}, {"Cd", 48},
        {"In", 49}, {"Sn", 50}, {"Sb", 51}, {"Te", 52}, {"I", 53},  {"Xe", 54},
        {"Cs", 55}, {"Ba", 56}, {"La", 57}, {"Ce", 58}, {"Pr", 59}, {"Nd", 60},
        {"Pm", 61}, {"Sm", 62}, {"Eu", 63}, {"Gd", 64}, {"Tb", 65}, {"Dy", 66},
        {"Ho", 67}, {"Er", 68}, {"Tm", 69}, {"Yb", 70}, {"Lu", 71}, {"Hf", 72},
        {"Ta", 73}, {"W", 74},  {"Re", 75}, {"Os", 76}, {"Ir", 77}, {"Pt", 78},
        {"Au", 79}, {"Hg", 80}, {"Tl", 81}, {"Pb", 82}, {"Bi", 83}, {"Po", 84},
        {"At", 85}, {"Rn", 86}, {"Fr", 87}, {"Ra", 88}, {"Ac", 89}, {"Th", 90},
        {"Pa", 91}, {"U", 92},  {"Np", 93}, {"Pu", 94}, {"Am", 95}, {"Cm", 96},
        {"Bk", 97}, {"Cf", 98}, {"Es", 99}, {"Fm", 100}, {"Md", 101}, {"No", 102},
        {"Lr", 103}, {"Rf", 104}, {"Db", 105}, {"Sg", 106}, {"Bh", 107}, {"Hs", 108},
        {"Mt", 109}, {"Ds", 110}, {"Rg", 111}, {"Cn", 112}, {"Nh", 113}, {"Fl", 114},
        {"Mc", 115}, {"Lv", 116}, {"Ts", 117}, {"Og", 118}};

    auto it = kAtomicNumbers.find(symbol);
    if (it == kAtomicNumbers.end()) {
        return 0;
    }
    return it->second;
}

std::string format_atom_for_canonical(const Atom &atom) {
    int number = atomic_number(atom.element);
    if (number > 0) {
        std::ostringstream oss;
        oss << 'Z' << std::setw(3) << std::setfill('0') << number;
        return oss.str();
    }
    return "EL:" + atom.element;
}

std::string canonical_subgraph_string(const MoleculeGraph &graph,
                                      const std::vector<int> &edge_indices) {
    std::vector<std::string> encoded_edges;
    encoded_edges.reserve(edge_indices.size());

    for (int edge_index : edge_indices) {
        if (edge_index <= 0 || static_cast<std::size_t>(edge_index) > graph.bonds.size()) {
            throw std::runtime_error("Edge index is out of range when encoding subgraph.");
        }

        const Bond &bond = graph.bonds[static_cast<std::size_t>(edge_index) - 1];
        if (bond.from <= 0 || bond.to <= 0 ||
            static_cast<std::size_t>(bond.from) > graph.atoms.size() ||
            static_cast<std::size_t>(bond.to) > graph.atoms.size()) {
            throw std::runtime_error("Bond references an atom index outside the allowed range.");
        }

        std::string left = format_atom_for_canonical(graph.atoms[static_cast<std::size_t>(bond.from) - 1]);
        std::string right = format_atom_for_canonical(graph.atoms[static_cast<std::size_t>(bond.to) - 1]);

        if (right < left) {
            std::swap(left, right);
        }

        encoded_edges.push_back(left + "-" + right + ":" + std::to_string(bond.order));
    }

    std::sort(encoded_edges.begin(), encoded_edges.end());

    std::ostringstream canonical;
    for (std::size_t i = 0; i < encoded_edges.size(); ++i) {
        if (i > 0) {
            canonical << '|';
        }
        canonical << encoded_edges[i];
    }

    return canonical.str();
}

std::vector<std::vector<int>> enumerate_connected_subgraphs(const MoleculeGraph &graph) {
    std::vector<std::vector<int>> results;
    if (graph.bonds.empty()) {
        return results;
    }

    const std::size_t atom_count = graph.atoms.size();

    auto dfs = [&](auto &&self, int last_edge, std::vector<int> &edge_set,
                   std::vector<bool> &active_atoms) -> void {
        for (std::size_t next = static_cast<std::size_t>(last_edge + 1);
             next < graph.bonds.size(); ++next) {
            const Bond &bond = graph.bonds[next];

            bool touches_existing = active_atoms[bond.from] || active_atoms[bond.to];
            if (!touches_existing) {
                continue; // adding this edge would create a disconnected component
            }

            bool added_from = !active_atoms[bond.from];
            bool added_to = !active_atoms[bond.to];

            active_atoms[bond.from] = true;
            active_atoms[bond.to] = true;
            edge_set.push_back(static_cast<int>(next) + 1); // store 1-based edge index

            results.push_back(edge_set);
            self(self, static_cast<int>(next), edge_set, active_atoms);

            edge_set.pop_back();
            if (added_from) {
                active_atoms[bond.from] = false;
            }
            if (added_to) {
                active_atoms[bond.to] = false;
            }
        }
    };

    for (std::size_t i = 0; i < graph.bonds.size(); ++i) {
        const Bond &start_bond = graph.bonds[i];
        std::vector<bool> active_atoms(atom_count + 1, false); // 1-based atoms
        active_atoms[start_bond.from] = true;
        active_atoms[start_bond.to] = true;

        std::vector<int> edge_set = {static_cast<int>(i) + 1};
        results.push_back(edge_set);
        dfs(dfs, static_cast<int>(i), edge_set, active_atoms);
    }

    return results;
}

std::vector<RepeatedFragment> collect_repeated_fragments(const MoleculeGraph &graph) {
    const auto subgraphs = enumerate_connected_subgraphs(graph);

    std::vector<std::string> canonical_signatures;
    canonical_signatures.reserve(subgraphs.size());

    std::unordered_map<std::string, int> multiplicity;
    for (const auto &edges : subgraphs) {
        std::string signature = canonical_subgraph_string(graph, edges);
        canonical_signatures.push_back(signature);
        ++multiplicity[signature];
    }

    std::vector<RepeatedFragment> fragments;
    std::unordered_set<std::string> emitted;
    for (std::size_t i = 0; i < subgraphs.size(); ++i) {
        const std::string &signature = canonical_signatures[i];
        auto count = multiplicity[signature];
        if (count < 2 || emitted.count(signature) > 0) {
            continue;
        }

        std::vector<int> edges = subgraphs[i];
        std::sort(edges.begin(), edges.end());

        fragments.push_back({signature, edges, count});
        emitted.insert(signature);
    }

    return fragments;
}

void print_repeated_fragments(const std::vector<RepeatedFragment> &fragments) {
    std::cout << "Fragments with multiplicity >= 2: " << fragments.size() << '\n';
    for (std::size_t i = 0; i < fragments.size(); ++i) {
        const auto &fragment = fragments[i];
        std::cout << "  Fragment " << (i + 1) << " (multiplicity " << fragment.multiplicity
                  << "):\n";
        std::cout << "    Canonical: " << fragment.canonical << '\n';
        std::cout << "    Example edges:";
        for (int edge_index : fragment.edge_indices) {
            std::cout << ' ' << edge_index;
        }
        std::cout << '\n';
    }
}

bool is_subset(const std::vector<int> &subset, const std::vector<int> &superset) {
    std::size_t i = 0;
    std::size_t j = 0;
    while (i < subset.size() && j < superset.size()) {
        if (subset[i] == superset[j]) {
            ++i;
            ++j;
        } else if (subset[i] > superset[j]) {
            ++j;
        } else {
            return false;
        }
    }
    return i == subset.size();
}

std::vector<int> subtract_edges(const std::vector<int> &source, const std::vector<int> &to_remove) {
    std::vector<int> remaining;
    remaining.reserve(source.size());

    std::size_t i = 0;
    std::size_t j = 0;
    while (i < source.size()) {
        if (j < to_remove.size() && source[i] == to_remove[j]) {
            ++i;
            ++j;
        } else if (j < to_remove.size() && source[i] > to_remove[j]) {
            ++j;
        } else {
            remaining.push_back(source[i]);
            ++i;
        }
    }

    return remaining;
}

int compute_min_ma(const MoleculeGraph &graph, const std::vector<int> &subgraph_edges,
                   const std::vector<RepeatedFragment> &fragments,
                   std::unordered_map<std::string, int> &memo) {
    std::vector<int> sorted_edges = subgraph_edges;
    std::sort(sorted_edges.begin(), sorted_edges.end());

    std::string key = canonical_subgraph_string(graph, sorted_edges);
    auto it = memo.find(key);
    if (it != memo.end()) {
        return it->second;
    }

    int naive_ma = static_cast<int>(sorted_edges.size());
    if (!sorted_edges.empty()) {
        naive_ma = static_cast<int>(sorted_edges.size()) - 1;
    }

    int best_ma = naive_ma;

    for (const auto &fragment : fragments) {
        if (fragment.edge_indices.size() >= sorted_edges.size()) {
            continue; // cannot split or would be identical
        }

        if (!is_subset(fragment.edge_indices, sorted_edges)) {
            continue;
        }

        int fragment_ma = compute_min_ma(graph, fragment.edge_indices, fragments, memo);

        std::vector<int> remaining = subtract_edges(sorted_edges, fragment.edge_indices);
        int remaining_ma = 0;
        int combined_cost = fragment_ma;
        if (!remaining.empty()) {
            remaining_ma = compute_min_ma(graph, remaining, fragments, memo);
            combined_cost = fragment_ma + remaining_ma + 1;
        }

        if (combined_cost < best_ma) {
            best_ma = combined_cost;
        }
    }

    memo[key] = best_ma;
    return best_ma;
}

int parse_fixed_width_int(const std::string &line, std::size_t offset, std::size_t width) {
    if (line.size() < offset + width) {
        throw std::runtime_error("Line is too short to parse expected integer field.");
    }
    std::string field = line.substr(offset, width);

    try {
        return std::stoi(field);
    } catch (const std::exception &) {
        throw std::runtime_error("Failed to parse integer value from field: '" + field + "'.");
    }
}

std::string parse_fixed_width_string(const std::string &line, std::size_t offset, std::size_t width) {
    if (line.size() < offset + width) {
        throw std::runtime_error("Line is too short to parse expected string field.");
    }
    std::string field = line.substr(offset, width);

    std::size_t first = field.find_first_not_of(' ');
    if (first == std::string::npos) {
        return "";
    }
    std::size_t last = field.find_last_not_of(' ');
    return field.substr(first, last - first + 1);
}

MoleculeGraph parse_mol_v2000(std::istream &input) {
    std::vector<std::string> lines;
    for (std::string line; std::getline(input, line);) {
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        lines.push_back(line);
    }

    if (lines.size() < 4) {
        throw std::runtime_error("Input does not contain enough lines for a valid molfile.");
    }

    const std::string &counts_line = lines[3];
    int atom_count = parse_fixed_width_int(counts_line, 0, 3);
    int bond_count = parse_fixed_width_int(counts_line, 3, 3);

    if (atom_count < 0 || bond_count < 0) {
        throw std::runtime_error("Atom or bond count cannot be negative.");
    }
    if (lines.size() < 4 + static_cast<std::size_t>(atom_count + bond_count)) {
        throw std::runtime_error("Input does not contain the declared number of atom and bond lines.");
    }

    MoleculeGraph graph;
    graph.atoms.reserve(static_cast<std::size_t>(atom_count));
    graph.bonds.reserve(static_cast<std::size_t>(bond_count));
    graph.adjacency.assign(static_cast<std::size_t>(atom_count) + 1, {}); // 1-based indexing

    // Atom block starts at line 4
    for (int i = 0; i < atom_count; ++i) {
        const std::string &atom_line = lines[4 + i];

        if (atom_line.size() < 34) {
            throw std::runtime_error("Atom line is too short: " + atom_line);
        }

        std::string element = parse_fixed_width_string(atom_line, 30, 4);
        if (element.empty()) {
            throw std::runtime_error("Atom symbol is missing in line: " + atom_line);
        }

        graph.atoms.push_back({i + 1, element});
    }

    // Bond block starts after atoms
    std::size_t bond_start = 4 + static_cast<std::size_t>(atom_count);
    for (int i = 0; i < bond_count; ++i) {
        const std::string &bond_line = lines[bond_start + static_cast<std::size_t>(i)];

        int from = parse_fixed_width_int(bond_line, 0, 3);
        int to = parse_fixed_width_int(bond_line, 3, 3);
        int order = parse_fixed_width_int(bond_line, 6, 3);

        if (from <= 0 || to <= 0 || from > atom_count || to > atom_count) {
            throw std::runtime_error("Bond references an atom index outside the allowed range.");
        }

        graph.bonds.push_back({from, to, order});
        graph.adjacency[static_cast<std::size_t>(from)].push_back(to);
        graph.adjacency[static_cast<std::size_t>(to)].push_back(from);
    }

    return graph;
}

void print_graph(const MoleculeGraph &graph) {
    std::cout << "Atom count: " << graph.atoms.size() << '\n';
    std::cout << "Bond count: " << graph.bonds.size() << '\n';

    std::cout << "Atoms:\n";
    for (const auto &atom : graph.atoms) {
        std::cout << "  (" << atom.index << ") " << atom.element << '\n';
    }

    std::cout << "Bonds:\n";
    for (const auto &bond : graph.bonds) {
        std::cout << "  (" << bond.from << ") -" << bond.order << "- (" << bond.to << ")\n";
    }

    std::cout << "Adjacency list (1-based indices):\n";
    for (std::size_t i = 1; i < graph.adjacency.size(); ++i) {
        std::cout << "  " << i << ":";
        for (int neighbor : graph.adjacency[i]) {
            std::cout << ' ' << neighbor;
        }
        std::cout << '\n';
    }
}

int main() {
    try {
        MoleculeGraph graph = parse_mol_v2000(std::cin);
        print_graph(graph);
        const auto fragments = collect_repeated_fragments(graph);
        print_repeated_fragments(fragments);

        std::vector<int> whole_graph_edges;
        whole_graph_edges.reserve(graph.bonds.size());
        for (std::size_t i = 0; i < graph.bonds.size(); ++i) {
            whole_graph_edges.push_back(static_cast<int>(i) + 1);
        }

        std::unordered_map<std::string, int> memo;
        int minimal_ma = compute_min_ma(graph, whole_graph_edges, fragments, memo);
        std::cout << "Minimal MA (split-branch): " << minimal_ma << '\n';
    } catch (const std::exception &ex) {
        std::cerr << "Failed to parse molfile: " << ex.what() << '\n';
        return 1;
    }
    return 0;
}

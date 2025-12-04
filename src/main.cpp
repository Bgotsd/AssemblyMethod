#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
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

void print_connected_subgraphs(const MoleculeGraph &graph) {
    const auto subgraphs = enumerate_connected_subgraphs(graph);
    std::cout << "Connected edge subgraphs: " << subgraphs.size() << '\n';
    for (std::size_t i = 0; i < subgraphs.size(); ++i) {
        std::cout << "  Subgraph " << (i + 1) << ':';
        for (int edge_index : subgraphs[i]) {
            std::cout << ' ' << edge_index;
        }
        std::cout << '\n';
    }
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
        print_connected_subgraphs(graph);
    } catch (const std::exception &ex) {
        std::cerr << "Failed to parse molfile: " << ex.what() << '\n';
        return 1;
    }
    return 0;
}

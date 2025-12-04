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
    } catch (const std::exception &ex) {
        std::cerr << "Failed to parse molfile: " << ex.what() << '\n';
        return 1;
    }
    return 0;
}

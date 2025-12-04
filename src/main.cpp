#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

struct Atom {
    int index{};            // 1-based
    std::string element;    // element symbol
};

struct Bond {
    int index{};            // 1-based
    int from{};             // 1-based atom index
    int to{};               // 1-based atom index
    int order{};            // bond order
};

struct MoleculeGraph {
    std::vector<Atom> atoms;                         // 1-based indices stored in index field
    std::vector<Bond> bonds;                         // 1-based indices stored in index field
    std::vector<std::vector<int>> adjacency;         // adjacency list of atom indices (1-based)
    std::unordered_map<long long, int> bond_lookup;  // key = (min<<32)|max -> bond index
};

// Motif abstractions
struct MotifType {
    std::string label;      // canonical label
    MoleculeGraph motif_graph;
    int local_ma{};
};

struct MotifEmbedding {
    const MotifType* type{};
    std::vector<int> atom_map;   // motif atom (1-based) -> host atom (1-based)
    std::vector<int> bond_indices; // host bond indices in this occurrence
};

struct MotifDatabase {
    std::vector<MotifType> types;
    std::vector<MotifEmbedding> embeddings;
    std::vector<std::vector<int>> type_to_embedding_indices; // parallel to types
};

int atomic_number(const std::string& symbol) {
    static const std::unordered_map<std::string, int> kAtomicNumbers = {
        {"H", 1},   {"He", 2},  {"Li", 3},  {"Be", 4},  {"B", 5},   {"C", 6},
        {"N", 7},   {"O", 8},   {"F", 9},   {"Ne", 10}, {"Na", 11}, {"Mg", 12},
        {"Al", 13}, {"Si", 14}, {"P", 15},  {"S", 16},  {"Cl", 17}, {"Ar", 18},
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

std::string format_atom_label(const std::string& element) {
    int number = atomic_number(element);
    if (number > 0) {
        std::ostringstream oss;
        oss << 'Z' << std::setw(3) << std::setfill('0') << number;
        return oss.str();
    }
    return "EL:" + element;
}

// Helper to build bond lookup key
long long bond_key(int a, int b) {
    int x = std::min(a, b);
    int y = std::max(a, b);
    return (static_cast<long long>(x) << 32) | static_cast<unsigned long long>(y);
}

int parse_fixed_width_int(const std::string& line, std::size_t offset, std::size_t width) {
    if (line.size() < offset + width) {
        throw std::runtime_error("Line is too short to parse expected integer field.");
    }
    std::string field = line.substr(offset, width);
    try {
        return std::stoi(field);
    } catch (const std::exception&) {
        throw std::runtime_error("Failed to parse integer value from field: '" + field + "'.");
    }
}

std::string parse_fixed_width_string(const std::string& line, std::size_t offset, std::size_t width) {
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

MoleculeGraph parse_mol_v2000(std::istream& input) {
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

    const std::string& counts_line = lines[3];
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
    graph.adjacency.assign(static_cast<std::size_t>(atom_count) + 1, {});

    for (int i = 0; i < atom_count; ++i) {
        const std::string& atom_line = lines[4 + i];
        if (atom_line.size() < 34) {
            throw std::runtime_error("Atom line is too short: " + atom_line);
        }
        std::string element = parse_fixed_width_string(atom_line, 30, 4);
        if (element.empty()) {
            throw std::runtime_error("Atom symbol is missing in line: " + atom_line);
        }
        graph.atoms.push_back({i + 1, element});
    }

    std::size_t bond_start = 4 + static_cast<std::size_t>(atom_count);
    for (int i = 0; i < bond_count; ++i) {
        const std::string& bond_line = lines[bond_start + static_cast<std::size_t>(i)];
        int from = parse_fixed_width_int(bond_line, 0, 3);
        int to = parse_fixed_width_int(bond_line, 3, 3);
        int order = parse_fixed_width_int(bond_line, 6, 3);
        if (from <= 0 || to <= 0 || from > atom_count || to > atom_count) {
            throw std::runtime_error("Bond references an atom index outside the allowed range.");
        }
        int bond_index = i + 1;
        graph.bonds.push_back({bond_index, from, to, order});
        graph.adjacency[static_cast<std::size_t>(from)].push_back(to);
        graph.adjacency[static_cast<std::size_t>(to)].push_back(from);
        graph.bond_lookup[bond_key(from, to)] = bond_index;
    }

    return graph;
}

void print_graph(const MoleculeGraph& graph) {
    std::cout << "Atom count: " << graph.atoms.size() << '\n';
    std::cout << "Bond count: " << graph.bonds.size() << '\n';

    std::cout << "Atoms:\n";
    for (const auto& atom : graph.atoms) {
        std::cout << "  (" << atom.index << ") " << atom.element << '\n';
    }

    std::cout << "Bonds:\n";
    for (const auto& bond : graph.bonds) {
        std::cout << "  " << bond.index << ": (" << bond.from << ") -" << bond.order << "- ("
                  << bond.to << ")\n";
    }
}

// Enumerate connected subgraphs by edge sets up to a maximum number of bonds
const int MAX_BONDS_PER_MOTIF = 6;

std::vector<std::vector<int>> enumerate_connected_subgraphs(const MoleculeGraph& graph) {
    std::set<std::vector<int>> unique_sets;

    auto dfs = [&](auto&& self, int last_edge_idx, std::vector<int>& edge_set,
                   std::vector<bool>& atom_included) -> void {
        unique_sets.insert(edge_set);
        if (static_cast<int>(edge_set.size()) >= MAX_BONDS_PER_MOTIF) {
            return;
        }

        for (std::size_t j = static_cast<std::size_t>(last_edge_idx + 1); j < graph.bonds.size(); ++j) {
            const Bond& b = graph.bonds[j];
            if (!(atom_included[b.from] || atom_included[b.to])) {
                continue; // would disconnect
            }
            bool added_from = !atom_included[b.from];
            bool added_to = !atom_included[b.to];
            atom_included[b.from] = true;
            atom_included[b.to] = true;
            edge_set.push_back(b.index);
            self(self, static_cast<int>(j), edge_set, atom_included);
            edge_set.pop_back();
            if (added_from) atom_included[b.from] = false;
            if (added_to) atom_included[b.to] = false;
        }
    };

    for (std::size_t i = 0; i < graph.bonds.size(); ++i) {
        const Bond& start = graph.bonds[i];
        std::vector<int> edge_set = {start.index};
        std::vector<bool> atom_included(graph.atoms.size() + 1, false);
        atom_included[start.from] = atom_included[start.to] = true;
        dfs(dfs, static_cast<int>(i), edge_set, atom_included);
    }

    return std::vector<std::vector<int>>(unique_sets.begin(), unique_sets.end());
}

MoleculeGraph build_motif_graph(const MoleculeGraph& host, const std::vector<int>& edge_indices) {
    std::unordered_map<int, int> atom_map; // host atom -> motif atom (1-based)
    MoleculeGraph motif;

    for (int edge_index : edge_indices) {
        const Bond& b = host.bonds[static_cast<std::size_t>(edge_index) - 1];
        if (!atom_map.count(b.from)) {
            int idx = static_cast<int>(atom_map.size()) + 1;
            atom_map[b.from] = idx;
            motif.atoms.push_back({idx, host.atoms[static_cast<std::size_t>(b.from) - 1].element});
        }
        if (!atom_map.count(b.to)) {
            int idx = static_cast<int>(atom_map.size()) + 1;
            atom_map[b.to] = idx;
            motif.atoms.push_back({idx, host.atoms[static_cast<std::size_t>(b.to) - 1].element});
        }
    }

    motif.adjacency.assign(motif.atoms.size() + 1, {});
    motif.bonds.reserve(edge_indices.size());
    motif.bond_lookup.clear();

    for (std::size_t i = 0; i < edge_indices.size(); ++i) {
        const Bond& b = host.bonds[static_cast<std::size_t>(edge_indices[i]) - 1];
        int from = atom_map[b.from];
        int to = atom_map[b.to];
        int order = b.order;
        int idx = static_cast<int>(i) + 1;
        motif.bonds.push_back({idx, from, to, order});
        motif.adjacency[static_cast<std::size_t>(from)].push_back(to);
        motif.adjacency[static_cast<std::size_t>(to)].push_back(from);
        motif.bond_lookup[bond_key(from, to)] = idx;
    }

    return motif;
}

std::string canonical_label_for_motif(const MoleculeGraph& motif) {
    int n = static_cast<int>(motif.atoms.size());
    if (n == 0) return "";

    // Build adjacency matrix of bond orders
    std::vector<std::vector<int>> adj(n, std::vector<int>(n, 0));
    for (const auto& b : motif.bonds) {
        adj[b.from - 1][b.to - 1] = b.order;
        adj[b.to - 1][b.from - 1] = b.order;
    }

    std::vector<int> indices(n);
    for (int i = 0; i < n; ++i) indices[i] = i;

    std::string best;
    bool first = true;
    do {
        std::ostringstream oss;
        oss << "A:";
        for (int i = 0; i < n; ++i) {
            if (i > 0) oss << ',';
            oss << format_atom_label(motif.atoms[static_cast<std::size_t>(indices[i])].element);
        }
        oss << "|B:";
        for (int i = 0; i < n; ++i) {
            if (i > 0) oss << ';';
            for (int j = 0; j < n; ++j) {
                if (j > 0) oss << ',';
                oss << adj[indices[i]][indices[j]];
            }
        }
        std::string candidate = oss.str();
        if (first || candidate < best) {
            best = candidate;
            first = false;
        }
    } while (std::next_permutation(indices.begin(), indices.end()));

    return best;
}

int compute_naive_ma(int bond_count) { return std::max(0, bond_count - 1); }

// Embedding search
int host_bond_index(const MoleculeGraph& host, int a, int b) {
    auto it = host.bond_lookup.find(bond_key(a, b));
    if (it == host.bond_lookup.end()) return -1;
    return it->second;
}

void backtrack_embedding(const MoleculeGraph& host, const MotifType& type, int depth,
                         std::vector<int>& motif_to_host, std::vector<bool>& host_used,
                         std::vector<MotifEmbedding>& out) {
    int n = static_cast<int>(type.motif_graph.atoms.size());
    if (depth == n) {
        MotifEmbedding emb;
        emb.type = &type;
        emb.atom_map = motif_to_host;
        // collect bond indices
        std::unordered_set<int> bond_set;
        for (const auto& b : type.motif_graph.bonds) {
            int ha = motif_to_host[b.from - 1];
            int hb = motif_to_host[b.to - 1];
            int idx = host_bond_index(host, ha, hb);
            if (idx > 0) bond_set.insert(idx);
        }
        emb.bond_indices.assign(bond_set.begin(), bond_set.end());
        std::sort(emb.bond_indices.begin(), emb.bond_indices.end());
        out.push_back(std::move(emb));
        return;
    }

    // choose next motif atom
    int next = -1;
    for (int i = 0; i < n; ++i) {
        if (motif_to_host[i] == 0) {
            next = i;
            break;
        }
    }
    if (next == -1) return;

    const Atom& motif_atom = type.motif_graph.atoms[static_cast<std::size_t>(next)];
    int motif_degree = static_cast<int>(type.motif_graph.adjacency[motif_atom.index].size());

    for (const auto& host_atom : host.atoms) {
        if (host_used[host_atom.index]) continue;
        if (host_atom.element != motif_atom.element) continue;
        if (static_cast<int>(host.adjacency[host_atom.index].size()) < motif_degree) continue;

        // check consistency with already mapped neighbors
        bool ok = true;
        for (int neighbor : type.motif_graph.adjacency[motif_atom.index]) {
            int mapped = motif_to_host[neighbor - 1];
            if (mapped == 0) continue;
            int hb = host_bond_index(host, host_atom.index, mapped);
            if (hb < 0) {
                ok = false;
                break;
            }
            // check bond order
            int motif_bond_idx = host_bond_index(type.motif_graph, motif_atom.index, neighbor);
            const Bond& host_b = host.bonds[static_cast<std::size_t>(hb) - 1];
            const Bond& motif_b = type.motif_graph.bonds[static_cast<std::size_t>(motif_bond_idx) - 1];
            if (host_b.order != motif_b.order) {
                ok = false;
                break;
            }
        }
        if (!ok) continue;

        motif_to_host[next] = host_atom.index;
        host_used[host_atom.index] = true;
        backtrack_embedding(host, type, depth + 1, motif_to_host, host_used, out);
        host_used[host_atom.index] = false;
        motif_to_host[next] = 0;
    }
}

std::vector<MotifEmbedding> find_embeddings(const MoleculeGraph& host, const MotifType& type) {
    int n = static_cast<int>(type.motif_graph.atoms.size());
    std::vector<int> motif_to_host(static_cast<std::size_t>(n), 0);
    std::vector<bool> host_used(host.atoms.size() + 1, false);
    std::vector<MotifEmbedding> embeddings;
    backtrack_embedding(host, type, 0, motif_to_host, host_used, embeddings);
    return embeddings;
}

MotifDatabase build_motif_database(const MoleculeGraph& host) {
    MotifDatabase db;
    auto subgraphs = enumerate_connected_subgraphs(host);

    std::unordered_map<std::string, MoleculeGraph> label_to_graph;
    for (const auto& edges : subgraphs) {
        if (edges.empty() || static_cast<int>(edges.size()) > MAX_BONDS_PER_MOTIF) continue;
        MoleculeGraph motif_graph = build_motif_graph(host, edges);
        if (motif_graph.bonds.empty()) continue;
        std::string label = canonical_label_for_motif(motif_graph);
        if (!label_to_graph.count(label)) {
            label_to_graph[label] = motif_graph;
        }
    }

    // Create motif types
    for (const auto& kv : label_to_graph) {
        MotifType type;
        type.label = kv.first;
        type.motif_graph = kv.second;
        type.local_ma = compute_naive_ma(static_cast<int>(type.motif_graph.bonds.size()));
        db.types.push_back(std::move(type));
    }

    db.type_to_embedding_indices.resize(db.types.size());

    // Enumerate embeddings and filter by multiplicity >= 2
    std::vector<MotifType> filtered_types;
    std::vector<std::vector<MotifEmbedding>> per_type_embeddings;

    for (std::size_t i = 0; i < db.types.size(); ++i) {
        auto embeddings = find_embeddings(host, db.types[i]);
        if (embeddings.size() < 2) continue;
        per_type_embeddings.push_back(std::move(embeddings));
        filtered_types.push_back(std::move(db.types[i]));
    }

    db.types = std::move(filtered_types);
    db.type_to_embedding_indices.resize(db.types.size());

    for (std::size_t i = 0; i < db.types.size(); ++i) {
        const auto& type = db.types[i];
        auto embeddings = find_embeddings(host, type);
        for (auto& emb : embeddings) {
            int idx = static_cast<int>(db.embeddings.size());
            db.type_to_embedding_indices[i].push_back(idx);
            db.embeddings.push_back(std::move(emb));
        }
    }

    return db;
}

struct MotifBenefitInfo {
    const MotifType* type{};
    int embedding_count{};
    int benefit{};
    std::vector<int> covered_bonds; // unique bond indices from all embeddings
};

std::vector<MotifBenefitInfo> compute_benefits(const MotifDatabase& db) {
    std::vector<MotifBenefitInfo> infos;
    for (std::size_t i = 0; i < db.types.size(); ++i) {
        const MotifType& type = db.types[i];
        std::unordered_set<int> bond_set;
        for (int emb_idx : db.type_to_embedding_indices[i]) {
            const auto& emb = db.embeddings[static_cast<std::size_t>(emb_idx)];
            bond_set.insert(emb.bond_indices.begin(), emb.bond_indices.end());
        }
        int k = static_cast<int>(db.type_to_embedding_indices[i].size());
        int E = static_cast<int>(type.motif_graph.bonds.size());
        int cost_naive = k * std::max(0, E - 1);
        int cost_ideal = type.local_ma + (k - 1);
        int benefit = cost_naive - cost_ideal;
        MotifBenefitInfo info;
        info.type = &type;
        info.embedding_count = k;
        info.benefit = benefit;
        info.covered_bonds.assign(bond_set.begin(), bond_set.end());
        infos.push_back(std::move(info));
    }
    return infos;
}

int compute_motif_based_ma(const MoleculeGraph& host, const MotifDatabase& db) {
    int naive = compute_naive_ma(static_cast<int>(host.bonds.size()));
    auto infos = compute_benefits(db);
    std::vector<bool> bond_covered(host.bonds.size() + 1, false);

    std::sort(infos.begin(), infos.end(), [](const MotifBenefitInfo& a, const MotifBenefitInfo& b) {
        if (a.benefit != b.benefit) return a.benefit > b.benefit;
        return a.type->motif_graph.bonds.size() > b.type->motif_graph.bonds.size();
    });

    double effective_savings = 0.0;
    for (const auto& info : infos) {
        if (info.benefit <= 0) continue;
        int new_cover = 0;
        for (int b : info.covered_bonds) {
            if (!bond_covered[b]) new_cover++;
        }
        if (new_cover == 0) continue;
        // accept motif
        for (int b : info.covered_bonds) {
            bond_covered[b] = true;
        }
        // scale benefit by fraction of newly covered bonds
        double scale = static_cast<double>(new_cover) /
                        static_cast<double>(info.covered_bonds.empty() ? 1 : info.covered_bonds.size());
        effective_savings += static_cast<double>(info.benefit) * scale;
    }

    int improved = static_cast<int>(std::round(static_cast<double>(naive) - effective_savings));
    if (improved < 0) improved = 0;
    return improved;
}

void print_motif_summary(const MotifDatabase& db) {
    auto infos = compute_benefits(db);
    std::cout << "Motif types discovered: " << infos.size() << '\n';
    for (std::size_t i = 0; i < infos.size(); ++i) {
        const auto& info = infos[i];
        const auto& motif_graph = info.type->motif_graph;
        std::cout << "  Type " << (i + 1) << " label=" << info.type->label << '\n';
        std::cout << "    Atoms: " << motif_graph.atoms.size() << ", Bonds: "
                  << motif_graph.bonds.size() << '\n';
        std::cout << "    Embeddings: " << info.embedding_count << '\n';
        std::cout << "    Local MA: " << info.type->local_ma << '\n';
        std::cout << "    Benefit: " << info.benefit << '\n';
    }
}

int main() {
    try {
        MoleculeGraph graph = parse_mol_v2000(std::cin);
        print_graph(graph);

        MotifDatabase db = build_motif_database(graph);
        print_motif_summary(db);

        int naive_ma = compute_naive_ma(static_cast<int>(graph.bonds.size()));
        int approx_ma = compute_motif_based_ma(graph, db);

        std::cout << "Naive MA: " << naive_ma << '\n';
        std::cout << "Approximate motif-based MA: " << approx_ma << '\n';
    } catch (const std::exception& ex) {
        std::cerr << "Failed to parse molfile: " << ex.what() << '\n';
        return 1;
    }
    return 0;
}

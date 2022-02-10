//
// Created by Bernardo Clavijo (EI) on 24/02/2018.
//

#include "WorkSpace.hpp"
#include <sstream>


const sdgVersion_t WorkSpace::min_compat = 0x0003;

void WorkSpace::dump_to_disk(std::string filename, bool graph_only) {
    sdglib::OutputLog()<<"Dumping workspace to "<<filename<<std::endl;
    std::ofstream of(filename);

    //Magic number
    of.write((char *) &SDG_MAGIC, sizeof(SDG_MAGIC));
    of.write((char *) &SDG_VN, sizeof(SDG_VN));
    SDG_FILETYPE type(WS_FT);
    of.write((char *) &type, sizeof(type));


    uint64_t count;

    //dump operations
    count = journal.size();
    of.write((char *) &count, sizeof(count));
    for (const auto &j:journal){
        sdglib::write_string(of,j.name);
        sdglib::write_string(of,j.detail);
        sdglib::write_string(of,j.tool);
        of.write((char *) &j.timestamp, sizeof(j.timestamp));
        count = j.entries.size();
        of.write((char *) &count, sizeof(count));
        for (const auto &e: j.entries) {
            sdglib::write_string(of, e.detail);
        }
    }

    //dump main graph
    sdg.write(of);

    count = distance_graphs.size();
    of.write((char *) &count, sizeof(count));
    for (auto i=0; i < count; ++i) {
        distance_graphs[i].write(of);
    }

    if (graph_only){
        count=0;
        of.write((char *) &count,sizeof(count));
        of.write((char *) &count,sizeof(count));
        of.write((char *) &count,sizeof(count));
        of.write((char *) &count,sizeof(count));
        return;
    }

    //paired read datastores
    count=paired_reads_datastores.size();
    of.write((char *) &count,sizeof(count));
    for (auto i=0;i<count;++i){
        paired_reads_datastores[i].write(of);
    }

    //linker read datastores
    count=linked_reads_datastores.size();
    of.write((char *) &count,sizeof(count));
    for (auto i=0;i<count;++i){
        linked_reads_datastores[i].write(of);
    }

    //long read datastores
    count=long_reads_datastores.size();
    of.write((char *) &count,sizeof(count));
    for (auto i=0;i<count;++i){
        long_reads_datastores[i].write(of);
    }

    // Kmer counts, keep these ones at the end of the workspace to make the workspaces editable by exteding the final part of the file
    count = kmer_counters.size();
    of.write((char *) &count,sizeof(count));
    for (auto i = 0; i < count; i++) {
        kmer_counters[i].write(of);
    }
}

void WorkSpace::load_from_disk(std::string filename, bool log_only) {
    std::ifstream wsfile(filename);
    if (!wsfile.good()) {
        std::cerr << filename << " opening error: " << strerror(errno) << std::endl;
        throw std::runtime_error("Error opening " + filename);
    }

    //==== Check magic number ====
    sdgMagic_t magic;
    sdgVersion_t version;
    SDG_FILETYPE type;
    wsfile.read((char *) &magic, sizeof(magic));
    wsfile.read((char *) &version, sizeof(version));
    wsfile.read((char *) &type, sizeof(type));

    if (magic != SDG_MAGIC) {
        throw std::runtime_error("Magic number not present in WorkSpace file");
    }

    if (version < min_compat) {
        throw std::runtime_error("WorkSpace file version: " + std::to_string(version) + " is not compatible with " + std::to_string(min_compat));
    }

    if (type != WS_FT) {
        throw std::runtime_error("File type supplied: " + std::to_string(type) + " is not compatible with WS_FT");
    }


    uint64_t count,count2;

    //read operations
    wsfile.read((char *) &count, sizeof(count));
    journal.resize(count);
    for (auto i=0; i < count; i++) {
        auto &j = journal[i];
        sdglib::read_string(wsfile,j.name);
        sdglib::read_string(wsfile,j.detail);
        sdglib::read_string(wsfile,j.tool);
        wsfile.read((char *) &j.timestamp, sizeof(j.timestamp));

        wsfile.read((char *) &count2, sizeof(count2));
        j.entries.resize(count2);
        for (auto e = 0; e < count2; ++e){
            sdglib::read_string(wsfile,j.entries[e].detail);
        }
    }


    if (log_only) return;

    //graph
    sdg.read(wsfile);
    sdglib::OutputLog() <<"Loaded graph with "<<sdg.nodes.size()-1<<" nodes" <<std::endl;

    //distance graphs
    wsfile.read((char *) &count,sizeof(count));
    //distance_graphs.reserve(count);
    for (auto i=0;i<count;++i) {
        distance_graphs.emplace_back(sdg);
        distance_graphs.back().read(wsfile);
    }

    //paired_reads_datastores
    wsfile.read((char *) &count,sizeof(count));
    //paired_reads_datastores.reserve(count);
    for (auto i=0;i<count;++i) {
        paired_reads_datastores.emplace_back(*this);
        paired_reads_datastores.back().read(wsfile);
    }

    //linked_reads_datastores
    wsfile.read((char *) &count,sizeof(count));
    //linked_reads_datastores.reserve(count);
    for (auto i=0;i<count;++i) {
        linked_reads_datastores.emplace_back(*this);
        linked_reads_datastores.back().read(wsfile);
    }

    //Long reads datastores
    wsfile.read((char *) &count,sizeof(count));
    //long_reads_datastores.reserve(count);
    for (auto i=0;i<count;++i) {
        long_reads_datastores.emplace_back(*this);
        long_reads_datastores.back().read(wsfile);
    }

    // Kmer counts datastore
    wsfile.read((char *) &count,sizeof(count));
    for (auto i = 0; i < count; i++) {
        kmer_counters.emplace_back(*this,wsfile);
    }
    sdglib::OutputLog() <<"WS loaded" <<std::endl;
}

std::string WorkSpace::ls(int level,bool recursive) const {
    std::stringstream ss;
    std::string spacer(2*level,' ');
    ss<<spacer<<"SDG Workspace"<<std::endl;
    if (recursive) ss<<sdg.ls(level+1,true);
    ss<<spacer<<"  Distance Graphs: "<<distance_graphs.size()<<std::endl;
    if (recursive){
        for (auto &d:distance_graphs) ss<<d.ls(level+1,true);
    }
    ss<<spacer<<"  Paired Read Datastores: "<<paired_reads_datastores.size()<<std::endl;
    if (recursive){
        for (auto &d:paired_reads_datastores) ss<<d.ls(level+2,true);
    }
    ss<<spacer<<"  Linked Read Datastores: "<<linked_reads_datastores.size()<<std::endl;
    if (recursive){
        for (auto &d:linked_reads_datastores) ss<<d.ls(level+2,true);
    }
    ss<<spacer<<"  Long Read Datastores: "<<long_reads_datastores.size()<<std::endl;
    if (recursive){
        for (auto &d:long_reads_datastores) ss<<d.ls(level+2,true);
    }
    ss<<spacer<<"  Kmer Count Datastores: "<<paired_reads_datastores.size()<<std::endl;
    if (recursive){
        for (auto &d:kmer_counters) ss<<d.ls(level+2,true);
    }
    return ss.str();
}

void WorkSpace::status() {
    for (const auto &j:journal){
        j.status();
    }
    //graph
    sdg.print_status();

    //PR datastores and mappings
    sdglib::OutputLog()<<"Workspace contains "<< paired_reads_datastores.size() << " paired reads datastores" <<std::endl;
    for (auto di=0;di<paired_reads_datastores.size();++di){
        paired_reads_datastores[di].print_status();

    }
    //10x datastores and mappings
    sdglib::OutputLog()<<"Workspace contains "<< linked_reads_datastores.size() << " linked reads datastores" <<std::endl;
    for (auto di=0;di<linked_reads_datastores.size();++di){
        linked_reads_datastores[di].print_status();
    }
    //LR datastores and mappings
    sdglib::OutputLog()<<"Workspace contains "<< long_reads_datastores.size() << " long reads datastores" <<std::endl;
    for (auto di=0;di<long_reads_datastores.size();++di){
        long_reads_datastores[di].print_status();
    }

}

std::vector<sgNodeID_t> WorkSpace::select_from_all_nodes(uint32_t min_size, uint32_t max_size, uint32_t min_tags, uint32_t max_tags,
                                         float min_ci, float max_ci) {
    std::vector<sgNodeID_t> nodes;
    sdglib::OutputLog()<<"Selecting nodes: " << min_size << "-" << max_size << " bp " <<
                      min_tags << "-" << max_tags << " tags " << min_ci << "-" << max_ci << " CI"<<std::endl;
    uint64_t tnodes=0,tbp=0;
    nodes.reserve(sdg.nodes.size());
#pragma omp parallel
    {
        std::vector<sgNodeID_t> thread_nodes;
        uint64_t ttbp=0;
#pragma omp for schedule(static, 100)
        for (auto n=1;n<sdg.nodes.size();++n) {
            if (sdg.nodes[n].sequence.size() < min_size) continue;
            if (sdg.nodes[n].sequence.size() > max_size) continue;
            if (!linked_reads_datastores.empty()) {
                auto ntags = linked_reads_datastores[0].mapper.get_node_tags(n);
                if (ntags.size() < min_tags or ntags.size() > max_tags) continue;
            }
            ++tnodes;
            thread_nodes.emplace_back(n);

            ttbp += sdg.nodes[n].sequence.size();
        }

#pragma omp critical(collect_selected_nodes)
        {
            nodes.insert(nodes.end(), thread_nodes.begin(), thread_nodes.end());
            tbp+=ttbp;
        }

    }
    sdglib::OutputLog()<< "Selected "<<tnodes<<" / "<<sdg.nodes.size()<<" with a total "<<tbp<<"bp"<< std::endl;
    return nodes;
}

PairedReadsDatastore &WorkSpace::add_paired_reads_datastore(const std::string &filename, const std::string &name) {
    if (paired_reads_datastores.size() > MAX_WORKSPACE_VECTOR_SIZE) {
        throw std::runtime_error("Maximum items exceeded, please increase MAX_WORKSPACE_VECTOR_SIZE compile option to add more items");
    }
    for (auto n:list_paired_reads_datastores()) if (n==name) throw std::runtime_error("Name is already in use");
    paired_reads_datastores.emplace_back(*this, filename);
    if (!name.empty()) paired_reads_datastores.back().name = name;
    return paired_reads_datastores.back();
}

LinkedReadsDatastore &WorkSpace::add_linked_reads_datastore(const std::string &filename, const std::string &name) {
    if (linked_reads_datastores.size() > MAX_WORKSPACE_VECTOR_SIZE) {
        throw std::runtime_error("Maximum items exceeded, please increase MAX_WORKSPACE_VECTOR_SIZE compile option to add more items");
    }
    for (auto n:list_linked_reads_datastores()) if (n==name) throw std::runtime_error("Name is already in use");
    linked_reads_datastores.emplace_back(*this, filename);
    if (!name.empty()) linked_reads_datastores.back().name = name;
    return linked_reads_datastores.back();
}

LongReadsDatastore &WorkSpace::add_long_reads_datastore(const std::string &filename, const std::string &name) {
    if (long_reads_datastores.size() > MAX_WORKSPACE_VECTOR_SIZE) {
        throw std::runtime_error("Maximum items exceeded, please increase MAX_WORKSPACE_VECTOR_SIZE compile option to add more items");
    }
    for (auto n:list_long_reads_datastores()) if (n==name) throw std::runtime_error("Name is already in use");
    long_reads_datastores.emplace_back(*this, filename);
    if (!name.empty()) long_reads_datastores.back().name = name;
    return long_reads_datastores.back();
}

DistanceGraph &WorkSpace::add_distance_graph(const DistanceGraph &dg, const std::string &name) {
    if (distance_graphs.size() > MAX_WORKSPACE_VECTOR_SIZE) {
        throw std::runtime_error("Maximum items exceeded, please increase MAX_WORKSPACE_VECTOR_SIZE compile option to add more items");
    }
    for (auto n:list_distance_graphs()) if (n==name) throw std::runtime_error("Name is already in use");
    distance_graphs.emplace_back(dg);
    if (!name.empty()) distance_graphs.back().name = name;
    return distance_graphs.back();
}

JournalOperation &WorkSpace::add_operation(const std::string &name, const std::string &tool, const std::string &detail) {
    journal.emplace_back(name, tool, detail);
    return journal.back();
}

KmerCounter &WorkSpace::add_kmer_counter(const std::string &filename, const std::string &name) {
    if (kmer_counters.size() > MAX_WORKSPACE_VECTOR_SIZE) {
        throw std::runtime_error("Maximum items exceeded, please increase MAX_WORKSPACE_VECTOR_SIZE compile option to add more items");
    }
    for (auto n:list_kmer_counters()) if (n==name) throw std::runtime_error("Name is already in use");
    kmer_counters.emplace_back(*this, filename);
    kmer_counters.back().name=name;
    // Should edit the WS here to include the new count
    return kmer_counters.back();
}

KmerCounter &WorkSpace::add_kmer_counter(const std::string &name, const uint8_t k, KmerCountMode count_mode) {
    if (kmer_counters.size() > MAX_WORKSPACE_VECTOR_SIZE) {
        throw std::runtime_error("Maximum items exceeded, please increase MAX_WORKSPACE_VECTOR_SIZE compile option to add more items");
    }
    for (auto n:list_kmer_counters()) if (n==name) throw std::runtime_error("Name is already in use");
    kmer_counters.emplace_back(*this, name, k, count_mode);
    // Should edit the WS here to include the new count

    return kmer_counters.back();
}

PairedReadsDatastore &WorkSpace::get_paired_reads_datastore(const std::string &name) {
    for (auto &ds : paired_reads_datastores){
        if (ds.name == name) return ds;
    }
    throw std::runtime_error("There are no PairedReadsDatastore named: " + name);
}

LinkedReadsDatastore &WorkSpace::get_linked_reads_datastore(const std::string &name) {
    for (auto &ds : linked_reads_datastores){
        if (ds.name == name) return ds;
    }
    throw std::runtime_error("There are no LinkedReadsDatastore named: " + name);
}

LongReadsDatastore &WorkSpace::get_long_reads_datastore(const std::string &name) {
    for (auto &ds : long_reads_datastores){
        if (ds.name == name) return ds;
    }
    throw std::runtime_error("There are no LongReadsDatastore named: " + name);

}

DistanceGraph &WorkSpace::get_distance_graph(const std::string &name) {
    for (auto &ds : distance_graphs){
        if (ds.name == name) return ds;
    }
    throw std::runtime_error("There are no DistanceGraphs named: " + name);
}

KmerCounter &WorkSpace::get_kmer_counter(const std::string &name) {
    for (auto &ds : kmer_counters) {
        if (ds.name == name) return ds;
    }
    throw std::runtime_error("Couldn't find a KmerCounter named: " + name);
}

WorkSpace::WorkSpace(const std::string &filename) : sdg(*this) {
    linked_reads_datastores.reserve(MAX_WORKSPACE_VECTOR_SIZE);
    paired_reads_datastores.reserve(MAX_WORKSPACE_VECTOR_SIZE);
    long_reads_datastores.reserve(MAX_WORKSPACE_VECTOR_SIZE);
    kmer_counters.reserve(MAX_WORKSPACE_VECTOR_SIZE);
    if (filename!="") load_from_disk(filename);
}

std::vector<std::string> WorkSpace::list_distance_graphs() {
    std::vector<std::string> names;
    for (const auto &d: distance_graphs) {
        names.emplace_back(d.name);
    }
    return names;
}

std::vector<std::string> WorkSpace::list_paired_reads_datastores() {
    std::vector<std::string> names;
    for (const auto &p: paired_reads_datastores) {
        names.emplace_back(p.name);
    }
    return names;
}

std::vector<std::string> WorkSpace::list_long_reads_datastores() {
    std::vector<std::string> names;
    for (const auto &l: long_reads_datastores) {
        names.emplace_back(l.name);
    }
    return names;
}

std::vector<std::string> WorkSpace::list_linked_reads_datastores() {
    std::vector<std::string> names;
    for (const auto &li: linked_reads_datastores) {
        names.emplace_back(li.name);
    }
    return names;
}

std::ostream &operator<<(std::ostream &os, const WorkSpace &ws) {
    os << "WorkSpace @ " << &ws;
    return os;
}

std::vector<std::string> WorkSpace::list_kmer_counters() {
    std::vector<std::string> names;
    for (const auto &kc: kmer_counters) {
        names.emplace_back(kc.name);
    }
    return names;
}

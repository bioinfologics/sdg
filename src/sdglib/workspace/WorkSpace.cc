//
// Created by Bernardo Clavijo (EI) on 24/02/2018.
//

#include "WorkSpace.hpp"


const bsgVersion_t WorkSpace::min_compat = 0x0001;

void WorkSpace::add_log_entry(std::string text) {
    log.emplace_back(std::time(0),std::string(GIT_COMMIT_HASH),text);
}

void WorkSpace::dump_to_disk(std::string filename) {
    sdglib::OutputLog()<<"Dumping workspace to "<<filename<<std::endl;
    std::ofstream of(filename);
    //dump log
    uint64_t count;
    count=log.size();
    of.write((char *) &BSG_MAGIC, sizeof(BSG_MAGIC));
    of.write((char *) &BSG_VN, sizeof(BSG_VN));
    BSG_FILETYPE type(WS_FT);
    of.write((char *) &type, sizeof(type));
    of.write((char *) &count,sizeof(count));
    for (auto &l:log){
        of.write((char *) &l.timestamp,sizeof(l.timestamp));
        count=l.bsg_version.size();
        of.write((char *) &count,sizeof(count));
        of.write((char *) l.bsg_version.data(), count);
        count=l.log_text.size();
        of.write((char *) &count,sizeof(count));
        of.write((char *) l.log_text.data(), count);

    }
    //dump main graph
    sdg.write(of);
    //dump KCI
    kci.write(of);

    //paired read datastores
    count=paired_read_datastores.size();
    of.write((char *) &count,sizeof(count));
    for (auto i=0;i<count;++i){
        paired_read_datastores[i].write(of);
        paired_read_datastores[i].mapper.write(of);
    }

    //linker read datastores
    count=linked_read_datastores.size();
    of.write((char *) &count,sizeof(count));
    for (auto i=0;i<count;++i){
        linked_read_datastores[i].write(of);
        linked_read_datastores[i].mapper.write(of);
    }

    //long read datastores
    count=long_read_datastores.size();
    of.write((char *) &count,sizeof(count));
    for (auto i=0;i<count;++i){
        long_read_datastores[i].write(of);
        long_read_datastores[i].mapper.write(of);
    }
    //dump element type then use that element's own dump to dump it to this file

//    //[GONZA]
//    int v = (int) read_counts_header.size();
//    of.write((const char *) &v, sizeof(v));
//    for (int i=0; i<read_counts_header.size(); ++i){
//        int hs = (int) read_counts_header[i].size();
//        of.write((const char *) &hs, sizeof(hs));
//        of.write((const char *) read_counts_header[i].data(), hs*sizeof(char));
//    }
}

void WorkSpace::load_from_disk(std::string filename, bool log_only) {
    std::ifstream wsfile(filename);
    if (!wsfile.good()) {
        std::cerr << filename << " opening error: " << strerror(errno) << std::endl;
        throw std::runtime_error("Error opening " + filename);
    }
    uint64_t count;
    //read log
    log.clear();
    std::string git_version,text;
    std::time_t timestamp;

    bsgMagic_t magic;
    bsgVersion_t version;
    BSG_FILETYPE type;
    wsfile.read((char *) &magic, sizeof(magic));
    wsfile.read((char *) &version, sizeof(version));
    wsfile.read((char *) &type, sizeof(type));

    if (magic != BSG_MAGIC) {
        throw std::runtime_error("Magic number not present in the kci file");
    }

    if (version < min_compat) {
        throw std::runtime_error("kci file version: " + std::to_string(version) + " is not compatible with " + std::to_string(min_compat));
    }

    if (type != WS_FT) {
        throw std::runtime_error("File type supplied: " + std::to_string(type) + " is not compatible with WS_FT");
    }


    wsfile.read((char *) &count,sizeof(count));
    for (auto i=0;i<count;++i){
        uint64_t ssize;
        wsfile.read((char *) &timestamp,sizeof(timestamp));
        wsfile.read((char *) &ssize,sizeof(ssize));
        git_version.resize(ssize);
        wsfile.read((char *) git_version.data(),ssize);
        wsfile.read((char *) &ssize,sizeof(ssize));
        text.resize(ssize);
        wsfile.read((char *) text.data(),ssize);
        log.emplace_back(timestamp,git_version,text);
    }
    if (log_only) return;
    //read element type, then use that element's read
    sdg.read(wsfile);
    sdglib::OutputLog() <<"Loaded graph with "<<sdg.nodes.size()-1<<" nodes" <<std::endl;
    kci.read(wsfile);
    //todo: read all datastores and mappers!!!
//    if (format>1) {
//        wsfile.read((char *) &count,sizeof(count));
//        for (auto i=0;i<count;++i){
//            paired_read_datastores.emplace_back();
//            paired_read_datastores.back().read(wsfile);
//            paired_read_mappers.emplace_back(sdg,linked_read_datastores.back());
//            paired_read_mappers.back().read(wsfile);
//        }
//    }
    wsfile.read((char *) &count,sizeof(count));
    paired_read_datastores.reserve(count);
    for (auto i=0;i<count;++i) {
        paired_read_datastores.emplace_back(*this, wsfile);
    }
    wsfile.read((char *) &count,sizeof(count));
    linked_read_datastores.reserve(count);

    for (auto i=0;i<count;++i) {
        linked_read_datastores.emplace_back(*this, wsfile);
    }

    wsfile.read((char *) &count,sizeof(count));
    long_read_datastores.reserve(count);
    for (auto i=0;i<count;++i) {
        long_read_datastores.emplace_back(*this, wsfile);
    }

}

void WorkSpace::print_log() {
    for (auto l:log){
        char buf[50];
        std::strftime(buf,50,"%F %R",std::localtime(&l.timestamp));
        std::cout<<buf<<" ["<<l.bsg_version<<"]: "<<l.log_text<<std::endl;
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
            if (!linked_read_datastores.empty()) {
                auto ntags = linked_read_datastores[0].mapper.get_node_tags(n);
                if (ntags.size() < min_tags or ntags.size() > max_tags) continue;
            }
            if (!kci.graph_kmers.empty()) {
                auto ci = kci.compute_compression_for_node(n, 1);
                if (ci < min_ci or ci > max_ci) continue;
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

void WorkSpace::remap_all() {
    sdglib::OutputLog()<<"Mapping reads..."<<std::endl;
    //auto pri=0;
    sdg.create_index();
    for (auto &ds:paired_read_datastores) {
        sdglib::OutputLog()<<"Mapping reads from paired library..."<<std::endl;
        ds.mapper.remap_all_reads();
        ds.print_status();
        sdglib::OutputLog()<<"Computing size distribution..."<<std::endl;
        //auto sdist=m.size_distribution();
        //std::ofstream df("prdist_"+std::to_string(pri++)+".csv");
        //for (auto i=0;i<sdist.size();i+=10){
        //    uint64_t t=0;
        //    for (auto j=i;j<i+10;++j) t+=sdist[j];
        //    if (t>0) df<<i<<", "<<t<<std::endl;
        //}
        add_log_entry("reads from "+ds.filename+" re-mapped to current graph");
        sdglib::OutputLog()<<"Mapping reads from paired library DONE."<<std::endl;
    }
    for (auto &ds:linked_read_datastores) {
        sdglib::OutputLog()<<"Mapping reads from linked library..."<<std::endl;
        ds.mapper.remap_all_reads();
        add_log_entry("reads from "+ds.filename+" re-mapped to current graph");
        sdglib::OutputLog()<<"Mapping reads from linked library DONE."<<std::endl;
    }
}

void WorkSpace::remap_all63() {
    sdglib::OutputLog()<<"Mapping reads..."<<std::endl;
    sdg.create_63mer_index();
    for (auto &ds:paired_read_datastores) {
        sdglib::OutputLog()<<"Mapping reads from paired library..."<<std::endl;
        ds.mapper.remap_all_reads63();
        ds.print_status();
        add_log_entry("reads from "+ds.filename+" re-mapped to current graph");
        sdglib::OutputLog()<<"Mapping reads from paired library DONE."<<std::endl;
    }
    for (auto &ds:linked_read_datastores) {
        sdglib::OutputLog()<<"Mapping reads from linked library..."<<std::endl;
        ds.mapper.remap_all_reads63();
        add_log_entry("reads from "+ds.filename+" re-mapped to current graph");
        sdglib::OutputLog()<<"Mapping reads from linked library DONE."<<std::endl;
    }
}


//
// Created by Bernardo Clavijo (EI) on 24/02/2018.
//

#include "WorkSpace.hpp"

void WorkSpace::add_log_entry(std::string text) {
    log.emplace_back(std::time(0),std::string(GIT_COMMIT_HASH),text);
}

void WorkSpace::dump_to_disk(std::string filename) {
    sglib::OutputLog()<<"Dumping workspace to "<<filename<<std::endl;
    std::ofstream of(filename);
    //dump log
    uint64_t count;
    count=log.size();
    of.write((char *) &count,sizeof(count));
    for (auto &l:log){
        of.write((char *) &l.timestamp,sizeof(l.timestamp));
        count=l.bsg_version.size();
        of.write((char *) &count,sizeof(count));
        of.write((char *) l.bsg_version.c_str(), count);
        count=l.log_text.size();
        of.write((char *) &count,sizeof(count));
        of.write((char *) l.log_text.c_str(), count);

    }
    //dump main graph
    sg.write(of);
    //dump KCI
    kci.write(of);

    //paired read datastores
    count=paired_read_datastores.size();
    of.write((char *) &count,sizeof(count));
    for (auto i=0;i<count;++i){
        paired_read_datastores[i].write(of);
        paired_read_mappers[i].write(of);
    }

    //linker read datastores
    count=linked_read_datastores.size();
    of.write((char *) &count,sizeof(count));
    for (auto i=0;i<count;++i){
        linked_read_datastores[i].write(of);
        linked_read_mappers[i].write(of);
    }


    count=path_datastores.size();
    of.write((char *) &count,sizeof(count));
    for (auto i=0;i<count;++i){
        path_datastores[i].write(of);
    }
    //dump element type then use that element's own dump to dump it to this file
}

void WorkSpace::load_from_disk(std::string filename, bool log_only) {
    std::ifstream wsfile(filename);
    uint64_t count;
    //read log
    log.clear();
    std::string version,text;
    std::time_t timestamp;
    wsfile.read((char *) &count,sizeof(count));
    for (auto i=0;i<count;++i){
        uint64_t ssize;
        wsfile.read((char *) &timestamp,sizeof(timestamp));
        wsfile.read((char *) &ssize,sizeof(ssize));
        version.resize(ssize);
        wsfile.read((char *) version.c_str(),ssize);
        wsfile.read((char *) &ssize,sizeof(ssize));
        text.resize(ssize);
        wsfile.read((char *) text.c_str(),ssize);
        log.emplace_back(timestamp,version,text);
    }
    if (log_only) return;
    //read element type, then use that element's read
    sg.read(wsfile);
    std::cout<<"Loaded graph with "<<sg.nodes.size()-1<<" nodes" <<std::endl;
    kci.read(wsfile);
    //todo: read all datastores and mappers!!!
//    if (format>1) {
//        wsfile.read((char *) &count,sizeof(count));
//        for (auto i=0;i<count;++i){
//            paired_read_datastores.emplace_back();
//            paired_read_datastores.back().read(wsfile);
//            paired_read_mappers.emplace_back(sg,linked_read_datastores.back());
//            paired_read_mappers.back().read(wsfile);
//        }
//    }
    wsfile.read((char *) &count,sizeof(count));
    for (auto i=0;i<count;++i){
        paired_read_datastores.emplace_back();
        paired_read_datastores.back().read(wsfile);
        paired_read_mappers.emplace_back(sg,paired_read_datastores.back());
        paired_read_mappers.back().read(wsfile);
    }
    wsfile.read((char *) &count,sizeof(count));
    for (auto i=0;i<count;++i){
        linked_read_datastores.emplace_back();
        linked_read_datastores.back().read(wsfile);
        linked_read_mappers.emplace_back(sg,linked_read_datastores.back());
        linked_read_mappers.back().read(wsfile);
    }
    if (!wsfile.eof()){
        wsfile.read((char *) &count,sizeof(count));
        for (auto i=0;i<count;++i){
            path_datastores.emplace_back(sg);
            path_datastores.back().read(wsfile);
        }
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
    sglib::OutputLog()<<"Selecting nodes: " << min_size << "-" << max_size << " bp " <<
                      min_tags << "-" << max_tags << " tags " << min_ci << "-" << max_ci << " CI"<<std::endl;
    uint64_t tnodes=0,tbp=0;
    nodes.reserve(sg.nodes.size());
#pragma omp parallel
    {
        std::vector<sgNodeID_t> thread_nodes;
#pragma omp for schedule(static, 100)
        for (auto n=1;n<sg.nodes.size();++n) {
            if (sg.nodes[n].sequence.size() < min_size) continue;
            if (sg.nodes[n].sequence.size() > max_size) continue;
            if (sg.get_fw_links(n).size() != 1 or sg.get_bw_links(n).size() != 1) continue;
            auto ntags = linked_read_mappers[0].get_node_tags(n);
            if (ntags.size() < min_tags or ntags.size() > max_tags) continue;
            auto ci = kci.compute_compression_for_node(n, 1);
            if (ci < min_ci or ci > max_ci) continue;
            ++tnodes;
            thread_nodes.emplace_back(n);
            tbp += sg.nodes[n].sequence.size();
        }

#pragma omp critical(collect_selected_nodes)
        nodes.insert(nodes.end(),thread_nodes.begin(),thread_nodes.end());

    }
    sglib::OutputLog()<< "Selected "<<tnodes<<" / "<<sg.nodes.size()<<" with a total "<<tbp<<"bp"<< std::endl;
    return nodes;
}
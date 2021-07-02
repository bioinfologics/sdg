//
// Created by Bernardo Clavijo (EI) on 30/06/2021.
//

#include "TotalSorter.hpp"

TotalSorter::TotalSorter(const ReadThreadsGraph &rtg, int min_thread_length, int min_node_threads) : rtg(rtg){
    //TODO: maybe copy the rtg and pop nodes/remove threads
    //init nodes with all nodes having enough threads of any kind
    for (auto nv:rtg.get_all_nodeviews(false,false)){
        if (rtg.node_threads(nv.node_id()).size()>=min_node_threads) nodes.insert(nv.node_id());
    }
    //init threads with all threads with enough nodes of any kind
    for (auto &ti:rtg.thread_info){
        if (ti.second.link_count+1>=min_thread_length) threads.insert(ti.first);
    }
    sdglib::OutputLog()<<"TotalSorter starting selection: "<<nodes.size()<<" nodes and "<<threads.size()<<" threads"<<std::endl;
    bool changed=true;
    while (changed) {
        changed=false;
        std::unordered_set<sgNodeID_t> new_nodes;
        for (auto nid:nodes){
            int valid_threads=0;
            for (auto tid:rtg.node_threads(nid)) if (threads.count(tid)) ++valid_threads;
            if (valid_threads>=min_node_threads) new_nodes.insert(nid);
            else changed=true;
        }
        std::swap(nodes,new_nodes);
        std::unordered_set<int64_t> new_threads;
        for (auto tid:threads){
            int valid_nodes=0;
            for (auto tn:rtg.get_thread(tid)) if (nodes.count(llabs(tn.node))) ++valid_nodes;
            if (valid_nodes>=min_thread_length) new_threads.insert(tid);
            else changed=true;
        }
        std::swap(threads,new_threads);
        if (changed)
            sdglib::OutputLog()<<"TotalSorter adjusted selection: "<<nodes.size()<<" nodes and "<<threads.size()<<" threads"<<std::endl;
    }
    sdglib::OutputLog()<<"TotalSorter final selection: "<<nodes.size()<<" nodes and "<<threads.size()<<" threads"<<std::endl;
}

TotalSorter::TotalSorter(const ReadThreadsGraph &rtg, std::string filename) : rtg(rtg){
    load(filename);
}

void TotalSorter::run_sorters_from_lines(int min_line_size) {
    //XXX: hardcoded parameters so far
    float line_occupancy=.7;
    int p=4,q=6;
    float node_happiness=.6;
    int node_threads=5;
    int end_size=30;
    int min_links=3;
    //END OF HARDCODE
    auto lines=rtg.sdg.get_all_lines(min_line_size);
    int finished_orders=sorters.size();
#pragma omp parallel for schedule(dynamic,1) shared(finished_orders)
    for (auto lineit=lines.cbegin();lineit<lines.cend();++lineit){
        auto &line=*lineit;
        int used_nodes=0;
        for (auto nid:line) if (node_sorters.count(llabs(nid)) and not node_sorters[llabs(nid)].empty()) ++used_nodes;
        if (used_nodes>=line_occupancy*line.size()) continue;
        try {
            int64_t si;
#pragma omp critical(orders)
            {
                si = next_sorter++;
                sorters.emplace(std::piecewise_construct,std::forward_as_tuple(si),std::forward_as_tuple(rtg, node_happiness, node_threads, end_size));
            }
            sorters.at(si).start_from_nodelist(line,min_links);
            sorters.at(si).recruit_all_happy_threads_q(p,q);
            while (sorters.at(si).grow(4,end_size,3,node_happiness)){
                sorters.at(si).recruit_all_happy_threads_q(p,q);
            }
#pragma omp critical(orders)
            {
                ++finished_orders;
                for (auto nid:sorters.at(si).order.as_signed_nodes()) if (nodes.count(llabs(nid))) node_sorters[llabs(nid)].emplace_back(nid>0?si:-si);
                for (auto tid:sorters.at(si).threads) if (threads.count(llabs(tid))) thread_sorters[llabs(tid)].emplace_back(tid>0?si:-si);

                int64_t nodes_in_sorter=0;
                for (auto &ns:node_sorters) if (not ns.second.empty()) ++nodes_in_sorter;
                int64_t threads_in_sorter=0;
                for (auto &ns:thread_sorters) if (not ns.second.empty()) ++threads_in_sorter;
                sdglib::OutputLog()<<finished_orders<<" sorters run, "<<nodes_in_sorter<<" / "<<nodes.size()<<" nodes and "<<threads_in_sorter<<" / "<<threads.size()<<" threads used"<<std::endl;

            }

        }
        catch (std::exception &e) {
            std::cout<<"Exception caught when trying to process line with node "<<line[0]<<": "<< typeid(e).name() <<": "<<e.what()<<", skipping"<<std::endl;
        }

    }
    sdglib::OutputLog()<<"TotalSorter run finished!!!"<<std::endl;

}

void TotalSorter::compute_sorter_classes() {

}

void TotalSorter::dump(std::string filename) {
    std::ofstream ofs(filename);
    ofs.write(reinterpret_cast<const char *>(&next_sorter),sizeof(next_sorter));
    sdglib::write_flat_unorderedset(ofs,nodes);
    sdglib::write_flat_unorderedset(ofs,threads);
    int64_t count=sorters.size();
    ofs.write(reinterpret_cast<const char *>(&count),sizeof(count));
    for (auto &s:sorters) {
        ofs.write(reinterpret_cast<const char *>(&s.first),sizeof(s.first));
        s.second.write(ofs);
    }
    sdglib::write_flat_unorderedmap(ofs,sorter_classes);
    count=node_sorters.size();
    ofs.write(reinterpret_cast<const char *>(&count),sizeof(count));
    for (auto &s:node_sorters) {
        ofs.write(reinterpret_cast<const char *>(&s.first),sizeof(s.first));
        sdglib::write_flat_vector(ofs,s.second);
    }
    count=thread_sorters.size();
    ofs.write(reinterpret_cast<const char *>(&count),sizeof(count));
    for (auto &s:thread_sorters) {
        ofs.write(reinterpret_cast<const char *>(&s.first),sizeof(s.first));
        sdglib::write_flat_vector(ofs,s.second);
    }
}

void TotalSorter::load(std::string filename) {
    std::ifstream ifs(filename);
    ifs.read(reinterpret_cast<char *>(&next_sorter),sizeof(next_sorter));
    sdglib::read_flat_unorderedset(ifs,nodes);
    sdglib::read_flat_unorderedset(ifs,threads);
    int64_t count,ix;
    ifs.read(reinterpret_cast<char *>(&count),sizeof(count));
    sorters.clear();
    sorters.reserve(count);
    for (auto i=count;i>0;--i) {
        ifs.read(reinterpret_cast<char *>(&ix),sizeof(ix));
        sorters.emplace(std::piecewise_construct,std::forward_as_tuple(ix),std::forward_as_tuple(rtg));
        sorters.at(ix).read(ifs);
    }
    sdglib::read_flat_unorderedmap(ifs,sorter_classes);
    ifs.read(reinterpret_cast<char *>(&count),sizeof(count));
    node_sorters.clear();
    node_sorters.reserve(count);
    for (auto i=count;i>0;--i) {
        ifs.read(reinterpret_cast<char *>(&ix),sizeof(ix));
        sdglib::read_flat_vector(ifs,node_sorters[ix]);
    }
    ifs.read(reinterpret_cast<char *>(&count),sizeof(count));
    thread_sorters.clear();
    thread_sorters.reserve(count);
    for (auto i=count;i>0;--i) {
        ifs.read(reinterpret_cast<char *>(&ix),sizeof(ix));
        sdglib::read_flat_vector(ifs,thread_sorters[ix]);
    }

}
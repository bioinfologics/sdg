//
// Created by Bernardo Clavijo (EI) on 30/06/2021.
//

#include "TotalSorter.hpp"

TotalSorter::TotalSorter(ReadThreadsGraph &rtg, int min_thread_length, int min_node_threads) : rtg(rtg){
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

TotalSorter::TotalSorter(ReadThreadsGraph &rtg, std::string filename) : rtg(rtg){
    load(filename);
}

void TotalSorter::prune_rtg() {
    sdglib::OutputLog()<<"Pruning unselected threads"<<std::endl;
    {
        std::vector<int64_t> to_delete;
        for (auto &tinf:rtg.thread_info) if (tinf.second.start!=0 and threads.count(tinf.first) == 0) to_delete.emplace_back(tinf.first);
        for (auto tid:to_delete) rtg.remove_thread(tid);
    }
    sdglib::OutputLog()<<"Pruning unselected nodes"<<std::endl;
    {
        std::vector<sgNodeID_t> to_delete;
        for (auto &nv:rtg.get_all_nodeviews(false, false)) if (nodes.count(nv.node_id()) == 0) to_delete.emplace_back(nv.node_id());
        for (auto &nid:to_delete) rtg.pop_node_from_all(nid);
    }
    sdglib::OutputLog()<<"Pruning done"<<std::endl;
}

void TotalSorter::run_sorters_from_lines(std::vector<std::vector<sgNodeID_t>> _lines, int min_line_size, int min_order_size,float line_occupancy,int p, int q, float node_happiness, int node_threads, int end_size, int min_links) {

    std::vector<std::vector<sgNodeID_t>> lines;
    if (not _lines.empty()) lines=_lines;
    else lines=rtg.sdg.get_all_lines(min_line_size);
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
            LocalOrder old_order=sorters.at(si).order;
            while (sorters.at(si).grow(4,end_size,3,node_happiness)){
                if (sorters.at(si).is_mixed()) {
                    sorters.at(si).order=old_order;
                    std::cout<<"Order was mixed, backtracking to last non-mixed order of size "<<old_order.size()<<std::endl;
                    break;
                }
                sorters.at(si).recruit_all_happy_threads_q(p,q);
                old_order=sorters.at(si).order;
            }
#pragma omp critical(orders)
            {
                if (sorters.at(si).order.size()<min_order_size or sorters.at(si).is_mixed()) sorters.erase(si);
                else {
                    ++finished_orders;
                    for (auto nid:sorters.at(si).order.as_signed_nodes())
                        if (nodes.count(llabs(nid)))
                            node_sorters[llabs(nid)].emplace_back(nid > 0 ? si : -si);
                    for (auto tid:sorters.at(si).threads)
                        if (threads.count(llabs(tid)))
                            thread_sorters[llabs(tid)].emplace_back(tid > 0 ? si : -si);

                    int64_t nodes_in_sorter = 0;
                    for (auto &ns:node_sorters) if (not ns.second.empty()) ++nodes_in_sorter;
                    int64_t threads_in_sorter = 0;
                    for (auto &ns:thread_sorters) if (not ns.second.empty()) ++threads_in_sorter;
                    sdglib::OutputLog() << finished_orders << " sorters run, " << nodes_in_sorter << " / "
                                        << nodes.size() << " nodes and " << threads_in_sorter << " / " << threads.size()
                                        << " threads used" << std::endl;
                }
            }

        }
        catch (std::exception &e) {
            std::cout<<"Exception caught when trying to process line with node "<<line[0]<<": "<< typeid(e).name() <<": "<<e.what()<<", skipping"<<std::endl;
        }

    }
    sdglib::OutputLog()<<"TotalSorter run finished!!!"<<std::endl;

}

void TotalSorter::update_usage() {
    node_sorters.clear();
    thread_sorters.clear();
    for (auto &si:sorters) {
        for (auto nid:si.second.order.as_signed_nodes())
            if (nodes.count(llabs(nid)))
                node_sorters[llabs(nid)].emplace_back(nid > 0 ? si.first : -si.first);
        for (auto tid:si.second.threads)
            if (threads.count(llabs(tid)))
                thread_sorters[llabs(tid)].emplace_back(tid > 0 ? si.first : -si.first);
    }
}

void TotalSorter::remove_mixed(int win, float fail) {
    std::vector<int64_t> to_delete;
    for (auto &s:sorters){
        if (s.second.is_mixed(win,fail)) to_delete.push_back(s.first);
    }
    for (auto &si:to_delete) sorters.erase(si);
}

void TotalSorter::compute_node_neighbours(int k, int max_f) {
    NKmerIndex nki(rtg.sdg,k,max_f);
    std::unordered_map<std::pair<sgNodeID_t,sgNodeID_t>,uint64_t> shared_counts;
    uint64_t last_kmer=-1;
    std::set<sgNodeID_t> kmer_nodes;
    for (auto &ke:nki){
        if (ke.kmer==last_kmer) kmer_nodes.insert(labs(ke.contigID));
        else {
            if (kmer_nodes.size()>1){
                for (auto n1:kmer_nodes)
                    for (auto &n2:kmer_nodes)
                        if (n1<n2) ++shared_counts[std::make_pair(n1,n2)];
            }
            last_kmer=ke.kmer;
            kmer_nodes.clear();
        }
    }
    for (auto &sc:shared_counts) {
        if ( rtg.sdg.are_connected(sc.first.first,sc.first.second) or
            rtg.sdg.are_connected(sc.first.first,-sc.first.second) or
            rtg.sdg.are_connected(-sc.first.first,sc.first.second) or
            rtg.sdg.are_connected(-sc.first.first,-sc.first.second) ) continue;
        node_neighbours[sc.first.first].emplace_back(sc.first.second,sc.second);
        node_neighbours[sc.first.second].emplace_back(sc.first.first,sc.second);
    }
}

void TotalSorter::compute_sorter_classes() {

}

int64_t TotalSorter::sorter_shared_nodes(int64_t oi1, int64_t oi2) {
    std::unordered_set<sgNodeID_t> n1;
    int64_t c=0;
    if (oi1>0) for (auto &nid:sorters.at(oi1).order.as_signed_nodes()) n1.insert(nid);
    else for (auto &nid:sorters.at(-oi1).order.as_signed_nodes()) n1.insert(-nid);
    if (oi2>0) {
        for (auto &nid:sorters.at(oi2).order.as_signed_nodes()) if (n1.count(nid)) ++c;
    }
    else {
        for (auto &nid:sorters.at(-oi2).order.as_signed_nodes()) if (n1.count(-nid)) ++c;
    }
    return c;
}

void TotalSorter::merge() {
    std::set<std::pair<int64_t,int64_t>> failed_merges;
    size_t last_sorters_count=sorters.size()+1;
    size_t last_fails_count=0;

    int64_t ls1,ls2;
    while (sorters.size()<last_sorters_count or failed_merges.size()>last_fails_count) {
        last_sorters_count=sorters.size();
        last_fails_count=failed_merges.size();
        std::cout<<"There are now "<<last_sorters_count<<" sorters, and "<<last_fails_count<<" tabu merges"<<std::endl;
        int64_t largest_shared=0;
        float largest_shared_percentage=0;
        //Finds the merge with the biggest shared nodes that has not been tried before, at least 100 shared nodes
        for (auto &ss1:sorters) {
            auto s1=ss1.first;
            auto size1=ss1.second.order.size();
            //if (ss1.second.order.size()<largest_shared) continue;
            for (auto &ss2:sorters) {
                auto size2=ss2.second.order.size();
                //if (ss2.second.order.size()<largest_shared) continue;
                auto s2=ss2.first;
                if (s1>=s2) continue;
                if (failed_merges.count({s1,s2})>0) continue;
                for (std::pair<int64_t,int64_t> sp : {std::make_pair(s1,s2),std::make_pair(s1,-s2)}) {
                    auto shared=sorter_shared_nodes(sp.first,sp.second);
                    float shared_percentage=((float) shared) / std::min(size1,size2);
                    //std::cout<<ss1.first<<" & "<<ss2.first<<": "<<shared<<" ("<<shared_percentage<<")"<<std::endl;
                    if (shared>100 and shared_percentage>largest_shared_percentage) {
                        largest_shared_percentage=shared_percentage;
                        largest_shared=shared;
                        ls1=sp.first;
                        ls2=sp.second;
                    }
                }
            }
        }
        if (largest_shared>100)
        {
            std::cout << "Merging orders " << ls1 << " and " << ls2 << ", which share " << largest_shared << " nodes ("<<(int)(100*largest_shared_percentage)<<"%)" << std::endl;
            auto ix=next_sorter++;
            std::cout << "New order to have id "<<ix << std::endl;
            sorters.emplace(std::piecewise_construct,std::forward_as_tuple(ix),std::forward_as_tuple(rtg));
            //copying the first order, which is always positive
            auto &s=sorters.at(ix);
            s.order=sorters.at(ls1).order;
            std::vector<sgNodeID_t> s2nodes=sorters.at(llabs(ls2)).order.as_signed_nodes();
            if (ls2<0) {
                std::vector<sgNodeID_t> rnodes;
                for (auto it=s2nodes.rbegin();it!=s2nodes.rend();++it) rnodes.emplace_back(-*it);
                std::swap(s2nodes,rnodes);
            }
            //now adding nodes from the other order one by one from the 10th shared up
            uint64_t shared=0;
            uint64_t placed=0,unplaced=0;

            for (auto nid:s2nodes){
                if (s.order.get_node_position(nid)<0) {
                    ++unplaced;
                }
                else if (s.order.get_node_position(nid)>0) ++shared;
                else if (shared>=10) {
                    try {
                        s.order.add_placed_nodes({{nid, s.order.place_node(rtg, nid)}});
                        ++placed;
                    }
                    catch (std::exception &e) {
                        ++unplaced;
                    }
                }
            }
            //std::cout<<"After FW loop: "<<placed<<" placed, and "<<unplaced<<" unplaced nodes"<<std::endl;
            shared=0;
            //now adding ramaining nodes from the other order one by one from the last shared up
            for (int64_t i=s2nodes.size()-1;i>=0;--i){
                auto nid=s2nodes[i];
                if (s.order.get_node_position(nid)<0) {
                    ++unplaced;
                }
                else if (s.order.get_node_position(nid)>0) ++shared;
                else if (shared>=10) {
                    try {
                        s.order.add_placed_nodes({{nid, s.order.place_node(rtg, nid)}});
                        ++placed;
                    }
                    catch (std::exception &e) {
                        ++unplaced;
                    }
                }
            }
            //std::cout<<"After BW loop: "<<placed<<" placed, and "<<unplaced<<" unplaced nodes"<<std::endl;
            //now recruiting threads, and checking if it is mixed and how many of each other order's node it has: 99% and not mixed to consider it a successful mix
            s.recruit_all_happy_threads_q(5,7);
            //if merge unsuccessful, delete order and add pair to failed merges
            if (sorters.at(ix).order.size()<.99*(sorters.at(ls1).order.size()+sorters.at(llabs(ls2)).order.size()-largest_shared) or sorters.at(ix).is_mixed()) {
                std::cout<<"merge failed with "<< sorters.at(ix).order.size()<<" / "<<(sorters.at(ls1).order.size()+sorters.at(llabs(ls2)).order.size()-largest_shared)<<" nodes, mixed: "<<sorters.at(ix).is_mixed()<<std::endl;
                sorters.erase(ix);
                if (sorters.at(ls1).order.size()>sorters.at(llabs(ls2)).order.size()) {
                    if(sorters.at(llabs(ls2)).order.size() * .99 <largest_shared ) {
                        sorters.erase(llabs(ls2));
                        std::cout << "deleting order " << llabs(ls2) << " as 99%+ contained" << std::endl;
                    }
                }
                else if(sorters.at(ls1).order.size() * .99 <largest_shared ){
                    sorters.erase(ls1);
                    std::cout<<"deleting order "<<ls1<<" as 99%+ contained"<<std::endl;
                }
                else {
                    failed_merges.insert({ls1, llabs(ls2)});
                    std::cout<<"adding merge to tabu"<<std::endl;
                }
            }
            else {
                sorters.erase(ls1);
                sorters.erase(llabs(ls2));
                std::cout<<"merge successful"<<std::endl;
            }
        }
    }
    std::cout<<"Merging finished"<<std::endl;
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
//
// Created by Bernardo Clavijo (EI) on 19/01/2018.
//

#include "GraphPartitioner.hpp"

std::vector<SequenceSubGraph> GraphPartitioner::partitions_as_subgraphs(
        const SequenceSubGraph subgraph,
        const std::vector<std::vector<bool>> partitions) {
    std::vector<SequenceSubGraph> subgraphs;
    for (auto &p:partitions){
        SequenceSubGraph ssg(sg);
        for (auto i=0;i<subgraph.nodes.size();++i) if (p[i]) ssg.nodes.push_back(subgraph.nodes[i]);
        subgraphs.push_back(ssg);
    }
    return subgraphs;
}

std::vector<std::vector<bool>> GraphPartitioner::tags_patterns(const SequenceSubGraph subgraph) {
    //Find all local tags
    std::set<uint32_t> tags;
    for (auto node:subgraph.nodes){
        auto n=(node>0?node:-node);
        for (auto r:rmappers[0].reads_in_node[n])
            if (rmappers[0].read_to_tag[r.read_id]>0) tags.insert(rmappers[0].read_to_tag[r.read_id]);
    }

    //create the node/tag read matrix.
    std::map<uint32_t,uint32_t> tags_to_local;
    std::vector<uint32_t>local_to_tag;
    {
        auto i=0;
        for (auto it = tags.begin(); it != tags.end(); ++it, ++i) {
            tags_to_local[*it] = i;
            local_to_tag.push_back(*it);
        }
    }

    std::vector<std::vector<uint64_t>> node_tag_readcount;
    for (auto node:subgraph.nodes) {
        node_tag_readcount.emplace_back(std::vector<uint64_t>(tags.size()));
        auto n = (node > 0 ? node : -node);
        for (auto r:rmappers[0].reads_in_node[n])
            if (rmappers[0].read_to_tag[r.read_id]>0)
                ++node_tag_readcount.back()[tags_to_local[rmappers[0].read_to_tag[r.read_id]]];
    }

    std::vector<std::vector<bool>> tp;
    for (auto xt=0;xt<node_tag_readcount.back().size();++xt){
        auto nodes=0;
        std::vector<bool> newpat;
        for (auto xn=0;xn<node_tag_readcount.size();++xn) {
            newpat.emplace_back((node_tag_readcount[xn][xt]>0));
            if (node_tag_readcount[xn][xt]>0) ++nodes;
        }
        if (nodes>1) {
            bool repeated=false;
            for (auto t:tp) if (newpat==t) {repeated=true;break;}
            if (!repeated) tp.push_back(newpat);
        }
    }

    return tp;
}

std::pair<float,float> GraphPartitioner::score_partition_set(const SequenceSubGraph subgraph,
                                                             const std::vector<std::vector<bool>> partitions,
                                                             const std::vector<std::vector<bool>> tag_patterns) {

    auto tp=tag_patterns;
    if (tag_patterns.size()==0){
        tp=tags_patterns(subgraph);
    }
    if (tp.size()==0) return {0,0};

    std::vector<std::vector<uint32_t>> tag_parts_all(tp.size()); //partitions containing ALL tag
    std::vector<std::vector<uint32_t>> tag_parts_some(tp.size()); //partitions containing SOME of the tag
    //TODO: save the tags into partitions (to use by further heuristics) also, print details
    //std::vector<std::vector<uint32_t>> part_tags(partitions.size());
    //std::vector<std::vector<uint32_t>> part_shared_tags(partitions.size());
    //std::vector<std::vector<uint32_t>> part_failed_tags(partitions.size());

    //go through partitions and tags, finding perfect and partial containment of tags in partitions
    for (auto xt=0;xt<tp.size();++xt){
        for (auto xp=0;xp<partitions.size();++xp){
            bool all_contained=true;
            bool none_contained=true;
            for (auto i=0;i<tp.back().size();++i){
                if (tp[xt][i]){
                    if(partitions[xp][i]) none_contained=false;
                    else all_contained=false;
                }
            }
            if (all_contained) tag_parts_all[xt].push_back(xp);
            else if (!none_contained) tag_parts_some[xt].push_back(xp);
        }

    }

    //Now for the satisfaction stats:
    uint32_t tags_ok=0,tags_unique=0,tags_fail=0;
    for (auto xt=0;xt<tp.size();++xt){
        if (!tag_parts_all[xt].empty()){
            ++tags_ok;
            if (tag_parts_all[xt].size()==1) { ++tags_unique; }
        }
        else if (!tag_parts_some[xt].empty()) ++tags_fail;
    }
    std::cout<<"Scoring tags in partitions:   OK: "<<tags_ok<<" ("<<tags_unique<<" single-partition)    Fail: "<<tags_fail<<std::endl;
    return{ ((float)tags_ok)/(tags_ok+tags_fail),((float)tags_unique)/tags_ok};
}


//def exclussions_from_patterns(tag_patterns):
//    node_count=len(tag_patterns[0])
//    #Compute the characteristic co-occurence string for each node
//    exclussion_mx=[[0]*node_count]*node_count
//    #print(co_occurence_mx)
//    for co in tag_patterns:
//        #print(co)
//        for n in range(node_count):
//            #print(n)
//            if co[n]==1:
//                exclussion_mx[n]=[x for x in map(sum,zip(exclussion_mx[n],co))]
//    #find "mutually exclussive" nodes
//    for n in range(node_count):
//        if exclussion_mx[n][n]:
//            exclussion_mx[n]=[float(x)/exclussion_mx[n][n]<1.0/(10*node_count) for x in exclussion_mx[n]]
//        else:
//            exclussion_mx[n]=[0.0 for x in exclussion_mx[n]]
//    return exclussion_mx

std::vector<std::vector<bool>> exclussions_from_patterns (const std::vector<std::vector<bool>> tag_patterns){
    auto node_count=tag_patterns.back().size();
    //first compute an exclussion mx based on tags that include one but not the other node for all vs. all
    std::vector<std::vector<float>> emxf(node_count);
    std::vector<std::vector<bool>> emxb(node_count);
    for (auto &e:emxf) e.resize(node_count);
    for (auto &e:emxb) e.resize(node_count);

    for (auto &co:tag_patterns){
        for (auto n=0;n<node_count;++n){
            if (co[n]) for (auto nn=0;nn<node_count;++nn) emxf[n][nn]+=co[nn];
        }
    }

    std::cout<<"Node exclusion matrix (float):"<<std::endl;
    for (auto stp:emxf) {
        for (auto p:stp) std::cout<<" "<<p;
        std::cout<<std::endl;
    }

    //If exclussion is more than 1/10*node_count in the number of nodes, then it is "real"?
    //TODO: the total number of reads/tags for a particular node influences its exclusion from others (size/map)
    for (auto n=0;n<node_count;++n){
        if (emxf[n][n]!=0) for (auto nn=0;nn<node_count;++nn)
                emxb[n][nn]=emxf[n][nn]/emxf[n][n]<1.0/(5*node_count);
    }
    return emxb;
}

std::vector<std::vector<bool>> GraphPartitioner::generate_partitions(const SequenceSubGraph subgraph,
                                                                    const std::vector<std::vector<bool>> tag_patterns,
                                                                    unsigned p) {

    auto tp=tag_patterns;
    if (tag_patterns.size()==0){
        tp=tags_patterns(subgraph);
    }
    if (tp.size()==0) return {};

    auto emxb=exclussions_from_patterns(tp);

    std::cout<<"Node exclusion matrix:"<<std::endl;
    for (auto stp:emxb) {
        for (auto p:stp) std::cout<<" "<<(p ? 1:0);
        std::cout<<std::endl;
    }

    std::vector<std::vector<bool>> partitions;
    std::vector<unsigned> node_exclussions(tp.back().size());

    for (auto xp=0;xp<p;++xp){
        //pick the node with the highest "exclussion score" (loneliest first)
        //exclussion score of n is the sum of 1/(1+exclussions[m]) for every node m that n excludes.
        float max_score=0;
        uint64_t max_score_node=0;
        for (auto n=0;n<node_exclussions.size();++n) {
            //compute the score
            float nscore=0;
            for (auto m=0;m<node_exclussions.size();++m){
                if (emxb[n][m]) nscore+=1.0/(1+node_exclussions[m]);
            }
            if (nscore>max_score){
                max_score=nscore;
                max_score_node=n;
            }
        }
        //create a partition with the node, and update exclussions list (+1 to exclussions to every node not in every partition)
        std::vector<bool> new_part(node_exclussions.size());
        for (auto i=0;i<node_exclussions.size();++i) {
            new_part[i]=!emxb[max_score_node][i];
            if (!new_part[i]) ++node_exclussions[i];
        }
        partitions.push_back(new_part);
    }

    //compute general exclussion score for every node from every partition (% of nodes in partition excluding it)
    std::vector<std::vector<float>> pes;
    for (auto n=0;n<node_exclussions.size();++n) {
        std::vector<float> spex;
        for (auto xp=0;xp<p;++xp){
            auto pnodes=0;
            auto pnodes_ex=0;
            for (auto i=0;i<node_exclussions.size();++i){
                if (partitions[xp][i]) {
                    ++pnodes;
                    if (emxb[i][n]) ++pnodes_ex;
                }
            }
            spex.push_back(((float)pnodes_ex)/pnodes);
        }
        pes.push_back(spex);
    }

    //try including unused nodes in their least exclusive partition
    for (auto n=0;n<node_exclussions.size();++n) {
        bool used=false;
        for (auto p:partitions) if (p[n]) used=true;
        if (used) continue;
        float min_ex=1;
        unsigned min_part=0;
        for (auto xp=0;xp<partitions.size();++xp) if (pes[n][xp]<min_ex) {min_ex=pes[n][xp];min_part=xp;};
        auto prev_score=score_partition_set(subgraph,partitions,tp);
        partitions[min_part][n]=1;
        auto new_score=score_partition_set(subgraph,partitions,tp);
        if (prev_score.first>new_score.first) partitions[min_part][n]=0;
    }


    //TODO: validate partitions as subgraphs with valid conectivity?

    return partitions;
}
//
// Created by Bernardo Clavijo (EI) on 18/07/2018.
//

#include "LocalHaplotypeAssembler.hpp"
#include "GraphMaker.hpp"

LocalHaplotypeAssembler::LocalHaplotypeAssembler(WorkSpace &_ws, std::vector<sgNodeID_t> _backbone): ws(_ws),backbone(_backbone) {
    std::cout<<"Creating a LocalHaplotypeAssembler instance from backbone"<<std::endl;
    std::cout<<"Backbone nodes:";
    for (auto n:backbone) std::cout<<" "<<n;
    std::cout<<std::endl;
    std::cout<<"Filling Candidate Tag Set..."<<std::endl;
    //Get tag reads in the nodes, and how many nodes with reads counts.
    std::map<bsg10xTag ,std::pair<uint32_t , uint32_t >> tagcounts; //tag -> nodes, reads
    for (auto &ln:backbone) {
        //partial counts in a single node
        std::map<bsg10xTag ,uint32_t> ntagcounts;
        for (auto rm:ws.linked_read_mappers[0].reads_in_node[llabs(ln)]){
            auto tag=ws.linked_read_datastores[0].get_read_tag(rm.read_id);
            if (tag==0) continue;
            ++ntagcounts[tag];
        }
        //added to totals
        for (auto ntc:ntagcounts) {
            ++tagcounts[ntc.first].first;
            tagcounts[ntc.first].second+=ntc.second;
        }
    }
    uint64_t total_reads=0;
    for (auto tc:tagcounts) {
        auto & tag=tc.first;
        auto & tag_backbone_nodes=tc.second.first;
        auto & tag_backbone_reads=tc.second.second;
        auto tag_total_reads=ws.linked_read_datastores[0].get_tag_reads(tc.first).size();//TODO: this can be optimised by not creating the vector
        if (tc.second.first>1 or tag_total_reads/50<tag_backbone_reads and tag_total_reads>10) {
            tagSet.insert(tc.first);
            total_reads+=tag_total_reads;
        }
    }
    std::cout<<"Local tags: "<<tagSet.size()<<" reads: "<<total_reads<<std::endl;
    std::cout<<"Creating set of relevant paired reads..."<<std::endl;
    for (auto prl=0;prl<ws.paired_read_mappers.size();++prl) {
        paired_reads.emplace_back(std::make_pair(prl,std::vector<uint64_t>()));
        for (auto &ln:backbone) {
            auto nreads=ws.paired_read_mappers[prl].get_node_readpairs_ids(ln);
            paired_reads.back().second.insert(paired_reads.back().second.end(),nreads.begin(),nreads.end());
        }
        std::cout<<paired_reads.back().second.size()<<" reads from "<<ws.paired_read_datastores[prl].filename<<std::endl;
    }
    std::cout<<"LocalHaplotypeAssembler created!"<<std::endl;

}

void LocalHaplotypeAssembler::assemble(int k, int min_cov, bool tag_cov){
    std::cout<<"Creating an uncleaned DBG"<<std::endl;
    //std::cout << "creating DBG for line #" << i << std::endl;
    std::ofstream anchf("local_dbg_" + std::to_string(backbone[0]) + "_anchors.fasta");
    for (auto n:backbone){
        anchf<<">seq"<<llabs(n)<<std::endl;
        anchf<<ws.sg.nodes[llabs(n)].sequence<<std::endl;

    }
    {
        BufferedLRSequenceGetter blrsg(ws.linked_read_datastores[0], 200000, 1000);
        auto ltkmers128 = ws.linked_read_datastores[0].get_tags_kmers128(k, min_cov, tagSet, blrsg, tag_cov);
        SequenceGraph dbg;
        GraphMaker gm(dbg);
        gm.new_graph_from_kmerset_trivial128(ltkmers128, k);
        dbg.write_to_gfa("local_dbg_" + std::to_string(backbone[0]) + "_uncleaned.gfa");
        gm.tip_clipping(200);
        gm.remove_small_unconnected(500);
        dbg.write_to_gfa("local_dbg_" + std::to_string(backbone[0]) + ".gfa");
    }
    {
        std::unordered_set<__uint128_t> allpr_kmers;
        for (auto pr:paired_reads) {
            auto lkmers128 = ws.paired_read_datastores[pr.first].get_reads_kmers128(k, 1, pr.second);
            allpr_kmers.insert(lkmers128.begin(),lkmers128.end());
        }
        SequenceGraph dbg;
        GraphMaker gm(dbg);
        gm.new_graph_from_kmerset_trivial128(allpr_kmers, k);
        dbg.write_to_gfa("local_dbg_" + std::to_string(backbone[0]) + "_pr_uncleaned.gfa");
        gm.tip_clipping(200);
        gm.remove_small_unconnected(500);
        dbg.write_to_gfa("local_dbg_" + std::to_string(backbone[0]) + "_pr.gfa");
    }
    std::cout<<"Analising junctions, one by one"<<std::endl;
    for (auto i=0;i<backbone.size()-1;++i){
        std::cout<<"Tring to joing "<<backbone[i]<<" (-) -> (+) "<<backbone[i+1]<<std::endl;
    }

}

void LocalHaplotypeAssembler::write_problem(std::string prefix) {
    std::ofstream output_file(prefix+".bsglhap");
    uint64_t count;

    //write down backbone;
    std::cout<<"writing backbone"<<std::endl;
    count=backbone.size();
    output_file.write((char *)&count,sizeof(count));
    output_file.write((char *)backbone.data(),count*sizeof(backbone[0]));

    //write down 10x tags;
    std::cout<<"writing tags"<<std::endl;
    count=tagSet.size();
    output_file.write((char *)&count,sizeof(count));
    for (auto &t:tagSet) output_file.write((char *)&t,sizeof(t));

    //write down paired read ids;
    std::cout<<"writing pair reads"<<std::endl;
    count=paired_reads.size();
    output_file.write((char *)&count,sizeof(count));
    for (auto &p:paired_reads){
        output_file.write((char *)&p.first,sizeof(p.first));
        count=p.second.size();
        output_file.write((char *)&count,sizeof(count));
        output_file.write((char *)p.second.data(),count*sizeof(p.second[0]));
    }

}

void LocalHaplotypeAssembler::write_full(std::string prefix, bool keep_full_sg=false) {
    //create a slimmed down bsgws
    std::set<sgNodeID_t> keep_nodes;
    for (auto n:backbone) keep_nodes.insert(llabs(n));
    WorkSpace sdws;
    for (auto n=1;n<ws.sg.nodes.size();++n) {
        if (keep_nodes.count(n)) sdws.sg.nodes.push_back(ws.sg.nodes[n]);
        else sdws.sg.nodes.emplace_back("",sgNodeDeleted);
    }

    //TODO: datastore:write_subdatastore_file

    std::ofstream(prefix+".bsglha");
}
//
// Created by Bernardo Clavijo (EI) on 18/07/2018.
//

#include "LocalHaplotypeAssembler.hpp"
#include "GraphMaker.hpp"

void LocalHaplotypeAssembler::init_from_backbone( std::vector<sgNodeID_t> _backbone) {
    backbone=_backbone;
    std::cout<<"Creating a LocalHaplotypeAssembler instance from backbone"<<std::endl;
    std::cout<<"Backbone nodes:";
    for (auto n:backbone) std::cout<<" "<<n;
    std::cout<<std::endl;
    for (auto n:backbone) backbone_nodes.emplace_back(ws.sg.nodes[llabs(n)].sequence);

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
        if (paired_reads.back().second.empty()) paired_reads.pop_back();
        else std::cout<<paired_reads.back().second.size()<<" reads from "<<ws.paired_read_datastores[prl].filename<<std::endl;
    }
    std::cout<<"LocalHaplotypeAssembler created!"<<std::endl;

}

void LocalHaplotypeAssembler::init_from_file(std::string problem_file) {
    std::cout<<"Loading LocalHaplotypeAssembler problem from file: "<<problem_file<<std::endl;

    std::ifstream input_file(problem_file);
    uint64_t count;

    //load backbone;
    input_file.read((char *)&count,sizeof(count));
    backbone.resize(count);
    input_file.read((char *)backbone.data(),count*sizeof(backbone[0]));
    for (auto n:backbone) backbone_nodes.emplace_back(ws.sg.nodes[llabs(n)].sequence);

    //load 10x tags;
    input_file.read((char *)&count,sizeof(count));
    bsg10xTag t;
    while (count--) {
        input_file.read((char *) &t, sizeof(t));
        tagSet.insert(t);
    }

    //load paired read ids;
    input_file.read((char *)&count,sizeof(count));
    paired_reads.resize(count);
    for (auto &p:paired_reads){
        input_file.read((char *)&p.first,sizeof(p.first));
        input_file.read((char *)&count,sizeof(count));
        p.second.resize(count);
        input_file.read((char *)p.second.data(),count*sizeof(p.second[0]));
    }

    std::cout<<"LocalHaplotypeAssembler created!"<<std::endl;

}

void LocalHaplotypeAssembler::assemble(int k, int min_cov, bool tag_cov, std::string output_prefix){
    if (output_prefix.empty()) output_prefix="local_dbg_" + std::to_string(backbone[0]);
    std::cout<<"Creating an uncleaned DBG"<<std::endl;
    //std::cout << "creating DBG for line #" << i << std::endl;
    std::ofstream anchf(output_prefix + "_anchors.fasta");
    for (auto n=0;n<backbone.size();++n){
        anchf<<">seq"<<llabs(backbone[n])<<std::endl;
        anchf<<backbone_nodes[n].sequence<<std::endl;

    }
    BufferedLRSequenceGetter blrsg(ws.linked_read_datastores[0], 200000, 1000);
    auto ltkmers128 = ws.linked_read_datastores[0].get_tags_kmers128(k, min_cov, tagSet, blrsg, tag_cov);
    GraphMaker gm(assembly);
    gm.new_graph_from_kmerset_trivial128(ltkmers128, k);
//        dbg.write_to_gfa(output_prefix + "_uncleaned.gfa");
    gm.tip_clipping(200);
    gm.remove_small_unconnected(500);
    path_all_reads();
    assembly.write_to_gfa(output_prefix + ".gfa");
    //std::cout<<"Analising junctions, one by one"<<std::endl;
    //for (auto i=0;i<backbone.size()-1;++i){
    //    std::cout<<"Tring to joing "<<backbone[i]<<" (-) -> (+) "<<backbone[i+1]<<std::endl;
    //}

}

void add_readkmer_nodes_lha(std::vector<sgNodeID_t> & kmernodes, std::vector<std::pair<uint64_t,bool>> & readkmers, std::unordered_map<uint64_t, graphPosition> & index, bool rev){
    //TODO allow for a minimum of kmers to count the hit?
    if (not rev) {
        for (auto rki=readkmers.begin();rki<readkmers.end();++rki) {
            auto kp=index.find(rki->first);
            if (kp==index.end()) continue; //kmer not found
            auto node=(rki->second ? -kp->second.node:kp->second.node);
            if (kmernodes.empty() or kmernodes.back()!=node) kmernodes.emplace_back(node);
        }
    }
    else {
        for (auto rki=readkmers.rbegin();rki<readkmers.rend();++rki) {
            auto kp=index.find(rki->first);
            if (kp==index.end()) continue; //kmer not found
            auto node=(rki->second ? kp->second.node:-kp->second.node);
            if (kmernodes.empty() or kmernodes.back()!=node) kmernodes.emplace_back(node);
        }
    }

}

void LocalHaplotypeAssembler::path_all_reads() {

    assembly.create_index();

    //now populate the linked read paths first

    CStringKMerFactory cskf(31);
    std::vector<std::pair<sgNodeID_t ,sgNodeID_t >> nodeproximity_thread;
    std::vector<std::pair<uint64_t,bool>> read1kmers,read2kmers;
    std::vector<sgNodeID_t> kmernodes;

    //BufferedPairedSequenceGetter bprsg(ws.paired_read_datastores[lib], 1000000, 1000);
    BufferedLRSequenceGetter blrsg(ws.linked_read_datastores[0],100000,1000);

    for (auto t:tagSet) {
        for (auto rid : ws.linked_read_datastores[0].get_tag_reads(t)) {
            //std::cout<<"analising reads "<<rid<<" and "<<rid+1<<std::endl;
            if (rid%2!=1) continue;

            read1kmers.clear();
            read2kmers.clear();
            kmernodes.clear();

            cskf.create_kmers_direction(read1kmers, blrsg.get_read_sequence(rid));
            cskf.create_kmers_direction(read2kmers, blrsg.get_read_sequence(rid + 1));
            //first put the kmers from read 1 in there;
            add_readkmer_nodes_lha(kmernodes, read1kmers, assembly.kmer_to_graphposition, false);
            //for (auto kn:kmernodes) std::cout<<" "<<kn; std::cout<<std::endl;
            add_readkmer_nodes_lha(kmernodes, read2kmers, assembly.kmer_to_graphposition, true);
            //for (auto kn:kmernodes) std::cout<<" "<<kn; std::cout<<std::endl;

            linkedread_paths.emplace_back(kmernodes);
        }
    }
    sglib::OutputLog()<<linkedread_paths.size()<<" linked-reads paths created!"<<std::endl;

    //now do the same for each paired library
    for (auto lpr:paired_reads) {
        BufferedPairedSequenceGetter bprsg(ws.paired_read_datastores[lpr.first],100000,1000);
        for (auto rid : lpr.second) {
            //std::cout<<"analising reads "<<rid<<" and "<<rid+1<<std::endl;
            if (rid%2!=1) continue;

            read1kmers.clear();
            read2kmers.clear();
            kmernodes.clear();

            cskf.create_kmers_direction(read1kmers, bprsg.get_read_sequence(rid));
            cskf.create_kmers_direction(read2kmers, bprsg.get_read_sequence(rid + 1));
            //first put the kmers from read 1 in there;
            add_readkmer_nodes_lha(kmernodes, read1kmers, assembly.kmer_to_graphposition, false);
            //for (auto kn:kmernodes) std::cout<<" "<<kn; std::cout<<std::endl;
            add_readkmer_nodes_lha(kmernodes, read2kmers, assembly.kmer_to_graphposition, true);
            //for (auto kn:kmernodes) std::cout<<" "<<kn; std::cout<<std::endl;

            pairedread_paths.emplace_back(kmernodes);
        }
    }
    sglib::OutputLog()<<pairedread_paths.size()<<" paired-reads paths created!"<<std::endl;

}

void LocalHaplotypeAssembler::write_problem(std::string prefix) {
    std::ofstream output_file(prefix+".bsglhap");
    uint64_t count;

    //write down backbone;
    count=backbone.size();
    output_file.write((char *)&count,sizeof(count));
    output_file.write((char *)backbone.data(),count*sizeof(backbone[0]));

    //write down 10x tags;
    count=tagSet.size();
    output_file.write((char *)&count,sizeof(count));
    for (auto &t:tagSet) output_file.write((char *)&t,sizeof(t));

    //write down paired read ids;
    count=paired_reads.size();
    output_file.write((char *)&count,sizeof(count));
    for (auto &p:paired_reads){
        output_file.write((char *)&p.first,sizeof(p.first));
        count=p.second.size();
        output_file.write((char *)&count,sizeof(count));
        output_file.write((char *)p.second.data(),count*sizeof(p.second[0]));
    }

}

void LocalHaplotypeAssembler::write_full(std::string prefix) {
    std::ofstream output_file(prefix+".bsglhapf");
    uint64_t count;
    //write down backbone;
    count=backbone.size();
    output_file.write((char *)&count,sizeof(count));
    output_file.write((char *)backbone.data(),count*sizeof(backbone[0]));
    //write  backbone node's sequences
    count=backbone.size();
    output_file.write((char *)&count,sizeof(count));
    for (auto n:backbone_nodes){
        count=n.sequence.size();
        output_file.write((char *)&count,sizeof(count));
        output_file.write(n.sequence.c_str(),count);
    }
    //write down 10x tags;
    count=tagSet.size();
    output_file.write((char *)&count,sizeof(count));
    for (auto &t:tagSet) output_file.write((char *)&t,sizeof(t));

    //write the condensation of the 10x workspace
    ws.linked_read_datastores[0].write_selection(output_file,tagSet);
    //write the condensation of the paired reads workspace.
    count=paired_reads.size();
    output_file.write((char *)&count,sizeof(count));
    for (auto &p:paired_reads){
        ws.paired_read_datastores[p.first].write_selection(output_file,p.second);
    }
}


void LocalHaplotypeAssembler::init_from_full_file(std::string full_file) {
    std::cout<<"Loading LocalHaplotypeAssembler problem from file: "<<full_file<<std::endl;

    std::ifstream input_file(full_file);
    uint64_t count;
    //load backbone;
    input_file.read((char *)&count,sizeof(count));
    backbone.resize(count);
    input_file.read((char *)backbone.data(),count*sizeof(backbone[0]));
    //write  backbone node's sequences
    std::cout<<"reading backbone sequences"<<std::endl;
    uint64_t bcount;
    input_file.read((char *)&bcount,sizeof(bcount));
    for (auto n=0;n<bcount;++n){
        std::string seq;
        input_file.read((char *)&count,sizeof(count));
        seq.resize(count);
        input_file.read((char *)seq.data(),count);
        backbone_nodes.emplace_back(seq);
    }
    //write down 10x tags;
    std::cout<<"writing tags"<<std::endl;
    input_file.read((char *)&count,sizeof(count));
    for (auto n=0;n<count;++n) {
        bsg10xTag t;
        input_file.read((char *) &t, sizeof(t));
        tagSet.insert(t);
    }

    //write the condensation of the 10x workspace
    std::cout<<"writing linked reads"<<std::endl;
    ws.linked_read_datastores.emplace_back();
    ws.linked_read_datastores[0].load_from_stream(full_file,input_file);
    //ws.linked_read_datastores[0].write_selection(output_file,tagSet);
    //write the condensation of the paired reads workspace.
    std::cout<<"writing paired reads"<<std::endl;
    input_file.read((char *)&count,sizeof(count));
    for (auto n=0;n<count;++n){
        ws.paired_read_datastores.emplace_back();
        ws.paired_read_datastores.back().load_from_stream(full_file,input_file);
        std::vector<uint64_t> read_ids;
        read_ids.reserve(ws.paired_read_datastores.back().size());
        for(uint64_t i=0;i<ws.paired_read_datastores.back().size();++i) read_ids.emplace_back(i);
        paired_reads.emplace_back(std::make_pair(n,read_ids));
    }

    std::cout<<"LocalHaplotypeAssembler created!"<<std::endl;

}
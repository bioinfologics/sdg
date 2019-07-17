//
// Created by Bernardo Clavijo (EI) on 10/02/2018.
//

#include "LinkedReadsDatastore.hpp"
#include <sdglib/utilities/OutputLog.hpp>
#include <sdglib/mappers/LinkedReadsMapper.hpp>
#include <sdglib/workspace/WorkSpace.hpp>
#include <fstream>
#include <strings.h>
#include <cstring>

const sdgVersion_t LinkedReadsDatastore::min_compat = 0x0003;

std::string bsg10xTag_to_seq(bsg10xTag tag, uint8_t k) {
    std::string seq;
    seq.reserve(k);
    for (int shift = (k - 1) * 2; shift >= 0; shift -= 2) {
//std::cout<<"kmer: "
        switch ((tag >> shift) % 4) {
            case 0:
                seq.push_back('A');
                break;
            case 1:
                seq.push_back('C');
                break;
            case 2:
                seq.push_back('G');
                break;
            case 3:
                seq.push_back('T');
                break;
        }
    }
    return seq;
}

void LinkedReadsDatastore::print_status() const {
    uint64_t tagcount=0;
    bsg10xTag prevtag=0;
    for (auto t:read_tag) {
        if (t != prevtag) {
            ++tagcount;
            t = prevtag;
        }
    }
    sdglib::OutputLog()<<"LinkedRead Datastore from "<<filename<<" contains "<<size()-1<<" reads, with a maximum size of " << readsize << " in "<< tagcount <<" tags."<<std::endl;
    mapper.print_status();
}

void LinkedReadsDatastore::build_from_fastq(std::string output_filename, std::string default_name, std::string read1_filename,
                                            std::string read2_filename,
                                            LinkedReadsFormat format, uint64_t readsize, size_t chunksize) {
    std::vector<uint32_t> read_tag;
    //std::cout<<"Memory used by every read's entry:"<< sizeof(LinkedRead)<<std::endl;
    //read each read, put it on the index and on the appropriate tag
    sdglib::OutputLog(sdglib::LogLevels::INFO)<<"Creating Datastore Index from "<<read1_filename<<" | "<<read2_filename<<std::endl;
    auto fd1=gzopen(read1_filename.c_str(),"r");
    if (!fd1) {
        std::cerr << "Failed to open " << read1_filename <<": " << strerror(errno);
        throw std::runtime_error("Could not open " + read1_filename);
    }
    auto fd2=gzopen(read2_filename.c_str(),"r");

    if (!fd2) {
        std::cerr << "Failed to open " << read2_filename <<": " << strerror(errno);
        throw std::runtime_error("Could not open " + read2_filename);
    }
    char readbuffer[1000];
    uint64_t r1offset,r2offset;
    uint64_t tagged_reads=0;
    //first, build an index of tags and offsets
    sdglib::OutputLog()<<"Building tag sorted chunks of "<<chunksize<<" pairs"<<std::endl;
    std::vector<LinkedReadData> readdatav;
    readdatav.reserve(chunksize);
    std::vector<std::ifstream> chunkfiles;
    std::vector<LinkedReadData> next_in_chunk;
    LinkedReadData currrent_read;
    //First, create the chunk files
    uint64_t pairs=0;
    while (!gzeof(fd1) and !gzeof(fd2)) {
        bsg10xTag newtag = 0;
        if (format==LinkedReadsFormat::UCDavis) {
            //LinkedRead r1,r2;
            if (NULL == gzgets(fd1, readbuffer, 999)) continue;
            if(!readbuffer[0] == '@') {
                throw std::runtime_error("Please check: " + read1_filename + ", it seems to be missing a header");
            }
            //Tag to number from r1's name
            for (auto i = 1; i < 17; ++i) {
                newtag <<= 2;
                if (readbuffer[i] == 'C') newtag += 1;
                else if (readbuffer[i] == 'G') newtag += 2;
                else if (readbuffer[i] == 'T') newtag += 3;
                else if (readbuffer[i] != 'A') {
                    newtag = 0;
                    break;
                }
            }
            currrent_read.tag=newtag;
            if (NULL == gzgets(fd1, readbuffer, 999)) continue;
            currrent_read.seq1=std::string(readbuffer);
            if (currrent_read.seq1.back()=='\n') currrent_read.seq1.resize(currrent_read.seq1.size()-1);
            if (NULL == gzgets(fd1, readbuffer, 999)) continue;
            if(!readbuffer[0] == '+') {
                throw std::runtime_error("Please check: " + read1_filename + ", it seems to be desynchronised a header");
            }
            if (NULL == gzgets(fd1, readbuffer, 999)) continue;

            if (NULL == gzgets(fd2, readbuffer, 999)) continue;
            if(!readbuffer[0] == '@') {
                throw std::runtime_error("Please check: " + read2_filename + ", it seems to be missing a header");
            }
            if (NULL == gzgets(fd2, readbuffer, 999)) continue;
            currrent_read.seq2=std::string(readbuffer);
            if (currrent_read.seq2.back()=='\n') currrent_read.seq2.resize(currrent_read.seq2.size()-1);
            if (NULL == gzgets(fd2, readbuffer, 999)) continue;
            if(!readbuffer[0] == '+') {
                throw std::runtime_error("Please check: " + read2_filename + ", it seems to be desynchronised a header");
            }
            if (NULL == gzgets(fd2, readbuffer, 999)) continue;
        }
        else if (format==LinkedReadsFormat::raw){
            if (NULL == gzgets(fd1, readbuffer, 999)) continue;
            if(!readbuffer[0] == '@') {
                throw std::runtime_error("Please check: " + read1_filename + ", it seems to be missing a header");
            }
            if (NULL == gzgets(fd1, readbuffer, 999)) continue;
            for (auto i = 0; i < 16; ++i) {
                newtag <<= 2;
                if (readbuffer[i] == 'C') newtag += 1;
                else if (readbuffer[i] == 'G') newtag += 2;
                else if (readbuffer[i] == 'T') newtag += 3;
                else if (readbuffer[i] != 'A') {
                    newtag = 0;
                    break;
                }
            }
            currrent_read.tag=newtag;
            currrent_read.seq1=std::string(readbuffer + 16 + 7);
            if (currrent_read.seq1.back()=='\n') currrent_read.seq1.resize(currrent_read.seq1.size()-1);
            if (NULL == gzgets(fd1, readbuffer, 999)) continue;
            if(!readbuffer[0] == '+') {
                throw std::runtime_error("Please check: " + read1_filename + ", it seems to be desynchronised a header");
            }
            if (NULL == gzgets(fd1, readbuffer, 999)) continue;

            if (NULL == gzgets(fd2, readbuffer, 999)) continue;
            if(!readbuffer[0] == '@') {
                throw std::runtime_error("Please check: " + read2_filename + ", it seems to be missing a header");
            }
            if (NULL == gzgets(fd2, readbuffer, 999)) continue;
            currrent_read.seq2=std::string(readbuffer);
            if (currrent_read.seq2.back()=='\n') currrent_read.seq2.resize(currrent_read.seq2.size()-1);
            if (NULL == gzgets(fd2, readbuffer, 999)) continue;
            if(!readbuffer[0] == '+') {
                throw std::runtime_error("Please check: " + read2_filename + ", it seems to be desynchronised a header");
            }
            if (NULL == gzgets(fd2, readbuffer, 999)) continue;
        }
        if (0 != newtag) tagged_reads += 2;
        ++pairs;
        readdatav.push_back(currrent_read);
        if (readdatav.size()==chunksize){
            //sort
            std::sort(readdatav.begin(),readdatav.end());
            //dump
            std::ofstream ofile("sorted_chunk_"+std::to_string(chunkfiles.size())+".data");
            sdglib::OutputLog()<<readdatav.size()<<" pairs dumping on chunk "<<chunkfiles.size()<<std::endl;
            if (!ofile) {
                std::cerr << "Failed to open " << ("sorted_chunk_" + std::to_string(chunkfiles.size()) + ".data") <<": " << strerror(errno);
                throw std::runtime_error("Could not open " + ("sorted_chunk_" + std::to_string(chunkfiles.size()) + ".data") );

            }

            sdglib::OutputLog()<<readdatav.size()<<" pairs dumping on chunk "<<chunkfiles.size()<<std::endl;
            //add file to vector of files
            char buffer[2*readsize+2];
            for (auto &r:readdatav){
                ofile.write((const char * ) &r.tag,sizeof(r.tag));
                bzero(buffer,2*readsize+2);
                memcpy(buffer,r.seq1.data(),(r.seq1.size()>readsize ? readsize : r.seq1.size()));
                memcpy(buffer+readsize+1,r.seq2.data(),(r.seq2.size()>readsize ? readsize : r.seq2.size()));
                ofile.write(buffer,2*readsize+2);
            }
            ofile.close();
            chunkfiles.emplace_back("sorted_chunk_"+std::to_string(chunkfiles.size())+".data");
            if (!chunkfiles.back()) {
                std::cerr << "Failed to open " << ("sorted_chunk_"+std::to_string(chunkfiles.size())+".data") <<": " << strerror(errno);
                throw std::runtime_error("Could not open " + ("sorted_chunk_"+std::to_string(chunkfiles.size())+".data") );
            }
            readdatav.clear();
            sdglib::OutputLog()<<"dumped!"<<std::endl;
        }
    }
    if (readdatav.size()>0) {
        //sort
        std::sort(readdatav.begin(), readdatav.end());
        //dump
        std::ofstream ofile("sorted_chunk_" + std::to_string(chunkfiles.size()) + ".data");
        sdglib::OutputLog() << readdatav.size() << " pairs dumping on chunk " << chunkfiles.size() << std::endl;
        if (!ofile) {
            std::cerr << "Failed to open " << ("sorted_chunk_" + std::to_string(chunkfiles.size()) + ".data") <<": " << strerror(errno);
            throw std::runtime_error("Could not open " + ("sorted_chunk_" + std::to_string(chunkfiles.size()) + ".data") );

        }
        sdglib::OutputLog() << readdatav.size() << " pairs dumping on chunk " << chunkfiles.size() << std::endl;
        //add file to vector of files
        char buffer[2 * readsize + 2];
        for (auto &r:readdatav) {
            ofile.write((const char *) &r.tag, sizeof(r.tag));
            bzero(buffer, 2 * readsize + 2);
            memcpy(buffer, r.seq1.data(), (r.seq1.size() > readsize ? readsize : r.seq1.size()));
            memcpy(buffer + readsize + 1, r.seq2.data(), (r.seq2.size() > readsize ? readsize : r.seq2.size()));
            ofile.write(buffer, 2 * readsize + 2);
        }
        ofile.close();
        chunkfiles.emplace_back("sorted_chunk_" + std::to_string(chunkfiles.size()) + ".data");
        if (!chunkfiles.back()) {
            std::cerr << "Failed to open " << ("sorted_chunk_" + std::to_string(chunkfiles.size()) + ".data") <<": " << strerror(errno);
            throw std::runtime_error("Could not open " + ("sorted_chunk_" + std::to_string(chunkfiles.size()) + ".data"));
        }
        readdatav.clear();
        sdglib::OutputLog() << "dumped!" << std::endl;
    }
    sdglib::OutputLog() << "performing merge from disk" << std::endl;
    //TODO: save space first for the tag index!!!
    std::ofstream output(output_filename.c_str());
    if (!output) {
        std::cerr << "Failed to open " << output_filename <<": " << strerror(errno);
        throw std::runtime_error("Could not open " + output_filename);
    }

    output.write((const char *) &SDG_MAGIC, sizeof(SDG_MAGIC));
    output.write((const char *) &SDG_VN, sizeof(SDG_VN));
    SDG_FILETYPE type(LinkedDS_FT);
    output.write((char *) &type, sizeof(type));

    sdglib::write_string(output, default_name);

    output.write((const char *) &readsize,sizeof(readsize));
    read_tag.resize(pairs);
    sdglib::OutputLog() << "leaving space for " <<pairs<<" read_tag entries"<< std::endl;

    sdglib::write_flat_vector(output,read_tag);

    //multi_way merge of the chunks
    int openfiles=chunkfiles.size();
    bsg10xTag next_tags[chunkfiles.size()];
    char buffer[2 * readsize + 2];
    //read a bit from each file
    for (auto i=0;i<chunkfiles.size();++i) chunkfiles[i].read((char *)&next_tags[i],sizeof(bsg10xTag));
    read_tag.clear();
    while(openfiles){
        bsg10xTag mintag=UINT32_MAX;
        //find the minimum tag
        for (auto &t:next_tags) if (t<mintag) mintag=t;
        //copy from each file till the next tag is > than minimum
        for (auto i=0;i<chunkfiles.size();++i){
            while (!chunkfiles[i].eof() and next_tags[i]==mintag){
                //read buffer from file and write to final file
                chunkfiles[i].read(buffer,2*readsize+2);
                output.write(buffer,2*readsize+2);
                read_tag.push_back(mintag);
                //read next tag... eof? tag=UINT32_MAX
                chunkfiles[i].read((char *)&next_tags[i],sizeof(bsg10xTag));
                if (chunkfiles[i].eof()) {
                    --openfiles;
                    next_tags[i]=UINT32_MAX;
                    sdglib::OutputLog() << "chunk "<<i<<" finished"<<std::endl;
                }
            }
        }
    }

    // Write empty mapper data
    output.write((char *) &SDG_MAGIC, sizeof(SDG_MAGIC));
    output.write((char *) &SDG_VN, sizeof(SDG_VN));
    type = PairedMap_FT;
    output.write((char *) &type, sizeof(type));

    std::vector<sgNodeID_t> read_to_node;
    sdglib::write_flat_vector(output, read_to_node);
    std::vector<std::vector<ReadMapping>> reads_in_node;
    sdglib::write_flat_vectorvector(output, reads_in_node);

    // Done

    //go back to the beginning of the file and write the read_tag part again
    output.seekp(sizeof(SDG_MAGIC)+sizeof(SDG_VN)+sizeof(type)+sizeof(readsize)+sizeof(uint64_t)+(default_name.size()*sizeof(char)));
    sdglib::OutputLog() << "writing down " <<pairs<<" read_tag entries"<< std::endl;
    sdglib::write_flat_vector(output, read_tag);

    output.close();
    //delete all temporary chunk files
    for (auto &c:chunkfiles) c.close();
    for (auto i=0;i<chunkfiles.size();++i) ::unlink(("sorted_chunk_"+std::to_string(i)+".data").c_str());
    //DONE!
    sdglib::OutputLog(sdglib::LogLevels::INFO)<<"Datastore with "<<(read_tag.size())*2<<" reads, "<<tagged_reads<<" reads with tags"<<std::endl; //and "<<reads_in_tag.size()<<"tags"<<std::endl;
    gzclose(fd1);
    gzclose(fd2);
}

void LinkedReadsDatastore::read(std::ifstream &input_file) {
    //read filename
    uint64_t s;
    sdglib::read_string(input_file, filename);
    sdglib::read_string(input_file, name);

    load_index(filename);
}

void LinkedReadsDatastore::load_index(std::string _filename){
    uint64_t s;
    filename=_filename;
    fd=fopen(filename.c_str(),"r");
    if (!fd) {
        std::cerr << "Failed to open " << filename <<": " << strerror(errno);
        throw std::runtime_error("Could not open " + filename);
    }
    sdgMagic_t magic;
    sdgVersion_t version;
    SDG_FILETYPE type;
    fread((char *) &magic, sizeof(magic),1,fd);
    fread((char *) &version, sizeof(version),1,fd);
    fread((char *) &type, sizeof(type),1,fd);

    if (magic != SDG_MAGIC) {
        throw std::runtime_error("Magic number not present in " + _filename);
    }

    if (version < min_compat) {
        throw std::runtime_error("LinkedReadsDS file version: " + std::to_string(version) + " is not compatible with " + std::to_string(min_compat));
    }

    if (type != LinkedDS_FT) {
        throw std::runtime_error("File type supplied: " + std::to_string(type) + " is not compatible with LinkedDS_FT");
    }

    uint64_t sname=0;
    fread( &sname, sizeof(sname), 1, fd);
    default_name.resize(sname);
    fread( (char *) default_name.data(), sizeof(char), sname, fd);

    fread( &readsize,sizeof(readsize),1,fd);
    fread(&s,sizeof(s),1,fd); read_tag.resize(s);
    fread(read_tag.data(),sizeof(read_tag[0]),read_tag.size(),fd);
    readpos_offset=ftell(fd);
    sdglib::OutputLog()<<"LinkedReadsDatastore open: "<<_filename<<"  max read length: "<<readsize<<" Total reads: " <<size()<<std::endl;
}

void LinkedReadsDatastore::write(std::ofstream &output_file) {
    //read filename
    sdglib::write_string(output_file, filename);
    sdglib::write_string(output_file, name);

    mapper.write(output_file);
}

void LinkedReadsDatastore::write_selection(std::ofstream &output_file, const std::set<bsg10xTag> &tagSet) {
    //write readsize
    output_file.write((const char *) &SDG_MAGIC, sizeof(SDG_MAGIC));
    output_file.write((const char *) &SDG_VN, sizeof(SDG_VN));
    SDG_FILETYPE type(LinkedDS_FT);
    output_file.write((char *) &type, sizeof(type));

    output_file.write((char *) &readsize, sizeof(readsize));
    //create a vector of read tags (including each tag as many times as reads it contains).
    //std::cout << "creating vector of tags" << std::endl;
    std::vector<bsg10xTag> readtags;
    //readtags.push_back(0);
    for (auto trc:get_tag_readcount()) {
        if (tagSet.count(trc.first)) readtags.resize(readtags.size() + trc.second, trc.first);
    }
    //readtags.push_back(0);

    // write tag count and tags.
    uint64_t rts = readtags.size();
    output_file.write((const char *) &rts, sizeof(rts));
    output_file.write((const char *) readtags.data(), sizeof(bsg10xTag) * readtags.size());

    //now for each read in each included tag, just copy the readsize+1 sequence into the file.
    uint64_t totaldata = 0;
    for (auto t:tagSet) {
        for (auto n = std::lower_bound(read_tag.begin(), read_tag.end(), t) - read_tag.begin(); read_tag[n] == t; ++n) {
            auto s = get_read_sequence(n * 2 + 1);
            output_file.write(s.c_str(), readsize + 1);
            s = get_read_sequence(n * 2 + 2);
            output_file.write(s.c_str(), readsize + 1);
            totaldata += 2 * readsize + 2;
        }
    }
}

std::string LinkedReadsDatastore::get_read_sequence(size_t readID) {
    if (readID==0 or readID>size()) return "";
    char buffer[readsize+1];
    size_t read_offset_in_file=readpos_offset+(readsize+1)*(readID-1);
    fseek(fd,read_offset_in_file,SEEK_SET);
    fread(buffer,readsize+1,1,fd);
    return std::string(buffer);
}

std::vector<uint64_t> LinkedReadsDatastore::get_tag_reads(bsg10xTag tag) const {
    std::vector<uint64_t> rids;
    rids.reserve(10000);
    for (auto n=std::lower_bound(read_tag.begin(),read_tag.end(),tag)-read_tag.begin();read_tag[n]==tag;++n) {
        rids.emplace_back( n * 2 + 1);
        rids.emplace_back( n * 2 + 2);
    }
    return rids;
}

bsg10xTag LinkedReadsDatastore::get_read_tag(size_t readID) {
    if (readID==0 or readID>size()) return 0;
    return read_tag[(readID-1)/2];
}

std::vector<std::pair<bsg10xTag, uint32_t>> LinkedReadsDatastore::get_tag_readcount() {
    std::vector<std::pair<bsg10xTag, uint32_t>> readcount;
    auto curr_tag=read_tag[0];
    uint32_t curr_count=0;
    for (auto &t:read_tag){
        if (t!=curr_tag){
            if (curr_tag!=0) readcount.emplace_back(curr_tag,curr_count);
            curr_count=0;
            curr_tag=t;
        }
        ++curr_count;
    }
    if (curr_tag!=0) readcount.emplace_back(curr_tag,curr_count);
    return readcount;
}

void LinkedReadsDatastore::dump_tag_occupancy_histogram(std::string filename) {
    std::ofstream tohf(filename);
    uint64_t oh[10001];
    for (auto &o:oh) o=0;
    for (auto &t:get_tag_readcount()){
        ++oh[(t.second>10000 ? 10000:t.second)];
    }
    for (auto i=0;i<10001;++i)
        if (oh[i]) tohf<<i<<","<<oh[i]<<std::endl;
}

std::unordered_set<uint64_t> LinkedReadsDatastore::get_tags_kmers(int k, int min_tag_cov, std::set<bsg10xTag> tags, ReadSequenceBuffer & blrsg) {
    class StreamKmerFactory : public  KMerFactory {
    public:
        explicit StreamKmerFactory(uint8_t k) : KMerFactory(k){}
        inline void produce_all_kmers(const char * seq, std::vector<uint64_t> &mers){
            // TODO: Adjust for when K is larger than what fits in uint64_t!
            last_unknown=0;
            fkmer=0;
            rkmer=0;
            auto s=seq;
            while (*s!='\0' and *s!='\n') {
                //fkmer: grows from the right (LSB)
                //rkmer: grows from the left (MSB)
                fillKBuf(*s, fkmer, rkmer, last_unknown);
                if (last_unknown >= K) {
                    if (fkmer <= rkmer) {
                        // Is fwd
                        mers.emplace_back(fkmer);
                    } else {
                        // Is bwd
                        mers.emplace_back(rkmer);
                    }
                }
                ++s;
            }
        }
    };
    StreamKmerFactory skf(k);

    //reserve space by counting reads first, save only the integer, do not merge just count and insert in the set
    std::vector<uint64_t> all_kmers;
    std::vector<uint64_t> read_ids;
    for (auto t:tags) {
        auto reads=get_tag_reads(t);
        read_ids.insert(read_ids.end(),reads.begin(),reads.end());
    }
    all_kmers.reserve(read_ids.size()*(readsize-k+1));
    for (auto rid:read_ids){
            skf.produce_all_kmers(blrsg.get_read_sequence(rid),all_kmers);
    }
    std::sort(all_kmers.begin(),all_kmers.end());
    std::unordered_set<uint64_t> kset;
    auto ri=all_kmers.begin();
    auto nri=all_kmers.begin();
    while (ri<all_kmers.end()){
        while (nri<all_kmers.end() and *nri==*ri) ++nri;
        if (nri-ri>=min_tag_cov) kset.insert(*ri);
        ri=nri;
    }
    return std::move(kset);
}

std::unordered_set<__uint128_t, int128_hash> LinkedReadsDatastore::get_tags_kmers128(int k, int min_tag_cov, std::set<bsg10xTag> tags, ReadSequenceBuffer & blrsg, bool count_tag_cvg) {
    class StreamKmerFactory128 : public  KMerFactory128 {
    public:
        explicit StreamKmerFactory128(uint8_t k) : KMerFactory128(k){}
        inline void produce_all_kmers(const char * seq, std::vector<__uint128_t> &mers){
            // TODO: Adjust for when K is larger than what fits in __uint128_t!
            last_unknown=0;
            fkmer=0;
            rkmer=0;
            auto s=seq;
            while (*s!='\0' and *s!='\n') {
                //fkmer: grows from the right (LSB)
                //rkmer: grows from the left (MSB)
                fillKBuf(*s, fkmer, rkmer, last_unknown);
                if (last_unknown >= K) {
                    if (fkmer <= rkmer) {
                        // Is fwd
                        mers.emplace_back(fkmer);
                    } else {
                        // Is bwd
                        mers.emplace_back(rkmer);
                    }
                }
                ++s;
            }
        }
    };
    StreamKmerFactory128 skf(k);

    //reserve space by counting reads first, save only the integer, do not merge just count and insert in the set
    std::vector<__uint128_t> all_kmers;
    std::vector<__uint128_t> read_ids;
    for (auto t:tags) {
        auto reads=get_tag_reads(t);
        read_ids.insert(read_ids.end(),reads.begin(),reads.end());
    }
    all_kmers.reserve(read_ids.size()*(readsize-k+1));
    if (count_tag_cvg) {
        for (auto t:tags) {
            auto tag_first_idx=all_kmers.size();
            for (auto rid:get_tag_reads(t)) {
                skf.produce_all_kmers(blrsg.get_read_sequence(rid), all_kmers);
            }
            std::sort(all_kmers.begin()+tag_first_idx,all_kmers.end());
            all_kmers.erase(std::unique(all_kmers.begin()+tag_first_idx,all_kmers.end()),all_kmers.end());
        }
    } else {
        for (auto rid:read_ids) {
            skf.produce_all_kmers(blrsg.get_read_sequence(rid), all_kmers);
        }
    }
    std::sort(all_kmers.begin(),all_kmers.end());
    std::unordered_set<__uint128_t, int128_hash> kset;
    auto ri=all_kmers.begin();
    auto nri=all_kmers.begin();
    while (ri<all_kmers.end()){
        while (nri<all_kmers.end() and *nri==*ri) ++nri;
        if (nri-ri>=min_tag_cov) kset.insert(*ri);
        ri=nri;
    }
    return std::move(kset);
}

LinkedReadsDatastore::LinkedReadsDatastore(WorkSpace &ws, std::string filename) : ws(ws), mapper(ws, *this){
    load_index(filename);
}

LinkedReadsDatastore::LinkedReadsDatastore(WorkSpace &ws, std::string read1_filename, std::string read2_filename,
                                           std::string output_filename, LinkedReadsFormat format,
                                           std::string default_name, int readsize) : ws(ws), mapper(ws, *this) {
    build_from_fastq(output_filename, default_name, read1_filename, read2_filename, format, readsize, 0);
}

LinkedReadsDatastore::LinkedReadsDatastore(WorkSpace &ws, std::ifstream &infile) : ws(ws), mapper(ws, *this) {
    read(infile);
    mapper.read(infile);
}

LinkedReadsDatastore::LinkedReadsDatastore(WorkSpace &ws, std::string filename, std::ifstream &infile) : ws(ws), mapper(ws, *this) {
    uint64_t s;
    filename=filename;
    fd=fopen(filename.c_str(),"r");
    sdgMagic_t magic;
    sdgVersion_t version;
    SDG_FILETYPE type;
    fread((char *) &magic, sizeof(magic),1,fd);
    fread((char *) &version, sizeof(version),1,fd);
    fread((char *) &type, sizeof(type),1,fd);

    if (magic != SDG_MAGIC) {
        throw std::runtime_error("Magic number not present in " + filename);
    }

    if (version < min_compat) {
        throw std::runtime_error("LinkedDS file version: " + std::to_string(version) + " is not compatible with " + std::to_string(min_compat));
    }

    if (type != LinkedDS_FT) {
        throw std::runtime_error("File type supplied: " + std::to_string(type) + " is not compatible with LinkedDS_FT");
    }

    infile.read( (char *) &readsize,sizeof(readsize));
    infile.read( (char *) &s,sizeof(s));
    read_tag.resize(s);
    infile.read( (char *) read_tag.data(),sizeof(read_tag[0])*s);
    readpos_offset=infile.tellg();
    fseek(fd,readpos_offset,SEEK_SET);
    infile.seekg(2*s*(readsize+1),std::ios_base::cur);
    sdglib::OutputLog()<<"LinkedReadsDatastore open: "<<filename<<"  max read length: "<<readsize<<" Total reads: " <<size()<<std::endl;
}

LinkedReadsDatastore::LinkedReadsDatastore(WorkSpace &ws, const LinkedReadsDatastore &o) : ws(ws), mapper(ws, *this) {
    mapper.read_to_node = o.mapper.read_to_node;
    mapper.reads_in_node = o.mapper.reads_in_node;
    mapper.tag_neighbours = o.mapper.tag_neighbours;
}

LinkedReadsDatastore &LinkedReadsDatastore::operator=(LinkedReadsDatastore const &o) {
    if ( &o == this) return *this;

    filename = o.filename;
    ws = o.ws;
    readsize = o.readsize;
    readpos_offset = o.readpos_offset;
    mapper = o.mapper;
    read_tag = o.read_tag;
    fd = fopen(filename.c_str(), "r");
}

std::ostream &operator<<(std::ostream &os, const LinkedReadsDatastore &lrds) {
    os << "LinkedReadsDatastore" << std::endl;
    os << "Name: " << (lrds.name.empty() ? lrds.default_name : lrds.name) << std::endl;
    lrds.print_status();
}

void BufferedTagKmerizer::get_tag_kmers(bsg10xTag tag) {
    auto read_ids=datastore.get_tag_reads(tag);
    //counts.reserve(counts.size()+read_ids.size()*(datastore.readsize-K+1));
    for (auto rid:read_ids){
        skf.produce_all_kmers(bprsg.get_read_sequence(rid),counts);
    }
}

std::unordered_set<uint64_t> BufferedTagKmerizer::get_tags_kmers(int min_tag_cov, std::set<bsg10xTag> tags) {
    counts.clear();
    for (auto t:tags) get_tag_kmers(t);
    std::cout<<"Sorting "<<counts.size()<<" kmers..."<<std::endl;
    std::sort(counts.begin(),counts.end());
    std::cout<<"...DONE"<<std::endl;
    uint64_t curr_k=UINT64_MAX;
    int curr_count=0;
    auto wi=counts.begin();
    for (auto &k:counts) {
        if (curr_k!=k){
            if (curr_count>=min_tag_cov) {
                *wi=curr_k;
                ++wi;
            }
            curr_count=1;
            curr_k=k;
        }
        else {
            ++curr_count;
        }
    }
    if (curr_count>=min_tag_cov) {
        *wi=counts.back();
        ++wi;
    }
    counts.resize(wi-counts.begin());
    std::unordered_set<uint64_t> r(counts.begin(),counts.end());
    return r;

}





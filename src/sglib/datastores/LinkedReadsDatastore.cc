//
// Created by Bernardo Clavijo (EI) on 10/02/2018.
//

#include "LinkedReadsDatastore.hpp"

void LinkedReadsDatastore::build_from_fastq(std::string read1_filename, std::string read2_filename,
                                            LinkedReadsFormat format) {
    if (format==LinkedReadsFormat::raw) {
        std::cout<<"Raw 10x format not supported yet, datastore not populated!!!"<<std::endl;
        return;
    }
    //std::cout<<"Memory used by every read's entry:"<< sizeof(LinkedRead)<<std::endl;
    //read each read, put it on the index and on the appropriate tag
    sglib::OutputLog(sglib::LogLevels::INFO)<<"Creating Datastore Index from "<<read1_filename<<" | "<<read2_filename<<std::endl;
    filename1=read1_filename;
    filename2=read2_filename;
    fd1=fopen(read1_filename.c_str(),"r");
    fd2=fopen(read2_filename.c_str(),"r");
    char readbuffer[1000];
    read_offset.resize(1);//leave read 0 empty.
    uint64_t r1offset,r2offset;
    uint64_t g1offset=0,g2offset=0;
    uint64_t readCountInFile=0;
    uint64_t tagged_reads=0;
    if (format==LinkedReadsFormat::UCDavis) {
        while (!feof(fd1) and !feof(fd2)) {
            //LinkedRead r1,r2;
            bsg10xTag newtag = 0;
            if (NULL == fgets(readbuffer, 999, fd1)) continue;
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
            r1offset = ftell(fd1);
            if (NULL == fgets(readbuffer, 999, fd1)) continue;
            if (NULL == fgets(readbuffer, 999, fd1)) continue;
            if (NULL == fgets(readbuffer, 999, fd1)) continue;

            if (NULL == fgets(readbuffer, 999, fd2)) continue;
            r2offset = ftell(fd2);
            if (NULL == fgets(readbuffer, 999, fd2)) continue;
            if (NULL == fgets(readbuffer, 999, fd2)) continue;
            if (NULL == fgets(readbuffer, 999, fd2)) continue;

            if (0 == readCountInFile % group_size) {
                g1offset = r1offset;
                group_offset1.push_back(r1offset);
                g2offset = r2offset;
                group_offset2.push_back(r2offset);
            }
            read_offset.push_back(r1offset - g1offset);
            read_offset.push_back(r2offset - g2offset);
            ++readCountInFile;
            read_tag.push_back(newtag);
            if (0 != newtag) tagged_reads += 2;
        }
    }
    else if (format==LinkedReadsFormat::seq){
        while (!feof(fd1) and !feof(fd2)) {
            //LinkedRead r1,r2;
            bsg10xTag newtag = 0;
            r1offset = ftell(fd1)+16;
            if (NULL == fgets(readbuffer, 999, fd1)) continue;
            //Tag to number from r1's name
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
            r2offset = ftell(fd2);
            if (NULL == fgets(readbuffer, 999, fd2)) continue;
            if (0 == readCountInFile % group_size) {
                g1offset = r1offset;
                group_offset1.push_back(r1offset);
                g2offset = r2offset;
                group_offset2.push_back(r2offset);
            }
            read_offset.push_back(r1offset - g1offset);
            read_offset.push_back(r2offset - g2offset);
            ++readCountInFile;
            read_tag.push_back(newtag);
            if (0 != newtag) tagged_reads += 2;
        }
    }
    sglib::OutputLog(sglib::LogLevels::INFO)<<"Datastore with "<<read_offset.size()-1<<" reads, "<<tagged_reads<<" reads with tags"<<std::endl; //and "<<reads_in_tag.size()<<"tags"<<std::endl;

}

void LinkedReadsDatastore::dump_index_to_disk(std::string filename) {
    std::ofstream f(filename);
    uint64_t s;
    s=filename1.size(); f.write((const char *) &s,sizeof(s));
    f<<filename1;
    s=filename2.size(); f.write((const char *) &s,sizeof(s));
    f<<filename2;
    s=group_size; f.write((const char *) &s,sizeof(s));
    s=group_offset1.size(); f.write((const char *) &s,sizeof(s));
    s=group_offset2.size(); f.write((const char *) &s,sizeof(s));
    s=read_tag.size(); f.write((const char *) &s,sizeof(s));
    s=read_offset.size(); f.write((const char *) &s,sizeof(s));
    f.write((const char *)group_offset1.data(),group_offset1.size()*sizeof(group_offset1[0]));
    f.write((const char *)group_offset2.data(),group_offset2.size()*sizeof(group_offset2[0]));
    f.write((const char *)read_tag.data(),read_tag.size()*sizeof(read_tag[0]));
    f.write((const char *)read_offset.data(),read_offset.size()*sizeof(read_offset[0]));
}

void LinkedReadsDatastore::load_index_from_disk(std::string filename) {
    std::ifstream f(filename);
    uint64_t s;
    f.read((char *) &s,sizeof(s)); filename1.resize(s);
    f.read((char *) filename1.data(),s*sizeof(filename1[0]));
    f.read((char *) &s,sizeof(s)); filename2.resize(s);
    f.read((char *) filename2.data(),s*sizeof(filename2[0]));
    f.read((char *) &s,sizeof(s)); group_size=s;
    f.read((char *) &s,sizeof(s)); group_offset1.resize(s);
    f.read((char *) &s,sizeof(s)); group_offset2.resize(s);
    f.read((char *) &s,sizeof(s)); read_tag.resize(s);
    f.read((char *) &s,sizeof(s)); read_offset.resize(s);
    f.read((char *)group_offset1.data(),group_offset1.size()*sizeof(group_offset1[0]));
    f.read((char *)group_offset2.data(),group_offset2.size()*sizeof(group_offset2[0]));
    f.read((char *)read_tag.data(),read_tag.size()*sizeof(read_tag[0]));
    f.read((char *)read_offset.data(),read_offset.size()*sizeof(read_offset[0]));
}


std::string LinkedReadsDatastore::get_read_sequence(size_t readID, FILE * file1, FILE * file2) {
#define SEQBUFFER_SIZE 260
    char buffer[SEQBUFFER_SIZE];
    if (0==readID%2){
        auto pos_in_file=readID/2-1;
//        std::cout<<"Read "<<readID<<" from file2, group "<<pos_in_file/group_size
//                 <<" group offset "<<group_offset2[pos_in_file/group_size]
//                 << " read_offset "<< read_offset[readID]<<std::endl;
        fseek(file2,group_offset2[pos_in_file/group_size]+read_offset[readID],SEEK_SET);
        fread(buffer,SEQBUFFER_SIZE,1,file2);

    } else {
        auto pos_in_file=readID/2;
//        std::cout<<"Read "<<readID<<" from file1, group "<<pos_in_file/group_size
//                 <<" group offset "<<group_offset1[pos_in_file/group_size]
//                 << " read_offset "<< read_offset[readID]<<std::endl;
        fseek(file1,group_offset1[pos_in_file/group_size]+read_offset[readID],SEEK_SET);
        fread(buffer,SEQBUFFER_SIZE,1,file1);
    }
    for (auto &c:buffer) if (c=='\n') {c='\0'; break;}
    return std::string(buffer);
}

std::string LinkedReadsDatastore::get_read_sequence_fd(size_t readID, int fd1, int fd2) {
#define SEQBUFFER_SIZE 260
    char buffer[SEQBUFFER_SIZE];
    if (0==readID%2){
        auto pos_in_file=readID/2-1;
//        std::cout<<"Read "<<readID<<" from file2, group "<<pos_in_file/group_size
//                 <<" group offset "<<group_offset2[pos_in_file/group_size]
//                 << " read_offset "<< read_offset[readID]<<std::endl;
        lseek(fd2,group_offset2[pos_in_file/group_size]+read_offset[readID],SEEK_SET);
        //fread(buffer,SEQBUFFER_SIZE,1,file2);
        read(fd2,buffer,SEQBUFFER_SIZE);

    } else {
        auto pos_in_file=readID/2;
//        std::cout<<"Read "<<readID<<" from file1, group "<<pos_in_file/group_size
//                 <<" group offset "<<group_offset1[pos_in_file/group_size]
//                 << " read_offset "<< read_offset[readID]<<std::endl;
        lseek(fd1,group_offset1[pos_in_file/group_size]+read_offset[readID],SEEK_SET);
        //fread(buffer,SEQBUFFER_SIZE,1,file1);
        read(fd1,buffer,SEQBUFFER_SIZE);
    }
    for (auto &c:buffer) if (c=='\n') {c='\0'; break;}
    return std::string(buffer);
}

void LinkedReadsDatastore::get_read_sequence_fd(size_t readID, int fd1, int fd2, char * dest) {
#define SEQBUFFER_SIZE 260
    if (0==readID%2){
        auto pos_in_file=readID/2-1;
        lseek(fd2,group_offset2[pos_in_file/group_size]+read_offset[readID],SEEK_SET);
        read(fd2,dest,SEQBUFFER_SIZE);

    } else {
        auto pos_in_file=readID/2;
        lseek(fd1,group_offset1[pos_in_file/group_size]+read_offset[readID],SEEK_SET);
        read(fd1,dest,SEQBUFFER_SIZE);
    }
}

bsg10xTag LinkedReadsDatastore::get_read_tag(size_t readID) {
    return read_tag[(readID-1)/2];
}

const char* BufferedLRSequenceGetter::get_read_sequence(uint64_t readID) {
    if (0==readID%2){
        auto pos_in_file=readID/2-1;
        size_t read_offset_in_file;
        read_offset_in_file=datastore.group_offset2[pos_in_file/datastore.group_size]+datastore.read_offset[readID];
        if (read_offset_in_file<buffer2_offset or read_offset_in_file+chunk_size>buffer2_offset+bufsize) {
            buffer2_offset=read_offset_in_file;
            lseek(fd2,read_offset_in_file,SEEK_SET);
            read(fd2,buffer2,bufsize);
        }
        return buffer2+(read_offset_in_file-buffer2_offset);

    } else {
        auto pos_in_file=readID/2;
        size_t read_offset_in_file;
        read_offset_in_file=datastore.group_offset1[pos_in_file/datastore.group_size]+datastore.read_offset[readID];
        if (read_offset_in_file<buffer1_offset or read_offset_in_file+chunk_size>buffer1_offset+bufsize) {
            buffer1_offset=read_offset_in_file;
            lseek(fd1,read_offset_in_file,SEEK_SET);
            read(fd1,buffer1,bufsize);
        }
        return buffer1+(read_offset_in_file-buffer1_offset);
    }
}

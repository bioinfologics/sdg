//
// Created by Luis Yanes (EI) on 14/09/2017.
//

#ifndef SEQSORTER_SMR_H
#define SEQSORTER_SMR_H

#include <fcntl.h>
#include <unistd.h>

#include <atomic>
#include <vector>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <string>
#include <iterator>

#include <cstdio>
#include <cstdlib>
#include <set>
#include <unordered_set>
#include <cmath>
#include <sglib/filesystem/check_or_create_directory.h>


static const std::uint64_t GB(1024*1024*1024);

/**
 * @brief
 * SMR is an external memory map reduce engine with configurable memory bounds
 * it requires a RecordFactory and a FileReader. The FileReader takes elements from the file
 * and sends them to the RecordFactory for parsing into single RecordType elements.
 *
 * @tparam RecordType
 * The RecordType class requires strict weak ordering and equivalence comparison functions.
 * @tparam RecordFactory
 * Generates a vector of RecordType elements to be merged into a single set.
 * @tparam FileReader
 * Gets a single FileRecord from the input file.
 * @tparam FileRecord
 * FileRecord containing zero or more RecordType elements.
 * @tparam ReaderParamStruct
 * Parameters for the FileReader constructor.
 * @tparam FactoryParamStruct
 * Parameters for the RecordFactory constructor.
 */
template <class RecordType, class RecordFactory, class FileReader, typename FileRecord, typename ReaderParamStruct, typename FactoryParamStruct >
class SMR {
public:

    /***
     * @brief
     * @param reader_parameters
     * The arguments to the FileReader, generally initialised using list initialization aka the {} initialization.
     * @param factory_parameters
     * The arguments to the FactoryParamStruct, generally initialised using list initialization aka the {} initialization.
     * @param maxMem
     * Maximum amount of memory available for the SMR in bytes, generally the units are multiplied on call.
     * I.E 4*GB, where GB is defined as (1024*1024*1024) on this header.
     * @param min
     * This refers to the *minimum* filter for the number of times a RecordType element has to be seen to be on the output.
     * @param max
     * This refers to the *maximum* filter for the number of times a RecordType element has to be seen to be on the output.
     * @param outdir
     * In order to reuse the results of the SMR, an output directory can be provided to store the final.kc file,
     * generally this directory is associated with the base working directory. The SMR will then compound the outdir with
     * another directory level "smr_files" and the "read_file" basename of the string passed to the read_from_file function call.
     * @param Otmp
     * To be able to execute multiple SMRs at the same time for different inputs, the Otmp parameter stores the path of the
     * directory to be used for the batch files, each input file will generate a new directory from the basename of the argument
     * passed to the read_from_file function call (read_file). The behaviour is equivalent to that of outdir.
     */
    SMR(ReaderParamStruct reader_parameters, FactoryParamStruct factory_parameters, uint64_t maxMem, unsigned int min = 0, unsigned int max = std::numeric_limits<unsigned int>::max(),
        const std::string outdir = "./", const std::string &Otmp = "./") : reader_parameters(reader_parameters),
                                                                        factory_parameters(factory_parameters),
                                                                        myBatches(0),
                                                                        totalFilteredRecords(0),
                                                                        totalRecordsGenerated(0), tmpBase(Otmp),
                                                                        outdir(outdir), minCount(min), maxCount(max),
                                                                        mergeCount(4), maxThreads((unsigned int) 1) {
        this->outdir = std::string(outdir+"smr_files/");
        sglib::check_or_create_directory(this->outdir);
        sglib::check_or_create_directory(this->tmpBase);
        numElementsPerBatch = (maxMem / sizeof(RecordType) / maxThreads);
        std::cout<<"SMR created with elements of size "<<sizeof(RecordType)<<" using "<<maxThreads<<" threads and "
                 <<maxMem<<" Bytes of memory -> batches of "<<numElementsPerBatch<<" elements"<<std::endl;

    }

    /**
     * @brief
     * Transforms the input file into a filtered RecordType vector using maxMem and merging subsets once every mergeCount times.
     * Removes the tmp directory used for the external memory reduction.
     * @param read_file
     * Path to the file to be read
     * @return
     * Filtered vector of RecordType elements
     */
    std::vector<RecordType> read_from_file(const std::string &read_file) {
        std::string readFileBasename(read_file.substr(0,read_file.find_last_of('.')).substr(read_file.find_last_of('/')+1));
        std::string finalFilePath(outdir+readFileBasename);
        tmpInstance = tmpBase+readFileBasename;
        sglib::check_or_create_directory(finalFilePath);
        outdir = finalFilePath;
        sglib::check_or_create_directory(tmpInstance);
        std::ifstream final_file(finalFilePath+"final.kc");
        if (final_file.is_open()){
            std::cout << "Using precomputed sum file at " << outdir << "final.kc" << std::endl;
            return readFinalkc(outdir+"final.kc");
        } else {
            uint64_t numFileRecords(0);

            std::chrono::time_point<std::chrono::system_clock> start, end;
            start = std::chrono::system_clock::now();

            std::cout << "Reading file: " << read_file << std::endl;
            FileReader myFileReader(reader_parameters, read_file);
            std::cout << "Begin reduction using " << numElementsPerBatch << " elements per batch (" << ceil(uint64_t((numElementsPerBatch*sizeof(RecordType)*maxThreads)) / (1.0f*1024*1024*1024)) << "GB)" << std::endl;
            mapElementsToBatches(myFileReader, numFileRecords);

            readerStatistics = myFileReader.getSummaryStatistics();

            end = std::chrono::system_clock::now();
            std::chrono::duration<double> elapsed_seconds = end - start;
            std::cout << "Done reduction in " << elapsed_seconds.count() << "s" << std::endl;
            sglib::remove_directory(tmpInstance);
            return getRecords();
        }
    };

    /**
     * @brief
     * Transforms the input from memory into a filtered RecordType vector using maxMem and merging subsets once every mergeCount times.
     * Removes the tmp directory used for the external memory reduction.
     * @param read_file
     * Path to the file to be read
     * @return
     * Filtered vector of RecordType elements
     */
    std::vector<RecordType> process_from_memory() {
        std::string finalFilePath(outdir+"SMR_mem");
        tmpInstance = tmpBase+"SMR_mem";
        sglib::check_or_create_directory(finalFilePath);
        outdir = finalFilePath;
        sglib::check_or_create_directory(tmpInstance);
        std::ifstream final_file(finalFilePath+"final.kc");

        uint64_t numFileRecords(0);

        std::chrono::time_point<std::chrono::system_clock> start, end;
        start = std::chrono::system_clock::now();

        std::cout << "Reading From memory"<< std::endl;
        FileReader myFileReader(reader_parameters);
        std::cout << "Begin reduction using " << numElementsPerBatch << " elements per batch (" << ceil(uint64_t((numElementsPerBatch*sizeof(RecordType)*maxThreads)) / (1.0f*1024*1024*1024)) << "GB)" << std::endl;
        mapElementsToBatches(myFileReader, numFileRecords);

        readerStatistics = myFileReader.getSummaryStatistics();

        end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = end - start;
        std::cout << "Done reduction in " << elapsed_seconds.count() << "s" << std::endl;
        sglib::remove_directory(tmpInstance);
        //TODO: remove instance / never create the final files?
        return getRecords();

    };

    /**
     * @brief
     * Provides summary statistics regarding the execution of the "read_from_file" function
     * @return
     * Returns a vector containing: total records generated by the factory,
     * total records after filtering,
     * total file records present in the file and the total number of records after filtering.
     */
    std::vector<uint64_t> summaryStatistics() {
        return {totalRecordsGenerated, totalFilteredRecords, readerStatistics.totalRecords, readerStatistics.filteredRecords};
    }

private:

    /***
     * @brief
     * Merges memory and disk batches. No filtering is applied at this stage.
     * @param tmpName
     * Name of the temporary output file to hold the batch
     * @param files
     * Vector of file names holding the disk batches to merge
     * @param elements
     * In memory batch of elements (represented as a vector)
     * @return
     * This function returns the number of elements merged
     */
    uint64_t merge(const std::string &tmpName, const std::vector<std::string> &files, std::vector<RecordType> &elements) {
        auto files_size(files.size());
        auto numMemoryElements(elements.size());
        uint64_t memoryElement(0);
        std::vector<int> in_fds(files_size );
        int out = ::open(tmpName.c_str(), O_CREAT|O_WRONLY|O_TRUNC, 0777);

        std::vector<bool> active_file(files_size,true);

        std::vector<RecordType*> next_element_from_file(files.size());
        std::vector<unsigned int> count_element_from_file(files_size,0);
        std::vector<unsigned int> size_element_from_file(files_size,0);
        std::vector<unsigned int> seen_element_from_file(files_size,0);
        uint64_t outCount(0), count0(0), bufferSize(10000000);
        uint64_t outBufferSize(50000000);

        ::write(out, &outCount, sizeof(outCount));

        auto *outfreqs=(RecordType*) malloc(outBufferSize* sizeof(RecordType));

        for (auto i = 0; i < files_size; i++) {
            in_fds[i] = ::open(files[i].c_str(), O_RDONLY);
            uint64_t size;
            ::read(in_fds[i], &size, sizeof(size));
            //std::cout << files[i] << " size " << size << " in disk " << sizeof(size) + size* sizeof(RecordType) << std::endl;
            next_element_from_file[i] = (RecordType*) malloc(bufferSize* sizeof(RecordType));
            count_element_from_file[i] = 0;
            size_element_from_file[i]=
                    ::read(in_fds[i], next_element_from_file[i], bufferSize* sizeof(RecordType)) / sizeof(RecordType);
            active_file[i] = true;
        }

        RecordType min_element;
        unsigned int min = 0;
        for (unsigned int i = 1; i < files_size; ++i) {
            if (next_element_from_file[i] < next_element_from_file[i-1]) min = i;
        }

        RecordType current_element;
        if (elements[memoryElement] < next_element_from_file[min][0]) {
            current_element = elements[memoryElement];
        } else {
            current_element = next_element_from_file[min][count_element_from_file[min]];
            count_element_from_file[min];
        }
        current_element.count=0;


        bool active;
        do {
            active=false;
            min_element = RecordType();

            for (int i=0; i<files_size; ++i) {
                if (size_element_from_file[i] <= count_element_from_file[i]) continue;

                if (current_element == next_element_from_file[i][count_element_from_file[i]]) {
                    current_element.merge(next_element_from_file[i][count_element_from_file[i]]);
                    count_element_from_file[i]++;
                    if (count_element_from_file[i] == size_element_from_file[i]) {
                        size_element_from_file[i] =
                                ::read(in_fds[i], (char *) next_element_from_file[i], bufferSize*sizeof(RecordType)) / sizeof(RecordType);
                        seen_element_from_file[i]+=count_element_from_file[i];
                        count_element_from_file[i] = 0;
                    }
                }
                if (count_element_from_file[i] < size_element_from_file[i]) {
                    if (!(next_element_from_file[i][count_element_from_file[i]] > min_element)){
                        active=true;
                        min_element = next_element_from_file[i][count_element_from_file[i]];
                    }
                }
            }
            if (memoryElement < numMemoryElements) {
                if (current_element == elements[memoryElement]) {
                    current_element.merge(elements[memoryElement]);
                    memoryElement++;
                }
            }
            if (memoryElement<numMemoryElements and !(elements[memoryElement] > min_element)) {
                active=true;
                min_element = elements[memoryElement];
            }

            outfreqs[count0] = current_element;
            ++count0;
            if (count0 == outBufferSize) {
                ::write(out, (char *) outfreqs, outBufferSize * sizeof(RecordType));
                outCount += count0;
                count0 = 0;
            }
            current_element = min_element;
            current_element.count=0;
        } while ( active );

        if (count0 > 0) {
            ::write(out, (char *) outfreqs, count0 * sizeof(RecordType));
            outCount += count0;
        }
        lseek(out, 0, SEEK_SET);
        ::write(out, (char *) &outCount, sizeof(outCount));
        free(outfreqs);
        close(out);

        for (auto i = 0; i < files_size; i++) {
            free(next_element_from_file[i]);
            std::remove(files[i].c_str());
        }
        return outCount;
    }

    /***
     * @brief
     * Merges disk batches and applies a min,max filter for the element counts.
     * Removes all the temporary files once it has finished using them.
     * @param tmpName
     * Temporary name of the output file
     * @param files
     * Vector of file names to merge
     * @return
     * Number of elements merged and filtered
     */
    uint64_t merge(const std::string &tmpName, const std::vector<std::string> &files) {
        auto files_size(files.size());
        std::vector<int> in_fds( files_size );
        int out = ::open(tmpName.c_str(), O_CREAT|O_WRONLY|O_TRUNC, 0777);
        if (out < 0) {
            perror(tmpName.c_str());
        }
        std::vector<bool> active_file(files_size,true);

        std::vector<RecordType*> next_element_from_file(files.size());
        std::vector<unsigned int> count_element_from_file(files_size,0);
        std::vector<unsigned int> size_element_from_file(files_size,0);
        std::vector<unsigned int> seen_element_from_file(files_size,0);
        uint64_t outCount(0), count0(0), bufferSize(10000000);
        uint64_t outBufferSize(50000000);

        ::write(out, &outCount, sizeof(outCount));

        auto *outfreqs=(RecordType*) malloc(outBufferSize* sizeof(RecordType));

        for (auto i = 0; i < files_size; i++) {
            in_fds[i] = ::open(files[i].c_str(), O_RDONLY);
            uint64_t size;
            ::read(in_fds[i], &size, sizeof(size));
            //std::cout << files[i] << " size " << size << " in disk " << sizeof(size) + size* sizeof(RecordType) << std::endl;
            next_element_from_file[i] = (RecordType*) malloc(bufferSize* sizeof(RecordType));
            count_element_from_file[i] = 0;
            size_element_from_file[i]=
                    ::read(in_fds[i], next_element_from_file[i], bufferSize* sizeof(RecordType)) / sizeof(RecordType);
            active_file[i] = true;
        }

        unsigned int min = 0;
        for (unsigned int i = 1; i < files_size; ++i) {
            if (next_element_from_file[i] < next_element_from_file[i-1]) min = i;
        }

        RecordType current_element(next_element_from_file[min][count_element_from_file[min]]);
        current_element.count=0;

        RecordType min_element;
        bool active;
        do {
            active=false;
            min_element = RecordType();

            for (int i=0; i<files_size; ++i) {
                if (size_element_from_file[i] <= count_element_from_file[i]) continue;

                if (current_element == next_element_from_file[i][count_element_from_file[i]]) {
                    current_element.merge(next_element_from_file[i][count_element_from_file[i]]);
                    count_element_from_file[i]++;
                    if (count_element_from_file[i] == size_element_from_file[i]) {
                        size_element_from_file[i] =
                                ::read(in_fds[i], (char *) next_element_from_file[i], bufferSize*sizeof(RecordType)) / sizeof(RecordType);
                        seen_element_from_file[i]+=count_element_from_file[i];
                        count_element_from_file[i] = 0;
                    }
                }
                if (count_element_from_file[i] < size_element_from_file[i]) {
                    if (!(next_element_from_file[i][count_element_from_file[i]] > min_element)){
                        active=true;
                        min_element = next_element_from_file[i][count_element_from_file[i]];
                    }
                }
            }

            if (minCount < current_element.count and current_element.count <= maxCount) {
                outfreqs[count0] = current_element;
                ++count0;
            }
            if (count0 == outBufferSize) {
                ::write(out, (char *) outfreqs, outBufferSize * sizeof(RecordType));
                outCount += count0;
                count0 = 0;
            }
            current_element = min_element;
            current_element.count=0;
        } while ( active );

        if (count0 > 0) {
            ::write(out, (char *) outfreqs, count0 * sizeof(RecordType));
            outCount += count0;
        }
        lseek(out, 0, SEEK_SET);
        ::write(out, (char *) &outCount, sizeof(outCount));
        free(outfreqs);
        close(out);

        for (auto i = 0; i < files_size; i++) {
            free(next_element_from_file[i]);
            std::remove(files[i].c_str());
        }
        return outCount;
    }

    /**
     * @brief
     * Loops over the records on myFileReader and generates a vector<RecordType> which is then appended to the _elements
     * vector, when the size of that vector is >= to the limit, a batch is sent to disk. If the number of batches on
     * the disk equals the mergeCount, a merge is done for the records in memory and the records on the disk.
     * This structure creates mergeCount batches and then reduces to one batch after the merge.
     * @param myFileReader
     * This is an already initialised FileReader class.
     * @param numReadsReduction
     * The number of file records generated by myFileReader.
     */
    void mapElementsToBatches(FileReader &myFileReader, uint64_t &numReadsReduction) {
        std::vector<RecordType> _elements;
        unsigned int totalBatches(0);
        unsigned int currentBatch(0);
        _elements.reserve(numElementsPerBatch);
        RecordFactory myRecordFactory(factory_parameters);
        FileRecord frecord;
        while (myFileReader.next_record(frecord)) {
            numReadsReduction++;
            std::vector<RecordType> record;
            myRecordFactory.setFileRecord(frecord);
            myRecordFactory.next_element(record);
            //while (!record.empty())
            {
                _elements.insert(_elements.end(), record.cbegin(), record.cend());
                if (_elements.size() >= numElementsPerBatch) {
                    totalBatches++;
                    totalRecordsGenerated+=_elements.size();
                    if (currentBatch==mergeCount) {
                        mergeBatches(currentBatch, _elements);
                        _elements.clear();
                        currentBatch=1;
                    } else {
                        dumpBatch(currentBatch, _elements);
                        _elements.clear();
                        currentBatch++;
                    }
                }
            }
        }
        totalRecordsGenerated+=_elements.size();

        if (totalBatches) mergeBatches(currentBatch, _elements);
        else {
            dumpBatch(currentBatch, _elements);
        }
        _elements.clear();
        rename(std::string(tmpInstance + "thread_0_batch_0.tmc").c_str(), std::string(tmpInstance + "thread_0.tmc").c_str());
    }

    /**
     * @brief
     * Collapses the _elements and writes the collapsed elements to the currentBatch file.
     * @param currentBatch
     * Number which is used on the file name to indicate batch count.
     * @param _elements
     * A list of non-collapsed elements to dump to disk.
     */
    void dumpBatch(const unsigned int currentBatch, std::vector<RecordType> &_elements) {
        collapse(_elements);
        // Simply dump the file
        std::string outBatchName(std::string(tmpInstance + "thread_0" + "_batch_" + std::to_string(currentBatch) + ".tmc"));
        // Simply dump the file
        std::ofstream outBatch(outBatchName.data(), std::ios_base::binary | std::ios_base::out | std::ios_base::trunc);
        auto size = _elements.size();
        outBatch.write((char *)&size, sizeof(size));
        outBatch.write((char *)_elements.data(), size*sizeof(RecordType));
        outBatch.close();
        _elements.clear();
    }

    void mergeBatches(unsigned int currentCount, std::vector<RecordType> &_elements) {
        collapse(_elements);
        unsigned int numFilesToMerge = currentCount;

        std::vector<std::string> inBatches(numFilesToMerge);
        for (auto batch=0; batch < numFilesToMerge; batch++) {
            std::string filename(tmpInstance+ "thread_0" + "_batch_" + std::to_string(batch)+ ".tmc");
            inBatches[batch] = filename;
        }
        std::string threadFile(tmpInstance+ "thread_0" + "_tmpBatch1.tmc");
        merge(threadFile, inBatches, _elements);
        rename(std::string(tmpInstance+ "thread_0" + "_tmpBatch1.tmc").c_str(), std::string(tmpInstance+ "thread_0" + "_batch_0.tmc").c_str());
    }

    /**
     * @brief
     * Sorts the _elements and then performs an inplace collapse of all equal elements using the merge function provided
     * by the RecordType, the _elements vector will be re-sized as to only contain the collapsed elements.
     * @param _elements
     * In/out parameter of sorted elements which are collapsed and resized.
     */
    void collapse(std::vector<RecordType> &_elements){
        if (_elements.size() == 0) return;
        std::sort(_elements.begin(), _elements.end());
        typename std::vector<RecordType>::iterator writePtr;
        writePtr = _elements.begin();
        typename std::vector<RecordType>::const_iterator endPtr;
        endPtr = _elements.end();

        // Keep the total in the first equal KMer
        typename std::vector<RecordType>::const_iterator readPtr;
        for (readPtr = _elements.begin(); readPtr != endPtr; ++writePtr) {
            *writePtr = *readPtr;
            ++readPtr;
            while (readPtr != endPtr and *writePtr == *readPtr) {
                writePtr->merge(*readPtr);
                ++readPtr;
            }
        }
        // After accumulating the values on the first KMer, remove all duplicates using unique
        // this operation is safe because the first copy is kept for all unique values
        // resize the object by erasing the duplicates
        _elements.resize(std::distance(_elements.begin(), writePtr));
    }


    std::vector<RecordType> readFinalkc(std::string finalFile) {
        std::ifstream outf(finalFile, std::ios_base::binary);

        outf.unsetf(std::ios::skipws);

        std::streampos fileSize;
        outf.seekg(0, std::ios::end);
        fileSize = outf.tellg();
        outf.seekg(sizeof(uint64_t), std::ios::beg);
        // reserve capacity
        std::vector<RecordType> vec;
        vec.reserve(fileSize);

        // read the data:
        vec.insert(vec.begin(), std::istream_iterator<RecordType>(outf), std::istream_iterator<RecordType>());

        return vec;

    }
    /**
     * @brief
     * Does a final merge of all thread's data and returns a collapsed vector along with writing the results to the final.kc file.
     * @return
     * The final vector of collapsed elements
     */
    std::vector<RecordType> getRecords() {

        std::vector<std::string> threadFiles(maxThreads);

        for (auto threadID = 0; threadID < maxThreads; threadID++) {
            std::string name(tmpInstance + "thread_" + std::to_string(threadID) + ".tmc");

            threadFiles[threadID] = name;
        }

        totalFilteredRecords = merge(outdir + "final.kc", threadFiles);

        return readFinalkc(outdir+"final.kc");
    }

    unsigned int minCount;                          /// Minimum times a record has to be seen to be accepted
    unsigned int maxCount;                          /// Maximum times a record has to been seen before being denied
    FactoryParamStruct factory_parameters;          /// Parameters for the RecordFactory class
    ReaderParamStruct reader_parameters;            /// Parameters for the FileReader class
    uint64_t numElementsPerBatch;                   /// Number of elements per batch (calculated using memory limit)
    std::atomic<uint64_t> myBatches;                /// Batches status variable

    ReaderStats readerStatistics;

    std::atomic<uint64_t> totalRecordsGenerated;    /// Total number of records
    uint64_t totalFilteredRecords;                  /// Total number of resulting elements after filtering and merging
    const int unsigned maxThreads;
    int mergeCount;                                 /// How many batches to keep rolling before merging
    std::string tmpBase;                            /// Directory to store the temporary files
    std::string tmpInstance;                        /// Directory to store the temporary files for the read_from_file argument
    std::string outdir;                             /// Output directory
};

#endif //SEQSORTER_SMR_H

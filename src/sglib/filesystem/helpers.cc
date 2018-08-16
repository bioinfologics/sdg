//
// Created by Luis Yanes (EI) on 23/11/2017.
//

#include <sys/param.h>
#include <cstring>
#include <array>
#include <fcntl.h>
#include <cassert>
#include "helpers.h"

bool sglib::check_file(std::string &filepath) {
    struct stat sb{};
    bool validate_dir(false);
    if (stat(filepath.c_str(), &sb) != 0) {
        if (stat(filepath.c_str(), &sb) != 0) {
            perror(filepath.c_str());
            validate_dir = false;
        }
    } else if (!S_ISDIR(sb.st_mode)) {
        validate_dir = true;
    } else {
        std::cout << filepath << " is not a file " << std::endl;
        validate_dir = false;
    }

    return validate_dir;
}
bool sglib::check_or_create_directory(std::string &output_prefix) {
        if (output_prefix.back() != '/') {
            output_prefix.push_back('/');
        }
        struct stat sb{};
        bool validate_dir(false);
        if (stat(output_prefix.c_str(), &sb) != 0) {
            if (errno == ENOENT) {
                mode_t mask = umask(0);
                umask(mask);
                std::cout<<"Creating: " << output_prefix << std::endl;
                mkdir(output_prefix.c_str(), mode_t(0777 - mask));
                validate_dir = true;
            }
            if (stat(output_prefix.c_str(), &sb) != 0) {
                perror(output_prefix.c_str());
                validate_dir = false;
            }
        } else if (!S_ISDIR(sb.st_mode)) {
            std::cout << output_prefix << " is not a directory " << std::endl;
        } else {
            validate_dir = true;
        }
        return validate_dir;
    }
void sglib::remove_directory(std::string path) {
    std::cout << "Removing: " << path << std::endl;
    ::rmdir(path.c_str());
}

std::string sglib::create_temp_directory(std::string prefix = "/tmp") {
    std::string tmplt (
            (prefix.empty()) ?
                std::string("/tmp/smr-tmp-XXXXXX").c_str() :
                std::string(prefix + "/smr-tmp-XXXXXX").c_str()
    );

    char buffer[MAXPATHLEN] = {0};
    strncpy(buffer, tmplt.c_str(), tmplt.size());
    auto result = mkdtemp(buffer);
    if (result ==  nullptr) {
        std::cerr << "Can't create " << buffer
                  << ", reason: " << strerror(errno) << "\n";
    } else {
        std::cout << "Created " << result << "\n";
    }

    return std::string(buffer);
}

bool sglib::copy_file(const std::string& from_p, const std::string& to_p, bool fail_if_exists)
{
    const std::size_t buf_sz = 32768;
    std::array<char, buf_sz> buf{};
    int infile=-1, outfile=-1;  // -1 means not open

    // bug fixed: code previously did a stat()on the from_file first, but that
    // introduced a gratuitous race condition; the stat()is now done after the open()

    if ((infile = ::open(from_p.c_str(), O_RDONLY))< 0)
    { return false; }

    struct stat from_stat{};
    if (::stat(from_p.c_str(), &from_stat)!= 0)
    {
        ::close(infile);
        return false;
    }

    int oflag = O_CREAT | O_WRONLY | O_TRUNC;
    if (fail_if_exists)
        oflag |= O_EXCL;
    if ((outfile = ::open(to_p.c_str(), oflag, from_stat.st_mode))< 0)
    {
        int open_errno = errno;
        assert(infile >= 0);
        ::close(infile);
        errno = open_errno;
        return false;
    }

    ssize_t sz, sz_read=1, sz_write;
    while (sz_read > 0
           && (sz_read = ::read(infile, buf.data(), buf_sz)) > 0)
    {
        // Allow for partial writes - see Advanced Unix Programming (2nd Ed.),
        // Marc Rochkind, Addison-Wesley, 2004, page 94
        sz_write = 0;
        do
        {
            assert(sz_read - sz_write > 0);  // #1
            // Possible infinite loop if write returns 0. My analysis
            // is that POSIX specifies 0 return only if 3rd arg is 0, and that will never
            // happen due to loop entry and continuation conditions. assert #1 above
            // and #2 below added to verify that analysis.
            if ((sz = ::write(outfile, buf.data() + sz_write,
                              sz_read - sz_write)) < 0)
            {
                sz_read = sz; // cause read loop termination
                break;        //  and error reported after closes
            }
            assert(sz > 0);                  // #2
            sz_write += sz;
        } while (sz_write < sz_read);
    }

    if (::close(infile)< 0)
        sz_read = -1;
    if (::close(outfile)< 0)
        sz_read = -1;

    return sz_read >= 0;
}
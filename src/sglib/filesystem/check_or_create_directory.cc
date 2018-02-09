//
// Created by Luis Yanes (EI) on 23/11/2017.
//

#include <sys/param.h>
#include <string.h>
#include "check_or_create_directory.h"

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
    const char * const tmplt (
            (prefix.empty()) ?
            std::string("/tmp/smr-tmp-XXXXXX").c_str() :
            std::string(prefix + "/smr-tmp-XXXXXX").c_str()
    );

    char buffer[MAXPATHLEN] = {0};
    strncpy(buffer, tmplt, strlen(tmplt));
    auto result = mkdtemp(buffer);
    if (result ==  nullptr) {
        std::cerr << "Can't create " << buffer
                  << ", reason: " << strerror(errno) << "\n";
    } else {
        std::cout << "Created " << result << "\n";
    }

    return std::string(buffer);
}
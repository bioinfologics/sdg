//
// Created by Luis Yanes (EI) on 23/11/2017.
//

#include "check_or_create_directory.h"
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
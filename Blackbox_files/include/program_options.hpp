#ifndef THAT_THIS_PROGRAM_OPTIONS_HEADER_FILE_IS_ALREADY_INCLUDED
#define THAT_THIS_PROGRAM_OPTIONS_HEADER_FILE_IS_ALREADY_INCLUDED

#include <vector>

namespace program_options {

    void parse(int argc, char* argv[]);

    const std::vector<double>&  variables();
    const std::vector<double>&  parameters();
    bool                        verbose();
    bool                        return_gradients();

};  // namespace program_options

#endif  // THAT_THIS_PROGRAM_OPTIONS_HEADER_FILE_IS_ALREADY_INCLUDED
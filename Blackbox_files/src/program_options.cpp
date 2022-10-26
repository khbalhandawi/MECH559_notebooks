#include "program_options.hpp"

#include <stdexcept>
#include <vector>
#include <string>

namespace {

static std::vector<double>  _variables;
static std::vector<double>  _parameters;
static bool                 _verbose = false;
static bool                 _return_gradients = false;
static std::string          _current_arg = "";

}  // namespace

void program_options::parse(int argc, char* argv[]) {
    const std::vector<std::string> args(argv + 1, argv + argc);

    for (const auto& arg : args) {

        if (arg == "-v" || arg == "--verbose") {
            _current_arg = arg;
            if (_verbose) {
                throw std::runtime_error("cannot use -v/--verbose parameter twice!");
            }
            _verbose = true;
            continue;
        }

        if (arg == "-G" || arg == "--gradients") {
            _current_arg = arg;
            if (_return_gradients) {
                throw std::runtime_error("cannot use -G/--gradients parameter twice!");
            }
            _return_gradients = true;
            continue;
        }

        if (_variables.empty()) {
            if (arg == "-x" || arg == "--variables") {
                _current_arg = arg;
                continue;
            }

        } else {
            if (arg == "-x" || arg == "--variables") {
                throw std::runtime_error("cannot use -x/--variables parameter twice!");
                continue;
            }
        }

        if (_parameters.empty()) {
            if (arg == "-p" || arg == "--parameters") {
                _current_arg = arg;
                continue;
            }

        } else {
            if (arg == "-p" || arg == "--parameters") {
                throw std::runtime_error("cannot use -p/--parameters parameter twice!");
                continue;
            }
        }

        if (_current_arg == "-x" || _current_arg == "--variables") {
            _variables.push_back(stod(arg));
            continue;
        }

        if (_current_arg == "-p" || _current_arg == "--parameters") {
            _parameters.push_back(stod(arg));
            continue;
        }

    }

    // default arguments
    if (_parameters.empty()) {
        _parameters = {8,-4,-3,-3,-3,3,2,0.8};
    } else {
        if (_parameters.size() != 8) {
            std::string error_msg = "exactly 8 parameters are needed! " + std::to_string(_parameters.size()) + " parameters are provided!";
            throw std::runtime_error(error_msg);
        }
    }

    // default arguments
    if (_variables.empty()) {
        std::string error_msg = "no variables x provided! Use the -x argument to provide them";
        throw std::runtime_error(error_msg);
    } else {
        if (_variables.size() != 2) {
            std::string error_msg = "exactly 2 variables are needed! " + std::to_string(_variables.size()) + " parameters are provided!";
            throw std::runtime_error(error_msg);
        }
    }
}

const std::vector<double>&
program_options::variables() {
    return _variables;
}

const std::vector<double>&
program_options::parameters() {
    return _parameters;
}

bool
program_options::verbose() {
    return _verbose;
}

bool
program_options::return_gradients() {
    return _return_gradients;
}

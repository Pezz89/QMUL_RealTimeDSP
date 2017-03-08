// Used for transfering CLI arguments from main to render function. Only way I
// could think of to dereference a void* used for command line arguments in
// render.cpp
struct UserOpts {
    bool bonus;
};

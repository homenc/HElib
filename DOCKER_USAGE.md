# Dockerized usage
To create an Ubuntu 18.04-based docker image with HElib and its dependencies
installed, use the following steps:

1. Clone HElib and `cd` into it.
2. Run `docker build .`, optionally tagging with `-t helib` or similar.  This
   will create an image with HElib installed globally in `/usr/local`.  Most
   build systems and compilers will look in here.  Furthermore, `cmake` can be
   used with `find_package(helib)` and invoked with no special command-line
   arguments.
3. The image is now ready to be used.  Some possibilities include:
   - Create a Dockerfile which begins with `FROM helib`, where `helib` may be
   any tag you used in step 2, and include your application logic in this
   Dockerfile.
   - Run a container from the HElib docker image directly to experiment, with
   `docker run -it helib bash`

FROM alpine:3.8

# Location of the src in the contrainer
ARG helibDir=/helib
# Location of the build in the contrainer
ARG buildDir=/build

COPY . ${helibDir}

# Currently fails when cmake is invoked
RUN apk add --update \
    g++              \
    make             \
    cmake            \
    libtool          \
  && rm -rf /var/cache/apk/* \
  && mkdir ${buildDir}       \
  && cd ${buildDir}          \
  && cmake ${helibDir}       \
  && make          

CMD [ "" ]


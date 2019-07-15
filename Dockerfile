FROM alpine:3.8

# Location of the src in the contrainer
ARG helibDir=/helib
# Location of the build in the contrainer
ARG buildDir=/build

# Currently fails when cmake is invoked
RUN apk add --update \
    g++              \
    make             \
    cmake            \
    libtool          \
    perl
RUN apk add gmp-dev \
    && apk add ca-certificates wget \
    && apk update
    
RUN wget https://www.shoup.net/ntl/ntl-11.3.2.tar.gz
RUN tar xf ntl-11.3.2.tar.gz \
    && cd ntl-11.3.2/src \
    && ./configure \
    && make \
    && make install

COPY . ${helibDir}

RUN  rm -rf /var/cache/apk/* \
  && mkdir ${buildDir}       \
  && cd ${buildDir}          \
  && cmake ${helibDir}       \
  && make          

CMD [ "" ]


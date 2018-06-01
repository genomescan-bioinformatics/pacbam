FROM debian:jessie-slim AS build

# ARG pkg_version

WORKDIR /build

RUN apt-get update                                               && \
    apt-get install -y gcc                                          \
                       g++                                          \
                       make                                         \
                       git                                          \
                       zlib1g-dev                                   \
                                                                 && \
    git clone https://CibioBCG@bitbucket.org/CibioBCG/pacbam.git && \
    cd pacbam                                                    && \
    make -f Makefile.linux

FROM debian:jessie-slim


RUN apt-get update                       && \
    apt-get install -y python2.7            \
                       python-matplotlib    \
                       python-numpy      && \
    rm -rf /var/lib/apt/lists/*

COPY --from=build /build/pacbam/pacbam /usr/bin/
COPY --from=build /build/pacbam/pacbam_report.py /usr/bin/

CMD [ "/usr/bin/pacbam" ]

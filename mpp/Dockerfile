FROM harbor.cluster.math.kit.edu/mpp_release/mpp_dev_ubuntu:22.04

ARG USER=mpi
ARG WORKDIR=./
ARG UPSTREAM_COMMIT=""
ARG CMAKE_ARGS=""

WORKDIR ${WORKDIR}
USER ${USER}

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV WORKDIR ${WORKDIR}
ENV USER ${USER}

COPY --chown=${USER}:${USER} . ${WORKDIR}

RUN echo ${UPSTREAM_COMMIT}

# Only needed if mpp is used as submodule
RUN if [ "${UPSTREAM_COMMIT}" = "" ]; then \
      echo "No UPSTREAM_COMMIT"; \
    else \
      cd mpp && \
      git config --get remote.origin.fetch && \
      git config remote.origin.fetch "+refs/heads/*:refs/remotes/origin/*" && \
      git config --get remote.origin.fetch && \
      git fetch --all && \
      git checkout ${UPSTREAM_COMMIT} && \
      cd ..; \
    fi

RUN echo ${CMAKE_ARGS}

RUN mkdir -p build && \
    cd build && \
    rm -f -- CMakeCache.txt && \
    cmake .. ${CMAKE_ARGS} && \
    make -j

# Default to a login shell
CMD ["bash", "-l"]

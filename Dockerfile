FROM ubuntu:18.04
CMD bash

ENV TZ=America/Los_Angeles

# Install Ubuntu packages.
# Please add packages in alphabetical order.
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -y update && \
    apt-get -y install \
      build-essential \
      clang-8 \
      clang-format-8 \
      clang-tidy-8 \
      cmake \
      doxygen \
      git \
      g++-7 \
      pkg-config \
      tzdata \
      valgrind \
      zlib1g-dev && \
    date
	
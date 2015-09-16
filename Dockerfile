FROM base/archlinux

# to avoid keyring problem
RUN sed -i "s/SigLevel    = Required DatabaseOptional/SigLevel = Never/" /etc/pacman.conf
RUN yes | pacman -Sy
RUN yes | pacman -S \
      gcc \
      glibc \
      libunistring \
      make \
      cmake \
      eigen \
      protobuf \
      gtest \
      gflags \
      google-glog \
      boost \
      yaml-cpp

ADD . /tmp/uvlm
RUN [ -d /tmp/uvlm/build ] || mkdir -p /tmp/uvlm/build
RUN cd /tmp/uvlm/build && cmake .. && make -j && make install
RUN cd /tmp/uvlm/build && make test

VOLUME ["/data"]

CMD ["/bin/bash"]

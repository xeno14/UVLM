FROM base/archlinux

RUN pacman -Syy
# to avoid keyring problem
RUN sed -i "s/SigLevel    = Required DatabaseOptional/SigLevel = Never/" /etc/pacman.conf
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
      google-glog

ADD . /tmp/uvlm
RUN cd /tmp/uvlm && cmake . && make -j && make install

VOLUME ["/data"]

CMD ["/bin/bash"]

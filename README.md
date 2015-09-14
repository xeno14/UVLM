UVLM
====

# Dependency
- [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
- [googletest](http://opencv.jp/googletestdocs/primer.html)
- [gflags](https://github.com/gflags/gflags)
- [glog](https://github.com/google/glog)
- [googletest](https://github.com/google/googletest)
- [protobuf](https://github.com/google/protobuf)
- [yaml-cpp](https://github.com/jbeder/yaml-cpp)

# Run in Docker

- Build the image

```
$ docker buikd -t uvlm .
```

- Run the image with specifing `--volume`

```
docker run -ti --volume /path/to/disk:/data uvlm
```

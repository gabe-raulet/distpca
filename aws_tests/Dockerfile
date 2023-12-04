FROM amazonlinux:latest

RUN yum -y groupinstall "Development tools"
RUN yum -y install gcc-c++ libcurl-devel cmake3 git
RUN cd /opt && curl "https://awscli.amazonaws.com/awscli-exe-linux-aarch64.zip" -o "awscliv2.zip" && unzip awscliv2.zip && ./aws/install && cd ..
RUN cd /opt && git clone https://github.com/awslabs/aws-lambda-cpp.git && cd aws-lambda-cpp && mkdir build && cd build && \
    cmake3 .. -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF -DCMAKE_INSTALL_PREFIX=/opt/out -DCMAKE_CXX_COMPILER=g++ && make && make install

# Use the official CentOS 7 image as the base
FROM centos:7

# Install necessary dependencies
RUN yum -y install epel-release && \
    yum -y update && \
    yum -y install R R-devel && \
    yum -y install libcurl-devel openssl-devel libxml2-devel && \
    yum -y install libgit2-devel libssh2-devel && \
    yum -y install pandoc pandoc-citeproc && \
    yum -y install bzip2-devel cmake cmake3 make git libcurl wget && \
    yum -y install gcc gcc-c++ && \
    yum clean all

# RUN useradd -m user
# USER user

# Dependencies for pkgdown
RUN \
    yum -y install harfbuzz-devel fribidi-devel && \
    yum clean all
RUN \
    yum -y install freetype-devel libpng-devel libtiff-devel libjpeg-turbo-devel &&\
    yum clean all

RUN \
  cd /root && \
  git clone https://github.com/pezmaster31/bamtools.git && \
  mkdir bamtools/build && \
  cd bamtools/build && \
  cmake3 .. && make && make install

RUN \
  mkdir ~/local && \
  wget --no-check-certificate http://sourceforge.net/projects/boost/files/boost/1.60.0/boost_1_60_0.tar.gz && \
  tar -xvzf boost_1_60_0.tar.gz && \
  cd boost_1_60_0 && \
  ./bootstrap.sh --with-libraries=filesystem,iostreams,log,system,thread,test && \
  ./b2 cxxflags=-std=c++11 link=shared threading=multi install --prefix=/root/local

# Set a working directory
WORKDIR /root

RUN \
  cd && \
  git clone https://github.com/hms-dbmi/dropEst.git

# Install devtools package using R's package installer
RUN \ 
    R -e "install.packages('devtools', repos='https://cran.asia/')" && \
    R -e 'library(devtools)'

# Install dependencies
RUN \
    R -e 'install.packages(c("RcppProgress", "Rcpp", "RInside", "Matrix", "optparse", "rmarkdown"), dependencies=T, repos = "https://cran.asia/")' && \
    R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/RcppEigen/RcppEigen_0.3.3.9.1.tar.gz", repos=NULL, type="source")' && \
    R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/kernlab/kernlab_0.9-25.tar.gz", repos=NULL, type="source")' && \
    R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/mvtnorm/mvtnorm_1.0-11.tar.gz", repos=NULL, type="source")' && \
    R -e 'install.packages("ks", repos = "https://cran.asia/")' && \
    R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/pcaPP/pcaPP_1.9-74.tar.gz", repos=NULL, type="source")'

# Install dropestr
RUN \
    # R -e 'devtools::install_local("~/dropEst/dropestr/", dependencies=T)'
    R -e 'install.packages("preseqR", repos = "https://cran.asia/")' && \
    R -e 'devtools::install_github("kharchenkolab/dropestr")'

# Install dropest
RUN \
  mkdir ~/dropEst/build && \
  cd ~/dropEst/build && \
  cmake -D BOOST_ROOT=~/local/ .. && \
  make

ENV PATH=/root/dropEst/build:$PATH

# Start a bash shell by default when the container is run
# CMD ["bash"]

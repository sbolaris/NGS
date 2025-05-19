FROM ubuntu:18.04

LABEL Bio-Rad Support <support@bio-rad.com>

RUN apt-get update && apt-get install -y \
    curl \
    unzip \
    perl \
    parallel \
    pigz \
    default-jdk \
    wget \
    samtools \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    python3

#############STAR##############

ENV STAR_VERSION 2.7.0f
# Same deal as above with FASTQC; using precompiled executable
RUN curl -SLO https://github.com/alexdobin/STAR/archive/${STAR_VERSION}.tar.gz && tar -zxvf ${STAR_VERSION}.tar.gz --directory /opt/ && rm ${STAR_VERSION}.tar.gz
ENV PATH /opt/STAR-${STAR_VERSION}/bin/Linux_x86_64:$PATH

###############################

####### Python ################
RUN apt-get install -y \
    python3-pip
RUN pip3 install --upgrade pip setuptools
RUN pip3 install pandas
RUN pip3 install Bio
RUN pip3 install ray

###############################

######## Timezone for R #######

RUN export DEBIAN_FRONTEND=noninteractive
RUN ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime
RUN apt-get install -y tzdata
RUN dpkg-reconfigure --frontend noninteractive tzdata
###############################

##########  R  ################
RUN echo "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/" > /etc/apt/sources.list.d/cran.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9

ENV R_BASE_VERSION 4.0.2
RUN apt-get update \
        && apt-get install -y --no-install-recommends \
                littler \
                libcurl4-openssl-dev \
                r-cran-littler \
                r-base=${R_BASE_VERSION}* \
                r-base-dev=${R_BASE_VERSION}* \
                r-recommended=${R_BASE_VERSION}* \
                && echo 'options(repos = c(CRAN = "https://cloud.r-project.org/"), download.file.method = "libcurl")' >> /etc/R/Rprofile.site \
        && echo 'source("/etc/R/Rprofile.site")' >> /etc/littler.r
RUN apt-get install -y libxml2-dev
RUN Rscript -e 'install.packages("XML", repos = "http://www.omegahat.net/R")'
RUN Rscript -e 'install.packages("restfulr", repos = "http://cran.r-project.org")'
RUN Rscript -e 'install.packages("BiocManager");BiocManager::install(c("rtracklayer"))'
###############################
WORKDIR /opt/biorad
COPY . .
# Pull in some ARGS for defining container name
ARG IMAGE_NAME
ARG SOURCE_BRANCH
ARG SOURCE_COMMIT
RUN printf "Container Name: ${IMAGE_NAME:-local}\n" > imageInfo.txt
RUN printf "Source Branch: ${SOURCE_BRANCH:-local}\n" >> imageInfo.txt
RUN printf "Source Commit: ${SOURCE_COMMIT:-local}\n" >> imageInfo.txt


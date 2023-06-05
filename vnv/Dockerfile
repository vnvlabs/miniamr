
ARG FROM_IMAGE=vnv_base
FROM ${FROM_IMAGE}

ARG COMMIT_HASH=unknown
LABEL vnv.version.miniamr=${COMMIT_HASH}



#Install MINIAMR with VnV
COPY . ${SOFTWARE_DIR}/miniAMR
RUN cd ${SOFTWARE_DIR}/miniAMR/ref \
    && bear make \
    && ${VNV_DIR}/bin/vnv-matcher --package MINIAMR \
                                  --output vnv_MINIAMR.c\
                         		  --fix-omp \
                                  --regFile vnv.__registration__\
				                  --targetFile ${PWD}/miniAMR.x \
				                  --cache vnv.__cache__\
                                  compile_commands.json\
    && make 

ENV VNV_CONFIG=${VNV_CONFIG}:${SOFTWARE_DIR}/miniAMR/ref/vnv.__registration__
ENV VNV_CONFIG=${VNV_CONFIG}:${SOFTWARE_DIR}/miniAMR/vnv/gui/vnv.__registration__
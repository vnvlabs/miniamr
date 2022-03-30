
ARG FROM_IMAGE=vnv_base
FROM ${FROM_IMAGE}:latest



#Install MINIAMR with VnV
COPY . ./miniAMR
RUN cd ./miniAMR/ref && bear make 
RUN cd ./miniAMR/ref && ${VNV_DIR}/bin/vnv-matcher --package MINIAMR --output vnv_MINIAMR.c --cache vnv.__cache__ compile_commands.json
RUN cd ./miniAMR/ref && make 
RUN cd ./miniAMR/ref && ./miniAMR.x

ARG COMMIT_HASH=unknown
LABEL vnv.version.miniamr=${COMMIT_HASH}

ENV VNV_CONFIG=${VNV_CONFIG}:/miniAMR/.vnv-config.json


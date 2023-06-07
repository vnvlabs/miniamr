
ARG FROM_IMAGE=vnv_base
FROM ${FROM_IMAGE}

#Install MINIAMR with VnV
COPY . ${SOFTWARE_DIR}/miniAMR
RUN cd ${SOFTWARE_DIR}/miniAMR/ref \
    && bear make \
    && ${VNV_DIR}/bin/vnv-matcher --package MINIAMR \
                                  --output vnv_MINIAMR.c\
                         		  --fix-omp \
				                  --targetFile ${PWD}/miniAMR.x \
				                  --cache vnv.__cache__\
                                  compile_commands.json\
    && make 
    
RUN ${VNV_REGISTER} gui miniAMR ${SOFTWARE_DIR}/miniAMR/vnv/gui/vnv.__registration__ 

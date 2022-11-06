FROM chaff800/cmake_pipeline:latest

COPY . /tmp

WORKDIR /tmp

RUN if [ -d build ]; then rm -rf build; fi; mkdir build
RUN cmake --preset GTEST -B build -S .
RUN cmake --build build -- -j

CMD ["./build/MPMDCC_exec"]
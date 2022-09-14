# syntax=docker/dockerfile:1

FROM python:3.7
WORKDIR /usr/src/app
RUN pip install sctriangulate==0.12.0 \
    && mkdir run \
    && cd run \
    && curl -O input.h5ad https://raw.githubusercontent.com/frankligy/scTriangulate/main/docker/input.h5ad \
    && cd .. \
    && chmod -R 777 ./
ENTRYPOINT ['python /usr/src/app/run_sctriangulate.py']

# enter the Dockerfile folder, docker build --no-cache -t sctriangulate .
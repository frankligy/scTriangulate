# syntax=docker/dockerfile:1

FROM python:3.7
WORKDIR /usr/src/app
RUN pip install sctriangulate==0.12.0 \
    && curl -o run_sctriangulate.py https://raw.githubusercontent.com/frankligy/scTriangulate/main/docker/run_sctriangulate.py \
    && mkdir run \
    && chmod -R 777 ./
ENTRYPOINT ['/usr/src/app/run_sctriangulate.py']


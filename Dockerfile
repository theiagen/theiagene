ARG THEIAGENE_VER="1.0.0"

### start of app stage ###
FROM python:3.12-slim AS app

ARG THEIAGENE_VER

LABEL base.image="python:3.12-slim"
LABEL dockerfile.version="1"
LABEL software="theiagene"
LABEL software.version="${THEIAGENE_VER}"
LABEL description="Gene coverage and variant annotation toolkit for Theiagen bioinformatics workflows"
LABEL website="https://github.com/theiagen/theiagene"
LABEL license="https://github.com/theiagen/theiagene/blob/main/LICENSE"
LABEL maintainer="Zachary Konkel"
LABEL maintainer.email="zachary.konkel@theiagen.com"

ENV LC_ALL=C

# install the theiagene package (pulls in pysam + biopython)
WORKDIR /theiagene
COPY pyproject.toml README.md ./
COPY src/ ./src/
RUN python3 -m pip install --no-cache-dir .

# final working directory is /data
WORKDIR /data

# default command prints the help menu
CMD ["theiagene", "--help"]

### start of test stage ###
FROM app AS test

WORKDIR /theiagene
COPY conftest.py ./
COPY tests/ ./tests/

RUN python3 -m pip install --no-cache-dir pytest

RUN theiagene --help && pytest -q

FROM continuumio/miniconda3:25.1.1-2
LABEL version="0.0.1" \
      description="Docker image for viralunity"

WORKDIR /app

COPY envs/vu_dependencies.yml /app/envs/vu_dependencies.yml
# Install the environment
RUN conda env create --quiet -f envs/vu_dependencies.yml && conda clean -a -y

ENV PATH=/opt/conda/envs/viralunity/bin:$PATH

COPY . /app/viralunity

WORKDIR /app/viralunity
RUN pip install . && rm -rf /root/.cache/pip
RUN viralunity_meta -v > viralunity-version.txt

WORKDIR /tmp/
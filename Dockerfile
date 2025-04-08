FROM continuumio/miniconda3:25.1.1-2
LABEL version="0.0.1" \
      description="Docker image for viralunity"

WORKDIR /app

COPY environment.yml /app/environment.yml
# Install the environment
RUN conda env create --quiet -f environment.yml && conda clean -a -y

ENV PATH=/opt/conda/envs/viralunity/bin:$PATH

COPY . /app/viralunity

WORKDIR /app/viralunity
RUN pip install . && rm -rf /root/.cache/pip
RUN viralunity -v > viralunity-version.txt

WORKDIR /tmp/
# FROM rocker/tidyverse

# RUN apt-get update && apt-get install -y \
#   build-essential \
#   libglpk40

# RUN R -e "install.packages('igraph', dependencies = T)"
# RUN R -e "install.packages('ergm', dependencies = T)"
# RUN R -e "install.packages('intergraph', dependencies = T)"
# RUN R -e "install.packages('sbm', dependencies = T)"
# RUN R -e "install.packages('languageserver', repos='https://p3m.dev/cran/__linux__/jammy/latest', dependencies = T)"



# CMD R -e "source('script/igraph.R')"

ARG VARIANT="4.3"
FROM ghcr.io/rocker-org/devcontainer/tidyverse:${VARIANT}

RUN apt-get update && apt-get install -y \
  build-essential \
  libglpk40

RUN install2.r --error --skipinstalled -n -1 \
        statip \
        patchwork \
        paletteer \
        here \
        doParallel \
        janitor \
        vip \
        ranger \
        palmerpenguins \
        skimr \
        nnet \
        kernlab \
        plotly \
        factoextra \
        cluster \
        tidymodels \
        markdown \
        ottr \
    && rm -rf /tmp/downloaded_packages \
    && R -q -e 'remotes::install_github("https://github.com/dcomtois/summarytools/tree/0-8-9")'

RUN R -e "install.packages('igraph', dependencies = T)" \
&& R -e "install.packages('ergm', dependencies = T)" \
&& R -e "install.packages('intergraph', dependencies = T)" \
&& R -e "install.packages('sbm', dependencies = T)" \
&& R -e "install.packages('languageserver', repos='https://p3m.dev/cran/__linux__/jammy/latest', dependencies = T)"

# Install Python packages
COPY requirements.txt /tmp/pip-tmp/
RUN python3 -m pip --disable-pip-version-check --no-cache-dir install -r /tmp/pip-tmp/requirements.txt
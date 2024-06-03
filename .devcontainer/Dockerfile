FROM rocker/tidyverse

RUN apt-get update && apt-get install -y \
  build-essential \
  libglpk40

RUN R -e "install.packages('igraph', dependencies = T)"
RUN R -e "install.packages('ergm', dependencies = T)"
RUN R -e "install.packages('intergraph', dependencies = T)"
RUN R -e "install.packages('sbm', dependencies = T)"
RUN R -e "install.packages('languageserver', repos='https://p3m.dev/cran/__linux__/jammy/latest', dependencies = T)"



CMD R -e "source('script/igraph.R')"
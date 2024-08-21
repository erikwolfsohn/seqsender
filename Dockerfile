# Create an argument to pull a particular version of micromamba image
ARG micromamba_version
ARG micromamba_version=${micromamba_version:-1.5.3}

############# base image ##################
FROM --platform=$BUILDPLATFORM google/cloud-sdk:489.0.0-stable AS base

# local apt mirror support
# start every stage with updated apt sources
ARG APT_MIRROR_NAME=
RUN if [ -n "$APT_MIRROR_NAME" ]; then sed -i.bak -E '/security/! s^https?://.+?/(debian|ubuntu)^http://'"$APT_MIRROR_NAME"'/\1^' /etc/apt/sources.list && grep '^deb' /etc/apt/sources.list; fi
RUN apt-get update --allow-releaseinfo-change --fix-missing

############# micromamba image ##################

FROM --platform=$BUILDPLATFORM mambaorg/micromamba:${micromamba_version} AS micromamba
RUN echo "Getting micromamba image"

############# Build Stage: Final ##################

# Build the final image 
FROM base AS final

# if image defaults to a non-root user, then we may want to make the
# next 3 ARG commands match the values in our image. 
ENV MAMBA_USER=$MAMBA_USER
ENV MAMBA_USER_ID=$MAMBA_USER_ID
ENV MAMBA_USER_GID=$MAMBA_USER_GID
ENV MAMBA_ROOT_PREFIX="/opt/conda"
ENV MAMBA_EXE="/bin/micromamba"

COPY --from=micromamba "$MAMBA_EXE" "$MAMBA_EXE"
COPY --from=micromamba /usr/local/bin/_activate_current_env.sh /usr/local/bin/_activate_current_env.sh
COPY --from=micromamba /usr/local/bin/_dockerfile_shell.sh /usr/local/bin/_dockerfile_shell.sh
COPY --from=micromamba /usr/local/bin/_entrypoint.sh /usr/local/bin/_entrypoint.sh
COPY --from=micromamba /usr/local/bin/_dockerfile_initialize_user_accounts.sh /usr/local/bin/_dockerfile_initialize_user_accounts.sh
COPY --from=micromamba /usr/local/bin/_dockerfile_setup_root_prefix.sh /usr/local/bin/_dockerfile_setup_root_prefix.sh

# Install system dependencies
ARG DEBIAN_FRONTEND=noninteractive

# Install system libraries of general use
RUN apt-get update --allow-releaseinfo-change --fix-missing \
  && apt-get install --no-install-recommends -y \
  dos2unix \
  ca-certificates \
  git \
  && apt clean autoclean \
  && apt autoremove --yes \
  && rm -rf /var/lib/{apt,dpkg,cache,log}/

ARG REPO_URL="https://github.com/erikwolfsohn/seqsender.git"
ARG REPO_BRANCH="seqsender-terra-submission"

RUN git clone --branch ${REPO_BRANCH} ${REPO_URL} /seqsender

# Create working directory 
ENV WORKDIR=/seqsender

# Set up volume directory
#VOLUME ${WORKDIR}

# Set up working directory
WORKDIR ${WORKDIR}

# Allow permission to read and write to working directory
RUN chmod -R a+rwx ${WORKDIR}

# Create a program variable
#ENV PROJECT_DIR=/seqsender

# Set up a volume directory 
#VOLUME ${PROJECT_DIR}

# Copy all files to project directory
#COPY . ${PROJECT_DIR}

############ Set-up micromamba environment ##################

# Copy requirement files to program directory
#COPY env.yaml "${PROJECT_DIR}/env.yaml"

# Set up environments
RUN micromamba install --yes --name base -f "${WORKDIR}/env.yaml" \
    && micromamba clean --all --yes

############# Launch PROGRAM ##################

# Copy bash script to run PROGRAM to docker image
# COPY seqsender-kickoff "${PROJECT_DIR}/seqsender-kickoff"

# Convert bash script from Windows style line endings to Unix-like control characters
RUN dos2unix "${WORKDIR}/seqsender-kickoff"

# Allow permission to excute the bash script
# RUN chmod a+x "${PROJECT_DIR}/seqsender-kickoff"

# Allow permission to read and write to program directory
RUN chmod a+rwx ${WORKDIR}

# Export bash script to path
ENV PATH="$PATH:${WORKDIR}"

# Activate conda environment
ENV PATH="$PATH:${MAMBA_ROOT_PREFIX}/bin"

# Execute the pipeline 
# ENTRYPOINT ["tail", "-f", "/dev/null"]

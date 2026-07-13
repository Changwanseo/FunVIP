# FunVIP in a self-contained conda environment (bundles the external tools).
# Useful on Windows/WSL or anywhere a native install is awkward.
#
#   docker build -t funvip .
#   docker run --rm -v "$PWD:/data" funvip --test terrei --email you@example.com --outdir /data/out
#
FROM condaforge/miniforge3:latest

# Build the environment from environment.yml (external tools via conda + FunVIP
# via pip). Clean caches to keep the image small.
COPY environment.yml /tmp/environment.yml
RUN mamba env create -f /tmp/environment.yml && mamba clean -a -y

# Headless Qt for ete4 tree rendering
ENV QT_QPA_PLATFORM=offscreen

# Run FunVIP inside the env by default; arguments passed to `docker run` go to FunVIP.
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "funvip", "FunVIP"]
CMD ["--help"]

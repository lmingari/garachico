# ---------- Stage 1: Build ----------
FROM fall3d-wgrib2:0.1 AS builder

WORKDIR /build

# Install building dependencies:
RUN apk add --no-cache \
    gcc \
    g++ \
    proj-util \
    netcdf-dev \
    hdf5-dev \
    musl-dev \
    gdal-dev \
    geos-dev \
    proj-dev \
    python3-dev \
    && rm -rf /var/cache/apk/*

COPY requirements.txt .

# Create virtual env
RUN python3 -m venv /build/venv

# Set PATH so executables are usable
ENV PATH="/build/venv/bin:${PATH}"

# Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# ---------- Stage 2: Executables ----------
FROM fall3d-wgrib2:0.1

WORKDIR /app

# Copy the virtual environment from the build stage
COPY --from=builder /build/venv /app/venv

ENV PATH="/app/venv/bin:${PATH}"

# Install runtime dependencies:
RUN apk add --no-cache \
    gdal \
    && rm -rf /var/cache/apk/*

# Copy static files
COPY cartopy_init.py .
COPY createCOG.py .
COPY fill_template.py .
COPY gfs.py .
COPY gfs_times.py .
COPY plotter.py .
COPY template.inp .
COPY GFS.tbl .
COPY logo.png .

RUN python cartopy_init.py

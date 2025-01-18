#!/bin/bash -e

# Define constants
CODECOV_BASE_URL="https://uploader.codecov.io/latest/linux"
PGP_KEYS_URL="https://keybase.io/codecovsecurity/pgp_keys.asc"

# Download and import the PGP keys
echo "Importing Codecov PGP keys..."
curl -fsSL "$PGP_KEYS_URL" | gpg --no-default-keyring --keyring trustedkeys.gpg --import

# Download Codecov uploader and verification files
echo "Downloading Codecov uploader and verification files..."
curl -fsSLO "$CODECOV_BASE_URL/codecov"
curl -fsSLO "$CODECOV_BASE_URL/codecov.SHA256SUM"
curl -fsSLO "$CODECOV_BASE_URL/codecov.SHA256SUM.sig"

# Verify the downloaded files
echo "Verifying signature of SHA256SUM file..."
gpgv codecov.SHA256SUM.sig codecov.SHA256SUM

echo "Checking SHA256 checksum..."
shasum -a 256 -c codecov.SHA256SUM

# Make the Codecov uploader executable
echo "Setting executable permissions for the Codecov uploader..."
chmod +x codecov

# Run the Codecov uploader with provided arguments
echo "Running the Codecov uploader..."
./codecov "$@"

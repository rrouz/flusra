#!/bin/bash

# Function to generate MD5 checksum for a given file using openssl
generate_md5_checksum() {
    openssl md5 -binary "$1" | openssl base64
}

# Function to compare checksums of files in expected and observed directories
compare_checksums() {
    expected_dir="$1"
    observed_dir="$2"

    while IFS= read -r -d '' file_name; do
        expected_file_path="$file_name"
        # remove the expected_dir prefix from the file path
        file_name="${file_name#$expected_dir/}"
        observed_file_path="$observed_dir/$file_name"

        if [ -f "$expected_file_path" ] && [ -f "$observed_file_path" ]; then
            expected_checksum=$(generate_md5_checksum "$expected_file_path")
            observed_checksum=$(generate_md5_checksum "$observed_file_path")

            if [ "$expected_checksum" != "$observed_checksum" ]; then
                echo "Checksum mismatch for file: $file_name"
                echo "Expected: $expected_checksum"
                echo "Observed: $observed_checksum"
            else
                echo "Checksum matched for file: $file_name"
            fi
        else
            echo "File not found: $file_name at $expected_file_path or $observed_file_path"
        fi
    done < <(find "$expected_dir" -type f -print0)
}

compare_checksums "test/output" "testing/output"
